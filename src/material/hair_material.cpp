#include "reflection.h"
#include "textures/constant.h"
#include "spectrum.h"
#include "texture.h"
#include "paramset.h"
#include "interaction.h"

#include <fstream>
#include <string>
#include <iostream>
#include <random>

//hair

#include <array>
#include <numeric>
#include "sampling.h"


#include"../table/brdf_c5.h" //キューティクル多層
//#include"../table/brdf_nc5.h" //キューティクルなし
//#include"../tablebrdf_hk5.h" //キューティクル剥落
//#include"../table/brdf_ck5.h"//キューティクル欠け
//#include"../table/brdf_hkck5.h"//キューティクル欠け+剥落

#include"morpho.h"
#include<typeinfo>

namespace pbrt {

// Hair Local Declarations

inline Float _I0(Float x), _LogI0(Float x);

// Hair Local Functions
static Float _Mp(Float cosThetaI, Float cosThetaO, Float sinThetaI,
                Float sinThetaO, Float v) {
    Float a = cosThetaI * cosThetaO / v;
    Float b = sinThetaI * sinThetaO / v;
    Float mp =
        (v <= .1)
            ? (std::exp(_LogI0(a) - b - 1 / v + 0.6931f + std::log(1 / (2 * v))))
            : (std::exp(-b) * _I0(a)) / (std::sinh(1 / v) * 2 * v);
    CHECK(!std::isinf(mp) && !std::isnan(mp));
    return mp;
}

 

inline Float _I0(Float x) {
    Float val = 0;
    Float x2i = 1;
    int64_t ifact = 1;
    int i4 = 1;
    // I0(x) \approx Sum_i x^(2i) / (4^i (i!)^2)
    for (int i = 0; i < 10; ++i) {
        if (i > 1) ifact *= i;
        val += x2i / (i4 * _Sqr(ifact));
        x2i *= x * x;
        i4 *= 4;
    }
    return val;
}

inline Float _LogI0(Float x) {
    if (x > 12)
        return x + 0.5 * (-std::log(2 * Pi) + std::log(1 / x) + 1 / (8 * x));
    else
        return std::log(_I0(x));
}

static std::array<Spectrum, _pMax + 1> _Ap(Float cosThetaO, Float eta, Float h,
                                         const RGBSpectrum &T) {
    std::array<Spectrum, _pMax + 1> ap;
    // Compute $p=0$ attenuation at initial cylinder intersection
    Float cosGammaO = _SafeSqrt(1 - h * h);
    Float cosTheta = cosThetaO * cosGammaO;
    Float f = FrDielectric(cosTheta, 1.f, eta); 
    ap[0] = f; //test16
    
    // Compute $p=1$ attenuation term
    ap[1] = _Sqr(1 - f) * T;
    //ap[1] = T;

    // Compute attenuation terms up to $p=__pMax_$
    for (int p = 2; p < _pMax; ++p) ap[p] =  ap[p - 1] * T * f;

    // Compute attenuation term accounting for remaining orders of scattering
    ap[_pMax] =  ap[_pMax - 1] * f * T / (Spectrum(1.f) - T * f); //T
    return ap;
}

inline Float _Phi(int p, Float gammaO, Float gammaT) {
    return 2 * p * gammaT - 2 * gammaO + p * Pi;
}

inline Float _Logistic(Float x, Float s) {
    x = std::abs(x);
    return std::exp(-x / s) / (s * _Sqr(1 + std::exp(-x / s)));
}

inline Float _LogisticCDF(Float x, Float s) {
    return 1 / (1 + std::exp(-x / s));
}

inline Float _TrimmedLogistic(Float x, Float s, Float a, Float b) {
    CHECK_LT(a, b);
    return _Logistic(x, s) / (_LogisticCDF(b, s) - _LogisticCDF(a, s));
}

inline Float _Np(Float phi, int p, Float s, Float gammaO, Float gammaT) {
    Float dphi = phi - _Phi(p, gammaO, gammaT);
    // Remap _dphi_ to $[-\pi,\pi]$
    while (dphi > Pi) dphi -= 2 * Pi;
    while (dphi < -Pi) dphi += 2 * Pi;
    return _TrimmedLogistic(dphi, s, -Pi, Pi);
}

static Float _SampleTrimmedLogistic(Float u, Float s, Float a, Float b) {
    CHECK_LT(a, b);
    Float k = _LogisticCDF(b, s) - _LogisticCDF(a, s);
    Float x = -s * std::log(1 / (u * k + _LogisticCDF(a, s)) - 1);
    CHECK(!std::isnan(x));
    return Clamp(x, a, b);
}


// HairMaterial Method Definitions
void MorphoMaterial::ComputeScatteringFunctions(SurfaceInteraction *si,
                                              MemoryArena &arena,
                                              TransportMode mode,
                                              bool allowMultipleLobes) const {
    Float bm = beta_m->Evaluate(*si);
    Float bn = beta_n->Evaluate(*si);
    Float a = alpha->Evaluate(*si);
    Float e = eta->Evaluate(*si);

    si->bsdf = ARENA_ALLOC(arena, BSDF)(*si, e);

    RGBSpectrum sig_a;
    if (sigma_a)
        sig_a = sigma_a->Evaluate(*si).Clamp();
    else if (color) {
        RGBSpectrum c = color->Evaluate(*si).Clamp();
        //c = sig_a; テストコード
        sig_a = MorphoBSDF::SigmaAFromReflectance(c, bn);
    } else {
        CHECK(eumelanin || pheomelanin);
        sig_a = MorphoBSDF::SigmaAFromConcentration(
            std::max(Float(0), eumelanin ? eumelanin->Evaluate(*si) : 0),
            std::max(Float(0), pheomelanin ? pheomelanin->Evaluate(*si) : 0));
    }

    // Offset along width
    Float h = -1 + 2 * si->uv[1];
    int wavelengthindex = 0;// = si->wavelengthindex; 
    si->bsdf->Add(ARENA_ALLOC(arena, MorphoBSDF)(h, e, sig_a, bm, bn, a, wavelengthindex));
}

MorphoMaterial *CreateMorphoMaterial(const TextureParams &mp) {


    std::shared_ptr<Texture<RGBSpectrum>> sigma_a =
        mp.GetSpectrumTextureOrNull("sigma_a");
    std::shared_ptr<Texture<RGBSpectrum>> color =
        mp.GetSpectrumTextureOrNull("color");
    std::shared_ptr<Texture<Float>> eumelanin =
        mp.GetFloatTextureOrNull("eumelanin");
    std::shared_ptr<Texture<Float>> pheomelanin =
        mp.GetFloatTextureOrNull("pheomelanin");
    if (sigma_a) {
        if (color)
            Warning(
                "Ignoring \"color\" parameter since \"sigma_a\" was provided.");
        if (eumelanin)
            Warning(
                "Ignoring \"eumelanin\" parameter since \"sigma_a\" was "
                "provided.");
        if (pheomelanin)
            Warning(
                "Ignoring \"pheomelanin\" parameter since \"sigma_a\" was "
                "provided.");
    } else if (color) {
        if (sigma_a)
            Warning(
                "Ignoring \"sigma_a\" parameter since \"color\" was provided.");
        if (eumelanin)
            Warning(
                "Ignoring \"eumelanin\" parameter since \"color\" was "
                "provided.");
        if (pheomelanin)
            Warning(
                "Ignoring \"pheomelanin\" parameter since \"color\" was "
                "provided.");
    } else if (eumelanin || pheomelanin) {
        if (sigma_a)
            Warning(
                "Ignoring \"sigma_a\" parameter since "
                "\"eumelanin\"/\"pheomelanin\" was provided.");
        if (color)
            Warning(
                "Ignoring \"color\" parameter since "
                "\"eumelanin\"/\"pheomelanin\" was provided.");
    } else {
        // Default: brown-ish hair.
        sigma_a = std::make_shared<ConstantTexture<RGBSpectrum>>(
            MorphoBSDF::SigmaAFromConcentration(1.3, 0.));
    }

    std::shared_ptr<Texture<Float>> eta = mp.GetFloatTexture("eta", 1.55f);
    std::shared_ptr<Texture<Float>> beta_m = mp.GetFloatTexture("beta_m", 0.3f);
    std::shared_ptr<Texture<Float>> beta_n = mp.GetFloatTexture("beta_n", 0.3f);
    std::shared_ptr<Texture<Float>> alpha = mp.GetFloatTexture("alpha", 2.f);

    return new MorphoMaterial(sigma_a, color, eumelanin, pheomelanin, eta, beta_m,
                            beta_n, alpha);
}

// HairBSDF Method Definitions
MorphoBSDF::MorphoBSDF(Float h, Float eta, const RGBSpectrum &sigma_a, Float beta_m,
                   Float beta_n, Float alpha, int wavelengthindex)
    : BxDF(BxDFType(BSDF_GLOSSY | BSDF_REFLECTION | BSDF_TRANSMISSION)),
      h(h),
      gammaO(_SafeASin(h)),
      eta(eta),
      sigma_a(sigma_a),
      beta_m(beta_m),
      beta_n(beta_n),
      wavelengthindex(wavelengthindex) {
    CHECK(h >= -1 && h <= 1);
    CHECK(beta_m >= 0 && beta_m <= 1);
    CHECK(beta_n >= 0 && beta_n <= 1);
    // Compute longitudinal variance from $\beta_m$
    static_assert(
        _pMax >= 1,
        "Longitudinal variance code must be updated to handle low _pMax");
    v[0] = _Sqr(0.726f * beta_m + 0.812f * _Sqr(beta_m) + 3.7f * _Pow<20>(beta_m));
    v[1] = .25 * v[0];
    v[2] = 4 * v[0];
    for (int p = 3; p <= _pMax; ++p)
        // TODO: is there anything better here?
        v[p] = v[2];

    // Compute azimuthal logistic scale factor from $\beta_n$
    s = _SqrtPiOver8 *
        (0.265f * beta_n + 1.194f * _Sqr(beta_n) + 5.372f * _Pow<22>(beta_n));
    CHECK(!std::isnan(s));

    // Compute $\alpha$ terms for hair scales
    sin2kAlpha[0] = std::sin(Radians(alpha));
    cos2kAlpha[0] = _SafeSqrt(1 - _Sqr(sin2kAlpha[0]));
    for (int i = 1; i < 3; ++i) {
        sin2kAlpha[i] = 2 * cos2kAlpha[i - 1] * sin2kAlpha[i - 1];
        cos2kAlpha[i] = _Sqr(cos2kAlpha[i - 1]) - _Sqr(sin2kAlpha[i - 1]);
    }
}

Spectrum MorphoBSDF::f(const Vector3f &wo, const Vector3f &wi) const {
    
    // Compute hair coordinate system terms related to _wo_
    Float sinThetaO = wo.x;
    Float cosThetaO = _SafeSqrt(1 - _Sqr(sinThetaO));
    Float phiO = std::atan2(wo.z, wo.y);
    Float thO = std::atan2(wo.z, wo.x);

    // Compute hair coordinate system terms related to _wi_
    Float sinThetaI = wi.x;
    Float cosThetaI = _SafeSqrt(1 - _Sqr(sinThetaI));
    Float phiI = std::atan2(wi.z, wi.y);
    Float thI = std::atan2(wi.z, wi.x);

    // Compute $\cos \thetat$ for refracted ray
    Float sinThetaT = sinThetaO / eta;
    Float cosThetaT = _SafeSqrt(1 - _Sqr(sinThetaT));

    // Compute $\gammat$ for refracted ray
    Float etap = std::sqrt(eta * eta - _Sqr(sinThetaO)) / cosThetaO;
    Float sinGammaT = h / etap;
    Float cosGammaT = _SafeSqrt(1 - _Sqr(sinGammaT));
    Float gammaT = _SafeASin(sinGammaT);
    float _sinThetaIp = sinThetaI * cos2kAlpha[1] + cosThetaI * sin2kAlpha[1];

    float _thetO = std::abs(std::asin(sinThetaO));
    float _thetI = std::abs(std::asin(sinThetaI));  
    float _phiI = std::abs(phiI);
    float _phiO = std::abs(phiO);
    float _gamT = std::abs(gammaT);
    float _thetIp = std::abs(std::asin(_sinThetaIp));

    int rphiI =int( _phiI * 180 / Pi);
    int rphiO = int(_phiO * 180 / Pi);
    int rgamT = int(_gamT * 180 / Pi);
    int rgamO = int(std::abs(gammaO) * 180 / Pi);
    int rthetO = int(_thetO * 180 / Pi);
    int rthetI = int(_thetI * 180 / Pi);


    //角度を算出 θi, θo
    auto wo_theta = std::atan2(wo.x, std::sqrt(wo.y * wo.y + wo.z * wo.z));
    auto wo_ang = wo_theta * 180 / Pi;
    int ot = std::abs(std::round(wo_ang));

    auto wi_theta = std::atan2(wi.x, std::sqrt(wi.y * wi.y + wi.z * wi.z));
    auto wi_ang = wi_theta * 180 / Pi;
    int it = std::abs(std::round(wi_ang));

    if(it > 90){
        it = it - 90;
    }

    if(ot < 0){
        ot = ot + 90;
    }

    
    if (rthetI > 90) {
        rthetI = int(abs(180 - rthetI));
    }
    if (rthetO > 90) {
        rthetO = int(abs(180 - rthetO));
    }

    if (rthetI < 0) {
        rthetI = int(abs(rthetI));
    }
    if (rthetO < 0) {
        rthetO = int(abs(rthetO));
    }

    //std::cout << "rt" << int(rthetI) << std::endl;
    int rthetIp = int(_thetIp * 180 / Pi);

    SampledSpectrum spectrum = 0;//data[std::round(rthetI)][std::round(rgamT)];
   
    auto ret = spectrum.ToRGBSpectrum();
    
    RGBSpectrum T1 = Exp(-sigma_a * (2 * cosGammaT / cosThetaT));


    /*
    int thetanum = 180;
    int thetanumm = 360;
    Float asino = _SafeASin(std::abs(sinThetaO));
    Float asinoi = _SafeASin(std::abs(sinThetaI));
    Float unit = Pi / 2.0 / (Float)thetanum;
    Float myU = 2*(Float)thetanum / Pi;
    int floortheta = Clamp(std::floor(asino / unit), -180, thetanum-1);
    int floorthetaI = Clamp(std::floor(asinoi / unit), -90, thetanum-1);
    int floorTheta = Clamp(std::floor(asino * myU), -180, thetanum-1);
    int floorThetaI = Clamp(std::floor(asinoi * myU), -90, thetanum-1);

    int Floortheta = Clamp(std::floor(thO * myU), 0, thetanum-1);
    int FloorthetaI = Clamp(std::floor(thI * myU), 0, thetanum-1);
    //std::cout << "floortheta" << floortheta << std::endl; 
    //std::cout << "floorthetaI" << floorthetaI << std::endl; 
*/

    
    

/   
    float sinThetaIp = sinThetaI * cos2kAlpha[1] + cosThetaI * sin2kAlpha[1];
    float cosThetaIp = cosThetaI * cos2kAlpha[1] - sinThetaI * sin2kAlpha[1];
    Float cosGammaO = _SafeSqrt(1 - h * h);
    Float cosTheta = cosThetaO * cosGammaO;
    Float fm = FrDielectric(cosTheta, 1.f, eta);
    
    cosThetaO = std::abs(cosThetaO);



    // Compute the transmittance _T_ of a single path through the cylinder
    RGBSpectrum T = Exp(-sigma_a * (2 * cosGammaT / cosThetaT));


    // Evaluate hair BSDF
    Float phi = phiI - phiO;
    std::array<Spectrum, _pMax + 1> ap = _Ap(cosThetaO, eta, h, T);
    
    //RGBSpectrum fsum_sub(0.);
    
    Spectrum fsum(0.);
    RGBSpectrum RN(0.);


    for(int i =0;i<60;++i){
                RN[i] = BRDFTABLE5[it][ot][i]/2.5;
                RN[i] *= SampledSpectrum::rgbIllum2SpectWhite[i];
        }
        fsum  = RN.ToRGBSpectrum();;

    
    for (int p = 0; p < _pMax; ++p) {
        // Compute $\sin \thetai$ and $\cos \thetai$ terms accounting for scales
        Float sinThetaIp, cosThetaIp;
        int floorthetaIp,floorThetaIp;
        int floorthetaT;
        

        if (p == 0) {
            sinThetaIp = sinThetaI * cos2kAlpha[1] + cosThetaI * sin2kAlpha[1];
            cosThetaIp = cosThetaI * cos2kAlpha[1] - sinThetaI * sin2kAlpha[1];
        }

        // Handle remainder of $p$ values for hair scale tilt
        else if (p == 1) {
            sinThetaIp = sinThetaI * cos2kAlpha[0] - cosThetaI * sin2kAlpha[0];
            cosThetaIp = cosThetaI * cos2kAlpha[0] + sinThetaI * sin2kAlpha[0];
        } else if (p == 2) {
            sinThetaIp = sinThetaI * cos2kAlpha[2] - cosThetaI * sin2kAlpha[2];
            cosThetaIp = cosThetaI * cos2kAlpha[2] + sinThetaI * sin2kAlpha[2];
        } else {
            sinThetaIp = sinThetaI;
            cosThetaIp = cosThetaI;
        }

        // Handle out-of-range $\cos \thetai$ from scale adjustment
        cosThetaIp = std::abs(cosThetaIp);
        //fsum += _Mp(cosThetaIp, cosThetaO, sinThetaIp, sinThetaO, v[p]) * ap[p] *_Np(phi, p, s, gammaO, gammaT);
    }
   
    // Compute contribution of remaining terms after __pMax_

    //fsum += _Mp(cosThetaI, cosThetaO, sinThetaI, sinThetaO, v[_pMax]) * ap[_pMax] / (2.f * Pi);

    //if (AbsCosTheta(wi) > 0) fsum /= AbsCosTheta(wi);
    //CHECK(!std::isinf(fsum.y()) && !std::isnan(fsum.y()));
    //fsum = RN;
    return fsum;
}




std::array<Float, _pMax + 1> MorphoBSDF::ComputeApPdf(Float cosThetaO, const Vector3f &wo) const {
    // Compute array of $A_p$ values for _cosThetaO_
    Float sinThetaO = _SafeSqrt(1 - cosThetaO * cosThetaO);

    // Compute $\cos \thetat$ for refracted ray
    Float sinThetaT = sinThetaO / eta;
    Float cosThetaT = _SafeSqrt(1 - _Sqr(sinThetaT));

    // Compute $\gammat$ for refracted ray
    Float etap = std::sqrt(eta * eta - _Sqr(sinThetaO)) / cosThetaO;
    Float sinGammaT = h / etap;
    Float cosGammaT = _SafeSqrt(1 - _Sqr(sinGammaT));

    Vector3f wm = Vector3f(-wo.x, -wo.y, wo.z);

    float th_O = std::asin(sinThetaO);
    float th_I = std::asin(sinThetaT);
    float rth_I = th_I * 180 / Pi;
    float rth_O = th_O * 180 / Pi;  


    int flag = 0;



    // Compute the transmittance _T_ of a single path through the cylinder

    //if(flag == 1) Spectrum T = Exp(-sigma_a * (2 * cosGammaT / cosThetaT));
    Spectrum T = Exp(-sigma_a * (2 * cosGammaT / cosThetaT));

   


    std::array<Spectrum, _pMax + 1> ap = _Ap(cosThetaO, eta, h, T);

    // Compute $A_p$ PDF from individual $A_p$ terms
    std::array<Float, _pMax + 1> apPdf;
    Float sumY =
        std::accumulate(ap.begin(), ap.end(), Float(0),
                        [](Float s, const Spectrum &ap) { return s + ap.y(); });
    for (int i = 0; i <= _pMax; ++i) apPdf[i] = ap[i].y() / sumY;
    return apPdf;
}

Spectrum MorphoBSDF::Sample_f(const Vector3f &wo, Vector3f *wi, const Point2f &u2,
                            Float *pdf, BxDFType *sampledType) const {
    // Compute hair coordinate system terms related to _wo_
    Float sinThetaO = wo.x;
    Float sinThetaT = sinThetaO / eta;
    Float cosThetaO = _SafeSqrt(1 - _Sqr(sinThetaO));
    Float phiO = std::atan2(wo.z, wo.y);

    // Derive four random samples from _u2_
    Point2f u[2] = {_DemuxFloat(u2[0]), _DemuxFloat(u2[1])};

    // Determine which term $p$ to sample for hair scattering
    std::array<Float, _pMax + 1> apPdf = ComputeApPdf(cosThetaO,wo);
    int p;
    for (p = 0; p < _pMax; ++p) {
        if (u[0][0] < apPdf[p]) break;
        u[0][0] -= apPdf[p];
    }

    // Sample $M_p$ to compute $\thetai$
    u[1][0] = std::max(u[1][0], Float(1e-5));
    Float cosTheta =
        1 + v[p] * std::log(u[1][0] + (1 - u[1][0]) * std::exp(-2 / v[p]));
    Float sinTheta = _SafeSqrt(1 - _Sqr(cosTheta));
    Float cosPhi = std::cos(2 * Pi * u[1][1]);
    Float sinThetaI = -cosTheta * sinThetaO + sinTheta * cosPhi * cosThetaO;
    Float cosThetaI = _SafeSqrt(1 - _Sqr(sinThetaI));

    // Update sampled $\sin \thetai$ and $\cos \thetai$ to account for scales
    Float sinThetaIp = sinThetaI, cosThetaIp = cosThetaI;
    if (p == 0) {
        sinThetaIp = sinThetaI * cos2kAlpha[1] - cosThetaI * sin2kAlpha[1];
        cosThetaIp = cosThetaI * cos2kAlpha[1] + sinThetaI * sin2kAlpha[1];
    } else if (p == 1) {
        sinThetaIp = sinThetaI * cos2kAlpha[0] + cosThetaI * sin2kAlpha[0];
        cosThetaIp = cosThetaI * cos2kAlpha[0] - sinThetaI * sin2kAlpha[0];
    } else if (p == 2) {
        sinThetaIp = sinThetaI * cos2kAlpha[2] + cosThetaI * sin2kAlpha[2];
        cosThetaIp = cosThetaI * cos2kAlpha[2] - sinThetaI * sin2kAlpha[2];
    }
    sinThetaI = sinThetaIp;
    cosThetaI = cosThetaIp;

    // Sample $N_p$ to compute $\Delta\phi$

    // Compute $\gammat$ for refracted ray
    Float etap = std::sqrt(eta * eta - _Sqr(sinThetaO)) / cosThetaO;
    Float sinGammaT = h / etap;
    Float gammaT = _SafeASin(sinGammaT);
    Float dphi;
    if (p < _pMax)
        dphi =
            _Phi(p, gammaO, gammaT) + _SampleTrimmedLogistic(u[0][1], s, -Pi, Pi);
    else
        dphi = 2 * Pi * u[0][1];

    // Compute _wi_ from sampled hair scattering angles
    Float phiI = phiO + dphi;
    *wi = Vector3f(sinThetaI, cosThetaI * std::cos(phiI),
                   cosThetaI * std::sin(phiI));

    // Compute PDF for sampled hair scattering direction _wi_
    *pdf = 0;
    
    

    
    for (int p = 0; p < _pMax; ++p) {
        // Compute $\sin \thetai$ and $\cos \thetai$ terms accounting for scales
        Float sinThetaIp, cosThetaIp;
        int floorthetaIp;
        int floorthetaT;
        if (p == 0) {
            sinThetaIp = sinThetaI * cos2kAlpha[1] + cosThetaI * sin2kAlpha[1];
            cosThetaIp = cosThetaI * cos2kAlpha[1] - sinThetaI * sin2kAlpha[1];
            Float asinoIp = _SafeASin(std::abs(sinThetaIp));
            floorthetaIp = Clamp(std::floor(asinoIp / unit), 0, thetanum-1);
        }

        // Handle remainder of $p$ values for hair scale tilt
        else if (p == 1) {
            sinThetaIp = sinThetaI * cos2kAlpha[0] - cosThetaI * sin2kAlpha[0];
            cosThetaIp = cosThetaI * cos2kAlpha[0] + sinThetaI * sin2kAlpha[0];

            *pdf += apPdf[p] * _Np(dphi, p, s, gammaO, gammaT);
        } else if (p == 2) {
            sinThetaIp = sinThetaI * cos2kAlpha[2] - cosThetaI * sin2kAlpha[2];
            cosThetaIp = cosThetaI * cos2kAlpha[2] + sinThetaI * sin2kAlpha[2];
        } else {
            sinThetaIp = sinThetaI;
            cosThetaIp = cosThetaI;
        }



        // Handle out-of-range $\cos \thetai$ from scale adjustment
        cosThetaIp = std::abs(cosThetaIp);
        *pdf += _Mp(cosThetaIp, cosThetaO, sinThetaIp, sinThetaO, v[p]) *apPdf[p] * _Np(dphi, p, s, gammaO, gammaT);
    }

    //*pdf = pp/60;

    *pdf += _Mp(cosThetaI, cosThetaO, sinThetaI, sinThetaO, v[_pMax]) *apPdf[_pMax] * (1 / (2 * Pi));

  
    *pdf = 1;//pdfについて,鏡面反射を想定してます。

    // if (std::abs(wi->x) < .9999) CHECK_NEAR(*pdf, Pdf(wo, *wi), .01);
    return f(wo, *wi);
}


Float MorphoBSDF::Pdf(const Vector3f &wo, const Vector3f &wi) const {
    // Compute hair coordinate system terms related to _wo_
    Float sinThetaO = wo.x;
    Float cosThetaO = _SafeSqrt(1 - _Sqr(sinThetaO));
    Float phiO = std::atan2(wo.z, wo.y);

    // Compute hair coordinate system terms related to _wi_
    Float sinThetaI = wi.x;
    Float cosThetaI = _SafeSqrt(1 - _Sqr(sinThetaI));
    Float phiI = std::atan2(wi.z, wi.y);

    // Compute $\gammat$ for refracted ray
    Float etap = std::sqrt(eta * eta - _Sqr(sinThetaO)) / cosThetaO;
    Float sinGammaT = h / etap;
    Float gammaT = _SafeASin(sinGammaT);

    // Compute PDF for $A_p$ terms
    std::array<Float, _pMax + 1> apPdf = ComputeApPdf(cosThetaO,wo);

    // Compute PDF sum for hair scattering events
    Float phi = phiI - phiO;
    Float pdf = 0;
    for (int p = 0; p < _pMax; ++p) {
        // Compute $\sin \thetai$ and $\cos \thetai$ terms accounting for scales
        Float sinThetaIp, cosThetaIp;
        if (p == 0) {
            sinThetaIp = sinThetaI * cos2kAlpha[1] + cosThetaI * sin2kAlpha[1];
            cosThetaIp = cosThetaI * cos2kAlpha[1] - sinThetaI * sin2kAlpha[1];
        }

        // Handle remainder of $p$ values for hair scale tilt
        else if (p == 1) {
            sinThetaIp = sinThetaI * cos2kAlpha[0] - cosThetaI * sin2kAlpha[0];
            cosThetaIp = cosThetaI * cos2kAlpha[0] + sinThetaI * sin2kAlpha[0];
        } else if (p == 2) {
            sinThetaIp = sinThetaI * cos2kAlpha[2] - cosThetaI * sin2kAlpha[2];
            cosThetaIp = cosThetaI * cos2kAlpha[2] + sinThetaI * sin2kAlpha[2];
        } else {
            sinThetaIp = sinThetaI;
            cosThetaIp = cosThetaI;
        }

        // Handle out-of-range $\cos \thetai$ from scale adjustment
        cosThetaIp = std::abs(cosThetaIp);
        pdf += _Mp(cosThetaIp, cosThetaO, sinThetaIp, sinThetaO, v[p]) *
               apPdf[p] * _Np(phi, p, s, gammaO, gammaT);
    }
    pdf += _Mp(cosThetaI, cosThetaO, sinThetaI, sinThetaO, v[_pMax]) *
           apPdf[_pMax] * (1 / (2 * Pi));
    return pdf;
}

std::string MorphoBSDF::ToString() const {
    return StringPrintf(
        "[ Hair h: %f gammaO: %f eta: %f beta_m: %f beta_n: %f "
        "v[0]: %f s: %f sigma_a: ", h, gammaO, eta, beta_m, beta_n,
        v[0], s) +
        sigma_a.ToString() +
        std::string("  ]");
}

RGBSpectrum MorphoBSDF::SigmaAFromConcentration(Float ce, Float cp) {
    Float sigma_a[3];
    Float eumelaninSigmaA[3] = {0.419f, 0.697f, 1.37f};
    Float pheomelaninSigmaA[3] = {0.187f, 0.4f, 1.05f};
    for (int i = 0; i < 3; ++i)
        sigma_a[i] = (ce * eumelaninSigmaA[i] + cp * pheomelaninSigmaA[i]);
    return RGBSpectrum::FromRGB(sigma_a);
}

RGBSpectrum MorphoBSDF::SigmaAFromReflectance(const RGBSpectrum &c, Float beta_n) {
    Spectrum sigma_a;
    for (int i = 0; i < Spectrum::nSamples; ++i)
        sigma_a[i] = _Sqr(std::log(c[i]) /
                         (5.969f - 0.215f * beta_n + 2.532f * _Sqr(beta_n) -
                          10.73f * _Pow<3>(beta_n) + 5.574f * _Pow<4>(beta_n) +
                          0.245f * _Pow<5>(beta_n)));
    return sigma_a;
}

}  // namespace pbrt