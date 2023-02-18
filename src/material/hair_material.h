#define NOMINMAX
#pragma once

// materials/mirror.h*
#include "pbrt.h"
#include "material.h"
#include "reflection.h"
#include "bssrdf.h"
#include <array>

namespace pbrt {


//hair

class MorphoMaterial : public Material {
  public:
    // HairMaterial Public Methods
    MorphoMaterial(const std::shared_ptr<Texture<RGBSpectrum>> &sigma_a,
                 const std::shared_ptr<Texture<RGBSpectrum>> &color,
                 const std::shared_ptr<Texture<Float>> &eumelanin,
                 const std::shared_ptr<Texture<Float>> &pheomelanin,
                 const std::shared_ptr<Texture<Float>> &eta,
                 const std::shared_ptr<Texture<Float>> &beta_m,
                 const std::shared_ptr<Texture<Float>> &beta_n,
                 const std::shared_ptr<Texture<Float>> &alpha)
        : sigma_a(sigma_a),
          color(color),
          eumelanin(eumelanin),
          pheomelanin(pheomelanin),
          eta(eta),
          beta_m(beta_m),
          beta_n(beta_n),
          alpha(alpha){}
    void ComputeScatteringFunctions(SurfaceInteraction *si, MemoryArena &arena,
                                    TransportMode mode,
                                    bool allowMultipleLobes) const;

  private:
    // HairMaterial Private Data
    std::shared_ptr<Texture<RGBSpectrum>> sigma_a, color;
    std::shared_ptr<Texture<Float>> eumelanin, pheomelanin, eta;
    std::shared_ptr<Texture<Float>> beta_m, beta_n, alpha;
};

MorphoMaterial *CreateMorphoMaterial(const TextureParams &mp);

// HairBSDF Constants

//static const int _pMax = 1;
static const int _pMax = 3;
static const Float _SqrtPiOver8 = 0.626657069f;


// HairBSDF Declarations
class MorphoBSDF : public BxDF {
  public:
    // HairBSDF Public Methods
    MorphoBSDF(Float h, Float eta, const RGBSpectrum &sigma_a, Float beta_m,
             Float beta_n, Float alpha,int wavelengthindex);
    Spectrum f(const Vector3f &wo, const Vector3f &wi) const;
    Spectrum Sample_f(const Vector3f &wo, Vector3f *wi, const Point2f &u,
                      Float *pdf, BxDFType *sampledType)const;
    Float Pdf(const Vector3f &wo, const Vector3f &wi) const;
    Float Pdf1(const Vector3f &wo, const Vector3f &wi) const;
    std::string ToString() const;
    static RGBSpectrum SigmaAFromConcentration(Float ce, Float cp);
    static RGBSpectrum SigmaAFromReflectance(const RGBSpectrum &c, Float beta_n);

    

  private:
    // HairBSDF Private Methods
    std::array<Float, _pMax + 1> ComputeApPdf(Float cosThetaO, const Vector3f &wo) const;

    // HairBSDF Private Data
    const Float h, gammaO, eta;
    const RGBSpectrum sigma_a;
    const Float beta_m, beta_n;
    Float v[_pMax + 1];
    Float s;
    Float sin2kAlpha[3], cos2kAlpha[3];
    int wavelengthindex;
};



// General Utility Functions
inline Float _Sqr(Float v) { return v * v; }
template <int n>
static Float _Pow(Float v) {
    static_assert(n > 0, "Power can't be negative");
    Float n2 = _Pow<n / 2>(v);
    return n2 * n2 * _Pow<n & 1>(v);
}

template <>
inline Float _Pow<1>(Float v) {
    return v;
}
template <>
inline Float _Pow<0>(Float v) {
    return 1;
}

inline Float _SafeASin(Float x) {//arcsin
    CHECK(x >= -1.0001 && x <= 1.0001);
    return std::asin(Clamp(x, -1, 1));
}



inline Float _SafeSqrt(Float x) {
    CHECK_GE(x, -1e-4);
    return std::sqrt(std::max(Float(0), x));
}

// https://fgiesen.wordpress.com/2009/12/13/decoding-morton-codes/
static uint32_t _Compact1By1(uint32_t x) {
    // TODO: as of Haswell, the PEXT instruction could do all this in a
    // single instruction.
    // x = -f-e -d-c -b-a -9-8 -7-6 -5-4 -3-2 -1-0
    x &= 0x55555555;
    // x = --fe --dc --ba --98 --76 --54 --32 --10
    x = (x ^ (x >> 1)) & 0x33333333;
    // x = ---- fedc ---- ba98 ---- 7654 ---- 3210
    x = (x ^ (x >> 2)) & 0x0f0f0f0f;
    // x = ---- ---- fedc ba98 ---- ---- 7654 3210
    x = (x ^ (x >> 4)) & 0x00ff00ff;
    // x = ---- ---- ---- ---- fedc ba98 7654 3210
    x = (x ^ (x >> 8)) & 0x0000ffff;
    return x;
}

static Point2f _DemuxFloat(Float f) {
    CHECK(f >= 0 && f < 1);
    uint64_t v = f * (1ull << 32);
    CHECK_LT(v, 0x100000000);
    uint32_t bits[2] = {_Compact1By1(v), _Compact1By1(v >> 1)};
    return {bits[0] / Float(1 << 16), bits[1] / Float(1 << 16)};
}


}  // namespace pbrt