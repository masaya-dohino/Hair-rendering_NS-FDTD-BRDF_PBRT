# Hair-rendering nano-structure optical BRDF
Hair rendering with BRDF by NS-FDTD on PBRT

## Over view
Windows

・Install git: https://git-scm.com/

・Install cmake: https://cmake.org/download/

・Install pbrt-v4： https://github.com/mmp/pbrt-v4

## PBRT 
https://pbr-book.org/

(reference)

Photorealistic computer graphics is ubiquitous today, with applications that include entertainment—notably, movies and video games; product design; and architecture. Over the past decade, physically-based approaches to rendering have become widely used, where an accurate modeling of the physics of light scattering is at the heart of image synthesis. These approaches offer both visual realism and predictability.


Physically Based Rendering describes both the mathematical theory behind a modern photorealistic rendering system as well as its practical implementation. A method known as “literate programming” combines human-readable documentation and source code into a single reference that is specifically designed to aid comprehension. The ideas and software in this book show the reader how to design and employ a full-featured rendering system capable of creating stunning imagery.

## How to use
### set up
・command prompt (administer), download.

```
cd %HOMEPATH%
rmdir /s /q pbrt-v4
git clone --recursive https://github.com/mmp/pbrt-v4
```

・Execute cmake

```
cd %HOMEPATH%
cd pbrt-v4
rmdir /s /q build
mkdir build
cd build
cmake .. -G "Visual Studio 17 2022" -A x64 -T host=x64
```

・Build

```
cd %HOMEPATH%
cd pbrt-v4
rmdir /s /q build
mkdir build
cd build
cmake .. -G "Visual Studio 17 2022" -A x64 -T host=x64
```

・Set PATH(c:\Program Files\PBRT-V4\bin). command prompt (administer)

```
call powershell -command "$oldpath = [System.Environment]::GetEnvironmentVariable(\"Path\", \"Machine\"); $oldpath += \";c:\Program Files\PBRT-V4\bin\"; [System.Environment]::SetEnvironmentVariable(\"Path\", $oldpath, \"Machine\")
```
・Sample rendering of SCENE file(.pbrt).

```
cd %HOMEPATH%
cd pbrt-v4
"c:\Program Files\PBRT-V4\bin\pbrt.exe" --outfile a.png ..\pbrt-v4-scenes\killeroos\killeroo-simple.pbrt
```


## Result
### Human hair
・Structure of human hair

This section describes the structure of human hair, which is important for this research. The human head contains approximately 100,000 hairs with diameters ranging from 50 to 100 µm, so the appearance of the hair reflects its optical properties as a fiber aggregate. Human hair basically consists of three layers. The cuticle, located on the top surface of the hair, is a transparent, plate-like cell about 0.5 µm thick made of non-keratin proteins, with a layered structure of about 5 or 6 overlapping cells. This layered structure is directional. Between the layers of cuticle cells, there is a cell membrane complex (CMC) about 30 nm thick, which serves as a pathway for substances to penetrate into the hair. And inside this is the cortex, which contains melanin granules that determine the hair's color, and the center is composed of the medulla, which is about 5% of the total.

One of the structures of damaged hair is the hollowing out of the cortex and medulla tissues, which leads to the loss of elasticity and luster. In addition, there is cuticle lifting on the surface. Cuticle lift-up is a condition in which the scaly cuticle layer shown in the figure on the right is missing or peeling off, and when the cuticle thins due to these factors, it becomes a cause of split ends and hair breakage. In this paper, we will focus on the structural changes of the cuticle as the structure of the hair surface.

![image](https://user-images.githubusercontent.com/57475794/219864271-3128cc04-cd74-4d81-9bbb-c928cae34529.png)

### 2D BRDF
![image](https://user-images.githubusercontent.com/57475794/219870855-3fe7801e-6888-43f8-b7ab-4fc6cbc13c00.png)


### Rendaring with nano-structure optical BRDF
Original BRDF express nano-structure(hair) optical scattering. Original BRDF access the table of wavelength spectrum for incident and reflected angles from calculation of NS-FDTD argolhythm.

![image](https://user-images.githubusercontent.com/57475794/219864552-4ce93d29-ccbc-4ef3-8f1e-2c673f64edfa.png)



## Learning
The following figure shows the BRDF results for each hair surface structure rendered with cylindrical geometry as fiber aggregates, using BRDF implemented from the values obtained by numerical simulation. A collimated light source (5800K) is used as a light source to simulate sunlight. The results for the multi-layered cuticle structure and the layer-less structure show weaker white highlights, giving the impression of a wet sheen to the hair as a whole. The exfoliated structure, the missing structure, and the exfoliated plus missing structure, which we defined as hair damage, have generally stronger white highlights and yellow reflections, and seem to have less effect in terms of gloss and transparency than the multi-layered cuticle structure. In particular, the increase of highlights in the peeled structure is remarkable and seems to be similar to the disordered state of cuticles in real life.

From the rendering results, it was found that through the simulation of the structure that imitates damaged hair surface, the white highlights became stronger, just like real hair, and that changes in the hair's microstructure affect the CG representation of the hair. Through these results, it can be inferred that the generation of cavities (gaps) due to cuticle shedding and lifting is one of the factors affecting the appearance of damaged hair.

## Next
‣Need to increase the number of samples for BRDF calculations and to increase the computational domain to reduce the influence of boundary conditions, thus increasing the speed of numerical calculations.

‣Extend NS-FDTD simulations to 3-D.

‣Need a rendering method that takes into account the entire cuticle, cortices, medulla, and melanin pigments.

‣Research is needed to determine why the width of the highlight effect (angel rings) differs, e.g. by increasing the number of hairs and performing ray tracing.
