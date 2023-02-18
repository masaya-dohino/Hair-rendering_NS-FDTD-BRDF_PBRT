# Hair-rendering_NS-FDTD-BRDF_PBRT
Hair rendering with BRDF by NS-FDTD on PBRT

## Over view
Windows

Install git:

Install cmake: 

Install pbrt-v4：

## PBRT 

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

## Result


## Learning
