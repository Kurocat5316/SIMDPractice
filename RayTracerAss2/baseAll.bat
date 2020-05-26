set runs=%1
cd x64
Release\RayTracerAss2.exe -runs %runs%
Release\RayTracerAss2.exe -runs %runs% -samples 4
Release\RayTracerAss2.exe -runs %runs% -size 500 300
Release\RayTracerAss2.exe -runs %runs% -input ../Scenes/cornell-256lights.txt -size 512 512
Release\RayTracerAss2.exe -runs %runs% -input ../Scenes/allmaterials.txt
Release\RayTracerAss2.exe -runs %runs% -input ../Scenes/5000spheres.txt -size 960 540
Release\RayTracerAss2.exe -runs %runs% -input ../Scenes/bunny500.txt
Release\RayTracerAss2.exe -runs %runs% -input ../Scenes/bunny10k.txt -size 256 256
cd ..
