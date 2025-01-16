cmake -G "MinGW Makefiles" -DCMAKE_CXX_STANDARD=17 -DCMAKE_CXX_FLAGS="-fexec-charset=UTF-8" ..
make


cd build
cmake -G "Visual Studio 17 2022" ..
cmake --build .

./build/Release/lajolla.exe ./scenes/cbox/cbox.xml

msbuild ./build/lajolla.sln /p:Configuration=Release
# HW1
./build/Release/lajolla.exe .\scenes\disney_bsdf_test\simple_sphere.xml
./build/Release/lajolla.exe .\scenes\disney_bsdf_test\disney_diffuse.xml
./build/Release/lajolla.exe .\scenes\disney_bsdf_test\disney_metal.xml
./build/Release/lajolla.exe .\scenes\disney_bsdf_test\disney_clearcoat.xml
./build/Release/lajolla.exe .\scenes\disney_bsdf_test\disney_glass.xml
./build/Release/lajolla.exe .\scenes\disney_bsdf_test\disney_sheen.xml
./build/Release/lajolla.exe .\scenes\disney_bsdf_test\disney_bsdf.xml
./build/Release/lajolla.exe .\scenes\disney_bsdf_test\disney_bsdf_array.xml

# HW2
./build/Release/lajolla.exe .\scenes\volpath_test\volpath_test1.xml
./build/Release/lajolla.exe .\scenes\volpath_test\volpath_test2.xml
