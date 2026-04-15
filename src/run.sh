rm -rf build
mkdir build
cd build
cmake ..
cmake --build . -j
./stump_demo
cd ../