#clang++ -O3 -march=native readH5YODA.cxx -I$PWD/../eigen-3.3.9 -I$PWD/../HighFive-2.2.2/include -I$PWD/../YODA-1.8.5/local_llvm/include \
    #-L$PWD/../YODA-1.8.5/local_llvm/lib -lYODA -lhdf5 \
    #-o readH5YODA

#clang++ -g -std=c++20 -march=native readAndWrite.cc -I$PWD/../eigen-3.3.9 -I$PWD/../HighFive-2.2.2/include -I$PWD/../YODA-1.8.5/local_llvm/include \
    #-L$PWD/../YODA-1.8.5/local_llvm/lib -lYODA -lhdf5 \
    #-o readAndWrite
clang++ -g  -march=native testxtensor.cxx -I$PWD/../xtl-0.6.23/include -I$PWD/../xtensor-0.21.10/include -I$PWD/../HighFive-2.2.2/include -lhdf5 \
    -o testxtensor
