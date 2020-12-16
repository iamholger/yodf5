#include <iostream>
#include <regex>
#include <unordered_map>
#include <functional>

#include <xtensor/xtensor.hpp>
#include <xtensor/xrandom.hpp>

#define H5_USE_XTENSOR
#include <highfive/H5Easy.hpp>
#include <highfive/H5DataSet.hpp>
#include <highfive/H5DataSpace.hpp>
#include <highfive/H5File.hpp>

using namespace std;
using namespace HighFive;


int main(int argc, const char** argv)
{
    int compression=atoi(argv[1]);

    H5Easy::File file("example.h5", H5Easy::File::ReadWrite | H5Easy::File::Create| H5Easy::File::Truncate);
    

    const size_t nbins = 10000;
    const size_t nvars = 250;
    const size_t nobjs = 10;
  
    xt::xtensor<double,3> E = xt::zeros<double>({nbins,nvars,nobjs});
    //xt::xtensor<double,3> E = xt::random::rand<double>({nbins,nvars,nobjs});


    H5Easy::dump(file, "E", E, H5Easy::DumpOptions(H5Easy::Compression(compression)));

    return 0;
}
