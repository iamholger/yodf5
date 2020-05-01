#include "YODA/ReaderYODA.h"
#include "YODA/Histo1D.h"
#include <iostream>
#include <regex>
#include <unordered_map>
#include <functional>


#include "boost/multi_array.hpp"
#define H5_USE_BOOST
#include <highfive/H5Easy.hpp>
#include <highfive/H5DataSet.hpp>
#include <highfive/H5DataSpace.hpp>
#include <highfive/H5File.hpp>

using namespace YODA;
using namespace std;
using namespace HighFive;


// pre C++20
inline bool ends_with(string const & value, string const & ending)
{
    if (ending.size() > value.size()) return false;
    return equal(ending.rbegin(), ending.rend(), value.rbegin());
}

inline vector<string> get_ao_names(const vector<AnalysisObject*>& aos, string const & aotype) {
    vector<string> aonames;
    for (auto ao : aos) 
    {
        if (ao->type()==aotype &! ao->path().ends_with( ']')) aonames.push_back(ao->path());
    }
    sort(aonames.begin(), aonames.end());
    return aonames;
}


inline vector<string> get_variations(const vector<AnalysisObject*>& aos, const string & hname) {
    vector<string> vars= {""};
    for (auto ao : aos) 
    {
        if (ao->path().starts_with(hname+"[") && ao->path().ends_with( ']'))
        {
            string s(ao->path());
            s.erase(s.begin(), s.begin()+hname.length()+1);
            s.erase(s.end()-1);
            vars.push_back(s);
        }
    }
    sort(vars.begin(), vars.end());
    return vars;
}

unordered_map<string, int> mk_nbinmap(unordered_map<string, AnalysisObject* > aomap, vector<string> const & hnames) {
    unordered_map<string, int> nbin_map;
    for (auto hn : hnames) {
        nbin_map.insert({hn, dynamic_cast<Histo1D*>(aomap[hn])->numBins()});
    }
    return nbin_map;
}

inline vector<string> mk_binids(vector<AnalysisObject*> const & aos, vector<string> const & hnames) 
{
    vector<string> binids;
    for (auto hn : hnames) {
        for (auto ao : aos) 
        {
            if (ao->path() == hn) {
                binids.push_back(hn+"#T");
                binids.push_back(hn+"#O");
                binids.push_back(hn+"#U");
                for (int nb=0;nb<dynamic_cast<Histo1D*>(ao)->numBins();++nb) binids.push_back(hn+"#" + to_string(nb));
                break;
            }
        }
    }
    return binids;
}

inline vector<string> mk_binids(unordered_map<string, AnalysisObject* > aomap, vector<string> const & hnames) 
{
    vector<string> binids;
    for (auto hn : hnames) 
    {
        binids.push_back(hn+"#T");
        binids.push_back(hn+"#O");
        binids.push_back(hn+"#U");
        for (int nb=0;nb<dynamic_cast<Histo1D*>(aomap[hn])->numBins();++nb) binids.push_back(hn+"#" + to_string(nb));
    }
    return binids;
}


inline vector<double> get_edge(unordered_map<string, AnalysisObject* > aomap, vector<string> const & hnames, double (HistoBin1D::*function)() const) 
{
    vector<double> edges;
    for (auto hn : hnames) 
    {
        edges.push_back(0.0);
        edges.push_back(0.0);
        edges.push_back(0.0);
        auto h = dynamic_cast<Histo1D*>(aomap[hn]);
        for (int nb=0;nb< h->numBins();++nb) edges.push_back((h->bin(nb).*function)());
    }
    return edges;
}

// TODO add variant that uses Eigen::ArrayXXd
boost::multi_array<double,2>  get_h1d_field(
        unordered_map<string, AnalysisObject* > aomap, 
        vector<string> const & hnames, vector<string> const & vars, size_t nbins,
        unordered_map<string, int>  nbinsMap,
        double (Dbn1D::*function)() const,
        double (HistoBin1D::*function2)() const
        )
{
    boost::multi_array<double, 2> field(boost::extents[nbins][vars.size()]);
    size_t offset(0);
    for (auto  hn : hnames) 
    {
        for (size_t j=0; j<vars.size();++j)
        {
            std::string varhname = (j>0) ? hn+"["+vars[j]+"]" : hn;
            if (aomap.count(varhname)>0)
            {
                auto h = dynamic_cast<Histo1D*>(aomap[varhname]);
                field[offset  ][j] = (h->totalDbn().*function)();
                field[offset+1][j] = (h->overflow().*function)();
                field[offset+2][j] = (h->underflow().*function)();
                for (int nb=0;nb< h->numBins();++nb) field[offset + 3 + nb][j] = (h->bin(nb).*function2)();
                if (j==vars.size()-1) offset += h->numBins() + 3;
            }
            else 
            {
                field[offset  ][j] = 0.0;
                field[offset+1][j] = 0.0;
                field[offset+2][j] = 0.0;
                for (int nb=0;nb< nbinsMap[hn];++nb) field[offset + 3 + nb][j] = 0.0;
                if (j==vars.size()-1) offset += nbinsMap[hn] + 3;
            }
        }
    }
    return field;
}

int main(int argc, const char** argv)
{
    int compression=atoi(argv[1]);
    std::cerr << argv[1] << "\n";
    size_t nFiles=argc-2;
    std::cerr << argv[2] << " "<< nFiles << "\n";
    auto aos = ReaderYODA::create().read(argv[2]);
    unordered_map<string, AnalysisObject*> aomap;
    for (auto ao : aos) aomap.insert({ao->path(), ao});

    // Determines sorting order --- in a robust application all hnames must be known
    auto const h1d_names = get_ao_names(aos, "Histo1D");
    auto const binids    = mk_binids(aomap, h1d_names);
    auto const vars      = get_variations(aos, h1d_names[0]);
    auto nbinmap = mk_nbinmap(aomap, h1d_names);

    H5Easy::File file("example.h5", H5Easy::File::ReadWrite | H5Easy::File::Create| H5Easy::File::Truncate);
    
    file.createGroup("/Histo1D");

    DataSetCreateProps props;
    props.add(Deflate(compression));
    // This is the same chunking strategy h5py applies seemingly
    props.add(Chunking(std::vector<hsize_t>{
                size_t(ceil(binids.size()/8.)),
                size_t(ceil(vars.size()/8.)),
                size_t(ceil(nFiles/8.))
                }));
    
    //DataSetAccessProps cacheConfig;
    //cacheConfig.add(Caching(1024, 1024, 0.5));

    // Initial DS extends
    size_t const NB = binids.size();
    size_t const NV = vars.size();
    size_t const NF = nFiles;


    // Dataset exttensible in 3rd dim, i.e. nfiles
    file.createDataSet<double>("/Histo1D/sumW",       DataSpace( { NB, NV, NF},{ NB, NV, DataSpace::UNLIMITED}), props);//, cacheConfig );
    file.createDataSet<double>("/Histo1D/sumW2",      DataSpace( { NB, NV, NF},{ NB, NV, DataSpace::UNLIMITED}), props);//, cacheConfig );
    file.createDataSet<double>("/Histo1D/sumWX",      DataSpace( { NB, NV, NF},{ NB, NV, DataSpace::UNLIMITED}), props);//, cacheConfig );
    file.createDataSet<double>("/Histo1D/sumWX2",     DataSpace( { NB, NV, NF},{ NB, NV, DataSpace::UNLIMITED}), props);//, cacheConfig );
    file.createDataSet<double>("/Histo1D/numEntries", DataSpace( { NB, NV, NF},{ NB, NV, DataSpace::UNLIMITED}), props);//, cacheConfig );

    // Dummy here, we write the data of the first input file as many times as there are command line arguments
    for (size_t fidx=0; fidx<NF; ++fidx)
    {
        file.getDataSet("/Histo1D/sumW").select(       {0, 0, fidx}, {NB, NV, 1}).write(get_h1d_field(aomap, h1d_names, vars, NB, nbinmap, &Dbn1D::sumW,       &HistoBin1D::sumW));
        file.getDataSet("/Histo1D/sumW2").select(      {0, 0, fidx}, {NB, NV, 1}).write(get_h1d_field(aomap, h1d_names, vars, NB, nbinmap, &Dbn1D::sumW2,      &HistoBin1D::sumW2));
        file.getDataSet("/Histo1D/sumWX").select(      {0, 0, fidx}, {NB, NV, 1}).write(get_h1d_field(aomap, h1d_names, vars, NB, nbinmap, &Dbn1D::sumWX,      &HistoBin1D::sumWX));
        file.getDataSet("/Histo1D/sumWX2").select(     {0, 0, fidx}, {NB, NV, 1}).write(get_h1d_field(aomap, h1d_names, vars, NB, nbinmap, &Dbn1D::sumWX2,     &HistoBin1D::sumWX2));
        file.getDataSet("/Histo1D/numEntries").select( {0, 0, fidx}, {NB, NV, 1}).write(get_h1d_field(aomap, h1d_names, vars, NB, nbinmap, &Dbn1D::numEntries, &HistoBin1D::numEntries));
    }
    
    file.createDataSet("Histo1D/xMin",   get_edge(aomap, h1d_names, &HistoBin1D::xMin));
    file.createDataSet("Histo1D/xMax",   get_edge(aomap, h1d_names, &HistoBin1D::xMax));
    
    H5Easy::dump(file, "Histo1D/binids", binids);
    H5Easy::dump(file, "Histo1D/names", h1d_names);
    

    return 0;
}
