#include "YODA/ReaderYODA.h"
#include "YODA/Histo1D.h"
#include <iostream>
#include <regex>
#include <unordered_map>
#include <functional>

#include "Eigen/Dense"

#include "boost/multi_array.hpp"
#define H5_USE_BOOST
#define H5_USE_EIGEN
#include <highfive/H5Easy.hpp>
#include <highfive/H5DataSet.hpp>
#include <highfive/H5DataSpace.hpp>
#include <highfive/H5File.hpp>

using namespace YODA;
using namespace std;
using namespace Eigen;
using namespace HighFive;


unordered_map< string, unordered_map< string, int> > myyodaindex(istream& stream_) {

    unordered_map<string, unordered_map< string, int>> hmap;
    hmap.insert({"Histo1D",   unordered_map< string, int>()});
    hmap.insert({"Histo2D",   unordered_map< string, int>()});
    hmap.insert({"Profile1D", unordered_map< string, int>()});
    hmap.insert({"Profile2D", unordered_map< string, int>()});
    hmap.insert({"Scatter1D", unordered_map< string, int>()});
    hmap.insert({"Scatter2D", unordered_map< string, int>()});
    hmap.insert({"Scatter3D", unordered_map< string, int>()});
    hmap.insert({"Counter",   unordered_map< string, int>()});

    //#ifdef HAVE_LIBZ
    //// NB. zstr auto-detects if file is deflated or plain-text
    //zstr::istream stream(stream_);
    //#else
    istream& stream = stream_;
    //#endif
    enum Context
    {
        NONE, SCATTER1D, SCATTER2D, SCATTER3D, COUNTER, HISTO1D, HISTO2D,
        PROFILE1D, PROFILE2D
    };


    /// State of the parser: line number, line, parser context, and pointer(s) to the object currently being assembled
    unsigned int nline = 0;
    string s;
    Context context = NONE;
    std::string curpath="";
    int nbins=0;

    // Loop over all lines of the input file
    bool in_anns = false;
    string fmt = "1";
    while (std::getline(stream, s))
    {
        nline += 1;
        if (!in_anns)
        {
            Utils::itrim(s);
            if (s.empty()) continue;
            if (s.find("#") == 0 && s.find("BEGIN") == string::npos && s.find("END") == string::npos) continue;
        }
        if (context == NONE)
        {
            if (s.find("BEGIN ") == string::npos)
            {
                stringstream ss;
                ss << "Unexpected line in YODA format parsing when BEGIN expected: '" << s << "' on line " << nline;
                throw ReadError(ss.str());
            }
            while (s.find("#") == 0) s = Utils::trim(s.substr(1));
            vector<string> parts;
            istringstream iss(s); string tmp;
            while (iss >> tmp) parts.push_back(tmp);

            if (parts.size() < 2 || parts[0] != "BEGIN")
            {
                stringstream ss;
                ss << "Unexpected BEGIN line structure when BEGIN expected: '" << s << "' on line " << nline;
                throw ReadError(ss.str());
            }
            const string ctxstr = parts[1];
            curpath = (parts.size() >= 3) ? parts[2] : "";
            nbins=0;
            if      (Utils::startswith(ctxstr, "YODA_COUNTER"))   {context = COUNTER;  }
            else if (Utils::startswith(ctxstr, "YODA_SCATTER1D")) {context = SCATTER1D;}
            else if (Utils::startswith(ctxstr, "YODA_SCATTER2D")) {context = SCATTER2D;}
            else if (Utils::startswith(ctxstr, "YODA_SCATTER3D")) {context = SCATTER3D;}
            else if (Utils::startswith(ctxstr, "YODA_HISTO1D"))   {context = HISTO1D;  }
            else if (Utils::startswith(ctxstr, "YODA_HISTO2D"))   {context = HISTO2D;  }
            else if (Utils::startswith(ctxstr, "YODA_PROFILE1D")) {context = PROFILE1D;}
            else if (Utils::startswith(ctxstr, "YODA_PROFILE2D")) {context = PROFILE2D;}
            const size_t vpos = ctxstr.find_last_of("V");
            fmt = vpos != string::npos ? ctxstr.substr(vpos+1) : "1";
            if (fmt != "1") in_anns = true;
        }
        else
        { //< not a BEGIN line
            if (s.find("BEGIN ") != string::npos) throw ReadError("Unexpected BEGIN line in YODA format parsing before ending current BEGIN..END block");
            // FINISHING THE CURRENT CONTEXT
            if (s.find("END ") != string::npos)
            {
                if      (context == HISTO1D)   {hmap["Histo1D"]  .insert({curpath, nbins});}
                else if (context == HISTO2D)   {hmap["Histo2D"]  .insert({curpath, nbins});}
                else if (context == PROFILE1D) {hmap["Profile1D"].insert({curpath, nbins});}
                else if (context == PROFILE2D) {hmap["Profile2D"].insert({curpath, nbins});}
                else if (context == SCATTER1D) {hmap["Scatter1D"].insert({curpath, nbins});}
                else if (context == SCATTER2D) {hmap["Scatter2D"].insert({curpath, nbins});}
                else if (context == SCATTER3D) {hmap["Scatter3D"].insert({curpath, nbins});}
                else if (context == COUNTER)   {hmap["Counter"]  .insert({curpath, nbins});}
                in_anns = false;
                context = NONE;
                continue;
            }
            // ANNOTATIONS PARSING
            if (fmt == "1")
            {
                // First convert to one-key-per-line YAML syntax
                const size_t ieq = s.find("=");
                if (ieq != string::npos) s.replace(ieq, 1, ": ");
                // Special-case treatment for syntax clashes
                const size_t icost = s.find(": *");
                if (icost != string::npos)
                {
                    s.replace(icost, 1, ": '*");
                    s += "'";
                }
                // Store reformatted annotation
                const size_t ico = s.find(":");
                if (ico != string::npos) continue;
            }
            else if (in_anns)
            {
                if (s == "---") in_anns = false;
                continue;
            }

            if ( (context == HISTO1D) || (context==HISTO2D) || (context==PROFILE1D) || (context==PROFILE2D))
            {
                if (s.find("Total") != string::npos || s.find("Underflow") != string::npos || s.find("Overflow") != string::npos) continue;
                nbins++;
            }
            else if ( (context==SCATTER1D)  || (context==SCATTER2D) || (context==SCATTER3D))
            {
                nbins++;
            }
        }
    }
    return hmap;
}



// pre C++20
inline bool ends_with(string const & value, string const & ending)
{
    if (ending.size() > value.size()) return false;
    return equal(ending.rbegin(), ending.rend(), value.rbegin());
}

inline vector<string> get_ao_names(const vector<AnalysisObject*>& aos, string const & aotype)
{
    vector<string> aonames;
    for (auto ao : aos) 
    {
        if (ao->type()==aotype &! ao->path().ends_with( ']')) aonames.push_back(ao->path());
    }
    sort(aonames.begin(), aonames.end());
    return aonames;
}


inline vector<string> get_variations(const vector<AnalysisObject*>& aos, const string & hname)
{
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

unordered_map<string, int> mk_nbinmap(unordered_map<string, AnalysisObject* > aomap, vector<string> const & hnames)
{
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


MatrixXd get_h1d_field_eig(
        unordered_map<string, AnalysisObject* > aomap,
        vector<string> const & hnames, vector<string> const & vars, size_t nbins,
        unordered_map<string, int>  nbinsMap,
        double (Dbn1D::*function)() const,
        double (HistoBin1D::*function2)() const
        )
{
    MatrixXd field(vars.size(), nbins);
    size_t offset(0);
    for (auto  hn : hnames)
    {
        for (size_t j=0; j<vars.size();++j)
        {
            std::string varhname = (j>0) ? hn+"["+vars[j]+"]" : hn;
            if (aomap.count(varhname)>0)
            {
                auto h = dynamic_cast<Histo1D*>(aomap[varhname]);
                field(j, offset  ) = (h->totalDbn().*function)();
                field(j, offset+1) = (h->overflow().*function)();
                field(j, offset+2) = (h->underflow().*function)();
                for (int nb=0;nb< h->numBins();++nb) field(j, offset + 3 + nb) = (h->bin(nb).*function2)();
                if (j==vars.size()-1) offset += h->numBins() + 3;
            }
            else
            {
                field(j, offset  ) = 0.0;
                field(j, offset+1) = 0.0;
                field(j, offset+2) = 0.0;
                for (int nb=0;nb< nbinsMap[hn];++nb) field(j, offset + 3 + nb) = 0.0;
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


    std::ifstream instream;
    instream.open(argv[2]);
    // This is a map of AOtype : {hname: nbins} --- runs 10x faster than read -- good for initial survey when
    // running over many files
    auto hmap = myyodaindex(instream);
    instream.close();
    auto nbinmap = hmap["Histo1D"];


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
        file.getDataSet("/Histo1D/sumW").select(       {0, 0, fidx}, {NB, NV, 1}).write(get_h1d_field_eig(aomap, h1d_names, vars, NB, nbinmap, &Dbn1D::sumW,       &HistoBin1D::sumW));
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
