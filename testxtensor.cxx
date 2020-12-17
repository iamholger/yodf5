#include <iostream>
#include <regex>
#include <unordered_map>
#include <functional>

#include <xtensor/xtensor.hpp>

#define H5_USE_XTENSOR
#include <highfive/H5Easy.hpp>
#include <highfive/H5DataSet.hpp>
#include <highfive/H5DataSpace.hpp>
#include <highfive/H5File.hpp>

#include "YODA/Histo1D.h"
#include "YODA/ReaderYODA.h"

using namespace std;
using namespace YODA;
using namespace HighFive;

// Pre C++20 this is needed
inline bool ends_with(string const & value, string const & ending)
{
    if (ending.size() > value.size()) return false;
    return equal(ending.rbegin(), ending.rend(), value.rbegin());
}

inline bool starts_with(string const & value, string const & pattern)
{
  if (value.rfind(pattern, 0) == 0)
  {
    return true;
  }
  return false;
}



inline vector<string> get_ao_names(const vector<AnalysisObject*>& aos, string const & aotype)
{
    vector<string> aonames;
    for (auto ao : aos) 
    {
        if (ao->type()==aotype &! ::ends_with(ao->path(), string(1, ']'))) aonames.push_back(ao->path());
    }
    sort(aonames.begin(), aonames.end());
    return aonames;
}

inline vector<string> get_variations(const vector<AnalysisObject*>& aos, const string & hname)
{
    vector<string> vars= {""};
    for (auto ao : aos) 
    {
        if ( ::starts_with(ao->path(), hname+"[") && ::ends_with(ao->path(), string(1, ']')))
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


// We need similar functions for the other AO types
// TODO can this be simplified?
void fill_h1d_field_xt(
        unordered_map<string, AnalysisObject* > aomap,
        vector<string> const & hnames, vector<string> const & vars, size_t nbins,
        unordered_map<string, int>  nbinsMap,
        double (Dbn1D::*function)() const,
        double (HistoBin1D::*function2)() const, xt::xtensor<double, 3> & target, size_t zindex
        )
{
    size_t offset(0);
    for (auto  hn : hnames)
    {
        for (size_t j=0; j<vars.size();++j)
        {
            std::string varhname = (j>0) ? hn+"["+vars[j]+"]" : hn;
            if (aomap.count(varhname)>0)
            {
                auto h = dynamic_cast<Histo1D*>(aomap[varhname]);
                target(j, offset  , zindex) = (h->totalDbn().*function)();
                target(j, offset+1, zindex) = (h->overflow().*function)();
                target(j, offset+2, zindex) = (h->underflow().*function)();
                for (int nb=0;nb< h->numBins();++nb) target(j, offset + 3 + nb, zindex) = (h->bin(nb).*function2)();
                if (j==vars.size()-1) offset += h->numBins() + 3;
            }
            else
            {
                target(j, offset  , zindex) = 0.0;
                target(j, offset+1, zindex) = 0.0;
                target(j, offset+2, zindex) = 0.0;
                for (int nb=0;nb< nbinsMap[hn];++nb) target(j, offset + 3 + nb, zindex) = 0.0;
                if (j==vars.size()-1) offset += nbinsMap[hn] + 3;
            }
        }
    }
}
        

void fillScaledBy(unordered_map<string, AnalysisObject* > aomap,  vector<string> const & hnames, vector<string> const & vars, xt::xtensor<double, 2> & target)
{
  for (size_t i=0; i<hnames.size();++i)
  {
    for (size_t j=0; j<vars.size();++j)
    {
      std::string varhname = (j>0) ? hnames[i]+"["+vars[j]+"]" : hnames[i];
      if (aomap[varhname]->hasAnnotation("ScaledBy"))
      {
        target(i,j) = std::stod(aomap[varhname]->annotation("ScaledBy"));
      }
      else
      {
        target(i,j) = -1.;
      }
    }
  }
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


int main(int argc, const char** argv)
{
    int compression=atoi(argv[1]);

    H5Easy::File file("example.h5", H5Easy::File::ReadWrite | H5Easy::File::Create| H5Easy::File::Truncate);
    
  
    // Read all objects into a map
    auto aos = ReaderYODA::create().read(argv[2]);
    unordered_map<string, AnalysisObject*> aomap;
    for (auto ao : aos) aomap.insert({ao->path(), ao});

    // yuck
    std::ifstream instream;
    instream.open(argv[2]);
    auto amap = ReaderYODA::create().mkIndex(instream);

    // TODO this should be a loop over all AnalysisObject types in the current yoda file --- or later what
    // is know to the rivet analysis handler in terms of AOs
    // We will have one big h5 dataset for each AO type

      // Determines sorting order --- in a robust application all hnames must be known
      auto const h1d_names = get_ao_names(aos, "Histo1D");
      auto const binids    = ::mk_binids(aomap, h1d_names);
      auto const vars      = get_variations(aos, h1d_names[0]);
    

      // TODO check what this does
      auto nbinmap = amap.getAOIndex()["Histo1D"];

      // Dataset extents
      size_t NB = binids.size();
      size_t NV = vars.size();


      xt::xtensor<double, 3>::shape_type sh1 = {NB, NV, 5};
      auto a1 = xt::empty<double>(sh1);
      fill_h1d_field_xt(aomap, h1d_names, vars, NB, nbinmap, &Dbn1D::sumW      ,    &HistoBin1D::sumW,        a1, 0);
      fill_h1d_field_xt(aomap, h1d_names, vars, NB, nbinmap, &Dbn1D::sumW2     ,    &HistoBin1D::sumW2,       a1, 1);
      fill_h1d_field_xt(aomap, h1d_names, vars, NB, nbinmap, &Dbn1D::sumWX     ,    &HistoBin1D::sumWX,       a1, 2);
      fill_h1d_field_xt(aomap, h1d_names, vars, NB, nbinmap, &Dbn1D::sumWX2    ,    &HistoBin1D::sumWX2,      a1, 3);
      fill_h1d_field_xt(aomap, h1d_names, vars, NB, nbinmap, &Dbn1D::numEntries,    &HistoBin1D::numEntries,  a1, 4);
      H5Easy::dump(file, "Histo1D", a1, H5Easy::DumpOptions(H5Easy::Compression(compression)));

      // TODO make the attribute creation a function    
      HighFive::DataSet ds = file.getDataSet("Histo1D");
      Attribute a = ds.createAttribute<std::string>("binids", DataSpace::From(binids));
      a.write(binids);
      Attribute b = ds.createAttribute<std::string>("histonames", DataSpace::From(h1d_names));
      b.write(h1d_names);
      Attribute c = ds.createAttribute<std::string>("weightnames", DataSpace::From(vars));
      c.write(vars);
      vector<string> zinfo = {"sumW", "sumW2", "sumWX", "sumWX2", "numEntries"};
      Attribute d = ds.createAttribute<std::string>("zinfo", DataSpace::From(zinfo));
      d.write(zinfo);

      
      // Scaling information
      xt::xtensor<double, 2>::shape_type sh2 = {h1d_names.size(), NV};
      auto sb = xt::empty<double>(sh2);
      fillScaledBy(aomap, h1d_names, vars, sb);

      vector<std::string> hasScaledBy;
      for (auto hn : h1d_names)
      {
        if (aomap[hn]->hasAnnotation("ScaledBy")) hasScaledBy.push_back("true");
        else                                      hasScaledBy.push_back("false");
      }

      H5Easy::dump(file, "Histo1D_ScaledBy", sb, H5Easy::DumpOptions(H5Easy::Compression(compression)));

      HighFive::DataSet ds2 = file.getDataSet("Histo1D_ScaledBy");
      Attribute e = ds2.createAttribute<std::string>("hasScaledBy", DataSpace::From(hasScaledBy));
      e.write(hasScaledBy);

      // Bin edges
      file.createDataSet("Histo1D_xMin",   get_edge(aomap, h1d_names, &HistoBin1D::xMin));
      file.createDataSet("Histo1D_xMax",   get_edge(aomap, h1d_names, &HistoBin1D::xMax));
      
    return 0;
}
