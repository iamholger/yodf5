#include <iostream>
#include <fstream>
#include <Eigen/Dense>
#include <highfive/H5Easy.hpp>
#include <highfive/H5Group.hpp>
#include <unordered_map>

#include "YODA/Histo1D.h"
#include "YODA/Dbn1D.h"
#include "YODA/WriterYODA.h"

using namespace Eigen;

void replace_all(std::string& input, const std::string& from, const std::string& to) {
  size_t pos = 0;
  while ((pos = input.find(from, pos)) != std::string::npos) {
    input.replace(pos, from.size(), to);
    pos += to.size();
  }
}

//AppEval AppEval::readHDF(const char* fname) {
    //H5Easy::File file(fname, H5Easy::File::ReadOnly);
    //auto PC        = H5Easy::load<MatrixXd>(file, "pcoeff"   );
    //auto QC        = H5Easy::load<MatrixXd>(file, "qcoeff"   );
    //auto structure = H5Easy::load<MatrixXi>(file, "structure");
    //auto offset    = H5Easy::load<MatrixXd>(file, "a"        );
    //auto xmin      = H5Easy::load<MatrixXd>(file, "xmin"     );
    //auto scaleterm = H5Easy::load<MatrixXd>(file, "scaleterm");
    //return {PC, QC, structure, offset, xmin, scaleterm};
//}
//

std::unordered_map<std::string, std::vector<int> > getListOfContents(std::string const & fname, std::string const & dtype="Histo1D") {
    H5Easy::File file(fname, H5Easy::File::ReadOnly);
    HighFive::Group g = file.getGroup(dtype);
    std::vector<std::string> const names = g.listObjectNames();

    std::unordered_map<std::string, std::vector<int> > index;
    for (auto s: names) index[s] = H5Easy::load<std::vector<int> >(file, dtype+"/"+s);

    return index;
}

std::vector<std::string> getListOfVariations(std::string const & fname) {
    H5Easy::File file(fname, H5Easy::File::ReadOnly);
    auto vars = H5Easy::load<std::vector<std::string> >(file, "variations");
    return vars;
}

int main(int argc, const char* argv[]) {
    H5Easy::File file(argv[1], H5Easy::File::ReadOnly);
    auto sumw       = H5Easy::load<MatrixXd>(file, "sumw"      );
    auto sumw2      = H5Easy::load<MatrixXd>(file, "sumw2"     );
    auto sumwx      = H5Easy::load<MatrixXd>(file, "sumwx"     );
    auto sumwx2     = H5Easy::load<MatrixXd>(file, "sumwx2"    );
    auto sumwy      = H5Easy::load<MatrixXd>(file, "sumwy"     );
    auto sumwy2     = H5Easy::load<MatrixXd>(file, "sumwy2"    );
    auto sumwxy     = H5Easy::load<MatrixXd>(file, "sumwxy"    );
    auto numEntries = H5Easy::load<MatrixXd>(file, "numEntries");
    auto xerrm      = H5Easy::load<MatrixXd>(file, "xerr-"     );
    auto xerrp      = H5Easy::load<MatrixXd>(file, "xerr+"     );
    auto xval       = H5Easy::load<MatrixXd>(file, "xval"      );
    auto yerrm      = H5Easy::load<MatrixXd>(file, "yerr-"     );
    auto yerrp      = H5Easy::load<MatrixXd>(file, "yerr+"     );
    auto yval       = H5Easy::load<MatrixXd>(file, "yval"      );
    auto xmin       = H5Easy::load<MatrixXd>(file, "xmin"      );
    auto xmax       = H5Easy::load<MatrixXd>(file, "xmax"      );
    auto ymin       = H5Easy::load<MatrixXd>(file, "ymin"      );
    auto ymax       = H5Easy::load<MatrixXd>(file, "ymax"      );

    auto h1idx = getListOfContents(argv[1], "Histo1D"  );
    auto h2idx = getListOfContents(argv[1], "Histo2D"  );
    auto p1idx = getListOfContents(argv[1], "Profile1D");
    auto s1idx = getListOfContents(argv[1], "Scatter1D");
    auto s2idx = getListOfContents(argv[1], "Scatter2D");
    auto ctidx = getListOfContents(argv[1], "Counter"  );


    auto vars = getListOfVariations(argv[1]);

    std::vector<YODA::Histo1D> allH1;
    for( const auto& h1 : h1idx )
    {
        for (unsigned int ivar=0; ivar<vars.size();++ivar) 
        {
            auto hname = h1.first;
            replace_all(hname, "|", "/");
            if (ivar>0)
            {
                hname+="[";
                hname+=vars[ivar];
                hname+="]";
            }
            YODA::Histo1D hcurr(hname);
            std::vector<YODA::HistoBin1D> h1binscurr;
            h1binscurr.reserve(h1.second.size()-3);

            int nbin=0;
            for (auto i : h1.second) 
            {
                const YODA::Dbn1D dbn(numEntries(i, ivar), sumw(i, ivar), sumw2(i, ivar), sumwx(i, ivar), sumwx2(i, ivar)); //i is the bin number, col is the variation
                if      (nbin==0) hcurr.setTotalDbn(dbn);
                else if (nbin==1) hcurr.setOverflow(dbn);
                else if (nbin==2) hcurr.setUnderflow(dbn);
                else h1binscurr.push_back(YODA::HistoBin1D(std::make_pair(xmin(i, ivar), xmax(i, ivar)), dbn));

                ++nbin;
            }
            hcurr.addBins(h1binscurr);
            allH1.push_back(hcurr);
        }
    }
    
    std::vector<YODA::Histo2D> allH2;
    for( const auto& h2 : h2idx )
    {
        for (unsigned int ivar=0; ivar<vars.size();++ivar) 
        {
            auto hname = h2.first;
            replace_all(hname, "|", "/");
            if (ivar>0)
            {
                hname+="[";
                hname+=vars[ivar];
                hname+="]";
            }
            YODA::Histo2D hcurr(hname);
            std::vector<YODA::HistoBin2D> h2binscurr;
            h2binscurr.reserve(h2.second.size()-3);

            int nbin=0;
            for (auto i : h2.second) 
            {
                const YODA::Dbn2D dbn(numEntries(i, ivar), sumw(i, ivar), sumw2(i, ivar), sumwx(i, ivar), sumwx2(i, ivar), sumwy(i, ivar), sumwy2(i, ivar), sumwxy(i, ivar)); //i is the bin number, col is the variation
                if      (nbin==0) hcurr.setTotalDbn(dbn);
                //else if (nbin==1) continue;//hcurr.setOverflow(dbn);
                //else if (nbin==2) continue;//hcurr.setUnderflow(dbn);
                if (nbin>2) h2binscurr.push_back(YODA::HistoBin2D(std::make_pair(xmin(i, ivar), xmax(i, ivar)),std::make_pair(ymin(i, ivar), ymax(i, ivar)), dbn));


                ++nbin;
            }
            hcurr.addBins(h2binscurr);
            allH2.push_back(hcurr);
        }
    }
    
    std::vector<YODA::Profile1D> allP1;
    for( const auto& p1 : p1idx )
    {
        for (unsigned int ivar=0; ivar<vars.size();++ivar) 
        {
            auto hname = p1.first;
            replace_all(hname, "|", "/");
            if (ivar>0)
            {
                hname+="[";
                hname+=vars[ivar];
                hname+="]";
            }
            YODA::Profile1D hcurr(hname);
            std::vector<YODA::ProfileBin1D> p1binscurr;
            p1binscurr.reserve(p1.second.size()-3);

            int nbin=0;
            for (auto i : p1.second) 
            {
                const YODA::Dbn2D dbn(numEntries(i, ivar), sumw(i, ivar), sumw2(i, ivar), sumwx(i, ivar), sumwx2(i, ivar), sumwy(i, ivar), sumwy2(i, ivar), 0); //i is the bin number, col is the variation
                if      (nbin==0) hcurr.setTotalDbn(dbn);
                else if (nbin==1) hcurr.setOverflow(dbn);
                else if (nbin==2) hcurr.setUnderflow(dbn);
                else p1binscurr.push_back(YODA::ProfileBin1D(std::make_pair(xmin(i, ivar), xmax(i, ivar)), dbn));

                ++nbin;
            }
            hcurr.addBins(p1binscurr);
            allP1.push_back(hcurr);
        }
    }

    
    std::vector<YODA::Scatter2D> allS2;
    for( const auto& s2 : s2idx )
    {
        for (unsigned int ivar=0; ivar<vars.size();++ivar) 
        {
            auto hname = s2.first;
            replace_all(hname, "|", "/");
            if (ivar>0)
            {
                hname+="[";
                hname+=vars[ivar];
                hname+="]";
            }
            YODA::Scatter2D scurr(hname);
            std::vector<YODA::Point2D> s2binscurr;
            s2binscurr.reserve(s2.second.size());

            for (auto i : s2.second) 
            {
                YODA::Point2D thispoint=YODA::Point2D(xval(i, ivar), yval(i, ivar), xerrm(i, ivar), xerrp(i, ivar), yerrm(i, ivar), yerrp(i, ivar));
                s2binscurr.push_back(thispoint);
            }
            scurr.addPoints(s2binscurr);
            allS2.push_back(scurr);
        }
    }

    if (argc>2) return 0;

    std::ofstream myfile;
    myfile.open("example.yoda");
    for (auto h1 : allH1) YODA::WriterYODA::write(myfile, h1);
    for (auto h2 : allH2) YODA::WriterYODA::write(myfile, h2);
    for (auto p1 : allP1) YODA::WriterYODA::write(myfile, p1);
    for (auto s2 : allS2) YODA::WriterYODA::write(myfile, s2);
    myfile.close();

    return 0;
}
