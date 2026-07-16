#include <string>
#include <filesystem>
#include <fstream>
#include <TROOT.h>
#include <TFile.h>
#include "/home/emile/GRSISort/include/TChannel.h"

void UpdateCalibration(int run_number_cal, int run_number_file) 
{
    std::string analysis_directory = "/home/emile/postdoc/data/s2232/sorted/analysis";
    std::string calibration_directory_files = "/home/emile/postdoc/analysis/s2232/sorter/calfiles/";

    char calfile[1024];

    std::string filename;
    std::string s1;

    char file_prefix[1024];

    snprintf(file_prefix, 1024, "analysis%d",run_number_file);
    snprintf(calfile, 1024, "%senergy_calibration%d.cal", calibration_directory_files.c_str(), run_number_cal);

    std::cout << calfile << std::endl;
    std::ifstream run_file(calfile);

    for(const auto & entry : std::filesystem::directory_iterator(analysis_directory))
    {
        s1 = entry.path();

        if(s1.find(file_prefix) != std::string::npos) 
        {
            filename = entry.path();
            std::cout << filename.c_str() << std::endl;

            auto file = TFile::Open(filename.c_str(), "update");

            TChannel::ReadCalFile(run_file);
            TChannel::WriteToRoot();

            file->Close();
        }
    }

    run_file.close();

    gROOT->ProcessLine(".q");
}
