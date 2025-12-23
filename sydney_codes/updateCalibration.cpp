#include <string>
#include <filesystem>
#include <fstream>
#include "/home/sydney/GRSISort/include/TChannel.h"
#include <TROOT.h>
#include <TFile.h>
#include "../include/GAOptions.h"

void updateCalibration(int run_number) {

    std::string analysis_directory = GAOptions::analysis_directory;
    std::string calibration_directory_files = GAOptions::calibration_directory_files;
    char calfile[1024];
    std::string filename;
    std::string s1;
    char file_prefix[1024];

	snprintf(file_prefix, 1024, "analysis%d",run_number);
    snprintf(calfile, 1024, "%senergy_calibration%d.cal", calibration_directory_files.c_str(), run_number);
    std::ifstream run_file(calfile);

    for (const auto & entry : std::filesystem::directory_iterator(analysis_directory)){
        s1 = entry.path();
        if (s1.find(file_prefix) != std::string::npos) {
			filename = entry.path();
            auto file = TFile::Open(filename.c_str(), "update");
            TChannel::ReadCalFile(run_file);
            TChannel::WriteToRoot();
            file->Close();
        }
    }

    run_file.close();
    gROOT->ProcessLine(".q");
}