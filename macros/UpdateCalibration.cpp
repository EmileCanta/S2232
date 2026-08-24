#include <string>
#include <filesystem>
#include <fstream>
#include <TROOT.h>
#include <TFile.h>
#include "/home/emile/GRSISort/include/TChannel.h"

void UpdateCalibration(int run_number_cal, int firstrun, int lastrun)
{
    std::string analysis_directory = "/home/emile/postdoc/data/s2232/sorted/analysis";
    std::string calibration_directory_files = "/home/emile/postdoc/analysis/s2232/sorter/calfiles/";

    char calfile[1024];
    snprintf(calfile, sizeof(calfile), "%senergy_calibration%d.cal",
              calibration_directory_files.c_str(), run_number_cal);
    std::cout << "Cal file: " << calfile << std::endl;

    for (int run_number_file = firstrun; run_number_file <= lastrun; run_number_file++)
    {
        char file_prefix[1024];
        snprintf(file_prefix, sizeof(file_prefix), "analysis%d", run_number_file);

        for (const auto & entry : std::filesystem::directory_iterator(analysis_directory))
        {
            std::string s1 = entry.path().string();
            if (s1.find(file_prefix) != std::string::npos)
            {
                std::string filename = entry.path().string();
                std::cout << "  Updating: " << filename << std::endl;

                // Re-open the cal file fresh for every target file
                std::ifstream run_file(calfile);
                if (!run_file.good()) {
                    std::cerr << "  Could not open cal file: " << calfile << std::endl;
                    continue;
                }

                auto file = TFile::Open(filename.c_str(), "update");
                if (!file || file->IsZombie()) {
                    std::cerr << "  Could not open ROOT file: " << filename << std::endl;
                    continue;
                }

                TChannel::ReadCalFile(run_file);
                TChannel::WriteToRoot();
                file->Close();

                run_file.close();
            }
        }
    }

    gROOT->ProcessLine(".q");
}
