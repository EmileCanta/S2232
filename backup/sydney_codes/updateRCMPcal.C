#include <iostream>
#include <string>
#include "TChannel.h"
#include "TFile.h"
#include <filesystem>

namespace fs = std::filesystem;

void updateRCMPcal() 
{
    cout << "boop" << endl;

    // Variables to set //
  
    int refRun = 28226;
    int copyRun = 28226;
    int firstRun = 28227;     
    int lastRun = 28227;     
    int NDET = 6;
    int NStrips = 32;

    std::string calfile_directory = "/home/emile/postdoc/analysis/s2232/sorter/calfiles/"; 
    std::string analysis_directory = "/home/emile/postdoc/data/s2232/sorted/analysis/"; 
    std::string filename = "energy_calibration";

    std::map<int, std::map<int, std::pair<double,double>>> coefMapP; 
    std::map<int, std::map<int, std::pair<double,double>>> coefMapN; 

    // Open reference calibration file // 

    std::string refFile = calfile_directory + filename + std::to_string(refRun) + ".cal";

    TChannel::ReadCalFile(refFile.c_str());

    std::cout << "Loaded calibration file: " << refFile << std::endl;

    std::string copyFile = calfile_directory + filename + std::to_string(copyRun) + ".cal";

    // Save the coefficients in the map // 
 
    for (int det = 1; det <= NDET; det++) 
    {
        for (int strip = 0; strip < NStrips; strip++)
        {
            std::string chNameP = Form("RCS%02dXP%02dX", det, strip); // Front Strips
                                                                     
            TChannel* chanP = TChannel::FindChannelByName(chNameP.c_str());

            if (!chanP) 
            {   
                std::cerr << "Warning: channel not found: " << chNameP << std::endl;
                continue;
            }

            std::vector<Float_t> coeffsP = chanP -> GetENGCoeff();

            coefMapP[det][strip] = {coeffsP[0], coeffsP[1]};

            std::string chNameN = Form("RCS%02dXN%02dX", det, strip); // Back Strips
                                                                     
            TChannel* chanN = TChannel::FindChannelByName(chNameN.c_str());

            if (!chanN) 
            {   
                std::cerr << "Warning: channel not found: " << chNameN << std::endl;
                continue;
            }

            std::vector<Float_t> coeffsN = chanN -> GetENGCoeff();

            coefMapN[det][strip] = {coeffsN[0], coeffsN[1]};
        }
    }


    // Update all of the cal files // 
    for (int run = firstRun; run <= lastRun; run++)
    {
        if (run == refRun) continue; // Skip run 28226

        std::string currentFile = calfile_directory + filename + std::to_string(run) + ".cal";

        if(!fs::exists(currentFile))
        {
            std::cerr << "Calibration file not found for run " << run << ". Creating a copy from " << copyFile << "\n";

            // Copy Calibration File
            fs::copy(copyFile, currentFile, fs::copy_options::overwrite_existing);
            std::cout << "Created calibration file: " << currentFile << " from " << copyFile << std::endl;

        }

        TChannel::ReadCalFile(currentFile.c_str());
        std::cout << "Opened calibration file: " << currentFile << std::endl;

        for (int det = 1; det <= NDET; det++ )
        {
            for (int strip = 0; strip < NStrips; strip++)
            {
                // P Channels
                std::string chNameP = Form("RCS%02dXP%02dX", det, strip); 
                TChannel* chanP = TChannel::FindChannelByName(chNameP.c_str());
                if (chanP)
                {
                    chanP -> DestroyENGCal();
                    chanP -> AddENGCoefficient(coefMapP[det][strip].first);
                    chanP -> AddENGCoefficient(coefMapP[det][strip].second);
                }

                // N Channels
                std::string chNameN = Form("RCS%02dXN%02dX", det, strip); 
                TChannel* chanN = TChannel::FindChannelByName(chNameN.c_str());
                if (chanN)
                {
                    chanN -> DestroyENGCal();
                    chanN -> AddENGCoefficient(coefMapN[det][strip].first);
                    chanN -> AddENGCoefficient(coefMapN[det][strip].second);
                }
            }
        }

        // Save the changes back to the cal file
        TChannel::WriteCalFile(currentFile);
        std::cout << "Updated and saved calibration file: " << currentFile << std::endl;


        // Update the ROOT files
        char file_prefix[1024];
        snprintf(file_prefix, 1024, "analysis%d",run);

        for (const auto& entry : fs::directory_iterator(analysis_directory))
        {
            std::string path = entry.path();
            if (path.find(file_prefix) != std::string::npos && path.find(".root") != std::string::npos)
            {
                auto file = TFile::Open(path.c_str(), "update");
                if (!file || file -> IsZombie())
                {
                    std::cerr << "Could not open ROOT file: " << path << std::endl;
                    continue;
                }

                TChannel::ReadCalFile(currentFile.c_str());
                TChannel::WriteToRoot();

                file -> Close();
                std::cout << "Updated ROOT file with calibration: " << path << std::endl;
            }
        }


    }
}
