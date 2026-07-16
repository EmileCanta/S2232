#include "/home/emile/GRSISort/include/TChannel.h"
#include "./FindHistogramFile.h"

#include <iostream>
#include <fstream>

#include <TFile.h>

void GetRoughCalibCoeffs(int run_number) 
{
	std::string histogram_directory = "/home/emile/postdoc/analysis/s2232/histograms/rcmp";
	std::string histogram_prefix = "RCMP";

    std::ofstream outfile("rcmp_calib_coeff_estimates_front.dat"); //Here decide front or back
    //std::ofstream outfile("rcmp_calib_coeff_estimates_back.dat"); //Here decide front or back

    if(!outfile) 
    { 
        std::cerr << "Error opening file.\n";
    }

    TFile* my_file = TFile::Open(findHistogramFile(run_number, histogram_directory, histogram_prefix).c_str(), "read");

    my_file->Get("Channel");

    outfile << "Detector" << '\t' << "Strip" << '\t' << "Offset" << '\t' << "Gain" <<'\n';
    
    for(int i = 1; i <= 6; i++) 
    { 
        for(int j = 0; j < 32; j++) 
        {

            std::string channel_name = Form("RCS%02dXP%02dX", i, j); //Here decide front or back
            //std::string channel_name = Form("RCS%02dXN%02dX", i, j); //Here decide front or back
            TChannel* detector_channel = TChannel::FindChannelByName(channel_name.c_str());

            if(!detector_channel) 
            {
                std::cerr << "Channel not found!" << std::endl;
                continue;            
            }

            double offset = detector_channel->GetENGCoeff()[0];
            double gain = detector_channel->GetENGCoeff()[1];

            if(i == 1) //This for the detector swap made in the calfile (because I didnt swap the coefficients, only the addresses)
            {
                outfile << "2" << '\t' << j+1 << '\t' << offset << '\t' << gain << '\n';
            }

            if(i == 2)
            {
                outfile << "1" << '\t' << j+1 << '\t' << offset << '\t' << gain << '\n';
            }

            if( 2 < i && i < 7)
            {
                outfile << i << '\t' << j+1 << '\t' << offset << '\t' << gain << '\n';
            }
        }
    }

    outfile.close();
}
