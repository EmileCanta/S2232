#include <string>
#include <chrono>
#include <iostream>
#include <fstream>
#include <filesystem>
#include <numeric>
#include <vector>
#include <map>
#include <limits>
#include <cmath>
#include <TROOT.h>
#include <TH1.h>
#include <TH2.h>
#include <TFile.h>
#include <TCanvas.h>
#include <TSpectrum.h>
#include <TGraph.h>
#include "/home/emile/GRSISort/include/TChannel.h"
#include "/home/emile/GRSISort/include/TRWPeak.h"
#include "./FindHistogramFile.h"

void RCMPEnergyCalProtonDet1_2(int run_number, int vis) 
{
    std::string log_directory = "./logs/";
    std::string histogram_directory = "/home/emile/postdoc/analysis/s2232/histograms/rcmp";
    std::string histogram_prefix = "RCMP";

    const int NStrips = 32;
    const int NDET = 6;
    int NPEAKS;

    std::vector<double> peak_energies;
    
    NPEAKS = 2;
    peak_energies.insert(peak_energies.end(), {741., 4088.}); //Proton

    int max_fit_attempts = 10;
    double peak_search_sigma = 0.7;

    double peak_search_threshold;
    peak_search_threshold = 0.55;

    double fit_success_tolerance = 15.;
    double max_plot_window = 100;
    double min_fit_radius = 20;
    double max_fit_radius = 40;
    double fit_radius;
    double charge_estimate_init = 0.;
    double minEntries = 1000;

    std::vector<std::vector<double>> peaksvec;

    std::ifstream gin;

    //gin.open("./charge_estimate_front_det1.dat"); 
    gin.open("./charge_estimate_front_det2.dat"); 

    if(!gin.is_open())
    {
        std::cout << "ERROR: Could not open rcmp_calib_coeff_estimates.dat\n";
    } 

    else 
    {
        double peak1 = 0.0;
        double peak2 = 0.0;

        int it = 0;

        while(gin >> peak1 >> peak2)
        {
            peaksvec.push_back({peak1, peak2});

            cout << peaksvec[it][0] << " " << peaksvec[it][1] << endl;

            it++;
        }

        gin.close();
    }

    TFile* my_file = TFile::Open(findHistogramFile(run_number, histogram_directory, histogram_prefix).c_str(), "read");

    my_file->Get("Channel");

    double peak_estimate;

    double min_search;
    double max_search;
    double *xpeaks;
    std::vector<double> xpeaks_full;
    xpeaks_full.assign(2,charge_estimate_init);   

    double xpeaks_fit;

    double strip_num[NStrips];

    bool fit_success;

    TCanvas *cFit = new TCanvas("cFit", "Peak Fits", 800, 600);
    
    TSpectrum* peak_search = new TSpectrum();

    for(int i = 2; i <= 2; i++) 
    { 
        double mean_chi_square;

        TH2D* py_charge;

        py_charge = (TH2D*)my_file->Get(Form("ChargeVSFrontStrip_MultTwo%d", i)); 

        for(int strip = 0; strip < NStrips; strip++) 
        { 
            TH1D* channel_charge = (TH1D*)py_charge->ProjectionY("channel_charge", strip+1, strip+1); 

            if(channel_charge->GetEntries() < minEntries)
			{
                std::cout << "Not enough statistics in det: " << i << " strip: " << strip+1 << std::endl;

				continue;
			}

            for(int j = 0; j < NPEAKS; j++) 
            {
                channel_charge->GetXaxis()->SetRangeUser(peaksvec[strip][j]-max_plot_window, peaksvec[strip][j]+max_plot_window);

                cout << "Range : " << peaksvec[strip][j]-max_plot_window << " - " << peaksvec[strip][j]+max_plot_window << endl;

                if(vis == 1)
                {
                    cFit->Update();
                    gSystem->ProcessEvents();
                    cFit->cd();
                    cFit->Update();
                    gSystem->ProcessEvents();
                    channel_charge->Draw();
                    cFit->Update();
                    gSystem->ProcessEvents();
                    gSystem->Sleep(2000);
                    cFit->Update();
                    gSystem->ProcessEvents();
                }
                
                peak_search->Search(channel_charge, peak_search_sigma, "nodraw", peak_search_threshold);

                xpeaks = peak_search->GetPositionX();

                int npeaks_found = peak_search->GetNPeaks();

                cout << "Number of peaks found: " << npeaks_found << " Charge 1: " << xpeaks[0] << " Charge 2: " << xpeaks[1] << endl;
                
                if(npeaks_found <= 0) 
                {
                    std::cout << "Not peaks found in det: " << i << " strip: " << strip+1 << " peak: " << j << std::endl;

                    continue;
				}
                
                peak_estimate = xpeaks[0];

                peaksvec[strip][j] = peak_estimate;

                channel_charge->GetXaxis()->SetRangeUser(peak_estimate-max_plot_window, peak_estimate+max_plot_window);

                if(channel_charge->Integral(channel_charge->FindBin(peak_estimate-max_plot_window), channel_charge->FindBin(peak_estimate+max_plot_window)) <= 100) 
                {
                    std::cout << "Not peaks found in det: " << i << " strip: " << strip+1 << " peak: " << j << std::endl;

                    continue;
                }

                int fit_attempts = 0;

                fit_success = false;
                fit_radius = max_fit_radius;

                while(!fit_success && fit_attempts < max_fit_attempts) 
                { 
                    fit_radius = fit_radius-(max_fit_radius-min_fit_radius)/(max_fit_attempts-1);
                    
                    min_search = peaksvec[strip][j]-fit_radius;
                    max_search = peaksvec[strip][j]+fit_radius;

                    TF1* peak_fitter;
                    
                    peak_fitter = new TF1("peak_fitter", "[0]*exp(-0.5*((x-[1])/[2])*((x-[1])/[2])) + [3]*x + [4]", min_search, max_search);

                    peak_fitter->SetParameters(channel_charge->GetMaximum(), channel_charge->GetBinCenter(channel_charge->GetMaximumBin()), channel_charge->GetRMS());
                    
                    TFitResultPtr fit_result = channel_charge->Fit(peak_fitter, "RQ0S");

                    if(vis == 1)
                    {
                        cFit->Update();
                        gSystem->ProcessEvents();
                        cFit->cd();
                        cFit->Update();
                        gSystem->ProcessEvents();
                        channel_charge->Draw();
                        cFit->Update();
                        gSystem->ProcessEvents();
                        peak_fitter->Draw("same");
                        cFit->Update();
                        gSystem->ProcessEvents();
                        cFit->Update();
                        gSystem->ProcessEvents();
                        gSystem->Sleep(20);
                        cFit->Update();
                        gSystem->ProcessEvents();
                    }

                    xpeaks_fit = fit_result->Parameter(1);

                    cout << "Fit charge: " << xpeaks_fit << endl << "\n";

                    if(abs(peak_estimate-xpeaks_fit) < fit_success_tolerance) 
                    {
                        fit_success = true;
                    }

                    else 
                    {
                        fit_attempts++;
                    }
                }

                if(fit_success)
                {
                    xpeaks_full[j] = xpeaks_fit;
                }

                else 
                {
                    xpeaks_full[j] = peaksvec[strip][j];
                }

                if(vis == 1)
                {
                    cFit->Update();
                    gSystem->ProcessEvents();
                    gSystem->Sleep(2000);
                    cFit->Update();
                    gSystem->ProcessEvents();
                }
            }

            TGraph* fit_graph;
            fit_graph = new TGraph(2, xpeaks_full.data(), peak_energies.data());

            TF1* linear_fit;
            linear_fit = new TF1("linear_fit", "[0]+[1]*x", xpeaks_full[0], xpeaks_full[1]);

            fit_graph->Fit(linear_fit, "qr");

            double offset = linear_fit->GetParameter(0), e_offset = linear_fit->GetParError(0);
            double gain = linear_fit->GetParameter(1), e_gain = linear_fit->GetParError(1);

            std::string channel_name;

            channel_name = Form("RCS%02dXP%02dX", i, strip); 
                                                                         
            TChannel* detector_channel = TChannel::FindChannelByName(channel_name.c_str());

            if(!detector_channel) 
            {
                std::cout << "Channel not found: " << channel_name << std::endl;

                continue; 
            }

            detector_channel->DestroyENGCal();
            detector_channel->AddENGCoefficient(offset);
            detector_channel->AddENGCoefficient(gain);

            std::cout << "Det : " << i << "    Strip : " << strip+1 << "    New gain : " << gain << "    New offset : " << offset << std::endl << "\n\n\n";
            
            std::string calibration_file = "/home/emile/postdoc/analysis/s2232/sorter/calfiles/energy_calibration" + std::to_string(run_number) + ".cal";

            TChannel::WriteCalFile(calibration_file);
        }
    }

    my_file->Close();

    gROOT->ProcessLine(".q");
}
