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

void RCMPEnergyCal(const std::string& run_type, int run_number, const std::string& side, int vis) 
{
    std::string log_directory = "./logs/";
    std::string histogram_directory = "/home/emile/postdoc/analysis/s2232/histograms/rcmp";
    std::string histogram_prefix = "RCMP";

    const int NStrips = 32;
    const int NDET = 6;
    int NPEAKS;

    std::vector<double> peak_energies;

    if(run_type == "source")
    {
        NPEAKS = 3;
        peak_energies.insert(peak_energies.end(), {5156.59, 5485.56, 5804.77}); 
    }
    
    if(run_type == "beam")
    {
        NPEAKS = 2;
        peak_energies.insert(peak_energies.end(), {2153.6, 4434.7}); //Alpha 
    }

    int max_fit_attempts = 100;
    double peak_search_sigma = 0.7;

    double peak_search_threshold;
    if(run_type == "source") peak_search_threshold = 0.75;
    if(run_type == "beam") peak_search_threshold = 0.55;

    double fit_success_tolerance = 15.;
    double max_plot_window = 150;
    double max_plot_window1 = 300;
    double max_plot_window2 = 500;
    double min_fit_radius = 20;
    double max_fit_radius = 40;
    double fit_radius;
    double charge_estimate_init = 0.;
    double minEntries = 1000;

    std::map<int, std::vector<double>> offsetMap;
    std::map<int, std::vector<double>> gainMap;

    gainMap[1] = std::vector<double>(NStrips, std::numeric_limits<double>::quiet_NaN());
    gainMap[2] = std::vector<double>(NStrips, std::numeric_limits<double>::quiet_NaN());
    gainMap[3] = std::vector<double>(NStrips, std::numeric_limits<double>::quiet_NaN());
    gainMap[4] = std::vector<double>(NStrips, std::numeric_limits<double>::quiet_NaN());
    gainMap[5] = std::vector<double>(NStrips, std::numeric_limits<double>::quiet_NaN());
    gainMap[6] = std::vector<double>(NStrips, std::numeric_limits<double>::quiet_NaN());

    offsetMap[1] = std::vector<double>(NStrips, std::numeric_limits<double>::quiet_NaN());
    offsetMap[2] = std::vector<double>(NStrips, std::numeric_limits<double>::quiet_NaN());
    offsetMap[3] = std::vector<double>(NStrips, std::numeric_limits<double>::quiet_NaN());
    offsetMap[4] = std::vector<double>(NStrips, std::numeric_limits<double>::quiet_NaN());
    offsetMap[5] = std::vector<double>(NStrips, std::numeric_limits<double>::quiet_NaN());
    offsetMap[6] = std::vector<double>(NStrips, std::numeric_limits<double>::quiet_NaN());

    std::ifstream gin;

    if(side == "front") gin.open("/home/emile/postdoc/analysis/s2232/macros/rcmp_calib_coeff_estimates_front.dat"); 
    if(side == "back") gin.open("/home/emile/postdoc/analysis/s2232/macros/rcmp_calib_coeff_estimates_back.dat");  

    if(!gin.is_open())
    {
        std::cerr << "ERROR: Could not open rcmp_calib_coeff_estimates.dat\n";
    } 

    else 
    {
        std::string header;
        std::getline(gin, header);

        int det = 0;
        int strip = 0;
        double offset = 0.0;
        double gain = 0.0;

        while(gin >> det >> strip >> offset >> gain)
        {
            if(strip >= 1 && strip <= NStrips) 
            {
                offsetMap[det][strip - 1] = offset;
                gainMap[det][strip - 1] = gain;
            }
        }

        gin.close();
    }

    TFile* my_file = TFile::Open(findHistogramFile(run_number, histogram_directory, histogram_prefix).c_str(), "read");

    my_file->Get("Channel");

    std::vector<double> charge_estimates;
    if(run_type == "source") charge_estimates.assign(3,charge_estimate_init);
    if(run_type == "beam") charge_estimates.assign(2,charge_estimate_init);

    double gain_estimate;
    double peak_estimate;

    double min_search;
    double max_search;
    double *xpeaks;
    std::vector<double> xpeaks_full;
    if(run_type == "source") xpeaks_full.assign(3,charge_estimate_init);   
    if(run_type == "beam") xpeaks_full.assign(2,charge_estimate_init);   

    double xpeaks_fit;

    double strip_num[NStrips];

    bool fit_success;

    TCanvas *cFit = new TCanvas("cFit", "Peak Fits", 800, 600);
    
    TSpectrum* peak_search = new TSpectrum();

    for(int i = 5; i <= NDET; i++) 
    { 
        double mean_chi_square;

        TH2D* py_charge;

        if(side == "front") py_charge = (TH2D*)my_file->Get(Form("ChargeVSFrontStrip%d", i)); 
        if(side == "back") py_charge = (TH2D*)my_file->Get(Form("ChargeVSBackStrip%d", i)); 

        for(int strip = 0; strip < NStrips; strip++) 
        { 
            double my_gain = gainMap[i][strip]; 
            double my_offset = offsetMap[i][strip];

            TH1D* channel_charge = (TH1D*)py_charge->ProjectionY("channel_charge", strip+3, strip+3); 

            if(channel_charge->GetEntries() < minEntries)
			{ 
				continue;
			}

            for(int j = 0; j < NPEAKS; j++) 
            {

                charge_estimates[j] = (peak_energies[j] - my_offset) / my_gain;

                cout << "Peak number: " << j << " Energy: " << peak_energies[j] <<  "    Gain: " << my_gain << "    Offset: " << my_offset << "    Estimated charge: " << charge_estimates[j] << endl;

                double this_max_plotwindow = 0.;

                if(j == 0) this_max_plotwindow = max_plot_window1;
                if(j == 1) this_max_plotwindow = max_plot_window2;

                if(run_type == "source") 
                {
                    channel_charge->GetXaxis()->SetRangeUser(charge_estimates[j]-max_plot_window, charge_estimates[j]+max_plot_window);

                    cout << "Range : " << charge_estimates[j]-max_plot_window << " - " << charge_estimates[j]+max_plot_window << endl;
                }
                
                if(run_type == "beam") 
                {
                    channel_charge->GetXaxis()->SetRangeUser(charge_estimates[j]-this_max_plotwindow, charge_estimates[j]+this_max_plotwindow);

                    cout << "Range : " << charge_estimates[j]-this_max_plotwindow << " - " << charge_estimates[j]+this_max_plotwindow << endl;
                }

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
                    gSystem->Sleep(20);
                    cFit->Update();
                    gSystem->ProcessEvents();
                }
                
                
                peak_search->Search(channel_charge, peak_search_sigma, "nodraw", peak_search_threshold);

                xpeaks = peak_search->GetPositionX();

                int npeaks_found = peak_search->GetNPeaks();

                cout << "Number of peaks found: " << npeaks_found << " Charge 1: " << xpeaks[0] << " Charge 2: " << xpeaks[1] << endl;
                
                if(npeaks_found <= 0) 
                {
                    continue;
				}
                
                peak_estimate = xpeaks[0];

                charge_estimates[j] = peak_estimate;

                if(run_type == "source")
                {
                    channel_charge->GetXaxis()->SetRangeUser(peak_estimate-max_plot_window, peak_estimate+max_plot_window);

                    if(channel_charge->Integral(channel_charge->FindBin(peak_estimate-max_plot_window), channel_charge->FindBin(peak_estimate+max_plot_window)) <= 100) continue;
                }

                if(run_type == "beam")
                {
                    channel_charge->GetXaxis()->SetRangeUser(peak_estimate-this_max_plotwindow, peak_estimate+this_max_plotwindow);

                    if(channel_charge->Integral(channel_charge->FindBin(peak_estimate-this_max_plotwindow), channel_charge->FindBin(peak_estimate+this_max_plotwindow)) <= 100) continue;
                }

                int fit_attempts = 0;

                fit_success = false;
                fit_radius = max_fit_radius;

                while(!fit_success && fit_attempts < max_fit_attempts) 
                { 
                    fit_radius = fit_radius-(max_fit_radius-min_fit_radius)/(max_fit_attempts-1);
                    
                    min_search = charge_estimates[j]-fit_radius;
                    max_search = charge_estimates[j]+fit_radius;

                    TF1* peak_fitter;
                    
                    if(run_type == "source") peak_fitter = new TF1("peak_fitter", "[0]*exp(-0.5*((x-[1])/[2])*((x-[1])/[2])) + [3]*x + [4]", min_search, max_search);
                    if(run_type == "beam") peak_fitter = new TF1("peak_fitter", "[0]*exp(-0.5*((x-[1])/[2])*((x-[1])/[2]))", min_search, max_search);

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
                    xpeaks_full[j] = xpeaks_fit;
                }

                if(vis == 1)
                {
                    cFit->Update();
                    gSystem->ProcessEvents();
                    gSystem->Sleep(20);
                    cFit->Update();
                    gSystem->ProcessEvents();
                }
            }

            TGraph* fit_graph;
            if(run_type == "source") fit_graph = new TGraph(3, xpeaks_full.data(), peak_energies.data());
            if(run_type == "beam") fit_graph = new TGraph(2, xpeaks_full.data(), peak_energies.data());

            TF1* linear_fit;
            if(run_type == "source") linear_fit = new TF1("linear_fit", "[0]+[1]*x", xpeaks_full[0], xpeaks_full[2]);
            if(run_type == "beam") linear_fit = new TF1("linear_fit", "[0]+[1]*x", xpeaks_full[0], xpeaks_full[1]);

            fit_graph->Fit(linear_fit, "qr");

            double offset = linear_fit->GetParameter(0), e_offset = linear_fit->GetParError(0);
            double gain = linear_fit->GetParameter(1), e_gain = linear_fit->GetParError(1);

            std::string channel_name;

            if(side == "front") channel_name = Form("RCS%02dXP%02dX", i, strip); 
            if(side == "back") channel_name = Form("RCS%02dXN%02dX", i, strip); 
                                                                         
            TChannel* detector_channel = TChannel::FindChannelByName(channel_name.c_str());

            if(!detector_channel) 
            {
                std::cerr << "Channel not found: " << channel_name << std::endl;

                continue; 
            }

            detector_channel->DestroyENGCal();
            detector_channel->AddENGCoefficient(offset);
            detector_channel->AddENGCoefficient(gain);

            std::cout << "Det : " << i << "    Strip : " << strip+1 << "    New gain : " << gain << "    New offset : " << offset << std::endl << "\n\n\n";
        }

        std::string calibration_file = "/home/emile/postdoc/analysis/s2232/sorter/calfiles/energy_calibration" + std::to_string(run_number) + ".cal";

        TChannel::WriteCalFile(calibration_file);
    }

    my_file->Close();

    gROOT->ProcessLine(".q");
}
