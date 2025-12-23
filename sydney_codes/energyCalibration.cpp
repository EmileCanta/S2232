#include <string>
#include <chrono>
#include <iostream>
#include <fstream>
#include <filesystem>
#include <numeric>
#include <vector>
#include <TROOT.h>
#include <TH1.h>
#include <TH2.h>
#include <TFile.h>
#include <TCanvas.h>
#include <TSpectrum.h>
#include <TGraph.h>
#include "/home/sydney/GRSISort/include/TChannel.h"
#include "/home/sydney/GRSISort/include/TRWPeak.h"
#include "/home/sydney/GRSISort/GRSIData/include/TGriffin.h"
#include "../include/FindHistogramFile.h"
#include "../include/GAOptions.h"

void energyCalibration(int run_number, std::string options = "f") {

	// Options
	bool auto_quit = options.find("q") != std::string::npos;
	bool hide_plots = options.find("h") != std::string::npos;
	bool do_fit = options.find("f") != std::string::npos;
	if (hide_plots || auto_quit){
		gROOT->SetBatch();
	}

	// Parameters to be set by user
	std::string log_directory = GAOptions::log_directory;
	std::string histogram_directory = GAOptions::histogram_directory;
	std::string histogram_prefix = GAOptions::alphanumeric_prefix;
	const int NDET = GAOptions::number_of_detectors;
	const int NPEAKS = 3;
	double peak_energies[NPEAKS]{238.27, 984.25, 1633.674}; //2614.522
	int max_fit_attempts = 4;
	double peak_search_sigma = 0.7;
	double peak_search_threshold = 0.75;
	double fit_success_tolerance = 1.5;
	double min_fit_radius = 2.5;
	double max_fit_radius = 36.0;
	double fit_radius;
	double charge_estimate_init = 186.0;

	// Create log file
	auto start_time = std::chrono::system_clock::now();
	std::time_t start_time_t = std::chrono::system_clock::to_time_t(start_time);
	std::stringstream str; 
	char log_file_name[1024];
	char time_string[32];
	strftime(time_string, 32, "%Y-%m-%d_%H:%M:%S", localtime(&start_time_t));
	snprintf(log_file_name, 1024, "%senergy_calibration%d_%s.txt", log_directory.c_str(), run_number, time_string);
	std::ofstream tout_log(log_file_name);
	tout_log << "Script executed at " << std::ctime(&start_time_t) << std::endl <<
	"Run = " << run_number << std::endl;

	// Get charge vs. channel histogram
	TFile* my_file = TFile::Open(findHistogramFile(run_number, histogram_directory, histogram_prefix).c_str(), "read");
	//my_file -> ls();
    TH2D* py_charge = (TH2D*)my_file->Get("ChargeVsChannel");
	my_file->Get("Channel");

	double charge_estimates[NPEAKS]{charge_estimate_init}; //Sets array of {186, 0, 0}
	double gain_estimate;
	double peak_estimate;

	// Variables for the peak finding algorithm
	double min_search;
	double max_search;
	double *xpeaks;
	double xpeaks_full[NPEAKS]; // make an array with numbers ot be filled in (shows up as zero initally)
	double xpeaks_fit;


	double crystal_num[NDET];
	double chisquares[NDET];
	double mean_chi_square;

	bool fit_success;
	int fit_success_counter = 0;
	int fit_success_peaks_counter[NPEAKS]{}; // sets values in arrays equal to zero

    for (int i = 0; i < NDET; i++) { //looping through all 64 crystals in data

	    TH1D* channel_charge = (TH1D*)py_charge->ProjectionY("channel_charge", i+1, i+1); //choses the 1D histogram to look at
		str << "\nFound peaks of channel " << i+1 << " at locations: "<<std::endl;

	    for (int j = 0; j < NPEAKS; j++) {

			// std::cout << "At j = " << j <<  "At i = " << i<< ", xpeaks_full[j] = " << xpeaks_full[j] << std::endl;
			if (j != 0) { // Figures out approximately where the peak is located
				gain_estimate = xpeaks_full[j-1]/peak_energies[j-1];
				charge_estimates[j] = gain_estimate*peak_energies[j];
			}
			else {
				charge_estimates[j] = charge_estimate_init;
			}
			//ask how the above loop works!
			//Lists the peaks from prominate height ex. xpeaks[0] would be the tallest (most likely 238)
	    
	    	channel_charge->GetXaxis()->SetRangeUser(charge_estimates[j]-max_fit_radius, charge_estimates[j]+max_fit_radius);
	    	TSpectrum* peak_search = new TSpectrum();
	    	peak_search->Search(channel_charge, peak_search_sigma, "nodraw", peak_search_threshold); // Cant find info about Search() online. appear Search()
			xpeaks = peak_search->GetPositionX();

			if (j!=0 && (abs(xpeaks[1]-charge_estimates[j]) < abs(xpeaks[0]-charge_estimates[j]))) { // if j doesn't = 0 AND the other condition is meet then .... 
				peak_estimate = xpeaks[1];
				cout << " Crystal: "<< i+1 << " triggered"; //doesnt seem like this is printing anywhere so is this condition never being met??
			}
			else if (j == 0 && (xpeaks[1] < xpeaks[0] && xpeaks[1] > charge_estimates[j]-max_fit_radius)) {
				peak_estimate = xpeaks[1];
				//cout << "this is working";
			}
			else {
				peak_estimate = xpeaks[0];
				// cout << "test2"; //only this one is working 
			}
			
			// if (j==0) {
			// 	cout << "Crystal " << i+1 << ":" << endl;
			// 	cout << xpeaks[0] << endl;
			// 	cout << xpeaks[1] << endl;
			// }
			
			charge_estimates[j] = peak_estimate;

			if (do_fit) {
				int fit_attempts = 0;
				fit_success = false;
				fit_radius = max_fit_radius;
				while (!fit_success && fit_attempts < max_fit_attempts) { // ! is logical NOT operator, ie. fit_sucess = true
					fit_radius = fit_radius-(max_fit_radius-min_fit_radius)/(max_fit_attempts-1);
					min_search = charge_estimates[j]-fit_radius;
					max_search = charge_estimates[j]+fit_radius;
					TRWPeak* peak1 = new TRWPeak(charge_estimates[j]); // RW - Raw Ware - basically gaussian but fancier
					TPeakFitter* peak_fitter = new TPeakFitter(min_search, max_search);
					peak_fitter->AddPeak(peak1);
					TFitResultPtr fit_result = peak_fitter->Fit(channel_charge, "q");
					xpeaks_fit = fit_result->Parameter(1);
					if (abs(peak_estimate-xpeaks_fit) < fit_success_tolerance) {
						fit_success = true;
					}
					else {
						fit_attempts++;
					}
				}
			}

			if (do_fit && fit_success){
				xpeaks_full[j] = xpeaks_fit;
				str << xpeaks_fit << "\t";
				fit_success_counter++;
				fit_success_peaks_counter[j]++;
				std::cout << "fit worked!" << std::endl;
			}
			else {
				xpeaks_full[j] = peak_estimate;
				str << peak_estimate << "\t";
				std::cout << "fit didnt work!" << std::endl;
			}
	    }
		
		// Perform linear fit to the peaks
		TGraph* fit_graph = new TGraph(NPEAKS, xpeaks_full, peak_energies);
	    TF1* linear_fit = new TF1("linear_fit", "[0]+[1]*x", xpeaks_full[0], xpeaks_full[NPEAKS-1]);
	    fit_graph->Fit(linear_fit, "qr");
		double offset = linear_fit->GetParameter(0), e_offset = linear_fit->GetParError(0);
	    double gain = linear_fit->GetParameter(1), e_gain = linear_fit->GetParError(1);
		double chi_square = linear_fit->GetChisquare();
		chisquares[i] = chi_square;

		// Write calibration coefficients
		int color_num = i/4+1;
		TChannel* detector_channel = TChannel::FindChannelByName(Form("GRG%02d%sN00A", color_num, TGriffin::GetColorFromNumber(i%4)));
		if (!detector_channel) {
				std::cerr << "Channel not found: " << Form("GRG%02d%sN00A", color_num, TGriffin::GetColorFromNumber(i%4)) << std::endl;
				continue; // skip to the next strip
			}
		detector_channel->DestroyENGCal();
 		detector_channel->AddENGCoefficient(offset);
 		detector_channel->AddENGCoefficient(gain);
	    str << "\n\nGain is: " << gain << " ± " << e_gain << std::endl <<
		"Offset is: " << offset << " ± " << e_offset << std::endl <<
		"Chi2 is: " << chi_square << "\n" << std::endl <<
		"~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~" << std::endl;
    }

	// Write new .cal file
	std::string calibration_file = GAOptions::calibration_directory_files + "energy_calibration" + std::to_string(run_number) + ".cal";
	TChannel::WriteCalFile(calibration_file);
	my_file->Close();

	// Plot reduced chi squared values for each crystal's fit
	std::iota(crystal_num, crystal_num+NDET, 1);
	TCanvas* c0 = new TCanvas("c0", "c0", 1000, 400);
	TGraph* chi_square_graph = new TGraph(NDET, crystal_num, chisquares);
	chi_square_graph->SetTitle("Chi squared values;Channel number;Chi^2");
	chi_square_graph->Draw("ap*");

	// Output statistics
	int chi_size = sizeof(chisquares) / sizeof(chisquares[0]); 
	std::vector<double> vectorchisquares(chisquares, chisquares+chi_size);
	mean_chi_square = TMath::Mean(vectorchisquares.begin(), vectorchisquares.end());
	str << "\nMean chi square: " << mean_chi_square << std::endl;
	if (do_fit) {
		double success_rate = fit_success_counter/(1.0*NPEAKS*NDET);
		str << "\nTotal fit success rate = " << 100*success_rate << "%" << std::endl;
		double success_rate_peaks[NPEAKS];
		for (int i = 0; i < 4; i++) {
			success_rate_peaks[i] = fit_success_peaks_counter[i]/(1.0*NDET);
			str << "Peak " << i+1 << " success rate = " << 100*success_rate_peaks[i] << "%" << std::endl;
		}
	}

	std::cout << str.str();
    tout_log << str.str();
	tout_log.close();

	if (auto_quit) {
		gROOT->ProcessLine(".q");
	}
	//gROOT -> ProcessLine(".q");
	
}