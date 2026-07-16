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
#include "/home/sydney/GRSISort/include/TChannel.h"
#include "/home/sydney/GRSISort/include/TRWPeak.h"
#include "/home/sydney/GRSISort/GRSIData/include/TGriffin.h"
#include "../include/FindHistogramFile.h"
#include "../include/RCMPOptions .h"

void energyCal_RCMP_Alpha(int run_number, std::string options = "f") {


	// Options
	bool auto_quit = options.find("q") != std::string::npos;
	bool hide_plots = options.find("h") != std::string::npos;
	bool do_fit = options.find("f") != std::string::npos;
	if (hide_plots || auto_quit){
		gROOT->SetBatch();
	}

	// Parameters to be set by user
	std::string log_directory = RCMPOptions::log_directory;
	std::string histogram_directory = RCMPOptions::histogram_directory;
	std::string histogram_prefix = RCMPOptions::cal_prefix;
	const int NStrips = RCMPOptions::number_of_strips;
	const int NDET = RCMPOptions::number_of_detectors;
	const int NPEAKS = 2;
	double peak_energies[NPEAKS]{2153.6, 4434.7}; //2614.522
	int max_fit_attempts = 4;
	double peak_search_sigma = 0.7;
	double peak_search_threshold = 0.75;
	double fit_success_tolerance = 1.5;
	double min_fit_radius = 60;
	double max_fit_radius = 140;
	double fit_radius;
	double charge_estimate_init = 985.0;
	double minEntries = 1000;

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

	// -------- Read gains from text file --------
	std::map<int, std::vector<double>> offsetMap;
	std::map<int, std::vector<double>> gainMap;
	// initialize detector 5 and 6
	gainMap[5] = std::vector<double>(NStrips, std::numeric_limits<double>::quiet_NaN());
	gainMap[6] = std::vector<double>(NStrips, std::numeric_limits<double>::quiet_NaN());

	offsetMap[5] = std::vector<double>(NStrips, std::numeric_limits<double>::quiet_NaN());
	offsetMap[6] = std::vector<double>(NStrips, std::numeric_limits<double>::quiet_NaN());

	std::ifstream gin("/home/sydney/S2232_analysis/cal_files/files/RCMP_gains_28267.txt");  
	if(!gin.is_open()){std::cerr << "ERROR: Could not open RCMP_gains_28267.txt\n";} 
	
	else {
		// skip first text/header line
		std::string header;
		std::getline(gin, header);

		int det = 0;
		int strip = 0;
		double offset = 0.0;
		double gain = 0.0;

		while(gin >> det >> strip >> offset >> gain) {

			if((det == 5 || det == 6) &&
			strip >= 1 && strip <= NStrips) {

				offsetMap[det][strip - 1] = offset;
				gainMap[det][strip - 1] = gain;
			}
		}

		gin.close();
	}
	

	// Optional: quick sanity check
	//std::cout << "Gain for crystal 64 (index 63) = " << gains[63] << std::endl;

	// Open ROOT file
	TFile* my_file = TFile::Open(findHistogramFile(run_number, histogram_directory, histogram_prefix).c_str(), "read");
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


	double strip_num[NStrips];
	//double chisquares[NStrips];
	// double chisquares_D5[NStrips];
	// double chisquares_D6[NStrips];
	//double mean_chi_square;
	std::map<int, std::vector<double>> chiSquareMap;
	chiSquareMap[5] = std::vector<double>(NStrips, 0.0);
	chiSquareMap[6] = std::vector<double>(NStrips, 0.0);

	bool fit_success;
	int fit_success_counter = 0;
	int fit_success_peaks_counter[NPEAKS]{}; // sets values in arrays equal to zero

    for (int i = 5; i <= NDET; i++) { 
		
		double mean_chi_square;
		
		// Get charge vs. channel histogram
		//my_file -> ls();
		TH2D* py_charge = (TH2D*)my_file->Get(Form("Det%dFECharge", i));
		//my_file->Get("Channel");
		
		// double chisquares[NStrips];
		// if (i == 5){
		// 	chi
		// }
		
		for (int strip = 0; strip < NStrips; strip++) { //looping through all 64 crystals in data

			double my_gain = gainMap[i][strip];  // i is 0..63 in your loop
			double my_offset = offsetMap[i][strip];
			//charge_estimate_init = (peak_energies[0]-my_offset) / my_gain; // initial guess for peak position based on gain and first peak energy

			TH1D* channel_charge = (TH1D*)py_charge->ProjectionY("channel_charge", strip+1, strip+1); //choses the 1D histogram to look at
			// === Check to make sure that the strip is active === // 
			if (channel_charge->GetEntries() < minEntries)
			{ 
				str << "Skipping strip " << strip+1 << " due to low statistics (" << channel_charge->GetEntries() << " entries)\n";
				continue;
			}

			str << "\nFound peaks of Detector " << i << " Strip " << strip+1 << " at locations: "<<std::endl;

			for (int j = 0; j < NPEAKS; j++) {

				charge_estimates[j] = (peak_energies[j] - my_offset) / my_gain;
				// std::cout << "At j = " << j <<  "At i = " << i<< ", xpeaks_full[j] = " << xpeaks_full[j] << std::endl;
				// if (j != 0) { // Figures out approximately where the peak is located
				// 	gain_estimate = xpeaks_full[j-1]/peak_energies[j-1];
				// 	charge_estimates[j] = gain_estimate*peak_energies[j];
				// }
				// else {
				// 	charge_estimates[j] = charge_estimate_init;
				// }
				

				channel_charge->GetXaxis()->SetRangeUser(charge_estimates[j]-max_fit_radius, charge_estimates[j]+max_fit_radius);
				TSpectrum* peak_search = new TSpectrum();
				peak_search->Search(channel_charge, peak_search_sigma, "nodraw", peak_search_threshold); // Cant find info about Search() online. appear Search()
				xpeaks = peak_search->GetPositionX();

				int npeaks_found = peak_search->GetNPeaks();

				if(npeaks_found <= 0) {
					str << "No peaks found for Detector " << i
						<< " Strip " << strip+1
						<< " Peak " << j+1 << "\n";
					continue;
				}

				peak_estimate = xpeaks[0];
				double best_distance = std::abs(xpeaks[0] - charge_estimates[j]);

				for(int p = 1; p < npeaks_found; p++) {

					// For second peak, only allow peaks higher than expected estimate
					if(j == 1 && xpeaks[p] < charge_estimates[j]) {
						continue;
					}

					double distance = std::abs(xpeaks[p] - charge_estimates[j]);

					if(distance < best_distance) {
						best_distance = distance;
						peak_estimate = xpeaks[p];
					}
				}
				// if (j!=0 && (abs(xpeaks[1]-charge_estimates[j]) < abs(xpeaks[0]-charge_estimates[j]))) { // if j doesn't = 0 AND the other condition is meet then .... 
				// 	peak_estimate = xpeaks[1];
				// 	//cout << " Crystal: "<< i+1 << " triggered"; //doesnt seem like this is printing anywhere so is this condition never being met??
				// }
				// else if (j == 0 && (xpeaks[1] < xpeaks[0] && xpeaks[1] > charge_estimates[j]-max_fit_radius)) {
				// 	peak_estimate = xpeaks[1];
				// 	//cout << "this is working";
				// }
				// else {
				// 	peak_estimate = xpeaks[0];
				// 	// cout << "test2"; //only this one is working 
				// }
				
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
						TF1 *peak_fitter = new TF1("peak_fitter", "gaus(0)+pol1(3)", min_search, max_search);
						//TFitResultPtr fit_result = peak_fitter->Fit(channel_charge, "RQ0S"); // R for fit range, Q for quiet mode, 0 to not store fit results, S to return fit result pointer
						TFitResultPtr fit_result = channel_charge->Fit(peak_fitter, "RQ0S");
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
					// std::cout << "fit worked!" << std::endl;
				}
				else {
					xpeaks_full[j] = peak_estimate;
					str << peak_estimate << "\t";
					// std::cout << "fit didnt work!" << std::endl;
				}
			}
			
			// Perform linear fit to the peaks
			TGraph* fit_graph = new TGraph(NPEAKS, xpeaks_full, peak_energies);
			TF1* linear_fit = new TF1("linear_fit", "[0]+[1]*x", xpeaks_full[0], xpeaks_full[NPEAKS-1]);
			fit_graph->Fit(linear_fit, "qr");
			double offset = linear_fit->GetParameter(0), e_offset = linear_fit->GetParError(0);
			double gain = linear_fit->GetParameter(1), e_gain = linear_fit->GetParError(1);
			double chi_square = linear_fit->GetChisquare();
			chiSquareMap[i][strip] = chi_square;

			// Write calibration coefficients
			std::string channel_name = Form("RCS%02dXP%02dX", i, strip); // P = Front, N = Back
			std::cout << "Trying to update: " << channel_name << std::endl;
			TChannel* detector_channel = TChannel::FindChannelByName(channel_name.c_str());
			if (!detector_channel) {
					std::cerr << "Channel not found: " << channel_name << std::endl;
					continue; // skip to the next strip
				}
				std::cout << channel_name << "  new offset = " << offset<< "  new gain = " << gain << std::endl;
			detector_channel->DestroyENGCal();
			detector_channel->AddENGCoefficient(offset);
			detector_channel->AddENGCoefficient(gain);
			str << "\n\nGain is: " << gain << " ± " << e_gain << std::endl <<
			"Offset is: " << offset << " ± " << e_offset << std::endl <<
			"Chi2 is: " << chi_square << "\n" << std::endl <<
			"~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~" << std::endl;
			std::cout << "UPDATED " << channel_name << std::endl;
		}

	// Write new .cal file
	std::string calibration_file = RCMPOptions::calibration_directory_files + "energy_calibration" + std::to_string(run_number) + ".cal";
	TChannel::WriteCalFile(calibration_file);
	std::cout << "Wrote calibration file: " << calibration_file << std::endl;

	std::cout << "Run Being Calibrated: " << run_number << std::endl;

	// Plot reduced chi squared values for each crystal's fit
	double strip_num[NStrips];
	std::iota(strip_num, strip_num + NStrips, 1);
	TCanvas* c0 = new TCanvas(Form("c_chi_D%d", i),Form("c_chi_D%d", i),1000,400);
	TGraph* chi_square_graph = new TGraph(NStrips, strip_num, chiSquareMap[i].data());
	chi_square_graph->SetTitle(Form("Detector %d Chi Square;Strip Number;Chi^{2}", i));
	chi_square_graph->Draw("ap*");

	// Output statistics
	// int chi_size = sizeof(chisquares) / sizeof(chisquares[0]); 
	// std::vector<double> vectorchisquares(chisquares, chisquares+chi_size);
	mean_chi_square = TMath::Mean(chiSquareMap[i].begin(),chiSquareMap[i].end());
	str << "\nMean chi square: " << mean_chi_square << std::endl;
	if (do_fit) {
		double success_rate = fit_success_counter/(1.0*NPEAKS*NStrips);
		str << "\nTotal fit success rate = " << 100*success_rate << "%" << std::endl;
		double success_rate_peaks[NPEAKS];
		for (int i = 0; i < NPEAKS; i++) {
			success_rate_peaks[i] = fit_success_peaks_counter[i]/(1.0*NStrips);
			str << "Peak " << i+1 << " success rate = " << 100*success_rate_peaks[i] << "%" << std::endl;
		}
	}

	}

	
	// std::string calibration_file = RCMPOptions::calibration_directory_files + "energy_calibration" + std::to_string(run_number) + ".cal";
	// TChannel::WriteCalFile(calibration_file);
	// std::cout << "Wrote calibration file: " << calibration_file << std::endl;
	my_file->Close();
	std::cout << str.str();
    tout_log << str.str();
	tout_log.close();

	if (auto_quit) {
		gROOT->ProcessLine(".q");
	}
	//gROOT -> ProcessLine(".q");
	
}