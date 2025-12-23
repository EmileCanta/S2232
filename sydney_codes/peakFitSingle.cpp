#include <string>
#include <iostream>
#include <filesystem>
#include <cmath>
#include <TFile.h>
#include <TH2.h>
#include <TRWPeak.h>
#include "../include/FindHistogramFile.h"
#include "../include/GAOptions.h"

void peakFitSingle(int run_number, std::string options = ""){

    bool auto_quit = options.find("q") != std::string::npos;

    // Parameters to be set by the user
    std::string histogram_directory = GAOptions::histogram_directory;
    std::string histogram_prefix = "Alphanumeric";
    double peak_energy = 238.27;
    double fit_radius = 10.0;

    // Find and open histogram file
    std::string file_name = findHistogramFile(run_number, histogram_directory, histogram_prefix);
    auto my_file = new TFile(file_name.c_str());

    // Get 1D gamma singles histogram
    TH2D* energy_2d = (TH2D*)my_file->Get("EnergyVsChannel");
    TH1D* energy_1d = (TH1D*)energy_2d->ProjectionY("py", 0, -1);

    // Fit peak
    TPeakFitter* peak_fitter = new TPeakFitter(peak_energy-fit_radius, peak_energy+fit_radius);
    TRWPeak* peak1 = new TRWPeak(peak_energy);
    peak_fitter->AddPeak(peak1);
    TFitResultPtr fit_result = peak_fitter->Fit(energy_1d, "q");
    double centroid = fit_result->Parameter(1);
    double std_dev = fit_result->Parameter(2);

    std::cout << "Centroid: " << centroid << " keV" << std::endl;
    std::cout << "Standard deviation: " << std_dev << " keV" << std::endl;
    std::cout << "FWHM: " << 2*sqrt(2*log(2))*std_dev << " keV" << std::endl;

    if (auto_quit) {
        gROOT->ProcessLine(".q");
    }
}