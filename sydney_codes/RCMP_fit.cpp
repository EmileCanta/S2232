

#include <iostream>
#include <string>
#include <TFile.h>
#include <TH2F.h>
#include <TH1F.h>
#include <TCanvas.h>
#include <TF1.h>
#include <TLegend.h>
#include <TGraphErrors.h>
#include <TLinearFitter.h>
#include <TCanvas.h>
#include <vector>
#include <TROOT.h>
#include "/home/sydney/GRSISort/include/RCMPOptions .h"

/*   
To do:
- Move the opening of the histogram to the top of the code !!done!!
- Get a way to store all of the centroids (change the void to something else?) !!done!!
- Run the linear fit from the centroids !!done!!
- Get the code to then open the cal file and update the gain and offset !!done!!
- Try to get the peak finder to find the rough location of the first peak and then set 
  that as the centre. !!nah!!
- Make sure this is working correctly for the Charge histograms!!! !!done!!

*/


/* 
Example for returining value in function:
double first_peak_fit(int detNum, int stripNum, int centre) {
    // ... do the fit ...
    double centroid = fDouble->GetParameter(1);
    return centroid;
}

*/
std::cout << "To do inital fitting: fit_all()" << std::endl;
std::cout << "To do calibration: calibration()" << std::endl;

// Store Global Results // 
struct FitResult {
    double doubleCentroid1, doubleError1;
    double doubleCentroid2, doubleError2;
    double doubleCentroid3, doubleError3;

    double singleCentroid1, singleError1;
    double singleCentroid2, singleError2;
    double singleCentroid3, singleError3;

    int detNum;
    int stripNum;
};

// Global variable to store the centroids for a given strip
FitResult storedCentroids;
bool centroidsStored = false;  // flag to know if centroids exist


// Open the histogram!! //
//static TFile* file = TFile::Open("/home/sydney/S2232_analysis/Histograms/RCMP/RCMP28226_000-014.root");
static TFile* file = TFile::Open("/home/sydney/S2232_analysis/Histograms/run28226/RCMP28226_009-014.root");

TH1F* GetStripHist(int detNum, int stripNum, const std::string& tag) {
    std::string histName = "Det" + std::to_string(detNum) + "FECharge";
    TH2F* hist = (TH2F*)file->Get(histName.c_str());
    if (!hist) {
        std::cerr << "Error: histogram " << histName << " not found in file!" << std::endl;
        return nullptr;
    }
    // Use unique projection name by appending tag
    return (TH1F*)hist->ProjectionY((histName + "_" + tag).c_str(), stripNum, stripNum);
}

// ===== Double and Triple Gaussians ===== //

// --- Function: Fit the first peak with two Gaussians (main + shoulder) ---
std::pair<double,double> first_peak_fit(int detNum, int stripNum, int centre, bool quiet = false) {
    static TCanvas* cFit_first = new TCanvas("cFit_first", "First Peak Fit", 800, 600);
    cFit_first->cd();  // make this canvas active

    TH1F* h1 = GetStripHist(detNum, stripNum, "first_dbl");

    double lowerBound = centre -115; //initally 5050
    double upperBound = centre +115; // initally 5230


    // --- Combined double Gaussian function: gaus(0)+gaus(3) ---
    static TF1 *fDouble_first = new TF1("fDouble_first", "gaus(0)+gaus(3)", lowerBound, upperBound);

    // Initial parameter guesses
    // [0] amplitude1, [1] mean1, [2] sigma1
    // [3] amplitude2, [4] mean2, [5] sigma2
    fDouble_first->SetParameters(
        5500, centre, 13,   // main peak
        1500,  centre - 40, 14    // shoulder 40
    );
    fDouble_first->SetParNames(
        "Peak 1 Amp", "Peak 1 Centroid", "Peak 1 Sigma",
        "Shoulder 1 Amp", "Shoulder 1 Centroid", "Shoulder 1 Sigma"
    );

    // Choose fit options based on quiet flag
    std::string fitOptions = quiet ? "RQ0+" : "R+";
    

    // Set colors and style
    fDouble_first->SetLineColor(kBlue);
    fDouble_first->SetLineWidth(2);

    // Fit in the range
    h1 -> Fit(fDouble_first, fitOptions.c_str());
    //h1->Fit(fDouble_first, "R+");

    // Get Centroid
    double centroid_1 = fDouble_first -> GetParameter(1);
    double error_1 = fDouble_first -> GetParError(1);

    // Expand view range ±50 beyond fit
    h1->GetXaxis()->SetRangeUser(lowerBound - 50, upperBound + 50);

    // Draw histogram
    h1->Draw("same");

    // To visualize the individual Gaussians, define them separately but don't fit
    static TF1 *fMain_first = new TF1("fMain_first", "gaus", lowerBound, upperBound);
    fMain_first->SetParameters(fDouble_first->GetParameter(0),
                         fDouble_first->GetParameter(1),
                         fDouble_first->GetParameter(2));
    fMain_first->SetLineColor(kRed);
    fMain_first->SetLineStyle(2);
    fMain_first->Draw("same");
    

    static TF1 *fShoulder_first = new TF1("fShoulder_first", "gaus", lowerBound, upperBound);
    fShoulder_first->SetParameters(fDouble_first->GetParameter(3),
                             fDouble_first->GetParameter(4),
                             fDouble_first->GetParameter(5));
    fShoulder_first->SetLineColor(kMagenta);
    fShoulder_first->SetLineStyle(2);
    fShoulder_first->Draw("same");

    // Legend
    // auto legend = new TLegend(0.6,0.7,0.9,0.9);
    // legend->AddEntry(fDouble_first,"Total Fit (Main+Shoulder)","l");
    // legend->AddEntry(fMain_first,"Main Peak","l");
    // legend->AddEntry(fShoulder_first,"Shoulder","l");
    // legend->Draw();

    return {centroid_1, error_1};
}

// --- Function: Fit the second peak with three Gaussians (main + 2 shoulders) ---
std::pair<double,double> second_peak_fit(int detNum, int stripNum, int centre, bool quiet = false) {
    static TCanvas* cFit = new TCanvas("cFit_second", "Second Peak Fit", 800, 600);
    cFit->cd();  // make this canvas active

    TH1F* h1 = GetStripHist(detNum, stripNum, "second_tbl");

    double new_centre = centre + 240;
    double lowerBound = new_centre - 140;   
    double upperBound = new_centre + 100; 

    // --- Combined triple Gaussian function: gaus(0)+gaus(3)+gaus(6) ---
    static TF1 *fTriple = new TF1("fTriple_second", "gaus(0)+gaus(3)+gaus(6)", lowerBound, upperBound);

    // Initial parameter guesses
    // [0] amplitude1, [1] mean1, [2] sigma1
    // [3] amplitude2, [4] mean2, [5] sigma2
    // [6] amplitude3, [7] mean3, [8] sigma3
    fTriple->SetParameters(
        5000, new_centre, 9,   // main peak 250 (D6), 200 (D3), 250
        4000,  new_centre - 30, 12,   // first shoulder 205 (D6), 160 (D3)
        200,   new_centre - 65, 12    // smallest shoulder 175 (D6), 145 (D3)
    );
    fTriple->SetParNames(
        "Peak 2 Amp", "Peak 2 Centroid", "Peak 2 Sigma",
        "Shoulder 2 Amp", "Shoulder 2 Centroid", "Shoulder 2 Sigma",
        "Small Shoulder 2 Amp", "Small Shoulder 2 Centroid", "Small Shoulder 2 Sigma"
    );
    fTriple->SetParLimits(2, 1, 50);  // sigma1 between 1 and 50
    fTriple->SetParLimits(5, 1, 50);  // sigma2 between 1 and 50
    fTriple->SetParLimits(8, 1, 50);

    // Choose fit options based on quiet flag
    std::string fitOptions = quiet ? "RQ0+" : "R+";

    // Set colors and style
    fTriple->SetLineColor(kBlue);
    fTriple->SetLineWidth(2);

    // Fit in the range
    h1 -> Fit(fTriple, fitOptions.c_str());
    //h1->Fit(fTriple, "R+");

    // Get Centroid
    double centroid_triple = fTriple -> GetParameter(1);
    double error_triple = fTriple -> GetParError(1);

    // Expand view range ±50 beyond fit
    h1->GetXaxis()->SetRangeUser(lowerBound - 50, upperBound + 50);

    // Draw histogram
    h1->Draw();

    // --- Visualize individual Gaussians ---
    static TF1 *fMain = new TF1("fMain_second", "gaus", lowerBound, upperBound);
    fMain->SetParameters(fTriple->GetParameter(0),
                         fTriple->GetParameter(1),
                         fTriple->GetParameter(2));
    fMain->SetLineColor(kRed);
    fMain->SetLineStyle(2);
    fMain->Draw("same");

    static TF1 *fShoulder1 = new TF1("fShoulder1_second", "gaus", lowerBound, upperBound);
    fShoulder1->SetParameters(fTriple->GetParameter(3),
                              fTriple->GetParameter(4),
                              fTriple->GetParameter(5));
    fShoulder1->SetLineColor(kMagenta);
    fShoulder1->SetLineStyle(2);
    fShoulder1->Draw("same");

    static TF1 *fShoulder2 = new TF1("fShoulder2_second", "gaus", lowerBound, upperBound);
    fShoulder2->SetParameters(fTriple->GetParameter(6),
                              fTriple->GetParameter(7),
                              fTriple->GetParameter(8));
    fShoulder2->SetLineColor(kGreen+2);
    fShoulder2->SetLineStyle(2);
    fShoulder2->Draw("same");

    // Legend
    // auto legend = new TLegend(0.6,0.7,0.9,0.9);
    // legend->AddEntry(fTriple,"Total Fit (Main+2 Shoulders)","l");
    // legend->AddEntry(fMain,"Main Peak (5497)","l");
    // legend->AddEntry(fShoulder1,"Shoulder (5450)","l");
    // legend->AddEntry(fShoulder2,"Small Shoulder (5403)","l");
    // legend->Draw();

    return {centroid_triple, error_triple};
}

// --- Function: Fit the second peak with two Gaussians (main + 1 shoulder) ---
void second_peak_simple_fit(int detNum, int stripNum, int centre) {
    static TCanvas* cFit = new TCanvas("cFit_second_simple", "Second Peak (Main + Shoulder)", 800, 600);
    cFit->cd();  // make this canvas active

    // Open file and get histogram
    TH1F* h1 = GetStripHist(detNum, stripNum, "second_dbl");
    
    double lowerBound = centre + 255; //5420
    double upperBound = centre + 385; // 5550

    // --- Combined double Gaussian function: gaus(0)+gaus(3) ---
    static TF1 *fDouble = new TF1("fDouble_second", "gaus(0)+gaus(3)", lowerBound, upperBound);

    // Initial parameter guesses
    // [0] amplitude1, [1] mean1, [2] sigma1
    // [3] amplitude2, [4] mean2, [5] sigma2
    fDouble->SetParameters(
        1398, centre + 330, 12,   // main peak
        263,  centre + 290, 16    // shoulder
    );
    
     fDouble->SetParNames(
        "Peak 2 Amp", "Peak 2 Centroid", "Peak 2 Sigma",
        "Shoulder 2 Amp", "Shoulder 2 Centroid", "Shoulder 2 Sigma"
    );

    // Set colors and style
    fDouble->SetLineColor(kBlue);
    fDouble->SetLineWidth(2);

    // Fit in the range
    h1->Fit(fDouble, "R+");

    // Get Centroid
    double centroid_2 = fDouble -> GetParameter(1);
    double error_2 = fDouble -> GetParError(1);

    // Expand view range ±50 beyond fit
    h1->GetXaxis()->SetRangeUser(lowerBound - 50, upperBound + 50);

    // Draw histogram
    h1->Draw();

    // --- Visualize individual Gaussians ---
    static TF1 *fMain = new TF1("fMain_second_simple", "gaus", lowerBound, upperBound);
    fMain->SetParameters(fDouble->GetParameter(0),
                         fDouble->GetParameter(1),
                         fDouble->GetParameter(2));
    fMain->SetLineColor(kRed);
    fMain->SetLineStyle(2);
    fMain->Draw("same");

    static TF1 *fShoulder = new TF1("fShoulder_second_simple", "gaus", lowerBound, upperBound);
    fShoulder->SetParameters(fDouble->GetParameter(3),
                             fDouble->GetParameter(4),
                             fDouble->GetParameter(5));
    fShoulder->SetLineColor(kMagenta);
    fShoulder->SetLineStyle(2);
    fShoulder->Draw("same");

    // Legend
    auto legend = new TLegend(0.6,0.7,0.9,0.9);
    legend->AddEntry(fDouble,"Total Fit (Main+Shoulder)","l");
    legend->AddEntry(fMain,"Main Peak (5497)","l");
    legend->AddEntry(fShoulder,"Shoulder (5450)","l");
    legend->Draw();
}

// --- Function: Fit the third peak with two Gaussians (main + shoulder) ---
std::pair<double,double> third_peak_fit(int detNum, int stripNum, int centre, bool quiet = false) {
    static TCanvas* cFit = new TCanvas("cFit_third", "Third Peak Fit", 800, 600);
    cFit->cd();  // make this canvas active

    // Open file and get histogram
    TH1F* h1 = GetStripHist(detNum, stripNum, "third_dbl");

    double new_centre = centre + 470;
    double lowerBound = new_centre - 120; // 380 (D6), 300 (D3)
    double upperBound = new_centre + 80; // 555 (D6), 460 (D3)

    // --- Combined double Gaussian function: gaus(0)+gaus(3) ---
    static TF1 *fDouble = new TF1("fDouble_third", "gaus(0)+gaus(3)", lowerBound, upperBound);

    // Initial parameter guesses
    // [0] amplitude1, [1] mean1, [2] sigma1
    // [3] amplitude2, [4] mean2, [5] sigma2
    fDouble->SetParameters(
        4000, new_centre, 10,   // main peak, 460 (D6), 485 (others), 415 (D3)
        1500, new_centre - 30, 12    // shoulder, 435 (D6), 455 (others), 390 (D3)
    );
     fDouble->SetParNames(
        "Peak 3 Amp", "Peak 3 Centroid", "Peak 3 Sigma",
        "Shoulder 3 Amp", "Shoulder 3 Centroid", "Shoulder 3 Sigma"
    );

    fDouble->SetParLimits(2, 1, 50);  // sigma1 between 1 and 50
    fDouble->SetParLimits(5, 1, 50);

    // Choose fit options based on quiet flag
    std::string fitOptions = quiet ? "RQ0+" : "R+";

    // Set colors and style
    fDouble->SetLineColor(kBlue);
    fDouble->SetLineWidth(2);

    // Fit in the range
    h1 -> Fit(fDouble, fitOptions.c_str());
    //h1->Fit(fDouble, "R+");

    // Get Centroid
    double centroid_3 = fDouble -> GetParameter(1);
    double error_3 = fDouble -> GetParError(1);

    // Expand view range ±50 beyond fit
    h1->GetXaxis()->SetRangeUser(lowerBound - 50, upperBound + 50);

    // Draw histogram
    h1->Draw();

    // To visualize the individual Gaussians, define them separately but don't fit
    static TF1 *fMain = new TF1("fMain_third", "gaus", lowerBound, upperBound);
    fMain->SetParameters(fDouble->GetParameter(0),
                         fDouble->GetParameter(1),
                         fDouble->GetParameter(2));
    fMain->SetLineColor(kRed);
    fMain->SetLineStyle(2);
    fMain->Draw("same");

    static TF1 *fShoulder = new TF1("fShoulder_third", "gaus", lowerBound, upperBound);
    fShoulder->SetParameters(fDouble->GetParameter(3),
                             fDouble->GetParameter(4),
                             fDouble->GetParameter(5));
    fShoulder->SetLineColor(kMagenta);
    fShoulder->SetLineStyle(2);
    fShoulder->Draw("same");

    // Legend
    // auto legend = new TLegend(0.6,0.7,0.9,0.9);
    // legend->AddEntry(fDouble,"Total Fit (Main+Shoulder)","l");
    // legend->AddEntry(fMain,"Main Peak","l");
    // legend->AddEntry(fShoulder,"Shoulder","l");
    // legend->Draw();

    return {centroid_3, error_3};
}


// ===== Single Gaussian ===== //

std::pair<double,double>  single_first_peak(int detNum, int stripNum, int centre, bool quiet = false) {
    static TCanvas* cFit = new TCanvas("cSingle_Fit_first", "First Peak Fit", 800, 600);
    cFit->cd();  // make this canvas active

    TH1F* h1 = GetStripHist(detNum, stripNum, "first");

    double lowerBound = centre -115; //initally 5050
    double upperBound = centre +115; // initally 5230

    // --- Combined double Gaussian function: gaus(0)+gaus(3) ---
    TF1 *fSingle = new TF1("fSingle", "gaus(0)", lowerBound, upperBound);

    // Initial parameter guesses
    // [0] amplitude1, [1] mean1, [2] sigma1
    fSingle->SetParameters(
        1900, centre, 13   // main peak
    );
    fSingle->SetParNames(
        "Peak 1 Amp", "Peak 1 Centroid", "Peak 1 Sigma"
    );

    // Choose fit options based on quiet flag
    std::string fitOptions = quiet ? "RQ0+" : "R+";

    // Set colors and style
    fSingle->SetLineColor(kBlue);
    fSingle->SetLineWidth(2);

    // Fit in the range
    h1 -> Fit(fSingle, fitOptions.c_str());
    //h1->Fit(fSingle, "R+");

    // Get Centroid
    double centroid_1S = fSingle -> GetParameter(1);
    double error_1S = fSingle -> GetParError(1);

    // Expand view range ±50 beyond fit
    h1->GetXaxis()->SetRangeUser(lowerBound - 50, upperBound + 50);

    // Draw histogram
    h1->Draw("same");



    return {centroid_1S, error_1S};
}

std::pair<double,double>  single_second_peak(int detNum, int stripNum, int centre, bool quiet = false) {
    TCanvas* cFit = new TCanvas("cSingle_Fit_second", "Second Peak Fit", 800, 600);
    cFit->cd();

    TH1F* h1 = GetStripHist(detNum, stripNum, "second");

    double new_centre = centre + 220;
    double lowerBound = new_centre - 140;   
    double upperBound = new_centre + 100;

    // Single Gaussian
    TF1 *fSingle = new TF1("fSingle_second", "gaus(0)", lowerBound, upperBound);
    fSingle->SetParameters(1400, new_centre, 12); // amplitude, mean, sigma
    fSingle->SetParNames("Peak 2 Amp", "Peak 2 Centroid", "Peak 2 Sigma");
    fSingle->SetLineColor(kBlue);
    fSingle->SetLineWidth(2);

    // Choose fit options based on quiet flag
    std::string fitOptions = quiet ? "RQ0+" : "R+";

    // Fit
    h1 -> Fit(fSingle, fitOptions.c_str());
    //h1->Fit(fSingle, "R+");

    // Get Centroid
    double centroid_2S = fSingle -> GetParameter(1);
    double error_2S = fSingle -> GetParError(1);

    // Expand view range ±50 beyond fit
    h1->GetXaxis()->SetRangeUser(lowerBound - 50, upperBound + 50);

    // Draw histogram and fit
    h1->Draw();
    fSingle->Draw("same");

    return {centroid_2S, error_2S};
}

std::pair<double,double>  single_third_peak(int detNum, int stripNum, int centre, bool quiet = false) {
    TCanvas* cFit = new TCanvas("cSingle_Fit_third", "Third Peak Fit", 800, 600);
    cFit->cd();

    TH1F* h1 = GetStripHist(detNum, stripNum, "third");

    double new_centre = centre + 430;
    double lowerBound = new_centre - 80; // 380 (D6), 300 (D3)
    double upperBound = new_centre + 80; // 555 (D6), 460 (D3)

    // Single Gaussian
    TF1 *fSingle = new TF1("fSingle_third", "gaus(0)", lowerBound, upperBound);
    fSingle->SetParameters(1080, new_centre, 11); // amplitude, mean, sigma
    fSingle->SetParNames("Peak 3 Amp", "Peak 3 Centroid", "Peak 3 Sigma");
    fSingle->SetLineColor(kBlue);
    fSingle->SetLineWidth(2);

    // Choose fit options based on quiet flag
    std::string fitOptions = quiet ? "RQ0+" : "R+";

    // Fit
    h1 -> Fit(fSingle, fitOptions.c_str());
    //h1->Fit(fSingle, "R+");

    // Get Centroid
    double centroid_3S = fSingle -> GetParameter(1);
    double error_3S = fSingle -> GetParError(1);

    // Expand view range ±50 beyond fit
    h1->GetXaxis()->SetRangeUser(lowerBound - 50, upperBound + 50);

    // Draw histogram and fit
    h1->Draw();
    fSingle->Draw("same");

    
    return {centroid_3S, error_3S};
}




void fit_all() {
    int detNum;
    int stripNum;
    int centre;
    std::cout << "Detector: ";
    std::cin >> detNum;
    std::cout << "Strip: ";
    std::cin >> stripNum;
    std::cout << "Centre of first peak: ";
    std::cin >> centre;

    // std::vector<std::pair<double,double>> double_fits;
    // double_fits.push_back(first_peak_fit(detNum, stripNum, centre));
    // double_fits.push_back(second_peak_fit(detNum, stripNum, centre));
    // double_fits.push_back(third_peak_fit(detNum, stripNum, centre));

    // --- Run double Gaussian fits ---
    auto first_double = first_peak_fit(detNum, stripNum, centre);
    auto second_double = second_peak_fit(detNum, stripNum, centre);
    auto third_double = third_peak_fit(detNum, stripNum, centre);

    // --- Run single Gaussian fits ---
    auto first_single = single_first_peak(detNum, stripNum, centre);
    auto second_single = single_second_peak(detNum, stripNum, centre);
    auto third_single = single_third_peak(detNum, stripNum, centre);

    // --- Store results in global struct ---
    storedCentroids.doubleCentroid1 = first_double.first;
    storedCentroids.doubleError1    = first_double.second;
    storedCentroids.doubleCentroid2 = second_double.first;
    storedCentroids.doubleError2    = second_double.second;
    storedCentroids.doubleCentroid3 = third_double.first;
    storedCentroids.doubleError3    = third_double.second;

    storedCentroids.singleCentroid1 = first_single.first;
    storedCentroids.singleError1    = first_single.second;
    storedCentroids.singleCentroid2 = second_single.first;
    storedCentroids.singleError2    = second_single.second;
    storedCentroids.singleCentroid3 = third_single.first;
    storedCentroids.singleError3    = third_single.second;

    storedCentroids.detNum = detNum;
    storedCentroids.stripNum = stripNum;

    centroidsStored = true; // mark that we have stored them

    // for(size_t i=0; i<double_fits.size(); i++) {
    //     std::cout << "Double fit " << i+1 
    //               << " = " << double_fits[i].first 
    //               << " ± " << double_fits[i].second << std::endl;
    // }
}




void calibration() {

    // Check if centroids were stored
    if (!centroidsStored) { 
        std::cerr << "Error: Please run fit_all() before calibration.\n";
        return;
    }

    // Turn off canvas
    bool oldBatch = gROOT ->IsBatch();
    gROOT -> SetBatch(kTRUE);

    int detNum = storedCentroids.detNum;
    int stripNum = storedCentroids.stripNum;
    // std::cout << "Detector: ";
    // std::cin >> detNum;
    // std::cout << "Strip: ";
    // std::cin >> stripNum;

    std::vector<std::pair<double,double>> reference = {
    {5156.59, 0.31},   // energy ± uncertainty
    {5485.56, 0.12},    // CHECK THE UNCERTAINTIES!!
    {5804.77, 0.06}
    };

    // Ask which fit to use
    std::cout << "Fit with double gaussian? (y/n): ";
    std::string answer;
    std::cin >> answer;

    // Select stored centroids
    std::vector<std::pair<double,double>> chosen;
    if(answer == "y" || answer == "Y") {
        chosen = {
            {storedCentroids.doubleCentroid1, storedCentroids.doubleError1},
            {storedCentroids.doubleCentroid2, storedCentroids.doubleError2},
            {storedCentroids.doubleCentroid3, storedCentroids.doubleError3}
        };
    } else {
        chosen = {
            {storedCentroids.singleCentroid1, storedCentroids.singleError1},
            {storedCentroids.singleCentroid2, storedCentroids.singleError2},
            {storedCentroids.singleCentroid3, storedCentroids.singleError3}
        };
    }

    // Add more uncertianty for bad fits
    // std::cout << "Add extra uncertainty to centroids? (y/n): ";
    // std::string addUnc;
    // std::cin >> addUnc;

    // double extra = 0.0;
    // if(addUnc == "y" || addUnc == "Y") {
    //     std::cout << "Enter extra uncertainty to add: ";
    //     std::cin >> extra;
    // }

    // for(auto& p : chosen) {
    //     p.second += extra;
    // }

    // Print out the chosen list
    for(const auto& p : chosen) {
    std::cout << p.first << " ± " << p.second << std::endl;
    }

    // Linear Fit //
    int NPEAKS = chosen.size();

    std::vector<double> xpeaks(NPEAKS);
    std::vector<double> xerrors(NPEAKS);
    std::vector<double> yenergies(NPEAKS);
    std::vector<double> yerrors(NPEAKS);

    for (int i = 0; i< NPEAKS; i++) {

        xpeaks[i] = chosen[i].first;
        xerrors[i] = chosen[i].second;
        yenergies[i] = reference[i].first;
        yerrors[i] = reference[i].second;

    }

    gROOT -> SetBatch(kFALSE);
    TCanvas* c1 = new TCanvas("c1", "Calibration Fit", 800, 600);

    // To take into account the errors of the x-axis and not the y-axis invert the 
    TGraphErrors* fit_graph = new TGraphErrors(NPEAKS, xpeaks.data(), yenergies.data(), xerrors.data(), yerrors.data());
    fit_graph -> SetTitle("Calibration Fit;Channel;Energy (keV)");
    fit_graph -> SetMarkerStyle(20);
    fit_graph -> SetMarkerSize(1.2);

    TF1* linear_fit = new TF1("linear_fit", "[0]+[1]*x", xpeaks.front(), xpeaks.back());
    fit_graph ->Fit(linear_fit, "Q");

    fit_graph -> Draw("AP"); // A = axes, P = points
    linear_fit -> Draw("same");
    c1 -> Update();

    double offset = linear_fit -> GetParameter(0);
    double e_offset = linear_fit -> GetParError(0);
    double gain = linear_fit -> GetParameter(1);
    double e_gain = linear_fit -> GetParError(1);
    double chi_square = linear_fit->GetChisquare();

    std::cout << "\n\nGain is: " << gain << " ± " << e_gain << std::endl <<
			"Offset is: " << offset << " ± " << e_offset << std::endl <<
			"Chi2 is: " << chi_square << "\n" << std::endl;

   

    
    // Comparing the results //
    std::cout << "\nCalibrated energies vs reference:\n";
    for (int i = 0; i < NPEAKS; i++) {
        double predictedEnergy = offset + gain * xpeaks[i];
        double residual = predictedEnergy - yenergies[i];
        std::cout << " Peak " << i+1 << ": predicted = " << predictedEnergy << " keV, reference = " << yenergies[i] << " keV, residual = " << residual << " keV" << std::endl;
    }


    // Writing Cal File //

    std::string channel_name = Form("RCS%02dXP%02dX", detNum, stripNum-1); // P = Front, N = Back

    std::string calibration_file = RCMPOptions::calibration_directory_files + "energy_calibration28226.cal";
	TChannel::ReadCalFile(calibration_file.c_str());

    TChannel* detector_channel = TChannel::FindChannelByName(channel_name.c_str());
    if (!detector_channel) 
        {
		std::cerr << "Channel not found: " << channel_name << std::endl;
		}
    else 
        {
        detector_channel->DestroyENGCal();
		detector_channel->AddENGCoefficient(offset);
		detector_channel->AddENGCoefficient(gain);

        //std::string calibration_file = RCMPOptions::calibration_directory_files + "energy_calibration28226.cal";
	    TChannel::WriteCalFile(calibration_file);

        std::cout << "Calibration updated for " << channel_name << " and written to: " << calibration_file << std::endl;

        }

}