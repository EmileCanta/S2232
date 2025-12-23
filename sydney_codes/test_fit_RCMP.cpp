

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
#include <vector>
#include <numeric>
#include <TROOT.h>
#include <TGraph.h>
#include "/home/sydney/GRSISort/include/RCMPOptions .h"




// Open the histogram!! //
static TFile* file = TFile::Open("/home/sydney/S2232_analysis/Histograms/RCMP/RCMP28226_000-014.root");
//static TFile* file = TFile::Open("/home/sydney/S2232_analysis/Histograms/run28226/RCMP28226_009-014.root");

TH1F* GetStripHist(int detNum, int firstBin, int lastBin, const std::string& tag, const std::string& plane) {
    std::string histName = "Det" + std::to_string(detNum) + plane;
    TH2F* hist = (TH2F*)file->Get(histName.c_str());
    if (!hist) {
        std::cerr << "Error: histogram " << histName << " not found in file!" << std::endl;
        return nullptr;
    }
    // Use unique projection name by appending tag
    return (TH1F*)hist->ProjectionY((histName + "_" + tag).c_str(), firstBin, lastBin);
}

// ===== Double and Triple Gaussians ===== //

// --- Function: Fit the first peak with two Gaussians (main + shoulder) ---
double first_peak_fit(int detNum, int firstXbin, int lastXbin, bool quiet, const std::string& plane) {
    
    TCanvas* cFit_first = new TCanvas(Form("cFit_first_%d", firstXbin), "First Peak Fit", 800, 600);
    cFit_first->cd();  // make this canvas active
    TH1F* h1 = GetStripHist(detNum, firstXbin, lastXbin, Form("first_dbl_%d", firstXbin), plane);

    double lowerBound = 5040; 
    double upperBound = 5270; 

    // --- Combined double Gaussian function: gaus(0)+gaus(3) ---
    TF1 *fDouble_first = new TF1(Form("fDouble_first_%d", firstXbin), "gaus(0)+gaus(3)", lowerBound, upperBound);

    // Initial parameter guesses
    // [0] amplitude1, [1] mean1, [2] sigma1
    // [3] amplitude2, [4] mean2, [5] sigma2
    fDouble_first->SetParameters(
        31000, 5156, 13,   // main peak
        7000,  5105, 14    // shoulder 40
    );
    fDouble_first->SetParNames(
        "Peak 1 Amp", "Peak 1 Centroid", "Peak 1 Sigma",
        "Shoulder 1 Amp", "Shoulder 1 Centroid", "Shoulder 1 Sigma"
    );

    // Set colors and style
    fDouble_first->SetLineColor(kBlue);
    fDouble_first->SetLineWidth(2);

    // Fit in the range
    std::string fitOptions = quiet ? "RQ0" : "R+";
    h1 -> Fit(fDouble_first, fitOptions.c_str());
    //h1->Fit(fDouble_first, "R+");

    // Get Centroid
    //double centroid_1 = fDouble_first -> GetParameter(1);
    //double error_1 = fDouble_first -> GetParError(1);

    // Get Sigma
    double sigma_1 = fDouble_first -> GetParameter(2);

    h1->GetXaxis()->SetRangeUser(lowerBound - 50, upperBound + 50); // Expand view range ±50 beyond fit

    h1->Draw();

    // To visualize the individual Gaussians, define them separately but don't fit
    TF1 *fMain_first = new TF1(Form("fMain_first_%d", firstXbin), "gaus", lowerBound, upperBound);
    fMain_first->SetParameters(fDouble_first->GetParameter(0),
                         fDouble_first->GetParameter(1),
                         fDouble_first->GetParameter(2));
    fMain_first->SetLineColor(kRed);
    fMain_first->SetLineStyle(2);
    fMain_first->Draw("same");
    

    TF1 *fShoulder_first = new TF1(Form("fShoulder_first_%d", firstXbin), "gaus", lowerBound, upperBound);
    fShoulder_first->SetParameters(fDouble_first->GetParameter(3),
                             fDouble_first->GetParameter(4),
                             fDouble_first->GetParameter(5));
    fShoulder_first->SetLineColor(kMagenta);
    fShoulder_first->SetLineStyle(2);
    fShoulder_first->Draw("same");
    


    return sigma_1;
}

// --- Function: Fit the second peak with three Gaussians (main + 2 shoulders) ---
double second_peak_fit(int detNum,int firstXbin, int lastXbin, bool quiet, const std::string& plane) {
    
    TCanvas* cFit = new TCanvas(Form("cFit_second_%d", firstXbin), "Second Peak Fit", 800, 600);
    cFit->cd();  // make this canvas active
    TH1F* h1 = GetStripHist(detNum, firstXbin, lastXbin, Form("second_tbl_%d", firstXbin), plane);

    double lowerBound = 5340;   
    double upperBound = 5580; 

    // --- Combined triple Gaussian function: gaus(0)+gaus(3)+gaus(6) ---
    TF1 *fTriple = new TF1(Form("fTriple_second_%d", firstXbin), "gaus(0)+gaus(3)+gaus(6)", lowerBound, upperBound);

    // Initial parameter guesses
    // [0] amplitude1, [1] mean1, [2] sigma1
    // [3] amplitude2, [4] mean2, [5] sigma2
    // [6] amplitude3, [7] mean3, [8] sigma3
    fTriple->SetParameters(
        35000, 5485, 9,   // main peak 250 (D6), 200 (D3), 250
        7000,  5442, 12,   // first shoulder 205 (D6), 160 (D3)
        5500,  5388, 12    // smallest shoulder 175 (D6), 145 (D3)
    );
    fTriple->SetParNames(
        "Peak 2 Amp", "Peak 2 Centroid", "Peak 2 Sigma",
        "Shoulder 2 Amp", "Shoulder 2 Centroid", "Shoulder 2 Sigma",
        "Small Shoulder 2 Amp", "Small Shoulder 2 Centroid", "Small Shoulder 2 Sigma"
    );
    fTriple->SetParLimits(2, 1, 50);  // sigma1 between 1 and 50
    fTriple->SetParLimits(5, 1, 50);  // sigma2 between 1 and 50
    fTriple->SetParLimits(8, 1, 50);


    // Set colors and style
    fTriple->SetLineColor(kBlue);
    fTriple->SetLineWidth(2);

    // Fit in the range
    std::string fitOptions = quiet ? "RQ0" : "R+";
    h1 -> Fit(fTriple, fitOptions.c_str());
    //h1->Fit(fTriple, "R+");

    // Get Centroid
    // double centroid_triple = fTriple -> GetParameter(1);
    // double error_triple = fTriple -> GetParError(1);

    // Get Sigma
    double sigma_2 = fTriple -> GetParameter(2);

    h1->GetXaxis()->SetRangeUser(lowerBound - 50, upperBound + 50); // Expand view range ±50 beyond fit
    h1->Draw();


    // --- Visualize individual Gaussians ---
    TF1 *fMain = new TF1(Form("fMain_second_%d", firstXbin), "gaus", lowerBound, upperBound);
    fMain->SetParameters(fTriple->GetParameter(0),
                         fTriple->GetParameter(1),
                         fTriple->GetParameter(2));
    fMain->SetLineColor(kRed);
    fMain->SetLineStyle(2);
    fMain->Draw("same");

    TF1 *fShoulder1 = new TF1(Form("fShoulder1_second_%d", firstXbin), "gaus", lowerBound, upperBound);
    fShoulder1->SetParameters(fTriple->GetParameter(3),
                              fTriple->GetParameter(4),
                              fTriple->GetParameter(5));
    fShoulder1->SetLineColor(kMagenta);
    fShoulder1->SetLineStyle(2);
    fShoulder1->Draw("same");
    

    TF1 *fShoulder2 = new TF1(Form("fShoulder2_second_%d", firstXbin), "gaus", lowerBound, upperBound);
    fShoulder2->SetParameters(fTriple->GetParameter(6),
                              fTriple->GetParameter(7),
                              fTriple->GetParameter(8));
    fShoulder2->SetLineColor(kGreen+2);
    fShoulder2->SetLineStyle(2);
    fShoulder2->Draw("same");

    return sigma_2;
}

// --- Function: Fit the third peak with two Gaussians (main + shoulder) ---
double third_peak_fit(int detNum, int firstXbin, int lastXbin, bool quiet, const std::string& plane) {

    TCanvas* cFit = new TCanvas(Form("cFit_third_%d", firstXbin), "Third Peak Fit", 800, 600);
    cFit->cd();  // make this canvas active
    // Open file and get histogram
    TH1F* h1 = GetStripHist(detNum, firstXbin, lastXbin, Form("third_dbl_%d", firstXbin), plane);

    
    double lowerBound = 5680; 
    double upperBound = 5890; 

    // --- Combined double Gaussian function: gaus(0)+gaus(3) ---
    TF1 *fDouble = new TF1(Form("fDouble_third_%d", firstXbin), "gaus(0)+gaus(3)", lowerBound, upperBound);

    // Initial parameter guesses
    // [0] amplitude1, [1] mean1, [2] sigma1
    // [3] amplitude2, [4] mean2, [5] sigma2
    fDouble->SetParameters(
        60000, 5804, 10,   // main peak, 460 (D6), 485 (others), 415 (D3)
        20000, 5750, 12    // shoulder, 435 (D6), 455 (others), 390 (D3)
    );
     fDouble->SetParNames(
        "Peak 3 Amp", "Peak 3 Centroid", "Peak 3 Sigma",
        "Shoulder 3 Amp", "Shoulder 3 Centroid", "Shoulder 3 Sigma"
    );

    fDouble->SetParLimits(2, 1, 50);  // sigma1 between 1 and 50
    fDouble->SetParLimits(5, 1, 50);


    // Set colors and style
    fDouble->SetLineColor(kBlue);
    fDouble->SetLineWidth(2);

    // Fit in the range
    std::string fitOptions = quiet ? "RQ0" : "R+";
    h1->Fit(fDouble, fitOptions.c_str());

    // Get Centroid
    // double centroid_3 = fDouble -> GetParameter(1);
    // double error_3 = fDouble -> GetParError(1);

    // Get Sigma
    double sigma_3 = fDouble -> GetParameter(5);

    h1->GetXaxis()->SetRangeUser(lowerBound - 50, upperBound + 50); // Expand view range ±50 beyond fit
    h1->Draw();


    // To visualize the individual Gaussians, define them separately but don't fit
    TF1 *fMain = new TF1(Form("fMain_third_%d", firstXbin), "gaus", lowerBound, upperBound);
    fMain->SetParameters(fDouble->GetParameter(0),
                         fDouble->GetParameter(1),
                         fDouble->GetParameter(2));
    fMain->SetLineColor(kRed);
    fMain->SetLineStyle(2);
    fMain->Draw("same");

    TF1 *fShoulder = new TF1(Form("fShoulder_third_%d", firstXbin), "gaus", lowerBound, upperBound);
    fShoulder->SetParameters(fDouble->GetParameter(3),
                             fDouble->GetParameter(4),
                             fDouble->GetParameter(5));
    fShoulder->SetLineColor(kMagenta);
    fShoulder->SetLineStyle(2);
    fShoulder->Draw("same");


    // if(draw)
    // {
    //     TCanvas* cFit = new TCanvas(Form("cFit_third_%d", firstXbin), "Third Peak Fit", 800, 600);
    //     cFit->cd();  // make this canvas active
    //     h1->GetXaxis()->SetRangeUser(lowerBound - 50, upperBound + 50); // Expand view range ±50 beyond fit
    //     h1->Draw();
    //     fMain->Draw("same");
    //     fShoulder->Draw("same");

    // }

    return sigma_3;
}

// ===== Single Gaussian ===== //

double single_first_peak(int detNum, int firstXbin, int lastXbin, bool quiet, const std::string& plane) {
    TCanvas* cFit = new TCanvas(Form("cSingle_Fit_first_%d", firstXbin), "First Peak Fit", 800, 600);
    cFit->cd();  // make this canvas active

    TH1F* h1 = GetStripHist(detNum, firstXbin, lastXbin, Form("first_%d", firstXbin), plane);

    double lowerBound = 5156 -30; //30
    double upperBound = 5156 +55; // 115

    // --- Combined double Gaussian function: gaus(0)+gaus(3) ---
    TF1 *fSingle = new TF1(Form("fSingle_%d", firstXbin), "gaus(0)", lowerBound, upperBound);

    // Initial parameter guesses
    // [0] amplitude1, [1] mean1, [2] sigma1
    fSingle->SetParameters(
        1500, 5156, 13   // main peak
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
    // double centroid_1S = fSingle -> GetParameter(1);
    // double error_1S = fSingle -> GetParError(1);
    double sigma1 = fSingle -> GetParameter(2);

    // Expand view range ±50 beyond fit
    h1->GetXaxis()->SetRangeUser(lowerBound - 50, upperBound + 50);

    // Draw histogram
    h1->Draw("same");



    return sigma1;
}

double single_second_peak(int detNum, int firstXbin, int lastXbin, bool quiet, const std::string& plane) {
    TCanvas* cFit = new TCanvas(Form("cSingle_Fit_second_%d", firstXbin), "Second Peak Fit", 800, 600);
    cFit->cd();

    TH1F* h1 = GetStripHist(detNum, firstXbin, lastXbin, Form("second_%d", firstXbin), plane);

    double lowerBound = 5485 - 30;  //30 
    double upperBound = 5485 + 100; // 100

    // Single Gaussian
    TF1 *fSingle = new TF1(Form("fSingle_second_%d", firstXbin), "gaus(0)", lowerBound, upperBound);
    fSingle->SetParameters(1400, 5485, 12); // amplitude, mean, sigma
    fSingle->SetParNames("Peak 2 Amp", "Peak 2 Centroid", "Peak 2 Sigma");
    fSingle->SetLineColor(kBlue);
    fSingle->SetLineWidth(2);

    // Choose fit options based on quiet flag
    std::string fitOptions = quiet ? "RQ0+" : "R+";

    // Fit
    h1 -> Fit(fSingle, fitOptions.c_str());
    //h1->Fit(fSingle, "R+");

    // Get Centroid
    //double centroid_2S = fSingle -> GetParameter(1);
    //double error_2S = fSingle -> GetParError(1);
    double sigma2 = fSingle -> GetParameter(2);

    // Expand view range ±50 beyond fit
    h1->GetXaxis()->SetRangeUser(lowerBound - 50, upperBound + 50);

    // Draw histogram and fit
    h1->Draw();
    fSingle->Draw("same");

    return sigma2;
}

double single_third_peak(int detNum, int firstXbin, int lastXbin, bool quiet, const std::string& plane) {
    TCanvas* cFit = new TCanvas(Form("cSingle_Fit_third_%d", firstXbin), "Third Peak Fit", 800, 600);
    cFit->cd();

    TH1F* h1 = GetStripHist(detNum, firstXbin, lastXbin, Form("third_%d", firstXbin), plane);

    double lowerBound = 5804 - 30; //25
    double upperBound = 5804 + 80; // 80

    // Single Gaussian
    TF1 *fSingle = new TF1(Form("fSingle_third_%d", firstXbin), "gaus(0)", lowerBound, upperBound);
    fSingle->SetParameters(1000, 5804, 11); // amplitude, mean, sigma
    fSingle->SetParNames("Peak 3 Amp", "Peak 3 Centroid", "Peak 3 Sigma");
    fSingle->SetLineColor(kBlue);
    fSingle->SetLineWidth(2);

    // Choose fit options based on quiet flag
    std::string fitOptions = quiet ? "RQ0+" : "R+";

    // Fit
    h1 -> Fit(fSingle, fitOptions.c_str());
    //h1->Fit(fSingle, "R+");

    // Get Centroid
    // double centroid_3S = fSingle -> GetParameter(1);
    // double error_3S = fSingle -> GetParError(1);
    double sigma3 = fSingle -> GetParameter(2);

    // Expand view range ±50 beyond fit
    h1->GetXaxis()->SetRangeUser(lowerBound - 50, upperBound + 50);

    // Draw histogram and fit
    h1->Draw();
    fSingle->Draw("same");

    
    return sigma3;
}



void one_fit() {
    int detNum;
    int strip;
    std::cout << "Detector: ";
    std::cin >> detNum;
    std::cout << "Strip: ";
    std::cin >> strip;
    std::string plane;
    std::cout << "Fit the front strips (y/n)? ";
    std::string answer;
    std::cin >> answer;

    if (answer =="y" || answer == "Y")
    {
        plane = "FE";
    }
    else 
    {
        plane = "BE";
    }
   

    // --- Run double Gaussian fits ---
    auto first_double = first_peak_fit(detNum, strip, strip, false, plane);
    auto second_double = second_peak_fit(detNum, strip, strip, false, plane);
    auto third_double = third_peak_fit(detNum, strip, strip, false, plane);

    std::cout << "FWHM from double Gaussians: " << std::endl;
    std::cout<< "FWHM of first peak: " << first_double * 2.35 << std::endl;
    std::cout<< "FWHM of second peak: " << second_double * 2.35 << std::endl;
    std::cout<< "FWHM of third peak: " << third_double * 2.35 << std::endl;

    auto first_single = single_first_peak(detNum, strip, strip, false, plane);
    auto second_single = single_second_peak(detNum, strip, strip, false, plane);
    auto third_single = single_third_peak(detNum, strip, strip, false, plane);

    std::cout << "FWHM from single Gaussians: " << std::endl;
    std::cout<< "FWHM of first peak: " << first_single * 2.35 << std::endl;
    std::cout<< "FWHM of second peak: " << second_single * 2.35 << std::endl;
    std::cout<< "FWHM of third peak: " << third_single * 2.35 << std::endl;
}

void all_fit() {
    int detNum;
    int strip;
    std::cout << "Detector: ";
    std::cin >> detNum;
    std::string plane;
    std::cout << "Fit the front strips (y/n)? ";
    std::string answer;
    std::cin >> answer;

    if (answer =="y" || answer == "Y")
    {
        plane = "FE";
    }
    else 
    {
        plane = "BE";
    }
   

    // --- Run double Gaussian fits ---
    auto first_double = first_peak_fit(detNum, 0, -1, false, plane);
    auto second_double = second_peak_fit(detNum, 0, -1, false, plane);
    auto third_double = third_peak_fit(detNum, 0, -1, false, plane);

    std::cout << "FWHM from double Gaussians: " << std::endl;
    std::cout<< "FWHM of first peak: " << first_double * 2.35 << std::endl;
    std::cout<< "FWHM of second peak: " << second_double * 2.35 << std::endl;
    std::cout<< "FWHM of third peak: " << third_double * 2.35 << std::endl;

    auto first_single = single_first_peak(detNum, 0, -1, false, plane);
    auto second_single = single_second_peak(detNum, 0, -1, false, plane);
    auto third_single = single_third_peak(detNum, 0, -1, false, plane);

    std::cout << "FWHM from single Gaussians: " << std::endl;
    std::cout<< "FWHM of first peak: " << first_single * 2.35 << std::endl;
    std::cout<< "FWHM of second peak: " << second_single * 2.35 << std::endl;
    std::cout<< "FWHM of third peak: " << third_single * 2.35 << std::endl;
}






// void fit_all() {
//     int detNum;
//     std::cout << "Detector: ";
//     std::cin >> detNum;
//     std::string plane;
//     std::cout << "Fit the front strips (y/n)? ";
//     std::string answer;
//     std::cin >> answer;

//     if (answer =="y" || answer == "Y")
//     {
//         plane = "FE";
//     }
//     else 
//     {
//         plane = "BE";
//     }
   
//     std::vector<double> fwhm_peak1, fwhm_peak2, fwhm_peak3;

//     // --- Overall Gaussian fits ---
//     auto first_double = first_peak_fit(detNum, 0, -1, false, plane, true);
//     auto second_double = second_peak_fit(detNum, 0, -1, false, plane, true);
//     auto third_double = third_peak_fit(detNum, 0, -1, false, plane, true);

//     // Individual strip fits
//     for (int strip = 1; strip <= 32; strip++) 
//     {
//         double sigma1 = first_peak_fit(detNum, strip, strip, true, plane, false);
//         fwhm_peak1.push_back(sigma1 * 2.35);

//         double sigma2 = second_peak_fit(detNum, strip, strip, true, plane, false);
//         fwhm_peak2.push_back(sigma2 * 2.35);

//         double sigma3 = third_peak_fit(detNum, strip, strip, true, plane, false);
//         fwhm_peak3.push_back(sigma3 * 2.35);
//     }
    
//     int nStrips = 32;
//     std::vector<double> x_vals;
//     std::vector<double> y_vals1, y_vals2, y_vals3;

//     for (int i = 0; i < nStrips; i++) 
//     {
//         x_vals.push_back(i+1);
//         y_vals1.push_back(fwhm_peak1[i]);
//         y_vals2.push_back(fwhm_peak2[i]);
//         y_vals3.push_back(fwhm_peak3[i]);
//     }

//     // Averages
//     double ave1 = std::accumulate(fwhm_peak1.begin(), fwhm_peak1.end(), 0.0) / fwhm_peak1.size();
//     double ave2 = std::accumulate(fwhm_peak2.begin(), fwhm_peak2.end(), 0.0) / fwhm_peak2.size();
//     double ave3 = std::accumulate(fwhm_peak3.begin(), fwhm_peak3.end(), 0.0) / fwhm_peak3.size();

//     std::cout << "Average FWHM (Peak 1): " << ave1 << std::endl;
//     std::cout << "Average FWHM (Peak 2): " << ave2 << std::endl;
//     std::cout << "Average FWHM (Peak 3): " << ave3 << std::endl;


//     // Peak 1
//     TGraph* gr1 = new TGraph(nStrips, x_vals.data(), y_vals1.data());
//     gr1 -> SetTitle("FWHM of First Peak;Strip Number;FWHM");
//     gr1 -> SetMarkerStyle(20);
//     gr1 -> SetMarkerColor(kBlue);
//     gr1 -> SetLineColor(kBlue);

//     TCanvas* cGraph1 = new TCanvas("cGraph1", "FWHM - Peak 1", 800, 600);
//     gr1 -> Draw("APL");

//     TLine* line1 = new TLine(1, ave1, nStrips, ave1);
//     line1 -> SetLineColor(kRed);
//     line1 -> SetLineStyle(2);
//     line1 -> SetLineWidth(2);

//     cGraph1 -> cd();
//     line1 -> Draw("same");



//     // Peak 2
//     TGraph* gr2 = new TGraph(nStrips, x_vals.data(), y_vals2.data());
//     gr2 -> SetTitle("FWHM of Second Peak;Strip Number;FWHM");
//     gr2 -> SetMarkerStyle(20);
//     gr2 -> SetMarkerColor(kBlue);
//     gr2 -> SetLineColor(kBlue);

//     TCanvas* cGraph2 = new TCanvas("cGraph2", "FWHM - Peak 2", 800, 600);
//     gr2 -> Draw("APL");

//     TLine* line2 = new TLine(1, ave2, nStrips, ave2);
//     line2 -> SetLineColor(kRed);
//     line2 -> SetLineStyle(2);
//     line2 -> SetLineWidth(2);

//     cGraph2 -> cd();
//     line2 -> Draw("same");

//     // Peak 3
//     TGraph* gr3 = new TGraph(nStrips, x_vals.data(), y_vals3.data());
//     gr3 -> SetTitle("FWHM of Third Peak;Strip Number;FWHM");
//     gr3 -> SetMarkerStyle(20);
//     gr3 -> SetMarkerColor(kBlue);
//     gr3 -> SetLineColor(kBlue);

//     TCanvas* cGraph3 = new TCanvas("cGraph3", "FWHM - Peak 3", 800, 600);
//     gr3 -> Draw("APL");

//     TLine* line3 = new TLine(1, ave3, nStrips, ave3);
//     line3 -> SetLineColor(kRed);
//     line3 -> SetLineStyle(2);
//     line3 -> SetLineWidth(2);

//     cGraph3 -> cd();
//     line3 -> Draw("same");

// }

