// Written by Dhruval Shah, 2024

#include "RCMPHelper.hh"
#include "/home/sydney/root/include/TStyle.h" // to use SetOptStat
#include <limits>

// Coincidences Gates (nanoseconds)
double ggLow  = -100.; //Set from grifGrifTimeDiff
double ggHigh = 100.;
double grLow  = 400.; //Set from rcmpGrifTimeDiff
double grHigh = 1100;
double rrLow = -175; //Set from rcmpRcmpTimeDiff
double rrHigh = 175;

// Griffin-Griffin gate - returns true if the time between two GRIFFIN hits is within the coincidence window
bool PromptCoincidence(TGriffinHit *h1, TGriffinHit *h2) 
{
	return ggLow < h1->GetTime() - h2->GetTime() && h1->GetTime() - h2->GetTime() < ggHigh; 
}

// Griffin-Rcmp gate - returns true if the time between a GRIFFIN hit and an RCMP hit is within the coincidence window
bool PromptCoincidence(TGriffinHit *h1, TRcmpHit *h2) 
{
	return grLow < h2->GetTime() - h1->GetTime() && h2->GetTime() - h1->GetTime() < grHigh; 
} 

// Rcmp-Rcmp gate - retunrs true if the time between two RCMP hits are within the coincidence window
bool PromptCoincidence(TRcmpHit *h1, TRcmpHit *h2)
{
   return rrLow < h1->GetTime() - h2->GetTime() && h1->GetTime() - h2->GetTime() < rrHigh;
}

/*PromptCoincidence is here to use to make a histogram based on coincidences. Can you it in a for loop, ex.
if (PromptCoincidence(grifHit1, grifHit2)){
}
or 
if (PromptCoincidence(griftHit, rcmpHit)) {
}
*/

// RCMP variables. !!The 0.5 is to make the bins centered!!
double xLowR  = 0;  //RCMP energy goes from 0-10,000 (0 - 10 MeV).
double xHighR = 10000; 
double xLowMul  = 0.5; // Number of hits the detector detects 
double xHighMul = 1100.5; //***** Why this number? *****
double detLow = 0; // Number of sides of the cube! originally set at -0.5 - 6.5 giving 7 total "detector bins", shifted it to 0.5 to only have 6.
double detHigh = 6; // Used in rcmpEDet
double timeLow = -10000.5; //Used in timeDiff histograms. Just sets the range of the x-axis
double timeHigh = 10000.5;
double frontStripLow = 0; // used in rcmpEFChannel. Looks at all of the front strips of all 6 detectors in one histogram
double frontStripHigh = 192; // 32 x 6 = 192
int stripBins = 32;
double stripLow = 0; // Number of strips per side
double stripHigh = 32; // used in the individual detector histograms
double energyBin = 10000;
double energyLow = 0;
double energyHigh = 10000;


// GRIFFIN variables
double xLowG  = -0.5; //Same as the energy limits for RCMP
double xHighG = 10000.5;

// 2D Hist variables for GRIFFIN

// ===  Channel Mapping  === //
// Detector 1
int frontMap1[32] = {31, 30, 29, 28, 27, 26, 25, 24, 23, 22, 21, 20, 19, 18, 17, 16, 15, 14, 13, 12, 11, 10, 9, 8, 7, 6, 5, 4, 3, 2, 1, 0};
int backMap1[32] = {1, 0, 3, 2, 5, 4, 7, 6, 9, 8, 11, 10, 13, 12, 15, 14, 17, 16, 19, 18, 21, 20, 23, 22, 25, 24, 27, 26, 29, 28, 31, 30}; 

// Detector 2
int frontMap2[32] = {31, 30, 29, 28, 27, 26, 25, 24, 23, 22, 21, 20, 19, 18, 17, 16, 15, 14, 13, 12, 11, 10, 9, 8, 7, 6, 5, 4, 3, 2, 1, 0};
int backMap2[32] = {31, 30, 29, 28, 27, 26, 25, 24, 23, 22, 21, 20, 19, 18, 17, 16, 15, 14, 13, 12, 11, 10, 9, 8, 7, 6, 5, 4, 3, 2, 1, 0};

// Detector 3
int backMap3[32] = {1, 0, 3, 2, 5, 4, 7, 6, 9, 8, 11, 10, 13, 12, 15, 14, 17, 16, 19, 18, 21, 20, 23, 22, 25, 24, 27, 26, 29, 28, 31, 30};
int frontMap3[32] = {1, 0, 3, 2, 5, 4, 7, 6, 9, 8, 11, 10, 13, 12, 15, 14, 17, 16, 19, 18, 21, 20, 23, 22, 25, 24, 27, 26, 29, 28, 31, 30};

// Detector 4
int frontMap4[32] = {0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31};
int backMap4[32] = {31, 30, 29, 28, 27, 26, 25, 24, 23, 22, 21, 20, 19, 18, 17, 16, 15, 14, 13, 12, 11, 10, 9, 8, 7, 6, 5, 4, 3, 2, 1, 0};

// Detector 5
int frontMap5[32] = {0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31};
int backMap5[32] = {30, 31, 28, 29, 26, 27, 24, 25, 22, 23, 20, 21, 18, 19, 16, 17, 14, 15, 12, 13, 10, 11, 8, 9, 6, 7, 4, 5, 2, 3, 0, 1};

// Detector 6
int frontMap6[32] = {0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31};       
int backMap6[32] = {0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31};


std::map<int, int*> frontMaps = {
   {1, frontMap1}, {2, frontMap2}, {3,frontMap3}, {4, frontMap4}, {5, frontMap5}, {6, frontMap6}
}; 

std::map<int, int*> backMaps = {
   {1, backMap1}, {2, backMap2}, {3, backMap3}, {4, backMap4}, {5, backMap5}, {6, backMap6}
}; 

void RCMPHelper::CreateHistograms(unsigned int slot)
{
   // try and get the cycle length if we have a PPG provided
   // only necessary for the first worker, this is shared with all other workers
   if(fPpg != nullptr) {
      // the ODB cycle length is in microseconds!!! and GetTime() is in nanoseconds
      fCycleLength = fPpg->OdbCycleLength();
      if(slot == 0) {
         std::stringstream str;
         str << "Got ODB cycle length " << fCycleLength << " us = " << fCycleLength / 1e6 << " s = " << fCycleLength * 1e3 << " ns" << std::endl;
         std::cerr << str.str();
      }
   } else if(slot == 0) {
      std::stringstream str;
      str << DRED << "No ppg provided, can't fill cycle spectra!" << RESET_COLOR << std::endl;
      std::cerr << str.str();
   }

   // FORMAT: fH1[slot]["histName"]  = new TH1F("histName", "Hist Title; X-axis Title; Y-axis Title", BinCount, xMin, xMax);

   // Rcmp Energy (1D)
   fH1[slot]["rcmpE"]  = new TH1F("rcmpE", "RCMP Pixel Energy; Front Energy (keV); Counts per keV", xHighR - xLowR, xLowR, xHighR);
   //fH1[slot]["rcmpEF"] = new TH1F("rcmpEF", "RCMP Front Energy; Front Energy (keV); Counts per keV", xHighR - xLowR, xLowR, xHighR);


   //FORMAT: fH2[slot]["histName"] = new TH2F("histName", "Hist Title; X-axis Title; Y-axis Title", BinCountX, xMin, xMax, BinCountY, yMin, yMax);

   // Rcmp 2D
   //fH2[slot]["frag"] = new TH2F("frag", "Fragments; Front; Back", detHigh - detLow, detLow, detHigh, detHigh - detLow, detLow, detHigh);
   fH2[slot]["rcmpEDet"] = new TH2F("rcmpEDet", "RCMP Pixel Energy and Detector Number; Energy (keV); Detector Number", xHighR - xLowR, xLowR, xHighR, detHigh - detLow, detLow, detHigh);
   // fH2[slot]["rcmpEFDet"] = new TH2F("rcmpEFDet", "RCMP Front Energy and Detector Number; Energy (keV); Detector Number", xHighR - xLowR, xLowR, xHighR, detHigh - detLow, detLow, detHigh);
   fH2[slot]["rcmpGrifE"] = new TH2F("rcmpGrifE", "RCMP GRIFFIN Energy; RCMP energy (keV); GRIFFIN energy (keV)", xHighR - xLowR, xLowR, xHighR, xHighG - xLowG, xLowG, xHighG);
   fH2[slot]["rcmpGrifE56"] = new TH2F("rcmpGrifE56", "RCMP GRIFFIN Energy; RCMP energy (keV); GRIFFIN energy (keV)", xHighR - xLowR, xLowR, xHighR, xHighG - xLowG, xLowG, xHighG);
   fH2[slot]["rcmpEMaxDet"] = new TH2F("rcmpEMaxDet", "Highest RCMP Pixel Energy and Detector Number; Energy (keV); Detector Number", xHighR - xLowR, xLowR, xHighR, detHigh - detLow, detLow, detHigh);
   fH2[slot]["rcmpENonMinDet"] = new TH2F("rcmpENonMinDet", "RCMP Pixel Energy (excluding lowest energy) and Detector Number; Energy (keV); Detector Number", xHighR - xLowR, xLowR, xHighR, detHigh - detLow, detLow, detHigh);
   fH2[slot]["rcmpEFChannel"] = new TH2F("rcmpEFChannel", "Energy vs Front Strip Number; Front Strips; Energy", frontStripHigh - frontStripLow, frontStripLow, frontStripHigh, xHighR - xLowR, xLowR, xHighR);

   for(int i = 1; i <= 6; ++i) 
   {

      bool flipAxes = (i == 2 || i == 4 || i == 6);
      if (!flipAxes)
      {
         fH2[slot][Form("Det%d", i)] = new TH2F(Form("Det%d", i), Form("RCMP Detector %d Using Pixel Mulitplicity; Front; Back", i), stripBins, stripLow, stripHigh, stripBins, stripLow, stripHigh);
         fH2[slot][Form("Det%dCheck", i)] = new TH2F(Form("Det%dCheck", i), Form("RCMP Detector %d Using Front/Back Multiplicity; Front; Back", i), stripBins, stripLow, stripHigh, stripBins, stripLow, stripHigh);
         //fH2[slot][Form("Det%dCut", i)] = new TH2F(Form("Det%dCut", i), Form("RCMP Detector %d, Energy > 4.9 MeV; Back; Front", i), stripBins, stripLow, stripHigh, stripBins, stripLow, stripHigh);
         // NOTE: The Det%d and Det%dCheck histograms are just ot check if two different ways of filling the histograms are valid! The Det%d histogram should probably have less counts!
      }  
      else 
      {
         fH2[slot][Form("Det%d", i)] = new TH2F(Form("Det%d", i), Form("RCMP Detector %d; Back; Front", i), stripBins, stripLow, stripHigh, stripBins, stripLow, stripHigh);  
         // fH2[slot][Form("Det%dCut", i)] = new TH2F(Form("Det%dCut", i), Form("RCMP Detector %d, E > 5 MeV; Front; Back", i), stripBins, stripLow, stripHigh, stripBins, stripLow, stripHigh); 
         fH2[slot][Form("Det%dCheck", i)] = new TH2F(Form("Det%dCheck", i), Form("RCMP Detector %d Using Front/Back Multiplicity; Back; Front", i), stripBins, stripLow, stripHigh, stripBins, stripLow, stripHigh);
         // fH2[slot][Form("Det%dFE_Old", i)] = new TH2F(Form("Det%dFE_Old", i), Form("RCMP Detector %d Using Pixel Multiplicity; Front; Front Energy", i), stripBins, stripLow, stripHigh, energyBin/4, energyLow, energyHigh); 
         // fH2[slot][Form("Det%dBE_Old", i)] = new TH2F(Form("Det%dBE_Old", i), Form("RCMP Detector %d Using Pixel Multiplicity; Back; Back Energy", i), stripBins, stripLow, stripHigh, energyBin/4, energyLow, energyHigh);
         // fH2[slot][Form("Det%dFECharge_Old", i)] = new TH2F(Form("Det%dFECharge_Old", i), Form("RCMP Detector %d Using Pixel Multiplicity; Front; Front Charge", i), stripBins, stripLow, stripHigh, energyBin, energyLow, energyHigh);
         // fH2[slot][Form("Det%dBECharge_Old", i)] = new TH2F(Form("Det%dBECharge_Old", i), Form("RCMP Detector %d Using Pixel Multiplicity; Back; Back Charge", i), stripBins, stripLow, stripHigh, energyBin, energyLow, energyHigh);
      }
      fH2[slot][Form("Det%dFE", i)] = new TH2F(Form("Det%dFE", i), Form("RCMP Detector %d; Front; Front Energy", i), stripBins, stripLow, stripHigh, energyBin/4, energyLow, energyHigh); 
      fH2[slot][Form("Det%dBE", i)] = new TH2F(Form("Det%dBE", i), Form("RCMP Detector %d; Back; Back Energy", i), stripBins, stripLow, stripHigh, energyBin/4, energyLow, energyHigh); 
      fH2[slot][Form("Det%dFECharge", i)] = new TH2F(Form("Det%dFECharge", i), Form("RCMP Detector %d; Front; Front Charge", i), stripBins, stripLow, stripHigh, energyBin/4, energyLow, energyHigh);
      fH2[slot][Form("Det%dBECharge", i)] = new TH2F(Form("Det%dBECharge", i), Form("RCMP Detector %d; Back; Back Charge", i), stripBins, stripLow, stripHigh, energyBin/4, energyLow, energyHigh);
   }

   // Rcmp sizes
   // fH1[slot]["rcmpFrF"] = new TH1I("rcmpFrF", "Front Fragments; Size of Front Frag Vector; Counts per bin", xHighMul - xLowMul, xLowMul, xHighMul);
   // fH1[slot]["rcmpFrB"] = new TH1I("rcmpFrB", "Back Fragments; Size of Back Frag Vector; Counts per bin", xHighMul - xLowMul, xLowMul, xHighMul);
   // fH1[slot]["rcmpMuF"] = new TH1I("rcmpMuF", "Front Multiplicity; Size of Front Multiplicity Vector; Counts per bin", xHighMul - xLowMul, xLowMul, xHighMul);
   // fH1[slot]["rcmpMuB"] = new TH1I("rcmpMuB", "Front Hits/Multiplicity; Size of Back Multiplicity Vector; Counts per bin", xHighMul - xLowMul, xLowMul, xHighMul);
   // fH1[slot]["rcmpMuPi"] = new TH1I("rcmpMuPi", "Front+Back (Pixel) Multiplicity; Size of Front+Back Multiplicity Vector; Counts per bin", xHighMul - xLowMul, xLowMul, xHighMul);
   // fH1[slot]["rcmpMuZ"] = new TH1I("rcmpMuZ", "Zero Multiplicity; Size of Pixel Multiplicity was zero; Counts", 10000, 0.5, 10000.5);

   // Griffin
   fH1[slot]["grifE"] = new TH1F("grifE", "GRIFFIN Energy; Energy (keV); Counts per keV", xHighG - xLowG, xLowG, xHighG);
   // fH1[slot]["grifEAddback"] = new TH1F("grifEAddback", "GRIFFIN Energy (Addback); Energy (keV); Counts per keV", xHighG - xLowG, xLowG, xHighG);
   // fH1[slot]["grifESupp"] = new TH1F("grifESupp", "GRIFFIN Energy (Suppressed); Energy (keV); Counts per keV", xHighG - xLowG, xLowG, xHighG);

   // timing spectra
   fH1[slot]["rcmpFrontBackTimeDiff"] = new TH1F("rcmpFrontBackTimeDiff", "RCMP-RCMP Time Diff;#DeltaT_{RCMP-RCMP} (ns)", timeHigh - timeLow, timeLow, timeHigh);
   fH1[slot]["rcmpFrontBackTimeStampDiff"] = new TH1F("rcmpFrontBackTimeStampDiff", "RCMP-RCMP Time Stamp Diff;#DeltaTS_{RCMP-RCMP} (ns)", timeHigh - timeLow, timeLow, timeHigh);
   fH1[slot]["grifGrifTimeDiff"] = new TH1F("grifGrifTimeDiff", "GRIFFIN-GRIFFIN Time Diff;#DeltaT_{GRIFFIN-GRIFFIN}", timeHigh - timeLow, timeLow, timeHigh);
   fH1[slot]["grifGrifTimeStampDiff"] = new TH1F("grifGrifTimeStampDiff", "GRIFFIN-GRIFFIN Time Stamp Diff;#DeltaTS_{GRIFFIN-GRIFFIN}", timeHigh - timeLow, timeLow, timeHigh);
   fH1[slot]["rcmpGriffinTimeDiff"] = new TH1F("rcmpGriffinTimeDiff", "RCMP-GRIFFIN Time Diff;#DeltaT_{RCMP-GRIFFIN} (ns)", timeHigh - timeLow, timeLow, timeHigh);
   fH1[slot]["rcmpGriffinTimeStampDiff"] = new TH1F("rcmpGriffinTimeStampDiff", "RCMP-GRIFFIN Time Stamp Diff;#DeltaTS_{RCMP-GRIFFIN} (ns)", timeHigh - timeLow, timeLow, timeHigh);
   // fH1[slot]["rcmpGriffinTimeDiff6"] = new TH1F("rcmpGriffinTimeDiff6", "RCMP-GRIFFIN Time Diff 6;#DeltaT_{RCMP-GRIFFIN} (ns)", timeHigh - timeLow, timeLow, timeHigh);
   // fH1[slot]["rcmpGriffinTimeStampDiff6"] = new TH1F("rcmpGriffinTimeStampDiff6", "RCMP-GRIFFIN Time Stamp Diff 6;#DeltaTS_{RCMP-GRIFFIN} (ns)", timeHigh - timeLow, timeLow, timeHigh);
   // fH2[slot]["timeDiff"]      = new TH2F("timeDiff", "Time Difference; Time (ns); Detector Number", 4000, -2000, 2000, 7, -0.5, 6.5);
   // fH2[slot]["griffinrcmpTime"]     = new TH2F("griffinrcmpTime", "GRIFFIN crystal vs. GRIFFIN-rcmp timing (suppressed addback);#Deltat_{GRIFFIN-rcmp}", 2000, -1000., 1000., 64, 0.5, 64.5);
   // fH2[slot]["griffinGriffinTS"]   = new TH2F("griffinGriffinTS", "GRIFFIN crystal vs. GRIFFIN-GRIFFIN timestamp difference (suppressed addback);#DeltaTS_{GRIFFIN-GRIFFIN}", 2000, -1000., 1000., 64, 0.5, 64.5);
   // fH2[slot]["griffinGriffinTime"] = new TH2F("griffinGriffinTime", "GRIFFIN crystal vs. GRIFFIN-GRIFFIN timing (suppressed addback);#Deltat_{GRIFFIN-GRIFFIN}", 2000, -1000., 1000., 64, 0.5, 64.5);

}

// TODO: Change the function arguments to match the detectors you want to use and the declaration in the header file!
void RCMPHelper::Exec(unsigned int slot, TGriffin& grif, TGriffinBgo& grifBgo, TRcmp& rcmp)
{
   /********** RCMP **********/

   // determine lowest and highest RCMP energy
   double highE = std::numeric_limits<double>::lowest(); //min returns the smallest positive value for floating point numbers
   double lowE = std::numeric_limits<double>::max();
   for(int i = 0; i < rcmp.GetPixelMultiplicity(); ++i)
   {
      auto* rcmp1 = rcmp.GetRcmpHit(i);
      if(rcmp1->GetEnergy() < lowE)  {  lowE = rcmp1->GetEnergy(); }
      if(rcmp1->GetEnergy() > highE) { highE = rcmp1->GetEnergy(); }
   }

   // Fill RCMP Pixel Multiplicity spectra
   for(int i = 0; i < rcmp.GetPixelMultiplicity(); ++i)
   {
      auto rcmp1 = rcmp.GetRcmpHit(i);

      fH1[slot].at("rcmpE")->Fill(rcmp1->GetEnergy());
      fH2[slot].at("rcmpEDet")->Fill(rcmp1->GetEnergy(), rcmp1->GetDetector());
      
      if(lowE < highE)
      {
	      // this means we have at least two hits in RCMP
	      // otherwise both lowE and highE would be set to the single energy found
	      if(std::abs(rcmp1->GetEnergy() - lowE) > 0.1)
         {
		      // not the lowest energy
		      fH2[slot].at("rcmpENonMinDet")->Fill(rcmp1->GetEnergy(), rcmp1->GetDetector());
	      } 
	      if(std::abs(rcmp1->GetEnergy() - highE) < 0.1)
         {
		      // highest energy only
		      fH2[slot].at("rcmpEMaxDet")->Fill(rcmp1->GetEnergy(), rcmp1->GetDetector());
	      }
      }

      // loop over GRIFFIN Multiplicity inside RCMP Pixel Multiplicity
      for(int j = 0; j < grif.GetMultiplicity(); ++j) 
      {
         auto grif1 = grif.GetGriffinHit(j);

      

         // Filling the following 2D histograms only for RCS05 and RCS06 -> Due to the two of them being upstream with no energy loss due to target
         if(rcmp1->GetDetector() == 5 || rcmp1->GetDetector() == 6) // the || mean detector 5 OR 6
         {
         	// fH2[slot].at("rcmpGrifE56")->Fill(rcmp1->GetEnergy(), grif1->GetEnergy());
         	fH1[slot].at("rcmpGriffinTimeDiff")->Fill(rcmp1->GetTime() - grif1->GetTime());
	         fH1[slot].at("rcmpGriffinTimeStampDiff")->Fill(rcmp1->GetTimeStampNs() - grif1->GetTimeStampNs());
            
         }

         // if(rcmp1->GetDetector() == 5) 
         // {
         // 	fH2[slot].at("rcmpGrifE56")->Fill(rcmp1->GetEnergy(), grif1->GetEnergy());
         // 	fH1[slot].at("rcmpGriffinTimeDiff5")->Fill(rcmp1->GetTime() - grif1->GetTime());
	      //    fH1[slot].at("rcmpGriffinTimeStampDiff5")->Fill(rcmp1->GetTimeStampNs() - grif1->GetTimeStampNs());
         //    fH1[slot].at("grifGrifTimeDiff")->Fill
         // }
         // if(rcmp1->GetDetector() == 6) 
         // {
         // 	fH2[slot].at("rcmpGrifE56")->Fill(rcmp1->GetEnergy(), grif1->GetEnergy());
         // 	fH1[slot].at("rcmpGriffinTimeDiff6")->Fill(rcmp1->GetTime() - grif1->GetTime());
	      //    fH1[slot].at("rcmpGriffinTimeStampDiff6")->Fill(rcmp1->GetTimeStampNs() - grif1->GetTimeStampNs());
         //    fH1[slot].at("grifGrifTimeDiff")->Fill
         // }

         
         // Looping at the grif-grif Time diff
         for(int k=j+1; k<grif.GetMultiplicity(); ++k){
            auto grif2 = grif.GetGriffinHit(k);
            double diffG = grif2->GetTime() - grif1->GetTime();
            // printf("Grif time 1; %lf   Grif time 2; %lf   Diff; %lf \n\n\n", grif1->GetTime(), grif2->GetTime(), diffG);
            // std::cout << "Grif time 1: " << grif1->GetTime() << "Grif time 2: " << grif2->GetTime()  << "Difference: "  << diffG << endl;
            fH1[slot].at("grifGrifTimeDiff")->Fill(diffG);
            // fH1[slot].at("grifGrifTimeDiff")->Fill(-diffG);

            double diffGStamp = grif2->GetTimeStampNs() - grif1->GetTimeStampNs();
            // printf("Grif time stamp 1; %lld   Grif time stamp 2; %lld   Diff; %lf \n\n\n", grif1->GetTimeStampNs(), grif2->GetTimeStampNs(), diffGStamp);
            fH1[slot].at("grifGrifTimeStampDiff")->Fill(diffGStamp);
            fH1[slot].at("grifGrifTimeStampDiff")->Fill(-diffGStamp);
         }
      }
      
   }
   // Coincidence in rcmp and grif
   for(int i = 0; i < rcmp.GetPixelMultiplicity(); ++i)
   {
      auto rcmp1 = rcmp.GetRcmpHit(i);

      for(int j = 0; j < grif.GetMultiplicity(); ++j)
      {
         auto grif1 = grif.GetGriffinHit(j);
         fH2[slot].at("rcmpGrifE")->Fill(rcmp1->GetEnergy(), grif1->GetEnergy());
      }
   }
   // ADD IN REQUIRMENT ABOUT THE FRONT AND BACK BEING IN COINCIDENCE

   // RCMP Heat map filling. Broken up into different loops based on how the detectors were installed

   for(int i = 0; i < rcmp.GetPixelMultiplicity(); ++i) 
   {
      auto rcmp1 = rcmp.GetRcmpHit(i);
      int det = rcmp1 -> GetDetector();
      int stripF = rcmp1 -> GetFrontStripNumber();
      int stripB = rcmp1 -> GetBackStripNumber();

      if (det == 2 || det == 4 || det == 6)
            {
               fH2[slot].at(Form("Det%d", det)) -> Fill(backMaps[det][stripB], frontMaps[det][stripF]);
            }

            else 
            {
               fH2[slot].at(Form("Det%d", det))->Fill(frontMaps[det][stripF], backMaps[det][stripB]);
            }
      
      //fH2[slot].at(Form("Det%d", det))->Fill(frontMaps[det][stripF], backMaps[det][stripB]);	  

      // fH2[slot].at(Form("Det%dFE_Old", det))->Fill(frontMaps[det][stripF], rcmp1->GetEnergy());	  
      // fH2[slot].at(Form("Det%dBE_Old", det))->Fill(backMaps[det][stripB], rcmp1->GetBackEnergy());	
         
      // fH2[slot].at(Form("Det%dFECharge", det)) -> Fill(stripF, rcmp1->GetCharge());
      // fH2[slot].at(Form("Det%dBECharge", det)) -> Fill(stripB, rcmp1->GetBackCharge());


      // if(rcmp1->GetEnergy() > 5000.) // Look at the hits that correspond the above 5MeV
      // {
      //    fH2[slot].at(Form("Det%dCut", rcmp1->GetDetector()))->Fill(rcmp1->GetFrontStripNumber(), rcmp1->GetBackStripNumber());	   
      // }
   }

   // ====== OLD WAY OF DOING THE FE AND BE HISTOGRAMS!! ====== // 
   // for (int i = 0; i < rcmp.GetPixelMultiplicity(); ++i) 
   // {
   //    auto rcmp1 = rcmp.GetRcmpHit(i);
   //    int det = rcmp1 -> GetDetector();

   //    fH2[slot].at(Form("Det%dFE_Old", det))->Fill(rcmp1 -> GetFrontStripNumber(), rcmp1->GetEnergy());	  
   //    fH2[slot].at(Form("Det%dBE_Old", det))->Fill(rcmp1 -> GetBackStripNumber(), rcmp1->GetEnergy());	

   //    fH2[slot].at(Form("Det%dFECharge_Old", det))->Fill(rcmp1 -> GetFrontStripNumber(), rcmp1->GetCharge());	  
   //    fH2[slot].at(Form("Det%dBECharge_Old", det))->Fill(rcmp1 -> GetBackStripNumber(), rcmp1->GetCharge());
   // }
   // PROBLEM: GetPixelMulitplicity only saves the front energy, time, etc. so no back strip information is saved!


   
   for (int i = 0; i < rcmp.GetFrontMultiplicity(); ++i)
   {
      auto frontHit = rcmp.GetRcmpFrontHit(i);
      int detF = frontHit -> GetDetector();
      int stripF = frontHit -> GetStripNumber();

      for (int j = 0; j < rcmp.GetBackMultiplicity(); ++j)
      {
         auto backHit = rcmp.GetRcmpBackHit(j);
         int detB = backHit -> GetDetector();
         int stripB = backHit -> GetStripNumber();
         if (detF != detB) continue;

         if (PromptCoincidence(frontHit, backHit))
         {
            if (detF == 2 || detF == 4 || detF == 6)
            {
               fH2[slot].at(Form("Det%dCheck", detF)) -> Fill(backMaps[detB][stripB], frontMaps[detF][stripF]);
            }

            else 
            {
               fH2[slot].at(Form("Det%dCheck", detF)) -> Fill(frontMaps[detF][stripF], backMaps[detB][stripB]);
            }
            //fH2[slot].at(Form("Det%dCheck", detF)) -> Fill(frontMaps[detF][stripF], backMaps[detB][stripB]);
         }
      }
   }

   // ====== NEW WAY OF DOING THE FE AND BE HISTOGRAMS!! ====== //
   for (int i = 0; i<rcmp.GetFrontMultiplicity(); ++i)
   {
      auto rcmp1 = rcmp.GetRcmpFrontHit(i);
      int det = rcmp1 -> GetDetector();
      int stripF = rcmp1 -> GetStripNumber();

      fH2[slot].at(Form("Det%dFE", det))->Fill(frontMaps[det][stripF], rcmp1->GetEnergy());
      // fH2[slot].at(Form("Det%dFE", det))->Fill(stripF, rcmp1->GetEnergy());
      fH2[slot].at(Form("Det%dFECharge", det)) -> Fill(stripF, rcmp1->GetCharge());
   }

   for (int i = 0; i< rcmp.GetBackMultiplicity(); ++i) 
   {
      auto rcmp1 = rcmp.GetRcmpBackHit(i);
      int det = rcmp1 -> GetDetector();
      int stripB = rcmp1 -> GetStripNumber();

      fH2[slot].at(Form("Det%dBE", det))->Fill(backMaps[det][stripB], rcmp1->GetEnergy());
      // fH2[slot].at(Form("Det%dBE", det))->Fill(stripB, rcmp1->GetEnergy());
      fH2[slot].at(Form("Det%dBECharge", det)) -> Fill(stripB, rcmp1->GetCharge());
   }
   // PROBLEM: Neighbours being affected by one event might be counted as separate events!


   // Fill RCMP Front Multiplicity spectra
   for(int i = 0; i < rcmp.GetFrontMultiplicity(); ++i)
   {
      auto rcmp1 = rcmp.GetRcmpFrontHit(i);
      
      fH2[slot].at("rcmpEFChannel")->Fill(rcmp1->GetArrayNumber(), rcmp1->GetEnergy());

      for(int j = 0; j < rcmp.GetBackMultiplicity(); ++j)
      {
         auto rcmp2 = rcmp.GetRcmpBackHit(j);
         
         if(rcmp1->GetDetector() == rcmp2->GetDetector())
         {
            double diffR = rcmp1->GetTime() - rcmp2->GetTime();
            // printf("RCMP time 1; %lf   RCMP time 2; %lf   Diff; %lf \n\n\n", rcmp1->GetTime(), rcmp2->GetTime(), diffR);
            // std::cout << "RCMP time 1: " << rcmp1->GetTime() << "RCMP time 2: " << rcmp2->GetTime()  << "Difference: "  << diffR << endl;
            fH1[slot].at("rcmpFrontBackTimeDiff")->Fill(diffR);
            double diffRStamp = rcmp1->GetTimeStampNs() - rcmp2->GetTimeStampNs();
            fH1[slot].at("rcmpFrontBackTimeStampDiff")->Fill(diffRStamp);
            // printf("RCMP time stamp 1; %lld   RCMP time stamp 2; %lld   Diff; %lf \n\n\n", rcmp1->GetTimeStampNs(), rcmp2->GetTimeStampNs(), diffRStamp);
         }
      }
   }

   /********** GRIFFIN **********/

   // loop over unsuppressed griffin hits
   for(int g = 0; g < grif.GetMultiplicity(); ++g)
   {
      auto grif1 = grif.GetGriffinHit(g);
      fH1[slot].at("grifE")->Fill(grif1->GetEnergy());
   }

   // // loop over suppressed griffin hits
   // for(int g = 0; g < grif.GetSuppressedMultiplicity(&grifBgo); ++g)
   // {
   //    auto grif1 = grif.GetSuppressedHit(g);
   //    fH1[slot].at("grifESupp")->Fill(grif1->GetEnergy());
   // }

   // // loop over addback griffin hits
   // for(int g = 0; g < grif.GetAddbackMultiplicity(); ++g)
   // {
   //    auto grif1 = grif.GetAddbackHit(g);
   //    fH1[slot].at("grifEAddback")->Fill(grif1->GetEnergy());
   // }

}

void RCMPHelper::EndOfSort(std::shared_ptr<std::map<std::string, TList>> list)
{

}




