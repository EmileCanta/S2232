// Written by Sydney Plante, 2026

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
bool PromptCoincidence(TGriffinHit *h1, TGriffinHit *h2) {return ggLow < h1->GetTime() - h2->GetTime() && h1->GetTime() - h2->GetTime() < ggHigh; }

// Griffin-Rcmp gate - returns true if the time between a GRIFFIN hit and an RCMP hit is within the coincidence window
bool PromptCoincidence(TGriffinHit *h1, TRcmpHit *h2) {return grLow < h2->GetTime() - h1->GetTime() && h2->GetTime() - h1->GetTime() < grHigh; } 

// Rcmp-Rcmp gate - retunrs true if the time between two RCMP hits are within the coincidence window
bool PromptCoincidence(TRcmpHit *h1, TRcmpHit *h2){return rrLow < h1->GetTime() - h2->GetTime() && h1->GetTime() - h2->GetTime() < rrHigh;}


// RCMP variables. !!The 0.5 is to make the bins centered!!
double xLowR  = 0;  //RCMP energy goes from 0-10,000 (0 - 10 MeV).
double xHighR = 10000; 
double xLowMul  = 0.5; // Number of hits the detector detects 
double xHighMul = 1100.5; 
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


   fH2[slot]["rcmpEF"] = new TH2F("rcmpEF", "Energy vs Front Strip Number; Front Strips; Energy", frontStripHigh - frontStripLow, frontStripLow, frontStripHigh, xHighR - xLowR, xLowR, xHighR);

   for(int i = 1; i <= 6; ++i) 
   {

      bool flipAxes = (i == 2 || i == 4 || i == 6);
      if (!flipAxes)
      {
         fH2[slot][Form("Det%d", i)] = new TH2F(Form("Det%d", i), Form("RCMP Detector %d; Front; Back", i), stripBins, stripLow, stripHigh, stripBins, stripLow, stripHigh);
      }  
      else 
      {
         fH2[slot][Form("Det%d", i)] = new TH2F(Form("Det%d", i), Form("RCMP Detector %d; Back; Front", i), stripBins, stripLow, stripHigh, stripBins, stripLow, stripHigh);  
      }

      fH2[slot][Form("Det%dFE", i)] = new TH2F(Form("Det%dFE", i), Form("RCMP Detector %d; Front; Front Energy", i), stripBins, stripLow, stripHigh, energyBin/2, energyLow, energyHigh); 
      fH2[slot][Form("Det%dBE", i)] = new TH2F(Form("Det%dBE", i), Form("RCMP Detector %d; Back; Back Energy", i), stripBins, stripLow, stripHigh, energyBin/2, energyLow, energyHigh); 
      fH2[slot][Form("Det%dFECharge", i)] = new TH2F(Form("Det%dFECharge", i), Form("RCMP Detector %d; Front Strips (OBE #); Front Charge", i), stripBins, stripLow, stripHigh, energyBin/4, energyLow, energyHigh);
      fH2[slot][Form("Det%dBECharge", i)] = new TH2F(Form("Det%dBECharge", i), Form("RCMP Detector %d; Back Strips (OBE #); Back Charge", i), stripBins, stripLow, stripHigh, energyBin/4, energyLow, energyHigh);
      
      // Cut out the 2 end strips on each side //
      fH2[slot][Form("Det%dFE_Cut", i)] = new TH2F(Form("Det%dFE_Cut", i), Form("RCMP Detector %d; Front; Front Energy", i), stripBins, stripLow, stripHigh, energyBin/2, energyLow, energyHigh); 
      fH2[slot][Form("Det%dBE_Cut", i)] = new TH2F(Form("Det%dBE_Cut", i), Form("RCMP Detector %d; Back; Back Energy", i), stripBins, stripLow, stripHigh, energyBin/2, energyLow, energyHigh);
      fH2[slot][Form("Det%dFECharge_Cut", i)] = new TH2F(Form("Det%dFECharge_Cut", i), Form("RCMP Detector %d; Front Strips (OBE #); Front Charge", i), stripBins, stripLow, stripHigh, energyBin/4, energyLow, energyHigh);
      fH2[slot][Form("Det%dBECharge_Cut", i)] = new TH2F(Form("Det%dBECharge_Cut", i), Form("RCMP Detector %d; Back Strips (OBE #); Back Charge", i), stripBins, stripLow, stripHigh, energyBin/4, energyLow, energyHigh);
   }

   // Griffin
   fH1[slot]["grifE"] = new TH1F("grifE", "GRIFFIN Energy; Energy (keV); Counts per keV", xHighG - xLowG, xLowG, xHighG);

   // timing spectra
   fH1[slot]["rcmpFrontBackTimeStampDiff"] = new TH1F("rcmpFrontBackTimeStampDiff", "RCMP-RCMP Time Stamp Diff;#DeltaTS_{RCMP-RCMP} (ns)", timeHigh - timeLow, timeLow, timeHigh);
   fH1[slot]["grifGrifTimeStampDiff"] = new TH1F("grifGrifTimeStampDiff", "GRIFFIN-GRIFFIN Time Stamp Diff;#DeltaTS_{GRIFFIN-GRIFFIN}", timeHigh - timeLow, timeLow, timeHigh);
   fH1[slot]["rcmpGriffinTimeStampDiff"] = new TH1F("rcmpGriffinTimeStampDiff", "RCMP-GRIFFIN Time Stamp Diff D5 & D6;#DeltaTS_{RCMP-GRIFFIN} (ns)", timeHigh - timeLow, timeLow, timeHigh);

}

// TODO: Change the function arguments to match the detectors you want to use and the declaration in the header file!
void RCMPHelper::Exec(unsigned int slot, TGriffin& grif, TGriffinBgo& grifBgo, TRcmp& rcmp)
{


   for (int i = 0; i < rcmp.GetFrontMultiplicity(); ++i) 
   {
      auto rcmp1 = rcmp.GetRcmpFrontHit(i);
      int detF = rcmp1 -> GetDetector();

      for (int j = 0; j < rcmp.GetBackMultiplicity(); ++j) 
      {

         auto rcmp2 = rcmp.GetRcmpBackHit(j);
         int detB = rcmp2 -> GetDetector();

         if ( detF == detB) 
         {
            if (!PromptCoincidence(rcmp1,rcmp2)){
               continue;
            }

            int stripF = rcmp1 -> GetStripNumber();
            int stripB = rcmp2 -> GetStripNumber();
            double energyF = rcmp1->GetEnergy();
            double chargeF = rcmp1->GetCharge();
            double energyB = rcmp2->GetEnergy();
            double chargeB = rcmp2->GetCharge();


            if (detF == 2 || detF == 4 || detF == 6 ) 
            {
               fH2[slot].at(Form("Det%d", detF)) -> Fill(backMaps[detF][stripB], frontMaps[detF][stripF]);
            }
            else
            {
               fH2[slot].at(Form("Det%d", detF))->Fill(frontMaps[detF][stripF], backMaps[detF][stripB]);
            }
            
            fH2[slot].at(Form("Det%dFE", detF))->Fill(frontMaps[detF][stripF], energyF);
            fH2[slot].at(Form("Det%dFECharge", detF)) -> Fill(stripF, chargeF);

            fH2[slot].at(Form("Det%dBE", detF))->Fill(backMaps[detB][stripB], energyB);
            fH2[slot].at(Form("Det%dBECharge", detF)) -> Fill(stripB, chargeB);

            fH2[slot].at("rcmpEF")->Fill(rcmp1->GetArrayNumber(), energyF);

            if (stripF >=2 && stripF <=29)
            {
               fH2[slot].at(Form("Det%dFE_Cut", detF))->Fill(frontMaps[detF][stripF], energyF);
               fH2[slot].at(Form("Det%dFECharge_Cut", detF)) -> Fill(stripF, chargeF);
            }

            if (stripB >=2 && stripB <=29)
            {
               fH2[slot].at(Form("Det%dBE_Cut", detF))->Fill(backMaps[detF][stripB], energyB);
               fH2[slot].at(Form("Det%dBECharge_Cut", detF)) -> Fill(stripB, chargeB);
            }

            // Timing Gate Histograms // 
            double diffRStamp = rcmp1->GetTimeStampNs() - rcmp2->GetTimeStampNs();
            fH1[slot].at("rcmpFrontBackTimeStampDiff")->Fill(diffRStamp);

            for (int k = 0; k < grif.GetMultiplicity(); ++k) 
            {
               auto grif1 = grif.GetGriffinHit(k);
               fH1[slot].at("grifE")->Fill(grif1->GetEnergy());

               if (detF == 5 || detF == 6) 
               {
                  fH1[slot].at("rcmpGriffinTimeStampDiff")->Fill(rcmp1->GetTimeStampNs() - grif1->GetTimeStampNs());
               }

            
               for(int l=k+1; l < grif.GetMultiplicity(); ++l)
               {
                  auto grif2 = grif.GetGriffinHit(l);
                  double diffG = grif2->GetTime() - grif1->GetTime();
                  double diffGStamp = grif2->GetTimeStampNs() - grif1->GetTimeStampNs();
         
                  fH1[slot].at("grifGrifTimeStampDiff")->Fill(diffGStamp);
               }
            
            }

         
         
         }

      }


   }

  
}

void RCMPHelper::EndOfSort(std::shared_ptr<std::map<std::string, TList>> list)
{

}




