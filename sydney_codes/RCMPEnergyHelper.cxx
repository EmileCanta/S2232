// Written by Sydney Plante, 2025

#include "RCMPEnergyHelper.hh"
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
// Rcmp-Rcmp gate - returns true if the time between two RCMP hits are within the coincidence window
bool PromptCoincidence(TRcmpHit *h1, TRcmpHit *h2)
{
   return rrLow < h1->GetTime() - h2->GetTime() && h1->GetTime() - h2->GetTime() < rrHigh;
}
// Neighbour check -> A general idea but needs work!
// bool NeighbourCoincidence(TRcmpHit* h1, TRcmpHit* h2, double timeWindow)
// {
//     // check same detector
//     if(h1->GetDetector() != h2->GetDetector()) return false;

//     // check neighbouring strips
//     if(std::abs(h1->GetFrontStrip() - h2->GetFrontStrip()) != 1) return false;

//     // check time coincidence
//     double dt = h1->GetTime() - h2->GetTime();
//     if(dt < -timeWindow || dt > timeWindow) return false;

//     return true;
// }



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
double mult = 50;


// GRIFFIN variables
double xLowG  = -0.5; //Same as the energy limits for RCMP
double xHighG = 10000.5;

// ===  Channel Mapping  === //
// Detector 1
//int frontMap1[32] = {0, 1, 2, 3, 4, 5, 6, 7, 8, 9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30,31};
int frontMap1[32] = {31, 30, 29, 28, 27, 26, 25, 24, 23, 22, 21, 20, 19, 18, 17, 16, 15, 14, 13, 12, 11, 10, 9, 8, 7, 6, 5, 4, 3, 2, 1, 0};

int backMap1[32] = {1, 0, 3, 2, 5, 4, 7, 6, 9, 8, 11, 10, 13, 12, 15, 14, 17, 16, 19, 18, 21, 20, 23, 22, 25, 24, 27, 26, 29, 28, 31, 30}; 

// Detector 2
//int frontMap2[32] = {0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31};
int frontMap2[32] = {31, 30, 29, 28, 27, 26, 25, 24, 23, 22, 21, 20, 19, 18, 17, 16, 15, 14, 13, 12, 11, 10, 9, 8, 7, 6, 5, 4, 3, 2, 1, 0};

//int backMap2[32] = {0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31};
int backMap2[32] = {31, 30, 29, 28, 27, 26, 25, 24, 23, 22, 21, 20, 19, 18, 17, 16, 15, 14, 13, 12, 11, 10, 9, 8, 7, 6, 5, 4, 3, 2, 1, 0};

// Detector 3
int backMap3[32] = {1, 0, 3, 2, 5, 4, 7, 6, 9, 8, 11, 10, 13, 12, 15, 14, 17, 16, 19, 18, 21, 20, 23, 22, 25, 24, 27, 26, 29, 28, 31, 30};

int frontMap3[32] = {1, 0, 3, 2, 5, 4, 7, 6, 9, 8, 11, 10, 13, 12, 15, 14, 17, 16, 19, 18, 21, 20, 23, 22, 25, 24, 27, 26, 29, 28, 31, 30};


int frontMap4[32] = {0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31};

//int backMap4[32] = {0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31};
int backMap4[32] = {31, 30, 29, 28, 27, 26, 25, 24, 23, 22, 21, 20, 19, 18, 17, 16, 15, 14, 13, 12, 11, 10, 9, 8, 7, 6, 5, 4, 3, 2, 1, 0};

// Detector 5
int frontMap5[32] = {0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31};
		
// int backMap5[32] =  {1, 0, 3, 2, 5, 4, 7, 6, 9, 8, 11, 10, 13, 12, 15, 14,  17, 16, 19, 18, 21, 20, 23, 22,  25, 24, 27, 26, 29, 28, 31, 30};
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


// 2D Hist variables for GRIFFIN

void RCMPEnergyHelper::CreateHistograms(unsigned int slot)
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

   fH1[slot]["PixelMult"]  = new TH1F("PixelMult", "RCMP Pixel Multiplicity; Pixel Multiplicity; Counts", mult, 0, mult);
   fH2[slot]["rcmpEF"] = new TH2F("rcmpEF", "Energy vs Front Strip Number; Front Strips; Front Energy", frontStripHigh - frontStripLow, frontStripLow, frontStripHigh, xHighR - xLowR, xLowR, xHighR);
   fH2[slot]["rcmpEB"] = new TH2F("rcmpEB", "Energy vs Back Strip Number; Back Strips; Back Energy", frontStripHigh - frontStripLow, frontStripLow, frontStripHigh, xHighR - xLowR, xLowR, xHighR);
    
   for(int i = 4; i <= 6; ++i) 
   {
      bool flipAxes = (i == 2 || i == 4 || i == 6); 

      if (!flipAxes)
      {
         fH2[slot][Form("Det%d", i)] = new TH2F(Form("Det%d", i), Form("RCMP Detector %d; Front; Back", i), stripBins, stripLow, stripHigh, stripBins, stripLow, stripHigh);
         fH2[slot][Form("Det%dCut", i)] = new TH2F(Form("Det%dCut", i), Form("RCMP Detector %d, Energy > 4.9 MeV; Back; Front", i), stripBins, stripLow, stripHigh, stripBins, stripLow, stripHigh);
      }  
      else 
      {
         fH2[slot][Form("Det%d", i)] = new TH2F(Form("Det%d", i), Form("RCMP Detector %d; Back; Front", i), stripBins, stripLow, stripHigh, stripBins, stripLow, stripHigh);
         fH2[slot][Form("Det%dCut", i)] = new TH2F(Form("Det%dCut", i), Form("RCMP Detector %d, Energy > 4.9 MeV; Back; Front", i), stripBins, stripLow, stripHigh, stripBins, stripLow, stripHigh);
      } 
        fH2[slot][Form("Det%dFE", i)] = new TH2F(Form("Det%dFE", i), Form("RCMP Detector %d; Front; Front Energy", i), stripBins, stripLow, stripHigh, energyBin/5, energyLow, energyHigh); 
        fH2[slot][Form("Det%dBE", i)] = new TH2F(Form("Det%dBE", i), Form("RCMP Detector %d; Back; Back Energy", i), stripBins, stripLow, stripHigh, energyBin/5, energyLow, energyHigh); 
        //fH2[slot][Form("Det%dFECharge", i)] = new TH2F(Form("Det%dFECharge", i), Form("RCMP Detector %d; Front; Front Charge", i), stripBins, stripLow, stripHigh, energyBin/4, energyLow, energyHigh);
        //fH2[slot][Form("Det%dBECharge", i)] = new TH2F(Form("Det%dBECharge", i), Form("RCMP Detector %d; Back; Back Charge", i), stripBins, stripLow, stripHigh, energyBin/4, energyLow, energyHigh);
        fH2[slot][Form("Det%dFECut", i)] = new TH2F(Form("Det%dFECut", i), Form("RCMP Detector %d, Energy > 4.9 MeV; Front; Front Energy", i), stripBins, stripLow, stripHigh, (energyBin-5000)/4, 3500, energyHigh - 1500); 
        fH2[slot][Form("Det%dBECut", i)] = new TH2F(Form("Det%dBECut", i), Form("RCMP Detector %d, Energy > 4.9 MeV; Back; Back Energy", i), stripBins, stripLow, stripHigh, (energyBin-5000)/4, 3500, energyHigh - 1500); 
        fH2[slot][Form("Det%dE", i)] = new TH2F(Form("Det%dE", i), Form("RCMP Detector %d; Front Energy; Back Energy", i), energyBin/4, energyLow, energyHigh, energyBin/4, energyLow, energyHigh);
        
        if (i == 5)
        {
         //fH2[slot][Form("Det%dF_B1", i)] = new TH2F(Form("Det%dF_B1", i), Form("RCMP Detector %d, Back Strip 1; Front; Front Energy", i), stripBins, stripLow, stripHigh, energyBin/10, energyLow, energyHigh);
         fH2[slot][Form("Det%dFE_10-20", i)] = new TH2F(Form("Det%dFE_10-20", i), Form("RCMP Detector %d, Back Strips 10-20; Front; Front Energy", i), stripBins, stripLow, stripHigh, energyBin/5, energyLow, energyHigh);
      }
        if (i == 6) 
        {
         //fH2[slot][Form("Det%dB_F32", i)] = new TH2F(Form("Det%dB_F32", i), Form("RCMP Detector %d, Front Strip 32; Back; Back Energy", i), stripBins, stripLow, stripHigh, energyBin/10, energyLow, energyHigh);
         fH2[slot][Form("Det%dBE_10-20", i)] = new TH2F(Form("Det%dBE_10-20", i), Form("RCMP Detector %d, Front Strips 10-20; Back; Back Energy", i), stripBins, stripLow, stripHigh, energyBin/5, energyLow, energyHigh);
        }
   }



   

}

void RCMPEnergyHelper::Exec(unsigned int slot, TGriffin& grif, TGriffinBgo& grifBgo, TRcmp& rcmp)
{
   int pixel = rcmp.GetPixelMultiplicity();
   fH1[slot].at("PixelMult") -> Fill(pixel);

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
         if(detF < 4 || detF > 6) continue; // Delete when looking at all detectors

         if (PromptCoincidence(frontHit, backHit))
         {
            //fH2[slot].at(Form("Det%d", detF)) -> Fill(frontMaps[detF][stripF], backMaps[detB][stripB]);
            fH2[slot].at(Form("Det%dFE", detF))->Fill(frontMaps[detF][stripF], frontHit->GetEnergy());
            //fH2[slot].at(Form("Det%dFECharge", detF)) -> Fill(stripF, frontHit->GetCharge());
            fH2[slot].at(Form("Det%dBE", detF))->Fill(backMaps[detF][stripB], backHit->GetEnergy());
            //fH2[slot].at(Form("Det%dBECharge", detF)) -> Fill(stripB, backHit->GetCharge());
            fH2[slot].at(Form("Det%dE", detF))->Fill(frontHit->GetEnergy(), backHit ->GetEnergy());

            fH2[slot].at("rcmpEF")->Fill(frontHit->GetArrayNumber(), frontHit->GetEnergy());
            fH2[slot].at("rcmpEB")->Fill(backHit->GetArrayNumber(), backHit->GetEnergy());

            if (detF == 6 && frontMaps[6][stripF] >= 10 && frontMaps[6][stripF] <= 20) //(detF == 5 && (frontMaps[5][stripF] == 31 || frontMaps[5][stripF] == 30 || frontMaps[5][stripF] == 29 ) )
            {
               fH2[slot].at(Form("Det%dBE_10-20", detF))->Fill(backMaps[detF][stripB], backHit->GetEnergy());
            }

            // if (detF == 5 && frontMaps[5][stripF] == 31) 
            // {
            //    fH2[slot].at(Form("Det%dB_F32", detF))->Fill(backMaps[detF][stripB], backHit->GetEnergy());
            // }



            if (detF == 5 && backMaps[5][stripB] >= 10 && backMaps[5][stripB] <= 20 ) //(detF == 6 && (backMaps[6][stripB] == 0 || backMaps[6][stripB] == 1 || backMaps[6][stripB] == 2 ))
            {
               fH2[slot].at(Form("Det%dFE_10-20", detF))->Fill(frontMaps[detF][stripF], frontHit->GetEnergy());
            }

            // if (detF == 6 && backMaps[6][stripB] == 0 ) 
            // {
            //    fH2[slot].at(Form("Det%dF_B1", detF))->Fill(frontMaps[detF][stripF], frontHit->GetEnergy());
            // }



            if (detF == 2 || detF == 4 || detF == 6)
            {
               fH2[slot].at(Form("Det%d", detF)) -> Fill(backMaps[detB][stripB], frontMaps[detF][stripF]);
            }

            else 
            {
               fH2[slot].at(Form("Det%d", detF)) -> Fill(frontMaps[detF][stripF], backMaps[detB][stripB]);
            }



            if(frontHit->GetEnergy() > 4900)
            {
               fH2[slot].at(Form("Det%dFECut", detF))->Fill(frontMaps[detF][stripF], frontHit->GetEnergy());
            }

            if(backHit->GetEnergy() > 4900)
            {
               fH2[slot].at(Form("Det%dBECut", detF))->Fill(backMaps[detF][stripB], backHit->GetEnergy());
            }

            if (frontHit->GetEnergy() > 4900 && backHit->GetEnergy() > 4900)
            {
               if (detF == 2 || detF == 4 || detF == 6)
               {
                  fH2[slot].at(Form("Det%dCut", detF)) -> Fill(backMaps[detB][stripB], frontMaps[detF][stripF]);
               }
               else 
               {
                  fH2[slot].at(Form("Det%dCut", detF)) -> Fill(frontMaps[detF][stripF], backMaps[detB][stripB]);
               }
            }

         }
      }
   }

   /*
   Idea for filtering out neighbours in the histogram:

   for (int i = 0; i < rcmp.GetFrontMultiplicity(); ++i) {
    auto frontHit = rcmp.GetRcmpFrontHit(i);
    bool hasNeighbour = false;

    for (int j = 0; j < rcmp.GetFrontMultiplicity(); ++j) {
        if(i == j) continue; // don't compare the hit to itself
        auto otherHit = rcmp.GetRcmpFrontHit(j);

        if(NeighbourCoincidence(frontHit, otherHit)) {
            hasNeighbour = true;
            break; // no need to check further, neighbour found
        }
    }

    if(!hasNeighbour) {
        // this hit has no neighbour: fill histograms normally
    } else {
        // hit has neighbour: optionally apply energy filter or take the highest energy
    }
}
   ====== More Filled In ======

   TRcmpHit* highestEnergyHit = frontHit;  // start with the current hit
    for (int j = 0; j < rcmp.GetFrontMultiplicity(); ++j) {
        if (i == j) continue; // skip self
        auto otherHit = rcmp.GetRcmpFrontHit(j);

        if (NeighbourCoincidence(frontHit, otherHit)) {
            if (otherHit->GetEnergy() > highestEnergyHit->GetEnergy()) {
                highestEnergyHit = otherHit; // update if higher energy neighbour
            }
        }
    }

    // Only fill histograms if this hit is the highest-energy among neighbours
    if(highestEnergyHit == frontHit) {
        for (int j = 0; j < rcmp.GetBackMultiplicity(); ++j) {
            auto backHit = rcmp.GetRcmpBackHit(j);
            int detB = backHit->GetDetector();
            int stripB = backHit->GetStripNumber();

            if (detF != detB) continue;
            if (!PromptCoincidence(frontHit, backHit)) continue;

            // Now fill your histograms exactly as before
            fH2[slot].at(Form("Det%dFE", detF))->Fill(frontMaps[detF][stripF], frontHit->GetEnergy());
            fH2[slot].at(Form("Det%dBE", detF))->Fill(backMaps[detF][stripB], backHit->GetEnergy());
            fH2[slot].at(Form("Det%dE", detF))->Fill(frontHit->GetEnergy(), backHit->GetEnergy());

            if(frontHit->GetEnergy() > 4900)
                fH2[slot].at(Form("Det%dFECut", detF))->Fill(frontMaps[detF][stripF], frontHit->GetEnergy());

            if(backHit->GetEnergy() > 4900)
                fH2[slot].at(Form("Det%dBECut", detF))->Fill(backMaps[detF][stripB], backHit->GetEnergy());

            if(frontHit->GetEnergy() > 4900 && backHit->GetEnergy() > 4900)
                fH2[slot].at(Form("Det%dCut", detF))->Fill(frontMaps[detF][stripF], backMaps[detB][stripB]);
        }
    }
   */
}



void RCMPEnergyHelper::EndOfSort(std::shared_ptr<std::map<std::string, TList>> list)
{

}
