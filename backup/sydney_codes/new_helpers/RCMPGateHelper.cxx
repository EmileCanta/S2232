// Written by Sydney Plante, 2025

#include "RCMPGateHelper.hh"
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


// RCMP variables. !!The 0.5 is to make the bins centered!!
double xLowR  = 0; 
double xHighR = 10000; 
double xLowMul  = 0.5; // Number of hits the detector detects 
double xHighMul = 1100.5; 
double detLow = 0.5; // Number of sides of the cube! originally set at -0.5 - 6.5 giving 7 total "detector bins", shifted it to 0.5 to only have 6.
double detHigh = 6.5; 
double timeLow = -10000.5;
double timeHigh = 10000.5;
double frontStripLow = -0.5; 
double frontStripHigh = 192.5; // 32 x 6 = 192
int stripBins = 33;
double stripLow = -0.5; // Number of strips per side
double stripHigh = 32.5; 
double energyBin = 6000;
double energyLow = 0;
double energyHigh = 6000;


// GRIFFIN variables
double xLowG  = -0.5; //Same as the energy limits for RCMP
double xHighG = 10000.5;

std::vector<int> energies = {238, 275, 984, 1634};

const double decayStart = 0.0;
const double decayEnd = 55.0;
const double daStart = 56.0;
const double daEnd = 60.0;



// 2D Hist variables for GRIFFIN

void RCMPGateHelper::CreateHistograms(unsigned int slot)
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

   // fH2[slot]["RCMPvsGrif"] = new TH2F("RCMPvsGrif", "RCMP vs GRIFFIN Energy; RCMP Energy, GRIFFIN Energy", 6000, 0, 6000, 6000, 0, 6000);
   // fH2[slot]["GrifvsGrif"] = new TH2F("GrifvsGrif", "GRIFFIN vs GRIFFIN Energy; GRIFFIN Energy, GRIFFIN Energy", 6000, 0, 6000, 6000, 0, 6000);
   fH1[slot]["BetaGate"] = new TH1F("BetaGate", "GRIFFIN Energy Gated on Beta Spectrum; GRIFFIN Energy [keV], Counts", 10000, 0, 10000);

   for (int energy : energies) {
        fH1[slot][Form("gate%dD56", energy)] = new TH1F(Form("gate%dD56", energy), Form("%d Peak Gate, Detectors 5 and 6; Energy (keV); Count", energy), (xHighR - xLowR)/2, xLowR, xHighR);
        fH1[slot][Form("gate%dLD56", energy)] = new TH1F(Form("gate%dLD56", energy), Form("%d Peak Gate Left, Detectors 5 and 6; Energy (keV); Count", energy), (xHighR - xLowR)/2, xLowR, xHighR);
        fH1[slot][Form("gate%dRD56", energy)] = new TH1F(Form("gate%dRD56", energy), Form("%d Peak Gate Right, Detectors 5 and 6; Energy (keV); Count", energy), (xHighR - xLowR)/2, xLowR, xHighR);
   }

}

void RCMPGateHelper::Exec(unsigned int slot, TGriffin& grif, TGriffinBgo& grifBgo, TRcmp& rcmp)
{

   for (int i = 0; i < rcmp.GetFrontMultiplicity(); ++i)
   {
      auto rcmp1 = rcmp.GetRcmpFrontHit(i);
      int detF = rcmp1 -> GetDetector();
      int stripF = rcmp1 -> GetStripNumber();

      if (stripF < 2 || stripF > 29) 
      {
         continue;
      }

      for (int j = 0; j < rcmp.GetBackMultiplicity(); ++j) 
      {
         auto rcmp2 = rcmp.GetRcmpBackHit(j);
         int detB = rcmp2 -> GetDetector();
         int stripB = rcmp2 -> GetStripNumber();

         if (stripB < 2 || stripB > 30 || detF != detB || (detF != 5 && detF !=6) || !PromptCoincidence(rcmp1,rcmp2))
         {
            continue;
         }

         double frontEnergy = rcmp1->GetEnergy();

         for (int k = 0; k < grif.GetMultiplicity(); ++k)
         {

            auto grif1 = grif.GetGriffinHit(k);
            double grifEnergy = grif1 -> GetEnergy();

            if (!PromptCoincidence(grif1, rcmp1))
            {
               continue;
            }

            if (frontEnergy <= 500) {
               fH1[slot].at("BetaGate") -> Fill(grifEnergy);
            }


            for (int energy : energies) 
            {
               // Middle Gate //
               if ((grifEnergy >= energy -3) && (grifEnergy <= energy +3))
               {
                  fH1[slot].at(Form("gate%dD56", energy))->Fill(frontEnergy);
               }
               
               // Left Gate //
               if ((grifEnergy >= energy -7) && (grifEnergy <= energy -4)) 
               {
                  fH1[slot].at(Form("gate%dLD56", energy))->Fill(frontEnergy);
               }
               
               // Right Gate // 
               if ((grifEnergy >= energy + 4) && (grifEnergy <= energy + 7)) 
               {
                  fH1[slot].at(Form("gate%dRD56", energy))->Fill(frontEnergy);
               }

            }
         }


      }

   }
   

}



void RCMPGateHelper::EndOfSort(std::shared_ptr<std::map<std::string, TList>> list)
{

}

