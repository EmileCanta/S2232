// Written by Dhruval Shah and Vinzenz Bildstein

#include "RCMPHelper.hh"
#include "TStyle.h" // to use SetOptStat
#include <limits>

// Coincidences Gates (nanoseconds)
double ggLow  = -300.;
double ggHigh = 300.;
double grLow  = -50.;
double grHigh = 50.;

// Griffin-Griffin gate - returns true if the time between two GRIFFIN hits is within the coincidence window

// RCMP variables
double xLowR  = 0.5;
double xHighR = 15000.5;
double xLowMul  = 0.5;
double xHighMul = 1100.5;
double detLow = -0.5;
double detHigh = 6.5;
double timeLow = -10000.5;
double timeHigh = 10000.5;

// GRIFFIN variables
double xLowG  = 0.5;
double xHighG = 16000.5;

// 2D Hist variables for GRIFFIN

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

   fH1[slot]["rcmpEF"] = new TH1F("rcmpEF", "RCMP Front Energy; Front Energy (keV); Counts per keV", xHighR - xLowR, xLowR, xHighR);
   fH1[slot]["rcmpEB"] = new TH1F("rcmpEB", "RCMP Back Energy; Back Energy (keV); Counts per keV", xHighR - xLowR, xLowR, xHighR);
   fH1[slot]["rcmpFB_Time"] = new TH1F("rcmpFB_Time", "RCMP Front Back Time difference; time(ns); Counts per ns", 400, -2000, 2000);

   fH2[slot]["frag"] = new TH2F("frag", "Fragments; Front; Back", detHigh - detLow, detLow, detHigh, detHigh - detLow, detLow, detHigh);
   fH2[slot]["rcmpEDet"] = new TH2F("rcmpEDet", "RCMP Pixel Energy and Detector Number; Energy (keV); Detector Number", xHighR - xLowR, xLowR, xHighR, detHigh - detLow, detLow, detHigh);
   //fH2[slot]["rcmpEFDet"] = new TH2F("rcmpEFDet", "RCMP Front Energy and Channel Number; Energy (keV); Detector Number", xHighR - xLowR, xLowR, xHighR, 32, -0.5, 31.5);
   //fH2[slot]["rcmpEBDet"] = new TH2F("rcmpEFDet", "RCMP Back Energy and Channel Number; Energy (keV); Detector Number", xHighR - xLowR, xLowR, xHighR, 32, -0.5, 31.5);
   fH2[slot]["rcmpFrontBack"] = new TH2F("rcmpFrontBack", "RCMP Front Back Hit Pattern;Front Channel; Back Channel", 32, 0.5, 32.5, 32, 0.5, 32.5);
   fH2[slot]["DetFE_Tagged"] = new TH2F("DetFE_Tagged", "RCMP Detector Front Hit Back tagged; Front; Front Energy" , 32, 0.5, 32.5, 10000, 0.5, 10000.5); 

    
   for(int i = 1; i <= 6; ++i) 
   {
	fH2[slot][Form("Det%d", i)] = new TH2F(Form("Det%d", i), Form("RCMP Detector %d; Front; Back", i), 33, -0.5, 32.5, 33, -0.5, 32.5);  

	fH2[slot][Form("Det%dFE", i)] = new TH2F(Form("Det%dFE", i), Form("RCMP Detector %d; Front; Front Energy", i), 32, 0.5, 32.5, 2000, 0.5, 2000.5); 
	fH2[slot][Form("Det%dBE", i)] = new TH2F(Form("Det%dBE", i), Form("RCMP Detector %d; Back; Front Energy", i), 32, 0.5, 32.5, 2000, 0.5, 2000.5); 
   }

   

 
}

// TODO: Change the function arguments to match the detectors you want to use and the declaration in the header file!
void RCMPHelper::Exec(unsigned int slot,TRcmp& rcmp)
{
   /********** RCMP **********/


   // Fill RCMP Front Multiplicity spectra
   for(int i = 0; i < rcmp.GetFrontMultiplicity(); ++i) 
   {
      auto rcmp1 = rcmp.GetRcmpFrontHit(i);

      fH1[slot].at("rcmpEF")->Fill(rcmp1->GetEnergy());
     // if(rcmp1->GetEnergy()>1000){
     //fill events for DSSSD1 
     //if(rcmp1->GetDetector()+1>0 && rcmp1->GetDetector()+1<=32){
     //std::cout<<rcmp1->GetDetector()+1<<std::endl;
     // if(rcmp1->GetDetector()==1){
      //std::cout<<rcmp1->GetSegment()+1<<std::endl;
      //fH2[slot].at("Det1FE")->Fill(rcmp1->GetSegment()+1, rcmp1->GetEnergy());}
      

	
	fH2[slot].at(Form("Det%dFE", rcmp1->GetDetector()))->Fill(rcmp1->GetSegment()+1, rcmp1->GetEnergy());

   
      
      
      
      
      
   
	for(int j = 0; j < rcmp.GetBackMultiplicity(); ++j) 
   {
      auto rcmp2 = rcmp.GetRcmpBackHit(j);
      double time_diff = rcmp1->GetTimeStampNs()-rcmp2->GetTimeStampNs();

    
      fH1[slot].at("rcmpFB_Time")->Fill(time_diff);
      if(time_diff<20){
      fH2[slot].at(Form("Det%d", rcmp2->GetDetector()))->Fill(rcmp1->GetSegment()+1, rcmp2->GetSegment()+1);}	 
      fH1[slot].at("rcmpEB")->Fill(rcmp2->GetEnergy());

      //std::cout<<rcmp2->GetSegment()+1<<std::endl;
      fH2[slot].at(Form("Det%dBE", rcmp2->GetDetector()))->Fill(rcmp2->GetSegment()+1, rcmp2->GetEnergy());

   
      
 
   }}
   

   
   
   

   
  }
   
  


void RCMPHelper::EndOfSort(std::shared_ptr<std::map<std::string, TList>> list)
{
}




