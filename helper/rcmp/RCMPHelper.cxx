#include "RCMPHelper.hh"

#include "ChannelMapping.h"

using namespace std;

void RCMPHelper::CreateHistograms(unsigned int slot)
{
    fH1[slot]["Energy"] = new TH1F("Energy", "Energy", 10000, 0, 10000);
    fH1[slot]["EnergyStripGatedCenterStrip"] = new TH1F("EnergyStripGatedCenterStrip", "EnergyStripGatedCenterStrip", 10000, 0, 10000);
    fH1[slot]["ChargeStripGatedCenterStrip"] = new TH1F("ChargeStripGatedCenterStrip", "ChargeStripGatedCenterStrip", 100000, 0, 100000);
    fH1[slot]["CalEnergyStripGatedCenterStrip"] = new TH1F("CalEnergyStripGatedCenterStrip", "CalEnergyStripGatedCenterStrip", 10000, 0, 10000);
    fH1[slot]["EnergyStripGatedOuterStrip"] = new TH1F("EnergyStripGatedOuterStrip", "EnergyStripGatedOuterStrip", 10000, 0, 10000);
    fH1[slot]["EnergyGriffin"] = new TH1F("EnergyGriffin", "EnergyGriffin", 10000, 0, 10000);
    fH1[slot]["Time"] = new TH1F("Time", "Time", 60000, 0, 60000e9);
    fH1[slot]["Multiplicity"] = new TH1F("Multiplicity", "Multiplicity", 100, 0, 100);

    fH2[slot]["FrontHitCorrelation"] = new TH2F("FrontHitCorrelation", "FrontHitCorrelation", 34, -1, 33, 34, -1, 33);
    fH2[slot]["BackHitCorrelation"] = new TH2F("BackHitCorrelation", "BackHitCorrelation", 34, -1, 33, 34, -1, 33);
    fH2[slot]["SideVSMultiplicity"] = new TH2F("SideVSMultiplicity", "SideVSMultiplicity", 100, 0, 100, 2, 0, 2);
    fH2[slot]["DetectorVSMultiplicity"] = new TH2F("DetectorVSMultiplicity", "DetectorVSMultiplicity", 8, 0, 8, 100, 0, 100);
    fH2[slot]["FrontVSBackEnergy"] = new TH2F("FrontVSBackEnergy", "FrontVSBackEnergy", 1000, 0, 10000, 1000, 0, 10000);

    for(int ndet = 1; ndet <= 6; ndet++)
    {
        fH2[slot][Form("EnergyVSFrontStrip%d", ndet)] = new TH2F(Form("EnergyVSFrontStrip%d", ndet), Form("EnergyVSFrontStrip%d", ndet), 36, -2, 34, 10000, 0, 10000);
        fH2[slot][Form("EnergyVSBackStrip%d", ndet)] = new TH2F(Form("EnergyVSBackStrip%d", ndet), Form("EnergyVSBackStrip%d", ndet), 36, -2, 34, 10000, 0, 10000);
        
        //Careful binning 1 bin = 10 keV
        fH2[slot][Form("ChargeVSFrontStrip%d", ndet)] = new TH2F(Form("ChargeVSFrontStrip%d", ndet), Form("ChargeVSFrontStrip%d", ndet), 36, -2, 34, 1000, 0, 10000);
        fH2[slot][Form("ChargeVSBackStrip%d", ndet)] = new TH2F(Form("ChargeVSBackStrip%d", ndet), Form("ChargeVSBackStrip%d", ndet), 36, -2, 34, 1000, 0, 10000);
        fH2[slot][Form("ChargeVSPixelFront%d", ndet)] = new TH2F(Form("ChargeVSPixelFront%d", ndet), Form("ChargeVSPixelFront%d", ndet), 1028, -2, 1026, 1000, 0, 10000);
        fH2[slot][Form("ChargeVSPixelBack%d", ndet)] = new TH2F(Form("ChargeVSPixelBack%d", ndet), Form("ChargeVSPixelBack%d", ndet), 1028, -2, 1026, 1000, 0, 10000);

        fH2[slot][Form("HitMapNonCorrected%d", ndet)] = new TH2F(Form("HitMapNonCorrected%d", ndet), Form("HitMapNonCorrected%d", ndet), 34, -1, 33, 34, -1, 33);
        fH2[slot][Form("HitMapCorrectedGood%d", ndet)] = new TH2F(Form("HitMapCorrectedGood%d", ndet), Form("HitMapCorrectedGood%d", ndet), 34, -1, 33, 34, -1, 33);
        fH2[slot][Form("HitMapCorrectedWeird%d", ndet)] = new TH2F(Form("HitMapCorrectedWeird%d", ndet), Form("HitMapCorrectedWeird%d", ndet), 34, -1, 33, 34, -1, 33);
    }
}

void RCMPHelper::Exec(unsigned int slot, TRcmp& rcmp, TGriffin& griffin, TGriffinBgo& griffinbgo)
{
    TRcmpHit* hit1;
    TRcmpHit* hit2;

    //double mult = rcmp.GetRcmpMultiplicity();
    double mult = rcmp.GetMultiplicity();
    double multGriffin = griffin.GetMultiplicity();

    if(mult == 2) 
    {
        if(rcmp.GetRcmpHit(0)->GetChannel()->GetMnemonic()->CollectedChargeString() == "P")
        {
            hit1 = rcmp.GetRcmpHit(0);
            hit2 = rcmp.GetRcmpHit(1);
        }

        if(rcmp.GetRcmpHit(0)->GetChannel()->GetMnemonic()->CollectedChargeString() == "N")
        {
            hit1 = rcmp.GetRcmpHit(1);
            hit2 = rcmp.GetRcmpHit(0);
        }

        int det1 = hit1->GetDetector();
        int det2 = hit2->GetDetector();

        int strip1 = hit1->GetSegment();
        int strip2 = hit2->GetSegment();

        int mappedstrip1 = frontMaps[det1][hit1->GetSegment()];
        int mappedstrip2 = backMaps[det2][hit2->GetSegment()];

        int mappedstrip1bis = frontMapsBis[det1][hit1->GetSegment()];
        int mappedstrip2bis = backMapsBis[det2][hit2->GetSegment()];

        string side1 = hit1->GetChannel()->GetMnemonic()->CollectedChargeString();
        string side2 = hit2->GetChannel()->GetMnemonic()->CollectedChargeString();

        int pixel = mappedstrip1 + 32 * mappedstrip2;

        fH1[slot].at("Energy")->Fill(hit1->GetEnergy());

        for(int ndet = 1; ndet <= 6; ndet++)
        {
            if((det1 == ndet) && (det1 == det2) && (side1 != side2))
            {
                fH2[slot].at(Form("EnergyVSFrontStrip%d", ndet))->Fill(mappedstrip1, hit1->GetEnergy());
                fH2[slot].at(Form("EnergyVSBackStrip%d", ndet))->Fill(mappedstrip2, hit2->GetEnergy());
                fH2[slot].at(Form("ChargeVSFrontStrip%d", ndet))->Fill(mappedstrip1, hit1->GetCharge());
                fH2[slot].at(Form("ChargeVSBackStrip%d", ndet))->Fill(mappedstrip2, hit2->GetCharge());

                fH2[slot].at(Form("HitMapNonCorrected%d", ndet))->Fill(strip1, strip2);
                fH2[slot].at(Form("HitMapCorrectedGood%d", ndet))->Fill(mappedstrip1, mappedstrip2);
                fH2[slot].at(Form("HitMapCorrectedWeird%d", ndet))->Fill(mappedstrip1bis, mappedstrip2bis);
                
                fH2[slot].at(Form("ChargeVSPixelFront%d", ndet))->Fill(pixel, hit1->GetCharge());
                fH2[slot].at(Form("ChargeVSPixelBack%d", ndet))->Fill(pixel, hit2->GetCharge());
            }
        }

        if((det1 == det2) && (side1 != side2))
        {
            if((det1 == 6) && (mappedstrip1 == 30) && (mappedstrip2 == 30))
            {
                fH1[slot].at("Time")->Fill(hit1->GetTime());
                
                fH2[slot].at("FrontVSBackEnergy")->Fill(hit1->GetEnergy(), hit2->GetEnergy());
                
                fH1[slot].at("EnergyStripGatedOuterStrip")->Fill(hit1->GetEnergy());
            }

            if((det1 == 6) && (mappedstrip1 == 16) && (mappedstrip2 == 16))
            {
                fH1[slot].at("EnergyStripGatedCenterStrip")->Fill(hit1->GetEnergy());
                fH1[slot].at("ChargeStripGatedCenterStrip")->Fill(hit1->GetCharge());
                fH1[slot].at("CalEnergyStripGatedCenterStrip")->Fill(hit1->GetCharge() * 1.354 + 191.083);
            }
        }

        if((det1 == det2) && (side1 == side2 && side1 == "P"))
        {
            fH2[slot].at("FrontHitCorrelation")->Fill(frontMaps[det1][hit1->GetSegment()], frontMaps[det2][hit2->GetSegment()]);
        }

        if((det1 == det2) && (side1 == side2 && side1 == "N"))
        {
            fH2[slot].at("BackHitCorrelation")->Fill(backMaps[det1][hit1->GetSegment()], backMaps[det2][hit2->GetSegment()]);
        }
    }

    /* for(int i = 0; i < mult; i++)
       {
       TRcmpHit* hit = rcmp.GetRcmpHit(i);

       int det = hit->GetDetector();
       string side = hit->GetChannel()->GetMnemonic()->CollectedChargeString();

       fH2[slot].at("DetectorVSMultiplicity")->Fill(det, mult, 1./mult); // Careful, weighted histogram here

       if(det == 6)
       {
       if(side == "P") 
       {
       fH2[slot].at("SideVSMultiplicity")->Fill(mult, 0);
       }

       if(side == "N") 
       {
       fH2[slot].at("SideVSMultiplicity")->Fill(mult, 1);
       }

       if(mult > 0) fH1[slot].at("Multiplicity")->Fill(mult, 1./mult); // Careful, weighted histogram here
       }
       }

       for(int i = 0; i < multGriffin; i++)
       {
       TGriffinHit* hit = griffin.GetGriffinHit(i);
       fH1[slot].at("EnergyGriffin")->Fill(hit->GetEnergy());
       } */ 
}

void RCMPHelper::EndOfSort(std::shared_ptr<std::map<std::string, TList>>& list)
{

}
