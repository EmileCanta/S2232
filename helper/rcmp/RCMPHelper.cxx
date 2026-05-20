#include "RCMPHelper.hh"

#include "ChannelMapping.h"

using namespace std;

void RCMPHelper::CreateHistograms(unsigned int slot)
{
    fH1[slot]["Energy"] = new TH1F("Energy", "Energy", 10000, 0, 10000);
    fH1[slot]["EnergyGriffin"] = new TH1F("EnergyGriffin", "EnergyGriffin", 10000, 0, 10000);
    fH1[slot]["Time"] = new TH1F("Time", "Time", 60000, 0, 60000e9);
    fH1[slot]["Multiplicity"] = new TH1F("Multiplicity", "Multiplicity", 100, 0, 100);

    fH2[slot]["HitMap"] = new TH2F("HitMap", "HitMap", 34, -1, 33, 34, -1, 33);
    fH2[slot]["FrontHitCorrelation"] = new TH2F("FrontHitCorrelation", "FrontHitCorrelation", 34, -1, 33, 34, -1, 33);
    fH2[slot]["BackHitCorrelation"] = new TH2F("BackHitCorrelation", "BackHitCorrelation", 34, -1, 33, 34, -1, 33);
    fH2[slot]["SideVSMultiplicity"] = new TH2F("SideVSMultiplicity", "SideVSMultiplicity", 100, 0, 100, 2, 0, 2);
    fH2[slot]["DetectorVSMultiplicity"] = new TH2F("DetectorVSMultiplicity", "DetectorVSMultiplicity", 8, 0, 8, 100, 0, 100);
    fH2[slot]["EnergyVSFrontStrip1"] = new TH2F("EnergyVSFrontStrip1", "EnergyVSFrontStrip1", 36, -2, 34, 10000, 0, 10000);
    fH2[slot]["EnergyVSFrontStrip2"] = new TH2F("EnergyVSFrontStrip2", "EnergyVSFrontStrip2", 36, -2, 34, 10000, 0, 10000);
    fH2[slot]["EnergyVSFrontStrip3"] = new TH2F("EnergyVSFrontStrip3", "EnergyVSFrontStrip3", 36, -2, 34, 10000, 0, 10000);
    fH2[slot]["EnergyVSFrontStrip4"] = new TH2F("EnergyVSFrontStrip4", "EnergyVSFrontStrip4", 36, -2, 34, 10000, 0, 10000);
    fH2[slot]["EnergyVSFrontStrip5"] = new TH2F("EnergyVSFrontStrip5", "EnergyVSFrontStrip5", 36, -2, 34, 10000, 0, 10000);
    fH2[slot]["EnergyVSFrontStrip6"] = new TH2F("EnergyVSFrontStrip6", "EnergyVSFrontStrip6", 36, -2, 34, 10000, 0, 10000);
    fH2[slot]["EnergyVSBackStrip1"] = new TH2F("EnergyVSBackStrip1", "EnergyVSBackStrip1", 36, -2, 34, 10000, 0, 10000);
    fH2[slot]["EnergyVSBackStrip2"] = new TH2F("EnergyVSBackStrip2", "EnergyVSBackStrip2", 36, -2, 34, 10000, 0, 10000);
    fH2[slot]["EnergyVSBackStrip3"] = new TH2F("EnergyVSBackStrip3", "EnergyVSBackStrip3", 36, -2, 34, 10000, 0, 10000);
    fH2[slot]["EnergyVSBackStrip4"] = new TH2F("EnergyVSBackStrip4", "EnergyVSBackStrip4", 36, -2, 34, 10000, 0, 10000);
    fH2[slot]["EnergyVSBackStrip5"] = new TH2F("EnergyVSBackStrip5", "EnergyVSBackStrip5", 36, -2, 34, 10000, 0, 10000);
    fH2[slot]["EnergyVSBackStrip6"] = new TH2F("EnergyVSBackStrip6", "EnergyVSBackStrip6", 36, -2, 34, 10000, 0, 10000);
    fH2[slot]["FrontVSBackEnergy"] = new TH2F("FrontVSBackEnergy", "FrontVSBackEnergy", 1000, 0, 10000, 1000, 0, 10000);
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

        else if(rcmp.GetRcmpHit(0)->GetChannel()->GetMnemonic()->CollectedChargeString() == "N")
        {
            hit1 = rcmp.GetRcmpHit(1);
            hit2 = rcmp.GetRcmpHit(0);
        }

        int det1 = hit1->GetDetector();
        int det2 = hit2->GetDetector();

        int strip1 = hit1->GetSegment();
        int strip2 = hit2->GetSegment();

        string side1 = hit1->GetChannel()->GetMnemonic()->CollectedChargeString();
        string side2 = hit2->GetChannel()->GetMnemonic()->CollectedChargeString();

        for(int ndet = 1; ndet <= 6; ndet++)
        {
            string nameFront = Form("EnergyVSFrontStrip%d", ndet);
            string nameBack = Form("EnergyVSBackStrip%d", ndet);
            
            if((det1 == ndet) && (det1 == det2) && (side1 != side2))
            {
                fH2[slot].at(nameFront)->Fill(strip1, hit1->GetEnergy());
                fH2[slot].at(nameBack)->Fill(strip2, hit2->GetEnergy());
            }
        }

        if((det1 == 6) && (det1 == det2) && (side1 != side2))
        {
            fH1[slot].at("Time")->Fill(hit1->GetTime());
            fH2[slot].at("FrontVSBackEnergy")->Fill(hit1->GetEnergy(), hit2->GetEnergy());
            fH2[slot].at("HitMap")->Fill(frontMaps[det1][hit1->GetSegment()], backMaps[det2][hit2->GetSegment()]);
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

    for(int i = 0; i < mult; i++)
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

    /* for(int i = 0; i < multGriffin; i++)
       {
       TGriffinHit* hit = griffin.GetGriffinHit(i);
       fH1[slot].at("EnergyGriffin")->Fill(hit->GetEnergy());
       } */
}

void RCMPHelper::EndOfSort(std::shared_ptr<std::map<std::string, TList>>& list)
{

}
