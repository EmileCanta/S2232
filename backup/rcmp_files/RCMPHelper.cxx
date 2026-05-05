#include "RCMPHelper.hh"

#include "ChannelMapping.h"

using namespace std;

void RCMPHelper::CreateHistograms(unsigned int slot)
{
    fH1[slot]["Energy"] = new TH1F("Energy", "Energy", 100000, 0, 100000);
    fH1[slot]["Time"] = new TH1F("Time", "Time", 60000, 0, 60000e9);
    fH1[slot]["Multiplicity"] = new TH1F("Multiplicity", "Multiplicity", 10, 0, 10);
    
    fH1[slot]["FrontOrBack"] = new TH1F("FrontOrBack", "FrontOrBack", 2, 0, 2);
    
    fH2[slot]["HitMap"] = new TH2F("HitMap", "HitMap", 34, -1, 33, 34, -1, 33);
    fH2[slot]["FrontHitCorrelation"] = new TH2F("FrontHitCorrelation", "FrontHitCorrelation", 34, -1, 33, 34, -1, 33);
    fH2[slot]["BackHitCorrelation"] = new TH2F("BackHitCorrelation", "BackHitCorrelation", 34, -1, 33, 34, -1, 33);
    fH2[slot]["SideVSMultiplicity"] = new TH2F("SideVSMultiplicity", "SideVSMultiplicity", 10, 0, 10, 2, 0, 2);
}

void RCMPHelper::Exec(unsigned int slot, TRcmp& rcmp)
{
    TRcmpHit* hit1 = nullptr;
    TRcmpHit* hit2 = nullptr;

    double mult = rcmp.GetMultiplicity();

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

        string side1 = hit1->GetChannel()->GetMnemonic()->CollectedChargeString();
        string side2 = hit2->GetChannel()->GetMnemonic()->CollectedChargeString();
        
        if((side1 != side2) && (det1 == 6) && (det1 == det2))
        {
            fH1[slot].at("Energy")->Fill(hit1->GetEnergy());
        }
        
        if((det1 == det2 && det1 == 6) && (side1 != side2))
        {
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

        if(det == 6)
        {
            if(side == "P") 
            {
                fH2[slot].at("SideVSMultiplicity")->Fill(mult, 0);

                //fH1[slot].at("Energy")->Fill(hit->GetEnergy());
                fH1[slot].at("Time")->Fill(hit->GetTime());
            }

            if(side == "N") 
            {
                fH2[slot].at("SideVSMultiplicity")->Fill(mult, 1);
            }

            if(mult > 0) fH1[slot].at("Multiplicity")->Fill(mult, 1./mult); // Careful, weighted histogram here
        }
    }
}

void RCMPHelper::EndOfSort(std::shared_ptr<std::map<std::string, TList>>& list)
{

}
