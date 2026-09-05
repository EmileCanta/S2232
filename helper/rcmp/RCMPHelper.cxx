#include "RCMPHelper.hh"

#include "ChannelMapping.h"

using namespace std;

void RCMPHelper::CreateHistograms(unsigned int slot)
{
    fH1[slot]["EnergyAll_MultTwo"] = new TH1F("EnergyAll_MultTwo", "EnergyAll_MultTwo", 10000, 0, 10000);
    fH1[slot]["TimeAll_MultTwo"] = new TH1F("TimeAll_MultTwo", "TimeAll_MultTwo", 6000, 0, 60e9); //1 bin = 0.01e9 ns = 10 ms

    fH1[slot]["EnergyGriffin_Singles"] = new TH1F("EnergyGriffin_Singles", "EnergyGriffin_Singles", 10000, 0, 10000);
    fH1[slot]["TimeGriffin_Singles"] = new TH1F("TimeGriffin_Singles", "TimeGriffin_Singles", 6000, 0, 60e9); //1 bin = 0.01e9 ns = 10 ms
    fH1[slot]["AbsTimeGriffin_Singles"] = new TH1F("AbsTimeGriffin_Singles", "AbsTimeGriffin_Singles", 600000, 0, 6000e9); //1 bin = 0.01e9 ns = 10 ms

    fH2[slot]["EnergyFront1VSTimeDiffFrontFront_Det5Or6"] = new TH2F("EnergyFront1VSTimeDiffFrontFront_Det5Or6", "EnergyFront1VSTimeDiffFrontFront_Det5Or6", 500, -2500, 2500, 1000, 0, 10000); //Less bins

    fH2[slot]["FrontHitCorrelation_MultTwo"] = new TH2F("FrontHitCorrelation_MultTwo", "FrontHitCorrelation_MultTwo", 32, 0, 32, 32, 0, 32);
    fH2[slot]["BackHitCorrelation_MultTwo"] = new TH2F("BackHitCorrelation_MultTwo", "BackHitCorrelation_MultTwo", 32, 0, 32, 32, 0, 32);
    fH2[slot]["DetRepartition_MultTwo"] = new TH2F("DetRepartition_MultTwo", "DetRepartition_MultTwo", 6, 1, 7, 6, 1, 7);
   
    fH2[slot]["DetectorVSMultiplicity_Singles"] = new TH2F("DetectorVSMultiplicity_Singles", "DetectorVSMultiplicity_Singles", 100, 0, 100, 6, 1, 7);

    fH2[slot]["EnergyFrontVSTimeFront_CondGRIF"] = new TH2F("EnergyFrontVSTimeFront_CondGRIF", "EnergyFrontVSTimeFront_CondGRIF", 6000, 0, 60e9, 1000, 0, 10000); //1 bin = 0.01e9 ns = 10 ms //Less bins
    fH2[slot]["EnergyFrontVSTimeFront_Singles"] = new TH2F("EnergyFrontVSTimeFront_Singles", "EnergyFrontVSTimeFront_Singles", 6000, 0, 60e9, 1000, 0, 10000); //1 bin = 0.01e9 ns = 10 ms //Less bins
    fH2[slot]["EnergyGriffinVSTimeGriffin_CondRCMP"] = new TH2F("EnergyGriffinVSTimeGriffin_CondRCMP", "EnergyGriffinVSTimeGriffin_CondRCMP", 6000, 0, 60e9, 10000, 0, 10000); //1 bin = 0.01e9 ns = 10 ms
    fH2[slot]["EnergyGriffinVSTimeGriffin_Singles"] = new TH2F("EnergyGriffinVSTimeGriffin_Singles", "EnergyGriffinVSTimeGriffin_Singles", 6000, 0, 60e9, 10000, 0, 10000); //1 bin = 0.01e9 ns = 10 ms

    for(int ndet = 1; ndet <= 6; ndet++)
    {
        fH2[slot][Form("EnergyVSFrontStrip_MultTwo%d", ndet)] = new TH2F(Form("EnergyVSFrontStrip_MultTwo%d", ndet), Form("EnergyVSFrontStrip_MultTwo%d", ndet), 32, 0, 32, 10000, 0, 10000);
        fH2[slot][Form("EnergyVSBackStrip_MultTwo%d", ndet)] = new TH2F(Form("EnergyVSBackStrip_MultTwo%d", ndet), Form("EnergyVSBackStrip_MultTwo%d", ndet), 32, 0, 32, 10000, 0, 10000);

        fH2[slot][Form("ChargeVSFrontStrip_MultTwo%d", ndet)] = new TH2F(Form("ChargeVSFrontStrip_MultTwo%d", ndet), Form("ChargeVSFrontStrip_MultTwo%d", ndet), 32, 0, 32, 2500, 0, 10000); //Less bins
        fH2[slot][Form("ChargeVSBackStrip_MultTwo%d", ndet)] = new TH2F(Form("ChargeVSBackStrip_MultTwo%d", ndet), Form("ChargeVSBackStrip_MultTwo%d", ndet), 32, 0, 32, 2500, 0, 10000); //Less bins
        
        fH2[slot][Form("ChargeVSPixelFront_MultTwo%d", ndet)] = new TH2F(Form("ChargeVSPixelFront_MultTwo%d", ndet), Form("ChargeVSPixelFront_MultTwo%d", ndet), 1023, 0, 1023, 1000, 0, 10000); //Less bins
        fH2[slot][Form("ChargeVSPixelBack_MultTwo%d", ndet)] = new TH2F(Form("ChargeVSPixelBack_MultTwo%d", ndet), Form("ChargeVSPixelBack_MultTwo%d", ndet), 1023, 0, 1023, 1000, 0, 10000); //Less bins

        fH2[slot][Form("HitMapNonCorrected_MultTwo%d", ndet)] = new TH2F(Form("HitMapNonCorrected_MultTwo%d", ndet), Form("HitMapNonCorrected_MultTwo%d", ndet), 32, 0, 32, 32, 0, 32);
        fH2[slot][Form("HitMapCorrectedGood_MultTwo%d", ndet)] = new TH2F(Form("HitMapCorrectedGood_MultTwo%d", ndet), Form("HitMapCorrectedGood_MultTwo%d", ndet), 32, 0, 32, 32, 0, 32);
        fH2[slot][Form("HitMapCorrectedWeird_MultTwo%d", ndet)] = new TH2F(Form("HitMapCorrectedWeird_MultTwo%d", ndet), Form("HitMapCorrectedWeird_MultTwo%d", ndet), 32, 0, 32, 32, 0, 32);

        fH2[slot][Form("HitMapCorrectedWeird_MultTwo%d", ndet)] = new TH2F(Form("HitMapCorrectedWeird_MultTwo%d", ndet), Form("HitMapCorrectedWeird_MultTwo%d", ndet), 32, 0, 32, 32, 0, 32);

        fH1[slot][Form("TimeDiffBackFront_MultTwo%d", ndet)] = new TH1F(Form("TimeDiffBackFront_MultTwo%d", ndet), Form("TimeDiffBackFront_MultTwo%d", ndet), 500, -2500, 2500);

        fH2[slot][Form("EnergyVSFrontStrip_Singles%d", ndet)] = new TH2F(Form("EnergyVSFrontStrip_Singles%d", ndet), Form("EnergyVSFrontStrip_Singles%d", ndet), 32, 0, 32, 10000, 0, 10000);
        fH2[slot][Form("EnergyVSBackStrip_Singles%d", ndet)] = new TH2F(Form("EnergyVSBackStrip_Singles%d", ndet), Form("EnergyVSBackStrip_Singles%d", ndet), 32, 0, 32, 10000, 0, 10000);

        fH2[slot][Form("EnergyFrontVSMult_Singles%d", ndet)] = new TH2F(Form("EnergyFrontVSMult_Singles%d", ndet), Form("EnergyFrontVSMult_Singles%d", ndet), 20, 0, 20, 1000, 0, 10000); //Less bins

        fH2[slot][Form("EnergyFrontVSEnergyGriffin_CondRCMP%d", ndet)] = new TH2F(Form("EnergyFrontVSEnergyGriffin_CondRCMP%d", ndet), Form("EnergyFrontVSEnergyGriffin_CondRCMP%d", ndet), 10000, 0, 10000, 1000, 0, 10000); //Less bins

        fH2[slot][Form("EnergyVSTimeDiffFrontGriffin_CondGRIF%d", ndet)] = new TH2F(Form("EnergyVSTimeDiffFrontGriffin_CondGRIF%d", ndet), Form("EnergyVSTimeDiffFrontGriffin_CondGRIF%d", ndet), 500, -2500, 2500, 1000, 0, 10000); //Less bins
        fH2[slot][Form("EnergyVSTimeDiffBackGriffin_CondGRIF%d", ndet)] = new TH2F(Form("EnergyVSTimeDiffBackGriffin_CondGRIF%d", ndet), Form("EnergyVSTimeDiffBackGriffin_CondGRIF%d", ndet), 500, -2500, 2500, 1000, 0, 10000); //Less bins

        fH2[slot][Form("EnergyBackVSEnergyFront_Gated%d", ndet)] = new TH2F(Form("EnergyBackVSEnergyFront_Gated%d", ndet), Form("EnergyBackVSEnergyFront_Gated%d", ndet), 1000, 0, 10000, 1000, 0, 10000); //Less bins
        
        fH1[slot][Form("TimeDiffBackFront_Gated%d", ndet)] = new TH1F(Form("TimeDiffBackFront_Gated%d", ndet), Form("TimeDiffBackFront_Gated%d", ndet), 500, -2500, 2500);
    }

    for(int multi = 1; multi <= 10; multi++)
    {
        fH2[slot][Form("EnergyFront1VSEnergyFront2Mult%d_Det5Or6", multi)] = new TH2F(Form("EnergyFront1VSEnergyFront2Mult%d_Det5Or6", multi), Form("EnergyFront1VSEnergyFront2Mult%d_Det5Or6", multi), 1000, 0, 10000, 1000, 0, 10000); //Less bins
        fH2[slot][Form("EnergyFront1VSEnergyFront2BgdMult%d_Det5Or6", multi)] = new TH2F(Form("EnergyFront1VSEnergyFront2BgdMult%d_Det5Or6", multi), Form("EnergyFront1VSEnergyFront2BgdMult%d_Det5Or6", multi), 1000, 0, 10000, 1000, 0, 10000); //Less bins
    }
}

void RCMPHelper::Exec(unsigned int slot, TRcmp& rcmp, TGriffin& griffin, TGriffinBgo& griffinbgo)
{
    TRcmpHit* hit1;
    TRcmpHit* hit2;

    int mult = rcmp.GetMultiplicity();
    int multGriffin = griffin.GetMultiplicity();

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

        if((det1 == det2) && (side1 != side2)) 
        {
            fH1[slot].at("EnergyAll_MultTwo")->Fill(hit1->GetEnergy());
            fH1[slot].at("TimeAll_MultTwo")->Fill(hit1->GetTimeStampNs()%60000000000LL);
        }

        if((det1 == det2) && (side1 != side2))
        {
            fH2[slot].at(Form("EnergyVSFrontStrip_MultTwo%d", det1))->Fill(mappedstrip1, hit1->GetEnergy());
            fH2[slot].at(Form("EnergyVSBackStrip_MultTwo%d", det1))->Fill(mappedstrip2, hit2->GetEnergy());

            fH2[slot].at(Form("ChargeVSFrontStrip_MultTwo%d", det1))->Fill(strip1, hit1->GetCharge()); //Careful to not use mapping for calibration code to work correctly!
            fH2[slot].at(Form("ChargeVSBackStrip_MultTwo%d", det1))->Fill(strip2, hit2->GetCharge());

            fH2[slot].at(Form("HitMapNonCorrected_MultTwo%d", det1))->Fill(strip1, strip2);
            fH2[slot].at(Form("HitMapCorrectedWeird_MultTwo%d", det1))->Fill(mappedstrip1bis, mappedstrip2bis); //This is for files calibrated with ODB previous to RUN28254

            fH2[slot].at(Form("ChargeVSPixelFront_MultTwo%d", det1))->Fill(pixel, hit1->GetCharge());
            fH2[slot].at(Form("ChargeVSPixelBack_MultTwo%d", det1))->Fill(pixel, hit2->GetCharge());

            if(det1 == 1) fH2[slot].at(Form("HitMapCorrectedGood_MultTwo%d", det1))->Fill(mappedstrip1, mappedstrip2); 
            if(det1 == 2) fH2[slot].at(Form("HitMapCorrectedGood_MultTwo%d", det1))->Fill(mappedstrip2, mappedstrip1);

            if(det1 == 3) fH2[slot].at(Form("HitMapCorrectedGood_MultTwo%d", det1))->Fill(mappedstrip1, mappedstrip2);
            if(det1 == 4) fH2[slot].at(Form("HitMapCorrectedGood_MultTwo%d", det1))->Fill(mappedstrip2, mappedstrip1);

            if(det1 == 5) fH2[slot].at(Form("HitMapCorrectedGood_MultTwo%d", det1))->Fill(mappedstrip1, mappedstrip2);
            if(det1 == 6) fH2[slot].at(Form("HitMapCorrectedGood_MultTwo%d", det1))->Fill(mappedstrip2, mappedstrip1);

            fH1[slot].at(Form("TimeDiffBackFront_MultTwo%d", det1))->Fill(hit2->GetTimeStampNs() - hit1->GetTimeStampNs());
        }

        if((det1 == det2) && (side1 == side2) && (side1 == "P"))
        {
            fH2[slot].at("FrontHitCorrelation_MultTwo")->Fill(frontMaps[det1][hit1->GetSegment()], frontMaps[det2][hit2->GetSegment()]);
        }

        if((det1 == det2) && (side1 == side2) && (side1 == "N"))
        {
            fH2[slot].at("BackHitCorrelation_MultTwo")->Fill(backMaps[det1][hit1->GetSegment()], backMaps[det2][hit2->GetSegment()]);
        }
        
        fH2[slot].at("DetRepartition_MultTwo")->Fill(det1, det2);
    }

    //cout << "MULT IS: " << mult << endl;

    for(int i = 0; i < mult; i++)
    {
        TRcmpHit* hit = rcmp.GetRcmpHit(i);

        int det = hit->GetDetector();
        
        double energy = hit->GetEnergy();
        double time = hit->GetTimeStampNs()%60000000000LL;

        string side = hit->GetChannel()->GetMnemonic()->CollectedChargeString();

        fH2[slot].at("DetectorVSMultiplicity_Singles")->Fill(mult, det, 1./((double)mult)); // Careful, weighted histogram here

        if(side == "P")
        {
            fH2[slot].at(Form("EnergyVSFrontStrip_Singles%d", det))->Fill(frontMaps[det][hit->GetSegment()], energy);
            fH2[slot].at(Form("EnergyFrontVSMult_Singles%d", det))->Fill(mult, energy);
            if(det != 3 && det != 4) fH2[slot].at("EnergyFrontVSTimeFront_Singles")->Fill(time, energy);
        }

        if(side == "N") 
        {
            fH2[slot].at(Form("EnergyVSBackStrip_Singles%d", det))->Fill(backMaps[det][hit->GetSegment()], energy);
        }

        for(int j = 0; j < multGriffin; j++)
        {
            TGriffinHit* hitgrif = griffin.GetGriffinHit(j);

            double energygrif = hitgrif->GetEnergy();
            double timegrif = hitgrif->GetTimeStampNs()%60000000000LL;

            double tdiffgrif = (hit->GetTimeStampNs() - hitgrif->GetTimeStampNs());

            if(side == "P")
            {
                fH2[slot].at(Form("EnergyVSTimeDiffFrontGriffin_CondGRIF%d", det))->Fill(tdiffgrif, energy);        

                if(tdiffgrif >= 400 && tdiffgrif <= 1200) 
                {
                    fH2[slot].at(Form("EnergyFrontVSEnergyGriffin_CondRCMP%d", det))->Fill(energygrif, energy);

                    if(det != 3 && det != 4) //Deal with these detectors later
                    {
                        fH2[slot].at("EnergyFrontVSTimeFront_CondGRIF")->Fill(time, energy);
                        fH2[slot].at("EnergyGriffinVSTimeGriffin_CondRCMP")->Fill(timegrif, energygrif);
                    }
                }
            }

            if(side == "N")
            {
                fH2[slot].at(Form("EnergyVSTimeDiffBackGriffin_CondGRIF%d", det))->Fill(tdiffgrif, energy);            
            }
        }

        for(int j = i+1; j < mult; j++)
        {
            TRcmpHit* hitbis = rcmp.GetRcmpHit(j);

            int detbis = hitbis->GetDetector();
            
            double energybis = hitbis->GetEnergy();

            string sidebis = hitbis->GetChannel()->GetMnemonic()->CollectedChargeString();

            if((mult <= 10) && (side == sidebis) && (side == "P") && (5 <= det <= 6) && (5 <= detbis <= 6)) 
            {
                double tdiff = (hitbis->GetTimeStampNs() - hit->GetTimeStampNs());

                fH2[slot].at("EnergyFront1VSTimeDiffFrontFront_Det5Or6")->Fill(tdiff, energy);            

                if((-60. <= tdiff) && (tdiff <= 60.)) 
                {
                    fH2[slot].at(Form("EnergyFront1VSEnergyFront2Mult%d_Det5Or6", mult))->Fill(energybis, energy);   
                }

                if((720. <= tdiff) && (tdiff <= 780.))
                {
                    fH2[slot].at(Form("EnergyFront1VSEnergyFront2BgdMult%d_Det5Or6", mult))->Fill(energybis, energy);
                }
            }

            //cout << det << " " << detbis << " " << side << " " << sidebis << " " << energy << " " << energybis << endl;
        }

        for(int j = 0; j < mult; j++)
        {
            TRcmpHit* hitbis = rcmp.GetRcmpHit(j);

            int detbis = hitbis->GetDetector();

            double energybis = hitbis->GetEnergy();

            string sidebis = hitbis->GetChannel()->GetMnemonic()->CollectedChargeString();

            double tdiff = (hitbis->GetTimeStampNs() - hit->GetTimeStampNs());

            if((side != sidebis) && (det == detbis) && TMath::Abs(tdiff) < 300. && TMath::Abs(energybis - energy) < energy/1.8)
            {
                fH1[slot].at(Form("TimeDiffBackFront_Gated%d", det))->Fill(tdiff);
                
                if(side == "P") fH2[slot].at(Form("EnergyBackVSEnergyFront_Gated%d", det))->Fill(energy, energybis);
            }
        }
    }

    for(int i = 0; i < multGriffin; i++)
    {
        TGriffinHit* hit = griffin.GetGriffinHit(i);
        double timegrif = hit->GetTimeStampNs()%60000000000LL;
        double timegrifabs = hit->GetTimeStampNs();
        double energygrif = hit->GetEnergy();

        fH1[slot].at("EnergyGriffin_Singles")->Fill(energygrif);
        fH1[slot].at("TimeGriffin_Singles")->Fill(timegrif);
        fH2[slot].at("EnergyGriffinVSTimeGriffin_Singles")->Fill(timegrif, energygrif);
        fH1[slot].at("AbsTimeGriffin_Singles")->Fill(timegrifabs);
    }
}

void RCMPHelper::EndOfSort(std::shared_ptr<std::map<std::string, TList>>& list)
{

}
