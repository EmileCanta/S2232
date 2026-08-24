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

    fH2[slot]["Mult2DetRepartition"] = new TH2F("Mult2DetRepartition", "Mult2DetRepartition", 6, 1, 7, 6, 1, 7);
    
    fH2[slot]["TDiffVSFrontStrip"] = new TH2F("TDiffVSFrontStrip", "TDiffVSFrontStrip", 36, -2, 34, 2000, -1000, 1000);
    fH2[slot]["TDiffVSFrontStripStamp"] = new TH2F("TDiffVSFrontStripStamp", "TDiffVSFrontStripStamp", 36, -2, 34, 2000, -1000, 1000);

    for(int ndet = 1; ndet <= 6; ndet++)
    {
        fH2[slot][Form("EnergyVSFrontStrip%d", ndet)] = new TH2F(Form("EnergyVSFrontStrip%d", ndet), Form("EnergyVSFrontStrip%d", ndet), 36, -2, 34, 10000, 0, 10000);
        fH2[slot][Form("EnergyVSFrontStripBis%d", ndet)] = new TH2F(Form("EnergyVSFrontStripBis%d", ndet), Form("EnergyVSFrontStripBis%d", ndet), 36, -2, 34, 10000, 0, 10000);
        fH1[slot][Form("EnergyAddback%d", ndet)] = new TH1F(Form("EnergyAddback%d", ndet), Form("EnergyAddback%d", ndet), 10000, 0, 10000);
        fH2[slot][Form("EnergyVSBackStrip%d", ndet)] = new TH2F(Form("EnergyVSBackStrip%d", ndet), Form("EnergyVSBackStrip%d", ndet), 36, -2, 34, 10000, 0, 10000);
        fH2[slot][Form("EnergyVSBackStripBis%d", ndet)] = new TH2F(Form("EnergyVSBackStripBis%d", ndet), Form("EnergyVSBackStripBis%d", ndet), 36, -2, 34, 10000, 0, 10000);
        fH2[slot][Form("TDiffVSFrontStrip%d", ndet)] = new TH2F(Form("TDiffVSFrontStrip%d", ndet), Form("TDiffVSFrontStrip%d", ndet), 36, -2, 34, 2000, -1000, 1000);
        fH2[slot][Form("TDiffVSFrontStripCfd%d", ndet)] = new TH2F(Form("TDiffVSFrontStripCfd%d", ndet), Form("TDiffVSFrontStripCfd%d", ndet), 36, -2, 34, 2000, -1000, 1000);
        
        //Careful binning 1 bin = 4 keV
        fH2[slot][Form("ChargeVSFrontStrip%d", ndet)] = new TH2F(Form("ChargeVSFrontStrip%d", ndet), Form("ChargeVSFrontStrip%d", ndet), 36, -2, 34, 2500, 0, 10000);
        fH2[slot][Form("ChargeVSBackStrip%d", ndet)] = new TH2F(Form("ChargeVSBackStrip%d", ndet), Form("ChargeVSBackStrip%d", ndet), 36, -2, 34, 2500, 0, 10000);

        //Careful binning 1 bin = 10 keV
        fH2[slot][Form("ChargeVSPixelFront%d", ndet)] = new TH2F(Form("ChargeVSPixelFront%d", ndet), Form("ChargeVSPixelFront%d", ndet), 1028, -2, 1026, 1000, 0, 10000);
        fH2[slot][Form("ChargeVSPixelBack%d", ndet)] = new TH2F(Form("ChargeVSPixelBack%d", ndet), Form("ChargeVSPixelBack%d", ndet), 1028, -2, 1026, 1000, 0, 10000);

        fH2[slot][Form("HitMapNonCorrected%d", ndet)] = new TH2F(Form("HitMapNonCorrected%d", ndet), Form("HitMapNonCorrected%d", ndet), 34, -1, 33, 34, -1, 33);
        fH2[slot][Form("HitMapCorrectedGood%d", ndet)] = new TH2F(Form("HitMapCorrectedGood%d", ndet), Form("HitMapCorrectedGood%d", ndet), 34, -1, 33, 34, -1, 33);
        fH2[slot][Form("HitMapCorrectedWeird%d", ndet)] = new TH2F(Form("HitMapCorrectedWeird%d", ndet), Form("HitMapCorrectedWeird%d", ndet), 34, -1, 33, 34, -1, 33);

        //Careful, hist limits smaller than usual here.
        fH2[slot][Form("EnergyGriffinVSEnergyFrontRCMP%d", ndet)] = new TH2F(Form("EnergyGriffinVSEnergyFrontRCMP%d", ndet), Form("EnergyGriffinVSEnergyFrontRCMP%d", ndet), 5000, 0, 5000, 3000, 0, 3000);
        fH2[slot][Form("EnergyGriffinVSEnergyFrontRCMPBis%d", ndet)] = new TH2F(Form("EnergyGriffinVSEnergyFrontRCMPBis%d", ndet), Form("EnergyGriffinVSEnergyFrontRCMPBis%d", ndet), 5000, 0, 5000, 3000, 0, 3000);
    }
}

void RCMPHelper::Exec(unsigned int slot, TRcmp& rcmp, TGriffin& griffin, TGriffinBgo& griffinbgo)
{
    frontstrip.clear();
    backstrip.clear();
    backenergy.clear();
    frontenergy.clear();
    frontdet.clear();
    backdet.clear();

    backaddbackvector.clear();
    frontaddbackvector.clear();

    firststripfront = 0;
    firstenergyfront = 0.;
    addbackenergyfront = 0.;

    firststripback = 0;
    firstenergyback = 0.;
    addbackenergyback = 0.;

    int countfront = 0;
    bool isaddbackfront = false;
    int countback = 0;
    bool isaddbackback = false;

    TRcmpHit* hit1;
    TRcmpHit* hit2;

    double mult = rcmp.GetMultiplicity();
    double multGriffin = griffin.GetMultiplicity();

    //MULT2 FILLER//

    /*if(mult == 2) 
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

        //cout << det1 << " " << strip1 << " " << mappedstrip1 << endl;

        string side1 = hit1->GetChannel()->GetMnemonic()->CollectedChargeString();
        string side2 = hit2->GetChannel()->GetMnemonic()->CollectedChargeString();

        int pixel = mappedstrip1 + 32 * mappedstrip2;

        if((det1 == det2) && (side1 != side2)) fH1[slot].at("Energy")->Fill(hit1->GetEnergy());

        if((det1 == det2) && (side1 != side2))
        {
            fH2[slot].at(Form("EnergyVSFrontStrip%d", det1))->Fill(mappedstrip1, hit1->GetEnergy());
            fH2[slot].at(Form("EnergyVSBackStrip%d", det1))->Fill(mappedstrip2, hit2->GetEnergy());
            fH2[slot].at(Form("ChargeVSFrontStrip%d", det1))->Fill(strip1, hit1->GetCharge()); //Careful to not use mapping for calibration code to work correctly!
            fH2[slot].at(Form("ChargeVSBackStrip%d", det1))->Fill(strip2, hit2->GetCharge());

            //fH2[slot].at(Form("HitMapNonCorrected%d", det1))->Fill(strip1, strip2);
            //fH2[slot].at(Form("HitMapCorrectedGood%d", det1))->Fill(mappedstrip1, mappedstrip2);
            //fH2[slot].at(Form("HitMapCorrectedWeird%d", det1))->Fill(mappedstrip1bis, mappedstrip2bis);

            fH2[slot].at(Form("ChargeVSPixelFront%d", det1))->Fill(pixel, hit1->GetCharge());
            fH2[slot].at(Form("ChargeVSPixelBack%d", det1))->Fill(pixel, hit2->GetCharge());

            if(det1 == 1) fH2[slot].at(Form("HitMapCorrectedGood%d", det1))->Fill(32 - mappedstrip2, mappedstrip1); 
            if(det1 == 2) fH2[slot].at(Form("HitMapCorrectedGood%d", det1))->Fill(mappedstrip1, 32 - mappedstrip2);

            if(det1 == 3) fH2[slot].at(Form("HitMapCorrectedGood%d", det1))->Fill(mappedstrip1, mappedstrip2);
            if(det1 == 4) fH2[slot].at(Form("HitMapCorrectedGood%d", det1))->Fill(mappedstrip2, mappedstrip1);

            if(det1 == 5) fH2[slot].at(Form("HitMapCorrectedGood%d", det1))->Fill(mappedstrip1, mappedstrip2);
            if(det1 == 6) fH2[slot].at(Form("HitMapCorrectedGood%d", det1))->Fill(mappedstrip2, mappedstrip1);

            for(int i = 0; i < multGriffin; i++)
            {
                TGriffinHit* hit = griffin.GetGriffinHit(i);
                fH2[slot].at(Form("EnergyGriffinVSEnergyFrontRCMP%d", det1))->Fill(hit->GetEnergy(), hit1->GetEnergy()); //Careful, no time condition here!
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
            
                fH2[slot].at("TDiffVSFrontStrip")->Fill(mappedstrip1, (hit1->GetTime() - hit2->GetTime()));
                fH2[slot].at("TDiffVSFrontStripStamp")->Fill(mappedstrip1, (hit1->GetTimeStamp() - hit2->GetTimeStamp()));
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
        
        fH2[slot].at("Mult2DetRepartition")->Fill(det1, det2);
    }*/
    
    //Filling temp front/back vectors

    for(int i = 0; i < mult; i++)
    {
        TRcmpHit* hit = rcmp.GetRcmpHit(i);

        int det = hit->GetDetector();
        string side = hit->GetChannel()->GetMnemonic()->CollectedChargeString();
        double energy = hit->GetEnergy();

        if(side == "P") 
        {
            int mappedstrip = frontMaps[det][hit->GetSegment()];
            
            frontenergy.push_back(energy);
            frontstrip.push_back(mappedstrip);
            frontdet.push_back(det);
        }

        if(side == "N") 
        {
            int mappedstrip = backMaps[det][hit->GetSegment()];
            
            backenergy.push_back(energy);
            backstrip.push_back(mappedstrip);
            backdet.push_back(det);
        }
    }

    //PRINTER//

    /*if(mult > 2)
    {
        cout << "mult : " << mult << endl;

        for(int i = 0; i < frontenergy.size(); i++)
        {
            cout << "front energy : " << frontenergy[i] << endl;
            cout << "front det : " << frontdet[i] << endl;
            cout << "front strip : " << frontstrip[i] << endl;
        }

        for(int i = 0; i < backenergy.size(); i++)
        {
            cout << "back energy : " << backenergy[i] << endl;
            cout << "back det : " << backdet[i] << endl;
            cout << "back strip : " << backstrip[i] << endl;
        }
    }*/

    //Physical hitter//

    //A physical hit is a hit that deposited a similar energy in front and back strips. It can be addback or non addback.

    /*for(int i = 0; i < frontenergy.size(); i++)
    {
        double en1 = frontenergy[i];
        double st1 = frontstrip[i];

        int counter = 0;

        for(int j = i+1; j < frontenergy.size(); j++)
        {
            double en2 = frontenergy[j];
            double st2 = frontstrip[j];

            if(TMath::Abs(st2-st1) == 1)
            {
                cout << "front addback energy : " << en1+en2 << endl;
                frontaddbackvector.push_back(en1+en2);
                counter++;
                break; 
            }

         }
        //this doesnt work, I count as non addback some addback events
        if(counter == 0) frontnonaddbackvector.push_back(en1);
    }

    for(int i = 0; i < backenergy.size(); i++)
    {
        double en1 = backenergy[i];
        double st1 = backstrip[i];

        for(int j = i+1; j < backenergy.size(); j++)
        {
            double en2 = backenergy[j];
            double st2 = backstrip[j];

            if(TMath::Abs(st2-st1) == 1)
            {
                cout << "back addback energy : " << en1+en2 << endl;
                backaddbackvector.push_back(en1+en2);
            }
        }
    }*/

    //ADDBACKER//

    for(int i = 0; i < mult; i++)
    {
        TRcmpHit* hit = rcmp.GetRcmpHit(i);

        int det = hit->GetDetector();

        string side = hit->GetChannel()->GetMnemonic()->CollectedChargeString();

        if(side == "P") //add condition on det has to exist in the elements of backdet 
        {
            if(countfront == 0) 
            {
                firstenergyfront = hit->GetEnergy();
                firststripfront = frontMaps[det][hit->GetSegment()];

                countfront++;
            }

            if(countfront > 0 && (frontMaps[det][hit->GetSegment()] == firststripfront + 1 || frontMaps[det][hit->GetSegment()] == firststripfront - 1))
            {
                addbackenergyfront = firstenergyfront+hit->GetEnergy();
                isaddbackfront = true;  
            }
        }

        if(side == "N" && mult == 3) //add condition on det has to exist in the elements of backdet 
        {
            if(countback == 0) 
            {
                firstenergyback = hit->GetEnergy();
                firststripback = backMaps[det][hit->GetSegment()];

                countback++;
            }

            if(countback > 0 && (backMaps[det][hit->GetSegment()] == firststripback + 1 || backMaps[det][hit->GetSegment()] == firststripback - 1))
            {
                addbackenergyback = firstenergyback+hit->GetEnergy();
                isaddbackback = true;  
            }
        }
    }

    //FILLER//

    for(int i = 0; i < mult; i++)
    {
        TRcmpHit* hit = rcmp.GetRcmpHit(i);

        int det = hit->GetDetector();

        string side = hit->GetChannel()->GetMnemonic()->CollectedChargeString();

        fH2[slot].at("DetectorVSMultiplicity")->Fill(det, mult, 1./mult); // Careful, weighted histogram here

        //cout << i << " " << backstrip.size() << " " << side << endl;

        if(side == "P" && backstrip.size() > 0) //add condition on det has to exist in the elements of backdet 
        {
            //fH2[slot].at(Form("EnergyVSFrontStripBis%d", det))->Fill(frontMaps[det][hit->GetSegment()], hit->GetEnergy());
            
            if(isaddbackfront)
            {
                if(isaddbackback && TMath::Abs(addbackenergyfront - addbackenergyback) < 60.)
                {
                    fH1[slot].at(Form("EnergyAddback%d", det))->Fill(addbackenergy);
                    break;
                }

                if(!isaddbackback && TMath::Abs(addbackenergyfront - addbackenergyback) < 60.)
                {
                    fH1[slot].at(Form("EnergyAddback%d", det))->Fill(addbackenergy);
                    break;
                }
            }
            
            if(!isaddback)
            {
                fH2[slot].at(Form("EnergyVSFrontStripBis%d", det))->Fill(frontMaps[det][hit->GetSegment()], hit->GetEnergy());
            }
        }

        if(side == "N") 
        {
            fH2[slot].at("SideVSMultiplicity")->Fill(mult, 1);
            fH2[slot].at(Form("EnergyVSBackStripBis%d", det))->Fill(backMaps[det][hit->GetSegment()], hit->GetEnergy());
        }

        if(mult > 0) fH1[slot].at("Multiplicity")->Fill(mult, 1./mult); // Careful, weighted histogram here

        //cout << hit->GetEnergy() << " " << side << " " << frontMaps[det][hit->GetSegment()] << endl;
    }

    //cout << "NEXT EVENT" << endl;

    //Beta spectrum attempt//

    /*for(int i = 0; i < mult; i++)
    {
        TRcmpHit* hit = rcmp.GetRcmpHit(i);

        int det = hit->GetDetector();

        string side = hit->GetChannel()->GetMnemonic()->CollectedChargeString();
        
        for(int j = 0; j < multGriffin; j++)
        {
            TGriffinHit* hitgrif = griffin.GetGriffinHit(j);

            if(side == "P") fH2[slot].at(Form("EnergyGriffinVSEnergyFrontRCMPBis%d", det))->Fill(hitgrif->GetEnergy(), hit->GetEnergy()); //Careful, no time condition here!
        }

        break;
    }

    for(int i = 0; i < multGriffin; i++)
    {
        TGriffinHit* hit = griffin.GetGriffinHit(i);
        fH1[slot].at("EnergyGriffin")->Fill(hit->GetEnergy());
    }*/
}

void RCMPHelper::EndOfSort(std::shared_ptr<std::map<std::string, TList>>& list)
{

}
