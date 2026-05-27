#include "/home/emile/GRSISort/include/TChannel.h"
#include "/home/emile/GRSISort/include/TRWPeak.h"

void UpdateRCMPEnergyCalibration() 
{
    TFile* analysis_file = new TFile("/home/emile/postdoc/data/s2232/sorted/analysis/analysis28264_000.root", "UPDATE");

    analysis_file->Get("Channel");

    vector<Float_t> OldCoefsP;
    vector<Float_t> OldCoefsN;

    double a = 1.00118;
    double b = -5.7388;

    for(int i = 1; i <= 6; i++)
    {
        for(int j = 0; j <= 31; j++)
        {
            TChannel* detector_channel_p = TChannel::FindChannelByName(Form("RCS%02dXP%02dX", i, j));
            TChannel* detector_channel_n = TChannel::FindChannelByName(Form("RCS%02dXN%02dX", i, j));

            if (!detector_channel_p || !detector_channel_n) 
            {
                std::cerr << "Channel not found" << std::endl;
            }

            vector<Float_t> OldCoeffsP = detector_channel_p->GetENGCoeff();

            detector_channel_p->DestroyENGCal();
            detector_channel_p->AddENGCoefficient(a*OldCoeffsP[0]+b);
            detector_channel_p->AddENGCoefficient(a*OldCoeffsP[1]);

            cout << "bip" << endl;   

            vector<Float_t> OldCoeffsN = detector_channel_n->GetENGCoeff();

            detector_channel_n->DestroyENGCal();
            detector_channel_n->AddENGCoefficient(a*OldCoeffsN[0]+b);
            detector_channel_n->AddENGCoefficient(a*OldCoeffsN[1]);

            cout << "boup" << endl;  
        }
    }

    TChannel::WriteCalFile("/home/emile/postdoc/analysis/s2232/sorter/calfiles/energy_calibration28264.cal");

    TChannel::WriteToRoot();

    analysis_file->Close();
}
