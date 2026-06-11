#include "TMath.h"

using namespace TMath;

void AngleHist()
{	
    TFile *filein = TFile::Open("/home/emile/postdoc/data/rcmp_sim/g4out.root", "READ");
    
    TFile *fileout = new TFile("./anglehist.root","RECREATE");
    
    TTree* tree = (TTree*)filein->Get("ntuple");
   
    ofstream out1("Det6.txt");
    ofstream out2("Det1.txt");
    ofstream out3("Det5.txt");
    ofstream out4("Det2.txt");
    ofstream out5("Det3.txt");
    ofstream out6("Det4.txt");

    int detnbr;
    double posx;
    double posy;
    double posz;
    
    double theta_ini;
    double phi_ini;

    double rad = 180./Pi();
   
    vector<vector<double>> thetaVec(6144);
    vector<vector<double>> phiVec(6144);

    tree->SetBranchAddress("detNumber", &detnbr);
    tree->SetBranchAddress("posx",&posx);
    tree->SetBranchAddress("posy",&posy);
    tree->SetBranchAddress("posz",&posz);
    
    tree->SetBranchAddress("primaryTheta",&theta_ini);
    tree->SetBranchAddress("primaryPhi",&phi_ini);

    TH2D *histsTheta[6]; 
    TH2D *histsPhi[6]; 
    TH2D *histsRho[6]; 
    TH2D *histsHits[6];
    
    TH2D *histsThetaIni[6]; 
    TH2D *histsPhiIni[6]; 

    vector<int> map = {5,0,4,1,2,3};

    TH3D *hist3D = new TH3D("hist", "hist", 20, -100., 100., 20, -100., 100., 20, -100., 100.);

    int fEntries = tree->GetEntries();
    
    for(int i =0; i<6; i++) 
    {
        TString nameHistTheta = Form("histTheta%d",i+1);
        TString nameHistPhi = Form("histPhi%d",i+1);
        TString nameHistRho = Form("histRho%d",i+1);
        TString nameHistHits = Form("histHits%d",i+1);
        
        TString nameHistPhiIni = Form("histPhiIni%d",i+1);
        TString nameHistThetaIni = Form("histThetaIni%d",i+1);
        
        histsTheta[i] = new TH2D(nameHistTheta, nameHistTheta, 6146, 0, 6146, 4000, -200, 200);
        histsPhi[i] = new TH2D(nameHistPhi, nameHistPhi, 6146, 0, 6146, 4000, -200, 200);
        histsRho[i] = new TH2D(nameHistRho, nameHistRho, 6146, 0, 6146, 1000, 0, 100);
        histsHits[i] = new TH2D(nameHistHits, nameHistHits, 34, -1, 33, 34, -1, 33);
        
        histsThetaIni[i] = new TH2D(nameHistThetaIni, nameHistThetaIni, 6146, 0, 6146, 4000, -200, 200);
        histsPhiIni[i] = new TH2D(nameHistPhiIni, nameHistPhiIni, 6146, 0, 6146, 4000, -200, 200);
    }
    
    for(int j = 0; j < fEntries; j++)
    {
        tree->GetEntry(j);

        hist3D->Fill(posx, posy, posz);

        double theta = ACos(posz/Sqrt(pow(posx,2)+pow(posy,2)+pow(posz,2))); 
        double phi = ATan2(posy,posx); 
        double rho = Sqrt(pow(posx,2)+pow(posy,2)+pow(posz,2)); 
        double newTheta = ACos(sin(theta)*sin(phi));
        double newPhi = ATan2(sin(theta)*cos(phi),cos(theta));
        double newTheta_ini = ACos(sin(theta_ini)*sin(phi_ini));
        double newPhi_ini = ATan2(sin(theta_ini)*cos(phi_ini),cos(theta_ini));
        
        for(int i =1; i<=6; i++)
        {
            if((1024*i-1023)<=detnbr && detnbr<=(i*1024))
            {
                histsTheta[map[i-1]]->Fill(detnbr, newTheta*rad);
                histsPhi[map[i-1]]->Fill(detnbr, newPhi*rad);
                histsRho[map[i-1]]->Fill(detnbr, rho);

                histsThetaIni[map[i-1]]->Fill(detnbr, newTheta_ini*rad);
                histsPhiIni[map[i-1]]->Fill(detnbr, newPhi_ini*rad);

                int modDet = detnbr - (i-1)*1024;
                int strip1 = (modDet - 1) % 32;
                int strip2 = (modDet - 1) / 32;

                if(map[i-1] <= 2) histsHits[map[i-1]]->Fill(31 - strip1, strip2); //Swapping order on x-axis to fit strip number definitions.
                if(map[i-1] >= 3) histsHits[map[i-1]]->Fill(strip1, strip2);

                thetaVec[detnbr-1].push_back(newTheta_ini*rad);
                phiVec[detnbr-1].push_back(newPhi_ini*rad);

            }
        }
    }

    for(int k = 1; k <= 6; k++)
    {
        string filename = "Det" + to_string(map[k-1]+1) + ".txt";
        ofstream outFile(filename);

        outFile << "phi " << "std_phi " << "theta " << "std_theta" << '\n';

        for(int i = (k*1024-1024); i<(k*1024); i++)
        { 
            double phiSum=0.;    
            double thetaSum=0.;    
            double phiSum2=0.;    
            double thetaSum2=0.;    

            double Nevents = phiVec[i].size();

            for(int j = 0; j<Nevents; j++)
            {
                phiSum += phiVec[i][j];
                thetaSum += thetaVec[i][j];

                phiSum2 += phiVec[i][j]*phiVec[i][j];
                thetaSum2 += thetaVec[i][j]*thetaVec[i][j];
            }

            if(i % 32 == 0)
            {
                outFile << "---" << '\n';
            }

            double meanPhi = phiSum/Nevents;
            double meanPhi2 = phiSum2/Nevents;
            double meanTheta = thetaSum/Nevents;
            double meanTheta2 = thetaSum2/Nevents;

            double stdPhi = sqrt(meanPhi2-(meanPhi*meanPhi));
            double stdTheta = sqrt(meanTheta2-(meanTheta*meanTheta));

            outFile << fixed << setprecision(3) << meanPhi << " " << stdPhi << " " << meanTheta << " " << stdTheta << '\n';
        }

        outFile.close();
    }

    fileout->Write();
}
