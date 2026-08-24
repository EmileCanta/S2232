// energy_vs_run.C
// Build a 2D histogram: (ProjectionY of EnergyTwo) vs Run number
//
// Usage:
//   root -l -b -q 'energy_vs_run.C("/path/to/data")'

#include <vector>
#include <algorithm>
#include <regex>
#include <string>

#include "TSystemDirectory.h"
#include "TSystemFile.h"
#include "TList.h"
#include "TFile.h"
#include "TH1.h"
#include "TH2.h"
#include "TCanvas.h"
#include "TString.h"

struct RunFile {
    int run;
    std::string path;
};

void EnergyVSRun(const char* dir = ".") {

    // ---- 1. Collect files matching RCMP<run>_... .root and extract run number ----
    TSystemDirectory sysdir("sysdir", dir);
    TList* files = sysdir.GetListOfFiles();

    std::vector<RunFile> runFiles;
    std::regex runRegex(R"(RCMP(\d+)_.*\.root)");

    if (files) {
        TSystemFile* file;
        TIter next(files);
        while ((file = (TSystemFile*)next())) {
            TString fname = file->GetName();
            if (file->IsDirectory()) continue;
            std::string sname = fname.Data();
            std::smatch m;
            if (std::regex_match(sname, m, runRegex)) {
                int run = std::stoi(m[1].str());
                runFiles.push_back({run, std::string(dir) + "/" + sname});
            }
        }
    }

    if (runFiles.empty()) {
        printf("No matching RCMP*.root files found in %s\n", dir);
        return;
    }

    // Sort by run number, drop duplicates
    std::sort(runFiles.begin(), runFiles.end(),
              [](const RunFile& a, const RunFile& b){ return a.run < b.run; });
    runFiles.erase(std::unique(runFiles.begin(), runFiles.end(),
                   [](const RunFile& a, const RunFile& b){ return a.run == b.run; }),
                   runFiles.end());

    int nRuns = runFiles.size();
    printf("Found %d runs, from %d to %d\n", nRuns, runFiles.front().run, runFiles.back().run);

    // ---- 2. Peek at the first file's "EnergyVSFrontStrip_Singles1" TH2D to get the Y-axis binning ----
    TFile* f0 = TFile::Open(runFiles[0].path.c_str());
    if (!f0 || f0->IsZombie()) {
        printf("Could not open %s\n", runFiles[0].path.c_str());
        return;
    }
    TH2* hE2D0 = (TH2*)f0->Get("EnergyVSFrontStrip_Singles1");
    if (!hE2D0) {
        printf("No 'EnergyVSFrontStrip_Singles1' histogram in %s\n", runFiles[0].path.c_str());
        return;
    }
    int nYBins = hE2D0->GetNbinsY();
    double yLow  = hE2D0->GetYaxis()->GetXmin();
    double yHigh = hE2D0->GetYaxis()->GetXmax();
    f0->Close();

    // ---- 3. Build the output 2D histogram: sequential run bins, labeled with run number ----
    TH2F* h2 = new TH2F("h2_EnergyVsRun", "Energy spectra (Y-projection of EnergyVSFrontStrip_Singles1) vs Run;Run number;Energy",
                         nRuns, 0, nRuns,
                         nYBins, yLow, yHigh);

    // ---- 4. Fill: for each run, project EnergyVSFrontStrip_Singles1 onto Y (full X range), copy into that column ----
    for (int i = 0; i < nRuns; ++i) {
        TFile* f = TFile::Open(runFiles[i].path.c_str());
        if (!f || f->IsZombie()) {
            printf("Warning: could not open %s, skipping\n", runFiles[i].path.c_str());
            continue;
        }
        TH2* hE2D = (TH2*)f->Get("EnergyVSFrontStrip_Singles1");
        if (!hE2D) {
            printf("Warning: no 'EnergyVSFrontStrip_Singles1' histogram in %s, skipping\n", runFiles[i].path.c_str());
            f->Close();
            continue;
        }

        if (hE2D->GetNbinsY() != nYBins) {
            printf("Warning: run %d has different Y binning (%d bins vs %d), skipping\n",
                   runFiles[i].run, hE2D->GetNbinsY(), nYBins);
            f->Close();
            continue;
        }

        // Full projection onto Y: all X bins included (bin 0 to nbinsX+1 to be safe with over/underflow excluded)
        TH1D* proj = hE2D->ProjectionY(Form("py_%d", runFiles[i].run), 1, hE2D->GetNbinsX());

        for (int yb = 1; yb <= nYBins; ++yb) {
            double content = proj->GetBinContent(yb);
            if (content != 0)
                h2->SetBinContent(i + 1, yb, content);
        }

        h2->GetXaxis()->SetBinLabel(i + 1, Form("%d", runFiles[i].run));

        delete proj;
        f->Close();
    }

    h2->GetXaxis()->LabelsOption("v"); // vertical labels, since there are many runs
    h2->SetStats(0);

    TCanvas* c = new TCanvas("c", "Energy vs Run", 1400, 700);
    c->SetLogz();
    h2->Draw("COLZ");

    TFile* fout = new TFile("../histograms/rcmp/20mg/EnergyVSRunDet1.root", "RECREATE");
    h2->Write();
    fout->Close();
}
