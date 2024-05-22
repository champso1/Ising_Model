#include "TROOT.h"
#include "TSystem.h"
#include "TFile.h"
#include "TNtuple.h"
#include "TCanvas.h"
#include "TFrame.h"
#include "TGraph.h"
#include "TH2F.h"
#include "TStyle.h"

#include <iostream>

const TString single_file_name = "./res/single_sim.root";

void plot(void) {
    TFile *tfile = nullptr;
    if (gSystem->AccessPathName(single_file_name, kFileExists)) {
        std::cout << "Single sim data does not exist!" << std::endl;
        return;
    }
    tfile = TFile::Open(single_file_name);

    TCanvas *c1 = new TCanvas("c1", "Single sim plots", 100,100, 1600,600);
    TPad *pad1 = new TPad("pad1", "Energies", 0.02, 0.02, 0.48, 0.98, 0);
    TPad *pad2 = new TPad("pad2", "Magnetizations", 0.52, 0.02, 0.98, 0.98, 0);
    pad1->Draw();
    pad2->Draw();
    TNtuple *data = (TNtuple *)tfile->Get("data");
    data->SetMarkerStyle(6); // larger dots

    pad1->cd();
    data->Draw("energy:N");
    pad2->cd();
    data->Draw("mag:N");

    // drawing any TTree will create a histogram called "htemp",
    // from which we can make some adjustments to the stuff
    // TH2F *htemp = (TH2F *)gPad->GetPrimitive("htemp");
    // htemp->SetTitle("Energy vs. Iterations (per lattice size)");
    // htemp->GetXaxis()->SetTitle("Number of iterations (per lattice site) N");
    // htemp->GetYaxis()->SetTitle("Energy (J/k_B)");
    // gStyle->SetTitleAlign(33);
    // gStyle->SetTitleX(0.99);

    c1->Update();
    c1->Print("./res/single_sim.pdf");

    delete tfile;
}