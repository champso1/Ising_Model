#include "TROOT.h"
#include "TSystem.h"
#include "TStyle.h"

#include "TFile.h"
#include "TNtuple.h"
#include "TCanvas.h"
#include "TH2F.h"
#include "TAxis.h"

#include <iostream>


const TString single_data_file_path = "./res/plot_data/single_sim.root";
const TString multi_data_file_path = "./res/plot_data/multi_sim.root";
const TString single_plot_output_path = "./res/plot_outputs/single_sim/";
const TString multi_plot_output_path = "./res/plot_outputs/multi_sim/";

void Plot_SingleData(void) {
    /* See if the file exists */
    if (gSystem->AccessPathName(single_data_file_path, kFileExists)) {
        std::cout << "Single sim data does not exist!" << std::endl;
        return;
    }
    TFile *file = new TFile(single_data_file_path);
    TNtuple *data = (TNtuple *)file->Get("data");

    TCanvas c("canvas", "Data from single simulation", 100, 100, 800, 600);

    data->SetMarkerStyle(kDot);
    data->SetMarkerColor(kBlack);

    data->Draw("energy:N>>energy_hist");
    TH2F *energy_hist = (TH2F *)gPad->GetPrimitive("energy_hist");
    energy_hist->SetTitle("Energy vs. Iterations");
    energy_hist->GetXaxis()->SetTitle("N (iterations per lattice site)");
    energy_hist->GetYaxis()->SetTitle("E (energy in J/k_B)");

    c.Update();
    c.Print(single_plot_output_path + "energy.pdf");

    data->Draw("mag:N>>mag_hist");
    TH2F *mag_hist = (TH2F *)gPad->GetPrimitive("mag_hist");
    mag_hist->SetTitle("Magnetization vs. Iterations");
    mag_hist->GetXaxis()->SetTitle("N (iterations per lattice site)");
    mag_hist->GetYaxis()->SetTitle("M (magnetization)");

    c.Update();
    c.Print(single_plot_output_path + "mag.pdf");

    delete file;
}

void Plot_MultiData(void) {
    /* See if the file exists */
    if (gSystem->AccessPathName(multi_data_file_path, kFileExists)) {
        std::cout << "Multi sim data does not exist!" << std::endl;
        return;
    }
    TFile *file = new TFile(multi_data_file_path);
    TNtuple *data = (TNtuple *)file->Get("data");

    TCanvas c("canvas", "Data from multi simulation", 100, 100, 800, 600);

    data->SetMarkerStyle(kFullDotLarge);
    data->SetMarkerColor(kBlack);

    data->Draw("energy:T>>energy_hist");
    TH2F *energy_hist = (TH2F *)gPad->GetPrimitive("energy_hist");
    energy_hist->SetTitle("Energy vs. Temperature");
    energy_hist->GetXaxis()->SetTitle("T (temperature)");
    energy_hist->GetYaxis()->SetTitle("E (energy in J/k_B)");

    c.Update();
    c.Print(multi_plot_output_path + "energy.pdf");

    data->Draw("mag:T>>mag_hist");
    TH2F *mag_hist = (TH2F *)gPad->GetPrimitive("mag_hist");
    mag_hist->SetTitle("Magnetization vs. Temperature");
    mag_hist->GetXaxis()->SetTitle("T (temperature)");
    mag_hist->GetYaxis()->SetTitle("M (magnetization)");

    c.Update();
    c.Print(multi_plot_output_path + "mag.pdf");



    c.Clear();
    c.SetWindowSize(1600, 600);

    TPad *pad1 = new TPad("pad1", "Susceptibilities", 0.02, 0.02, 0.48, 0.98);
    TPad *pad2 = new TPad("pad2", "Specific Heats", 0.52, 0.02, 0.98, 0.98);
    pad1->Draw();
    pad2->Draw();

    pad1->cd();
    data->Draw("suscept:T>>suscept_hist");
    TH2F *suscept_hist = (TH2F *)gPad->GetPrimitive("suscept_hist");
    suscept_hist->SetTitle("Susceptibility vs. Temperature;T (temperature);chi (susceptibility)");

    pad2->cd();
    data->Draw("heat:T>>heat_hist");
    TH2F *heat_hist = (TH2F *)gPad->GetPrimitive("heat_hist");
    heat_hist->SetTitle("Specific Heat vs. Temperature;T (temperature);c (specific heat)");

    c.Update();
    c.Print(multi_plot_output_path + "suscept_spec_heats.pdf");

    delete file;
}



void plot(void) {
    Plot_SingleData();
    Plot_MultiData();
}