
#include<iostream>
#include<fstream>
#include<vector>
#include "TGraph.h"
#include "TMarker.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TROOT.h"
#include "TH1.h"
#include "TH2.h"
#include "TProfile.h"
#include "TFile.h"


using namespace std;

int main() {
  double OV[5] = {8,6,4,3,2};
  double MIP_pin[5] = {18.21,12.27,7.588,5.319,2.746};
  double MIP_flex[5] = {14.30,10.21,6.485,4.363,1.478};
  double RMS_pin[5] = {0.041,0.044,0.056,0.068,0.10};
  double RMS_flex[5] = {0.053,0.059,0.077,0.12,0.17};
  double RMS_sub8_pin[5] = {0,0.01597,0.03814,0.05425,0.091209};
  double RMS_sub8_flex[5] = {0,0.025923,0.05586,0.10766,0.16153}; 

  TGraph *MIPpin_vs_OV = new TGraph (5, OV, MIP_pin);
  TGraph *MIPflex_vs_OV = new TGraph (5, OV, MIP_flex);
  TGraph *RMSpin_vs_OV = new TGraph (5, OV, RMS_pin);
  TGraph *RMSflex_vs_OV = new TGraph (5,OV, RMS_flex);
  TGraph *RMSpin_sub8_vs_OV = new TGraph (5, OV, RMS_sub8_pin);
  TGraph *RMSflex_sub8_vs_OV = new TGraph (5,OV, RMS_sub8_flex);
  TLegend *leg;

  TCanvas * cMIP_peak = new TCanvas ("cMIP_peak", "cMIP_peak", 600, 500);
  cMIP_peak->cd();
  MIPpin_vs_OV->Draw("APE");
  MIPpin_vs_OV->GetYaxis()->SetRangeUser(0, 20);
  MIPpin_vs_OV->SetTitle("MIP Peak vs. OV");
  MIPpin_vs_OV->GetYaxis()->SetTitle("MIP peak [a.u.]");
  MIPpin_vs_OV->GetXaxis()->SetTitle("Overvoltage (V)");
  MIPpin_vs_OV->SetMarkerStyle(20);
  MIPpin_vs_OV->SetMarkerColor(kBlue+1);
  MIPpin_vs_OV->SetLineColor(kBlue+1);
  MIPflex_vs_OV->SetMarkerStyle(21);
  MIPflex_vs_OV->SetMarkerColor(kRed+1);
  MIPflex_vs_OV->SetLineColor(kRed+1);
  MIPflex_vs_OV->Draw("same PE");
  leg = new TLegend(0.15,0.75,0.5,0.88,NULL,"brNDC");
  leg->AddEntry(MIPpin_vs_OV, "Pin connected array", "lp" );
  leg->AddEntry(MIPflex_vs_OV, "Flex connected array", "lp" );
  leg->Draw();
  gPad->SetGridy();
  cMIP_peak->SaveAs(Form("MIP Peak Values vs OV.pdf"));

  TCanvas * cIC_RMS = new TCanvas ("cIC_RMS", "cIC_RMS", 600, 500);
  cIC_RMS->cd();
  RMSpin_vs_OV->Draw("APE");
  RMSpin_vs_OV->GetYaxis()->SetRangeUser(0, 0.23);
  RMSpin_vs_OV->SetTitle("RMS of Intercalibration Coefficient vs. OV");
  RMSpin_vs_OV->GetYaxis()->SetTitle("RMS of IC Values [a.u.]");
  RMSpin_vs_OV->GetXaxis()->SetTitle("Overvoltage (V)");
  RMSpin_vs_OV->SetMarkerStyle(20);
  RMSpin_vs_OV->SetMarkerColor(kBlue+1);
  RMSpin_vs_OV->SetLineColor(kBlue+1);
  RMSflex_vs_OV->SetMarkerStyle(21);
  RMSflex_vs_OV->SetMarkerColor(kRed+1);
  RMSflex_vs_OV->SetLineColor(kRed+1);
  RMSflex_vs_OV->Draw("same PE");
  leg = new TLegend(0.15,0.75,0.5,0.88,NULL,"brNDC");
  leg->AddEntry(RMSpin_vs_OV, "Pin connected array", "lp" );
  leg->AddEntry(RMSflex_vs_OV, "Flex connected array", "lp" );
  leg->Draw();
  gPad->SetGridy();
  cIC_RMS->SaveAs(Form("RMS of IC Values vs OV.pdf"));

  TCanvas * cIC_RMS_sub8 = new TCanvas ("cIC_RMS_sub8", "cIC_RMS_sub8", 600, 500);
  cIC_RMS_sub8->cd();
  RMSpin_sub8_vs_OV->Draw("APE");
  RMSpin_sub8_vs_OV->GetYaxis()->SetRangeUser(0, 0.23);
  RMSpin_sub8_vs_OV->SetTitle("#sqrt{RMS_{OV} - RMS_{8}} of Intercalibration Coefficient vs. OV");
  RMSpin_sub8_vs_OV->GetYaxis()->SetTitle("#sqrt{RMS_{OV} - RMS_{8}} of IC Values [a.u.]");
  RMSpin_sub8_vs_OV->GetXaxis()->SetTitle("Overvoltage (V)");
  RMSpin_sub8_vs_OV->SetMarkerStyle(20);
  RMSpin_sub8_vs_OV->SetMarkerColor(kBlue+1);
  RMSpin_sub8_vs_OV->SetLineColor(kBlue+1);
  RMSflex_sub8_vs_OV->SetMarkerStyle(21);
  RMSflex_sub8_vs_OV->SetMarkerColor(kRed+1);
  RMSflex_sub8_vs_OV->SetLineColor(kRed+1);
  RMSflex_sub8_vs_OV->Draw("same PE");
  leg = new TLegend(0.15,0.75,0.5,0.88,NULL,"brNDC");
  leg->AddEntry(RMSpin_sub8_vs_OV, "Pin connected array", "lp" );
  leg->AddEntry(RMSflex_sub8_vs_OV, "Flex connected array", "lp" );
  leg->Draw();
  gPad->SetGridy();
  cIC_RMS_sub8->SaveAs(Form("RMS sub8 of IC Values vs OV.pdf"));

  return 0;
}
