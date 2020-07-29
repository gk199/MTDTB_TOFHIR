
#include<iostream>
#include<fstream>
#include<vector>
#include "TGraph.h"
#include "TMarker.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TROOT.h"
#include "TH1.h"
#include "TF1.h"
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

  // fitting function for the MIP peak vs OV plot
  TFile * inputFile = new TFile ("outFileFitHDR2.root", "READ");
  TF1 * myFit_pin = (TF1*) inputFile->Get("fitHDR2_signal");
  TF1 * myFit_flex = (TF1*) inputFile->Get("fitHDR2_signal");
  for (int i=2; i<5; i++) myFit_pin->FixParameter(i, myFit_pin->GetParameter(i));
  for (int i=2; i<5; i++) myFit_flex->FixParameter(i, myFit_flex->GetParameter(i));

  TCanvas * cMIP_peak = new TCanvas ("cMIP_peak", "cMIP_peak", 600, 500);
  cMIP_peak->cd();
  MIPpin_vs_OV->Draw("APE");
  MIPpin_vs_OV->GetYaxis()->SetRangeUser(0, 20);
  MIPpin_vs_OV->GetXaxis()->SetLimits(0,9);
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
  cMIP_peak->SaveAs(Form("MIP_Peak_Values_vs_OV.pdf"));
  cMIP_peak->SaveAs(Form("/eos/user/g/gkopp/www/BTL_TB/MIP_Peak_Values_vs_OV.pdf"));

  // fitting function, linear
  TF1 * linear_pin = new TF1("linear_pin","[0]*(x-[1])", 0,9);//"pol1");
  TF1 * linear_flex = new TF1("linear_flex","[0]*(x-[1])", 0,9);//,"pol1");
  MIPpin_vs_OV->Fit(linear_pin);
  MIPflex_vs_OV->Fit(linear_flex);
  linear_pin->Draw("same");
  linear_flex->Draw("same");
  std::cout << linear_pin->GetParameter(0) << " = offset for pin, and " << linear_flex->GetParameter(0) << " = offset for flex " << std::endl;
  std::cout << linear_pin->GetParameter(1) << " = slope for pin, and "<< linear_flex->GetParameter(1) << " = slope for flex " << std::endl;
  cMIP_peak->SaveAs(Form("LinearFit_MIP_Peak_Values_vs_OV.pdf"));
  cMIP_peak->SaveAs(Form("/eos/user/g/gkopp/www/BTL_TB/LinearFit_MIP_Peak_Values_vs_OV.pdf"));
  MIPpin_vs_OV->Fit(myFit_pin,"","",0,5);
  MIPflex_vs_OV->Fit(myFit_flex,"","",0,5);
  myFit_pin->Draw("same");
  myFit_flex->Draw("same");
  std::cout << myFit_pin->GetParameter(0) << " = offset for pin, and " << myFit_flex->GetParameter(0) << " = offset for flex, fancy fit" << std::endl;
  std::cout << myFit_pin->GetParameter(1) << " = slope for pin, and "<< myFit_flex->GetParameter(1) << " = slope for flex, fancy fit" << std::endl;
  cMIP_peak->SaveAs(Form("FancyFit_MIP_Peak_Values_vs_OV.pdf"));
  cMIP_peak->SaveAs(Form("/eos/user/g/gkopp/www/BTL_TB/FancyFit_MIP_Peak_Values_vs_OV.pdf"));

  TCanvas * cIC_RMS = new TCanvas ("cIC_RMS", "cIC_RMS", 600, 500);
  cIC_RMS->cd();
  RMSpin_vs_OV->Draw("LAPE");
  RMSpin_vs_OV->GetYaxis()->SetRangeUser(0, 0.23);
  RMSpin_vs_OV->GetXaxis()->SetLimits(0,9);
  RMSpin_vs_OV->SetTitle("RMS of Intercalibration Coefficient vs. OV");
  RMSpin_vs_OV->GetYaxis()->SetTitle("RMS of IC Values [a.u.]");
  RMSpin_vs_OV->GetXaxis()->SetTitle("Overvoltage (V)");
  RMSpin_vs_OV->SetMarkerStyle(20);
  RMSpin_vs_OV->SetMarkerColor(kBlue+1);
  RMSpin_vs_OV->SetLineColor(kBlue+1);
  RMSflex_vs_OV->SetMarkerStyle(21);
  RMSflex_vs_OV->SetMarkerColor(kRed+1);
  RMSflex_vs_OV->SetLineColor(kRed+1);
  RMSflex_vs_OV->Draw("same LPE");
  leg = new TLegend(0.15,0.75,0.5,0.88,NULL,"brNDC");
  leg->AddEntry(RMSpin_vs_OV, "Pin connected array", "lp" );
  leg->AddEntry(RMSflex_vs_OV, "Flex connected array", "lp" );
  leg->Draw();
  gPad->SetGridy();
  cIC_RMS->SaveAs(Form("RMS_of_IC_Values_vs_OV.pdf"));
  cIC_RMS->SaveAs(Form("/eos/user/g/gkopp/www/BTL_TB/RMS_of_IC_Values_vs_OV.pdf"));

  TCanvas * cIC_RMS_sub8 = new TCanvas ("cIC_RMS_sub8", "cIC_RMS_sub8", 600, 500);
  cIC_RMS_sub8->cd();
  RMSpin_sub8_vs_OV->Draw("LAPE");
  RMSpin_sub8_vs_OV->GetYaxis()->SetRangeUser(0, 0.23);
  RMSpin_sub8_vs_OV->GetXaxis()->SetLimits(0,9);
  RMSpin_sub8_vs_OV->SetTitle("#sqrt{RMS^{2}_{OV} - RMS^{2}_{8}} of Intercalibration Coefficient vs. OV");
  RMSpin_sub8_vs_OV->GetYaxis()->SetTitle("#sqrt{RMS_{OV} - RMS_{8}} of IC Values [a.u.]");
  RMSpin_sub8_vs_OV->GetXaxis()->SetTitle("Overvoltage (V)");
  RMSpin_sub8_vs_OV->SetMarkerStyle(20);
  RMSpin_sub8_vs_OV->SetMarkerColor(kBlue+1);
  RMSpin_sub8_vs_OV->SetLineColor(kBlue+1);
  RMSflex_sub8_vs_OV->SetMarkerStyle(21);
  RMSflex_sub8_vs_OV->SetMarkerColor(kRed+1);
  RMSflex_sub8_vs_OV->SetLineColor(kRed+1);
  RMSflex_sub8_vs_OV->Draw("same LPE");
  leg = new TLegend(0.15,0.75,0.5,0.88,NULL,"brNDC");
  leg->AddEntry(RMSpin_sub8_vs_OV, "Pin connected array", "lp" );
  leg->AddEntry(RMSflex_sub8_vs_OV, "Flex connected array", "lp" );
  leg->Draw();
  gPad->SetGridy();
  cIC_RMS_sub8->SaveAs(Form("RMS_sub8_of_IC_Values_vs_OV.pdf"));
  cIC_RMS_sub8->SaveAs(Form("/eos/user/g/gkopp/www/BTL_TB/RMS_sub8_of_IC_Values_vs_OV.pdf"));

  TCanvas * cIC_RMS_overlay = new TCanvas ("cIC_RMS_overlay", "cIC_RMS_overlay", 600, 500);
  cIC_RMS_overlay->cd();
  RMSpin_vs_OV->Draw("LAPE");
  RMSpin_vs_OV->GetYaxis()->SetRangeUser(0, 0.23);
  RMSpin_vs_OV->GetXaxis()->SetLimits(0,9);
  RMSpin_vs_OV->SetTitle("RMS and #sqrt{RMS^{2}_{OV} - RMS^{2}_{8}} of Intercalibration Coefficient vs. OV");
  RMSpin_vs_OV->GetYaxis()->SetTitle("RMS of IC Values [a.u.]");
  RMSpin_vs_OV->GetXaxis()->SetTitle("Overvoltage (V)");
  RMSpin_vs_OV->SetMarkerStyle(20);
  RMSpin_vs_OV->SetMarkerColor(kBlue+1);
  RMSpin_vs_OV->SetLineColor(kBlue+1);
  RMSflex_vs_OV->SetMarkerStyle(21);
  RMSflex_vs_OV->SetMarkerColor(kRed+1);
  RMSflex_vs_OV->SetLineColor(kRed+1);
  RMSflex_vs_OV->Draw("same LPE");
  RMSpin_sub8_vs_OV->SetMarkerStyle(20);
  RMSpin_sub8_vs_OV->SetMarkerColor(kBlue+1);
  RMSpin_sub8_vs_OV->SetLineColor(kBlue+1);
  RMSpin_sub8_vs_OV->SetLineStyle(2);
  RMSpin_sub8_vs_OV->Draw("same LPE");
  RMSflex_sub8_vs_OV->SetMarkerStyle(21);
  RMSflex_sub8_vs_OV->SetMarkerColor(kRed+1);
  RMSflex_sub8_vs_OV->SetLineColor(kRed+1);
  RMSflex_sub8_vs_OV->SetLineStyle(2);
  RMSflex_sub8_vs_OV->Draw("same LPE");
  leg = new TLegend(0.15,0.71,0.5,0.89,NULL,"brNDC");
  leg->AddEntry(RMSpin_vs_OV, "Pin connected array, RMS_{OV}", "lp" );
  leg->AddEntry(RMSflex_vs_OV, "Flex connected array, RMS_{OV}", "lp" );
  leg->AddEntry(RMSpin_sub8_vs_OV, "Pin connected array,  #sqrt{RMS^{2}_{OV} - RMS^{2}_{8}}", "lp" );
  leg->AddEntry(RMSflex_sub8_vs_OV, "Flex connected array,  #sqrt{RMS^{2}_{OV} - RMS^{2}_{8}}", "lp" );
  leg->Draw();
  gPad->SetGridy();
  cIC_RMS_overlay->SaveAs(Form("RMS_of_IC_Values_vs_OV_overlay.pdf"));
  cIC_RMS_overlay->SaveAs(Form("/eos/user/g/gkopp/www/BTL_TB/RMS_of_IC_Values_vs_OV_overlay.pdf"));

  return 0;
}
