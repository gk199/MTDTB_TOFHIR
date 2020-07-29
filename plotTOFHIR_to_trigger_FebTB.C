// to compile: g++ -Wall -o plotTOFHIR_to_trigger plotTOFHIR_to_trigger_Feb.C `root-config --cflags --glibs` -lSpectrum

#include <cstdlib>
#include <iostream>
#include <fstream>
#include <iomanip>
#include <map>
#include <vector>
#include <string>

#include "TFile.h"
#include "TTree.h"
#include "TBranch.h"
#include "TChain.h"
#include "TGraph.h"
#include "TGraphErrors.h"
#include "TProfile2D.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TProfile.h"
#include "TF1.h"
#include "TCanvas.h"
#include "TLatex.h"
#include "TLegend.h"
#include "TLine.h"
#include "TApplication.h"
#include "TFormula.h"

#include "TFile.h"
#include "TStyle.h"
#include "TSystem.h"
#include "TKey.h"

#include "TString.h"
#include "TTree.h"
#include "TBranch.h"

#include "TSpline.h"
#include "TCanvas.h"
#include "TSpectrum.h"
// #include "TDirectory.h." 
#include "TObject.h"
#include <algorithm>
#include "functions.hh"

int main(int argc, char** argv)
{
  // open the file that important things are written to
  std::ofstream myfile;
  myfile.open ("example.txt");
  myfile << "Saving bar position and MPV from the lambda fit." << std::endl;

    gROOT->SetBatch(true);
    TApplication* theApp = new TApplication("App", &argc, argv);
    //TLegend *leg;

    //read input files - based on first and last run number listed in the command line arguments
    int firstRun = 1;    
    
    if (argc > 1) 
    {
        firstRun= atoi(argv[1]);
    }
    std::cout << "plotting from run = " << firstRun << std::endl;
    
    int lastRun = firstRun;
    if (argc > 2) 
    {
        lastRun= atoi(argv[2]);
    }
    
    std::string data_path = "/eos/cms/store/group/dpg_mtd/comm_mtd/TB/MTDTB_FNAL_Feb2020/TOFHIR/RecoData/v1/RecoWithTracks/";  
    if (argc > 3) 
    {
      std::cout << "using data path from input line" << std::endl;
      data_path = std::string(argv[3]);
    }
    std::cout << "data_path = " << data_path << std::endl;
    
    
    //define tree
    float	step1;
    float	step2;
    float	x_dut;
    float	y_dut;
    int 	ntracks;        
    long long   time[400];    
    float	tot[256];
    float       qfine[256];
    float       energy[256];
            
    TChain * tree = new TChain("data", "data");
    
    // going over each of the runs listed as command line arguments
    for (int iRun = firstRun; iRun <= lastRun; iRun++)
      {
	//if ( std::find(bad_runs.begin(), bad_runs.end(), iRun) != bad_runs.end()) continue;       
	tree->Add( Form("%s/run%.5d_events.root", data_path.c_str(), iRun) );  
	std::cout << "adding run: " << iRun << std::endl;
      }

    // extract information from the ROOT tree
    tree->SetBranchStatus("*",0);
    tree->SetBranchAddress("step1", &step1); // overvoltage
    tree->SetBranchAddress("step2", &step2); // (threshold 1 * 10000 + 1) + (threshold 2 * 100 + 1) + (threshold E + 1), first 2 digits are vth1, second 2 are vth2, last 2 are energy 
    tree->SetBranchAddress("x_dut", &x_dut);
    tree->SetBranchAddress("y_dut", &y_dut);
    tree->SetBranchAddress("ntracks", &ntracks);    
    tree->SetBranchAddress("tot", &tot);
    tree->SetBranchAddress("time", &time);
    tree->SetBranchAddress("qfine", &qfine); // qfine is integral of amplified signal in range defined by tot
    tree->SetBranchAddress("energy", &energy); // energy found from calibration plot of qfine vs tot, and then difference between calibration tot and measured
    
    //count how many steps are contained in the file
    std::cout << "loop to count steps in file" << std::endl;
    int NSTEP1 = 0;
    int NSTEP2 = 0;        
    std::vector<float> step1_vct;
    std::vector<float> step2_vct;
    step1_vct.clear();
    step2_vct.clear();
        
    // starting the EVENT loop now
    Long64_t NEVENTS = tree->GetEntries();
    std::cout << "nEvents = " << NEVENTS << std::endl;    
    for (Int_t iEvt= 0; iEvt < NEVENTS; iEvt++) 
      {
        tree->GetEntry(iEvt);
        //std::cout << "step1 = " << step1 << " :: step2 = " << step2 <<std::endl;
	// getting vth1 vth2 and vthe from step2 for each event
	float vth1 = float(int(step2/10000)-1);;
	float vth2 = float(int((step2-10000*vth1)/100)-1);
	float vthe = float(int((step2-10000*vth1-step2-100*vth2)/1)-1);

	if (ntracks != 1) continue;
        if (std::find(step1_vct.begin(), step1_vct.end(), step1) == step1_vct.end() )
	  {             
	    if (step1 > 1.e-2)
	      {
		step1_vct.push_back(step1);
		NSTEP1++;
	      }
	  }                
        if (std::find(step2_vct.begin(), step2_vct.end(), step2) == step2_vct.end() )
	  {            
            step2_vct.push_back(step2);
            NSTEP2++;
	  }
      }
     
    // step 1 is overvoltage                                                                                                                                           
    for (int iStep1 = 0; iStep1<NSTEP1; iStep1++)
      {
        std::cout << "step1 (" << iStep1+1 << "/" << NSTEP1 << "): " << step1_vct.at(iStep1)  << std::endl;        
      }
    // step 2 is vth1, vth2, vthe in each two digits
    for (int iStep2 = 0; iStep2<NSTEP2; iStep2++)
      {
        std::cout << "step2 (" << iStep2+1 << "/" << NSTEP2 << "): " << step2_vct.at(iStep2)  << std::endl;        
      }    

    // define histos    
    double minTot = 0;
    double maxTot = 700;
    double minEnergy = 0;
    double maxEnergy = 50;
    double minTime = -200000;
    double maxTime = 200000;
    double minXpos = -10;
    double maxXpos = 40;
    double minYpos = -10;
    double maxYpos = 50;
    int REBIN_COEFF = 32;
    
    #define NCH 400

    // this is for configuration 2, 3 (connector was rotated wrong), 4, or 5
    // conf 6 and 8 have the horizontal array as the low channel numbers. Config 9 has both arrays vertical - mostly focusing on this one
    //    int myChList[] = {0,1,2,3,4,5,6,7,8,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30,31,32}; // from channelMapping1 in mtd drawMatrices.cfg, VERTICAL
    //    int myChList1[] = {3,10,0,1,7,14,5,12,21,20,23,22,16,18,17,19}; // one side from channelMapping1
    //    int myChList2[] = {4,6,15,8,13,2,11,27,32,31,30,29,28,26,24,25}; // one side from channelMapping1
    //    int myChList[] = {57,63,50,60,59,55,61,56,58,53,62,54,9,51,38,52,34,43,33,42,36,44,35,46,37,45,39,47,41,48,40,49}; // channelMapping2 HORIZONTAL
    //    int myChList1[] = {57,50,59,61,58,62,9,38,34,33,36,35,37,39,41,40}; // one side from channelMapping2, totalEnergy[0]
    //    int myChList2[] = {63,60,55,56,53,54,51,52,43,42,44,46,45,47,48,49}; // one side from channelMapping2, totalEnergy[1]
    // mapping for the pin connector array, caltech array
    int myChList[] = {63,57,60,50,55,59,56,61,53,58,54,62,51,9,52,38,40,49,41,48,39,47,37,45,35,46,36,44,33,42,34,43}; // VERTICAL
    int myChList1[] = {63,57,60,50,55,59,56,61,53,58,54,62,51,9,52,38}; // one side from channelMapping2, totalEnergy[0]
    int myChList2[] = {40,49,41,48,39,47,37,45,35,46,36,44,33,42,34,43}; // one side from channelMapping2, totalEnergy[1]
    int NBARS = 16; // full array has 16 bars, 32 SiPM readouts
    double center[NCH] = {0};
    double MIP[NCH] = {0};
    double IC[NCH] = {0};
    double MIP_corr[NCH] = {0}; // initialize these arrays to 0 and then only define values for the channels we are working with
    double IC_corr[NCH] = {0};

    // center of bar x position, conf 9.1
    center[63] = center[40] = 0;
    center[57] = center[49] = 1.2;
    center[60] = center[41] = 1.2;
    center[50] = center[48] = 1.1;
    center[55] = center[39] = 2.4;
    center[59] = center[47] = 6;
    center[56] = center[37] = 8.4;
    center[61] = center[45] = 11.9;
    center[53] = center[35] = 15; 
    center[58] = center[46] = 18.4;
    center[54] = center[36] = 21.6;
    center[62] = center[44] = 24.8;
    center[51] = center[33] = 28.2;
    center[9] = center[42] = 31.2; 
    center[52] = center[34] = 34.3;
    center[38] = center[43] = 37.9;

    // list MIP peak energies based off of Landau fit, this is for v1 reconstruction
    MIP[45] = MIP[56] = MIP[61] = 290;
    MIP[41] = MIP[33] = MIP[34] = MIP[38] = MIP[39] = MIP[42] = MIP[46] = MIP[47] = MIP[57] = MIP[58] = MIP[59] = MIP[60] = MIP[51] = 300;
    MIP[35] = MIP[44] = MIP[49] = 305;
    MIP[36] = MIP[37] = 295;
    MIP[48] = MIP[52] = MIP[53] = MIP[54] = MIP[55] = MIP[62] = 310;
    MIP[9] = MIP[50] = 315;
    MIP[40] = MIP[43] = MIP[63] = 0; // these channels don't work

// list MIP peak energies based off of Landau fit, this is for v5 recomstruction - using Landau to CORRECTED TOT
    MIP_corr[9] = 12.2;
    MIP_corr[33] = 10.4;
    MIP_corr[34] = 11.6;
    MIP_corr[35] = 11.7;
    MIP_corr[36] = 10.6;
    MIP_corr[37] = 11.6;
    MIP_corr[38] = 11; // no events when cut around this bar?? 
    MIP_corr[39] = 11.1;
    MIP_corr[41] = 10.9;
    MIP_corr[42] = 11.1;
    MIP_corr[44] = 11.6;
    MIP_corr[45] = 11.1;
    MIP_corr[46] = 11.6;
    MIP_corr[47] = 11.9;
    MIP_corr[48] = 11.2;
    MIP_corr[49] = 11.9;
    MIP_corr[50] = 12.0;
    MIP_corr[51] = 12.0;
    MIP_corr[52] = 12.5;
    MIP_corr[53] = 12.7;
    MIP_corr[54] = 12.6;
    MIP_corr[55] = 12.5;
    MIP_corr[56] = 11.2;
    MIP_corr[57] = 12.1;
    MIP_corr[58] = 12.2;
    MIP_corr[59] = 12.1;
    MIP_corr[60] = 10.3;
    MIP_corr[61] = 11.1;
    MIP_corr[62] = 12.8;

    MIP_corr[40] = MIP_corr[43] = MIP_corr[63] = 0; // these channels don't work

    // calculate intercallibration coefficients
    double avgMIP = 0;
    double avgMIP_corr = 0;
    for (int iCh = 0; iCh < NCH; iCh++) 
      {
	avgMIP += MIP[iCh];
	avgMIP_corr += MIP_corr[iCh];
      }
    avgMIP = avgMIP / 29; // we have 29 active channels!!
    avgMIP_corr = avgMIP_corr / 29;

    for (int iCh= 0; iCh < NCH; iCh++)
      { 
	IC[iCh] = MIP[iCh] / avgMIP;
	IC_corr[iCh] = MIP_corr[iCh] / avgMIP_corr;
      }
 
 // declare the histograms, these will be filled in the channel loop    
 std::map<float, std::map<float, std::map<int, TH1F * > > > hTot;
 std::map<float, std::map<float, std::map<int, TH1F * > > > hTot_correction;
 std::map<float, std::map<float, std::map<int, TH1F * > > > hTot_cut;
 std::map<float, std::map<float, std::map<int, TH1F * > > > hTot_cut_correction;
 std::map<float, std::map<float, std::map<int, TH1F * > > > hTime;
 std::map<float, std::map<float, std::map<int, TProfile * > > > pTot_vs_Xpos;
 std::map<float, std::map<float, std::map<int, TProfile * > > > pTot_vs_Ypos;
 std::map<float, std::map<float, std::map<int, TH1F * > > > pEff_vs_Xpos; 
 std::map<float, std::map<float, std::map<int, TH1F * > > > pEff_vs_Ypos;
 std::map<float, std::map<float, std::map<int, TH2F * > > > pXpos_Ypos_Tot;
 std::map<float, std::map<float, std::map<int, TH1F * > > > hCTR_UD;
 std::map<float, std::map<float, TProfile * > > pTot_vs_Xpos_overlay;
 std::map<float, std::map<float, TProfile * > > pTot_vs_Ypos_overlay;
 std::map<float, std::map<float, std::map<int, TProfile2D * > > > pXY_Edep;
 std::map<float, std::map<float, std::map<int, TH1F * > > > pCrossTalk;
 std::map<float, std::map<float, std::map<int, TH1F * > > > pCrossTalkBar;
 std::map<float, std::map<float, std::map<int, TH1F * > > > pCrossTalkBar1away;
 std::map<float, std::map<float, std::map<int, TH1F * > > > pCrossTalkBar2away;
 std::map<float, std::map<float, std::map<int, TH1F * > > > pCrossTalk_corr;
 std::map<float, std::map<float, std::map<int, TH1F * > > > pCrossTalkBar_corr;
 std::map<float, std::map<float, std::map<int, TH1F * > > > pCrossTalkBar1away_corr;
 std::map<float, std::map<float, std::map<int, TH1F * > > > pCrossTalkBar2away_corr;
 
 //    std::map<float, std::map<float, std::map<int, TH1F * > > > pCrossTalkOverlay;

 TH1F * hPosX = new TH1F ("hPosX", "hPosX", 200, minXpos, maxXpos);
 TH1F * hPosY = new TH1F ("hPosY", "hPosY", 200, minYpos, maxYpos);
 
 // step 1, step 2, and channel listing loops. Histograms are defined inside the loops
 for (int iStep1 = 0; iStep1< NSTEP1; iStep1++)
   {
     for (int iStep2 = 0; iStep2< NSTEP2; iStep2++)
       {
	 for (int iCh = 0; iCh < NCH; iCh++)
	   {       
	     hTot[step1_vct.at(iStep1)][step2_vct.at(iStep2)][iCh] = new TH1F (Form("hTot_ch%.3d_step1vct_%.1f_step2vct_%.1f_step1_%i_step2_%i", iCh, step1_vct.at(iStep1), step2_vct.at(iStep2), iStep1, iStep2), Form("Time over threshold, ch%.3d, step1_%.1f, step2_%.1f", iCh, step1_vct.at(iStep1), step2_vct.at(iStep2)), 4000, minTot, maxTot );
	     hTot_correction[step1_vct.at(iStep1)][step2_vct.at(iStep2)][iCh] = new TH1F (Form("hTot_corr_ch%.3d_step1vct_%.1f_step2vct_%.1f_step1_%i_step2_%i", iCh, step1_vct.at(iStep1), step2_vct.at(iStep2), iStep1, iStep2), Form("Energy, ch%.3d, step1_%.1f, step2_%.1f", iCh, step1_vct.at(iStep1), step2_vct.at(iStep2)), 4000, minEnergy, maxEnergy);
	     hTot_cut[step1_vct.at(iStep1)][step2_vct.at(iStep2)][iCh] = new TH1F (Form("hTot_cut_ch%.3d_step1vct_%.1f_step2vct_%.1f_step1_%i_step2_%i", iCh, step1_vct.at(iStep1), step2_vct.at(iStep2), iStep1, iStep2), Form("Time over threshold with x cut, ch%.3d, step1_%.1f, step2_%.1f", iCh, step1_vct.at(iStep1), step2_vct.at(iStep2)), 4000, minTot, maxTot );
	     hTot_cut_correction[step1_vct.at(iStep1)][step2_vct.at(iStep2)][iCh] = new TH1F (Form("hTot_cut_corr_ch%.3d_step1vct_%.1f_step2vct_%.1f_step1_%i_step2_%i", iCh, step1_vct.at(iStep1), step2_vct.at(iStep2), iStep1, iStep2), Form("Energy with x cut, ch%.3d, step1_%.1f, step2_%.1f", iCh, step1_vct.at(iStep1), step2_vct.at(iStep2)), 4000, minEnergy, maxEnergy );
	     hTime[step1_vct.at(iStep1)][step2_vct.at(iStep2)][iCh] = new TH1F (Form("hTime_ch%.3d_step1vct_%.1f_step2vct_%.1f_step1_%i_step2_%i", iCh, step1_vct.at(iStep1), step2_vct.at(iStep2), iStep1, iStep2), Form("Time Hist, ch%.3d, step1_%.1f, step2_%.1f", iCh, step1_vct.at(iStep1), step2_vct.at(iStep2)), 4000, minTime, maxTime );
	     pTot_vs_Xpos[step1_vct.at(iStep1)][step2_vct.at(iStep2)][iCh] = new TProfile (Form("pTot_vs_Xpos_ch%.3d_step1vct_%.1f_step2vct_%.1f_step1_%i_step2_%i", iCh, step1_vct.at(iStep1), step2_vct.at(iStep2), iStep1, iStep2), Form("ToT vs. X pos, ch%.3d, step1_%.1f, step2_%.1f", iCh, step1_vct.at(iStep1), step2_vct.at(iStep2)), 4000, minXpos, maxXpos );
	     pTot_vs_Ypos[step1_vct.at(iStep1)][step2_vct.at(iStep2)][iCh] = new TProfile (Form("pTot_vs_Ypos_ch%.3d_step1vct_%.1f_step2vct_%.1f_step1_%i_step2_%i", iCh, step1_vct.at(iStep1), step2_vct.at(iStep2), iStep1, iStep2), Form("ToT vs. Y pos, ch%.3d, step1_%.1f, step2_%.1f", iCh, step1_vct.at(iStep1), step2_vct.at(iStep2)), 4000, minYpos, maxYpos );
	     pEff_vs_Xpos[step1_vct.at(iStep1)][step2_vct.at(iStep2)][iCh] = new TH1F (Form("pEff_vs_Xpos_ch%.3d_step1vct_%.1f_step2vct_%.1f_step1_%i_step2_%i", iCh, step1_vct.at(iStep1), step2_vct.at(iStep2), iStep1, iStep2), Form("Efficiency vs. X pos, ch%.3d, step1_%.1f, step2_%.1f", iCh, step1_vct.at(iStep1), step2_vct.at(iStep2)), 400, 0, maxXpos);
	     pEff_vs_Ypos[step1_vct.at(iStep1)][step2_vct.at(iStep2)][iCh] = new TH1F (Form("pEff_vs_Ypos_ch%.3d_step1vct_%.1f_step2vct_%.1f_step1_%i_step2_%i", iCh, step1_vct.at(iStep1), step2_vct.at(iStep2), iStep1, iStep2), Form("Efficiency vs. Y pos, ch%.3d, step1_%.1f, step2_%.1f", iCh, step1_vct.at(iStep1), step2_vct.at(iStep2)), 400, 0, maxYpos);
	     pXpos_Ypos_Tot[step1_vct.at(iStep1)][step2_vct.at(iStep2)][iCh] = new TH2F (Form("pXpos_Ypos_Tot_ch%.3d_step1vct_%.1f_step2vct_%.1f_step1_%i_step2_%i", iCh, step1_vct.at(iStep1), step2_vct.at(iStep2), iStep1, iStep2), Form("X and Y pos vs. ToT, ch%.3d, step1_%.1f, step2_%.1f", iCh, step1_vct.at(iStep1), step2_vct.at(iStep2)), 40, minXpos, maxXpos, 40, minYpos, maxYpos );	     
	     pXY_Edep[step1_vct.at(iStep1)][step2_vct.at(iStep2)][iCh] = new TProfile2D (Form("pXY_Edep_ch%.3d_step1_%.1f_step2_%.1f", iCh, step1_vct.at(iStep1), step2_vct.at(iStep2)), Form("XY Energy dep ch%.3d_step1_%.1f_step2_%.1f", iCh, step1_vct.at(iStep1), step2_vct.at(iStep2)), 400, minXpos, maxXpos, 400, minXpos, maxXpos );

	     pCrossTalk[step1_vct.at(iStep1)][step2_vct.at(iStep2)][iCh] = new TH1F (Form("pCrossTalk_ch%.3d_step1_%.1f_step2_%.1f", iCh, step1_vct.at(iStep1), step2_vct.at(iStep2)), Form("Cross Talk ch%.3d_step1_%.1f_step2_%.1f", iCh, step1_vct.at(iStep1), step2_vct.at(iStep2)), 400, 0, 1);
	     pCrossTalkBar[step1_vct.at(iStep1)][step2_vct.at(iStep2)][iCh] = new TH1F (Form("pCrossTalkBar_ch%.3d_step1_%.1f_step2_%.1f", iCh, step1_vct.at(iStep1), step2_vct.at(iStep2)), Form("Cross Talk (with bar cut) ch%.3d_step1_%.1f_step2_%.1f", iCh, step1_vct.at(iStep1), step2_vct.at(iStep2)), 400, 0, 1);
	     pCrossTalkBar1away[step1_vct.at(iStep1)][step2_vct.at(iStep2)][iCh] = new TH1F (Form("pCrossTalkBar1away_ch%.3d_step1_%.1f_step2_%.1f", iCh, step1_vct.at(iStep1), step2_vct.at(iStep2)), Form("Cross Talk 1 away (with bar cut) ch%.3d_step1_%.1f_step2_%.1f", iCh, step1_vct.at(iStep1), step2_vct.at(iStep2)), 400, 0, 1);
	     pCrossTalkBar2away[step1_vct.at(iStep1)][step2_vct.at(iStep2)][iCh] = new TH1F (Form("pCrossTalkBar2away_ch%.3d_step1_%.1f_step2_%.1f", iCh, step1_vct.at(iStep1), step2_vct.at(iStep2)), Form("Cross Talk 2 away (with bar cut) ch%.3d_step1_%.1f_step2_%.1f", iCh, step1_vct.at(iStep1), step2_vct.at(iStep2)), 400, 0, 1);
	     
	     //cross talk with corrected tot
	     pCrossTalk_corr[step1_vct.at(iStep1)][step2_vct.at(iStep2)][iCh] = new TH1F (Form("pCrossTalk_corr_ch%.3d_step1_%.1f_step2_%.1f", iCh, step1_vct.at(iStep1), step2_vct.at(iStep2)), Form("Cross Talk (energy) ch%.3d_step1_%.1f_step2_%.1f", iCh, step1_vct.at(iStep1), step2_vct.at(iStep2)), 400, 0, 1);
	     pCrossTalkBar_corr[step1_vct.at(iStep1)][step2_vct.at(iStep2)][iCh] = new TH1F (Form("pCrossTalkBar_corr_ch%.3d_step1_%.1f_step2_%.1f", iCh, step1_vct.at(iStep1), step2_vct.at(iStep2)), Form("Cross Talk (with bar cut, energy) ch%.3d_step1_%.1f_step2_%.1f", iCh, step1_vct.at(iStep1), step2_vct.at(iStep2)), 400, 0, 1);
	     pCrossTalkBar1away_corr[step1_vct.at(iStep1)][step2_vct.at(iStep2)][iCh] = new TH1F (Form("pCrossTalkBar1away_corr_ch%.3d_step1_%.1f_step2_%.1f", iCh, step1_vct.at(iStep1), step2_vct.at(iStep2)), Form("Cross Talk 1 away (with bar cut, energy) ch%.3d_step1_%.1f_step2_%.1f", iCh, step1_vct.at(iStep1), step2_vct.at(iStep2)), 400, 0, 1);
	     pCrossTalkBar2away_corr[step1_vct.at(iStep1)][step2_vct.at(iStep2)][iCh] = new TH1F (Form("pCrossTalkBar2away_corr_ch%.3d_step1_%.1f_step2_%.1f", iCh, step1_vct.at(iStep1), step2_vct.at(iStep2)), Form("Cross Talk 2 away (with bar cut, energy) ch%.3d_step1_%.1f_step2_%.1f", iCh, step1_vct.at(iStep1), step2_vct.at(iStep2)), 400, 0, 1);
	   }

	 pTot_vs_Xpos_overlay[step1_vct.at(iStep1)][step2_vct.at(iStep2)] = new TProfile (Form("pTot_vs_Xpos_over_step1vct_%.1f_step2vct_%.1f_step1_%i_step2_%i", step1_vct.at(iStep1), step2_vct.at(iStep2), iStep1, iStep2), Form("ToT vs. X pos overlay, step1_%.1f, step2_%.1f", step1_vct.at(iStep1), step2_vct.at(iStep2)), 4000, minXpos, maxXpos );
	 pTot_vs_Ypos_overlay[step1_vct.at(iStep1)][step2_vct.at(iStep2)] = new TProfile (Form("pTot_vs_Ypos_over_step1vct_%.1f_step2vct_%.1f_step1_%i_step2_%i", step1_vct.at(iStep1), step2_vct.at(iStep2), iStep1, iStep2), Form("ToT vs. Y pos overlay, step1_%.1f, step2_%.1f", step1_vct.at(iStep1), step2_vct.at(iStep2)), 4000, minYpos, maxYpos );

	 // bar loop (outside of channel loop) to define time resolution for a single bar
	 for (int iBar = 0; iBar < NBARS; iBar++)
	   {
	     hCTR_UD[step1_vct.at(iStep1)][step2_vct.at(iStep2)][iBar] = new TH1F (Form("hCTR_UD_ch%.3d_step1vct_%.1f_step2vct_%.1f_step1_%i_step2_%i", iBar, step1_vct.at(iStep1), step2_vct.at(iStep2), iStep1, iStep2), Form("Bar Time Res, UD_ch%.3d, step1_%.1f, step2_%.1f", iBar, step1_vct.at(iStep1), step2_vct.at(iStep2)), 4000, -1000, 1000 );
	   }                       
       }
   }
    
 //************************************************************************************//
 //              loop 0
 //************************************************************************************//
 
 std::cout << "(0) looping over events to get MIP peak position" << std::endl;
 for (Int_t iEvt= 0; iEvt < NEVENTS; iEvt++) 
   {
     //std::cout << "processing event before 1 MIP req "<< iEvt<< std::endl;
     tree->GetEntry(iEvt);
     if (ntracks != 1) continue; // require 1 MIP per event
     if (iEvt%1000 == 0) std::cout << "processing event: " << iEvt << "\r" << std::flush;
     //	if (step1!=6) continue;
     //hCTR[step1][step2]->Fill(time2-time1);
     long long time_ref = time[384];
     hPosX->Fill(x_dut);
     hPosY->Fill(y_dut);
     // channel loop to calculate total energy deposit
     double TotalEnergy[2] = {0};
     double TotalEnergy_corr[2] = {0};
     double corrected_tot[NCH] = {0}; // initialize the corrected_tot for each channel to 0 at the start of each event loop
     for (int iCh = 0; iCh<NCH; iCh++) // channel loop to calculate total energy per event, per side of the array (top half or bottom half)
       {
	 // if not one of channels that we care about, skip it
	 if (std::find(std::begin(myChList), std::end(myChList), iCh) == std::end(myChList) ) continue;
	 // want to find the total energy in the event, in any channel
	 if ( (qfine[iCh] <= 14) || (qfine[iCh] > 500) || (energy[iCh] == 0) || (energy[iCh] == 20) || (energy[iCh] == 40) || (energy[iCh] >= 60)) continue; // skip any thing where the energy isnt a normal value (should always be postive, default value is -9999) for both energy and tot. Excluse high and low qfine regions, and energy values where there is a spike (multiplies of 20)
	 corrected_tot[iCh] = 14.13 * (exp(0.01562 * (tot[iCh]/1.e3))-1); // corrected energy given the correction to the tot linearization         

	 // define total energy separately for the top and bottom half of the bars. TotalEnergy[0] = even, TotalEnergy[1] = odd channels
	 if (std::find(std::begin(myChList1), std::end(myChList1), iCh) == std::end(myChList1) ) {
	   TotalEnergy[0] += (tot[iCh]/1.e3) / IC[iCh];
	   TotalEnergy_corr[0] += energy[iCh] / IC_corr[iCh];
	 }
	 if (std::find(std::begin(myChList2), std::end(myChList2), iCh) == std::end(myChList2) ) {
	   TotalEnergy[1] += (tot[iCh]/1.e3) / IC[iCh];
	   TotalEnergy_corr[1] += energy[iCh] / IC_corr[iCh];
	 }
       }

     for (int iCh = 0; iCh<NCH; iCh++) // channel loop in the event loop
       {
	 if (std::find(std::begin(myChList), std::end(myChList), iCh) == std::end(myChList) ) continue;
         if ( (qfine[iCh] <= 14) || (qfine[iCh] > 500) || (energy[iCh] == 0) || (energy[iCh] == 20) || (energy[iCh] == 40) || (energy[iCh] >= 60)) continue; // skip any thing where the energy isnt a normal value (should always be postive, default value is -9999) for both energy and tot. Excluse high and low qfine regions, and energy values where there is a spike (multiplies of 20) 
	 //if (iCh%2 == 1 ) continue;
	 //if (step1 != 6) continue;
	 //if (step2 != 0) continue;
	 
	 hTot[step1][step2][iCh]->Fill(tot[iCh]/1.e3); // raw ToT information, no cuts
	 hTot_correction[step1][step2][iCh]->Fill(energy[iCh]); // cut out low region where energy and tot are non-linear with each other
	 pXpos_Ypos_Tot[step1][step2][iCh]->Fill(x_dut, y_dut);
	 hTime[step1][step2][iCh]->Fill(time[iCh] - time_ref);
	 
	 // calculate the cross talk, and plot this
	 // plot Tot, normalized by MIP peak energy, and then as a fraction of the total energy deposited in all channels
	 if (std::find(std::begin(myChList1), std::end(myChList1), iCh) == std::end(myChList1) ) { // then use TotalEnergy[0]
	   if ((TotalEnergy[0] < 4000) && (tot[iCh]/1.e3 > 5) && (tot[iCh]/1.e3 < 500) && (energy[iCh] > -1)) // tot>5 for a zero supression from low energy deposits, no cut around MIP peak, this is blue on the cross talk plots. Energy > -1 to make sure event is good since used for cross talk plots
	     {
	       pCrossTalk[step1][step2][iCh]->Fill(((tot[iCh]/1.e3) / IC[iCh] )  / TotalEnergy[0] );
	     }
	   if ((TotalEnergy_corr[0] < 4000) && (qfine[iCh]>13) && (energy[iCh] > -1) )
	     {
	       pCrossTalk_corr[step1][step2][iCh]->Fill((energy[iCh] / IC_corr[iCh] )  / TotalEnergy_corr[0] );
	     }
	 }
	 if (std::find(std::begin(myChList2), std::end(myChList2), iCh) == std::end(myChList2) ) { // then use TotalEnergy[1]
	   if ((TotalEnergy[1] < 4000) && (tot[iCh]/1.e3 > 5) && (tot[iCh]/1.e3 < 500) && (energy[iCh] > -1)) // tot>5 for a zero supression from low energy deposits, no cut around MIP peak, this is blue on the cross talk plots
	     {
	       pCrossTalk[step1][step2][iCh]->Fill(((tot[iCh]/1.e3) / IC[iCh] )  / TotalEnergy[1] );
	     }
	   if ((TotalEnergy_corr[1] < 4000) && (qfine[iCh]>13)&& (energy[iCh] > -1) )
	     {
	       pCrossTalk_corr[step1][step2][iCh]->Fill((energy[iCh] / IC_corr[iCh] )  / TotalEnergy_corr[1] );
	     }
	 }
	 
	 // cut on a specific bar to see how this affects MIP peak (expect to pick out Landau peak for one bar) 
	 // do this for each channel, based off of the stats found from the fit to the efficiency plots
	 // for the cross talk (green plot) cut around the MIP position requring signal in the central bar to be 0.85*MIP - 4*MIP
	 // MIP and x_dut cut should be redundant actually
	 if (x_dut < center[iCh]+1 && x_dut > center[iCh]-1 )
	   {
	     hTot_cut[step1][step2][iCh]->Fill(tot[iCh]/1.e3);
	     hTot_cut_correction[step1][step2][iCh]->Fill(energy[iCh]);
	     
	     if (std::find(std::begin(myChList1), std::end(myChList1), iCh) == std::end(myChList1) ) { // then use TotalEnergy[0]
	       if ((TotalEnergy[0] < 4000) && (tot[iCh]/1.e3 > 0.85*MIP[iCh]) && (tot[iCh]/1.e3<4*MIP[iCh])) {
		 pCrossTalkBar[step1][step2][iCh]->Fill(((tot[iCh]/1.e3) / IC[iCh] )  / TotalEnergy[0] );
	       }
	       if ((energy[iCh] > 0.85*MIP_corr[iCh]) && (energy[iCh] < 4*MIP_corr[iCh]) && (qfine[iCh] > 13)) { // energy MIP peak cut, but exclusing low tot region where tot vs. energy is inversely related
		 pCrossTalkBar_corr[step1][step2][iCh]->Fill((energy[iCh] / IC_corr[iCh] )  / TotalEnergy_corr[0] );
	       }
	     }
	     if (std::find(std::begin(myChList2), std::end(myChList2), iCh) == std::end(myChList2) ) { // then use TotalEnergy[1]
	       if ((TotalEnergy[1] < 4000) && (tot[iCh]/1.e3 > 0.85*MIP[iCh]) && (tot[iCh]/1.e3<4*MIP[iCh])) {
		 pCrossTalkBar[step1][step2][iCh]->Fill(((tot[iCh]/1.e3) / IC[iCh] )  / TotalEnergy[1] );
	       }
	       if ((energy[iCh] > 0.85*MIP_corr[iCh]) && (energy[iCh] < 4*MIP_corr[iCh]) && (qfine[iCh] > 13)) { // energy MIP peak cut, but exclusing low tot region where tot vs. energy is inversely related 
		 pCrossTalkBar_corr[step1][step2][iCh]->Fill((energy[iCh] / IC_corr[iCh] )  / TotalEnergy_corr[1] ); 
	       }
	     }
	   } // closing loop based on position cut

	 // fill the 1 bar away cross talk for channels. Restrict to only consider ones that are not right on the edge
	 int xtalkchannel = 0;
	 int ChList1Pos = -100;
	 int ChList2Pos = -100;
	 // for cross talk analysis, iCh is where cross talk hit is originating from, and xtalkch = myChList1[ChList1Pos+-1] is the channel we care about xtalk in (this is what is plotted on final overlay distributions)
	 for (int ListPos = 0; ListPos < 16; ListPos++ ) // find which position in myChList we are in
	   {
	     if (iCh == myChList1[ListPos]) ChList1Pos = ListPos;
	     if (iCh == myChList2[ListPos]) ChList2Pos = ListPos;
	   }
	 if  (x_dut < center[iCh]+1 && x_dut > center[iCh]-1 ) 
	   {
	     // put a MIP cut on the bar +-1 away from iCh of interest
	     if (tot[iCh]/1.e3 > 0.85*MIP[iCh] && tot[iCh]/1.e3<4*MIP[iCh] ) // no totalEnergy cut since already require non-zero tot on this side of array
	       {
		 if ( ChList1Pos > 1 && ChList1Pos < 14 ) {
		   pCrossTalkBar1away[step1][step2][myChList1[ChList1Pos-1]]->Fill(((tot[myChList1[ChList1Pos-1]]/1.e3) / IC[myChList1[ChList1Pos-1]] ) / TotalEnergy[0] );
		   pCrossTalkBar1away[step1][step2][myChList1[ChList1Pos+1]]->Fill(((tot[myChList1[ChList1Pos+1]]/1.e3) / IC[myChList1[ChList1Pos+1]] ) / TotalEnergy[0] );
		 }
		 if ( ChList2Pos > 1 && ChList2Pos < 14 ) {
                   pCrossTalkBar1away[step1][step2][myChList2[ChList2Pos-1]]->Fill(((tot[myChList2[ChList2Pos-1]]/1.e3) / IC[myChList2[ChList2Pos-1]] ) / TotalEnergy[1] );
                   pCrossTalkBar1away[step1][step2][myChList2[ChList2Pos+1]]->Fill(((tot[myChList2[ChList2Pos+1]]/1.e3) / IC[myChList2[ChList2Pos+1]] ) / TotalEnergy[1] );
		 }
	       } // ToT MIP peak cut closing

	     // Energy MIP peak cut: separate out for channel +- 1 away to see if these have different distributions. Using third index as iCh = myChList1[ChList1Pos] instead of myChList1[ChList1Pos-1] inorder to separate out x talk from neighboring +-1 bars. When plotting need to associate to which channel this is x talking to
	     if ((energy[iCh] > 0.85*MIP_corr[iCh]) && (energy[iCh] < 4*MIP_corr[iCh]) && (qfine[iCh] > 13)) { // energy MIP peak cut, but exclusing low tot region where tot vs. energy is inversely related 
	       if (ChList1Pos > 1 && ChList1Pos < 14 ) {
		 pCrossTalkBar1away_corr[step1][step2][iCh]->Fill((energy[myChList1[ChList1Pos-1]] / IC_corr[myChList1[ChList1Pos-1]] )  / TotalEnergy_corr[0] );
		 pCrossTalkBar1away_corr[step1][step2][iCh]->Fill((energy[myChList1[ChList1Pos+1]] / IC_corr[myChList1[ChList1Pos+1]] )  / TotalEnergy_corr[0] );
		 pCrossTalkBar2away_corr[step1][step2][iCh]->Fill((energy[myChList1[ChList1Pos-2]] / IC_corr[myChList1[ChList1Pos-2]] )  / TotalEnergy_corr[0] );
		 pCrossTalkBar2away_corr[step1][step2][iCh]->Fill((energy[myChList1[ChList1Pos+2]] / IC_corr[myChList1[ChList1Pos+2]] )  / TotalEnergy_corr[0] );
	       }
	       if ( ChList2Pos > 1 && ChList2Pos < 14 ) {
		 pCrossTalkBar1away_corr[step1][step2][iCh]->Fill((energy[myChList2[ChList2Pos-1]] / IC_corr[myChList2[ChList2Pos-1]] )  / TotalEnergy_corr[1] );
		 pCrossTalkBar1away_corr[step1][step2][iCh]->Fill((energy[myChList2[ChList2Pos+1]] / IC_corr[myChList2[ChList2Pos+1]] )  / TotalEnergy_corr[1] );
		 pCrossTalkBar2away_corr[step1][step2][iCh]->Fill((energy[myChList2[ChList2Pos-2]] / IC_corr[myChList2[ChList2Pos-2]] )  / TotalEnergy_corr[1] );
		 pCrossTalkBar2away_corr[step1][step2][iCh]->Fill((energy[myChList2[ChList2Pos+2]] / IC_corr[myChList2[ChList2Pos+2]] )  / TotalEnergy_corr[1] );
	       }
	     }
	   } // closing loop of x_dut over center iCh
	 
	 //	    if (tot[iCh] > 0. ) // needs more troubleshooting for this, why does it make the distributions so broad
	 pTot_vs_Xpos[step1][step2][iCh]->Fill(x_dut, tot[iCh]/1.e3);
	 pTot_vs_Ypos[step1][step2][iCh]->Fill(y_dut, tot[iCh]/1.e3);
	 
	 // make efficiency plots - 1 if 0.9 - 3 * MIP energy, 0 otherwise
	 if (x_dut != 0 && y_dut != 0)
	   {
	     pXY_Edep[step1][step2][iCh]->Fill(x_dut, y_dut, tot[iCh]/1.e3);
	     if ((tot[iCh]/1.e3) >= 0.9 * MIP[iCh] && (tot[iCh]/1.e3) <= 3 * MIP[iCh] && x_dut > 0)
	       {
		 pEff_vs_Xpos[step1][step2][iCh]->Fill(x_dut, 1);
		 pEff_vs_Ypos[step1][step2][iCh]->Fill(y_dut, 1);
	       }
	   }
	 if (step1 > 1.e-2)
	   {
	     // fill the overlay plot with the x position from each channel   
	     pTot_vs_Xpos_overlay[step1][step2]->SetMaximum(60);
	     pTot_vs_Xpos_overlay[step1][step2]->SetMinimum(-5);
	     pTot_vs_Xpos_overlay[step1][step2]->Fill(x_dut, tot[iCh]/1.e3);
	     
	     pTot_vs_Ypos_overlay[step1][step2]->SetMaximum(60);
	     pTot_vs_Ypos_overlay[step1][step2]->SetMinimum(-5);
	     pTot_vs_Ypos_overlay[step1][step2]->Fill(y_dut, tot[iCh]/1.e3);
	   }
       } // closing channel loop
     
     for (int iBar = 0; iBar<NBARS; iBar++) // TO DO fix this for the new channel numbers
       {
	 int chBarUp = iBar*2+128;
	 int chBarDown = iBar*2+129;
         
	 if (tot[chBarUp]>130 && tot[chBarDown]>130 ) 
	   {
	     hCTR_UD[step1][step2][iBar]->Fill(time[chBarUp] - time[chBarDown]);                        
	   }
       }
   }// closing NEVENT loop
 
    TCanvas *cTots_scan[NSTEP1][NSTEP2][NCH];
    TCanvas *cTots_scan_correction[NSTEP1][NSTEP2][NCH];
    TCanvas *cXpos_scan[NSTEP1][NSTEP2][NCH];
    TCanvas *cYpos_scan[NSTEP1][NSTEP2][NCH];
    TCanvas *cEff_scanX[NSTEP1][NSTEP2][NCH];
    TCanvas *cEff_scanY[NSTEP1][NSTEP2][NCH];
    TCanvas *cCrossTalkOverlay[NSTEP1][NSTEP2][NCH];
    TCanvas *cCrossTalkOverlay_corr[NSTEP1][NSTEP2][NCH];
    TCanvas *cCrossTalk12away_corr[NSTEP1][NSTEP2][NCH];
    TCanvas *cXpos_over_scan;
    TCanvas *cYpos_over_scan;

    cXpos_over_scan = new TCanvas (Form("cXpos_over"), Form("cXpos_over"), 800, 400);
    cYpos_over_scan = new TCanvas (Form("cYpos_over"), Form("cYpos_over"), 800, 400);

    // now draw all the histograms, again have step 1, step 2, channel loop - step 1 is the overvoltage, step 2 is the threshold
    for (int iStep1 = 0; iStep1< NSTEP1; iStep1++)
    {
        for (int iStep2 = 0; iStep2< NSTEP2; iStep2++)
        {
            for (int iCh = 0; iCh< NCH; iCh++)
            {
	      // if not one of channels that we care about, skip it                                                                                                  
                if (std::find(std::begin(myChList), std::end(myChList), iCh) == std::end(myChList) ) continue;
      		//if (step1_vct.at(iStep1) !=6) continue;
		//if (step2_vct.at(iStep2) !=211606) continue;

                cTots_scan[iStep1][iStep2][iCh] = new TCanvas (Form("cTots_ch%.3d_step1_%.1f_step2_%.1f", iCh, step1_vct.at(iStep1), step2_vct.at(iStep2)), Form("cTots_ch%.3d_step1_%.1f_step2_%.1f", iCh, step1_vct.at(iStep1), step2_vct.at(iStep2)), 800, 400);
                cTots_scan[iStep1][iStep2][iCh]->cd();            
                hTot[step1_vct.at(iStep1)][step2_vct.at(iStep2)][iCh]->Rebin(REBIN_COEFF/4);
		hTot[step1_vct.at(iStep1)][step2_vct.at(iStep2)][iCh]->SetAxisRange(0,600,"X");
		hTot[step1_vct.at(iStep1)][step2_vct.at(iStep2)][iCh]->SetAxisRange(1,1500,"Y");
                hTot[step1_vct.at(iStep1)][step2_vct.at(iStep2)][iCh]->Draw();
                hTot[step1_vct.at(iStep1)][step2_vct.at(iStep2)][iCh]->GetXaxis()->SetTitle("tot [ns]");
                hTot[step1_vct.at(iStep1)][step2_vct.at(iStep2)][iCh]->GetYaxis()->SetTitle("Counts");
                gPad->SetLogy();
		// hTot_cut is for selecting the position of one bar and plotting the Landau peak
		// use the same binning for the plot before and after the cut on a single bar
		// different line color for the MIP peak after the bar cut (this is landau peak after choosing each bar)
		hTot_cut[step1_vct.at(iStep1)][step2_vct.at(iStep2)][iCh]->Rebin(REBIN_COEFF/4);
		hTot_cut[step1_vct.at(iStep1)][step2_vct.at(iStep2)][iCh]->SetLineColor(kGreen+2);
		hTot_cut[step1_vct.at(iStep1)][step2_vct.at(iStep2)][iCh]->SetAxisRange(0,600,"X");
                hTot_cut[step1_vct.at(iStep1)][step2_vct.at(iStep2)][iCh]->SetAxisRange(1,1500,"Y");
		hTot_cut[step1_vct.at(iStep1)][step2_vct.at(iStep2)][iCh]->Draw("same");
		// do the Landau fit, and set the range over which the fit is done (reduced from 120-400 to 130-250
		TF1 * fitLandau = new TF1 ("fitLandau","landau",MIP[iCh]-20,MIP[iCh]+80);
		hTot_cut[step1_vct.at(iStep1)][step2_vct.at(iStep2)][iCh]->Fit(fitLandau, "QRL");
		std::cout << "Ch # " << iCh << " Landau fit normalization coeff: " << fitLandau->GetParameter(0) << " most  probable value: " << fitLandau->GetParameter(1) << " Lambda value: " << fitLandau->GetParameter(2) << std::endl;
		cTots_scan[iStep1][iStep2][iCh]->SaveAs(Form("hTot_ch%.3d_step1_%.1f_step2_%.1f.pdf", iCh, step1_vct.at(iStep1), step2_vct.at(iStep2)));


		cTots_scan_correction[iStep1][iStep2][iCh] = new TCanvas (Form("cTots_corr_ch%.3d_step1_%.1f_step2_%.1f", iCh, step1_vct.at(iStep1), step2_vct.at(iStep2)), Form("cTots_corr_ch%.3d_step1_%.1f_step2_%.1f", iCh, step1_vct.at(iStep1), step2_vct.at(iStep2)), 800, 400);
                cTots_scan_correction[iStep1][iStep2][iCh]->cd();
                hTot_correction[step1_vct.at(iStep1)][step2_vct.at(iStep2)][iCh]->Rebin(REBIN_COEFF/4);
		hTot_correction[step1_vct.at(iStep1)][step2_vct.at(iStep2)][iCh]->SetAxisRange(0,600,"X");
                hTot_correction[step1_vct.at(iStep1)][step2_vct.at(iStep2)][iCh]->SetAxisRange(1,1500,"Y");
                hTot_correction[step1_vct.at(iStep1)][step2_vct.at(iStep2)][iCh]->Draw();
                hTot_correction[step1_vct.at(iStep1)][step2_vct.at(iStep2)][iCh]->GetXaxis()->SetTitle("Energy");
                hTot_correction[step1_vct.at(iStep1)][step2_vct.at(iStep2)][iCh]->GetYaxis()->SetTitle("Counts");
                gPad->SetLogy();
                // hTot_cut is for selecting the position of one bar and plotting the Landau peak
                // use the same binning for the plot before and after the cut on a single bar                                                         
                // different line color for the MIP peak after the bar cut (this is landau peak after choosing each bar)                                                   
                hTot_cut_correction[step1_vct.at(iStep1)][step2_vct.at(iStep2)][iCh]->Rebin(REBIN_COEFF/4);
                hTot_cut_correction[step1_vct.at(iStep1)][step2_vct.at(iStep2)][iCh]->SetLineColor(kGreen+2);
		hTot_cut_correction[step1_vct.at(iStep1)][step2_vct.at(iStep2)][iCh]->SetAxisRange(0,600,"X");
                hTot_cut_correction[step1_vct.at(iStep1)][step2_vct.at(iStep2)][iCh]->SetAxisRange(1,1500,"Y");
                hTot_cut_correction[step1_vct.at(iStep1)][step2_vct.at(iStep2)][iCh]->Draw("same");
                // do the Landau fit, and set the range over which the fit is done (reduced from 120-400 to 130-220                              
		TF1 * fitLandau_corr = new TF1 ("fitLandau_corr","landau",MIP_corr[iCh]-3,MIP_corr[iCh]+10);
                hTot_cut_correction[step1_vct.at(iStep1)][step2_vct.at(iStep2)][iCh]->Fit(fitLandau_corr, "QRL");
		std::cout << "Ch # " << iCh << " Landau fit on corrected tot normalization coeff: " << fitLandau_corr->GetParameter(0) << " most  probable value: " << fitLandau_corr->GetParameter(1) << " Lambda value: " << fitLandau_corr->GetParameter(2) << std::endl;
		myfile << "Ch # " << iCh << " Landau fit on corrected tot normalization coeff: " << fitLandau_corr->GetParameter(0) << " most  probable value: " << fitLandau_corr->GetParameter(1) << " Lambda value: " << fitLandau_corr->GetParameter(2) << std::endl;

                cTots_scan_correction[iStep1][iStep2][iCh]->SaveAs(Form("hTot_corr_ch%.3d_step1_%.1f_step2_%.1f.pdf", iCh, step1_vct.at(iStep1), step2_vct.at(iStep2)));

		// overlaying the cross talk with the cut on a certain bar
		cCrossTalkOverlay[iStep1][iStep2][iCh] = new TCanvas (Form("cCrossTalkOverlay_ch%.3d_step1_%.1f_step2_%.1f", iCh, step1_vct.at(iStep1), step2_vct.at(iStep2)), Form("cCrossTalkOverlay_ch%.3d_step1_%.1f_step2_%.1f", iCh, step1_vct.at(iStep1), step2_vct.at(iStep2)), 800, 400);
		cCrossTalkOverlay[iStep1][iStep2][iCh]->cd();
		pCrossTalk[step1_vct.at(iStep1)][step2_vct.at(iStep2)][iCh]->Draw();
		pCrossTalkBar[step1_vct.at(iStep1)][step2_vct.at(iStep2)][iCh]->SetLineColor(kGreen+2);
		pCrossTalkBar[step1_vct.at(iStep1)][step2_vct.at(iStep2)][iCh]->Draw("same");
		pCrossTalk[step1_vct.at(iStep1)][step2_vct.at(iStep2)][iCh]->GetXaxis()->SetTitle("Fractional energy deposit");
		pCrossTalk[step1_vct.at(iStep1)][step2_vct.at(iStep2)][iCh]->GetYaxis()->SetTitle("Events");
		gPad->SetLogy();
		cCrossTalkOverlay[iStep1][iStep2][iCh]->SaveAs(Form("crosstalkOverlay_ch%.3d_step1_%.1f_step2_%.1f.pdf", iCh, step1_vct.at(iStep1), step2_vct.at(iStep2)));

		//cross talk overlay with corrected tot and energy
		cCrossTalkOverlay_corr[iStep1][iStep2][iCh] = new TCanvas (Form("cCrossTalkOverlay_corr_ch%.3d_step1_%.1f_step2_%.1f", iCh, step1_vct.at(iStep1), step2_vct.at(iStep2)), Form("cCrossTalkOverlay_corr_ch%.3d_step1_%.1f_step2_%.1f", iCh, step1_vct.at(iStep1), step2_vct.at(iStep2)), 800, 400);
                cCrossTalkOverlay_corr[iStep1][iStep2][iCh]->cd();
                pCrossTalk_corr[step1_vct.at(iStep1)][step2_vct.at(iStep2)][iCh]->Draw();
                pCrossTalkBar_corr[step1_vct.at(iStep1)][step2_vct.at(iStep2)][iCh]->SetLineColor(kGreen+2);
                pCrossTalkBar_corr[step1_vct.at(iStep1)][step2_vct.at(iStep2)][iCh]->Draw("same");
                pCrossTalk_corr[step1_vct.at(iStep1)][step2_vct.at(iStep2)][iCh]->GetXaxis()->SetTitle("Fractional energy deposit (corrected)");
		pCrossTalk_corr[step1_vct.at(iStep1)][step2_vct.at(iStep2)][iCh]->GetYaxis()->SetTitle("Events");
                gPad->SetLogy();
                cCrossTalkOverlay_corr[iStep1][iStep2][iCh]->SaveAs(Form("crosstalkOverlay_corr_ch%.3d_step1_%.1f_step2_%.1f.pdf", iCh, step1_vct.at(iStep1), step2_vct.at(iStep2)));

		int ChList1Pos = -100;
		int ChList2Pos = -100;
		for (int ListPos = 0; ListPos < 16; ListPos++ )
		  {
		    if (iCh == myChList1[ListPos]) ChList1Pos = ListPos;
		    if (iCh == myChList2[ListPos]) ChList2Pos = ListPos;
		  }
		if ((ChList1Pos > 1 && ChList1Pos < 14) || (ChList2Pos > 1 && ChList2Pos < 14)) { // make sure cross talk can be defined for at least one side of the array
		  cCrossTalk12away_corr[iStep1][iStep2][iCh] = new  TCanvas (Form("cCrossTalk12away%.3d_corr_step1_%.1f_step2_%.1f", iCh, step1_vct.at(iStep1), step2_vct.at(iStep2)), Form("cCrossTalk12away%.3d_corr_step1_%.1f_step2_%.1f", iCh, step1_vct.at(iStep1), step2_vct.at(iStep2)), 800, 400);
		  cCrossTalk12away_corr[iStep1][iStep2][iCh]->cd();
		  pCrossTalk_corr[step1_vct.at(iStep1)][step2_vct.at(iStep2)][iCh]->Draw();
		  pCrossTalkBar_corr[step1_vct.at(iStep1)][step2_vct.at(iStep2)][iCh]->SetLineColor(kGreen+2);
		  pCrossTalkBar_corr[step1_vct.at(iStep1)][step2_vct.at(iStep2)][iCh]->Draw("same");
		  if (ChList1Pos > -1) {
		    std::cout << "making plot with cross talk 1 and 2 away" << std::endl;
		    pCrossTalkBar1away_corr[step1_vct.at(iStep1)][step2_vct.at(iStep2)][myChList1[ChList1Pos-1]]->SetLineColor(kRed+1); // left                         
		    pCrossTalkBar1away_corr[step1_vct.at(iStep1)][step2_vct.at(iStep2)][myChList1[ChList1Pos+1]]->SetLineColor(kRed-7); // right
		    pCrossTalkBar1away_corr[step1_vct.at(iStep1)][step2_vct.at(iStep2)][myChList1[ChList1Pos-1]]->Draw("same");
		    pCrossTalkBar1away_corr[step1_vct.at(iStep1)][step2_vct.at(iStep2)][myChList1[ChList1Pos+1]]->Draw("same");
		    pCrossTalkBar2away_corr[step1_vct.at(iStep1)][step2_vct.at(iStep2)][myChList1[ChList1Pos-2]]->SetLineColor(kMagenta+3);
		    pCrossTalkBar2away_corr[step1_vct.at(iStep1)][step2_vct.at(iStep2)][myChList1[ChList1Pos+2]]->SetLineColor(kMagenta-4);
		    pCrossTalkBar2away_corr[step1_vct.at(iStep1)][step2_vct.at(iStep2)][myChList1[ChList1Pos-2]]->Draw("same");
		    pCrossTalkBar2away_corr[step1_vct.at(iStep1)][step2_vct.at(iStep2)][myChList1[ChList1Pos+2]]->Draw("same");
		  }
		  if (ChList2Pos > -1 ) {
		    pCrossTalkBar1away_corr[step1_vct.at(iStep1)][step2_vct.at(iStep2)][myChList2[ChList2Pos-1]]->SetLineColor(kRed+1); // left
		    pCrossTalkBar1away_corr[step1_vct.at(iStep1)][step2_vct.at(iStep2)][myChList2[ChList2Pos+1]]->SetLineColor(kRed-7); // right         
                    pCrossTalkBar1away_corr[step1_vct.at(iStep1)][step2_vct.at(iStep2)][myChList2[ChList2Pos-1]]->Draw("same");
                    pCrossTalkBar1away_corr[step1_vct.at(iStep1)][step2_vct.at(iStep2)][myChList2[ChList2Pos+1]]->Draw("same");
                    pCrossTalkBar2away_corr[step1_vct.at(iStep1)][step2_vct.at(iStep2)][myChList2[ChList2Pos-2]]->SetLineColor(kMagenta+3);
		    pCrossTalkBar2away_corr[step1_vct.at(iStep1)][step2_vct.at(iStep2)][myChList2[ChList2Pos+2]]->SetLineColor(kMagenta-4);
                    pCrossTalkBar2away_corr[step1_vct.at(iStep1)][step2_vct.at(iStep2)][myChList2[ChList2Pos-2]]->Draw("same");
                    pCrossTalkBar2away_corr[step1_vct.at(iStep1)][step2_vct.at(iStep2)][myChList2[ChList2Pos+2]]->Draw("same");
                  }
		  pCrossTalk_corr[step1_vct.at(iStep1)][step2_vct.at(iStep2)][iCh]->GetXaxis()->SetTitle("Fractional energy deposit (corrected)");
		  pCrossTalk_corr[step1_vct.at(iStep1)][step2_vct.at(iStep2)][iCh]->GetYaxis()->SetTitle("Events");
		  gPad->SetLogy();
		  cCrossTalk12away_corr[iStep1][iStep2][iCh]->SaveAs(Form("crosstalk12away%.3d_corr_step1_%.1f_step2_%.1f.pdf", iCh, step1_vct.at(iStep1), step2_vct.at(iStep2)));
		  TF1 * fitGausGreen = new TF1 ("fitGausGreen", "gaus", 0.6, 0.84 ); // for array 2
		  TF1 * fitGausRedr = new TF1 ("fitGausRedr", "gaus", 0.03, 0.25 );
		  TF1 * fitGausRedl = new TF1 ("fitGausRedl", "gaus", 0.03, 0.25 );
		  TF1 * fitGausPurpler = new TF1 ("fitGausPurpler", "gaus", 0.01, 0.08 );
		  TF1 * fitGausPurplel = new TF1 ("fitGausPurplel", "gaus", 0.01, 0.08 );
		  pCrossTalkBar_corr[step1_vct.at(iStep1)][step2_vct.at(iStep2)][iCh]->Fit(fitGausGreen, "QRL");
		  //		  pCrossTalkBar1away_corr[step1_vct.at(iStep1)][step2_vct.at(iStep2)][iChr1]->Fit(fitGausRedr, "QRL");
		  //		  pCrossTalkBar1away_corr[step1_vct.at(iStep1)][step2_vct.at(iStep2)][iChl1]->Fit(fitGausRedl, "QRL");
		  //		  pCrossTalkBar2away_corr[step1_vct.at(iStep1)][step2_vct.at(iStep2)][iChr2]->Fit(fitGausPurpler, "QRL");
		  //		  pCrossTalkBar2away_corr[step1_vct.at(iStep1)][step2_vct.at(iStep2)][iChl2]->Fit(fitGausPurplel, "QRL");
		  //		std::cout << "fractional energy ch " << iCh << " :" << fitGausGreen->GetParameter(1) << " fractional energy 1 away (l,r): " << fitGausRedl->GetParameter(1) << " and " << fitGausRedr->GetParameter(1) << " fractional energy 2 away (l,r): " << fitGausPurplel->GetParameter(1) << " and " << fitGausPurpler->GetParameter(1) << std::endl;
		}

		// making plots for the x position of the device under test (x_dut)
		cXpos_scan[iStep1][iStep2][iCh] = new TCanvas (Form("cXpos_ch%.3d_step1_%.1f_step2_%.1f", iCh, step1_vct.at(iStep1), step2_vct.at(iStep2)), Form("cXpos_ch%.3d_step1_%.1f_step2_%.1f", iCh, step1_vct.at(iStep1), step2_vct.at(iStep2)), 800, 400);
                cXpos_scan[iStep1][iStep2][iCh]->cd();
                pTot_vs_Xpos[step1_vct.at(iStep1)][step2_vct.at(iStep2)][iCh]->Rebin(REBIN_COEFF);
                pTot_vs_Xpos[step1_vct.at(iStep1)][step2_vct.at(iStep2)][iCh]->Draw();
                pTot_vs_Xpos[step1_vct.at(iStep1)][step2_vct.at(iStep2)][iCh]->GetXaxis()->SetTitle("X position");
                pTot_vs_Xpos[step1_vct.at(iStep1)][step2_vct.at(iStep2)][iCh]->GetYaxis()->SetTitle("tot [ns]");
		TF1 * fitGaus = new TF1 ("fitGaus", "gaus", center[iCh]-2.5 , center[iCh]+2.5 ); 
		pTot_vs_Xpos[step1_vct.at(iStep1)][step2_vct.at(iStep2)][iCh]->Fit(fitGaus, "QRL");
		std::cout << "xPos ch[" << iCh << "] = " << fitGaus->GetParameter(1) << " x position centered " << std::endl;
		myfile << "xPos ch[" << iCh << "] = " << fitGaus->GetParameter(1) << " x position centered" << std::endl;

                cXpos_scan[iStep1][iStep2][iCh]->SaveAs(Form("pTot_vs_Xpos_ch%.3d_step1_%.1f_step2_%.1f.pdf", iCh, step1_vct.at(iStep1), step2_vct.at(iStep2)));

		cYpos_scan[iStep1][iStep2][iCh] = new TCanvas (Form("cYpos_ch%.3d_step1_%.1f_step2_%.1f", iCh, step1_vct.at(iStep1), step2_vct.at(iStep2)), Form("cYpos_ch%.3d_step1_%.1f_step2_%.1f", iCh, step1_vct.at(iStep1), step2_vct.at(iStep2)), 800, 400);
                cYpos_scan[iStep1][iStep2][iCh]->cd();
                pTot_vs_Ypos[step1_vct.at(iStep1)][step2_vct.at(iStep2)][iCh]->Rebin(REBIN_COEFF);
		pTot_vs_Ypos[step1_vct.at(iStep1)][step2_vct.at(iStep2)][iCh]->Draw();
                pTot_vs_Ypos[step1_vct.at(iStep1)][step2_vct.at(iStep2)][iCh]->GetXaxis()->SetTitle("Y position");
                pTot_vs_Ypos[step1_vct.at(iStep1)][step2_vct.at(iStep2)][iCh]->GetYaxis()->SetTitle("tot [ns]");
                pTot_vs_Ypos[step1_vct.at(iStep1)][step2_vct.at(iStep2)][iCh]->Fit(fitGaus, "QRL");
		std::cout << "yPos ch[" << iCh << "] = " << fitGaus->GetParameter(1) << " y position centered " << std::endl;
		//                cYpos_scan[iStep1][iStep2][iCh]->SaveAs(Form("pTot_vs_Ypos_ch%.3d_step1_%.1f_step2_%.1f.pdf", iCh, step1_vct.at(iStep1), step2_vct.at(iStep2)));

		// efficiency plots
		cEff_scanX[iStep1][iStep2][iCh] = new TCanvas (Form("cEffx_ch%.3d_step1_%.1f_step2_%.1f", iCh, step1_vct.at(iStep1), step2_vct.at(iStep2)), Form("cEffx__ch%.3d_step1_%.1f_step2_%.1f", iCh, step1_vct.at(iStep1), step2_vct.at(iStep2)), 800, 400);
		cEff_scanX[iStep1][iStep2][iCh]->cd();
		pEff_vs_Xpos[step1_vct.at(iStep1)][step2_vct.at(iStep2)][iCh]->Rebin(2);
		pEff_vs_Xpos[step1_vct.at(iStep1)][step2_vct.at(iStep2)][iCh]->Draw();
		pEff_vs_Xpos[step1_vct.at(iStep1)][step2_vct.at(iStep2)][iCh]->GetXaxis()->SetTitle("X position");
		pEff_vs_Xpos[step1_vct.at(iStep1)][step2_vct.at(iStep2)][iCh]->GetYaxis()->SetTitle("Within MIP Peak Energy");
		cEff_scanX[iStep1][iStep2][iCh]->SaveAs(Form("cEffx_ch%.3d_step1_%.1f_step2_%.1f.pdf", iCh, step1_vct.at(iStep1), step2_vct.at(iStep2)));

		cEff_scanY[iStep1][iStep2][iCh] = new TCanvas (Form("cEffy_ch%.3d_step1_%.1f_step2_%.1f", iCh, step1_vct.at(iStep1), step2_vct.at(iStep2)), Form("cEffy__ch%.3d_step1_%.1f_step2_%.1f", iCh, step1_vct.at(iStep1), step2_vct.at(iStep2)), 800, 400);
                cEff_scanY[iStep1][iStep2][iCh]->cd();
                pEff_vs_Ypos[step1_vct.at(iStep1)][step2_vct.at(iStep2)][iCh]->Rebin(2);
                pEff_vs_Ypos[step1_vct.at(iStep1)][step2_vct.at(iStep2)][iCh]->Draw();
                pEff_vs_Ypos[step1_vct.at(iStep1)][step2_vct.at(iStep2)][iCh]->GetXaxis()->SetTitle("Y position");
                pEff_vs_Ypos[step1_vct.at(iStep1)][step2_vct.at(iStep2)][iCh]->GetYaxis()->SetTitle("Within MIP Peak Energy");
		//                cEff_scanY[iStep1][iStep2][iCh]->SaveAs(Form("cEffy_ch%.3d_step1_%.1f_step2_%.1f.pdf", iCh, step1_vct.at(iStep1), step2_vct.at(iStep2)));
		// need to find the center of the efficiency plot for one channel - use a specific fit function - this is done later in overlay plots in channel / bar loop

		// plots for overlay of x position 
		cXpos_over_scan->cd();
		pTot_vs_Xpos[step1_vct.at(iStep1)][step2_vct.at(iStep2)][iCh]->Draw("same");
		pTot_vs_Xpos[step1_vct.at(iStep1)][step2_vct.at(iStep2)][iCh]->GetYaxis()->SetRange(-5,70);
		pTot_vs_Xpos[step1_vct.at(iStep1)][step2_vct.at(iStep2)][iCh]->SetMarkerColor(6);
		pTot_vs_Xpos[step1_vct.at(iStep1)][step2_vct.at(iStep2)][iCh]->SetTitle("X position");
		pTot_vs_Xpos[step1_vct.at(iStep1)][step2_vct.at(iStep2)][iCh]->SetTitle("tot [ns]");

		cYpos_over_scan->cd();
                pTot_vs_Ypos[step1_vct.at(iStep1)][step2_vct.at(iStep2)][iCh]->Draw("same");
                pTot_vs_Ypos[step1_vct.at(iStep1)][step2_vct.at(iStep2)][iCh]->GetYaxis()->SetRange(-5,70);
                pTot_vs_Ypos[step1_vct.at(iStep1)][step2_vct.at(iStep2)][iCh]->SetMarkerColor(6);
                pTot_vs_Ypos[step1_vct.at(iStep1)][step2_vct.at(iStep2)][iCh]->SetTitle("Y position");
                pTot_vs_Ypos[step1_vct.at(iStep1)][step2_vct.at(iStep2)][iCh]->SetTitle("tot [ns]");
            }
          cXpos_over_scan->SaveAs(Form("pTot_vs_Xpos_overlay_step1_%.1f_step2_%.1f.pdf", step1_vct.at(iStep1), step2_vct.at(iStep2)));
	  cYpos_over_scan->SaveAs(Form("pTot_vs_Ypos_overlay_step1_%.1f_step2_%.1f.pdf", step1_vct.at(iStep1), step2_vct.at(iStep2)));
        }
    }
    
    
    int selStep1 = 0;
    int selStep2 = 0;
    
    std::cout<< "before making array plots" << std::endl;

    TCanvas * cBeamPosX = new TCanvas ("cBeamPosX", "cBeamPosX", 500, 500);
    cBeamPosX->cd();
    hPosX->Rebin(2);
    hPosX->Draw();
    hPosX->GetXaxis()->SetTitle("beam X pos [mm]");
    hPosX->GetYaxis()->SetTitle("Counts");        
    gPad->SetLogy();

    TCanvas * cBeamPosY = new TCanvas ("cBeamPosY", "cBeamPosY", 500, 500);
    cBeamPosY->cd();
    hPosY->Rebin(2);
    hPosY->Draw();
    hPosY->GetXaxis()->SetTitle("beam Y pos [mm]");
    hPosY->GetYaxis()->SetTitle("Counts");
    gPad->SetLogy();

    TCanvas * cArrayTots = new TCanvas ("cArrayTots", "cArrayTots", 1600, 300);
    std::cout<< "canvas created" << std::endl;
    cArrayTots->Divide(NBARS, 2);
    std::cout<< "canvas divided" << std::endl;
    for (int chId = 0; chId<NBARS*2; chId++)
      {
	std::cout<< "in channel loop with " << chId << " and ch# " << myChList[chId]  << std::endl;
        cArrayTots->cd(chId+1);
        hTot[step1_vct.at(selStep1)][step2_vct.at(selStep2)][myChList[chId]]->Draw();
	hTot[step1_vct.at(selStep1)][step2_vct.at(selStep2)][myChList[chId]]->GetXaxis()->SetRange(0,400);
        hTot[step1_vct.at(selStep1)][step2_vct.at(selStep2)][myChList[chId]]->GetXaxis()->SetTitle("tot [ns]");
        hTot[step1_vct.at(selStep1)][step2_vct.at(selStep2)][myChList[chId]]->GetYaxis()->SetTitle("Counts");
	if ( myChList[chId] == 142 )
	  {
	    cArrayTots->SaveAs(Form("hTot_array_chId%.3d.pdf", myChList[chId]));
	    break;
	  }
      }
    
    TCanvas * cArrayTimes = new TCanvas ("cArrayTimes", "cArrayTimes", 1600, 300);
    cArrayTimes->Divide(NBARS, 2);
    for (int chId = 0; chId<NBARS*2; chId++)
      {
        cArrayTimes->cd(chId+1);
        hTime[step1_vct.at(selStep1)][step2_vct.at(selStep2)][myChList[chId]]->Draw();
	hTime[step1_vct.at(selStep1)][step2_vct.at(selStep2)][myChList[chId]]->GetXaxis()->SetRange(-250,-100);
        hTime[step1_vct.at(selStep1)][step2_vct.at(selStep2)][myChList[chId]]->GetXaxis()->SetTitle("Time [ps]");
        hTime[step1_vct.at(selStep1)][step2_vct.at(selStep2)][myChList[chId]]->GetYaxis()->SetTitle("Counts");
	if ( myChList[chId] == 142 )
	  {
	    cArrayTimes->SaveAs(Form("hTime_array_chId%.3d.pdf", myChList[chId]));
	    break;
	  }
      }
    
    
    TCanvas * cArrayCTR_UD = new TCanvas ("cArrayCTR_UD", "cArrayCTR_UD", 1600, 300);
    cArrayCTR_UD->Divide(NBARS, 1);
    for (int barId = 0; barId<NBARS; barId++)
      {
        cArrayCTR_UD->cd(barId+1);
        hCTR_UD[step1_vct.at(selStep1)][step2_vct.at(selStep2)][barId]->Rebin(REBIN_COEFF);
        hCTR_UD[step1_vct.at(selStep1)][step2_vct.at(selStep2)][barId]->Draw();
        hCTR_UD[step1_vct.at(selStep1)][step2_vct.at(selStep2)][barId]->GetXaxis()->SetTitle("t_{up} - t_{down} [ps]");
        hCTR_UD[step1_vct.at(selStep1)][step2_vct.at(selStep2)][barId]->GetYaxis()->SetTitle("Counts");
        
        float mean = hCTR_UD[step1_vct.at(selStep1)][step2_vct.at(selStep2)][barId]->GetMean();
        float rms  = hCTR_UD[step1_vct.at(selStep1)][step2_vct.at(selStep2)][barId]->GetRMS();
        
        TF1 * fitGaus = new TF1 ("fitGaus", "gaus", mean-rms , mean+rms);
        hCTR_UD[step1_vct.at(selStep1)][step2_vct.at(selStep2)][barId]->Fit(fitGaus, "QRL");
	if (barId == NBARS - 1)
	  {
	    cArrayCTR_UD->SaveAs(Form("hCTR_UD_array_barId%.3d.pdf", barId));
	    break;
	  }
        std::cout << "CTR_UD [" << barId << "] = " << fitGaus->GetParameter(2) << " ps --> sigma_bar = " << fitGaus->GetParameter(2)/2 << " ps " << std::endl;
      }

    TCanvas * cArrayEdepXY = new TCanvas ("cArrayEdepXY", "cArrayEdepXY", 1600, 300);
    cArrayEdepXY->Divide(NBARS, 2);
    for (int chId = 0; chId<NBARS*2; chId++)
      {
	cArrayEdepXY->cd(chId+1);
        pXY_Edep[step1_vct.at(selStep1)][step2_vct.at(selStep2)][myChList[chId]]->Draw("COLZ");
        pXY_Edep[step1_vct.at(selStep1)][step2_vct.at(selStep2)][myChList[chId]]->GetXaxis()->SetTitle("beam X pos [mm]");
        pXY_Edep[step1_vct.at(selStep1)][step2_vct.at(selStep2)][myChList[chId]]->GetYaxis()->SetTitle("beam Y pos [mm]");
        pXY_Edep[step1_vct.at(selStep1)][step2_vct.at(selStep2)][myChList[chId]]->GetZaxis()->SetTitle("Mean Tot [ns]");
      }
    cArrayEdepXY->SaveAs(Form("XY_scatter_all.pdf"));

    TCanvas * cArrayEffx = new TCanvas("cArrayEffx","cArrayEffx", 1600, 300);
    TCanvas * cArrayEffy = new TCanvas("cArrayEffy","cArrayEffy", 1600, 300);
    cArrayEffx->Divide(NBARS,2);
    cArrayEffy->Divide(NBARS,2);
    for (int chId = 0; chId<NBARS*2; chId++)
      {
	// normalize the efficiency plots
	pEff_vs_Xpos[step1_vct.at(selStep1)][step2_vct.at(selStep2)][myChList[chId]]->Rebin(2);
	pEff_vs_Ypos[step1_vct.at(selStep1)][step2_vct.at(selStep2)][myChList[chId]]->Rebin(2);
        for (int iBin = 0; iBin < 200; iBin++)
	  {
            if (hPosX->GetBinContent(iBin+1)>0) pEff_vs_Xpos[step1_vct.at(selStep1)][step2_vct.at(selStep2)][myChList[chId]]->SetBinContent(iBin+1, pEff_vs_Xpos[step1_vct.at(selStep1)][step2_vct.at(selStep2)][myChList[chId]]->GetBinContent(iBin+1)/hPosX->GetBinContent(iBin+1) );
	    if (hPosY->GetBinContent(iBin+1)>0) pEff_vs_Ypos[step1_vct.at(selStep1)][step2_vct.at(selStep2)][myChList[chId]]->SetBinContent(iBin+1, pEff_vs_Ypos[step1_vct.at(selStep1)][step2_vct.at(selStep2)][myChList[chId]]->GetBinContent(iBin+1)/hPosY->GetBinContent(iBin+1) );
	  }
	cArrayEffx->cd(chId+1);
        pEff_vs_Xpos[step1_vct.at(selStep1)][step2_vct.at(selStep2)][myChList[chId]]->SetAxisRange(0,1,"Y"); // set Y range
	pEff_vs_Xpos[step1_vct.at(selStep1)][step2_vct.at(selStep2)][myChList[chId]]->SetAxisRange(0,40,"X"); // set X range
        pEff_vs_Xpos[step1_vct.at(selStep1)][step2_vct.at(selStep2)][myChList[chId]]->Draw();
	pEff_vs_Xpos[step1_vct.at(selStep1)][step2_vct.at(selStep2)][myChList[chId]]->SetStats(0); // no stats box for the efficiency plots
        pEff_vs_Xpos[step1_vct.at(selStep1)][step2_vct.at(selStep2)][myChList[chId]]->SetTitle("X position");
        pEff_vs_Xpos[step1_vct.at(selStep1)][step2_vct.at(selStep2)][myChList[chId]]->SetTitle("Efficiency within MIP Peak Energy");

	cArrayEffy->cd(chId+1);
	pEff_vs_Ypos[step1_vct.at(selStep1)][step2_vct.at(selStep2)][myChList[chId]]->SetAxisRange(0,1,"Y"); // set Y range
        pEff_vs_Ypos[step1_vct.at(selStep1)][step2_vct.at(selStep2)][myChList[chId]]->SetAxisRange(0,30,"X"); // set X range
        pEff_vs_Ypos[step1_vct.at(selStep1)][step2_vct.at(selStep2)][myChList[chId]]->Draw();
        pEff_vs_Ypos[step1_vct.at(selStep1)][step2_vct.at(selStep2)][myChList[chId]]->SetStats(0); // no stats box for the efficiency plots  
        pEff_vs_Ypos[step1_vct.at(selStep1)][step2_vct.at(selStep2)][myChList[chId]]->SetTitle("Y position");
        pEff_vs_Ypos[step1_vct.at(selStep1)][step2_vct.at(selStep2)][myChList[chId]]->SetTitle("Efficiency within MIP Peak Energy");
      }

    TF1 * fitBarPos = new TF1 ("fitBarPos", fitBarEffErr, 2, 32, 5);
    fitBarPos->SetParameters(5., 3., 0.01, 0.8, 0.2);
    fitBarPos->SetNpx(5000);
    fitBarPos->SetParLimits(0, minXpos, maxXpos);
    fitBarPos->SetParLimits(1, 2, 4);
    fitBarPos->SetParLimits(3, 0.2, 1.);

    TCanvas *cArrayEffOverlayX = new TCanvas("cArrayEffOverlayX","cArrayEffOverlayX",800,400);
    cArrayEffOverlayX->cd();
    // look at just one channel first
    pEff_vs_Xpos[step1_vct.at(selStep1)][step2_vct.at(selStep2)][myChList[0]]->Draw();
    pEff_vs_Xpos[step1_vct.at(selStep1)][step2_vct.at(selStep2)][myChList[0]]->SetTitle("X position");
    pEff_vs_Xpos[step1_vct.at(selStep1)][step2_vct.at(selStep2)][myChList[0]]->SetTitle("Efficiency within MIP Peak Energy");

    // then go over all the bars
    for (int chId = 0; chId<NBARS; chId++)
      {
	pEff_vs_Xpos[step1_vct.at(selStep1)][step2_vct.at(selStep2)][myChList[chId]]->SetLineColor(chId+1);
	fitBarPos->SetLineColor(chId+1);
	pEff_vs_Xpos[step1_vct.at(selStep1)][step2_vct.at(selStep2)][myChList[chId]]->Draw("same");
	fitBarPos->SetParameter(0, center[myChList[chId]]);
	pEff_vs_Xpos[step1_vct.at(selStep1)][step2_vct.at(selStep2)][myChList[chId]]->Fit(fitBarPos,"QR");

        float posBar   = fitBarPos->GetParameter(0);
        float widthBar = fitBarPos->GetParameter(1);
        float effBar   = fitBarPos->GetParameter(3);
        float trackRes = fitBarPos->GetParameter(4);
	std::cout << "posBar[" << chId << "] = " << posBar << " :: width = " << widthBar << " :: effBar = " << effBar << " :: trackRes = " << trackRes << std::endl; // printing for x position of bars currently
	if (chId == NBARS - 1)
	  {
	    cArrayEffx->SaveAs(Form("EfficiencyX_array_barId%.3d.pdf", chId));
	    cArrayEffOverlayX->SaveAs(Form("EfficiencyX_array_overlay_barId%.3d.pdf", chId));
	  }
      }

    TCanvas *cArrayEffOverlayY = new TCanvas("cArrayEffOverlayY","cArrayEffOverlayY",800,400);
    cArrayEffOverlayY->cd();
    pEff_vs_Ypos[step1_vct.at(selStep1)][step2_vct.at(selStep2)][myChList[0]]->Draw();
    pEff_vs_Ypos[step1_vct.at(selStep1)][step2_vct.at(selStep2)][myChList[0]]->SetTitle("Y position");
    pEff_vs_Ypos[step1_vct.at(selStep1)][step2_vct.at(selStep2)][myChList[0]]->SetTitle("Efficiency within MIP Peak Energy");
    for (int chId = 0; chId<NBARS; chId++)
      {
        pEff_vs_Ypos[step1_vct.at(selStep1)][step2_vct.at(selStep2)][myChList[chId]]->SetLineColor(chId+1);
	fitBarPos->SetLineColor(chId+1);
        pEff_vs_Ypos[step1_vct.at(selStep1)][step2_vct.at(selStep2)][myChList[chId]]->Draw("same");
        fitBarPos->SetParameter(0, 4 + chId*3.);
        for (int i = 0; i< 3; i++)
          {
            pEff_vs_Ypos[step1_vct.at(selStep1)][step2_vct.at(selStep2)][myChList[chId]]->Fit(fitBarPos,"QR");                                                                                         
          }
        float posBar   = fitBarPos->GetParameter(0);
        float widthBar = fitBarPos->GetParameter(1);
        float effBar   = fitBarPos->GetParameter(3);
        float trackRes = fitBarPos->GetParameter(4);
	std::cout << "posBar[" << chId << "] = " << posBar << " :: width = " << widthBar << " :: effBar = " << effBar << " :: trackRes = " << trackRes << std::endl; // printing for x position of bars currently                                                                                                                                                                                                           
	if (chId == NBARS - 1)
	  {
	    cArrayEffy->SaveAs(Form("EfficiencyY_array_barId%.3d.pdf", chId));
	    cArrayEffOverlayY->SaveAs(Form("EfficiencyY_array_overlay_barId%.3d.pdf", chId));
	  }
      }

    /*
    int selBar = 5;
    int selBarCh = selBar*2+128;

    TCanvas * cSingleScatterEdep = new TCanvas ("cSingleScatterEdep", "cSingleScatterEdep", 500, 500);
    cSingleScatterEdep->cd();
    pXY_Edep[step1_vct.at(selStep1)][step2_vct.at(selStep2)][selBarCh]->SetStats(0);
    pXY_Edep[step1_vct.at(selStep1)][step2_vct.at(selStep2)][selBarCh]->Draw("COLZ");
    pXY_Edep[step1_vct.at(selStep1)][step2_vct.at(selStep2)][selBarCh]->GetXaxis()->SetTitle("beam X pos [mm]");
    pXY_Edep[step1_vct.at(selStep1)][step2_vct.at(selStep2)][selBarCh]->GetYaxis()->SetTitle("beam Y pos [mm]");
    pXY_Edep[step1_vct.at(selStep1)][step2_vct.at(selStep2)][selBarCh]->GetZaxis()->SetTitle("Mean Tot [ns]");
    pXY_Edep[step1_vct.at(selStep1)][step2_vct.at(selStep2)][selBarCh]->GetXaxis()->SetRangeUser(0, 30);
    pXY_Edep[step1_vct.at(selStep1)][step2_vct.at(selStep2)][selBarCh]->GetYaxis()->SetRangeUser(10, 40);
    cSingleScatterEdep->SaveAs(Form("XY_scatter.pdf"));
    */

    TF1 * fitTestMyFunc = new TF1 ("fitTestMyFunc", fitBarEffErr, -20, 20, 5);
    fitTestMyFunc->SetParameters(10., 3., 0.0, 0.8, 0.2);
    fitTestMyFunc->SetLineColor(kBlack);
    fitTestMyFunc->SetNpx(10000);
    fitTestMyFunc->Draw("same");
    
    /*    
    //save plots or histos/graphs to output file
    TFile * outputFile = new TFile(Form("./output/output_run_%.3d.root", firstRun), "RECREATE");
    outputFile->cd();    
    
    outputFile->Write();
    outputFile->Close();    
    */
    
    std::cout << "end of program" << std::endl;
    theApp->Run();
    myfile.close();

}
