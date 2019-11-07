// to compile: g++ -Wall -o plotTOFHIR_to_trigger plotTOFHIR_to_trigger.C `root-config --cflags --glibs` -lSpectrum

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
    
    //    std::string data_path = "../data/RecoData/v1/RecoWithTracks";
    std::string data_path = "/eos/cms/store/group/dpg_mtd/comm_mtd/TB/MTDTB_FNAL_Jun2019/TOFHIR/RecoData/RecoWithTracks/v2";
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
    long long   chTime[400];    
    float	chtot[400];
            
    TChain * tree = new TChain("data", "data");
    
    // going over each of the runs listed as command line arguments
    for (int iRun = firstRun; iRun <= lastRun; iRun++)
      {
	//if ( std::find(bad_runs.begin(), bad_runs.end(), iRun) != bad_runs.end()) continue;       
	tree->Add( Form("%s/run%.5d_events_withTrack.root", data_path.c_str(), iRun) );
	//tree->Add( Form("%s/run%.5d_events.root", data_path.c_str(), iRun) ); 
	std::cout << "adding run: " << iRun << std::endl;
      }

    // extract information from the ROOT tree
    tree->SetBranchStatus("*",0);
    tree->SetBranchAddress("step1", &step1);
    tree->SetBranchAddress("step2", &step2);    
    tree->SetBranchAddress("x_dut", &x_dut);
    tree->SetBranchAddress("y_dut", &y_dut);
    tree->SetBranchAddress("ntracks", &ntracks);    
    tree->SetBranchAddress("chtot", &chtot);
    tree->SetBranchAddress("chTime", &chTime);
    
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
     
    for (int iStep1 = 0; iStep1<NSTEP1; iStep1++)
      {
        std::cout << "step1 (" << iStep1+1 << "/" << NSTEP1 << "): " << step1_vct.at(iStep1)  << std::endl;        
      }
    for (int iStep2 = 0; iStep2<NSTEP2; iStep2++)
      {
        std::cout << "step2 (" << iStep2+1 << "/" << NSTEP2 << "): " << step2_vct.at(iStep2)  << std::endl;        
      }    

    // define histos    
    double minTot = 0;
    double maxTot = 700;
    
    double minTime = -200000;
    double maxTime = 200000;
    
    double minXpos = 0;
    double maxXpos = 40;

    double minYpos = -10;
    double maxYpos = 50;

    int REBIN_COEFF = 32;
    
    int NCH = 400;

    // these are from the board mapping on the google sheet tab 
    int myChList[] = {128, 130, 132, 134, 136, 138, 140, 142, 129, 131, 133, 135, 137, 139, 141, 143}; // this is for the first array
    //int myChList[] = {128, 130, 132, 134, 136, 138, 140, 142}; // this is for the second array
    int NBARS = 8;
    double center[NCH];
    double MIP[NCH];
    double IC[NCH];
    double MIP_corr[NCH];
    double IC_corr[NCH];

    // centers for first array when both sides of the bars are read out
    /*
    center[128] = 5;
    center[129] = 5;
    center[130] = 0;
    center[131] = 8;
    center[132] = 11;
    center[133] = 11;
    center[134] = 14;
    center[135] = 14;
    center[136] = 17;
    center[137] = 17;
    center[138] = 20;
    center[139] = 20;
    center[140] = 23;
    center[141] = 23;
    center[142] = 26;
    center[143] = 26;
    */
    // centers for second array when only one side of the bar is read out    
    // /*
    center[128] = 7;           
    center[132] = 11.18;                                                                         
    center[134] = 14.32;                                        
    center[136] = 17.49;                                    
    center[138] = 20.62;                              
    center[140] = 23.95;                                             
    center[142] = 26.79;                              
    // */

    // list MIP peak energies based off of Landau fit, this is for v2 recomstruction
    MIP[128] = 144.993;
    MIP[132] = 141.019;
    MIP[134] = 144.835;
    MIP[136] = 140.924;
    MIP[138] = 152.04;
    MIP[140] = 169.621;
    MIP[142] = 152.091;
    
    // calculate intercallibration coefficients
    double avgMIP = (MIP[128] + MIP[132] + MIP[134] + MIP[136] + MIP[138] + MIP[140] + MIP[142]) / 7;
    IC[128] = MIP[128] / avgMIP;
    IC[132] = MIP[132] / avgMIP;
    IC[134] = MIP[134] / avgMIP;
    IC[136] = MIP[136] / avgMIP;
    IC[138] = MIP[138] / avgMIP;
    IC[140] = MIP[140] / avgMIP;
    IC[142] = MIP[142] / avgMIP;
    std::cout << "128: " << IC[128] << " 132: " << IC[132] << " 134: " << IC[134] << " 136: " << IC[136] << " 138: " << IC[138] << " 140: " << IC[140] << " 142: " << IC[142] << std::endl;

    // list MIP peak energies based off of Landau fit, this is for v2 recomstruction - using Landau to CORRECTED TOT
    MIP_corr[128] = 120.565;
    MIP_corr[132] = 113.45;
    MIP_corr[134] = 121.412;
    MIP_corr[136] = 113.446;
    MIP_corr[138] = 137.645;
    MIP_corr[140] = 185.344;
    MIP_corr[142] = 137.809;

    // calculate intercallibration coefficients   
    double avgMIP_corr = (MIP_corr[128] + MIP_corr[132] + MIP_corr[134] + MIP_corr[136] + MIP_corr[138] + MIP_corr[140] + MIP_corr[142]) / 7;
    IC_corr[128] = MIP_corr[128] / avgMIP_corr;
    IC_corr[132] = MIP_corr[132] / avgMIP_corr;
    IC_corr[134] = MIP_corr[134] / avgMIP_corr;
    IC_corr[136] = MIP_corr[136] / avgMIP_corr;
    IC_corr[138] = MIP_corr[138] / avgMIP_corr;
    IC_corr[140] = MIP_corr[140] / avgMIP_corr;
    IC_corr[142] = MIP_corr[142] / avgMIP_corr;
    std::cout << "corrected 128: " << IC_corr[128] << " 132: " << IC_corr[132] << " 134: " << IC_corr[134] << " 136: " << IC_corr[136] << " 138: " << IC_corr[138] << " 140: " << IC_corr[140] << " 142: " << IC_corr[142] << std::endl;

    // declare the histograms, these will be filled in the channel loop    
    std::map<float, std::map<float, std::map<int, TH1F * > > > hTot;
    std::map<float, std::map<float, std::map<int, TH1F * > > > hTot_correction;
    std::map<float, std::map<float, std::map<int, TH1F * > > > hTot_cut;
    std::map<float, std::map<float, std::map<int, TH1F * > > > hTot_cut_correction;
    std::map<float, std::map<float, std::map<int, TH1F * > > > hTime;
    std::map<float, std::map<float, std::map<int, TProfile * > > > pTot_vs_Xpos;
    std::map<float, std::map<float, std::map<int, TProfile * > > > pTot_vs_Ypos;
    std::map<float, std::map<float, std::map<int, TH1F * > > > pEff_vs_Xpos; 
    std::map<float, std::map<float, std::map<int, TH2F * > > > pXpos_Ypos_Tot;
    std::map<float, std::map<float, std::map<int, TH1F * > > > hCTR_UD;
    std::map<float, std::map<float, TProfile * > > pTot_vs_Xpos_overlay;
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
    
    // step 1, step 2, and channel listing loops. Histograms are defined inside the loops
    for (int iStep1 = 0; iStep1< NSTEP1; iStep1++)
      {
        for (int iStep2 = 0; iStep2< NSTEP2; iStep2++)
	  {
            for (int iCh = 0; iCh < NCH; iCh++)
	      {       
		hTot[step1_vct.at(iStep1)][step2_vct.at(iStep2)][iCh] = new TH1F (Form("hTot_ch%.3d_step1vct_%.1f_step2vct_%.1f_step1_%i_step2_%i", iCh, step1_vct.at(iStep1), step2_vct.at(iStep2), iStep1, iStep2), Form("Time over threshold, ch%.3d, step1_%.1f, step2_%.1f", iCh, step1_vct.at(iStep1), step2_vct.at(iStep2)), 4000, minTot, maxTot );

		hTot_correction[step1_vct.at(iStep1)][step2_vct.at(iStep2)][iCh] = new TH1F (Form("hTot_corr_ch%.3d_step1vct_%.1f_step2vct_%.1f_step1_%i_step2_%i", iCh, step1_vct.at(iStep1), step2_vct.at(iStep2), iStep1, iStep2), Form("Time over threshold corrected, ch%.3d, step1_%.1f, step2_%.1f", iCh, step1_vct.at(iStep1), step2_vct.at(iStep2)), 4000, minTot, maxTot );

		hTot_cut[step1_vct.at(iStep1)][step2_vct.at(iStep2)][iCh] = new TH1F (Form("hTot_cut_ch%.3d_step1vct_%.1f_step2vct_%.1f_step1_%i_step2_%i", iCh, step1_vct.at(iStep1), step2_vct.at(iStep2), iStep1, iStep2), Form("Time over threshold with x cut, ch%.3d, step1_%.1f, step2_%.1f", iCh, step1_vct.at(iStep1), step2_vct.at(iStep2)), 4000, minTot, maxTot );
                
		hTot_cut_correction[step1_vct.at(iStep1)][step2_vct.at(iStep2)][iCh] = new TH1F (Form("hTot_cut_corr_ch%.3d_step1vct_%.1f_step2vct_%.1f_step1_%i_step2_%i", iCh, step1_vct.at(iStep1), step2_vct.at(iStep2), iStep1, iStep2), Form("Time over threshold corrected with x cut, ch%.3d, step1_%.1f, step2_%.1f", iCh, step1_vct.at(iStep1), step2_vct.at(iStep2)), 4000, minTot, maxTot );

		hTime[step1_vct.at(iStep1)][step2_vct.at(iStep2)][iCh] = new TH1F (Form("hTime_ch%.3d_step1vct_%.1f_step2vct_%.1f_step1_%i_step2_%i", iCh, step1_vct.at(iStep1), step2_vct.at(iStep2), iStep1, iStep2), Form("Time Hist, ch%.3d, step1_%.1f, step2_%.1f", iCh, step1_vct.at(iStep1), step2_vct.at(iStep2)), 4000, minTime, maxTime );

		pTot_vs_Xpos[step1_vct.at(iStep1)][step2_vct.at(iStep2)][iCh] = new TProfile (Form("pTot_vs_Xpos_ch%.3d_step1vct_%.1f_step2vct_%.1f_step1_%i_step2_%i", iCh, step1_vct.at(iStep1), step2_vct.at(iStep2), iStep1, iStep2), Form("ToT vs. X pos, ch%.3d, step1_%.1f, step2_%.1f", iCh, step1_vct.at(iStep1), step2_vct.at(iStep2)), 4000, minXpos, maxXpos );

		pEff_vs_Xpos[step1_vct.at(iStep1)][step2_vct.at(iStep2)][iCh] = new TH1F (Form("pEff_vs_Xpos_ch%.3d_step1vct_%.1f_step2vct_%.1f_step1_%i_step2_%i", iCh, step1_vct.at(iStep1), step2_vct.at(iStep2), iStep1, iStep2), Form("Efficiency vs. X pos, ch%.3d, step1_%.1f, step2_%.1f", iCh, step1_vct.at(iStep1), step2_vct.at(iStep2)), 400, 0, maxXpos);

		pXpos_Ypos_Tot[step1_vct.at(iStep1)][step2_vct.at(iStep2)][iCh] = new TH2F (Form("pXpos_Ypos_Tot_ch%.3d_step1vct_%.1f_step2vct_%.1f_step1_%i_step2_%i", iCh, step1_vct.at(iStep1), step2_vct.at(iStep2), iStep1, iStep2), Form("X and Y pos vs. ToT, ch%.3d, step1_%.1f, step2_%.1f", iCh, step1_vct.at(iStep1), step2_vct.at(iStep2)), 40, minXpos, maxXpos, 40, minYpos, maxYpos );

		pXY_Edep[step1_vct.at(iStep1)][step2_vct.at(iStep2)][iCh] = new TProfile2D (Form("pXY_Edep_ch%.3d_step1_%.1f_step2_%.1f", iCh, step1_vct.at(iStep1), step2_vct.at(iStep2)), Form("XY Energy dep ch%.3d_step1_%.1f_step2_%.1f", iCh, step1_vct.at(iStep1), step2_vct.at(iStep2)), 400, minXpos, maxXpos, 400, minXpos, maxXpos );

		pCrossTalk[step1_vct.at(iStep1)][step2_vct.at(iStep2)][iCh] = new TH1F (Form("pCrossTalk_ch%.3d_step1_%.1f_step2_%.1f", iCh, step1_vct.at(iStep1), step2_vct.at(iStep2)), Form("Cross Talk ch%.3d_step1_%.1f_step2_%.1f", iCh, step1_vct.at(iStep1), step2_vct.at(iStep2)), 400, 0, 1.1);

		pCrossTalkBar[step1_vct.at(iStep1)][step2_vct.at(iStep2)][iCh] = new TH1F (Form("pCrossTalkBar_ch%.3d_step1_%.1f_step2_%.1f", iCh, step1_vct.at(iStep1), step2_vct.at(iStep2)), Form("Cross Talk (with bar cut) ch%.3d_step1_%.1f_step2_%.1f", iCh, step1_vct.at(iStep1), step2_vct.at(iStep2)), 400, 0, 1.1);

		pCrossTalkBar1away[step1_vct.at(iStep1)][step2_vct.at(iStep2)][iCh] = new TH1F (Form("pCrossTalkBar1away_ch%.3d_step1_%.1f_step2_%.1f", iCh, step1_vct.at(iStep1), step2_vct.at(iStep2)), Form("Cross Talk 1 away (with bar cut) ch%.3d_step1_%.1f_step2_%.1f", iCh, step1_vct.at(iStep1), step2_vct.at(iStep2)), 400, 0, 1.1);

		pCrossTalkBar2away[step1_vct.at(iStep1)][step2_vct.at(iStep2)][iCh] = new TH1F (Form("pCrossTalkBar2away_ch%.3d_step1_%.1f_step2_%.1f", iCh, step1_vct.at(iStep1), step2_vct.at(iStep2)), Form("Cross Talk 2 away (with bar cut) ch%.3d_step1_%.1f_step2_%.1f", iCh, step1_vct.at(iStep1), step2_vct.at(iStep2)), 400, 0, 1.1);


		//cross talk with corrected tot
		pCrossTalk_corr[step1_vct.at(iStep1)][step2_vct.at(iStep2)][iCh] = new TH1F (Form("pCrossTalk_corr_ch%.3d_step1_%.1f_step2_%.1f", iCh, step1_vct.at(iStep1), step2_vct.at(iStep2)), Form("Cross Talk (corrected tot) ch%.3d_step1_%.1f_step2_%.1f", iCh, step1_vct.at(iStep1), step2_vct.at(iStep2)), 400, 0, 1.1);

		pCrossTalkBar_corr[step1_vct.at(iStep1)][step2_vct.at(iStep2)][iCh] = new TH1F (Form("pCrossTalkBar_corr_ch%.3d_step1_%.1f_step2_%.1f", iCh, step1_vct.at(iStep1), step2_vct.at(iStep2)), Form("Cross Talk (with bar cut, corrected tot) ch%.3d_step1_%.1f_step2_%.1f", iCh, step1_vct.at(iStep1), step2_vct.at(iStep2)), 400, 0, 1.1);

		pCrossTalkBar1away_corr[step1_vct.at(iStep1)][step2_vct.at(iStep2)][iCh] = new TH1F (Form("pCrossTalkBar1away_corr_ch%.3d_step1_%.1f_step2_%.1f", iCh, step1_vct.at(iStep1), step2_vct.at(iStep2)), Form("Cross Talk 1 away (with bar cut, corrected tot) ch%.3d_step1_%.1f_step2_%.1f", iCh, step1_vct.at(iStep1), step2_vct.at(iStep2)), 400, 0, 1.1);

                pCrossTalkBar2away_corr[step1_vct.at(iStep1)][step2_vct.at(iStep2)][iCh] = new TH1F (Form("pCrossTalkBar2away_corr_ch%.3d_step1_%.1f_step2_%.1f", iCh, step1_vct.at(iStep1), step2_vct.at(iStep2)), Form("Cross Talk 2 away (with bar cut, corrected tot) ch%.3d_step1_%.1f_step2_%.1f", iCh, step1_vct.at(iStep1), step2_vct.at(iStep2)), 400, 0, 1.1);

		//pCrossTalkOverlay[step1_vct.at(iStep1)][step2_vct.at(iStep2)][iCh] = new TH1F (Form("pCrossTalkOverlay_ch%.3d_step1_%.1f_step2_%.1f", iCh, step1_vct.at(iStep1), step2_vct.at(iStep2)), Form("Cross Talk Overlayed ch%.3d_step1_%.1f_step2_%.1f", iCh, step1_vct.at(iStep1), step2_vct.at(iStep2)), 400, 0, 1.1);
	      }

	    pTot_vs_Xpos_overlay[step1_vct.at(iStep1)][step2_vct.at(iStep2)] = new TProfile (Form("pTot_vs_Xpos_over_step1vct_%.1f_step2vct_%.1f_step1_%i_step2_%i", step1_vct.at(iStep1), step2_vct.at(iStep2), iStep1, iStep2), Form("ToT vs. X pos overlay, step1_%.1f, step2_%.1f", step1_vct.at(iStep1), step2_vct.at(iStep2)), 4000, minXpos, maxXpos );

	    // bar loop (outside of channel loop) to define time resolution for a single bar
            for (int iBar = 0; iBar < NBARS; iBar++)
	      {
		hCTR_UD[step1_vct.at(iStep1)][step2_vct.at(iStep2)][iBar] = new TH1F (Form("hCTR_UD_ch%.3d_step1vct_%.1f_step2vct_%.1f_step1_%i_step2_%i", iBar, step1_vct.at(iStep1), step2_vct.at(iStep2), iStep1, iStep2), Form("Bar Time Res, UD_ch%.3d, step1_%.1f, step2_%.1f", iBar, step1_vct.at(iStep1), step2_vct.at(iStep2)), 4000, -1000, 1000 );
	      }                       
        }
    }
    
    //define more histos...               
    //TH2F * hTotTimeWalk  = new TH2F ("hAmpTimeWalk", "hAmpTimeWalk", 1000, -5, 5, 4000, -6000, 6000);
    //TProfile * pTotTimeWalk  = new TProfile ("pAmpTimeWalk", "pAmpTimeWalk", 1000, -5, 5);
    
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
	if (step1!=6) continue;
	//hCTR[step1][step2]->Fill(time2-time1);
	long long time_ref = chTime[384];
	hPosX->Fill(x_dut);
	// channel loop to calculate total energy deposit
	double TotalEnergy = 0;
	double TotalEnergy_corr = 0;
	for (int iCh = 0; iCh<NCH; iCh++)
	  {
	    // if not one of channels that we care about, skip it
	    if (std::find(std::begin(myChList), std::end(myChList), iCh) == std::end(myChList) ) continue;
	    // want to find the total energy in the event, in any channel
	    //std::cout << "energy " << chtot[iCh] << " for channel " << iCh << std::endl;
	    if ( chtot[iCh] < -100 ) continue; // skip any thing where the energy isnt a normal value (should always be postive, default value is -9999)

	    double ICcoeff = MIP[iCh] / avgMIP;
	    TotalEnergy += (chtot[iCh]/1.e3) / ICcoeff;

	    // corrected energy given the correction to the tot linearization
	    //double corrected_tot = (0.3639 * chtot[iCh]/1.e3)  -  (0.001911 * pow((chtot[iCh]/1.e3),2))  +  (0.00003503 * pow((chtot[iCh]/1.e3),3)) ;
	    double corrected_tot = 14.13 * (exp(0.01562 * chtot[iCh]/1.e3)-1);
	    double ICcoeff_corr = MIP_corr[iCh] / avgMIP_corr;
	    TotalEnergy_corr += corrected_tot / ICcoeff_corr;

	  }
	//std::cout << "total energy: " << TotalEnergy << std::endl;
	// channel loop in the event loop
	double corrected_tot[NCH];
	for (int iCh = 0; iCh<NCH; iCh++)
	  {
	    if (std::find(std::begin(myChList), std::end(myChList), iCh) == std::end(myChList) ) continue;
	    //if (step1 != 6) continue;
	    //if (step2 != 0) continue;

	    //if (chtot[iCh]>0) std::cout << "filling ch[" << iCh << "] with tot = " << chtot[iCh]/1.e3 << " ns :: and time-t_ref = " << chTime[iCh] << "  - " << time_ref << " = " << chTime[iCh]-time_ref << std::endl;            
	    //std::cout << Form("filling histograms for Ch_%i, step1_%f, step2_%f",iCh, step1, step2) << std::endl;
	    //std::cout << Form("channel time over threshold_%f",chtot[iCh]) << std::endl;

	    //corrected_tot[iCh] = (0.3639 * chtot[iCh]/1.e3)  -  (0.001911 * pow((chtot[iCh]/1.e3),2))  +  (0.00003503 * pow((chtot[iCh]/1.e3),3)) ;
	    corrected_tot[iCh] = 14.13 * (exp(0.01562 * chtot[iCh]/1.e3)-1);

	    hTot[step1][step2][iCh]->Fill(chtot[iCh]/1.e3);
	    hTot_correction[step1][step2][iCh]->Fill(corrected_tot[iCh]);
	    pXpos_Ypos_Tot[step1][step2][iCh]->Fill(x_dut, y_dut);
	    hTime[step1][step2][iCh]->Fill(chTime[iCh] - time_ref);

	    // calculate the cross talk, and plot this
	    // plot Tot, normalized by MIP peak energy, and then as a fraction of the total energy deposited in all channels
	    // if no energy recorded, TotalEnergy = 0, ignore this case for the cross talk calculation

	    if (TotalEnergy > 0.00001 && chtot[iCh]/1.e3 > 5 && chtot[iCh]/1.e3< 400) // tot>5 for a zero supression from low energy deposits, no cut around MIP peak, this is blue on the cross talk plots
	      {
		double ICcoeff = MIP[iCh] /avgMIP;
		pCrossTalk[step1][step2][iCh]->Fill(((chtot[iCh]/1.e3) / ICcoeff )  / TotalEnergy );

		double ICcoeff_corr = MIP_corr[iCh] /avgMIP_corr;
                pCrossTalk_corr[step1][step2][iCh]->Fill((corrected_tot[iCh] / ICcoeff_corr )  / TotalEnergy_corr );
	      }

	    //      std::cout << "htot" << std::endl;
	    //      std::cout << "ypos"<< std::endl;
	    //	    std::cout << "time"<< std::endl;

	    // try and cut on a specific bar to see how this affects MIP peak (expect to pick out Landau peak for one bar) 
	    // do this for each channel, based off of the stats found from the fit to the efficiency plots
	    // for the cross talk (green plot) cut around the MIP position requring signal in the central bar to be 0.85*MIP - 4*MIP
	    // MIP and x_dut cut should be redundant actually
	    if (x_dut < center[iCh]+1 && x_dut > center[iCh]-1 )
	      {
		hTot_cut[step1][step2][iCh]->Fill(chtot[iCh]/1.e3);
		hTot_cut_correction[step1][step2][iCh]->Fill(corrected_tot[iCh]);
		if (TotalEnergy > 0.00001 && chtot[iCh]/1.e3 > 0.85*MIP[iCh] && chtot[iCh]/1.e3<4*MIP[iCh]) 
		  {
		    double ICcoeff = MIP[iCh] /avgMIP;
		    pCrossTalkBar[step1][step2][iCh]->Fill(((chtot[iCh]/1.e3) / ICcoeff )  / TotalEnergy );

		    double ICcoeff_corr = MIP_corr[iCh] /avgMIP_corr;
		    pCrossTalkBar_corr[step1][step2][iCh]->Fill((corrected_tot[iCh] / ICcoeff_corr )  / TotalEnergy_corr );
		  }
	      }

	    // fill the 1 bar away cross talk for ch 138
	    if ( (x_dut < center[140]+1 && x_dut > center[140]-1) || (x_dut < center[136]+1 && x_dut > center[136]-1) )
	      {
		// put a MIP cut on the bar +-1 away from 138
		if (TotalEnergy > 0.00001 && chtot[iCh]/1.e3 > 0.85*MIP[iCh] && chtot[iCh]/1.e3<4*MIP[iCh] && ( iCh == 140 || iCh == 136))
		  {
		    double ICcoeff = MIP[138] /avgMIP;
		    pCrossTalkBar1away[step1][step2][138]->Fill(((chtot[138]/1.e3) / ICcoeff ) / TotalEnergy );

		    double ICcoeff_corr = MIP_corr[138] /avgMIP_corr;
                    pCrossTalkBar1away_corr[step1][step2][138]->Fill((corrected_tot[138] / ICcoeff_corr )  / TotalEnergy_corr );
		  }
	      }
	    // fill the 2 bars away cross talk for ch 138
	    if ( (x_dut < center[142]+1 && x_dut > center[142]-1 ) || (x_dut < center[134]+1 && x_dut > center[134]-1) )
	      {
		if (TotalEnergy > 0.00001 && chtot[iCh]/1.e3 > 0.85*MIP[iCh] && chtot[iCh]/1.e3<4*MIP[iCh] && (iCh == 142 || iCh == 134) )
		  {
		    double ICcoeff = MIP[138] /avgMIP;
		    pCrossTalkBar2away[step1][step2][138]->Fill(((chtot[138]/1.e3) / ICcoeff ) / TotalEnergy );

		    double ICcoeff_corr = MIP_corr[138] /avgMIP_corr;
                    pCrossTalkBar2away_corr[step1][step2][138]->Fill((corrected_tot[138] / ICcoeff_corr )  / TotalEnergy_corr );
		  }
	      }

            // fill the 1 bar away cross talk for ch 136
            if ( (x_dut < center[138]+1 && x_dut > center[138]-1) || (x_dut < center[134]+1 && x_dut > center[134]-1) )
              {
                if (TotalEnergy > 0.00001 && chtot[iCh]/1.e3 > 0.85*MIP[iCh] && chtot[iCh]/1.e3<4*MIP[iCh] && (iCh == 138 || iCh == 134))
                  {
                    double ICcoeff = MIP[136] /avgMIP;
                    pCrossTalkBar1away[step1][step2][136]->Fill(((chtot[136]/1.e3) / ICcoeff ) / TotalEnergy );

		    double ICcoeff_corr = MIP_corr[136] /avgMIP_corr;
                    pCrossTalkBar1away_corr[step1][step2][136]->Fill((corrected_tot[136] / ICcoeff_corr )  / TotalEnergy_corr );
                  }
              }
            // fill the 2 bars away cross talk for ch 136                                                  
            if ( (x_dut < center[140]+1 && x_dut > center[140]-1 ) || (x_dut < center[132]+1 && x_dut > center[132]-1) )
              {
                if (TotalEnergy > 0.00001 && chtot[iCh]/1.e3 > 0.85*MIP[iCh] && chtot[iCh]/1.e3<4*MIP[iCh] && (iCh == 140 || iCh == 132))
                  {
                    double ICcoeff = MIP[136] /avgMIP;
                    pCrossTalkBar2away[step1][step2][136]->Fill(((chtot[136]/1.e3) / ICcoeff ) / TotalEnergy );

		    double ICcoeff_corr = MIP_corr[136] /avgMIP_corr;
                    pCrossTalkBar2away_corr[step1][step2][136]->Fill((corrected_tot[136] / ICcoeff_corr )  / TotalEnergy_corr );
                  }
              }

	    //	    if (chtot[iCh] > 0. ) // needs more troubleshooting for this, why does it make the distributions so broad
	    pTot_vs_Xpos[step1][step2][iCh]->Fill(x_dut, chtot[iCh]/1.e3);
	    //pTot_vs_Ypos[step1][step2][iCh]->Fill(y_dut, chtot[iCh]/1.e3);
		
	    // make efficiency plots - 1 if 0.9 - 3 * MIP energy, 0 otherwise
	    if (x_dut != 0 && y_dut != 0)
	      {
		pXY_Edep[step1][step2][iCh]->Fill(x_dut, y_dut, chtot[iCh]/1.e3);
		if ((chtot[iCh]/1.e3) >= 0.9 * MIP[iCh] && (chtot[iCh]/1.e3) <= 3 * MIP[iCh] && x_dut > 0)
		  {
		    //std::cout<< "energy = " << chtot[iCh]/1.e3 << std::endl;
		    //std::cout<< "MIP peak value = " << MIP[iCh] << std::endl;
		    pEff_vs_Xpos[step1][step2][iCh]->Fill(x_dut, 1);
		  }
	      }
	    if (step1 > 1.e-2)
	      {
		// fill the overlay plot with the x position from each channel   
	        pTot_vs_Xpos_overlay[step1][step2]->SetMaximum(60);
	        pTot_vs_Xpos_overlay[step1][step2]->SetMinimum(-5);
		pTot_vs_Xpos_overlay[step1][step2]->Fill(x_dut, chtot[iCh]/1.e3);
	      }
	  }
        
        for (int iBar = 0; iBar<NBARS; iBar++)
        {
            int chBarUp = iBar*2+128;
            int chBarDown = iBar*2+129;
            
            if (chtot[chBarUp]>130 && chtot[chBarDown]>130 ) 
            {
                hCTR_UD[step1][step2][iBar]->Fill(chTime[chBarUp] - chTime[chBarDown]);                        
//                 std::cout << "CT_UD [" << iBar << "] = " << chTime[chBarUp] - chTime[chBarDown] << std::endl;
            }
        }
    }
    
    TCanvas *cTots_scan[NSTEP1][NSTEP2][NCH];
    TCanvas *cTots_scan_correction[NSTEP1][NSTEP2][NCH];
    TCanvas *cXpos_scan[NSTEP1][NSTEP2][NCH];
    TCanvas *cEff_scan[NSTEP1][NSTEP2][NCH];
    TCanvas *cCrossTalk[NSTEP1][NSTEP2][NCH];
    TCanvas *cCrossTalkBar[NSTEP1][NSTEP2][NCH];
    TCanvas *cCrossTalkOverlay[NSTEP1][NSTEP2][NCH];
    TCanvas *cCrossTalk138[NSTEP1][NSTEP2][NCH];
    TCanvas *cCrossTalk136[NSTEP1][NSTEP2][NCH];
    TCanvas *cCrossTalk_corr[NSTEP1][NSTEP2][NCH];
    TCanvas *cCrossTalkBar_corr[NSTEP1][NSTEP2][NCH];
    TCanvas *cCrossTalkOverlay_corr[NSTEP1][NSTEP2][NCH];
    TCanvas *cCrossTalk138_corr[NSTEP1][NSTEP2][NCH];
    TCanvas *cCrossTalk136_corr[NSTEP1][NSTEP2][NCH];
    TCanvas *cXpos_over_scan;
    //    TCanvas *cXposYpos_scan[NSTEP1][NSTEP2][NCH];
    //    TCanvas *cXpos_over_scan[NSTEP1][NSTEP2];
    //    TCanvas *cYpos_scan[NSTEP1][NSTEP2][NCH];
    //float minTotForPeakSearch = 40;

    cXpos_over_scan = new TCanvas (Form("cXpos_over"), Form("cXpos_over"), 800, 400);

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
		//if (step2_vct.at(iStep2) !=0) continue;
		std::cout<< "step1, step2 " << step1_vct.at(iStep1) << step2_vct.at(iStep2) << std::endl;

                cTots_scan[iStep1][iStep2][iCh] = new TCanvas (Form("cTots_ch%.3d_step1_%.1f_step2_%.1f", iCh, step1_vct.at(iStep1), step2_vct.at(iStep2)), Form("cTots_ch%.3d_step1_%.1f_step2_%.1f", iCh, step1_vct.at(iStep1), step2_vct.at(iStep2)), 800, 400);
                cTots_scan[iStep1][iStep2][iCh]->cd();            
                hTot[step1_vct.at(iStep1)][step2_vct.at(iStep2)][iCh]->Rebin(REBIN_COEFF/8);
                hTot[step1_vct.at(iStep1)][step2_vct.at(iStep2)][iCh]->Draw();
                hTot[step1_vct.at(iStep1)][step2_vct.at(iStep2)][iCh]->GetXaxis()->SetTitle("tot [ns]");
                hTot[step1_vct.at(iStep1)][step2_vct.at(iStep2)][iCh]->GetYaxis()->SetTitle("Counts");
                gPad->SetLogy();
		// hTot_cut is for selecting the position of one bar and plotting the Landau peak
		// use the same binning for the plot before and after the cut on a single bar
		// different line color for the MIP peak after the bar cut (this is landau peak after choosing each bar)
		hTot_cut[step1_vct.at(iStep1)][step2_vct.at(iStep2)][iCh]->Rebin(REBIN_COEFF/8);
		hTot_cut[step1_vct.at(iStep1)][step2_vct.at(iStep2)][iCh]->SetLineColor(kGreen+2);
		hTot_cut[step1_vct.at(iStep1)][step2_vct.at(iStep2)][iCh]->Draw("same");
		// do the Landau fit, and set the range over which the fit is done (reduced from 120-400 to 130-250
		TF1 * fitLandau = new TF1 ("fitLandau","landau",130,250);
		hTot_cut[step1_vct.at(iStep1)][step2_vct.at(iStep2)][iCh]->Fit(fitLandau, "QRL");
		std::cout << "Ch # " << iCh << " Landau fit normalization coeff: " << fitLandau->GetParameter(0) << " most  probable value: " << fitLandau->GetParameter(1) << " Lambda value: " << fitLandau->GetParameter(2) << std::endl;
		cTots_scan[iStep1][iStep2][iCh]->SaveAs(Form("hTot_ch%.3d_step1_%.1f_step2_%.1f.pdf", iCh, step1_vct.at(iStep1), step2_vct.at(iStep2)));


		cTots_scan_correction[iStep1][iStep2][iCh] = new TCanvas (Form("cTots_corr_ch%.3d_step1_%.1f_step2_%.1f", iCh, step1_vct.at(iStep1), step2_vct.at(iStep2)), Form("cTots_corr_ch%.3d_step1_%.1f_step2_%.1f", iCh, step1_vct.at(iStep1), step2_vct.at(iStep2)), 800, 400);
                cTots_scan_correction[iStep1][iStep2][iCh]->cd();
                hTot_correction[step1_vct.at(iStep1)][step2_vct.at(iStep2)][iCh]->Rebin(REBIN_COEFF/8);
                hTot_correction[step1_vct.at(iStep1)][step2_vct.at(iStep2)][iCh]->Draw();
                hTot_correction[step1_vct.at(iStep1)][step2_vct.at(iStep2)][iCh]->GetXaxis()->SetTitle("Corrected tot [ns]");
                hTot_correction[step1_vct.at(iStep1)][step2_vct.at(iStep2)][iCh]->GetYaxis()->SetTitle("Counts");
                gPad->SetLogy();
                // hTot_cut is for selecting the position of one bar and plotting the Landau peak
                // use the same binning for the plot before and after the cut on a single bar                                                         
                // different line color for the MIP peak after the bar cut (this is landau peak after choosing each bar)                                                   
                hTot_cut_correction[step1_vct.at(iStep1)][step2_vct.at(iStep2)][iCh]->Rebin(REBIN_COEFF/8);
                hTot_cut_correction[step1_vct.at(iStep1)][step2_vct.at(iStep2)][iCh]->SetLineColor(kGreen+2);
                hTot_cut_correction[step1_vct.at(iStep1)][step2_vct.at(iStep2)][iCh]->Draw("same");
                // do the Landau fit, and set the range over which the fit is done (reduced from 120-400 to 130-250                              
		TF1 * fitLandau_corr = new TF1 ("fitLandau_corr","landau",92,519);
                hTot_cut_correction[step1_vct.at(iStep1)][step2_vct.at(iStep2)][iCh]->Fit(fitLandau_corr, "QRL");
		std::cout << "Ch # " << iCh << " Landau fit on corrected tot normalization coeff: " << fitLandau_corr->GetParameter(0) << " most  probable value: " << fitLandau_corr->GetParameter(1) << " Lambda value: " << fitLandau_corr->GetParameter(2) << std::endl;
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

		if (iCh == 138 )
		  {
		    cCrossTalk138[iStep1][iStep2][138] = new  TCanvas (Form("cCrossTalk138_step1_%.1f_step2_%.1f", step1_vct.at(iStep1), step2_vct.at(iStep2)), Form("cCrossTalk138_step1_%.1f_step2_%.1f", step1_vct.at(iStep1), step2_vct.at(iStep2)), 800, 400);
		    cCrossTalk138[iStep1][iStep2][138]->cd();
		    pCrossTalk[step1_vct.at(iStep1)][step2_vct.at(iStep2)][138]->Draw();
		    pCrossTalkBar[step1_vct.at(iStep1)][step2_vct.at(iStep2)][138]->SetLineColor(kGreen+2);
		    pCrossTalkBar[step1_vct.at(iStep1)][step2_vct.at(iStep2)][138]->Draw("same");
		    pCrossTalkBar1away[step1_vct.at(iStep1)][step2_vct.at(iStep2)][138]->SetLineColor(kRed+2);
		    pCrossTalkBar1away[step1_vct.at(iStep1)][step2_vct.at(iStep2)][138]->Draw("same");
		    pCrossTalkBar2away[step1_vct.at(iStep1)][step2_vct.at(iStep2)][138]->SetLineColor(kMagenta+2);
		    pCrossTalkBar2away[step1_vct.at(iStep1)][step2_vct.at(iStep2)][138]->Draw("same");
		    pCrossTalk[step1_vct.at(iStep1)][step2_vct.at(iStep2)][138]->GetXaxis()->SetTitle("Fractional energy deposit");
		    pCrossTalk[step1_vct.at(iStep1)][step2_vct.at(iStep2)][138]->GetYaxis()->SetTitle("Events");
		    gPad->SetLogy();
		    cCrossTalk138[iStep1][iStep2][138]->SaveAs(Form("crosstalk138_step1_%.1f_step2_%.1f.pdf", step1_vct.at(iStep1), step2_vct.at(iStep2)));

		    cCrossTalk138_corr[iStep1][iStep2][138] = new  TCanvas (Form("cCrossTalk138_corr_step1_%.1f_step2_%.1f", step1_vct.at(iStep1), step2_vct.at(iStep2)), Form("cCrossTalk138_corr_step1_%.1f_step2_%.1f", step1_vct.at(iStep1), step2_vct.at(iStep2)), 800, 400);
                    cCrossTalk138_corr[iStep1][iStep2][138]->cd();
                    pCrossTalk_corr[step1_vct.at(iStep1)][step2_vct.at(iStep2)][138]->Draw();
                    pCrossTalkBar_corr[step1_vct.at(iStep1)][step2_vct.at(iStep2)][138]->SetLineColor(kGreen+2);
                    pCrossTalkBar_corr[step1_vct.at(iStep1)][step2_vct.at(iStep2)][138]->Draw("same");
                    pCrossTalkBar1away_corr[step1_vct.at(iStep1)][step2_vct.at(iStep2)][138]->SetLineColor(kRed+2);
                    pCrossTalkBar1away_corr[step1_vct.at(iStep1)][step2_vct.at(iStep2)][138]->Draw("same");
                    pCrossTalkBar2away_corr[step1_vct.at(iStep1)][step2_vct.at(iStep2)][138]->SetLineColor(kMagenta+2);
                    pCrossTalkBar2away_corr[step1_vct.at(iStep1)][step2_vct.at(iStep2)][138]->Draw("same");
                    pCrossTalk_corr[step1_vct.at(iStep1)][step2_vct.at(iStep2)][138]->GetXaxis()->SetTitle("Fractional energy deposit (corrected)");
                    pCrossTalk_corr[step1_vct.at(iStep1)][step2_vct.at(iStep2)][138]->GetYaxis()->SetTitle("Events");
                    gPad->SetLogy();
                    cCrossTalk138_corr[iStep1][iStep2][138]->SaveAs(Form("crosstalk138_corr_step1_%.1f_step2_%.1f.pdf", step1_vct.at(iStep1), step2_vct.at(iStep2)));
		  }

                if (iCh == 136 )
                  {
                    cCrossTalk136[iStep1][iStep2][136] = new  TCanvas (Form("cCrossTalk136_step1_%.1f_step2_%.1f", step1_vct.at(iStep1), step2_vct.at(iStep2)), Form("cCrossTalk136_step1_%.1f_step2_%.1f", step1_vct.at(iStep1), step2_vct.at(iStep2)),800, 400);
                    cCrossTalk136[iStep1][iStep2][136]->cd();
                    pCrossTalk[step1_vct.at(iStep1)][step2_vct.at(iStep2)][136]->Draw();
                    pCrossTalkBar[step1_vct.at(iStep1)][step2_vct.at(iStep2)][136]->SetLineColor(kGreen+2);
                    pCrossTalkBar[step1_vct.at(iStep1)][step2_vct.at(iStep2)][136]->Draw("same");
                    pCrossTalkBar1away[step1_vct.at(iStep1)][step2_vct.at(iStep2)][136]->SetLineColor(kRed+2);
                    pCrossTalkBar1away[step1_vct.at(iStep1)][step2_vct.at(iStep2)][136]->Draw("same");
                    pCrossTalkBar2away[step1_vct.at(iStep1)][step2_vct.at(iStep2)][136]->SetLineColor(kMagenta+2);
                    pCrossTalkBar2away[step1_vct.at(iStep1)][step2_vct.at(iStep2)][136]->Draw("same");
                    pCrossTalk[step1_vct.at(iStep1)][step2_vct.at(iStep2)][136]->GetXaxis()->SetTitle("Fractional energy deposit");
                    pCrossTalk[step1_vct.at(iStep1)][step2_vct.at(iStep2)][136]->GetYaxis()->SetTitle("Events");
                    gPad->SetLogy();
                    cCrossTalk136[iStep1][iStep2][136]->SaveAs(Form("crosstalk136_step1_%.1f_step2_%.1f.pdf", step1_vct.at(iStep1), step2_vct.at(iStep2)));

		    cCrossTalk136_corr[iStep1][iStep2][136] = new  TCanvas (Form("cCrossTalk136_corr_step1_%.1f_step2_%.1f", step1_vct.at(iStep1), step2_vct.at(iStep2)), Form("cCrossTalk136_corr_step1_%.1f_step2_%.1f", step1_vct.at(iStep1), step2_vct.at(iStep2)),800, 400);
                    cCrossTalk136_corr[iStep1][iStep2][136]->cd();
                    pCrossTalk_corr[step1_vct.at(iStep1)][step2_vct.at(iStep2)][136]->Draw();
                    pCrossTalkBar_corr[step1_vct.at(iStep1)][step2_vct.at(iStep2)][136]->SetLineColor(kGreen+2);
                    pCrossTalkBar_corr[step1_vct.at(iStep1)][step2_vct.at(iStep2)][136]->Draw("same");
                    pCrossTalkBar1away_corr[step1_vct.at(iStep1)][step2_vct.at(iStep2)][136]->SetLineColor(kRed+2);
                    pCrossTalkBar1away_corr[step1_vct.at(iStep1)][step2_vct.at(iStep2)][136]->Draw("same");
                    pCrossTalkBar2away_corr[step1_vct.at(iStep1)][step2_vct.at(iStep2)][136]->SetLineColor(kMagenta+2);
                    pCrossTalkBar2away_corr[step1_vct.at(iStep1)][step2_vct.at(iStep2)][136]->Draw("same");
                    pCrossTalk_corr[step1_vct.at(iStep1)][step2_vct.at(iStep2)][136]->GetXaxis()->SetTitle("Fractional energy deposit (corrected)");
                    pCrossTalk_corr[step1_vct.at(iStep1)][step2_vct.at(iStep2)][136]->GetYaxis()->SetTitle("Events");
                    gPad->SetLogy();
                    cCrossTalk136_corr[iStep1][iStep2][136]->SaveAs(Form("crosstalk136_corr_step1_%.1f_step2_%.1f.pdf", step1_vct.at(iStep1), step2_vct.at(iStep2)));
                  }

		/*// plotting beam profile
		cXposYpos_scan[iStep1][iStep2][iCh] = new TCanvas (Form("cXpos_ch%.3d_step1_%.1f_step2_.%1f", iCh, step1_vct.at(iStep1), step1_vct.at(iStep1)), Form("cXpos_ch%.3d_step1_%.1f_step2_.%1f", iCh, step1_vct.at(iStep1), step1_vct.at(iStep1)), 800, 400);
                cXposYpos_scan[iStep1][iStep2][iCh]->cd();
		pXpos_Ypos_Tot[step1_vct.at(iStep1)][step2_vct.at(iStep2)][iCh]->Draw("colz");
		pXpos_Ypos_Tot[step1_vct.at(iStep1)][step2_vct.at(iStep2)][iCh]->GetXaxis()->SetTitle("x position");
		pXpos_Ypos_Tot[step1_vct.at(iStep1)][step2_vct.at(iStep2)][iCh]->GetYaxis()->SetTitle("y position");
		cXposYpos_scan[iStep1][iStep2][iCh]->SaveAs(Form("pXposYpos_ch%.3d.pdf", iCh));
		*/

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
                cXpos_scan[iStep1][iStep2][iCh]->SaveAs(Form("pTot_vs_Xpos_ch%.3d_step1_%.1f_step2_%.1f.pdf", iCh, step1_vct.at(iStep1), step2_vct.at(iStep2)));

		// efficiency plots
		cEff_scan[iStep1][iStep2][iCh] = new TCanvas (Form("cEff_ch%.3d_step1_%.1f_step2_%.1f", iCh, step1_vct.at(iStep1), step2_vct.at(iStep2)), Form("cEff__ch%.3d_step1_%.1f_step2_%.1f", iCh, step1_vct.at(iStep1), step2_vct.at(iStep2)), 800, 400);
		cEff_scan[iStep1][iStep2][iCh]->cd();
		pEff_vs_Xpos[step1_vct.at(iStep1)][step2_vct.at(iStep2)][iCh]->Rebin(2);
		pEff_vs_Xpos[step1_vct.at(iStep1)][step2_vct.at(iStep2)][iCh]->Draw();
		pEff_vs_Xpos[step1_vct.at(iStep1)][step2_vct.at(iStep2)][iCh]->GetXaxis()->SetTitle("X position");
		pEff_vs_Xpos[step1_vct.at(iStep1)][step2_vct.at(iStep2)][iCh]->GetYaxis()->SetTitle("Within MIP Peak Energy");
		cEff_scan[iStep1][iStep2][iCh]->SaveAs(Form("cEff_ch%.3d_step1_%.1f_step2_%.1f.pdf", iCh, step1_vct.at(iStep1), step2_vct.at(iStep2)));
		// need to find the center of the efficiency plot for one channel - use a specific fit function - this is done later in overlay plots in channel / bar loop

		// plots for overlay of x position 
		cXpos_over_scan->cd();
		pTot_vs_Xpos[step1_vct.at(iStep1)][step2_vct.at(iStep2)][iCh]->Draw("same");
		pTot_vs_Xpos[step1_vct.at(iStep1)][step2_vct.at(iStep2)][iCh]->GetYaxis()->SetRange(-5,70);
		pTot_vs_Xpos[step1_vct.at(iStep1)][step2_vct.at(iStep2)][iCh]->SetMarkerColor(6);
		pTot_vs_Xpos[step1_vct.at(iStep1)][step2_vct.at(iStep2)][iCh]->SetTitle("X position");
		pTot_vs_Xpos[step1_vct.at(iStep1)][step2_vct.at(iStep2)][iCh]->SetTitle("tot [ns]");
            }
          cXpos_over_scan->SaveAs(Form("pTot_vs_Xpos_overlay_step1_%.1f_step2_%.1f.pdf", step1_vct.at(iStep1), step2_vct.at(iStep2)));
            
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
	if ( myChList[chId] == 143 )
	  {
	    cArrayTots->SaveAs(Form("hTot_array_chId%.3d.pdf", myChList[chId]));
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
	if ( myChList[chId] == 143 )
	  {
	    cArrayTimes->SaveAs(Form("hTime_array_chId%.3d.pdf", myChList[chId]));
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

    TCanvas * cArrayEff = new TCanvas("cArrayEff","cArrayEff", 1600, 300);
    cArrayEff->Divide(NBARS,2);
    for (int chId = 0; chId<NBARS*2; chId++)
      {
	// normalize the efficiency plots
	pEff_vs_Xpos[step1_vct.at(selStep1)][step2_vct.at(selStep2)][myChList[chId]]->Rebin(2);
        for (int iBin = 0; iBin < 200; iBin++)
	  {
            if (hPosX->GetBinContent(iBin+1)>0) pEff_vs_Xpos[step1_vct.at(selStep1)][step2_vct.at(selStep2)][myChList[chId]]->SetBinContent(iBin+1, pEff_vs_Xpos[step1_vct.at(selStep1)][step2_vct.at(selStep2)][myChList[chId]]->GetBinContent(iBin+1)/hPosX->GetBinContent(iBin+1) );
	  }
	cArrayEff->cd(chId+1);
        //pEff_vs_Xpos[step1_vct.at(selStep1)][step2_vct.at(selStep2)][myChList[chId]]->Rebin(2);
        pEff_vs_Xpos[step1_vct.at(selStep1)][step2_vct.at(selStep2)][myChList[chId]]->SetAxisRange(0,1.2,"Y");
        pEff_vs_Xpos[step1_vct.at(selStep1)][step2_vct.at(selStep2)][myChList[chId]]->Draw();
        pEff_vs_Xpos[step1_vct.at(selStep1)][step2_vct.at(selStep2)][myChList[chId]]->SetTitle("X position");
        pEff_vs_Xpos[step1_vct.at(selStep1)][step2_vct.at(selStep2)][myChList[chId]]->SetTitle("Efficiency within MIP Peak Energy");
      }

    TF1 * fitBarPos = new TF1 ("fitBarPos", fitBarEffErr, 2, 32, 5);
    fitBarPos->SetParameters(5., 3., 0.01, 0.8, 0.2);
    fitBarPos->SetNpx(5000);
    fitBarPos->SetParLimits(0, minXpos, maxXpos);
    fitBarPos->SetParLimits(1, 2, 4);
    fitBarPos->SetParLimits(3, 0.2, 1.);

    TCanvas *cArrayEffOverlay = new TCanvas("cArrayEffOverlay","cArrayEffOverlay",800,400);
    cArrayEffOverlay->cd();
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
	fitBarPos->SetParameter(0, 4 + chId*3.);
	for (int i = 0; i< 3; i++) 
	  {
	    pEff_vs_Xpos[step1_vct.at(selStep1)][step2_vct.at(selStep2)][myChList[chId]]->Fit(fitBarPos,"QR");
	  }
        float posBar   = fitBarPos->GetParameter(0);
        float widthBar = fitBarPos->GetParameter(1);
        float effBar   = fitBarPos->GetParameter(3);
        float trackRes = fitBarPos->GetParameter(4);
	std::cout << "posBar[" << chId << "] = " << posBar << " :: width = " << widthBar << " :: effBar = " << effBar << " :: trackRes = " << trackRes << std::endl;
	if (chId == NBARS - 1)
	  {
	    cArrayEff->SaveAs(Form("Efficiency_array_barId%.3d.pdf", chId));
	    cArrayEffOverlay->SaveAs(Form("Efficiency_array_overlay_barId%.3d.pdf", chId));
	  }
      }

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
    cSingleScatterEdep->SaveAs(Form("XY_xcatter.pdf"));

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

}
