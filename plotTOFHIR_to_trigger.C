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

int main(int argc, char** argv)
{
    TApplication* theApp = new TApplication("App", &argc, argv);
    TLegend *leg;

    //read input files
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
      data_path = std::string(argv[3]);
    }
    std::cout << "data_path = " << data_path << std::endl;
    
    
    //define tree
    float	step1;
    float	step2;
    float	x_dut;
    float	y_dut;
    int 	ntracks;        
    long long	 chTime[400];    
    float	 chtot[400];
            
    TChain * tree = new TChain("data", "data");
    
    for (int iRun = firstRun; iRun <= lastRun; iRun++)
    {
//        if ( std::find(bad_runs.begin(), bad_runs.end(), iRun) != bad_runs.end()) continue;       
       tree->Add( Form("%s/run%.5d_events_withTrack.root", data_path.c_str(), iRun) );
       //       tree->Add( Form("%s/run%.5d_events.root", data_path.c_str(), iRun) ); 
       std::cout << "adding run: " << iRun << std::endl;
    }
    
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
        
    Long64_t NEVENTS = tree->GetEntries();
    std::cout << "nEvents = " << NEVENTS << std::endl;
    
    for (Int_t iEvt= 0; iEvt < NEVENTS; iEvt++) 
    {
        tree->GetEntry(iEvt);
//         std::cout << "step1 = " << step1 << " :: step2 = " << step2 <<std::endl;
	if (ntracks != 1) continue;
        if (std::find(step1_vct.begin(), step1_vct.end(), step1) == step1_vct.end() )
        {             
            step1_vct.push_back(step1);
            NSTEP1++;
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
    
    double minXpos = -5;
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
    // centers for first array when both sides of the bars are read out
    //    /*
    center[128] = 5.5;
    center[129] = 5.5;
    center[130] = 0;
    center[131] = 8;
    center[132] = 11;
    center[133] = 12;
    center[134] = 14;
    center[135] = 15;
    center[136] = 17;
    center[137] = 18;
    center[138] = 20;
    center[139] = 20;
    center[140] = 24;
    center[141] = 23;
    center[142] = 26;
    center[143] = 26;
    //    */
    // centers for second array when only one side of the bar is read out    
    /*
    center[128] = 5;           
    center[132] = 12;                                                                         
    center[134] = 15;                                                                                        
    center[136] = 18;                                                                                                                    
    center[138] = 21;                                                                        
    center[140] = 23;                                                                                                         
    center[142] = 27;                                                                          
    */

    // declare the histograms, these will be filled in the channel loop    
    std::map<float, std::map<float, std::map<int, TH1F * > > > hTot;
    std::map<float, std::map<float, std::map<int, TH1F * > > > hTot_cut;
    std::map<float, std::map<float, std::map<int, TH1F * > > > hTime;
    std::map<float, std::map<float, std::map<int, TProfile * > > > pTot_vs_Xpos;
    std::map<float, std::map<float, std::map<int, TProfile * > > > pTot_vs_Ypos;
    std::map<float, std::map<float, std::map<int, TH2F * > > > pXpos_Ypos_Tot;
    std::map<float, std::map<float, std::map<int, TH1F * > > > hCTR_UD;
    std::map<float, std::map<float, TProfile * > > pTot_vs_Xpos_overlay;
    
//     std::map<float, std::map<float, TProfile2D * > > hTot_XY[NCH];
//     std::map<float, std::map<float, float > > tot_mean;    
//     std::map<float, std::map<float, float > > time_mean;
        
    // step 1, step 2, and channel listing loops. Histograms are defined inside the loops
    for (int iStep1 = 0; iStep1< NSTEP1; iStep1++)
    {
        for (int iStep2 = 0; iStep2< NSTEP2; iStep2++)
        {
            for (int iCh = 0; iCh < NCH; iCh++)
            {       
	      //std::cout << "in histogram loop " << std::endl;
	      //std::cout << Form("hTot_ch%.3d_step1_%.1f_step2_%.1f", iCh, step1_vct.at(iStep1), step2_vct.at(iStep2)) << std::endl;
	      //std::cout << Form("step1_%i, step2_%i",iStep1, iStep2) << std::endl;

	      hTot[step1_vct.at(iStep1)][step2_vct.at(iStep2)][iCh] = new TH1F (Form("hTot_ch%.3d_step1vct_%.1f_step2vct_%.1f_step1_%i_step2_%i", iCh, step1_vct.at(iStep1), step2_vct.at(iStep2), iStep1, iStep2), Form("Time over threshold, ch%.3d, step1_%.1f, step2_%.1f", iCh, step1_vct.at(iStep1), step2_vct.at(iStep2)), 4000, minTot, maxTot );

	      hTot_cut[step1_vct.at(iStep1)][step2_vct.at(iStep2)][iCh] = new TH1F (Form("hTot_cut_ch%.3d_step1vct_%.1f_step2vct_%.1f_step1_%i_step2_%i", iCh, step1_vct.at(iStep1), step2_vct.at(iStep2), iStep1, iStep2), Form("Time over threshold with x cut, ch%.3d, step1_%.1f, step2_%.1f", iCh, step1_vct.at(iStep1), step2_vct.at(iStep2)), 4000, minTot, maxTot );
                
	      hTime[step1_vct.at(iStep1)][step2_vct.at(iStep2)][iCh] = new TH1F (Form("hTime_ch%.3d_step1vct_%.1f_step2vct_%.1f_step1_%i_step2_%i", iCh, step1_vct.at(iStep1), step2_vct.at(iStep2), iStep1, iStep2), Form("Time Hist, ch%.3d, step1_%.1f, step2_%.1f", iCh, step1_vct.at(iStep1), step2_vct.at(iStep2)), 4000, minTime, maxTime );

	      pTot_vs_Xpos[step1_vct.at(iStep1)][step2_vct.at(iStep2)][iCh] = new TProfile (Form("pTot_vs_Xpos_ch%.3d_step1vct_%.1f_step2vct_%.1f_step1_%i_step2_%i", iCh, step1_vct.at(iStep1), step2_vct.at(iStep2), iStep1, iStep2), Form("ToT vs X pos, ch%.3d, step1_%.1f, step2_%.1f", iCh, step1_vct.at(iStep1), step2_vct.at(iStep2)), 4000, minXpos, maxXpos );

	      pTot_vs_Ypos[step1_vct.at(iStep1)][step2_vct.at(iStep2)][iCh] = new TProfile (Form("pTot_vs_Ypos_ch%.3d_step1vct_%.1f_step2vct_%.1f_step1_%i_step2_%i", iCh, step1_vct.at(iStep1), step2_vct.at(iStep2), iStep1, iStep2), Form("ToT vs Y pos, ch%.3d, step1_%.1f, step2_%.1f", iCh, step1_vct.at(iStep1), step2_vct.at(iStep2)), 4000, minYpos, maxYpos );

	      pXpos_Ypos_Tot[step1_vct.at(iStep1)][step2_vct.at(iStep2)][iCh] = new TH2F (Form("pXpos_Ypos_Tot_ch%.3d_step1vct_%.1f_step2vct_%.1f_step1_%i_step2_%i", iCh, step1_vct.at(iStep1), step2_vct.at(iStep2), iStep1, iStep2), Form("X and Y pos vs. ToT, ch%.3d, step1_%.1f, step2_%.1f", iCh, step1_vct.at(iStep1), step2_vct.at(iStep2)), 40, minXpos, maxXpos, 40, minYpos, maxYpos );

            }

	    pTot_vs_Xpos_overlay[step1_vct.at(iStep1)][step2_vct.at(iStep2)] = new TProfile (Form("pTot_vs_Xpos_over_step1vct_%.1f_step2vct_%.1f_step1_%i_step2_%i", step1_vct.at(iStep1), step2_vct.at(iStep2), iStep1, iStep2), Form("ToT vs X pos overlay, step1_%.1f, step2_%.1f", step1_vct.at(iStep1), step2_vct.at(iStep2)), 4000, minXpos, maxXpos );

	    // bar loop (outside of channel loop) to define time resolution for a single bar
            for (int iBar = 0; iBar < NBARS; iBar++)
            {
	      hCTR_UD[step1_vct.at(iStep1)][step2_vct.at(iStep2)][iBar] = new TH1F (Form("hCTR_UD_ch%.3d_step1vct_%.1f_step2vct_%.1f_step1_%i_step2_%i", iBar, step1_vct.at(iStep1), step2_vct.at(iStep2), iStep1, iStep2), Form("Bar Time Res, UD_ch%.3d, step1_%.1f, step2_%.1f", iBar, step1_vct.at(iStep1), step2_vct.at(iStep2)), 4000, -1000, 1000 );
            }                       
        }
    }

    
    //define more histos...               

//     TH2F * hTotTimeWalk  = new TH2F ("hAmpTimeWalk", "hAmpTimeWalk", 1000, -5, 5, 4000, -6000, 6000);
//     TProfile * pTotTimeWalk  = new TProfile ("pAmpTimeWalk", "pAmpTimeWalk", 1000, -5, 5);
    
    

    //************************************************************************************//
    //              loop 0
    //************************************************************************************//
    
    std::cout << "(0) looping over events to get MIP peak position" << std::endl;
    for (Int_t iEvt= 0; iEvt < NEVENTS; iEvt++) 
    {
        tree->GetEntry(iEvt);
	if (ntracks != 1) continue; // require 1 MIP per event

        if (iEvt%1000 == 0) std::cout << "processing event: " << iEvt << "\r" << std::flush;
//         if (step1!=0) continue;
//         hCTR[step1][step2]->Fill(time2-time1);
        long long time_ref = chTime[384];
	
	// channel loop in the event loop
        for (int iCh = 0; iCh<NCH; iCh++)
        {
            if (std::find(std::begin(myChList), std::end(myChList), iCh) == std::end(myChList) ) continue;
            
//             if (chtot[iCh]>0) std::cout << "filling ch[" << iCh << "] with tot = " << chtot[iCh]/1.e3 << " ns :: and time-t_ref = " << chTime[iCh] << "  - " << time_ref << " = " << chTime[iCh]-time_ref << std::endl;
            
	    //std::cout << Form("filling histograms for Ch_%i, step1_%f, step2_%f",iCh, step1, step2) << std::endl;
	    //	    std::cout << Form("channel time over threshold_%f",chtot[iCh]) << std::endl;

	    hTot[step1][step2][iCh]->Fill(chtot[iCh]/1.e3);
	    //	    std::cout << "htot" << std::endl;
	    pXpos_Ypos_Tot[step1][step2][iCh]->Fill(x_dut, y_dut);
	    //	    std::cout << "ypos"<< std::endl;
	    hTime[step1][step2][iCh]->Fill(chTime[iCh] - time_ref);
	    //	    std::cout << "time"<< std::endl;

	    // try and cut on a specific bar to see how this affects MIP peak (expect to pick out Landau peak for one bar) 
	    if (x_dut < 4.59+1 && x_dut > 4.59-1 )
	      {
		hTot_cut[step1][step2][iCh]->Fill(chtot[iCh]/1.e3);
	      }

	    //	    if (chtot[iCh] > 0. ) // needs more troubleshooting for this, why does it make the distributions so broad
	      {
		pTot_vs_Xpos[step1][step2][iCh]->Fill(x_dut, chtot[iCh]/1.e3);
		pTot_vs_Ypos[step1][step2][iCh]->Fill(y_dut, chtot[iCh]/1.e3);
		// fill the overlay plot with the x position from each channel   
	        // pTot_vs_Xpos_overlay[step1][step2]->GetYaxis()->SetRange(-5,60);
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
    
    
    float minTotForPeakSearch = 40;
    
    TCanvas *cTots_scan[NSTEP1][NSTEP2][NCH];
    TCanvas *cXpos_scan[NSTEP1][NSTEP2][NCH];
    TCanvas *cYpos_scan[NSTEP1][NSTEP2][NCH];
    TCanvas *cXposYpos_scan[NSTEP1][NSTEP2][NCH];
    TCanvas *cXpos_over_scan;
    //    TCanvas *cXpos_over_scan[NSTEP1][NSTEP2];

    cXpos_over_scan = new TCanvas (Form("cXpos_over"), Form("cXpos_over"), 800, 400);

    // now draw all the histograms, again have step 1, step 2, channel loop
    for (int iStep1 = 0; iStep1< NSTEP1; iStep1++)
    {
        for (int iStep2 = 0; iStep2< NSTEP2; iStep2++)
        {
            for (int iCh = 0; iCh< NCH; iCh++)
            {
                if (std::find(std::begin(myChList), std::end(myChList), iCh) == std::end(myChList) ) continue;
                
                cTots_scan[iStep1][iStep2][iCh] = new TCanvas (Form("cTots_ch%.3d_step1_%.1f_step2_%.1f", iCh, step1_vct.at(iStep1), step2_vct.at(iStep2)), Form("cTots_ch%.3d_step1_%.1f_step2_%.1f", iCh, step1_vct.at(iStep1), step2_vct.at(iStep2)), 800, 400);
//                 cTots_scan[iStep1][iStep2][iCh]->Divide(2,1);
                
                cTots_scan[iStep1][iStep2][iCh]->cd();            
                hTot[step1_vct.at(iStep1)][step2_vct.at(iStep2)][iCh]->Rebin(REBIN_COEFF);
                hTot[step1_vct.at(iStep1)][step2_vct.at(iStep2)][iCh]->Draw();
                hTot[step1_vct.at(iStep1)][step2_vct.at(iStep2)][iCh]->GetXaxis()->SetTitle("tot [ns]");
                hTot[step1_vct.at(iStep1)][step2_vct.at(iStep2)][iCh]->GetYaxis()->SetTitle("Counts");
                gPad->SetLogy();
		// hTot_cut is for selecting the position of one bar and plotting the Landau peak
		hTot_cut[step1_vct.at(iStep1)][step2_vct.at(iStep2)][iCh]->Draw("same");
		cTots_scan[iStep1][iStep2][iCh]->SaveAs(Form("hTot_ch%.3d_step1_%.1f_step2_%.1f.pdf", iCh, step1_vct.at(iStep1), step2_vct.at(iStep2)));

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

		//float mean = pTot_vs_Xpos[step1_vct.at(iStep1)][step2_vct.at(iStep2)][iCh]->GetMean();
		//float rms  = pTot_vs_Xpos[step1_vct.at(iStep1)][step2_vct.at(iStep2)][iCh]->GetRMS();
		//TF1 * fitGaus = new TF1 ("fitGaus", "gaus", mean-rms , mean+rms);
		//pTot_vs_Xpos[step1_vct.at(iStep1)][step2_vct.at(iStep2)][iCh]->Fit(fitGaus, "QRL");
		//std::cout << "Xpos ch [" << iCh << "] = " << fitGaus->GetParameter(2) << " ps --> sigma_bar = " << fitGaus->GetParameter(2)/2 << " ps " << std::endl;
		std::cout << "xPos ch[" << iCh << "] = " << fitGaus->GetParameter(1) << " x position centered " << std::endl;
                cXpos_scan[iStep1][iStep2][iCh]->SaveAs(Form("pTot_vs_Xpos_ch%.3d_step1_%.1f_step2_%.1f.pdf", iCh, step1_vct.at(iStep1), step2_vct.at(iStep2)));

                // making plots for the y position of the device under test (y_dut)
                cYpos_scan[iStep1][iStep2][iCh] = new TCanvas (Form("cYpos_ch%.3d_step1_%.1f_step2_%.1f", iCh, step1_vct.at(iStep1), step2_vct.at(iStep2)), Form("cYpos_ch%.3d_step1_%.1f_step2_%.1f", iCh, step1_vct.at(iStep1), step2_vct.at(iStep2)), 800, 400);
                cYpos_scan[iStep1][iStep2][iCh]->cd();
                pTot_vs_Ypos[step1_vct.at(iStep1)][step2_vct.at(iStep2)][iCh]->Rebin(REBIN_COEFF);
                pTot_vs_Ypos[step1_vct.at(iStep1)][step2_vct.at(iStep2)][iCh]->Draw();
                pTot_vs_Ypos[step1_vct.at(iStep1)][step2_vct.at(iStep2)][iCh]->GetXaxis()->SetTitle("Y position");
                pTot_vs_Ypos[step1_vct.at(iStep1)][step2_vct.at(iStep2)][iCh]->GetYaxis()->SetTitle("tot [ns]");
		//TF1 * fitGausY = new TF1 ("fitGausY", "gaus", -100 , 150 );
		//pTot_vs_Ypos[step1_vct.at(iStep1)][step2_vct.at(iStep2)][iCh]->Fit(fitGausY, "QRL");
		//std::cout << "yPos ch[" << iCh << "] = " << fitGausY->GetParameter(1) << " y position centered " << std::endl;
                cYpos_scan[iStep1][iStep2][iCh]->SaveAs(Form("pTot_vs_Ypos_ch%.3d_step1_%.1f_step2_%.1f.pdf", iCh, step1_vct.at(iStep1), step2_vct.at(iStep2)));

		// plots for overlay of x position 
		//cXpos_over_scan[iStep1][iStep2] = new TCanvas (Form("cXpos_over_step1_%.1f_step2_.%1f", step1_vct.at(iStep1), step2_vct.at(iStep2)), Form("cXpos_over_step1_%.1f_step2_.%1f", step1_vct.at(iStep1), step2_vct.at(iStep2)), 800, 400);
		cXpos_over_scan->cd();
		//pTot_vs_Xpos[step1_vct.at(iStep1)][step2_vct.at(iStep2)][iCh]->Rebin(REBIN_COEFF);
		pTot_vs_Xpos[step1_vct.at(iStep1)][step2_vct.at(iStep2)][iCh]->Draw("same");
		pTot_vs_Xpos[step1_vct.at(iStep1)][step2_vct.at(iStep2)][iCh]->GetYaxis()->SetRange(-5,70);
		pTot_vs_Xpos[step1_vct.at(iStep1)][step2_vct.at(iStep2)][iCh]->SetMarkerColor(6);
		pTot_vs_Xpos[step1_vct.at(iStep1)][step2_vct.at(iStep2)][iCh]->SetTitle("X position");
		pTot_vs_Xpos[step1_vct.at(iStep1)][step2_vct.at(iStep2)][iCh]->SetTitle("tot [ns]");
		//cXpos_over_scan->SaveAs(Form("pTot_vs_Xpos_overlay.pdf"));
            }
            cXpos_over_scan->SaveAs(Form("pTot_vs_Xpos_overlay_step1_%.1f_step2_%.1f.pdf", step1_vct.at(iStep1), step2_vct.at(iStep2)));
            
            /*
            hTot1[step1_vct.at(iStep1)][step2_vct.at(iStep2)]->GetXaxis()->SetRange(hTot1[step1_vct.at(iStep1)][step2_vct.at(iStep2)]->GetXaxis()->FindBin(minTotForPeakSearch) , hTot1[step1_vct.at(iStep1)][step2_vct.at(iStep2)]->GetXaxis()->GetNbins());
//             float x_max = hTot1[step1_vct.at(iStep1)][step2_vct.at(iStep2)]->GetBinCenter(hTot1[step1_vct.at(iStep1)][step2_vct.at(iStep2)]->GetMaximumBin());
            
            int npeaks = 3;
            TSpectrum *s = new TSpectrum(npeaks);
            Int_t nfound = s->Search(hTot1[step1_vct.at(iStep1)][step2_vct.at(iStep2)], 20, "nobackground", 0.1);    //sigma, option, threshold
            printf("Found %d candidate peaks in tot1 to fit\n",nfound);
            Double_t *xpeaks = s->GetPositionX();               
            Double_t xp = xpeaks[0];
            for (int ip = 0; ip<npeaks; ip++)
            {
               if (xp<xpeaks[ip]) xp = xpeaks[ip];
            }
            float x_max = xp;
        
            TF1 * fitTots = new TF1 ("fitTots", "gaus", x_max*0.9, x_max*1.15);
            hTot1[step1_vct.at(iStep1)][step2_vct.at(iStep2)]->Fit(fitTots, "QR");
            tot1_mean[step1_vct.at(iStep1)][step2_vct.at(iStep2)] = fitTots->GetParameter(1);
            tot1_sigma[step1_vct.at(iStep1)][step2_vct.at(iStep2)] = fitTots->GetParameter(2);
            
            std::cout << "tot1[" << step1_vct.at(iStep1)<<"]["<<step2_vct.at(iStep2)<<"] = " << tot1_mean[step1_vct.at(iStep1)][step2_vct.at(iStep2)] << " +/- " << tot1_sigma[step1_vct.at(iStep1)][step2_vct.at(iStep2)] << std::endl; // << " :: res = " << tot1_sigma/tot1_mean << std::endl;    
            gStep1Tot1[iStep1]->SetPoint(iStep2, step2_vct.at(iStep2), tot1_mean[step1_vct.at(iStep1)][step2_vct.at(iStep2)]);
            gStep2Tot1[iStep2]->SetPoint(iStep1, step1_vct.at(iStep1), tot1_mean[step1_vct.at(iStep1)][step2_vct.at(iStep2)]);
            */
        }
    }
    
    
    int selStep1 = 0;
    int selStep2 = 0;
    
    TCanvas * cArrayTots = new TCanvas ("cArrayTots", "cArrayTots", 1600, 300);
    cArrayTots->Divide(8, 2);
    for (int chId = 0; chId<16; chId++)
    {
        cArrayTots->cd(chId+1);
        hTot[step1_vct.at(selStep1)][step2_vct.at(selStep2)][myChList[chId]]->Draw();
        hTot[step1_vct.at(selStep1)][step2_vct.at(selStep2)][myChList[chId]]->GetXaxis()->SetTitle("tot [ns]");
        hTot[step1_vct.at(selStep1)][step2_vct.at(selStep2)][myChList[chId]]->GetYaxis()->SetTitle("Counts");
	if ( myChList[chId] == 143 )
	  {
	    cArrayTots->SaveAs(Form("hTot_array_chId%.3d.pdf", myChList[chId]));
	  }
    }
    
    
    TCanvas * cArrayTimes = new TCanvas ("cArrayTimes", "cArrayTimes", 1600, 300);
    cArrayTimes->Divide(8, 2);
    for (int chId = 0; chId<16; chId++)
    {
        cArrayTimes->cd(chId+1);
        hTime[step1_vct.at(selStep1)][step2_vct.at(selStep2)][myChList[chId]]->Draw();
        hTime[step1_vct.at(selStep1)][step2_vct.at(selStep2)][myChList[chId]]->GetXaxis()->SetTitle("Time [ps]");
        hTime[step1_vct.at(selStep1)][step2_vct.at(selStep2)][myChList[chId]]->GetYaxis()->SetTitle("Counts");
	if ( myChList[chId] == 143 )
	  {
	    cArrayTimes->SaveAs(Form("hTime_array_chId%.3d.pdf", myChList[chId]));
	  }
    }
    
    
    TCanvas * cArrayCTR_UD = new TCanvas ("cArrayCTR_UD", "cArrayCTR_UD", 1600, 300);
    cArrayCTR_UD->Divide(8, 1);
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
