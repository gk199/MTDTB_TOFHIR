# MTDTB_TOFHIR

Analysis code for the SiPM array FNAL testbeams. To compile,

    cmsenv
    g++ -Wall -o plotTOFHIR_to_trigger plotTOFHIR_to_trigger.C functions.hh `root-config --cflags --glibs` -lSpectrum
    g++ -Wall -o plotTOFHIR_to_trigger_FebTB plotTOFHIR_to_trigger_FebTB.C functions.hh `root-config --cflags --glibs` -lSpectrum
    g++ -Wall -o plotTOFHIR_to_triggerFebTBv2 plotTOFHIR_to_trigger_FebTB_v2.C functions_v2.cc `root-config --cflags --glibs` -lSpectrum

and to run,

    ./plotTOFHIR_to_trigger <first run number> <last run number> <path to data>
    ./plotTOFHIR_to_trigger_FebTB <first run number> <last run number> <path to data>
    ./plotTOFHIR_to_triggerFebTBv2 <first run number> <last run number>
    ./plotTOFHIR_to_triggerFebTBv2 24702 24712

## June Testbeam Details
Data from the June 2019 FNAL Testbeam is stored in 

    /eos/cms/store/group/dpg_mtd/comm_mtd/TB/MTDTB_FNAL_Jun2019/TOFHIR/RecoData/RecoWithTracks/v2/
        
For the array with both sides read out, the range of run numbers is 15554 to 15673. (This array is 3x3x50 mm^3 with both sides readout with 3.2 mm SiPMs). For the array matched to the SiPMs (3x3.2x57 mm^3 with one side readout), the run numbers are 15674 to 15716. The [run logbook](https://docs.google.com/spreadsheets/d/1ilOMxOy2Qlut1EweUQIKJa14cDbHkf-NVMCfdacFLVk/edit#gid=555033895 "June 2019 Run Logbook") has the channel mappings and specifics about the configurations for each run.

For the second array (where only even channel numbers are read out), runs 15676-15710 have step 1 (overvoltage) = 6, and step 2 (threshold) = 0. Run 15692 is very well reconstructed and good to run over.

For Run 15692, this is with vth1 = 0 and vth2 = 30 (these are the starting time and falling edge of how the ToT is determined), and therefore the correction function for this is applied. The fit used is either:

    pol3 with TF1 *fCorr = new TF1("fCorr","[0]*x + [1]*x*x + [2]*x*x*x", 0, 1000);

    expo with TF1 *fCorr = new TF1("fCorr","[0] * ( exp([1]*x) - 1 )", 0, 1000);

with [pol3 plotted here](https://malberti.web.cern.ch/malberti/MTD/Lab5015/TOFHIR/NonLinearityToT/pol3/vth1_0_vth2_30/c_correction.pdf "Martina Malberti Tot non linearity pol3") and [expo plotted here](https://malberti.web.cern.ch/malberti/MTD/Lab5015/TOFHIR/NonLinearityToT/vth1_0_vth2_30/c_correction.png "Martina Malberti Tot non linearity exponential"). This ToT linearization is applied for the Landau fit and the cross talk studies. 

## February Testbeam Details
Data from February 2020 FNAL Testbeam is stored in
     
     /eos/cms/store/group/dpg_mtd/comm_mtd/TB/MTDTB_FNAL_Feb2020/TOFHIR/RecoData/v1/RecoWithTracks/

[Run logbook](https://docs.google.com/spreadsheets/d/1BTyPMHOUD97ctFpJEPvZQs3BoDRCWzZvoSHltiQdLSg/edit#gid=1048771982 "February 2020 Run Logbook") has details about configuration for each run. Additional code for plots in 

     /afs/cern.ch/user/m/mtd/Lab5015Analysis/main

with plots [here](http://miptimingdetector.web.cern.ch/miptimingdetector/). Mostly run in QDC mode where tot, qfine, and energy are read out. Details for [TOFPET](https://drive.google.com/file/d/15z_Hjv814W3Wo-l9Emn5No1tnk_T7Np1/view).

## Plotting Cross Talk
Cross talk: blue plot shows fractional energy in a bar when a MIP is deposited somewhere in the array, green shows when the MIP is in the central bar, red shows when MIP is in neighboring bar, and purple when the MIP is in a second to neighboring bar.