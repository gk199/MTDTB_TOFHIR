# MTDTB_TOFHIR

Analysis code for the SiPM array FNAL testbeams. To compile,

    g++ -Wall -o plotTOFHIR_to_trigger plotTOFHIR_to_trigger.C functions.hh `root-config --cflags --glibs` -lSpectrum

and to run,

    ./plotTOFHIR_to_trigger <first run number> <last run number> <path to data>

Data from the June 2019 FNAL Testbeam is stored in 

    /eos/cms/store/group/dpg_mtd/comm_mtd/TB/MTDTB_FNAL_Jun2019/TOFHIR/RecoData/RecoWithTracks/v2/
        
For the array with both sides read out, the range of run numbers is 15554 to 15673. (This array is 3x3x50 mm^3 with both sides readout with 3.2 mm SiPMs). For the array matched to the SiPMs (3x3.2x57 mm^3 with one side readout), the run numbers are 15674 to 15716. The [run logbook](https://docs.google.com/spreadsheets/d/1ilOMxOy2Qlut1EweUQIKJa14cDbHkf-NVMCfdacFLVk/edit#gid=555033895 "June 2019 Run Logbook") has the channel mappings and specifics about the configurations for each run.

For the second array (where only even channel numbers are read out), runs 15676-15710 have step 1 (overvoltage) = 6, and step 2 (threshold) = 0. Run 15692 is very well reconstructed and good to run over.

For Run 15692, this is with vth1 = 0 and vth2 = 30 (these are the starting time and falling edge of how the ToT is determined), and therefore the correction function for this is applied. The fit used is pol3 with TF1 *fCorr = new TF1("fCorr","[0]*x + [1]*x*x + [2]*x*x*x", 0, 1000); and [plotted here](https://malberti.web.cern.ch/malberti/MTD/Lab5015/TOFHIR/NonLinearityToT/pol3/vth1_0_vth2_30/c_correction.pdf "Martina Malberti Tot non linearity"). This ToT linearization is applied for the Landau fit and the cross talk studies. 