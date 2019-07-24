# MTDTB_TOFHIR

Analysis code for the SiPM array FNAL testbeams. To compile,

    g++ -Wall -o plotTOFHIR_to_trigger plotTOFHIR_to_trigger.C `root-config --cflags --glibs` -lSpectrum

and to run,

    ./plotTOFHIR_to_trigger <first run number> <last run number> <path to data>

Data from the June 2019 FNAL Testbeam is stored in 

    /eos/cms/store/group/dpg_mtd/comm_mtd/TB/MTDTB_FNAL_Jun2019/TOFHIR/RecoData/v1/RecoWithTracks/
        
For the array with both sides read out, the range of run numbers is 15554 to 15673. (This array is 3x3x50 mm^3 with both sides readout with 3.2 mm SiPMs). For the array matched to the SiPMs (3x3.2x57 mm^3 with one side readout), the run numbers are 15674 to 15716. The [run logbook](https://docs.google.com/spreadsheets/d/1ilOMxOy2Qlut1EweUQIKJa14cDbHkf-NVMCfdacFLVk/edit#gid=555033895 "June 2019 Run Logbook") has the channel mappings and specifics about the configurations for each run.
