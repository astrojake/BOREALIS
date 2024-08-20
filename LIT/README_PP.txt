README for Post-Processing Code in LOFAR Beam-formed Pipeline 

This README describes two codes: 
 - lofar_processing.pro: The main PP code
 - lofar_postprocess_textfile: Program to create/read the quantitative txt file created by the main program. 

-------------------------------------------------------------------------------------------------------------
Code name: lofar_processing.pro
Code location: LOFAR_Bf_PipelineV2/Main

INPUTS:  
;         save_file_name ; save file name used for the sav files (e.g. files named rebindata_"save_file_name"_beamXX.sav)
;	  Date           ; Date of observation (e.g. In LOFAR format: L570725)  
;         Target_Beam    ; Beam Number of Target (e.g. 0 for ON BEAM in 55Cnc data)
;         Sky_Beam       ; Beam Number of Sky (e.g. 2 for OFF BEAM in 55Cnc data)
;         S              ; Polarization (see below for more info; 0 = I, 1 = Abs(Q), 2 = Abs(U), 3 = Abs(V), 4 = Vprime, 5 = L)
;         rebin_time     ; Orginal rebin time in sav file (secs)
;         rebin_freq     ; Orginal rebin frequency in sav file (MHz)
;         rebin_time2    ; New rebin time (secs) --- Will be timescale used in Q2 and Q3 
;         rebin_freq2    ; New rebin freq (MHz)  --- Frequency interval for int spectrum Q1b
;         fmaxs          ; Max Freqs to search (e.g. [50,40,30]) (MHz)
;         fmins          ; Min Freqs to search (e.g. [40,30,20]) (MHz)
;         threshold      ; Threshold values for Q4 (in sigmas) (e.g. array from 1-6 in steps of 0.1)
;
;REQUIRED FLAGS:
;         Q1             ; Run and plot Q1 (0=NO, 1 = YES)
;         Q2             ; Run and plot Q2 (0=NO, 1 = YES)
;         Q3             ; Run and plot Q3 (0=NO, 1 = YES)
;         Q4             ; Run and plot Q4 (0=NO, 1 = YES)
;         threshold_p    ; Mask Theshold for indiviual pixels (e.g. 0.10, anything below 0.10 will be set to zero and not used)
;         threshold_t    ; Mask threshold for integrated time series (e.g. x = 0.10, then 1-x = 90%, everything above 90% is kept)
;         threshold_f    ; Mask threshold for integrated spectrum mask (e.g. x = 0.10, then 1-x = 90%, everything above 90% is kept)
;         VERBOSE        ; Verbose 
;         PS 		 ; Create PS (pdf) files
;
;OPTIONAL FLAGS:
;         TI             ; Interval for searching for a signal in Q3 in seconds (e.g. 120 seconds) (Default: 120)
;         slope          ; Slope value comparing beams for Q3f/Q4f (Default: 2) 
;         max_threshold  ; Max sigma value for Q2 (sigma) (Default: 10)
;         tau_Qplot      ; Sigma value for Q3 plot vs. time (sigma) (Default: 2)
;         Qwin           ; Window size for high pass filter to create Q2 (units of time steps) (Default: 10)
;         Default_Time   ; Run Qs through all time (Default: Yes)
;         QTmin          ; If Default_Time is not set, Min time to run Qs (days) 
;         QTmax          ; If Default_Time is not set, Max time to run Qs (days)
;         READ_GAUSS     ; Keyword READ_GAUSS is NOT set, program creates Gaussian sav set (See more details below)
;         File_G         ; Filename of Gaussian sav set (see more details below)
;         gauss_steps    ; Number of Gaussian steps to take for Gaussian Q4 comparison (Default: 10,000)
;         OnBeam_Label   ; Label for On Beam (e.g. 55Cnc Beam) (Default: ON Beam)
;         OffBeam_Label  ; Label for Off Beam (e.g. OFF-Beam 2) (Default: OFF Beam) 
;         NOELLCORR      ; Run non-ellipical correction on Q2 (0=NO, 1 = YES) (Default: No)
;         save_Q         ; Save Q values in save set (saves Q1,Q2,Q3,Q4)
;         UT 		 ; Plot time axis in UT time (default time axis is in days) -useful for papers
;         TARGET_ONLY    ; Only plot the target (no OFF Beam): Only plots Q1 for target beam
;	  Extra_Beam1    ; Beam # of extra beam 1 to plot (e.g 2 for 55Cnc data) -- Plots Only for Q1
;         Extra_Beam2    ; Beam # of extra beam 2 to plot (e.g. 3 for 55Cnc data)-- Plots Only for Q1
;
;FLAGS UNDER CONSTRUCTION:
;         LIST           ; This is the elements of each date, so that sevearl dates can be analyzed  at same time (NOT CURRENTLY USED)
;         PHASE 	 ; plot phase instead of time (NOT USED yet, will be useful for exoplanets) 
;         Planet_Name    ; Name of the planet  (NOT USED)
;         PERIOD         ; Period of planet (days)                                    --only used if phase is to plotted 
;         T0             ; Time of 0 phase for the exoplanet from first transit (MJD) --only used if phase is to plotted 
;
;Flags used in specials cases:
;         RUN_JUPITER    ; Run Jupiter tests (NOTE, not general)
;
;File Input:
;         The only files input to the code are the rebined data files. 
;         The name for these files are rebindata_"save_file_name"_beamXX.sav,
;            where "save_file_name" is input into the code and XX is the beam number (0,1,2,3)
;         This rebinned save set should contain the following variables:
;              xt_rebin     ;Time array in seconds since the start of the observation (rebinned) 
;              MJD_rebin    ;Time array of the observations in MJD (rebinned)
;              UT_rebin     ;UT time array of the observations (rebinned)
;              xf_rebin     ;frequency array of the observations (rebinned) 
;              data_rebin   ;Rebinned Dynamic Spectra (t,f) 
;              p2_rebin     ;Rebinned RFI mask (t,f) -same dimenions as data_rebin
;
;Notes on format for save_file_name:
;         The format for the save_file_name can be anytime. In LOFAR pipeline, they are normally long (polarization, time of obs, RFI info included) 
;         
;Polarization Input (S): The polarization input into the code will perform certain functions
;              S = 0; Stokes-I; Input dyn spec is I, I is ran through PP 
;	       S = 1; Stokes-Q; Input dyn spec is Abs(Q), Abs(Q) is ran through PP
;              S = 2; Stokes-U; Input dyn spec is Abs(U), Abs(U) is ran through PP
;	       S = 3; Stokes-V; Input dyn spec is Abs(V), Abs(V) is ran through PP
;              S = 4; Stokes-V; Input dyn spec is Vprime, Abs(V), V+, V- is ran through PP 
;              S = 5; L; Input dyn spec is L, L is ran through PP 
;
;GAUSSIAN NOTES:
;	       The gaussian Q refernce curve is always ran. Default is to run 10000 or the keyword gauss_steps. 
;              3 options are available: 
;               (1) Keyword READ_GAUSS is NOT set (or equal to 0). Program creates the new Gaussian file. 
;                   - Takes about 30-40 mins to perform this task. Other options takes less than 1 min. 
;                   - If multiple frequency sections are chosen, the Guassian is only calculated on the 1st run and then used in the following sections
;               (2) Keyword READ_Gauss is set (or equal to 1) and FILE_G is not set. Automatically creates the name and reads in a saved Gaussian sav file 
;               (3) Keyword READ_Gauss is set (or equal to 1) and FILE_G is set. Reads in a saved Gaussian sav file called FILE_G. 
;                   - The format of FILE_G is as follows: Q4_Gauss_"Ntime"_"rebin_time2".sav, where "NTime" is the number of time elements in Q2 and 
;                       "rebin_time2" is the rebinned time in secs. (eg. Q4_Gauss_10787_1.00000.sav)
;              
;DEFAULT OUTPUTS:
;       Panel PDF of PP: "save_file_name"_"Fmin"to"Fmax"MHz_thresT"threshold_t"_thresF"threshold_f"_postprocessing_"S"_EllC_Output.pdf
;                        (e.g. L570725_50.0to60.0MHz_thresT0.10_thresF0.10_postprocessing_I_EllC_Output.pdf) 
;       Text File with PP values: "save_file_name"_"Fmin"to"Fmax"MHz_thresT"threshold_t"_thresF"threshold_f"_postprocessing_"S"_EllC.txt 
;                        (e.g. L570725_50.0to60.0MHz_thresT0.10_thresF0.10_postprocessing_I_EllC.txt)
;
;OPTIONAL OUTPUTS (If keyword save_Q is set):
;       Sav file with Qs (File will contain xt,xf and the Qs): "save_file_name"_"Fmin"to"Fmax"MHz_thresT"threshold_t"_thresF"threshold_f"_postprocessing_"S"_EllC_Qs.sav
;                         (e.g. (e.g. L570725_50.0to60.0MHz_thresT0.10_thresF0.10_postprocessing_I_EllC_Qs.sav)
------------------------------------------------------------------------------------------------------------------------------------------------------------------------

Code name: lofar_postprocess_textfile
Code location: LOFAR_Bf_PipelineV2/Main/Secondary

INPUTS: 
       filename: Filename of txt file of interest (e.g. 570725_50.0to60.0MHz_thresT0.10_thresF0.10_postprocessing_I_EllC.txt) 

FLAGS:        
       Date: Date 
       S   : Polarization 
       fmin: Min Frequency
       fmax: Max Frequency 
       BT  ; Beam Number of Target Beam (e.g. 0) 
       BSky; Beam Number of Sky Beam (e.g. 2)
       Q1aV: Q1a value
       Q1bV: Q1b value
       Q4aV: Q4a value
       Q4bV: Q4b value
       Q4cV: Q4c value
       Q4dV: Q4d value
       Q4eV: Q4e value 
       Q4fV: Q4f value 
       Detect: If set, code tells you if there is a detection in Q4f (if the ratio to the Gaussian is gt 1)
       WRITE: Write the txt file 
       READ:  Read the txt file

See PP_Quantitative_Overview.pdf for a description of each Q quantitative values. Located in LOFAR_Bf_PipelineV2/LIT/ 