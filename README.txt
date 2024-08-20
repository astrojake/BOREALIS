BeamfOrmed Radio Emission AnaLysIS (BOREALIS) pipeline
By    Jake D. Turner (University of Virginia/Cornell University)
      Jean-Mathias Griessmeier (CNRS)
      Philippe Zarka (OBSPM)
Beta Version 1.0 (November 30, 2017)
Beta Version 2.0 (June 2, 2018)   
Version 3.0 (September 12, 2018) 
Beta Version 4.0 (March 23, 2021) - Version 4 (named changed to BOREALIS); 
				   Updated with NENUFAR (new fits files) and DynspecMS
GNU General Public License v3

FOLDER ORGANIZATION: 
	-BOREALISv4
          - README.text 				 --> This file 
	  - README_InputFile.dat			 --> Explanation of the Input file (line by line)
	  - Input_File.dat				 --> Example Input file for the code
	  - Input_File_NenuFAR.dat		         --> Example Input file for the code on NenuFAR data
          - BOREALISv4/Main                    		 - All the main codes for the pipeline
		- BOREALISv4/Main/IDL_pro_util  	 - IDL codes by our team
                - BOREALISv4/Main/Secondary    		 - IDL codes created for this code by our team
          - BOREALISv4/LIT                      	- Literature (with refs in papers) describing the code     
          - BOREALISv4/sub		         	- Folder for non-standard IDL packages
		- BOREALISv4/sub/astron.dir     	- IDL astro lib
		- BOREALISv4/sub/catalyst       	- Used by IDL coyote lib
		- BOREALISv4/sub/coyote         	- IDL coyote lib
	  - BOREALISv4/IMAGING_BF		 	- Test directory of codes for BF data from imaging

README FILES:
	The download of the code has two REAME files (make sure to read both):
	  -README.txt           ——>This file 
	  -README_InputFile.dat ——>Explanation of the Input file (line by line)
          - Read the LIT folder 

LIT FOLDER (located in BOREALISv4/LIT) 
		- Turner et al. 2017 (PRE 8) and refs within
		- Turner et al. 2019 (Jupiter paper)
		- Turner et al. 2021 (1st results and tentative detections on Tau Bootis)
                - LOFAR-I-V-processing-FINAL.docx : Final (near) Overview of Processing
		- PP_Quantitative_Overview.pdf: Overview of the PP Quantitative txt file
		- README_PP.txt: Description of how to run the stand alone PP code.
                                 It also includes how to read in the text file created by the PP. 

WHAT IT DOES: 
	     The code does the following on beam-formed data from LOFAR: 
		- Creates a Quicklook
 		- Finds an RFI mask
		- Finds the time-frequency (t-f) response function 
		  - Method 1: Average of the spectrum - described in Turner et al 2017 (PRE 8)
		  - Method 2: Average spectrum scaled from the 10% quantile (or average) over each time slice (described in Turner et al. 2019)
		- Applies the RFI mask and normalizes the data by time-frequency response function
		- Rebins the data and mask	
		- Runs post-processing in order to detect bursty signals 

		The parts of the code are as follows:
		- Quicklook
		- Processing 
		  - RFI mask, find t-f response function, apply mask & t-f to data, and rebin 
		- Post-processing

		The codes runs in the following order:
		 - QuickLook
		 - RFI for all beams
		 - Combine masks from beams if set
		 - Processing for all beams
		   - if de-dispersion is set, this runs parallel to the normal processing
		   - FFT code runs after processing code is done 
		 - Post-processing 

HOW TO USE:
          - The code runs by editing an Input file and then running the code. No other changes besides the input file are needed.
     LOFAR:     
	- Hdf5 files with LOFAR headers are used for the code
	- The data should be located in ROOT_FOLDER/LOFARFORMAT/ (both the ROOT_FOLDER and LOFARFORMAT are defined in Input file)
		- e.g. ROOT_FOLDER = /databf/LOFAR/EXO/LC5_DDT_002/ and LOFARFORMAT = L429868
	  - You have to have a symbolic link to the data located at ROOT_FOLDER/LOFARFORMAT in whatever directory you are working in. 
	    For hdf5 files you have to link both the .h5 and .raw files
	 	- You can create symbolic link using in a terminal > ln -s ROOT_FOLDER/LOFARFORMAT/FILENAME 
		   (e.g. >> ln -s /databf/LOFAR/EXO/LC5_DDT_002/L433872/raw/L433872_SAP000_B003_S0_P000_bf.raw .)
    NenuFAR: 
	- Runs on fits files created after initial processing by NenuFAR team (L1 data)
     	- You have to have a symbolic link to the data located at ROOT_FOLDER in whatever directory you are working in. 
		e.g: ROOT_FOLDER: /databf2/nenufar-tf/ES02/2021/02/20210210_060000_20210210_100000_CR_DRA_TRACKING/L1/
	- In the input file for the "Filename to run" only input the name and not the time:
 		e.g. files are TAU_BOOTIS_JOINT_CAMPAIGN_TAU_BOOTIS_TRACKING_20200411_195938_#.spectra.fits
		     only put TAU_BOOTIS_JOINT_CAMPAIGN_TAU_BOOTIS_TRACKING_20200411 for the "Filename to run"
   DynspecMS: 
	- To be updated -- in progress 

TO COMPILE:
	Put the procedures in your IDL path
	Example (procedures located at /home/codes/): 
          - Type from command line in IDL: PREF_SET, 'IDL_PATH', '+/data/jake.turner/BOREALISv4:<IDL_DEFAULT>', /COMMIT
          - Doing this once updates the IDL_PATH, with the requested directory + the standard libraries ('<IDL_DEFAULT>').
          - This setting then stays forever; no need to do this after every start of IDL. 
          - The .bashrc should NOT contain the definitions of the type IDL_PATH.

TO RUN: 
 	- Edit the input file (e.g. Input.dat)
	- In IDL type:
	 IDL> .r BOREALISv4
	 IDL> BOREALISv4,’Input.dat’

INPUT REQUIREMENTS:
	The only input into the code is the input file. 
	An example input file (Input_File.dat) can be found in the main downloaded directory Loar_BF_PipelineV1. 
        A description of each line in the input file can found in README_InputFile.txt 

POLARIZATION INPUT (S)
 LOFAR:
 The polarization input into the code will perform certain functions:
 S = 0; Stokes-I; Input dyn spec is I; I is ran through PP 
 S = 1; Stokes-Q; Input dyn spec is Abs(Q), Abs(Q) is ran through PP
 S = 2; Stokes-U; Input dyn spec is Abs(U); Abs(U) is ran through PP
 S = 3; Stokes-V; Input dyn spec is Abs(V); Abs(V) is ran through PP
 S = 4; Stokes-V; Input dyn spec is Vprime; Abs(V), V+, V- is ran through PP 
 S = 5; L; Input dyn spec is L, L is ran through PP 

 For NenuFAR: only S= 0, S = 3; S = 4; and S=5 work   

Special notes about NenuFAR data:
 - The mask from the L1 fits in in bytes. Therefore, you can't run or save the RFI in BOREALISv# in bits like you can for LOFAR
 - The pointing jumps make RFI more tricky. Therefore, it is not recommended to use the 10% quantile for normalization. Also, make sure the steps are shorter than the jump time (more in Input file)
 - Default RFI setting will be different for NenuFAR as they were tuned for the LOFAR data (time and freq resolution and environment)   
 - See Input_File_NenuFAR.dat for more recommendations 

IDL REQUIREMENTS:
	IDL v 8.4 or later
	May be compatible with earlier versions, but they have not been tested

TIME REQUIREMENTS:
    LOFAR Data:
	- Running the code on 4 hours of observations takes:
	  - QuickLook: ~4 hours each beam 
          - RFI: ~20 hours each beam
	  - Processing: ~5 hours each beam
	  - Post-processing: ~30 minutes if running gaussian, ~5 minutes if reading gaussian from a file —> (both with Ell Correction & all wavelength ranges)
                              30 secs: reading gaussian from a file and one wavelength range
    NenuFAR Data:
	 - Takes less than an 1 hour for all beams 


PROGRAM DESCRIPTIONS AND REQUIREMENTS:
	-All the needed programs and subfiles needed for the code can be found in BOREALISv#/Main/. 
	-All the programs developed by our team are located in BOREALISv#/Main and BOREALISv#/Main/IDL_PRO_UTIL
        -The following external libraries are used (located in BOREALISv#/sub):
	 - coyote (idl coyote libraries)
         - catalyst (used with idl coyote libraries) 
	 - idlastro lib
	-IDL Libs used: 
	  - Default packages 
	  - hdf5 library (in default IDL packages) 
	-MAIN CODES w/ descriptions located at LOFAR_Bf_PipelineV1/Main (Each code has a lengthy description of inputs/outputs in the header):
	   For outputs described below (e.g. SAVE_FILE_NAME =   L429868_pol0_beam0)
	 - BOREALISv4.pro
	    - This is only file needed to run the program. All other programs are called as subroutines in this program. 
	    - The only input into this file is the name of the Input file
	 - BOREALIS_rfi.pro: 
	   - Runs RFI mitigation on the data. 
	   - The outputs are the RFI mask and the 10% quantile at each step needed for the t-f (Method 2).
	   - OUTPUTS: 
		- rfi plots mask for each beam (e.g. SAVE_FILE_NAME_RFI.pdf ) 
		- rfi mask for each beam (e.g. mask_full_bit_steps_SAVE_FILE_NAME.sav)
		- If method 2 for t-f curve: freq response curves (e.g. Freq_Response_SAVE_FILE_NAME.sav)
	 - BOREALIS_processing.pro 
 	   - Finds t-f normalization curve 
	   - Applies t-f norm cure and RFI mask to the data 
	   - Rebins the data and mask 
	   - OUTPUTS: 
		- Plots for various aspects of processing
		- Plots for creation of t-f curve (SAVE_FILE_NAME_processing_TF.pdf)
		- Plots for rebinning (SAVE_FILE_NAME_processing_rebin.pdf)
		- Plots of final rfi-masked, normalized, and rebinned data (SAVE_FILE_NAME_processing_final.pdf)
		- Save file of t-f curve (SAVE_FILE_NAM_beam#_TimeFreqCorrection.sav) 
		- If de-dispersion is set a sav file with data (SAVE_FILE_NAME_disp.sav)
	 - BOREALIS_postprocessing.pro 
	   - Evaluates the observable quantities (Qs) 
	   - OUTPUTS: 
		- (e.g. SAVE_FILE_NAME2 =   L429868_pol0)
		- Plots for post-processing (e.g. SAVE_FILE_NAME2_26.1to73.7MHz_thresT0.10_thresF0.10_postprocessing.pdf) 
		- If El Correction (e.g. SAVE_FILE_NAME2_26.1to73.7MHz_thresT0.10_thresF0.10_postprocessing_EllC.pdf) 
                - Save file of Qs (e.g. SAVE_FILE_NAME2_26.1to73.7MHz_thresT0.10_thresF0.10_postprocessing_PP_Qs.sav
		- Gaussian file (e.g. Q4_Gauss_15572_1.00000.sav)

	-****Will add more codes later*** 	

COMMON RUNNING MODES:
	Run Quicklook: To run a quicklook, just edit the Input_file with YES for ‘Run QUICKLOOK’ and edit the section ‘QuickLook Settings’
	Run All      : To run everything for all beams, edit in the the Input_file with YES for 
                       ‘Run RFI’, ‘RUN Processing’, ‘Run Post-processing’ and edit the main settings

USING POST-PROCESSING CODE AS STAND-ALONE PROGRAM:
 

IMAGING DATA CODES (NOT COMPLETE!):
	- The codes to run the imaging data dynamic spectrum are not general yet. 
	- These codes can be found in BOREALISv# under /IMAGING_BF
    TO RUN: 
	- Convert fits files to sav file format (edit image_createsav.pro)
	- edit dynspec_imaging.pro
	- IDL> .r dynspec_imaging.pro 
        - IDL> dynspec_imaging,Target_Beam=Target_Beam,Sky_Beam=Sky_Beam
   
    MAIN CODES:
		imagine_createsav.pro   - creates the sav files from the fits files
		dynspec_imaging.pro 	- runs post-processing on the sav files (uses LOFAR_postprocessing.pro)

Please cite the following papers if you use this code for your research: 
1.) Turner J.D., Griessmeier J.-M., Zarka P., Vasylieva I. 2017. The search for radio emission from exoplanets using LOFAR low-frequency beam-formed observations: Data pipeline and preliminary results for the 55 Cnc system. Planetary Radio Emissions VIII. arXiv:1710.04997.
2.) Turner J.D., Griessmeier J.-M., Zarka P., Vasylieva I. 2019. The search for radio emission from exoplanets using LOFAR beam-formed observations: Jupiter as an exoplanet. A&A. 624A. 40.
3.) Turner J.D., Zarka P., Griessmeier J.-M., et al. 2021. The search for radio emission from the exoplanetary systems 55 Cancri,Andromedae, and τ Boötis using LOFAR beam-formed observations. A&A. 645. 59.

CONTACT INFO: 
If you find a bug please let me know: 
Jake D. Turner 
astrojaketurner@gmail.com

* Copyright 2017 Jake Turner <astrojaketurner@gmail.com>
 * These programs are free 
*  software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, see <http://www.gnu.org/licenses/>
}
