This README file describes all the outputs from BOREALIS
Version 1: May 17, 2021
It is recommended that you read the BOREALIS README, Input files, and BOREALIS papers before reading this README. 

Table of Contents:
	- Quicklook
	- RFI
	- Processing
	- Post-processing
	- TO DO: De-dispersion and FFT 

------------
Quicklook: 
------------
 	1.) This commands outputs a pdf of the dynamic spectra of the observation. 
    		- Each page is one step of N time (specified in the Input file). 
    		- If plotting the raw data, normalized data, or the data divide by a slope, this will be described in the headers 
    		- The last page is the full quicklook for the entire file 
    		-- Example name of output for LOFAR: L776035_SAP000_B001_S0_QuickLook.pdf           
    		-- Example name of output for NENUFAR: 20200411_TauBoo_QuickLook_beam0_Pol3.pdf     *Pol3=Stokes-V
		-- Example name of output for DYNSPECMS: L352758_Quicklook_beam0_Pol3.pdf

------------
   RFI: 
------------
  This command outputs several files:
 	 1.) PDF showing steps of the RFI processing
		- For the first time step of N spectra, every processing step of the RFI is plotted. 
           This includes: raw data and normalized data. Integrated spectra and time series for both.  
           Each RFI step in order (Patrol, LE_SIG, SUM, PEX)
         - For every other time step, only the normalized data with RFI set as the median value and RFI mask are shown
         -- EXAMPLE output for LOFAR: L776035_pol4_beam0_RFI.pdf
         -- EXAMPLE output for NENUFAR: 20200411_TauBoo_pol4_beam1_RFI.pdf   
         -- EXAMPLE output for DYNSPECMS: L352758_pol0_beam1_RFI.pdf
	2.) RFI mask
		- This file is saved as an idl save file 
		- The format of this file depends on whether the file is saved in bits or bytes. Saving in bits saves space
          - For NENUFAR, the default is to save in bytes because the mask from the pre-processing steps is in bytes. 
          
         If files is in bytes: 
		 -- EXAMPLE output for NENUFAR: mask_full_L545211_TauBoo_pol0_beam0.sav
		 Contents of output: 
			- JD0    		; Reference JD time at start of observation 
			- P2_FULL		; Full RFI mask in bytes [nt,n f]      
			- XF 		; freq ramp [nf]
			- XT			; time ramp [nt]

	    If mask is in bits (ONLY for LOFAR and DYNSPECMS):
          There are two options to save and use mask in bits. 
		 These options are specified under General Settings in the input file:"
			Read/Save Mask in Bits (0=NO, 1=YES):       - Read/save the mask in bits (saves disk space)
			Manipulate the mask in Bits (0=NO, 1=YES):  - Always keep in the mask in bits (this will save memory) 

		If the option to manipulate the mask in Bits is set to 1=YES: 
		 -- EXAMPLE output for LOFAR:	  mask_full_bit_steps_pol0_beam0.sav
  		 -- EXAMPLE output for DYNSPECMS: mask_full_bit_steps_L352758_pol0_beam0.sav
          - To convert the pointer to a bit array use the command pointer_to_array
          - To convert the bit mask into bytes use the command BITARRAY_TO_BYTEARRAY 
		 Contents of output: 
		   	- P2_FULL_BIT_P     ; pointer of the bit mask 
		   	- XF 			 ; freq ramp
			- XT				 ; time ramp
			- XSIZE_STEP		 ; size of original byte array for each step (used in BITARRAY_TO_BYTEARRAY) 
		
		If the option to manipulate the mask in Bits is set to 0=NO: 
		 -- EXAMPLE output for LOFAR:	  mask_full_bit_pol0_beam0.sav
  		 -- EXAMPLE output for DYNSPECMS: mask_full_bit_L352758_pol0_beam0.sav
          - To convert the bit mask into bytes use the command BITARRAY_TO_BYTEARRAY 
		 Contents of output: 
		   	- P2_FULL_BIT       ; full bit mask 
		   	- XF 			 ; freq ramp
			- XT			; time ramp
			- XSIZE.  		 ; size of original byte array (used in BITARRAY_TO_BYTEARRAY) 

	3.)	If RFI is ran with method 2 (default in input file under 'Processing Settings') option for processing, 
		some additional files are created that are used in the processing step (makes the code run faster)
		- For Stokes-V:
		 -- EXAMPLE output for LOFAR: DataMean_L776035_pol4_beam0.sav
		 -- EXAMPLE output for NENUFAR: DataMean_L545211_TauBoo_pol4_beam0.sav
		  Contents of output:
			- DATA_MEAN_SAVE  - Average background at each step calculated using the mean [nt,nf] 
			- PF_FULL         - mask for each step  [nt,nf] 
			- XF              - freq ramp [nf]
			- XT_REBIN        - rebinned time ramp for processing step (nt = number of steps)        
		  
		- For Stokes-I*, 
		  -- EXAMPLE output for LOFAR: Freq_Response_L776035_pol0_beam0.sav 
 		  -- EXAMPLE output for NENUFAR: Freq_Response_L545211_TauBoo_pol0_beam0.sav
		 Contents of output:
             - DATA_MEAN_SAVE   - Average background at each step calculated using the 10% quantile [nt,nf] 
			- PF_FULL         - mask for each step [nt,nf] 
			- XF              - freq ramp [nf]
			- XT_REBIN        - rebinned time ramp for processing step (nt = number of steps)        
		*This file might be created for Stokes-V as well but it not used in the processing step 

------------
Processing:
------------
	1.) PDF showing the creation of the time-freq curve 
		- The 1st pg shows the freq-response averaged over all time
	    	- The next 5 pgs show the time-response for several select freqs fit with a quadratic function
		- The next pgs show the 2-d surface created from the t-f fit, 
            including the dynamic spectra, integrated spectra, and time series
	    -- EXAMPLE output for LOFAR: L779361_pol0_beam0_processing_TF.pdf
         -- EXAMPLE output for NENUFAR: 20200411_TauBoo_pol4_beam1_processing_TF.pdf   
         -- EXAMPLE output for DYNSPECMS: L352758_pol0_beam1_processing_TF.pdf
	2.) PDF showing each step of the rebinning process 
		- Each step is shown with the rebinned data and mask
	    -- EXAMPLE output for LOFAR: L779361_pol0_beam0_processing_rebin.pdf
         -- EXAMPLE output for NENUFAR: 20200411_TauBoo_pol4_beam1_processing_rebin.pdf   
         -- EXAMPLE output for DYNSPECMS: L352758_pol0_beam1_processing_rebin.pdf
	3.) PDF of final processed data
		- Dynamic spectra, integrated spectra, and time series of the mask and final processed data
	4.) sav file of t-f curve
	    -- EXAMPLE output for LOFAR: L776035_pol4_beam0_TimeFreqCorrection.sav
         -- EXAMPLE output for NENUFAR: 20200411_TauBoo_pol4_beam1_TimeFreqCorrection.sav  
	   Contents of output:
		- TIMEFREQCORRECT_ARRAY			-- [4, nf]   where surface(t,f) = a(f) + b(f)*t + c(f)*t^2 
		  TIMEFREQCORRECT_ARRAY[0,*]  	-- freq ramp at raw resolution [nf]
		  TIMEFREQCORRECT_ARRAY[1,*]  	-- a [nf]
		  TIMEFREQCORRECT_ARRAY[2,*]  	-- b [nf]
		  TIMEFREQCORRECT_ARRAY[3,*]  	-- c [nf]
	    - xt                           	-- time ramp at raw resolution [nt]         
	    - For Stokes-V: VMean 		  	-- Average mean across all time [nf]  
	5.) sav file of normalized and rebinned data
        - This file is input into the post-processing part of the code. 
	   - If you would like to run post-processing without running the processing code, 
		the format of the sav file needs to be exactly the same as this file. 
		-- EXAMPLE output for LOFAR:  rebindata_L776035_pol4_beam0.sav
         -- EXAMPLE output for NENUFAR:  rebindata_20200411_TauBoo_pol4_beam1.sav 
         -- EXAMPLE output for DYNSPECMS: rebindata_L352758_pol0_beam1.sav
 	   Contents of output:
		- DATA_REBIN      	- rebinned and processed data [nt,nf]  
	    	- MJD_REBIN       	- rebinned MJD time (only for LOFAR) [nf]
	    	- P2_REBIN        	- rebinned mask [nt,nf] 
	     - UT_REBIN       	- rebinned UT time (only for LOFAR) [nt]
		- XF_REBIN        	- rebinned freq (MHz) [nf]
		- XT_REBIN        	- rebinned time (s)  [nt] 

----------------
Postprocessing:
----------------
	1.) PDF for post-processing plots 
 		- These plots include the dynamic spectra, Q1, Q2, Q3, and Q4 for ON, OFF, and ON-OFF. 
		- The code creates outputs for the polarization, thresholds, 
		  and freq range covered (by default: full freq, 1st half, 2nd half) 
		- For Stokes-V, |V|, V+, and V- are ran separately. 
		-- EXAMPLE output for LOFAR (Elliptical Correction) : L779361_pol4_14.7to26.7MHz_thresT0.10_thresF0.10_postprocessing_AbsV_EllC_Output.pdf
		-- EXAMPLE output for LOFAR (No El Correction): L779361_pol4_14.7to26.7MHz_thresT0.10_thresF0.10_postprocessing_AbsV_Output.pdf
	2.) Txt file of PP summary numbers 
		- A text file including all the summary numbers for Q1 and Q4
		-- EXAMPLE output for LOFAR (Elliptical Correction): L779361_pol4_14.7to26.7MHz_thresT0.10_thresF0.10_postprocessing_AbsV_EllC.txt
		-- EXAMPLE output for LOFAR (No Elliptical Correction): L779361_pol4_14.7to26.7MHz_thresT0.10_thresF0.10_postprocessing_AbsV.txt
	3.) Save file of Q numbers: 
		- save file with all the Q numbers to reproduce the pdf plots
	    -- EXAMPLE output for LOFAR (Elliptical Correction): L779361_pol4_26.7to38.6MHz_thresT0.10_thresF0.10_postprocessing_AbsV_EllC_PP_Qs.sav
 	   Contents of output:
		Q1a             ; [beams, nt]  beams= ON, OFF1, OFF 											
		Q2              ; [beams, nnt]  nnt= rebinned time from processing step					
		Q4              ; [Beams, Nthres, Q4#] 									
		Q4_GAUSS        ; Gaussian curves for Q4  [Beams, Nthres, Q4#]
		THRESHOLD       ; Threshold values for Q4 [Nthres]
		XF              ; Rebinned Freq ramp [nf] 								
		XT              ; Q1 (usually 2 mins) Time ramp [nt] 
	3.) Save file of Q1b numbers (not standardized): 
		- save file with all the Q1a numbers to reproduce the pdf plots
		- TO DO: Need to make consistent and more readable outputs 
	    -- EXAMPLE output for LOFAR: L779361_pol4_26.7to38.6MHz_thresT0.10_thresF0.10_postprocessing_AbsV_EllCFREQPLOT.sav
	   Contents of output:
		ERROR           ; error in Q1a (ON or OFF)     
		ERROR_DIFF      ; error in Q1a diff (ON - OFF)  
		ERROR_I         ; error in Q1b (ON or OFF)     
		ERROR_I_DIFF    ; error in Q1b diff (ON - OFF) 
		UT_REBIN        ; rebin time in UT 	[nt]
		X2_Q1           ; Q1a OFF			[nt] 
		X_Q1            ; Q1a ON			[nt] 
		XD_Q1           ; Q1a Diff (ON-OFF) 	[nt]	
		XFQ1_DIFF       ; Q1b Diff (ON-OFF)	[nf] 
		XF_Q1           ; Q1b ON 			[nf]         
		X2F_Q1          ; Q1b OF			[nf]			
		XF_REBIN_Q1     ;Freq ramp for Q1b 	[nf] 
	4.) TO DO: Combine two save files together and consolidate 

-----------
TO DO: De-dispersion and FFT  
-----------

BeamfOrmed Radio Emission AnaLysIS (BOREALIS) pipeline
By    Jake D. Turner (University of Virginia/Cornell University)
      Jean-Mathias Griessmeier (CNRS)
      Philippe Zarka (OBSPM)
Beta Version 4.0 (March 23, 2021) - BOREALIS  Version 4 