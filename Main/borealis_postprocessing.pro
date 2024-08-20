;***************************************************************
;      Post-Processing Lofar exoplanet data
;**************************************************************
;***************************************************************
;Author: Jake Turner (UVA)
;        J.M. Griessmeier (CNRS) and Philippe Zarka (Obspm)
;        
;INPUTS:  
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
;         READ_GAUSS     ; Keyword READ_GAUSS is NOT set, program creates Gaussian sav set (See more details below)  --NOT USED ANYMORE
;         File_G         ; Filename of Gaussian sav set (see more details below)                                     --Still used
;         gauss_steps    ; Number of Gaussian steps to take for Gaussian Q4 comparison (Default: 10,000)             -NOT USED
;         OnBeam_Label   ; Label for On Beam (e.g. 55Cnc Beam) (Default: ON Beam)
;         OffBeam_Label  ; Label for Off Beam (e.g. OFF-Beam 2) (Default: OFF Beam) 
;         NOELLCORR      ; Run non-ellipical correction on Q2 (0=NO, 1 = YES) (Default: No)
;         save_Q         ; Save Q values in save set (saves Q1,Q2,Q3,Q4)
;         UT 		 ; Plot time axis in UT time (default time axis is in days) -useful for papers
;         TARGET_ONLY    ; Only plot the target (no OFF Beam): Only plots Q1 for target beam
;	        Extra_Beam1    ; Beam # of extra beam 1 to plot (e.g 2 for 55Cnc data) -- Plots Only for Q1
;         Extra_Beam2    ; Beam # of extra beam 2 to plot (e.g. 3 for 55Cnc data)-- Plots Only for Q1
;         MINT           ; Min threshold vallue used for false-postiive calc on Q4
;         MAXT           ; Max threshold value used for false-postiive calc on Q4

;FLAGS UNDER CONSTRUCTION:
;         LIST           ; This is the elements of each date, so that sevearl dates can be analyzed  at same time (NOT CURRENTLY USED)
;         PHASE 	       ; plot phase instead of time (NOT USED yet, will be useful for exoplanets) 
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
;               (2) Keyword READ_Gauss is set (or equal to 1) and FILE_G is no set. Automatically creates the name and reads in a saved Gaussian sav file 
;               (3) Keyword READ_Gauss is set (or equal to 1) and FILE_G is set. Reads in a saved Gaussian sav file called FILE_G. 
;                   - The format of FILE_G is as follows: Q4_Gauss_"Ntime"_"rebin_time2".sav, where "NTime" is the number of time elements in Q2 and 
;                       "rebin_time2" is the rebinned time in secs. (eg. Q4_Gauss_10787_1.00000.sav)
;              
;DEFAULT OUTPUTS:
;       Panel PDF of PP: "save_file_name"_"Fmin"to"Fmax"MHz_thresT"threshold_t"_thresF"threshold_f"_postprocessing_"S"_EllC_Output.pdf
;                        (e.g. L570725_50.0to60.0MHz_thresT0.10_thresF0.10_postprocessing_I_EllC_Output.pdf) 
;       Text File with PP Values: "save_file_name"_"Fmin"to"Fmax"MHz_thresT"threshold_t"_thresF"threshold_f"_postprocessing_"S"_EllC.txt 
;                        (e.g. L570725_50.0to60.0MHz_thresT0.10_thresF0.10_postprocessing_I_EllC.txt)
;
;OPTIONAL OUTPUTS (If keyword save_Q is set):
;       Sav file with Qs (File will contain xt,xf and the Qs): "save_file_name"_"Fmin"to"Fmax"MHz_thresT"threshold_t"_thresF"threshold_f"_postprocessing_"S"_EllC_Qs.sav
;                         (e.g. (e.g. L570725_50.0to60.0MHz_thresT0.10_thresF0.10_postprocessing_I_EllC_Qs.sav)
;
;Observables: 
;      Q1: total power of each time-interval (TI), itegrated over all frequencies
;      Q2: Normalized time series (high-pass filtered and divided by standard deviation)
;      Q4a: number of peaks per TI > threshold tau, compared to number of peaks < -tau 
;      Q4b: sum of power of peaks per TI > threshold tau, compared to sum of abs(Power) of peaks < -tau 
;      Q4c: Number of peaks per TI > threshold tau and exceeding the cooresponding Off Values by a factor >= 2 
;      Q4d: Sum of power of peaks per TI above a threshold tau and exceeding the corresponding off values by a factor of > 2 
;
;History;  May 23, 2017: Combine dates and old Qs (on OSURC)
;          May 24, 2017: (no LIST) will update when needed
;          May 31, 2018: Include polarization for Vprime
;          Sept 12, 2018:Update Gaussian to run only once and then scale Q values
;
pro BOREALIS_postprocessing,save_file_name,Date,Target_Beam,Sky_Beam,S,rebin_time,rebin_time2,rebin_freq,rebin_freq2,fmaxs,fmins,threshold,$
                          Q1=Q1,Q2=Q2,Q3=Q3,Q4=Q4_run,$
                          threshold_p=threshold_p,threshold_t=threshold_t,threshold_f=threshold_f,$
                          VERBOSE=VERBOSE,PS=PS,$ 
                          TI=TI,slope=slope,max_threshold=max_threshold,tau_Qplot=tau_Qplot,Qwin=Qwin,$
                          Default_Time=Default_Time,QTmin=QTmin,QTmax=QTmax,$
                          READ_GAUSS=READ_GAUSS,FILE_G=file_g,gauss_steps=gauss_steps,$
			                    OnBeam_Label=OnBeam_Label,OffBeam_Label=OffBeam_Label,$
                          ELLCORR=ELLCORR,$
                          save_Q = save_Q,$
                          UT=UT,$
                          TARGET_ONLY=TARGET_ONLY,$
                          Extra_Beam1=Extra_Beam1,Extra_Beam2=Extra_Beam2,$
                          PHASE=PHASE,PERIOD=P,T0=T0,LIST=LIST,Planet_Name=Planet_Name,$
                          RUN_JUPITER=RUN_JUPITER,$
                          MINT=MINT,MAXT=MAXT,$
                          DYNSPECMS=DYNSPECMS,$
                          REMOVE=REMOVE,PAPERPLOTS=PAPERPLOTS,$
                          UTR2=UTR2,XTMIN=XTMIN,XTMAX=XTMAX
  stp          =  SYSTIME(1)        ;start time
  startmem_stp =  MEMORY(/CURRENT)  ;start memory
 
  ;---------------------------------------------------------------
  ;----------------------------------------------------------------
  ;				Inputs
  ;----------------------------------------------------------------
  ;---------------------------------------------------------------- 
  If not(keyword_set(REMOVE)) then REMOVE= 1   ;(1=yes, 0=NO) Remove pdf of separated plots (used to make panel plots)
  If not(keyword_set(PAPERPLOTS)) then PAPERPLOTS= 1   ;(1=yes, 0=NO) creat plot for papers (no labels and jargon) 
    
;----------------------------------------------------------------------------------------------------------------
;  				Don't need to change anything below this line
;----------------------------------------------------------------------------------------------------------------
  panel_ps  =1         ;Always plot ps file in panels
  square = 0
  char_size = 2 
  
  ;Check Flags
  If not(keyword_set(Extra_Beam1)) then ExBeam1_info = 0 
  If not(keyword_set(Extra_Beam1)) then ExBeam2_info = 0 
  
  print,'***********************************************'
  if n_elements(Extra_Beam1) ne 0 OR arg_present(Extra_Beam1) then $
    print, 'Plot ExBeam1: set' else $
    print, 'Plot ExBeam1: not set'
  if n_elements(Extra_Beam2) ne 0 OR arg_present(Extra_Beam2) then $
    print, 'Plot ExBeam2 Source: set' else $
    print, 'Plot ExBeam2 Source: not set'  
  if n_elements(threshold_t) ne 0 OR arg_present(threshold_t) then $
    print, 'Threshold t: set', threshold_t  else $
    print, 'Threshold t: not set'
  if n_elements(threshold_f) ne 0 OR arg_present(threshold_f) then $
    print, 'Threshold f: set', threshold_f else $
    print, 'Threshold f: not set'   
  if keyword_set(Q1) then $
    print, 'Q1: set' else $
    print, 'Q1: not set' 
  if keyword_set(Q2)  then $
    print, 'Q2: set' else $
    print, 'Q2: not set' 
  if keyword_set(Q3) then $
    print, 'Q3: set' else $
    print, 'Q3: not set'
  if keyword_set(Q4_run) then $
    print, 'Q4: set' else $
    print, 'Q4: not set' 
  ;If keyword_set(read_gauss) then begin 
  ;   skip_gauss = 1 ;skip running the gaussian
  ;   print, 'Read gaussian file: Set'
  ;Endif  
  ;If not(keyword_set(read_gauss)) then print,'Save gaussian file (this takes a long time!)'  
  if n_elements(LIST) gt 2 then $
    print, 'LIST: set' else $
    print, 'LIST: not set'   
  print, 'Reduce Time (0=NO, 1=Yes)',Default_time
  print, 'QTmin, QTmax (secs): ', QTmin, QTmax  
  print,'Rebin Time 2 (secs)', rebin_time2
  print,'Rebin Freq 2 (MHz)', rebin_freq2
  If not(keyword_set(tau_Qplot)) then tau_Qplot = 2.0 
  if not(keyword_set(slope)) then slope = 2 
  If not(keyword_set(Qwin)) then Qwin = 10
  If not(keyword_set(root_gauss)) then root_gauss =''
  If not(keyword_set(gauss_steps)) then gauss_steps = 10000.
  If not(keyword_set(MINT)) then MINT = 1.0 
  If not(keyword_set(MAXT)) then MAXT = 5.9
  If not(keyword_set(OffBeam_Label)) then OffBeam_Label= 'Off Beam' 
    If not(keyword_set(OnBeam_Label)) then OnBeam_Label= 'On Beam' 
  print,'Tau for Q3 plot',tau_Qplot
  print,'Slope for Q3f/Q4f: ', slope
  print,'Tau for Q3 plot',tau_Qplot
  
  ps_dynspec = PS                      ;plot dynspec (1 = YES, 0 = NO)
  Target_Beam_l = 'beam'+strtrim(string(UINT(Target_Beam)),1)
  Sky_Beam_l    = 'beam'+strtrim(string(UINT(Sky_Beam)),1)

   ;---------------------
   ;Polarization Setup 
   ;----------------------
   If valid_num(S,/integer) eq 0 then begin
       pol_label    = strtrim(string(UINT(S)),1)
   endif
   If valid_num(S,/integer) eq 1 then begin
       If S eq 0 then begin    ;I
          pol_label    = 'I'
          pol_start  = 0.
          pol_end    = 0. 
       endif 
       If S eq 1 then begin    ;Q = Abs(Q)
          pol_label    = 'AbsQ'
          pol_start  = 0.
          pol_end    = 0. 
       endif 
       If S eq 2 then begin   ;U = Abs(U)
          pol_label    = 'AbsU'
          pol_start  = 0.
          pol_end    = 0. 
       endif 
       If S eq 3 then begin  ;V = Abs(V)
          ;loop labels
          pol_label = 'AbsV
          pol_start = 0 ;0 = Abs(V) 
          pol_end   = 0 
          
          If keyword_set(DynspecMS) then begin
            Vprime_label = strarr(3)
            Vprime_label[0] = 'AbsV'
            ;Vprime_label[1] = 'V2'
            Vprime_label[1] = 'V+'
            Vprime_label[2] = 'V-'
            ;loop labels
            pol_start = 0;0 = abs(V), 1=V+, 2 = V-
            pol_end   = 2  ;2 = V-
          Endif
       endif
        If S eq 4 then begin ; V= Vprime 
          Vprime_label = strarr(3)
          Vprime_label[0] = 'AbsV'
          ;Vprime_label[1] = 'V2'
          Vprime_label[1] = 'V+'
          Vprime_label[2] = 'V-'
          ;loop labels
          pol_start = 0;0 = abs(V), 1=V+, 2 = V- 
          pol_end   = 2  ;2 = V-
          
       endif 
       If S eq 5 then begin 
         pol_label      = 'L'
         pol_start = 0  ;0=L
         pol_end   = 0 
       endif 
       If S eq 6 then begin     ; this is for beams from imaging data!!!
        pol_label = 'V+'
        pol_start = 0
        pol_end  = 0  
       Endif
       If S eq 7 then begin     ; this is for beams from imaging data!!!
         pol_label = 'V-'
         pol_start = 0
         pol_end  = 0
       Endif
   endif ;end pol

  ;----------
  ;Q2 window
  ;----------
  Qwin = long(Qwin)       ; 10 steps; 10 second window if rebin_time = 1 secs
  If keyword_set(VERBOSE) then print, 'Window time (times x rebin time):', Qwin
  
  ;find Index for tau
  tau_index = where(threshold eq tau_Qplot)
  
  ;------------------------------
  ;Time Setup (not default time)
  ;------------------------------
  If not(keyword_set(Default_Time)) then begin
    Print,'----Choose Subset of Q2 time-------'
    ;QTmin = 0.079 ;Move Jupiter
    ;QTmax = 0.11  ;first half hour
    ;QTmax = 0.13   ;All of second day emission
    print,'Tmin Q2 (days):',QTmin
    print,'Tmax Q2 (days):',QTmax
  endif
  
  ;----------------------
  ;Ellipical Correction
  ;----------------------
  ;ElCorrect_or = ElCorrect
  ;If keyword_set(VERBOSE) then print,'Ellipical Correction Q2 (0 = NO, 1 = YES): ',ElCorrect
  
  If keyword_set(RUN_JUPITER) then begin
    input_file_J = 'Input_Jupiter.dat'  ;flag! (hardcoded)
    read_jupiter_input,input_file_J,scale=scale
    alpha = 1.0d/scale
    print,'Alpha: ', alpha
  Endif
  
  ;--------------------
  ;   Rebin time 
  ;--------------------
  Num_rebin_t2 = long(rebin_time2/rebin_time) ;Number of bins to combine for time
  
  If ELLCORR eq 0 then begin  ;run ell but no ELL
   STARTELL = 0 
   ENDELL   = 0
  Endif
  
  If ELLCORR eq 2 then begin  ;run both 
    STARTELL = 0
    ENDELL  = 1
  Endif
  
  If ELLCORR eq 3 then begin
    STARTELL = 1  ;start at 1 to run non-elliptical only  (i=0 is EL, i=1 is non el)
    ENDELL  = 1  
  Endif
 
  ;***************
  ;   Freq loop  
  ;***************
  for iii=0,n_elements(fmaxs)-1 do begin
    ;********************
    ;mask threshold loop
    ;********************
    for iij=0,n_elements(threshold_t)-1 do begin
     ;**************************
     ;Elliptical Correction Loop
     ;***************************
     for iil=STARTELL,ENDELL do begin    ;for loop for Elliptical Correction 
      ;***************************
      ; Polarization Outcome Loop 
      ;***************************
      ;loop for each polarization over its possible outcomes (0=Normal or Abs Value, 1=V-, 2=V+)
      for iih=pol_start,pol_end do begin
       If S eq 4 then pol_label =  Vprime_label[iih]     ;Vprime labels
       If S eq 3 and keyword_set(DYNSPECMS) then pol_label =  Vprime_label[iih]
        If iil eq 0 then begin
         print,'-------------------------------------------------'
         print,'------Running elliptical Correction--------------'
         print,'-------------------------------------------------'
         ElCorrect = 1  ;note 
        Endif
        If iil eq 1 then begin
         print,'-------------------------------------------------'
         Print,'----Running Regular (no Elliptical Correction)---'
         print,'-------------------------------------------------'
         ElCorrect = 0
        Endif
        count = 0.0d
        If keyword_set(VERBOSE) then print,'Freq: ',fmins[iii],fmaxs[iii],', Threshold: ', threshold_t[iij], ' Polarization:', S, ' Pol Type:', pol_label

        ;**********************
        ;restore Target Beam
        ;**********************
        If not(keyword_set(Run_Jupiter)) then restore,filename='rebindata_'+save_file_name+Target_Beam_l+'.sav'
        If keyword_set(Run_Jupiter) then restore,filename='rebindata_'+save_file_name+'_'+Target_Beam_l+'.sav'
          ;inputs: xt_rebin, MJD_rebin, xf_rebin, data_rebin, p2_rebin,UT_rebin
        If keyword_set(VERBOSE) then print,'Filename for Target Beam: ', 'rebindata_'+save_file_name+Target_Beam_l+'.sav'
         
        wf        = where(xf_rebin ge fmins[iii] and xf_rebin le fmaxs[iii])
        If not(keyword_set(Default_Time)) then wt = where(xt_rebin ge Qtmin and xt_rebin le Qtmax) else wt = where(xt_rebin eq xt_rebin)
        XTMIN = min(xt_rebin)/3600.0d ;hours
        XTMAX = max(xt_rebin)/3600.0d ;hours
        xf_rebin1   = xf_rebin(wf)
        x           = data_rebin(*,wf)
        px          = p2_rebin(*,wf)
        x           = x(wt,*)
        px          = px(wt,*) 
        xt_rebin1   = xt_rebin(wt)/3600.0d ;hours
        MJD_rebin   = MJD_rebin(wt)
        If keyword_set(LOFAR) then UT_Rebin     = UT_Rebin(wt) 
        phase_rebin  = phase_rebin(wt) 
        
        nf = n_elements(xf_Rebin1)

      ;----------------------------------
      ;Mask Threshold for individual  pixels
      ;---------------------------------- 
      count_bpix = 0 ;counter of bad pixels
      for i=0,n_elements(xf_rebin1)-1 do begin
        w_pf = where(px(*,i) lt threshold_p,count)
        If count gt 0 then begin 
         px(w_pf,i) = 0.0
         x(w_pf,i)  = 0.0
        endif 
        count_bpix = count_bpix + count
      endfor
      count = 0.0d ;reset count
      If keyword_set(VERBOSE) then print,'% Bad Individual Pixels (Target):', count_bpix/n_elements(px)   
      
      ;-----------------
      ;    Rebin time
      ;-----------------
      If keyword_set(VERBOSE) then print, 'Number of Time Bins to rebin',Num_rebin_t2
      reduce_array,x,[Num_rebin_t2,1],x,/dbl
      reduce_array,px*1.,[Num_rebin_t2,1],px,/dbl
      reduce_array,xt_rebin1,[Num_rebin_t2],xt_rebin1,/dbl
      If keyword_SET(LOFAR) then reduce_array,MJD_rebin1,[Num_rebin_t2],MJD_rebin1,/dbl
      If not(keyword_set(IMAGING)) and keyword_set(LOFAR) then reduce_array,UT_Rebin,[Num_rebin_t2],UT_rebin,/dbl
      nf = n_elements(xf_rebin1)
      nt = n_elements(xt_rebin1)
    
      ;--------------------
      ;  Normalization
      ;--------------------
      ;  x = x - mean(x)
      
      ;**********************
      ;restore Sky Beam
      ;**********************
      If not(keyword_set(Run_Jupiter)) then restore,filename='rebindata_'+save_file_name+Sky_Beam_l+'.sav'
      print,'Filename for Sky Beam: ', 'rebindata_'+save_file_name+Sky_Beam_l+'.sav'
      If keyword_set(Run_Jupiter) then begin
        l = STRPOS(save_file_name, 'Jupiter')
        save_file_name2 = STRMID(save_file_name, 0, l)
        print, 'New save filename for Sky: ', save_file_name2
        restore,filename='rebindata_'+save_file_name2+'_'+Sky_Beam_l+'.sav'
      Endif
        ;Input: xt_rebin, MJD_rebin, xf_rebin, data_rebin, p2_rebin, UT_rebin
      wf          = where(xf_rebin ge fmins[iii] and xf_rebin le fmaxs[iii])
      If not(keyword_set(Default_Time)) then wt = where(xt_rebin ge Qtmin and xt_rebin le Qtmax) else wt = where(xt_rebin eq xt_rebin)
      xf_rebin2     = xf_rebin(wf)
      x2            = data_rebin(*,wf)
      px2           = p2_rebin(*,wf)
      x2            = x2(wt,*)
      px2           = px2(wt,*)
      MJD_rebin2    = MJD_rebin(wt)
      If keyword_set(LOFAR) then UT_Rebin2     = UT_Rebin(wt) 
      phase_rebin2  = phase_rebin(wt) 
      xt_rebin2     = xt_rebin(wt)/3600.0d ;hours
      
      ;----------------------------------
      ;Mask Threshold for individual  pixels
      ;----------------------------------
      for i=0,n_elements(xf_rebin2)-1 do begin
        w_pf = where(px2(*,i) lt threshold_p)
        px2(w_pf,i) = 0.0
        x2(w_pf,i)  = 0.0
      endfor
    
      ;------------ 
      ;Rebin time
      ;------------
       reduce_array,x2,[Num_rebin_t2,1],x2,/dbl
       reduce_array,px2*1.,[Num_rebin_t2,1],px2,/dbl
       reduce_array,xt_rebin2,[Num_rebin_t2],xt_rebin2,/dbl

      ;--------------------
      ;  Normalization
      ;--------------------
      ;  x2 = x2 - mean(x2)      

      ;********************************************************************************
      ;                         Extra Beam 1: ExBeam1 Beam
      ;********************************************************************************
      If keyword_set(Extra_Beam1) then begin
        print,'Extra Beam 1'
        restore,filename='rebindata_'+save_file_name+$
          'beam'+strtrim(string(cgnumber_formatter(Extra_Beam1[0],decimals=0)),1)+'.sav'
            ;Input: xt_rebin, MJD_rebin, xf_rebin, data_rebin, p2_rebin
        wf = where(xf_rebin ge fmins[iii] and xf_rebin le fmaxs[iii])
        If not(keyword_set(Default_Time)) then wt = where(xt_rebin ge Qtmin and xt_rebin le Qtmax) else wt = where(xt_rebin eq xt_rebin)
        xf_rebin      = xf_rebin(wf)
        x3            = data_rebin(*,wf)
        px3           = p2_rebin(*,wf)
        x3            = x3(wt,*)
        px3           = px3(wt,*)
        MJD_rebin3    = MJD_rebin(wt)
        If keyword_set(LOFAR) then UT_Rebin3     = UT_Rebin(wt) 
        phase_rebin3  = phase_rebin(wt) 
        xt_rebin3     = xt_rebin(wt)/3600.0d ;hours
   
        ;----------------------------------
        ;Mask Threshold for individual  pixels
        ;----------------------------------
        for i=0,n_elements(xf_rebin)-1 do begin
          w_pf = where(px3(*,i) lt threshold_p)
          px3(w_pf,i) = 0.0
          x3(w_pf,i)  = 0.0
        endfor
        
        ;Rebin time
        reduce_array,x3,[Num_rebin_t2,1],x3,/dbl
        reduce_array,px3*1.,[Num_rebin_t2,1],px3,/dbl
        reduce_array,xt_rebin3,[Num_rebin_t2],xt_rebin3,/dbl
        reduce_array,MJD_rebin3,[Num_rebin_t2],MJD_rebin3,/dbl

        ;--------------------
        ;  Normalization
        ;--------------------
        ;x3 = x3 - mean(x3)
      Endif ;extra beam 1 input
      
      ;********************************************************************************
      ;                         Extra Beam 2: ExBeam2 Beam
      ;********************************************************************************
      If keyword_set(Extra_Beam2) then begin
        restore,filename='rebindata_'+save_file_name+$
          'beam'+strtrim(string(cgnumber_formatter(Extra_Beam2[0],decimals=0)),1)+'.sav'
           ;xt_rebin, MJD_rebin, xf_rebin, data_rebin, p2_rebin
        wf         = where(xf_rebin ge fmins[iii] and xf_rebin le fmaxs[iii])
        If not(keyword_set(Default_Time)) then wt = where(xt_rebin ge Qtmin and xt_rebin le Qtmax) else wt = where(xt_rebin eq xt_rebin)
        xf_rebin      = xf_rebin(wf)
        x4            = data_rebin(*,wf)
        px4           = p2_rebin(*,wf)
        x4            = x4(wt,*)
        px4           = px4(wt,*)
        xt_rebin4     = xt_rebin(wt)/3600.0d ;hours
        MJD_rebin4    = MJD_rebin(wt)
        If keyword_set(LOFAR) then UT_Rebin4     = UT_Rebin(wt) 
        phase_rebin4  = phase_rebin(wt) 
 
        ;Rebin time
        reduce_array,x4,[Num_rebin_t2,1],x4,/dbl
        reduce_array,px4*1.,[Num_rebin_t2,1],px4,/dbl
        reduce_array,xt_rebin4,[Num_rebin_t2],xt_rebin4,/dbl
        reduce_array,MJD_rebin4,[Num_rebin_t2],MJD_rebin4,/dbl

        ;--------------------
        ;  Normalization
        ;--------------------
        ;x4 = x4 - mean(x4)
      Endif ;extra beam 2 input 
     
        ;strings
        thres_t_string = 'thresT'+strtrim(string(cgnumber_formatter(threshold_t[iij],decimals=2)),1)
        thres_f_string = 'thresF'+strtrim(string(cgnumber_formatter(threshold_f[iij],decimals=2)),1)
        upper_string = strtrim(string(cgnumber_formatter(fmaxs[iii],decimals=1)),1)+'MHz'
        lower_string = strtrim(string(cgnumber_formatter(fmins[iii],decimals=1)),1)
        
        ;------------
        ;Setup Plot
        ;-----------
        filename_ps = save_file_name+lower_string+'to'+$
                     upper_string+'_'+thres_t_string+'_'+thres_f_string+'_postprocessing_'+pol_label
        If keyword_set(ElCorrect) then filename_ps = filename_ps+'_EllC'
        If not(keyword_set(Default_time)) then begin
          xtmin_label = strtrim(string(round(QTmin)))
          xtmax_label = strtrim(string(round(QTmax)))
          filename_ps =  filename_ps+'_Time'+xtmin_label+'-'+xtmax_label
          endif
        If keyword_set(VERBOSE) then print,'Filename Ps:', filename_ps
        ;device,filename=filename_ps+'.ps',/landscape
	       cgPS_Open,filename_ps+'.ps'
        
        ;Labels for plots
        freq_label   = strtrim(string(cgnumber_formatter(min(xf_rebin1),decimals=1)),0)+'-'+$
                     strtrim(string(cgnumber_formatter(max(xf_rebin1),decimals=1)),0)+'MHz'           
        thres_label  = strtrim(string(cgnumber_formatter(long((1. - threshold_t[iij])*100.))),0)+'%'
        label_main11 = pol_label+' '+freq_label+' '+thres_label  ;Header for Plots 
        If keyword_set(RUN_JUPITER) then begin
          Jup_label = strtrim(string(alpha,format='(E0.1)'),0)
          l = strlen(Jup_label)
          exp = strmid(Jup_label,l-3,l) 
          ;label_main11 = label_main11+' '+cgGreek('alpha')+':10!u'+exp+'!N'
          label_main11 = label_main11+' '+cgGreek('alpha')+':'+Jup_label
        Endif
        
        ;*********************
        ;   Time Threshold
        ;**********************
        If keyword_set(threshold_t) then begin
            ;*********************************************
            ; Combine Masks (in case not already combined)
            ;*********************************************
            REDUCE_ARRAY, px, [1,nf], px1_t
            REDUCE_ARRAY, px2, [1,nf], px2_t
            w12   = (px1_t ge 1. - threshold_t[iij])*(px2_t ge 1. - threshold_t[iij])
            w12t  =  where(w12 eq 1)  
           
            ;Apply to Target
            x           = x(w12t,*)
            px          = px(w12t,*)
            If keyword_Set(LOFAR) then MJD_rebin1  = MJD_rebin1(w12t)
            xt_rebin1   = xt_rebin1(w12t)
            nt1 = n_elements(xt_rebin1)
            
            ;Apply to Sky
            x2          = x2(w12t,*)
            px2         = px2(w12t,*)
            If keyword_Set(LOFAR) then MJD_rebin2  = MJD_rebin2(w12t)
            xt_rebin2   = xt_rebin2(w12t)
            nt2 = n_elements(xt_rebin2)             

            print,'Orginal Size of t:', n_elements(px1_t)
            print,'Size of t after threshold:', nt1 
            
            If nt1 le 1 then begin
              print,'*****************'
              print,'*****************'
              print,'Skipping Calc'
              print,'*****************'
              print,'*****************' 
              GOTO, skip
             endif 
           Endif ; end threshold_t 
          
          ;*******************
          ; Freq Threshold
          ;*******************
          If keyword_set(threshold_f) then begin
             REDUCE_ARRAY, px, [nt1,1], px1_f
             REDUCE_ARRAY, px2, [nt1,1], px2_f
             
             w12   = (px1_f ge 1. - threshold_f[iij])*(px2_f ge 1. - threshold_f[iij])
             ww12f  =  where(w12 eq 1)
             
             x        = x(*,ww12f)
             px       = px(*,ww12f)
             x2       = x2(*,ww12f)
             px2      = px2(*,ww12f)
             xf_rebin1 = xf_rebin1(ww12f)
             nf1 = n_elements(ww12f)
             nf2 = n_elements(ww12f)  
             print,'Orginal Size of f:', n_elements(px1_f)
             print,'Size of f after threshold:', n_elements(ww12f)
             
             If n_elements(ww12f) le 1 then begin
              print,'*****************'
              print,'*****************'
              print,'Skipping Calc'
              print,'*****************'
              print,'*****************'
              GOTO, skip
             endif 
        Endif
  
        If not(keyword_set(threshold_f)) then begin
          nf1 = n_elements(xf_rebin1)
          nf2 = n_elements(xf_rebin2)
        Endif
        
        ;-----------------------------------
        ; Polarization Setup
        ;-----------------------------------
        If S eq 4 or S eq 3 then begin 
           If iih eq 0 then begin  ; Abs(V)
            x  = abs(x)
            x2 = abs(x2)
           endif 
           If iih eq 1 then begin  ;V' > 0; V+
              wneg1   = where(x lt 0)   ;negative values
              wneg2   = where(x2 lt 0,count1)   ;negative values
              x(wneg1)   = 0.0   ;set all to zero V' < 0 
              px(wneg1)  = 0.0 
              x = abs(x) 
              x2(wneg2)   = 0.0   ;set all to zero V' < 0 
              px2(wneg2)  = 0.0 
              x2 = abs(x2) 
              count = count1 + count
           endif  ;iih eq 1
           xb =x
           If iih eq 2 then begin  ;V' < 0; V-
             wpos1     = where(x gt 0,count1)  ;positive values
             x(wpos1)  = 0.0   ;set all to zero V' > 0 
             px(wpos1) = 0.0
             x        = abs(x) 
             wpos2     = where(x2 gt 0,count1)  ;positive values
             x2(wpos2) = 0.0   ;set all to zero V' > 0 
             px2(wpos2)= 0.0
             x2       = abs(x2) 
             count    = count1 + count
           endif ; iih eq 2 
            
           If keyword_set(VERBOSE) then print, 'Polarization Count (All)',count, ' % after cut: ', count/n_elements(x)*100.0d      
         endif  ;S=4 or 3 
        
        ;---------------------------
        ;Thresholds; Extra Beam1 
        ;--------------------------- 
        If keyword_set(Extra_Beam1) then begin
          REDUCE_ARRAY, px3, [1,nf], px3_t
          REDUCE_ARRAY, px3, [nt,1], px3_f
          
          ;*******************
          ; Time Threshold
          ;*******************
           If keyword_set(threshold_t) then begin
             w12   = (px1_t ge 1. - threshold_t[iij])*(px3_t ge 1. - threshold_t[iij])*(px2_t ge 1. - threshold_t[iij])
             ww12  =  where(w12 eq 1)
             x3 = x3(ww12,*)
             px3 = px3(ww12,*)      
             If keyword_Set(LOFAR) then MJD_rebin3 = MJD_rebin(ww12)
             xt_rebin3 = xt_rebin3(ww12)
             nt3 = n_elements(xt_rebin3)
           Endif
          
          ;*******************
          ; Freq Threshold
          ;*******************
          If keyword_set(threshold_f) then begin
            w12   = (px1_f ge 1. - threshold_f[iij])*(px3_f ge 1. - threshold_f[iij])*(px2_f ge 1. - threshold_f[iij])
            w12  =  where(w12 eq 1)
            
            x3 = x3(*,w12)
            px3 = px3(*,w12)   
            nf3 = n_elements(xf_rebin1(w12))
          Endif
          
          If not(keyword_set(threshold_f)) then begin
            nf3 = n_elements(xf_rebin1)
          Endif
          
        Endif
        
        ;---------------------------
        ;Thresholds, Extra Beam 2
        ;---------------------------
        If keyword_set(Extra_Beam2) then begin 
          REDUCE_ARRAY, px4, [1,nf], px4_t
          REDUCE_ARRAY, px4, [nt,1], px4_f
          
          ;*******************
          ; Time Threshold
          ;*******************
           If keyword_set(threshold_t) then begin
             w12   = (px1_t ge 1. - threshold_t[iij])*(px4_t ge 1. - threshold_t[iij])*(px2_t ge 1. - threshold_t[iij])
             ww12  =  where(w12 eq 1)
             
             ;ww12 =  w12t  ;target
             x4 = x4(ww12,*)
             px4 = px4(ww12,*)
             If keyword_Set(LOFAR) then MJD_rebin4 = MJD_rebin(ww12)
             xt_rebin4 = xt_rebin4(ww12)
             nt4 = n_elements(xt_rebin3(ww12))
           Endif
          
         ;*******************
         ; Freq Threshold
         ;******************* 
          If keyword_set(threshold_f) then begin
            w12   = (px1_f ge 1. - threshold_f[iij])*(px4_f ge 1. - threshold_f[iij])*(px2_f ge 1. - threshold_f[iij])
            ww12  =  where(w12 eq 1)
            
           ; ww12 = ww12f 
            x4       = x4(*,ww12)
            px4      = px4(*,ww12)
            xf_rebin = xf_rebin(ww12)
            nf4      = n_elements(xf_rebin)
          Endif
          
          If not(keyword_set(threshold_f)) then begin
           nf4 = n_elements(xf_rebin)
          Endif
        Endif ;ExBeam2source 
       
       ;---------------
       ;Plot DynSpecs
       ;---------------
       If keyword_set(ps_dynspec) then begin
          If not(keyword_set(PAPERPLOTS)) then label   = label_main11 Else label ='Dynamic Spectrum'
          If not(keyword_set(UT)) then begin 
           STANDARD_PLOTS,x,px,xt_rebin1,xf_rebin1,label,/ONLY_PLOT,/NO_ZOOM, xunit = 'Hour',ytitle_bar='Intensity (SEFD)'
           If not(keyword_set(PAPERPLOTS)) then label   = label_main11 Else label =''
           STANDARD_PLOTS,x2,px,xt_rebin1,xf_rebin1,label,/ONLY_PLOT,/NO_ZOOM, xunit = 'Hour',ytitle_bar='Intensity (SEFD)'
          endif
          If keyword_set(UT) then begin
             STANDARD_PLOTS,x,px,UT_rebin,xf_rebin1,'Dynamic Spectrum',/ONLY_PLOT,/NO_ZOOM,skyplot=x2, xlabel='UT',xunit = 'hour',ytitle_bar='Intensity (SEFD)',/panel_ps
             !p.multi=[0,1,1]
             STANDARD_PLOTS,x2,px,UT_rebin,xf_rebin1,'',/ONLY_PLOT,/NO_ZOOM,skyplot=x2, xlabel='UT',xunit = 'hour',ytitle_bar='Intensity (SEFD)',/panel_ps
          endif 
       Endif
       
  ;***********************************************************
  ;           Q1: total power of each TI over all frequencies
  ;***********************************************************
      If keyword_set(Q1) then begin
        Num_rebin_t = long(TI/rebin_time2)         ;Number of bins to combine for time interval
        Num_rebin_f = long(rebin_freq2/(rebin_freq*1d-3)) ;Number of bins to combine for freq interval (rebin_freq in kHz and rebin_freq2 in MHz)
        ;Guassian Errors 
        Error = 1./sqrt(TI*(fmaxs[iii]-fmins[iii])*1e6)  ;error on flux; 1/sqrt(b*t)
        Error_diff = sqrt(2)*Error

      ; If not(keyword_set(LIST)) then begin
        ;print,'******NO LIST******'
        ;--------------------------------------------------------
        ;                      Target
        ;--------------------------------------------------------
         If keyword_Set(LOFAR) then REDUCE_ARRAY, MJD_rebin1,[Num_rebin_t],MJD_rebin1_q1,/dbl  ;rebin time
         REDUCE_ARRAY, xt_rebin1,[Num_rebin_t],xt_rebin1_q1,/dbl    ;rebin time 
         REDUCE_ARRAY, x, [Num_rebin_t,nf1], x_q1,/dbl              ;rebin data for time series
         REDUCE_ARRAY, px*1., [Num_rebin_t,nf1], px_q1,/dbl         ;rebin mask
         If keyword_Set(LOFAR) then begin 
           If not(keyword_set(IMAGING)) then REDUCE_ARRAY, UT_rebin,[Num_rebin_t],UT_rebin_q1,/dbl    ;rebin time 
         endif 
         x_q1 = x_q1/(px_q1 + (px_q1 eq 0))                         ;apply mask ElCorrectly
       
         ;Freq 
         REDUCE_ARRAY, x, [nt1,Num_rebin_f], xf_q1,/dbl              ;rebin data for Intensity (SEFD)
         REDUCE_ARRAY, px*1., [nt1,Num_rebin_f], pxf_q1,/dbl         ;rebin mask
         REDUCE_ARRAY,xf_rebin1,[Num_rebin_f],xf_rebin_q1,/dbl        ;rebin freq
         xf_q1 = xf_q1/(pxf_q1 + (pxf_q1 eq 0))
        
        ;--------------------------------------------------------
        ;                      Sky
        ;--------------------------------------------------------
         If keyword_Set(LOFAR) then REDUCE_ARRAY, MJD_rebin2,[Num_rebin_t],MJD_rebin2_q1,/dbl
         REDUCE_ARRAY, xt_rebin2,[Num_rebin_t],xt_rebin2_q1,/dbl
         REDUCE_ARRAY, x2, [Num_rebin_t,nf2], x2_q1,/dbl
         REDUCE_ARRAY, px2*1., [Num_rebin_t,nf2], px2_q1,/dbl
         x2_q1 = x2_q1/(px2_q1 + (px2_q1 eq 0))
         
         ;Freq
         REDUCE_ARRAY, x2, [nt1,Num_rebin_f], x2f_q1,/dbl              ;rebin data for time series
         REDUCE_ARRAY, px2*1., [nt1,Num_rebin_f], px2f_q1,/dbl         ;rebin mask
         x2f_q1 = x2f_q1/(px2f_q1 + (px2f_q1 eq 0))
        
         ;Save Q1
         Q1 = dblarr(2,n_elements(xt_rebin2_q1))
         Q1[0,*] = x_q1[*]
         Q1[1,*] = x2_q1[*]
         
         ;--------------------------------------------------------
         ;                      Extra Beam 1
         ;--------------------------------------------------------
         If keyword_set(Extra_Beam1) then begin
           If keyword_Set(LOFAR) then REDUCE_ARRAY, MJD_rebin3,[Num_rebin_t],MJD_rebin3_q1,/dbl
           REDUCE_ARRAY, xt_rebin3,[Num_rebin_t],xt_rebin3_q1,/dbl
           REDUCE_ARRAY, x3, [Num_rebin_t,nf3], x3_q1,/dbl
           REDUCE_ARRAY, px3*1., [Num_rebin_t,nf3], px3_q1,/dbl
           x3_q1 = x3_q1/(px3_q1 + (px3_q1 eq 0))
           ExBeam1_info      = dblarr(3,n_elements(x3_q1))
           If keyword_Set(LOFAR) then ExBeam1_info[0,*] = MJD_rebin3_q1
           ExBeam1_info[1,*] = x3_q1
           ExBeam1_info[2,*] = xt_rebin3_q1
         Endif
  
         ;--------------------------------------------------------
         ;                      Extra Beam 2
         ;--------------------------------------------------------
         If keyword_set(Extra_Beam2) then begin
           If keyword_Set(LOFAR) then REDUCE_ARRAY, MJD_rebin4,[Num_rebin_t],MJD_rebin4_q1,/dbl
           REDUCE_ARRAY, xt_rebin4,[Num_rebin_t],xt_rebin4_q1,/dbl
           REDUCE_ARRAY, x4, [Num_rebin_t,nf4], x4_q1,/dbl
           REDUCE_ARRAY, px4*1., [Num_rebin_t,nf4], px4_q1,/dbl
           x4_q1 = x4_q1/(px4_q1 + (px4_q1 eq 0))
           ExBeam2_info = dblarr(3,n_elements(x4_q1))
           If keyword_Set(LOFAR) then ExBeam2_info[0,*] = MJD_rebin4_q1
           ExBeam2_info[1,*] = x4_q1
           ExBeam2_info[2,*] = xt_rebin4_q1
         Endif  
      ;endif ; not list 
     Endif ;end Q1
      
      ;--------------------------------------------
      ;Derivative of the Beams (NOT used anymore)
      ;--------------------------------------------
      If keyword_set(plot_derv) then begin
       If not(keyword_set(TARGET_ONLY)) then begin
         ;-------------------------------------
         ;Beam Comparison: Derivative   
         ;-------------------------------------
         dx_q1    = deriv(xt_rebin1_q1,x_q1)   
         dx2_q1   = deriv(xt_rebin1_q1,x2_q1)
         error_dx = derivsig(xt_rebin1_q1,xt_rebin1_q1,0,error)
         
         If not(keyword_set(PAPERPLOTS)) then label='Q1: Derivative '+label_main11 else label='Q1 Derivative' 
         ytitle = 'Derivative'
          diagnostic_plots,dx_q1,dx2_q1,MJD_rebin1_q1,xt_rebin1_q1,P,T0,label,ytitle,t2=MJD_rebin1_q1,$
           LIST=LIST_REBIN,/single_plot,error =error_dx,TARGET_ONLY=TARGET_ONLY,diff=(dx_q1-dx2_q1);,UT=UT_q1
       endif ;if not target only 
      endif 

      ;-----------------Q1 plots ------------------------
      ;Plot Target, sky, and difference
      ;--------------------------------------------------
      If not(keyword_set(TARGET_ONLY)) then begin
         xd_q1 = x_q1 - x2_q1   ;Q1 Diff
        
         ;---------------------------
         ;Values for Output text file 
         ;---------------------------
         calc_ppextrav2,xt_rebin1_q1,x_q1,/RunQ1,/Q1a,y2=x2_q1,out=Q1aV,gauss= error_diff,file_g=file_g    
         
         print,'************************************'
         print,'Prob of Q1a being false detection:', Q1av[3]         
         print,'************************************'

         ;------------------------------
         ;Plot of Q1 with diff included 
         ;;-----------------------------          
         If not(keyword_set(PAPERPLOTS)) then label='Q1a: '+label_main11 else label = 'Q1a'
         ytitle = 'Intensity (SEFD)'
         diagnostic_plots,x_q1,x2_q1,MJD_rebin1_q1,xt_rebin1_q1,P,T0,label,ytitle,t2=MJD_rebin2_q1,$
           TARGET_ONLY=TARGET_ONLY,LIST=LIST_REBIN,$
           ExBeam1=ExBeam1_info,ExBeam2=ExBeam2_info,$
           /single_plot,error=error,ONBEAM_LABEL=ONBEAM_LABEL,OFFBEAM_LABEL=OFFBEAM_LABEL,max_scale = 1.19;diff=xd_q1
          
          ;----------------
          ;Only diff plot 
          ;----------------
          If not(keyword_set(PAPERPLOTS)) then label='Q1a (ON-OFF): '+label_main11 else label = 'Q1a (ON-OFF)'       
	        diagnostic_plots,xd_q1,x2_q1,MJD_rebin1_q1,xt_rebin1_q1,P,T0,label,'Intensity (SEFD)',t2=MJD_rebin2_q1,$
             /TARGET_ONLY,LIST=LIST_REBIN,$
             ExBeam1=0,ExBeam2=0,line=1,$
             /single_plot,error=error_diff;,UT=UT_q1
           
           If keyword_set(VERBOSE) then print,'Q1: Target STDEV:', stddev(x_q1)
           If keyword_set(VERBOSE) then print, 'Q1: Sky STDEV:', stddev(x2_q1)
           If keyword_set(VERBOSE) then print,'Q1: Difference:',stddev(xd_q1) 
      endif        

      ;*****************************
      ;     Integrated Spectrum
      ;*****************************
      If not(keyword_set(TARGET_ONLY)) then begin
         ;error in integrated spectrum
         Error_i = 1./sqrt(rebin_freq2*1e6*(max(xt_rebin1)-min(xt_rebin1))*86400.d)
         error_i_diff  = sqrt(2)*Error_i
         
         ;Difference
         xfq1_diff = xf_q1 - x2f_q1       
    
         ;---------------------
         ;Values for text file 
         ;---------------------
         calc_ppextrav2,xf_rebin_q1,xf_q1,/RunQ1,/Q1b,y2=x2f_q1,out=Q1bV,gauss= error_i_diff,file_g=file_g 

         print,'************************************'
         print,'Prob of Q1b being false detection:', Q1bV[3]
         print,'************************************'

        If not(keyword_set(PAPERPLOTS)) then label='Q1b: '+label_main11 else label='Q1b' 
        ytitle = 'Intensity (SEFD)'
        cgplot,xf_rebin_q1,xf_q1,err_yhigh=error_i,err_ylow=error_i,title=label,ytitle=ytitle,$
          xtitle='Frequency (MHz)',color=cgcolor('Black'),psym=46,$
          yrange=[(min([xf_q1,x2f_q1])-error[0]*1.2),(max([xf_q1,x2f_q1])+error[0])*1.2],/xstyle,/ystyle,charsize=char_size
        cgplot,xf_rebin_q1,x2f_q1,err_yhigh=error_i,err_ylow=error_i,/overplot,color=cgcolor('Red'),psym=46
        cgplot,[-1000,1000],[0,0],LineStyle=0,/overplot
        al_legend,[ONBEAM_LABEL,OFFBEAM_LABEL],$
         colors=['Black','red'],psym=[46,46],/left,/top,charsize=2.0,symsize=2      

        ;-----------
        ;diff plot
        ;-----------
        maxs = max(xfq1_diff) + error_i_diff[0]*1.02
        mins = min(xfq1_diff) - error_i_diff[0]*1.02
        If not(keyword_set(PAPERPLOTS)) then label='Q1b (ON-OFF): '+label_main11 else label='Q1b (ON-OFF)' 
        !p.multi=0
        cgplot,xf_rebin_q1,xfq1_diff,title=label,ytitle='Intensity (SEFD)',xtitle='Frequency (MHz)',color=cgcolor('Black'),psym=-46,$
          /xstyle,/ystyle,err_yhigh=error_i_diff,err_ylow=error_i_diff,yrange=[mins,maxs],charsize=char_size
        cgplot,[-1000,1000],[0,0],LineStyle=0,/overplot
      endif
    
     save,xf_rebin_q1,xf_q1,xfq1_diff,x2f_q1,x_q1,x2_q1,xd_q1,error_i_diff,error_i,error,error_diff,MJD_rebin,phase_rebin,filename=filename_ps+'FREQPLOT.sav'

    ;***********************************************************
    ;           Q2: Normalized time series (high-pass filtered and divided by standard deviation)
    ;***********************************************************
    If keyword_set(Q2) then begin
      ;Create high pass filter 
      print, 'Window Size (steps of rebin_time2)',Qwin
      ;If not(keyword_set(LIST)) then begin
        
          ;***********
          ;   Target
          ;***********
          REDUCE_ARRAY, x, [1,nf1], x_q2,/dbl
          REDUCE_ARRAY, px*1., [1,nf1], px_q2,/dbl
          x_q2    = x_q2/(px_q2 + (px_q2 eq 0))
          x_q1_HP = smooth(x_q2,Qwin,/EDGE_TRUNCATE)   
          x_q2    = x_q2 - x_q1_HP
          stddev_x1 = stddev(x_q2)
          y       = (x_q2 - mean(x_q2))/stddev(x_q2) 
          print,'Sigma of Q2 (Target): ',stddev(x_q2) 
          
          ;dynspec
          y_dyn = x - smooth(x,[Qwin,1])   
          print, 'High-pass Dyn Spec Target Sigma',stddev(y_dyn)
          y_dyn = (y_dyn - mean(y_dyn))/stddev(y_dyn)
          label = 'High-Pass Filter (Target): '+label_main11
          
          ;***********
          ;   Sky
          ;***********
          REDUCE_ARRAY, x2, [1,nf1], x2_q2,/dbl
          REDUCE_ARRAY, px2*1., [1,nf1], px2_q2,/dbl
          x2_q2    = x2_q2/(px2_q2 + (px2_q2 eq 0))
          x2_q2_HP = smooth(x2_q2,Qwin,/EDGE_TRUNCATE) 
          x2_q2    = x2_q2 - x2_q2_HP 
          y2       = (x2_q2 - mean(x2_q2))/stddev(x2_q2) 
          print, 'Sigma for Q2 (Sky): ', stddev(x2_q2) 
          
          ;dynspec 
          y2_dyn = (x2 - smooth(x2,[Qwin,1]))
          stddev_y2dyn = stddev(y2_dyn)
          print, 'High-pass Dyn Spec Sky Sigma',stddev(y2_dyn)
          y2_dyn = (y2_dyn - mean(y2_dyn))/stddev(y2_dyn)
          label = 'High-Pass Filter: '+label_main11
          
          Q2 = dblarr(2,n_elements(y))
          Q2[0,*] = y[*]
          Q2[1,*] = y2[*]
          
          ;----------------
          ;Max Threshold
          ;--------------
          If keyword_set(max_threshold) then begin 
            w  = where(y lt max_threshold and y gt -1.0*max_threshold and y2 lt max_threshold and y2 gt -1.0*max_threshold)
            y = y[w]
            y2 = y2[w]
          endif 
          
          ;-----------------
          ;Cut Time
          ;-----------------
          ;If not(keyword_set(Default_Time)) then begin
          ;  ww = where(xt_rebin1 gt QTmin and xt_rebin1 le QTmax) 
          ;  xt_rebin1  = xt_rebin1[ww] 
          ;  MJD_rebin1 = MJD_rebin1[ww] 
           ; y          = y[ww]
          ;  y2         = y2[ww]
          ;  print, 'Number of Elements: ', n_elements(y)
          ;endif
          
          If not(keyword_set(PAPERPLOTS)) then label = 'Q2: '+label_main11 else label = 'Q2'
          ytitle = 'Intensity (Sigma)'
          
          If keyword_set(Extra_Beam1) then begin
            REDUCE_ARRAY, x3, [1,nf3], x3_q2,/dbl
            REDUCE_ARRAY, px3*1., [1,nf3], px3_q2,/dbl
            x3_q2 = x3_q2/(px3_q2 + (px3_q2 eq 0))
            x3_q2_HP = smooth(x3_q2,Qwin,/EDGE_TRUNCATE) 
            x3_q2 = x3_q2 - x3_q2_HP
            y3 = (x3_q2 - mean(x3_q2))/stddev(x3_q2)
            
            ExBeam1_info = dblarr(3,n_elements(y3))
            ExBeam1_info[0,*] = MJD_rebin3
            ExBeam1_info[1,*] = y3
            ExBeam1_info[2,*] = xt_rebin3
          endif
  
          If keyword_set(Extra_Beam2) then begin
            REDUCE_ARRAY, x4, [1,nf4], x4_q2,/dbl
            REDUCE_ARRAY, px4*1., [1,nf4], px4_q2,/dbl
            x4_q2 = x4_q2/(px4_q2 + (px4_q2 eq 0))
            x4_q2_HP = smooth(x4_q2,Qwin,/EDGE_TRUNCATE) 
            x4_q2 = x4_q2 - x4_q2_HP
            y4 = (x4_q2 - mean(x4_q2))/stddev(x4_q2)
            
            ExBeam2_info = dblarr(3,n_elements(y))
            ExBeam2_info[0,*] = MJD_rebin4
            ExBeam2_info[1,*] = y4
            ExBeam2_info[2,*] = xt_rebin4
          endif  
      ;endif ;not LIST
      
      ;--------------
      ;Density plot
      ;--------------
      If keyword_Set(density) then begin
        xrange = [Min(y), Max(y)]
        yrange = [Min(y2), Max(y2)]
        xbinsize = 0.10
        ybinsize = 0.10
        cgLoadCT, 33
        TVLCT, cgColor('gray', /Triple), 0
        TVLCT, r, g, b, /Get
        palette = [ [r], [g], [b] ]
        density = Hist_2D(y, y2,Bin1=xbinsize,Bin2=ybinsize)
        maxDensity = Ceil(Max(density))
        scaledDensity = BytScl(density, Min=0, Max=maxDensity)
        cgImage, scaledDensity, XRange=xrange, YRange=yrange, /Axes, Palette=palette, $
        XTitle='Target', YTitle='Sky', $
        Position=[0.15, 0.25, 0.95, 0.90],title='Q2: '+label_main11
        cgplot,[0,0],[-1e6,1e6],/overplot
        cgplot,[-1e6,1e6],[0,0],/overplot
        thick = (!D.Name EQ 'PS') ? 6 : 2
        cgColorbar, Position=[0.15, 0.05, 0.95, 0.10], Title='Density', $
        Range=[0, maxDensity], NColors=254, Bottom=1, OOB_Low='gray', $
        TLocation='Top'
        cgLoadCT, 0
        TVLCT, cgColor('gray', /Triple), 0
       endif ;end density
      endif ;end Q2
     
     ;------------------------------------------- 
     ;ElCorrected; Sym Scatter symmterical 
     ;-------------------------------------------
      ;renormalize the sigmas (Jupiter paper explains) is now in the elliptical_correction procedure
     If keyword_set(ElCorrect) then begin
       ynew = y 
       y2new = y2
       elliptical_correction,ynew,y2new,xout=yel,yout=y2el
     endif 
        
      ;---------------------
      ;2D scatter plot
      ;--------------------
      If keyword_set(ElCorrect) then range_max = max( [abs(y),abs(y2),abs(yel),abs(y2el)])+0.5  else range_max = max( [abs(y),abs(y2)])+0.5
      print, 'Range Max: ', range_max
      
      If not(keyword_set(PAPERPLOTS)) then label = 'Q2: '+label_main11 else label = 'Q2'
      cgplot,y,y2,title=label,xtitle=OnBeam_Label+' (sigma)',ytitle=OffBeam_Label+' (sigma)',$
            yrange=[-range_max,range_max],xrange=[-range_max,range_max],psym='star',$
            charsize=char_size,/ISOTROPIC,XGridStyle=1, YGridStyle=1;,aspect=1;,XTicklen=1.0, YTicklen=1.0
      cgplot,[0,0],[-1e6,1e6],/overplot
      cgplot,[-1e6,1e6],[0,0],/overplot
      for i=-(round(range_max)),round(range_max) do begin
	    cgplot,[-100,100],[i,i],/overplot,thick=0.5,linestyle=1
	    cgplot,[i,i],[-100,100],/overplot,thick=0.5,linestyle=1
      endfor

      If not(keyword_set(ElCorrect)) then begin 
       cgplot,y,y2,title=label,xtitle=OnBeam_Label+' (sigma)',ytitle=OffBeam_Label+' (sigma)',$
            yrange=[-range_max,range_max],xrange=[-range_max,range_max],psym='star',$
            charsize=char_size,/ISOTROPIC,XGridStyle=1, YGridStyle=1;,aspect=1;,XTicklen=1.0, YTicklen=1.0
       cgplot,[0,0],[-1e6,1e6],/overplot
       cgplot,[-1e6,1e6],[0,0],/overplot
       for i=-(round(range_max)),round(range_max) do begin
	     cgplot,[-100,100],[i,i],/overplot,thick=0.5,linestyle=1
	     cgplot,[i,i],[-100,100],/overplot,thick=0.5,linestyle=1
       endfor
      endif 
 
      ;--------- EL Correction -----------
      If keyword_set(ElCorrect) then begin 
       ;2D scatter plot
       If not(keyword_set(PAPERPLOTS)) then label='Q2 (ElCorrected): '+label_main11 else label = 'Q2 (Elliptical Correction)'

       ;update 
       y = yel
       y2 = y2el           
       cgplot,y,y2,title=label,xtitle=OnBeam_Label+' (sigma)',ytitle=OffBeam_Label+' (sigma)',$
         yrange=[-range_max,range_max],xrange=[-range_max,range_max],psym='star',charsize=char_size,/ISOTROPIC,$
         XGridStyle=1, YGridStyle=1,aspect=1
       cgplot,[0,0],[-100,100],/overplot
       cgplot,[-100,100],[0,0],/overplot 
        for i=-(round(range_max)),round(range_max) do begin
	       cgplot,[-100,100],[i,i],/overplot,thick=0.5,linestyle=1
	       cgplot,[i,i],[-100,100],/overplot,thick=0.5,linestyle=1
        endfor
     endif ; El Correction

      ;---------------------
      ;dynspec plot; normal
      ;---------------------
       label = 'High-Pass Filter: '+label_main11
      If not(keyword_set(UT)) then begin ;Flagg
        STANDARD_PLOTS,y_dyn,px,xt_rebin1,xf_rebin1,label,mask=mask,/VERBOSE,/ONLY_PLOT,$
	     /NO_ZOOM,xunit = 'Hour',ytitle_bar='Intensity (Sigma)'
        STANDARD_PLOTS,y2_dyn,px,xt_rebin1,xf_rebin1,'Sky: '+label,mask=mask,/VERBOSE,/ONLY_PLOT,$
	     /NO_ZOOM,xunit = 'Hour',ytitle_bar='Intensity (Sigma)'
       endif 
       If keyword_set(UT) then begin  ;Flagg
         STANDARD_PLOTS,y_dyn,px,UT_rebin,xf_rebin1,label,mask=mask,/VERBOSE,/ONLY_PLOT,skyplot=y2_dyn,$
            /NO_ZOOM,xlabel='UT', xunit = 'hour',ytitle_bar='Intensity (Sigma)'
       endif 

       If keyword_set(TARGET_ONLY) then begin 
         label='Q1a: '+label_main11
         ytitle = 'Intensity (SEFD)'
         diagnostic_plots,x_q1,x2_q1,MJD_rebin1_q1,xt_rebin1_q1,P,T0,label,ytitle,t2=MJD_rebin2_q1,$
            TARGET_ONLY=TARGET_ONLY,LIST=LIST_REBIN,$
            ExBeam1=ExBeam1_info,ExBeam2=ExBeam2_info,$
            /single_plot,error=error,UT=UT_q1
       endif    
          
          ;---------------------------------------------------------------------
          ;dyspec plot w/ high sigma highlighted (NOT Working; work in progress)
          ;---------------------------------------------------------------------
           y_dynr = y_dyn 
           y2_dynr = y2_dyn
           new_bandwidth = (xf_rebin_q1[1] - xf_rebin_q1[0])  ;MHz
           bandwidth     = max(xf_rebin_q1) - min(xf_rebin_q1) ;MHz
           new_sigma = tau_Qplot*1./(sqrt(bandwidth/new_bandwidth))*(stddev_x1/stddev_y2dyn)^(-1.0)   ;sqrt(old/new) - convert between sigmas, (sdev/stdev2) -- units are different 
           wr = where(y_dynr ge new_sigma and y_dynr ge 2.0*y2_dynr, count_y)
           w2r = where(y2_dynr ge new_sigma and y2_dynr gt 2.0*y_dynr,count_y2)
           print, 'Number of pixels above new tau in dynspec (Target)',count_y
           print, 'Number of pixels above new tau in dynspec (Sky)', count_y2
           nt_r = n_elements(y_dynr[*,0])
           nf_r = n_elements(y_dynr[0,*])
           print,'Sigma in Dynamic plot to highlight:', new_sigma 
           b      = fltarr(nt_r,nf_r) ;filler arrays
           b2     = fltarr(nt_r,nf_r) ;filler
           If count_y gt 0 then b[wr]   = 1.
           If count_y2 gt 0 then b2[w2r] = 1. 
           b      = dilate(b,fltarr(3,3)+1) gt 0
           b2     = dilate(b2,fltarr(3,3)+1) gt 0     
           ;wr     = where(b eq 1)
           ;wr2    = where(b2 eq 1)
           y_dynr  = b*30.+randomn(seed,nt_r,nf_r)    ;final array
           y2_dynr = b2*30.+randomn(seed,nt_r,nf_r)   ;final array 

      ;-----------------------------------------------------------
      ;			 Q2 vs time
      ;-----------------------------------------------------------
      If not(keyword_set(PAPERPLOTS)) then label = 'Q2: '+label_main11 else label='Q2'
      ytitle = 'Intensity (Sigma)'
	    ;multiplot
      diagnostic_plots,y,y2,MJD_rebin1,xt_rebin1,P,T0,label,ytitle,t2=MJD_rebin2,/single_plot,ONBEAM_LABEL=ONBEAM_LABEL,OFFBEAM_LABEL=OFFBEAM_LABEL,max_scale=1.20
      
      If keyword_set(UT) then begin
        diagnostic_plots,y,y2,MJD_rebin1,xt_rebin1,P,T0,label,ytitle,t2=MJD_rebin,$
         /single_plot,ONBEAM_LABEL=ONBEAM_LABEL,OFFBEAM_LABEL=OFFBEAM_LABEL,max_scale=1.20,UT=UT_rebin
      endif 

      If keyword_set(Extra_Beam1) or keyword_set(Extra_Beam2) then begin
        diagnostic_plots,y,y2,MJD_rebin1,xt_rebin1,P,T0,label,ytitle,ExBeam1=ExBeam1_info,t2=MJD_rebin2,$
          ExBeam2=ExBeam2_info
      endif
       
      y_diff = y - y2  
      
      ;Difference Plot
      label = 'Q2 (ON-OFF): '+label_main11
      ytitle = 'Intensity (Sigma)'	
      diagnostic_plots,y_diff,y_diff,MJD_rebin1,xt_rebin1,P,T0,label,ytitle,t2=MJD_rebin2,$
        /TARGET_ONLY,/single_plot;,UT=UT_rebin
         
      axis_scale  = 1.05
      
      ;---------------------------
      ;    Gaussian Comparison Function
      ;---------------------------
       n_thres = n_elements(threshold)
       theshold_sav = threshold ;save threshold
        ;--------------------------
        ;Restore Gaussian Q values
        ;--------------------------
         Ny = n_elements(y)
         print, 'Number of Elements used for Gaussian: ',ny
         load_gaussianv2,theshold_sav,Ny,Q=q4_gauss,file=File_G
         
       ;---------------------
       ;Q2 Histogram 
       ;---------------------
      If keyword_set(Q2_histogram) then begin 
       yhist = cgHistogram(y,nbins = 50,BINSIZE=0.5)
       yhist2 = cgHistogram(y2,nbins = 50,BINSIZE=0.5)
       cghistoplot,y,xrange= [-10,10],yrange= [0.,max([yhist,yhist2])],nbins = 50,BINSIZE=0.5,Title='Histogram of Q2',THICK=2,color='black'
       cghistoplot,y2,xrange= [-10,10],yrange= [0.,max([yhist,yhist2])],nbins = 50,BINSIZE=0.5,/oplot,color='Red'
       al_legend,['On Beam','Off Beam'],textcolor=['Black','Red'],/right,/top,charsize=2.0,symsize=2
      
       cghistoplot,y,xrange= [-10,10],yrange= [0.,max(yhist)/10],nbins = 50,BINSIZE=0.5,Title='Histogram of Q2',THICK=2,color='black'
       cghistoplot,y2,xrange= [-10,10],yrange= [0.,max(yhist)/10],nbins = 50,BINSIZE=0.5,/oplot,color='Red'
       al_legend,['On Beam','Off Beam'],textcolor=['Black','Red'],/right,/top,charsize=2.0,symsize=2
      endif 

    ;------------
    ;create q4
    ;------------
    q4  = dblarr(2,n_thres,6)  ;ON and OFF(2), nthres, (3a,b,c,d,e,f)
  
    ;---------------------------
    ;Loop over threshold values 
    ;---------------------------  
    for iik=0,n_thres-1 do begin      
      N_label     = strtrim(string(cgnumber_formatter(tau_Qplot,decimal=1)),0)
      label_main12 = freq_label+' '+thres_label+' '+cgGreek('tau')+':'+N_label
      If keyword_set(RUN_JUPITER) then begin
        label_main12 = label_main12+' '+cgGreek('alpha')+':10!u'+exp+'!N'
      Endif
             
      ;---------------------------------    
      ;---------------------------------
      ; Q4: Sum over all time 
      ;---------------------------------
      ;----------------------------------
       If keyword_set(Q4_run) then begin
         ;--------------------
         ;Find Peaks and Power
         ;--------------------
         postprocessing_calc,y,y2,threshold[iik],slope=slope,$
                             NxA=NxA,PNxA=Power_NxA,$
                             Nx2A=Nx2A,PNx2A=Power_Nx2A,$
                             Nx3Ae=NxA_3e,PNAx3f=Power_NxA_3f,$
                             Nx3Be=NxB_3e,PNBx3f=Power_NxB_3f,$
                             NxB=NxB,PNxB=Power_NxB,$
                             Nx2B=Nx2B,PN2xB=Power_Nx2B
         Qc_x    = 0.  & Q3c_x2 = 0. 
         Q3d_x   = 0.  & Q3d_x2  = 0.     
         Q3c_x  = NXA  - NxB               ;[(N > tau) - (N < - tau)] - On
         Q3c_x2 = Nx2A - Nx2B              ;[(N > tau) - (N < - tau)] - Off  
         Q3d_x  = Power_NXA  - Power_NxB   ;[Power(N > tau) - Power(N < - tau)] - On
         Q3d_x2 = Power_Nx2A - Power_Nx2B  ;[Power(N > tau) - Power(N < - tau)] - Off
         
        ;-----------------------------------------------------------
        ;                       gaussian
        ;-------------------------------------------------------------
;        If not(keyword_set(read_gauss)) then begin
;          ;print,'Start Gaussian'
;          gauss_lofar_postprocessing,gauss_steps,n_elements(y),slope=slope,threshold[iik],Q3a=Q3a_g,dQ3a=Q3a_diff_g,Q3b=Q3b_g,dQ3b=Q3b_diff_g,$
;            Q3c=Q3c_g,dQ3c=Q3c_diff_g,Q3d=Q3d_g,dQ3d=Q3d_diff_g,Q3e=Q3e_g,dQ3e=Q3e_diff_g,Q3f=Q3f_g,dQ3f=Q3f_diff_g,square=square
;           If iil eq 0 then skip_gauss = 1 ;skip next round
;        endif 
         
      
        ;------------------------------------------------------------------------
        ;                 Save Q4 
        ;------------------------------------------------------------------------
          q4[0,iik,0]       = NxA    ;sum of 3a for ON BEAM  
          q4[1,iik,0]       = Nx2A   ;sum of 3a for OFF BEAM
;          If not(keyword_set(read_gauss)) then  q4_gauss[0,iik,0] = Q3a_g          ;total(NA_g)   ;sum of 3a for Guass 
;          If not(keyword_set(read_gauss)) then q4_gauss[1,iik,0] = Q3a_diff_g           ;Difference 3a for Gauss        

          q4[0,iik,1]       = Power_NxA    ;sum of 3b for ON BEAM
          q4[1,iik,1]       = Power_Nx2A    ;sum of 3b for OFF BEAM
;          If not(keyword_set(read_gauss)) then  q4_gauss[0,iik,1] = Q3b_g
;          If not(keyword_set(read_gauss)) then  q4_gauss[1,iik,1] = Q3b_diff_g
          
          q4[0,iik,2]       = Q3c_x    ;sum of 3c for ON BEAM
          q4[1,iik,2]       = Q3c_x2    ;sum of 3c for OFF BEAM
;          If not(keyword_set(read_gauss)) then  q4_gauss[0,iik,2] = Q3c_g
;          If not(keyword_set(read_gauss)) then q4_gauss[1,iik,2] = Q3c_diff_g
          
          q4[0,iik,3]       = Q3d_x   ;sum of 3d for ON BEAM
          q4[1,iik,3]       = Q3d_x2   ;sum of 3d for OFF BEAM
;          If not(keyword_set(read_gauss)) then  q4_gauss[0,iik,3] = Q3d_g
;          If not(keyword_set(read_gauss)) then  q4_gauss[1,iik,3] = Q3d_diff_g
  
          q4[0,iik,4]       = NxA_3e  ;sum of 3e for ON BEAM
          q4[1,iik,4]       = NxB_3e   ;sum of 3e for OFF BEAM (FLIP!!)
;          If not(keyword_set(read_gauss)) then  q4_gauss[0,iik,4] = Q3e_g
;          If not(keyword_set(read_gauss)) then  q4_gauss[1,iik,4] = Q3e_diff_g
          
          q4[0,iik,5]       = Power_NxA_3f    ;sum of 3f for ON BEAM
          q4[1,iik,5]       = Power_NxB_3f   ;sum of 3f for OFF BEAM (FLIP!!)
;          If not(keyword_set(read_gauss)) then   q4_gauss[0,iik,5] = Q3f_g
;          If not(keyword_set(read_gauss)) then  q4_gauss[1,iik,5] = Q3f_diff_g
      endif ;End Q4     
     endfor ;end Thresholds sigma values  
     
     ;----------------------------------
     ;           Gaussian
     ;----------------------------------
     ;Save the Gausssian
;     If not(keyword_set(read_gauss)) then begin
;       n = n_elements(y)
;       n_label = strtrim(string(n),1)
;       rebin_str = strtrim(string(rebin_time2),1)
;       ;save Gaussian
;       save,q4_gauss,n,threshold,filename='Q4_Gauss_'+n_label+'_'+rebin_str+'.sav'
;       file_g = 'Q4_Gauss_'+n_label+'_'+rebin_str+'.sav'
;     endif
;     
;     ;Read the Gaussian
;     If keyword_set(read_gauss) and keyword_set(file_G) then begin
;       print,'************************'
;       print,'Restore gaussian File: ', file_g
;       print,'Number of elements in y: ',n_elements(y)
;       print,'***********************'
;       restore,filename=file_g
;       ;Q4_gauss,n,threshold
;     Endif
;
;     ;Read the Gaussian
;     If keyword_set(read_gauss) and not(keyword_set(file_G)) then begin
;       n = n_elements(y)
;       n_label = strtrim(string(n),1)
;       rebin_str = strtrim(string(rebin_time2),1)
;       file_g = 'Q4_Gauss_'+n_label+'_'+rebin_str+'.sav'
;       print,'************************'
;       print,'Create File G Name'
;       print,'Restore gaussian File: ', file_g
;       print,'Number of elements in y: ',n_elements(y)
;       print,'***********************'
;       restore,filename=file_g
;       ;Q4_gauss,n,threshold
;     Endif
;
;     If not(keyword_set(read_gauss)) and skip_gauss eq 1 then read_gauss = 1   ;override (read in file next time)
     
     ;---------
     ;Q4 Plots
     ;--------
     If keyword_set(Q4_run) then begin 
       ;----------------------------
       ;           Q4
       ;----------------------------
       label_Q4 = 'S'+pol_label+' '+freq_label  
        ;------------
        ;Q4a (for 3a)
        ;------------  
        ;multiplot       
        ymax = max([q4[0,*,0],q4[1,*,0],q4_gauss[0,*,0]])*axis_scale
        ymin = min([q4[0,*,0],q4[1,*,0],q4_gauss[0,*,0]])*axis_scale
        If not(keyword_set(PAPERPLOTS)) then label = 'Q4a: '+label_main11 else label = 'Q4a' 
        cgplot,threshold,q4[0,*,0],xtitle='Threshold (units of sigma)',ytitle='Sum of Number of Peaks (#)',$
          title=label,/ystyle,/xstyle,yrange=[ymin,ymax],charsize=char_size
        cgplot,threshold,q4[1,*,0],/overplot,color='red'  ;OFF  
        cgplot,threshold,q4_gauss[0,*,0],/overplot,LineStyle=2
        al_legend,[OnBeam_Label,OffBeam_Label],textcolor=['Black','red'],/right,/top,charsize=2.0,symsize=2
        
          ;Q4a Diff
          Q4a_diff = q4[0,*,0] - q4[1,*,0]
          
          Q4a_diff_g = dblarr(n_thres)
          Q4a_diff_g(*) = q4_gauss[1,*,0]
          ymin = min([reform(Q4a_diff),-Q4a_diff_g])*axis_scale
          ymax = max([reform(Q4a_diff),Q4a_diff_g])*axis_scale
	        If not(keyword_set(PAPERPLOTS)) then label = 'Q4a (ON-OFF): '+label_main11 else label = 'Q4a (ON-OFF)'
          cgplot,threshold,Q4a_diff,xtitle='Threshold (units of sigma)',ytitle='Sum of Number of Peaks (#)',$
            title=label,/ystyle,/xstyle,yrange=[ymin,ymax],charsize=char_size
          cgplot,threshold,Q4a_diff_g,/overplot,LineStyle=2  
          cgplot,threshold,Q4a_diff_g*2.0,/overplot,LineStyle=2 
          cgplot,threshold,Q4a_diff_g*3.0,/overplot,LineStyle=2
          cgplot,threshold,-Q4a_diff_g,/overplot,LineStyle=2  
          cgplot,threshold,-Q4a_diff_g*2.0,/overplot,LineStyle=2
          cgplot,threshold,-Q4a_diff_g*3.0,/overplot,LineStyle=2
          cgplot,[0,100],[0,0],/overplot
          al_legend,[OnBeam_Label+' - '+OffBeam_Label],textcolor=['Black'],/right,/top,charsize=1.9

          ;-----------------
          ;Numerical Values
          ;-----------------
          calc_ppextrav2,threshold,Q4a_diff,gauss=Q4a_diff_g,/RunQ4,/Q4a,out=Q4aV,file_g=file_g,MINT=MINT,MAXT=MAXT

        ;------------
        ;Q4b (for 3b)
        ;------------
        ymin = min([q4[0,*,1],q4[1,*,1]])
        ymax = max([q4[0,*,1],q4[1,*,1]])
        If not(keyword_set(PAPERPLOTS)) then label = 'Q4b: '+label_main11 else label ='Q4b'
        cgplot,threshold,q4[0,*,1],xtitle='Threshold (units of sigma)',ytitle='Sum of the Power of Peaks (sigma)',title=label,$
          /ystyle,/xstyle,yrange=[ymin,ymax],charsize=char_size
        cgplot,threshold,q4[1,*,1],/overplot,color='red'  ;OFF
        cgplot,threshold,q4_gauss[0,*,1],/overplot,LineStyle=2
        al_legend,[OnBeam_Label,OffBeam_Label],textcolor=['Black','red'],/right,/top,charsize=2.0,symsize=2
        
          ;Q4b Diff
          Q4b_diff = q4[0,*,1] - q4[1,*,1]
          Q4b_diff_g = dblarr(n_thres)  
          Q4b_diff_g(*) =   q4_gauss[1,*,1]
          ymin = min([reform(Q4b_diff),-Q4b_diff_g])*axis_scale
          ymax = max([reform(Q4b_diff),Q4b_diff_g])*axis_scale
    
           If not(keyword_set(PAPERPLOTS)) then label ='Q4b (ON-OFF): '+label_main11 else label = 'Q4b (ON-OFF)'
          cgplot,threshold,Q4b_diff,xtitle='Threshold (units of sigma)',ytitle='Sum of the Power of Peaks (sigma)',$
            title=label,/ystyle,/xstyle,yrange=[ymin,ymax],charsize=char_size
          cgplot,threshold,Q4b_diff_g,/overplot,LineStyle=2  
          cgplot,threshold,Q4b_diff_g*2.0,/overplot,LineStyle=2
          cgplot,threshold,Q4b_diff_g*3.0,/overplot,LineStyle=2
          cgplot,threshold,-Q4b_diff_g,/overplot,LineStyle=2  
          cgplot,threshold,-Q4b_diff_g*2.0,/overplot,LineStyle=2 
          cgplot,threshold,-Q4b_diff_g*3.0,/overplot,LineStyle=2 
          cgplot,[0,100],[0,0],/overplot
          al_legend,[OnBeam_Label+' - '+OffBeam_Label],textcolor=['Black'],/right,/top,charsize=1.9     
                
	        ;-----------------
          ;Numerical Values
          ;-----------------
          calc_ppextrav2,threshold,Q4b_diff,gauss=Q4b_diff_g,/RunQ4,/Q4b,out=Q4bV,file_g=file_g,MINT=MINT,MAXT=MAXT
         
        ;------------
        ;Q4c (for 3c)
        ;------------
        ymin = min([q4[0,*,2],q4[1,*,2]])
        ymax = max([q4[0,*,2],q4[1,*,2]])
        If not(keyword_set(PAPERPLOTS)) then label='Q4c: '+label_main11 else label = 'Q4c'
        cgplot,threshold,q4[0,*,2],xtitle='Threshold (units of sigma)',ytitle='Sum of the Peak Asymmetry (#)',title=label,/ystyle,/xstyle,$
          yrange =[ymin,ymax],charsize=char_size 
        cgplot,threshold,q4[1,*,2],/overplot,color='red'  ;OFF
        cgplot,threshold,q4_gauss[0,*,2],/overplot,LineStyle=2
         cgplot,[0,100],[0,0],/overplot
        al_legend,[OnBeam_Label,OffBeam_Label],textcolor=['Black','red'],/right,/top,charsize=2.0,symsize=2
    
           ;Q4c Diff
          Q4c_diff = q4[0,*,2] - q4[1,*,2]
          Q4c_diff_g = dblarr(n_thres)
          Q4c_diff_g(*) =   q4_gauss[1,*,2]
          ymin = min([reform(Q4c_diff),-Q4c_diff_g])*axis_scale
          ymax = max([reform(Q4c_diff),Q4c_diff_g])*axis_scale
           If not(keyword_set(PAPERPLOTS)) then label ='Q4c (ON-OFF): '+label_main11 else label = 'Q4c (ON-OFF)'
          cgplot,threshold,Q4c_diff,xtitle='Threshold (units of sigma)',ytitle='Sum of the Peak Asymmetry (#)',$
            title=label,/ystyle,/xstyle,yrange=[ymin,ymax],charsize=char_size
          cgplot,threshold,Q4c_diff_g,/overplot,LineStyle=2
          cgplot,threshold,Q4c_diff_g*2.0,/overplot,LineStyle=2
          cgplot,threshold,Q4c_diff_g*3.0,/overplot,LineStyle=2
          cgplot,threshold,-Q4c_diff_g,/overplot,LineStyle=2
          cgplot,threshold,-Q4c_diff_g*2.0,/overplot,LineStyle=2
          cgplot,threshold,-Q4c_diff_g*3.0,/overplot,LineStyle=2
          cgplot,[0,100],[0,0],/overplot
          al_legend,[OnBeam_Label+' - '+OffBeam_Label],textcolor=['Black'],/right,/top,charsize=1.9
          
          ;-----------------
          ;Numerical Values
          ;-----------------
          calc_ppextrav2,threshold,Q4c_diff,gauss=Q4c_diff_g,/RunQ4,/Q4c,out=Q4cV,file_g=file_g,MINT=MINT,MAXT=MAXT

        ;---------------
        ;Q4d (for 3d)
        ;---------------
        ymin = min([q4[0,*,3],q4[1,*,3]])
        ymax = max([q4[0,*,3],q4[1,*,3]])
        If not(keyword_set(PAPERPLOTS)) then label = 'Q4d: '+label_main11 else label = 'Q4d'
        cgplot,threshold,q4[0,*,3],xtitle='Threshold (units of sigma)',ytitle='Sum of the Power Asymmetry (sigma)',title=label,/ystyle,/xstyle,$
          yrange=[ymin,ymax],charsize=char_size
        cgplot,threshold,q4[1,*,3],/overplot,color='red'  ;OFF
        cgplot,threshold,q4_gauss[0,*,3],/overplot,LineStyle=2
         cgplot,[0,100],[0,0],/overplot
        al_legend,[OnBeam_Label,OffBeam_Label],textcolor=['Black','red'],/right,/top,charsize=2.0,symsize=2
        
          ;Q4d Diff
          Q4d_diff = q4[0,*,3] - q4[1,*,3]
          Q4d_diff_g = dblarr(n_thres)
          Q4d_diff_g(*) =   q4_gauss[1,*,3]
          ;print,'Q4d Diff G',Q4d_diff_g[0]
          ymin = min([reform(Q4d_diff),-Q4d_diff_g])*axis_scale
          ymax = max([reform(Q4d_diff),Q4d_diff_g])*axis_scale
           If not(keyword_set(PAPERPLOTS)) then label ='Q4d (ON-OFF): '+label_main11 else label = 'Q4d (ON-OFF)'
          cgplot,threshold,Q4d_diff,xtitle='Threshold (units of sigma)',ytitle='Sum of the Power Asymmetry (sigma)',$
            title=label,/ystyle,/xstyle,yrange=[ymin,ymax],charsize=char_size
          cgplot,threshold,Q4d_diff_g,/overplot,LineStyle=2
          cgplot,threshold,Q4d_diff_g*2,/overplot,LineStyle=2
          cgplot,threshold,Q4d_diff_g*3,/overplot,LineStyle=2
          cgplot,threshold,-Q4d_diff_g,/overplot,LineStyle=2
          cgplot,threshold,-Q4d_diff_g*2.0,/overplot,LineStyle=2
          cgplot,threshold,-Q4d_diff_g*3.0,/overplot,LineStyle=2
          cgplot,[0,100],[0,0],/overplot
          al_legend,[OnBeam_Label+' - '+OffBeam_Label],textcolor=['Black'],/right,/top,charsize=1.9
          
          ;-----------------
          ;Numerical Values
          ;-----------------
          calc_ppextrav2,threshold,Q4d_diff,gauss=Q4d_diff_g,/RunQ4,/Q4d,out=Q4dV,file_g=file_g,MINT=MINT,MAXT=MAXT
                    
        ;------------
        ;Q4e (for 3e)
        ;------------
        ymin = min([q4[0,*,4],q4[1,*,4]])
        ymax = max([q4[0,*,4],q4[1,*,4]])
	      If not(keyword_set(PAPERPLOTS)) then label ='Q4e: '+label_main11 else label = 'Q4e'
        cgplot,threshold,q4[0,*,4],xtitle='Threshold (units of sigma)',ytitle='Sum of the Peak Offset (#)',title=label,$
          /ystyle,/xstyle,yrange=[ymin,ymax],charsize=char_size
        cgplot,threshold,q4[1,*,4],/overplot,color='red'  ;OFF
        cgplot,threshold,q4_gauss[0,*,4],/overplot,LineStyle=2
        al_legend,[OnBeam_Label,OffBeam_Label],textcolor=['Black','red'],/right,/top,charsize=2.0,symsize=2
     
          ;Q4e Diff 
          Q4e_diff = q4[0,*,4] - q4[1,*,4]
          Q4e_diff_g    = dblarr(n_thres)
          Q4e_diff_g(*) =   q4_gauss[1,*,4]
          ymin = min([reform(Q4e_diff),-Q4e_diff_g])*axis_scale
          ymax = max([reform(Q4e_diff),Q4e_diff_g])*axis_scale
           If not(keyword_set(PAPERPLOTS)) then label ='Q4e (ON-OFF): '+label_main11 else label = 'Q4e (ON-OFF)'
          cgplot,threshold,Q4e_diff,xtitle='Threshold (units of sigma)',ytitle='Sum of the Peak Offset (#)',$
            title=label,/ystyle,/xstyle,yrange=[ymin,ymax],charsize=char_size
          cgplot,threshold,Q4e_diff_g,/overplot,LineStyle=2
          cgplot,threshold,Q4e_diff_g*2,/overplot,LineStyle=2
          cgplot,threshold,Q4e_diff_g*3,/overplot,LineStyle=2
          cgplot,threshold,-Q4e_diff_g,/overplot,LineStyle=2
          cgplot,threshold,-Q4e_diff_g*2.0,/overplot,LineStyle=2
          cgplot,threshold,-Q4e_diff_g*3.0,/overplot,LineStyle=2
          cgplot,[0,100],[0,0],/overplot
          al_legend,[OnBeam_Label+' - '+OffBeam_Label],textcolor=['Black'],/right,/top,charsize=1.9

          ;-----------------
          ;Numerical Values
          ;-----------------
          calc_ppextrav2,threshold,Q4e_diff,gauss=Q4e_diff_g,/RunQ4,/Q4e,out=Q4eV,file_g=file_g,MINT=MINT,MAXT=MAXT
         
        ;------------- 
        ;Q4f (for 3f)
        ;-------------
        print,'--------- Q4f -------------'
        ymin = min([q4[0,*,5],q4[1,*,5],q4_gauss[0,*,5]])
        ymax = max([q4[0,*,5],q4[1,*,5],q4_gauss[0,*,5]])
 	      If not(keyword_set(PAPERPLOTS)) then label ='Q4f: '+label_main11 else label = 'Q4f'
        cgplot,threshold,q4[0,*,5],xtitle='Threshold (units of sigma)',ytitle='Sum of the Power Offset (sigma)',title=label,/ystyle,/xstyle,$
           yrange=[ymin,ymax],charsize=char_size
        cgplot,threshold,q4[1,*,5],/overplot,color='red'  ;OFF
        cgplot,threshold,q4_gauss[0,*,5],/overplot,LineStyle=2
        al_legend,[OnBeam_Label,OffBeam_Label],textcolor=['Black','red'],/right,/top,charsize=2.0,symsize=2
    
        ;Q4f Diff
        Q4f_diff = q4[0,*,5] - q4[1,*,5]
        Q4f_diff_g    = dblarr(n_thres)
        Q4f_diff_g(*) =   q4_gauss[1,*,5]
        ymin = min([reform(Q4f_diff),-Q4f_diff_g])*axis_scale
        ymax = max([reform(Q4f_diff),Q4f_diff_g])*axis_scale   
        If not(keyword_set(PAPERPLOTS)) then label ='Q4f (ON-OFF): '+label_main11 else label = 'Q4f (ON-OFF)'
        cgplot,threshold,q4f_diff,xtitle='Threshold (units of sigma)',ytitle='Sum of the Power Offset (sigma)',title=label,$
            /ystyle,/xstyle,yrange=[ymin,ymax],charsize=char_size,xrange=[min(threshold),max(threshold)]
        cgplot,threshold,Q4f_diff_g,/overplot,LineStyle=2
        cgplot,threshold,Q4f_diff_g*2.0,/overplot,LineStyle=2
        cgplot,threshold,Q4f_diff_g*3.0,/overplot,LineStyle=2
        cgplot,threshold,-Q4f_diff_g,/overplot,LineStyle=2
        cgplot,threshold,-Q4f_diff_g*2,/overplot,LineStyle=2
        cgplot,threshold,-Q4f_diff_g*3,/overplot,LineStyle=2
        cgplot,[0,100],[0,0],/overplot
	      al_legend,[OnBeam_Label+' - '+OffBeam_Label],textcolor=['Black'],/right,/top,charsize=1.9
	      
          ;-----------------
          ;Numerical Values
          ;-----------------
          ;test_detection,Q4f_diff,THRESHOLD,Q4f_diff_g,/save_gaus,/VERBOSE,TEST_AVG=TEST_AVG_Q4f,PROP_FALSE=PROP_FALSE_Q4f 
          calc_ppextrav2,threshold,Q4f_diff,gauss=Q4f_diff_g,/RunQ4,/Q4f,out=Q4fV,file_g=file_g,MINT=MINT,MAXT=MAXT
          print,'----------------------'
      endif ;End Q4 
      
      ;********************************************************************************
      ; 			   Q3 (Qi vs time)
      ;********************************************************************************
      If keyword_set(Q3) then begin
        bins       = n_elements(x_q1)
        NxA        = dblarr(bins) ;3a (on) ;number of peaks for on beam above theshold
        Nx2A       = dblarr(bins) ;3a (off);number of peaks for off beam above theshold
        Q3a_diff   = dblarr(bins) ;3a diff
        NA_g       = dblarr(bins) ;3a gaussian
        NA_g2      = dblarr(bins) ;3a gaussian 2

        Power_NxA  = dblarr(bins) ;3b (on) ;power of peaks for on beam
        Power_Nx2A = dblarr(bins) ;3b (off);power of peaks for off beam
        Power_NA_g = dblarr(bins) ;3b gaussian
        Power_NA_g2 = dblarr(bins) ;3b gaussian 2

        NxB        = dblarr(bins) ;3c (Number negative peaks; ON)
        Nx2B       = dblarr(bins) ;3c (Number negative peaks; Off)
        NB_g       = dblarr(bins) ;3c gaussian
        NB_g2      = dblarr(bins) ;3b gaussian 2

        Power_NxB  = dblarr(bins) ;3c (Power negative peaks; ON)
        Power_Nx2B = dblarr(bins) ;3c (Power negative peaks; OFF)
        Power_NB_g = dblarr(bins) ;3c gaussian
        Power_NB_g2 = dblarr(bins) ;3c gaussian 2

        NxA_3e       = dblarr(bins) ;3e (number of peaks above threshold and twice off values)
        NxB_3e       = dblarr(bins) ;3e (number of peaks above threshold for OFF and twice ON values)
        NA_3e_g    = dblarr(bins)
        Power_NxA_3f = dblarr(bins) ;3f (power of 3e: ON)
        Power_NxB_3f = dblarr(bins) ;3f (power of 3e OFF)
        Power_NA_g_3f = dblarr(bins)
        TI_win    = long(TI/rebin_time2)  ;Amount of elements that TI covers

        ;-----------------------------------------------------------
        ;                       gaussian
        ;-----------------------------------------------------------
        gauss_lofar_postprocessing,gauss_steps,TI_win,tau_Qplot,slope=slope,$
          Q3a=Q3a_g,dQ3a=Q3a_g_diff,$
          Q3b=Q3b_g,dQ3b=Q3b_g_diff,$
          Q3c=Q3c_g,dQ3c=Q3c_g_diff,$
          Q3d=Q3d_g,dQ3d=Q3d_g_diff,$
          Q3e=Q3e_g,dQ3e=Q3e_g_diff,$
          Q3f=Q3f_g,dQ3f=Q3f_g_diff,square=square

        ;---------------------------------------------------
        ;               Start Data
        ;--------------------------------------------------
        for i=0,bins-1 do begin
          end_array = (i+1)*TI_win
          If i eq bins-1. then begin
            If end_array gt n_elements(y) then end_array = n_elements(y)  ;in case of not perfect coverage
          Endif
          ytemp  = y[i*TI_win:end_array-1.]              ;ON
          ytemp2 = y2[i*TI_win:end_array-1.]             ;OFF

          postprocessing_calc,ytemp,ytemp2,tau_Qplot,slope=slope,$
                             NxA=NxAo,PNxA=Power_NxAo,$
                             Nx2A=Nx2Ao,PNx2A=Power_Nx2Ao,$
                             Nx3Ae=NxA_3eo,PNAx3f=Power_NxA_3fo,$
                             Nx3Be=NxB_3eo,PNBx3f=Power_NxB_3fo,$
                             NxB=NxBo,PNxB=Power_NxBo,$
                             Nx2B=Nx2Bo,PN2xB=Power_Nx2Bo
           
           NxA(i)          =NxAo
           Power_NxA(i)    =Power_NxAo 
           Nx2A(i)         =Nx2Ao
           Power_Nx2A(i)   =Power_Nx2Ao
           NxA_3e(i)       =NxA_3eo
           Power_NxA_3f(i) =Power_NxA_3fo
           NxB_3e(i)       =NxB_3eo                          
           Power_NxB_3f(i) =Power_NxB_3fo
           NxB(i)          =NxBo
           Nx2B(i)         =Nx2Bo
           Power_NxB(i)    =Power_NxBo
           Power_Nx2B(i)   =Power_Nx2Bo       
        endfor ;end bins

        ;------------------------------------------------
        ;            Q3a: Number of peaks above tau
        ;------------------------------------------------
        Q3a_diff = 0 
        Q3a_diff = NxA - Nx2A

        If not(keyword_set(PAPERPLOTS)) then label='Q3a: '+label_main12 else label='Q3a'
        ytitle = 'Number of Peaks > '+cgGreek('tau')+'='+N_label
          diagnostic_plots,NxA,Nx2A,MJD_rebin1_q1,xt_rebin1_q1,P,T0,label,ytitle,t2=MJD_rebin2_q1,$
            /single_plot,Q_Gauss = Q3a_g,/line,OnBeam_Label=OnBeam_Label,OffBeam_Label=OffBeam_Label;,UT=UT_q1
        
        If keyword_set(UT) then begin
          diagnostic_plots,NxA,Nx2A,MJD_rebin1_q1,xt_rebin1_q1,P,T0,label,ytitle,t2=MJD_rebin2_q1,$
            /single_plot,Q_Gauss = Q3a_g,/line,OnBeam_Label=OnBeam_Label,OffBeam_Label=OffBeam_Label,UT=UT_rebin_q1
        endif 

          label ='Q3a (ON-OFF): '+label_main12 
          diagnostic_plots,Q3a_diff,Q3a_diff,MJD_rebin1_q1,xt_rebin1_q1,P,T0,label,ytitle,t2=MJD_rebin2_q1,$
            /single_plot,/TARGET_ONLY,Q_Gauss =Q3a_g_diff,/NEG_GAUSS,/line,OnBeam_Label=OnBeam_Label,OffBeam_Label=OffBeam_Label,/ALL_Gauss;,UT=UT_q1

        ;------------------------------------------------
        ;            Q3b: Power of sum of peaks above tau
        ;------------------------------------------------
        Q3b_diff = 0 
        Q3b_diff = Power_NxA - Power_Nx2A
          label = 'Q3b: '+label_main12
          ytitle = 'Sum of the Power of Peaks > '+cgGreek('tau')+'='+N_label
          diagnostic_plots,Power_NxA,Power_Nx2A,MJD_rebin1_q1,xt_rebin1_q1,P,T0,label,ytitle,t2=MJD_rebin2_q1,$
            /single_plot,Q_Gauss = Q3b_g,/line

          label='Q3b (ON-OFF): '+label_main12
          diagnostic_plots,Q3b_diff,Q3b_diff,MJD_rebin1_q1,xt_rebin1_q1,P,T0,label,ytitle,t2=MJD_rebin2_q1,$
            /single_plot,/TARGET_ONLY, Q_Gauss =  Q3b_g_diff,/NEG_GAUSS,/line,/ALL_Gauss

        ;----------------------------------------------------------
        ; Q3c: Difference of N of peaks above tau and below - tau
        ;----------------------------------------------------------
        Q3c_x = 0 
        Q3c_x2 = 0 
        Q3c_x  = NXA  - NxB   ;[(N > tau) - (N < - tau)](time) - On
        Q3c_x2 = Nx2A - Nx2B  ;[(N > tau) - (N < - tau)](time) - Off
        Q3c_diff = Q3c_x - Q3c_x2
          label = 'Q3c: '+label_main12
          ytitle = 'Sum of the Peak Asymmetry (#)'
          diagnostic_plots,Q3c_x,Q3c_x2,MJD_rebin1_q1,xt_rebin1_q1,P,T0,label,ytitle,t2=MJD_rebin2_q1,$
            /single_plot,Q_Gauss = Q3c_g,/NEG_GAUSS,/line

           label = 'Q3c (ON-OFF): '+label_main12
          diagnostic_plots,Q3c_diff,Q3c_diff,MJD_rebin1_q1,xt_rebin1_q1,P,T0,label,ytitle,t2=MJD_rebin2_q1,$
            /single_plot,/TARGET_ONLY, Q_Gauss =  Q3c_g_diff,/NEG_GAUSS,/line,/ALL_Gauss


        ;----------------------------------------------------------------------
        ; Q3d: Difference of Power of peaks above tau and below - tau
        ;----------------------------------------------------------------------
        Q3d_x = 0 
        Q3d_x2 = 0 
        Q3d_x  = Power_NXA  - Power_NxB   ;[Power(N > tau) - Power(N < - tau)](time) - On
        Q3d_x2 = Power_Nx2A - Power_Nx2B  ;[Power(N > tau) - Power(N < - tau)](time) - Off
        Q3d_diff=Q3d_x - Q3d_x2 
          label = 'Q3d: '+label_main12
          ;ytitle = 'Power Asymmetry [Power (N > '+cgGreek('tau')+') - (N < - '+cgGreek('tau')+')]'
           ytitle = 'Sum of the Power Asymmetry (#)'
          diagnostic_plots,Q3d_x,Q3d_x2,MJD_rebin1_q1,xt_rebin1_q1,P,T0,label,ytitle,t2=MJD_rebin2_q1,$
            /single_plot,Q_Gauss =  Q3d_g,/NEG_GAUSS,/line

            label = 'Q3d (ON-OFF): '+label_main12
          diagnostic_plots,Q3d_diff,Q3d_diff,MJD_rebin1_q1,xt_rebin1_q1,P,T0,label,ytitle,t2=MJD_rebin2_q1,$
            /single_plot,/TARGET_ONLY, Q_Gauss =  Q3d_g_diff,/NEG_GAUSS,/line,/ALL_Gauss

        ;-------------------------------------------------------------------------
        ;Q3e:  Number of peaks above tau and twice the off Values
        ;-------------------------------------------------------------------------
        Q3e_diff = NxA_3e - NxB_3e
        label = 'Q3e: '+label_main12
        ytitle = 'Sum of the Peak Offset (#)'
          diagnostic_plots,NxA_3e,NxB_3e,MJD_rebin1_q1,xt_rebin1_q1,P,T0,label,ytitle,t2=MJD_rebin2_q1,$
            /single_plot,Q_Gauss = Q3e_g,/line

            label = 'Q3e (ON-OFF): '+label_main12
          diagnostic_plots,Q3e_diff,Q3e_diff,MJD_rebin1_q1,xt_rebin1_q1,P,T0,label,ytitle,t2=MJD_rebin2_q1,$
            /single_plot,/TARGET_ONLY, Q_Gauss =  Q3e_g_diff,/NEG_GAUSS,/line,/ALL_Gauss  
        
        ;-------------------------------------------------------------------------
        ;Q3f:  Power of the Number of peaks above tau and twice the off Values
        ;-------------------------------------------------------------------------
        Q3f_diff = Power_NxA_3f-Power_NxB_3f
        label = 'Q3f: '+label_main12
        ytitle = 'Sum of the Power Offset (#)'
          diagnostic_plots,Power_NxA_3f,Power_NxB_3f,MJD_rebin1_q1,xt_rebin1_q1,P,T0,label,ytitle,t2=MJD_rebin2_q1,$
            /single_plot,Q_Gauss =  Q3f_g,/line
         
            label = 'Q3f (ON-OFF): '+label_main12
          diagnostic_plots,Q3f_diff,Q3f_diff,MJD_rebin1_q1,xt_rebin1_q1,P,T0,label,ytitle,t2=MJD_rebin2_q1,$
            /single_plot,/TARGET_ONLY, Q_Gauss =  Q3f_g_diff,/NEG_GAUSS,/line,/ALL_Gauss 
    
        ;-----------------------------------------------------------------
 	; Dyn Spec high pass filter  (tau) 
	;-----------------------------------------------------------------
            ;STANDARD_PLOTS,y_dynr,px,xt_rebin1,xf_rebin_q1,label,/VERBOSE,/ONLY_PLOT,skyplot=y2_dynr,/NO_ZOOM,xunit = 'Days',$
             ; fracmax=1.00,fracmin=0.05,ytitle_bar='Intensity (Arbitrary units)'
             label = 'High-Pass Filter; '+cgGreek('tau')+':'+strtrim(string(cgnumber_formatter(new_sigma,decimals=1)),1)+' '+label_main11
             If not(keyword_set(UT))then begin 
               STANDARD_PLOTS,y_dynr,px,xt_rebin1,xf_rebin_q1,label,mask=mask,/VERBOSE,/ONLY_PLOT,skyplot=y2_dynr,$
	       /NO_ZOOM,xunit = 'Hour',fracmax=1.00,fracmin=0.05,ytitle_bar='Intensity (Arbitrary units)',/panel_ps
             endif 
             If keyword_set(UT) then begin 
               STANDARD_PLOTS,y_dynr,px,UT_rebin,xf_rebin_q1,label,mask=mask,/VERBOSE,/ONLY_PLOT,skyplot=y2_dynr,$
               /NO_ZOOM,xlabel='UT', xunit = 'hour',fracmax=1.00,fracmin=0.05,ytitle_bar='Intensity (Arbitrary units)',/panel_ps
             endif 
      endif ;End Q3
      
      !p.multi=[0,1,1]
      cgPS_Close
      cgfixps,filename_ps+'.ps'    
         
      ;Create panels pdfs
      spawn,'ps2pdf -dPDFSETTINGS=/ebook '+filename_ps+'.ps'
      panel_pdf,filename_ps,Q3=Q3
      spawn,'rm '+filename_ps+'.ps'
     If keyword_set(REMOVE) then spawn,'rm '+filename_ps+'.pdf'
      
      ;------------------
      ;Save Output file
      ;------------------ 
      If not(keyword_set(date)) then date = strmid(filename_ps,0,7)
      lofar_postprocess_textfilev2,filename_ps,Date=Date,S=pol_label,fmin=fmins[iii],fmax=fmaxs[iii],BT=Target_Beam,BSky=Sky_Beam,$
    	Q1aV=Q1aV,Q1bV=Q1bV,Q4aV=Q4aV,Q4bV=Q4bV,Q4cV=Q4CV,Q4dV=Q4dV,Q4eV=Q4eV,Q4fV=Q4fV,$
        /WRITE        
       
      ;Save Q1,2,4 values
      If keyword_set(save_Q) then begin
        If keyword_set(Q4_run) and keyword_set(Q2) and keyword_set(Q1) then begin
          xt = xt_rebin1_q1
          xf = xf_rebin_q1
          save,threshold,xt,xf,q1,q2,q4,q4_gauss,phase_rebin,MJD_rebin,filename=filename_ps+'_PP_Qs.sav'
        endif
        If keyword_set(Q4_run) and keyword_set(Q2) and keyword_set(Q1) and keyword_set(Q3) then begin
          xt = xt_rebin1
          xf = xf_rebin_q1
          save,threshold,xt,xf,q1,q2,q3,q4,q4_gauss,phase_rebin,MJD_rebin,filename=filename_ps+'_PP_Qs.sav'
        endif
        
        skip:
      endif ;saveQ
   endfor ;end Ellptical Correction Loop  
  endfor ;pol loop (differnet ways to search)
 endfor ;end threshold(mask) loop
endfor ;end freq loop
 
  print,'SYSTIME for Post-Processing =',SYSTIME(1) - stp, ' sec'
  PRINT, 'Memory required for Post-Processing: ', (MEMORY(/HIGHWATER) - startmem_stp)/1d9 ,' Gb'
  print,'The End of Post-Processing Program'
  return 
end
