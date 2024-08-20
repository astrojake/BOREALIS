;************************************************************************************************
;************************************************************************************************
;         Main Program for the BeamfOrmed Radio Emission AnaLysIS (BOREALIS) pipeline 
;************************************************************************************************
;************************************************************************************************
;Authors: Jake Turner (UVA/Cornell)
;        J.M. Griessmeier (CNRS) and Philippe Zarka (Obspm)
;
;Inputs: 
;       Input_file  ; Name of the Input file (See README_InputFile.dat) ;       
;
;Purpose: 
;         Analyze beamdformed data. 
;         BOREALIS works on LOFAR, NenuFAR, and dynspecMS (beams extracted from data.
;         It can be generalized to work with any beamformed data
;       
;Citations: 
;         If you use this program, please cite: 
;          - Turner et al. 2017 (https://arxiv.org/abs/1710.04997)
;          - Turner et al. 2019 (https://arxiv.org/abs/1802.07316)
;          - Turner et al. 2021 (https://arxiv.org/abs/2012.07926)      
;       
;Uses: 
;  Internal (our group)
;       lofar_rfi.pro           :Runs Pre-Procesing and runs RFI 
;       lofar_processing.pro    :Processes Data (Finds time-frequency norm curve, normalizes, apply RFI, and rebins data)
;       lofar_fft.pro           :Performs an FFT
;       lofar_postprocessing.pro :Preforms Post-Processing 
;       lofar_combinemask.pro    ;Combines masks from different beams 
;       lofar_combine.pro        ;Combines data from multiple dates for post-processing 
;       standard_plots.pro      :Stadard Plotting routine for dynspec
;      -Programs from Phillip (rfi_mitigate and subprograms)
;        - rfi_mitigate.pro: has been updated to include plots
;          -(Programs within)  
;      -Programs from Ian's pipeline     
;  External Programs:          
;      -hdf5 programs (are downloaded on nancep machines)
;      -cgcoyate programs
;      -idlastro lib       
;
;Outputs (Examples):
;       Log file                                                  (e.g. log_L429868_beam_1_pol0_3523sec_start2_steps84_patrol5.5_lesig3.5_pex_sum_RFI0_P1_2016-11-03T12:19:07Z.dat)
;       RFI Plots                                                 (e.g. L429868_pol0_3523sec_RFI_patrol5.5_pex_lesig3.5_sum_beam1.ps)   
;       RFI mask (if in bits)                                     (e.g. mask_full_bit_L429868_pol0_15770sec_RFI_patrol5.5_pex_lesig3.5_sum_beam0.sav)
;       RFI mask (if in bytes)                                    (e.g. mask_full_byte_L429868_pol0_15770sec_RFI_patrol5.5_pex_lesig3.5_sum_beam0.sav)
;       TimeFreq-Correction Save file                             (e.g. 3523sec_RFI_patrol5.5_pex_lesig3.5_sum_beam1_TimeFreqCorrection.sav)
;       Processing File 1: Finding Time-F Normal                  (e.g. L429868_pol0_3523sec_RFI_patrol5.5_pex_lesig3.5_sum_beam1_processing_TF.ps)
;       Processing File 2: Rebin                                  (e.g. L429868_pol0_3523sec_RFI_patrol5.5_pex_lesig3.5_sum_beam1_processing_rebin.ps)
;       Processing File 3: Final mask & rebinned normalized data  (e.g. L429868_pol0_3523sec_RFI_patrol5.5_pex_lesig3.5_sum_beam1_processing_final.ps)
;       FFT-Processing Plots                                      (e.g. L429868_pol0_3523sec_RFI_patrol5.5_pex_lesig3.5_sum_beam1_fft_pulsar.ps)
;       Post-Processing Plots                                     (e.g. L429868_pol0_3523sec_RFI_patrol5.5_pex_lesig3.5_sum_beam1_26.1to40.0MHz_thresT0.10_thresF0.10_postprocessing.ps)
;       Post-Processing: Thresholds                               (e.g. L429868_pol0_3523sec_RFI_patrol5.5_pex_lesig3.5_sum_thresholds.dat)
; 
; - All Outputs will be written out in local directory
;
;TO RUN PROGRAM:
; Compile:
;  > .r LOFAR_BF_PipelineV4.pro
;  > LOFAR_BF_PipelineV4,'Inputv5.dat' 
;             
;-----------------
;Lofar File Format
;-----------------
;   L248553_SAP000_B000_S0_P000_bf.h5', L is the date
;   B000 - B003: 4 beams: 4 coherent (55cnc-exoplanet, pulsar, sky, bright source)
;   S0-S3:  4 polarizations: IQUV
;   
;Polarization INPUT (in Input File): 
; S = 0; Stokes-I; I is ran through code 
; S = 1; Stokes-Q; Abs(Q) is ran through code
; S = 2; Stokes-U; Abs(U) is ran through code
; S = 3; Stokes-V; Abs(V) is ran through code
; S = 4; Stokes-V; Vprime is ran through code, Abs(V), V+, V- is ran through PP 
; S = 5; L; L is ran through through code 
;
;MODIFCATION HISTORY: 
;                2017/Nov/30      Jake Turner (UVA) 
;                                 Beta version 1; V1  
;                2018/June/02     Jake Turner (UVA) 
;                                 Beta version 3; V2 (Updated with V Polarization)
;                2019/March/22    Jake Turner (Cornell) 
;                                 Version 4 (Beta); Updated with NENUFAR and DYNSPECMS 
;                                            Update Q, U, L polarization for LOFAR
;               2021/March/23    Jake Turner (Cornell)
;                                 Version 4 (named changed to BOREALIS); 
;                                           Updated with NENUFAR (new fits files) and DynspecMS
;---------------------------------------------------------------------------------------------------------------------------   
;---------------------------------------------------------------------------------------------------------------------------
pro BOREALISv4,input_file
       
;---------------------------------------------------------------------------------------
;                               Read Input File 
;---------------------------------------------------------------------------------------
read_inputv5,input_file,planet=planet,period=Period,TO=TO,date=date,fileroot=file_root,$
               RQUICKLOOK=run_QUICKLOOK,RRFI=Run_RFI,RProcessing=Run_processing,RFFT=run_FFT,RDISP=run_DISP,RPostProcessing=run_postprocessing,$
               PLOT_RAW=PLOT_RAW,DATA_NORM=DATA_NORM,DATA_SLOPE=DATA_SLOPE,fminql=fmin_ql,fmaxql=fmax_ql,qltmin=tmin_ql,N_timeql=N_time_ql,$
               qlsteps=steps_ql,bminql=bmin_ql,bmaxql=bmax_ql,Sminql=Smin_ql,Smaxql=Smax_ql,$
               start_beam=start_beam,end_beam=end_beam,$
               Polarization=S,tmin=tmin,nspectra=nspectra,steps=steps,$
               PS=PS,VERBOSE=VERBOSE,mask_bit=mask_bit,bit=bit,send_email=send_email,email_address=email_address,$
               sigma_patrol=patrol_value,Si_le_sig=le_sig_value,patrol=patrol,LE_SIG=LE_SIG,SUM=SUM,PEX=PEX,$
               TF_correction=TimeFreq_correction,Apply_Corrections=Apply_Corrections,$
               Drebin_data=rebin_data,Frebin_freq=rebin_freq,Trebin_time=rebin_time,$
               beam_DISP=beam_DISP,pulsarname=pulsarname,dm=dm,size_disp=size_disp,$
               beam_FFT=beam_FFT,PFFT=Period_FFT,pmin=pmin,fmaxp=fmaxp,fminp=fminp,$
               Target_Beam = Target_Beam, Sky_Beam=Sky_Beam,Q1=Q1,Q2=Q2,Q3=Q3,Q4=Q4,$
               REMOVE_BADF=REMOVE_BADF,$
               SELECT_Freq,fminimum=fmin,fmaximum=fmax,$
               quantile_norm,quantile=quantile,mask_save=mask_save,$
               NFREQ_PATROL=NFREQ_PATROL,NTIME_PATROL=NTIME_PATROL,NFREQ_LESIG=NFREQ_LESIG,NTIME_LESIG=NTIME_LESIG,$
               MF=MF,MT=MT,THRESHOLD=THRESHOLD,FPMIN=PMIN_F,ExpansionF1=EXPF,TPMIN=PMIN_T,ExpansionT1=EXPT,F2PMIN=PMIN_F2,EXPF2=EXPF2,T2PMIN=PMIN_T2,EXPT2=EXPT2,$
               RFI=RFI,full_clean_save=full_clean_save,save_datagain=save_datagain,$
               RUN_DISP_TIME=RUN_DISP_TIME, rebin_time_disp=rebin_time_disp,$
               ddm=ddm,$
               create_combine_mask=create_combine_mask,combine_mask=combine_mask,NUM_MASKS=NUM_MASKS,COMBO_BEAMS=COMBO_BEAMS,$
               Smooth_Surface=Smooth_Surface,N_smooth=N_smooth,Process_METHOD2=Process_METHOD2,$
               skip_rebinloop=skip_rebinloop,save_RFIMaskData=save_RFIMaskData,$
               TI=TI,Default_MASK=Default_MASK,Run_Mask=Run_Mask,p_threshold=threshold_mask,Default_FreqRange=Default_FreqRange,Run_FreqRange=Run_FreqRange,$
               freq_minextra1=freq_minextra1,freq_maxextra2=freq_maxextra2,Default_Time=Default_Time,QTmin=QTmin,QTmax=QTmax,Default_Threshold=Default_Threshold,Exthreshold_extra=threshold_extra,$
               max_threshold=max_threshold,pthreshold=threshold_p,tau_Qplot=tau_Qplot,slope=slope,Qwin=Qwin,ElC=ELLCORR,save_gauss=save_gauss,gauss_steps=gauss_steps,file_g=file_g,$
               root_gauss =root_gauss,save_Q=save_Q,OnBeam_Label=OnBeam_Label,OffBeam_Label=OffBeam_Label,ExB1=Extra_Beam1,ExB2=Extra_Beam2,$
               Run_Rebintime2=Run_Rebintime2,rebin_time2=rebin_time2,Run_RebinFreq2=Run_RebinFreq2,rebin_freq2=rebin_freq2,$
               file_root_combine=file_root_combine,filenname_combinedates=filenname_combinedates,$
               create_dates_combine=create_dates_combine,Combine_NBeams=Combine_NBeams,dates_combine,Combine_Beams=Combine_Beams,$
               RUN_JUPITER=RUN_JUPITER,input_file_J=input_file_J,$
               PREMOVE=REMOVE,PAPERPLOTS=PAPERPLOTS,MINT=MINT,MAXT=MAXT,$
               LOFAR=LOFAR,NENUFAR=NENUFAR,DYNSPECMS=DYNSPECMS,filename=filename,$
               COMPLEX_NAME=COMPLEX_NAME,PCA=PCA,Full_RFI=Full_RFI,use_combinemask=use_combinemask
  
  loop_polar = 0                  
  If S eq 6 then loop_polar = 1   ; loop both Stokes I and Stokes-V
  If loop_polar eq 1 then begin
    print,'------------------------------------'
    print,'Loop Stokes-I and Stokes-V polarizations'
    print,'------------------------------------'
  Endif
  
  for i=0,loop_polar do begin
    
    print,'I for loop',i
    ;------------------------
    ;Staring memory and time
    ;------------------------
    start_mem = MEMORY(/CURRENT) ;start memory
    st= SYSTIME(1)               ;start time
    binTime = BIN_DATE(systime(/ut))  ;
    timestamp_string = TIMESTAMP(YEAR = binTime[0], MONTH = binTime[1],$
      DAY = binTime[2], HOUR = binTime[3], MINUTE = binTime[4], SECOND =binTime[5])  ;2012-09-04T11:25:15Z
    print,'------------------------------------'
    print,'Program Start at ', timestamp_string
    print,'------------------------------------'
    
    ;---------------
    ;Read Header Info
    ;---------------
    If loop_polar eq 1 then begin 
      If i eq 0 then begin
       print,'------------------------------------'
       print,' ----- Running Stokes-I  -----------'
       print,'------------------------------------'
       S = 0    ;Stokes I 
      Endif
      If i eq 1 then begin
        print,'------------------------------------'
        print,' ----- Running Stokes-V   -----------'
        print,'------------------------------------'
        S = 4    ;Stokes I   Stokes-V prime 
      Endif
    Endif ; end loop polar   
    If keyword_set(LOFAR) then begin 
      If S lt 4 then S_in = S   ; S = S (I, U, Q, V)
      If S eq 4 then S_in = 3   ; V^2
      If S eq 5 then S_in = 3   ; V'
      If S eq 6 then S_in = 1   ; L = sqrt(Q^2 + U^2) 
      file = file_root+date+'/'+date+'_SAP000_B00'+strtrim(string(uint(start_beam)),1)+'_S'+strtrim(string(uint(S_in)),1)+'_P000_bf'
      READ_H5_DATA,file,0,1,fmin,fmax,start_beam,S_in,x,t,xf,n_t,nf,Time=tf,/NSPECTRA,$
        NOF_SAMPLES=NOF_SAMPLES,START_MJD=START_MJD,MJD_sampling_time=MJD_sampling_time,$
        ANTENNA_SET=ANTENNA_SET,SAMPLING_TIME=SAMPLING_TIME,NOF_stations=NOF_stations,$
        channel_width=channel_width,/READ_HEADER
    endif 
    If keyword_set(NENUFAR) then begin
      If keyword_set(RUN_QUICKLOOK) then beam_in = bmin_ql ELSE beam_in = start_beam
      If S eq 0 then Stokes = 0 ;Stokes-I
      If S eq 1 then Stokes = 1 ;Stokes-Q
      If S eq 2 then Stokes = 2 ;Stokes-U
      If S eq 3 then Stokes = 3 ;Stokes-V
      If S eq 4 then Stokes = 3 ;Stokes-V^2
      If S eq 5 then Stokes = 3 ;Stokes-V prime
      If S eq 0 then nstokes = 1 ELSE nstokes = 4
      temp       = filename+'*_'+strtrim(string(start_beam),1)+'.spectra*'
      file_input = findfile(temp)
      READ_NU_FITS, file_input, command,param,variab, nt,dt,nf,df,ns,jd0,h0,fref,$
        temporary(data3),xt,xf,beam,ndata,corrt,corrf,/nodata
        sampling_time = dt 
      ;read_nenufar, file_root+filename, temporary(data2),tt,ff,beam_in,ndata, ntt,sampling_time,nff,df,ns, Stokes,nstokes=nstokes,tmin=0,tmax=1,fmin=9,fmax=70,/VERBOSE,MINFREQ=min_freq,MAXFREQ=max_freq
    Endif
    If keyword_set(DYNSPECMS) then begin
      read_dynspec_fits,filename+'_beam'+strtrim(string(uint(start_beam)),1),$
        0,1,fmin,fmax,0,temporary(data2),xt,xf,nt,nf,dt=dt,df=df,/PRINT_HEAD
      sampling_time =dt
    Endif
   
    ;------------
    ;Time info
    ;------------
    tmin_org = tmin  ;save
    If nspectra gt 0 then N_time = sampling_time*(nspectra)
    total_time = N_time*(steps)
    
    ;----------
    ;Send Email 
    ;-----------
    email,email_address,'Program Started ',send=send_email
    
   ;--------------------------------------------------------------------------------------
   ;---------------------------------------------------------------------------------------
   ;         Don't need to change any free parameters below here 
   ;---------------------------------------------------------------------------------------
   ;---------------------------------------------------------------------------------------
   
   ;------------------
   ;    Quicklook
   ;------------------
   If keyword_set(RUN_QUICKLOOK) then begin
      print,'---Runing Quicklook----'
      BOREALIS_quicklook,planet,file_root,date,fmin_ql,fmax_ql,tmin_ql,N_time_ql,steps_ql,bmin_ql,$
                      bmax_ql,Smin_ql,Smax_ql,PLOT_RAW=PLOT_RAW,DATA_NORM=DATA_NORM,DATA_SLOPE=DATA_SLOPE,$
                      FULL_QUICKLOOK=FULL_QUICKLOOK,PAPERPLOTS=PAPERPLOTS,$
                      LOFAR=LOFAR,NENUFAR=NENUFAR,DYNSPECMS=DYNSPECMS,filename=filename
   
   endif 
   
  If keyword_set(run_RFI) or keyword_set(run_processing) or keyword_set(run_FFT) then begin 
   for i=start_beam,end_beam do begin ;beam for
      beam = i
      print,'-------------------------------------------------'
      print,'-------- Starting New Beam-----------------------'
      print,'-------------------------------------------------'
      print,'--------------Beam: ',strtrim(string(cgnumber_formatter(beam,decimals=0)),1),'---------------------'
      patrol_string =strtrim(string(patrol_value),1)
      lesig_string =strtrim(string(le_sig_value),1)
      
      ;-------------------------------------
      ;           Setup Log
      ;------------------------------------
       log_name = 'log_'+date+'_beam_'+strtrim(string(cgnumber_formatter(beam,decimals=0)),1)+'_pol'+strtrim(string(cgnumber_formatter(S,decimals=0)),1)+$
                   '_'+timestamp_string+'_RFI'+strtrim(string(cgnumber_formatter(run_RFI,decimals=0)),1)+'_P'+$
                   strtrim(string(cgnumber_formatter(run_processing,decimals=0)),1)+'.dat'
       If keyword_set(COMPLEX_NAME) then begin 
         log_name = 'log_'+date+'_beam_'+strtrim(string(cgnumber_formatter(beam,decimals=0)),1)+'_pol'+strtrim(string(cgnumber_formatter(S,decimals=0)),1)+'_'+$
                   strtrim(string(fix(total_time)),1)+$
                  'sec_start'+strtrim(string(fix(tmin_org)),1)+'_steps'+strtrim(string(fix(steps)),1)+'_patrol'+$
                   strtrim(string(cgnumber_formatter(patrol_value,decimals=1)),1)+$
                   '_lesig'+strtrim(string(cgnumber_formatter(le_sig_value,decimals=1)),1)+'_pex_sum_'+$
                   'RFI'+strtrim(string(cgnumber_formatter(run_RFI,decimals=0)),1)+'_P'+$
                   strtrim(string(cgnumber_formatter(run_processing,decimals=0)),1)+'_'+timestamp_string+'.dat'
       endif            
       journal,log_name
    
       print,'**********************************************************************'
       print,'******************************Inputs**********************************'
       print,'**********************************************************************'
       print,'Start Beam', start_beam
       print,'End Beam', end_beam
       print,'Beam: ',beam 
       print,'Polarization: ', S
       print,'Run RFI (1 = Yes, 0 = NO) ', run_RFI
       print,'Run Processing (1 = Yes, 0 = NO) ', run_processing
       print,'Run Post-Processing (1 = Yes, 0 = NO) ', run_postprocessing
       print,'Run FFT (1 = Yes, 0 = NO) ', run_FFT
       print,'Start time (s): ', tmin
       print,'# Spectra: ', nspectra
       print, 'Time of each slice (s)', N_time
       print, 'Total Time', total_time
       print,'# of Slices: ', steps
       print,'Start Frequency (MHz): ',fmin 
       print,'End Frequency (MHz): ', fmax
       print,'Save Mask in Bits (1 = Yes, 0 = NO):', mask_bit
       print,'Manipulate the mask in bits (1 = Yes, 0 = NO):',bit
       If keyword_Set(run_processing) then print,'Use combined mask (1 = YES, 0=NO):', dates_combine
       print,'**********************************************************************'
       print,'**********************************************************************'
       
   ;-----------------------------------------
   ;              Run Pre-Processing/RFI
   ;-----------------------------------------
       If run_RFI eq 1 then begin
          BOREALIS_RFI,tmin,N_time,steps,fmin,fmax,patrol_value,le_sig_value,beam,S,file_root,date,$
            PS=PS,VERBOSE=VERBOSE,$
            RFI=RFI,patrol=patrol,le_sig=le_sig,sum=sum,PEX=PEX,q_norm=quantile_norm,$
            save_file_name2=save_file_name2,$
            mask_save=mask_save,mask_bit=mask_bit,bit=bit,$
            full_clean_save=full_clean_save,save_datagain=save_datagain,$
            RUN_JUPITER=RUN_JUPITER,REMOVE_BADF=REMOVE_BADF,quantile=quantile,$
            Process_METHOD2=Process_METHOD2,$
            LOFAR=LOFAR,NENUFAR=NENUFAR,DYNSPECMS=DYNSPECMS,filename=filename,$
            NFREQ_PATROL=NFREQ_PATROL,NTIME_PATROL=NTIME_PATROL,NFREQ_LESIG=NFREQ_LESIG,NTIME_LESIG=NTIME_LESIG,$
            PLANET=PLANET,Full_RFI=Full_RFI
  
           print,'Save File Name for Processing/Post-Processing: ', save_file_name2
          
          ;----Email------
          email,email_address,save_file_name2+': RFI Done ',send=send_email
         
          ;-------------------
          ;update start time 
          ;-------------------
          tmin = tmin_org 
       Endif
   endfor ; end beams 
   
   ;-----------------------------------------
   ;           Create Combine Mask
   ;-----------------------------------------
    If keyword_set(create_combine_mask) then begin
       save_create,N_time,steps,patrol_value,le_sig_value,beam,RFI=RFI,patrol=patrol,le_sig=le_sig,sum=sum,pex=pex,$
          save_file_name=save_file_name,run_JUPITER=run_JUPITER,COMPLEX_NAME=COMPLEX_NAME,/NOBEAM
       
       If keyword_set(NENUFAR) then save_file_name2 = date+'_'+planet+'_pol'+strtrim(string(S),1)
       If keyword_set(LOFAR) then   save_file_name2 = date+'_pol'+strtrim(string(S),1)+'_'+save_file_name
      
      ;If S eq 4 then save_file_name2 = date+'_pol'+strtrim(string(3),1)+'_'+save_file_name
          
      If keyword_set(mask_bit) then lofar_combinemask,save_file_name2,mask_bit=mask_bit,VERBOSE=VERBOSE,PS=PS,BEAMS=COMBO_BEAMS
      ;think about 
      If not(keyword_set(mask_bit)) then lofar_combinemask,save_file_name2,mask_bit=1,VERBOSE=VERBOSE,PS=PS,BEAMS=COMBO_BEAMS
    endif 
  
    ;Save DISP
     ;  ;------------------------------------------
     ;  ;De-dispersion and rebin
     ;  ;-------------------------------------------
     disp_values    = fltarr(4)
     disp_values[0] = run_DISP                 ;Do you want to run De-dispersion; (0=NO, 1=YES)
     ;size_disp        = 1.0d                ;rebin after de-dispersion in MHz
     ;size_disp        = 3051.7578d/1e6*256. ;MHz (~0.7MHz)
     ;rebin_time_disp  = 0                   ;secs (0 = don't resample)
     disp_values[1] = size_disp          ;size of de-dispersion
     disp_values[2] = rebin_time_disp    ;if 0 then don't resample time
     disp_values[3] = dm
     disp_values_org = disp_values
  ;-----------------------------------------
  ;              Run Processing
  ;----------------------------------------- 
    i = 0 
    for i=start_beam,end_beam do begin
     If run_processing eq 1 then begin
      beam = i
       
      print,'----------------------------------'
      print,'--------Start Processing----------'
      print,'----------------------------------'
      print,'Start Beam:' , start_beam, 'End Beam: ', end_beam, 'Beam: ', i
      save_create,N_time,steps,patrol_value,le_sig_value,beam,RFI=RFI,patrol=patrol,le_sig=le_sig,sum=sum,pex=pex,$
          save_file_name=save_file_name,run_JUPITER=run_JUPITER,COMPLEX_NAME=COMPLEX_NAME
          
      save_file_name2 = date+'_pol'+strtrim(string(S),1)+'_'+save_file_name 
      If keyword_set(DynspecMS) then begin
        If S eq 4 then begin 
          Sin = 3
          save_file_name2 = date+'_pol'+strtrim(string(Sin),1)+'_'+save_file_name
        endif
      endif
      If keyword_set(NENUFAR) then save_file_name2 = date+'_'+planet+'_pol'+strtrim(string(S),1)+'_'+save_file_name
  
      ;-------------------------
      ;Update de-disp paramters 
      ;   (more comments)
      ;-------------------------
      If beam eq beam_DISP and disp_values_org[0] eq 1 then disp_values[0] = 1 Else disp_values[0] = 0
      
      ;---------------
      ;Run processing
      ;---------------
      If keyword_set(PCA) then begin
            print,'*****************************************************'
            print,'********************* PCA Sysrem*********************'
            print,'*****************************************************'
      endif       
        BOREALIS_processing,tmin,N_time,steps,fmin,fmax,beam,S,file_root,date,save_file_name2,$
             disp_values=disp_values,$
             PS=PS,VERBOSE=VERBOSE,$
             rebin_time=rebin_time,rebin_freq=rebin_freq,rebin_data=rebin_data,$
             skip_rebinloop = skip_rebinloop,$
             save_RFIMaskData = save_RFIMaskData,$
             mask_bit=mask_bit,bit=bit,$
             use_combinemask=use_combinemask,$
             TimeFreq_correction=TimeFreq_correction,$
             Apply_Corrections=Apply_Corrections,Smooth_Surface=Smooth_Surface,$
             Run_JUPITER=Run_JUPITER,REMOVE_BADF=REMOVE_BADF,N_smooth=N_smooth,$
             Process_METHOD2=Process_METHOD2,$
             LOFAR=LOFAR,NENUFAR=NENUFAR,DYNSPECMS=DYNSPECMS,$
             filename=filename,COMPLEX_NAME=COMPLEX_NAME,$
             TARGET=PLANET,FULL=Full_RFI,TO=TO,Period=Period,PCA=PCA
  ;    endif 
      
  
  ;       lofar_processing_pca,tmin,N_time,steps,fmin,fmax,beam,S,file_root,date,save_file_name2,$
  ;       disp_values=disp_values,$
  ;       PS=PS,VERBOSE=VERBOSE,$
  ;       rebin_time=rebin_time,rebin_freq=rebin_freq,rebin_data=rebin_data,$
  ;       skip_rebinloop = skip_rebinloop,$
  ;       save_RFIMaskData = save_RFIMaskData,$
  ;       mask_bit=mask_bit,bit=bit,$
  ;       use_combinemask=use_combinemask,$
  ;       TimeFreq_correction=TimeFreq_correction,$
  ;       Apply_Corrections=Apply_Corrections,Smooth_Surface=Smooth_Surface,$
  ;       Run_JUPITER=Run_JUPITER,REMOVE_BADF=REMOVE_BADF,N_smooth=N_smooth,$
  ;       Process_METHOD2=Process_METHOD2
  ;     Endif
   
       ;----Email------
       email,email_address,save_file_name2+': Processsing Done ',send=send_email
     
       tmin = tmin_org    ;restore tmin      
     Endif ;end processing
        
     ;-----------------------------------------
     ;              Run FFT
     ;-----------------------------------------
     ;Update FFT to run 
     If run_FFT eq 1 then begin
       If beam eq beam_FFT then begin 
         If not(keyword_set(run_processing)) then begin
           save_create,N_time,steps,patrol_value,le_sig_value,beam,RFI=RFI,patrol=patrol,le_sig=le_sig,sum=sum,pex=pex,$
           save_file_name=save_file_name
           save_file_name2 = date+'_pol'+strtrim(string(S),1)+'_'+save_file_name
           print, 'Save File Name for FFT:', save_file_name2
           ;If S eq 4 then save_file_name2 = date+'_pol'+strtrim(string(3),1)+'_'+save_file_name
         Endif
      
         lofar_FFT,pulsarname,dm,Period_FFT,pmin,fminp,fmaxp,ddm,save_file_name2
        
         ;----Email------
         email,email_address,save_file_name2+': FFT Done ',send=send_email
       endif ;if beam= beam FFT 
     Endif ; end fft 
   endfor ;beamloop
  endif ;if running rfi,processing,or fft
  
  ;-----------------------------------------
  ;              Run Post-Processing
  ;-----------------------------------------
   If run_postprocessing eq 1 then begin
      ;-------------------------------------
     ;           Setup Log
     ;------------------------------------
     log_name = 'log_'+date+'_pol0_'+'PostProcess_'+timestamp_string+'.dat'
     ;log_name = 'log_'+date+'_pol0_'+strtrim(string(fix(total_time)),1)+$
     ;  'sec_start'+strtrim(string(fix(tmin_org)),1)+'_steps'+strtrim(string(fix(steps)),1)+'_patrol'+$
     ;  strtrim(string(cgnumber_formatter(patrol_value,decimals=1)),1)+$
     ;  '_lesig'+strtrim(string(cgnumber_formatter(le_sig_value,decimals=1)),1)+'_pex_sum_'+$
     ;   'PostProcess_'+timestamp_string+'.dat'
     journal,log_name
  
     print,'**********************************************************************'
     print,'******************************Inputs**********************************'
     print,'**********************************************************************'
     print,'Polarization: ', S
     print,'Run Post-Processing (1 = Yes, 0 = NO) ', run_postprocessing
     print,'Start time (s): ', tmin
     print,'# Spectra: ', nspectra
     print, 'Time of each slice (s)', N_time
     print, 'Total Time', total_time
     print,'# of Slices: ', steps
     print,'Start Frequency (MHz): ',fmin
     print,'End Frequency (MHz): ', fmax
     print,'**********************************************************************'
     print,'**********************************************************************'
     
      print,'----------------------------------'
      print,'-----Start Post-Processing--------'
      print,'----------------------------------'
    
      ;---------------------------
      ;       Combine Data
      ;---------------------------  
      If keyword_set(create_dates_combine) then begin
        print,'**************'
        print, 'combine dates'
        print,'**************'
        
        readcol,filenname_combinedates,dates,dates_steps,format='(A,F)'
  
        
        create_filenames,file_root,dates,dates_steps,N_time,patrol_value,le_sig_value,S,COMPLEX_NAME=COMPLEX_NAME,$
          RFI=RFI,patrol=patrol,le_sig=le_sig,sum=sum,pex=pex,dates_filenames=dates_filenames,save_file_name=save_file_name
        
        for i=0,Combine_NBeams-1 do begin
          BOREALIS_combine,file_root,save_file_name,dates,dates_filenames,Combine_Beams[i],LIST=LIST,P=Period,TO=TO
        endfor    
       
        save_file_name3 = 'combine_'+save_file_name  
        save,LIST,save_file_name3,filename=save_file_name3+'_LIST.sav'  
        save_file_name2 = save_file_name3
       Endif
      
      If not(keyword_set(dates_combine)) and not(keyword_set(RUN_JUPITER)) then begin
           beam = 0.
           save_create,N_time,steps,patrol_value,le_sig_value,beam,RFI=RFI,patrol=patrol,le_sig=le_sig,sum=sum,pex=pex,$
              save_file_name=save_file_name,/NOBEAM,run_JUPITER=run_JUPITER
            save_file_name2 = date+'_pol'+strtrim(string(S),1)+'_'+save_file_name
            If keyword_set(DynspecMS) then begin 
               If S eq 4 or S eq 3 then begin 
                Sin = 3   
                save_file_name2 = date+'_pol'+strtrim(string(Sin),1)+'_'+save_file_name
               endif
            endif
            If keyword_set(NENUFAR) then save_file_name2 = date+'_'+planet+'_pol'+strtrim(string(S),1)+'_'+save_file_name
           ;If S eq 4 then save_file_name2 = date+'_pol'+strtrim(string(3),1)+'_'+save_file_name
     endif 
         
  ;   If keyword_set(dates_combine) then begin
  ;     LIST =  1
  ;     save_create,N_time,steps,patrol_value,le_sig_value,beam,RFI=RFI,patrol=patrol,le_sig=le_sig,sum=sum,pex=pex,$
  ;       save_file_name=save_file_name,COMPLEX_NAME=COMPLEX_NAME,/LIST
  ;     save_file_name2 = date+'_pol'+strtrim(string(S),1)+'_'+save_file_name
  ;     ;If S eq 4 then save_file_name2 = date+'_pol'+strtrim(string(3),1)+'_'+save_file_name
  ;   endif
         
     If keyword_set(RUN_JUPITER) then begin
       beam = Target_beam
       save_create,N_time,steps,patrol_value,le_sig_value,beam,RFI=RFI,patrol=patrol,le_sig=le_sig,sum=sum,pex=pex,run_JUPITER=run_JUPITER,$
         save_file_name=save_file_name,/NOBEAM
        save_file_name2 = date+'_pol'+strtrim(string(S),1)+'_'+save_file_name
     Endif
     
     print, save_file_name2
     print,'Target Beam: ', Target_Beam
     print,'Sky Beam: ', Sky_Beam
     
      ;------------------------
      ;Threshold of observables (Q3-Q7)
      ;-------------------------
      If keyword_Set(Default_Threshold) then begin
        threshold = dindgen(50,start=1,increment=0.1)  ; 1 - 5.9 sigma
      endif
      
      ;-----------------------
      ;Thresholds of the mask (loop)
      ;----------------------
      If keyword_set(Default_MASK) then begin
        threshold_t = [0.10] ;threshold = 1 - x
        threshold_f = [0.10] ;0.05 = 95%
        If keyword_set(Run_Mask) then begin 
          threshold_t = [threshold_t,threshold_mask]
          threshold_f = [threshold_f,threshold_mask]
        endif 
      endif
  
      If not(keyword_set(Default_MASK)) then begin
        threshold_t = [1.00] ;threshold = 1 - x
        threshold_f = [1.00] ;0.05 = 95%
        If keyword_set(Run_Mask) then begin 
          threshold_t = [threshold_t,threshold_mask]
          threshold_f = [threshold_f,threshold_mask]
        endif 
      endif
      
      ;------------------------
      ;     Freq ranges 
      ;------------------------
      ;(10 MHz ranges)
      If keyword_set(Default_FreqRange) and keyword_set(Freq_10MHz) then begin
        ssize = 10.0 ;10MHz step size
        ss = floor((fmax - fmin)/ssize)    ;steps
        for i=0,ss do begin
         If i eq 0 then begin
          fmins = fmin
          fmaxs = fmax
         endif
        If i gt 0 then begin
         fmaxn  = fmax - ssize*(i-1)
         fminn  = fmaxn - ssize
         fmins   = [fmins,fminn]
         fmaxs   = [fmaxs,fmaxn]
        endif
        endfor
      endif
      
      ;Half Freq range
      If keyword_set(Default_FreqRange) then begin
          fhalf = (fmax - fmin)/2.0
          fmins = [fmin,fmin,fmax-fhalf]
          fmaxs = [fmax,fmin+fhalf,fmax]
          If keyword_set(Run_FreqRange) then begin 
            fmins = [freq_minextra1,fmins]
            fmaxs = [freq_maxextra2,fmaxs]
          endif 
      endif
      
      If not(keyword_set(Default_FreqRange)) then begin
        fmins = fmin
        fmaxs = fmax 
        If keyword_set(Run_FreqRange) then begin 
          fmins = [freq_minextra1]
          fmaxs = [freq_maxextra2]
        endif 
      endif  
       
        If keyword_set(save_gauss) then read_gauss = 0      ;save Gauss
        If not(keyword_set(save_gauss)) then read_gauss = 1 ;read Gauss
    ;Run with default QTmin and QTmin    
        BOREALIS_postprocessing,save_file_name2,Date,Target_Beam,Sky_Beam,S,rebin_time,rebin_time2,rebin_freq,rebin_freq2,fmaxs,fmins,threshold,$
          threshold_p=threshold_p,threshold_t=threshold_t,threshold_f=threshold_f,max_threshold=max_threshold,$
          Extra_Beam1=Extra_Beam1,Extra_Beam2=Extra_Beam2,$
          LIST=LIST,Planet_Name=Planet_Name,$
          TI=TI,tau_Qplot=tau_Qplot,slope=slope,Qwin=Qwin,$
          Q1=Q1,Q2=Q2,Q3=Q3,Q4=Q4,save_Q=save_Q,TARGET_ONLY=TARGET_ONLY,$
          ELLCORR=ELLCORR,$
          read_gauss=read_gauss,gauss_steps=gauss_steps,file_g=root_gauss,$
          OnBeam_Label=OnBeam_Label,OffBeam_Label=OffBeam_Label,$
          Default_Time=Default_Time,QTmin=QTmin,QTmax=QTmax,$
          PERIOD=Period,T0=T0,$
          VERBOSE=VERBOSE,PS=PS,$
          RUN_JUPITER=RUN_JUPITER,$
          REMOVE=REMOVE,PAPERPLOTS=PAPERPLOTS,MINT=MINT,MAXT=MAXT,$
          DYNSPECMS=DYNSPECMS,XTMIN=XTMIN,XTMAX=XTMAX
   
  ;      ;For loop over time
  ;      for t=0,2 do begin
  ;        XTMIN = XTMIN*3600.0d ; secs 
  ;        XTMAX = XTMAX*3600.0d ; secs
  ;        Default_Time= 0 
  ;        ;T=0 ; 1st third
  ;        ;T=1 ; 2nd third
  ;        ;T=2 ; last third  
  ;        xt_step = (XTMAX - XTMIN)/3.0d
  ;        If t eq 0 then begin 
  ;          QTmin = XTMIN   
  ;          QTmax =  XTMIN + xt_step
  ;        endif
  ;         If t eq 1 then begin
  ;          QTmin = XTMIN + xt_step
  ;          QTmax =  XTMIN + xt_step*2.0d
  ;        endif 
  ;        If t eq 2 then begin
  ;          QTmin = XTMIN + xt_step*2.0d
  ;          QTmax =  XTMAX
  ;        endif
  ;        BOREALIS_postprocessing,save_file_name2,Date,Target_Beam,Sky_Beam,S,rebin_time,rebin_time2,rebin_freq,rebin_freq2,fmaxs,fmins,threshold,$
  ;          threshold_p=threshold_p,threshold_t=threshold_t,threshold_f=threshold_f,max_threshold=max_threshold,$
  ;          Extra_Beam1=Extra_Beam1,Extra_Beam2=Extra_Beam2,$
  ;          LIST=LIST,Planet_Name=Planet_Name,$
  ;          TI=TI,tau_Qplot=tau_Qplot,slope=slope,Qwin=Qwin,$
  ;          Q1=Q1,Q2=Q2,Q3=Q3,Q4=Q4,save_Q=save_Q,TARGET_ONLY=TARGET_ONLY,$
  ;          ELLCORR=ELLCORR,$
  ;          read_gauss=read_gauss,gauss_steps=gauss_steps,file_g=root_gauss,$
  ;          OnBeam_Label=OnBeam_Label,OffBeam_Label=OffBeam_Label,$
  ;          Default_Time=Default_Time,QTmin=QTmin,QTmax=QTmax,$
  ;          PERIOD=Period,T0=T0,$
  ;          VERBOSE=VERBOSE,PS=PS,$
  ;          RUN_JUPITER=RUN_JUPITER,$
  ;          REMOVE=REMOVE,PAPERPLOTS=PAPERPLOTS,MINT=MINT,MAXT=MAXT,$
  ;          DYNSPECMS=DYNSPECMS
  ;      endfor
        endif ;end postprocessing
        
    binTime = BIN_DATE(systime(/ut))  ;
    timestamp_string = TIMESTAMP(YEAR = binTime[0],MONTH = binTime[1], $
       DAY = binTime[2], HOUR = binTime[3], MINUTE = binTime[4], SECOND =binTime[5])  ;2012-09-04T11:25:15Z
  
   ps_pdf,/REMOVE  ;convert all ps to pdf
   print,'-----------------------------------------------------------------------'
   print,'SYSTIME entire program = ',SYSTIME(1) - st, ' sec'
   print,'Program Done at ', timestamp_string
   print,'-----------------------------------------------------------------------'
   print, 'End of Main program: Goodbye'
      
   ;----Email------
   email,email_address,'Program Done ',send=send_email
   journal ;close log
   If S eq 0 then i = 0 
 endfor ; End loop for polarizations  
end ;end lofar pipeline
