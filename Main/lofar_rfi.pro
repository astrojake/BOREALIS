;***************************************************************
;      RFI on lofar data for exoplanet search
;**************************************************************
;***************************************************************
;Author: Jake Turner (UVA)
;        J.M. Griessmeier (CNRS) and Philippe Zarka (Obspm)
;        
;Inputs:
;      tmin     ;Starting time (seconds)
;      N_time         ;Size of steps in seconds
;      steps          ;Number of steps to take in N_time
;      fmin     ;Starting Frequency (MHZ)
;      fmax       ;Ending Frequency (MHZ)
;      Beam           ;Beam (Beam 0 - Beam 4; 55cnc, pulsar, sky, bright source)
;      S              ;Polarization
;      patrol_value   ;# of STDEV above which a value is considered as a spike
;      le_sig_value   ;# of STDEV above which a value is considered as RFI (above and below)
;      file_root      ;root directory for raw data (eg. file_root = '/data/jake.turner/exoplanet/LC5_DDT_002/')
;      date           ;Dates in an array (eg. date = 'L429868'; folder must exist in root directory)
;      
;Keywords:      
;      RFI             ;Run RFI 
;      patrol          ;Run patrol RFI
;      le_Sig          ;Run Leg Sig RFI
;      sum             ;Run Sum RFI
;      hex             ;Run hex RFI
;      allsave         ;Save all of the variables (makes very large file but good for debugging)
;      full_clean_save ;Save the full non-binned cleaned data (only useful for debugging)
;      PS              ;Make ps file of plots 
;      save_datagain   ;Save each time-frequency curve for each time step /save_tfall
;      save_file_name   ;save file name for processing and fft 
;      mask_save        ;Save the mask created from the RFI mitigation (called 'mask_full_'+save_file_name+'.sav')
;      mask_bit          ; Read the mask from the save file in bits (then it is converted to bytes)
;      bit               ; Manipulate the mask in bits
;      Run_JUPITER     ; Run Jupiter analysis
;      quantile        ; What quantile to use 
;      Process_METHOD2 ;Run method 2 in processing
;      
;Uses: -hdf5 programs
;      -Programs from Phillip (rfi_mitigate and subprograms)
;        - rfi_mitigate.pro: has been updated
;          -(Programs within)
;      -cgcoyate programs
;      -idlastro lib
;      
;Important Notes:
;          1.) Removed 0 and 64 element of each subband; this is a flag
;          2.) Normalized around 0, so that bad pixels are set as 0 in fft and other steps
;          3.) Only one free parameter not input into code (quantile)
;          4.) Free parameters in rfi_mitigate.pro are hard coded (I can change this)
;          5.) The Jupiter emission tests are hardcoded and not general 
;      
;File Format:
;  L248553_SAP000_B000_S0_P000_bf.h5'
;    L is the date
;    B000 - B003: 4 beams: (exoplanet, pulsar, blanck sky, bright source)
;    S0-S3:  4 polizations: IQUV
;
;To Run: 
; Put rfi_mitigate.pro and subprograms in path
; compile lofar_rfi.pro
; e.g. lofar_RFI,tmin,N_time,steps,fmin,fmax,patrol_value,le_sig_value,Beam,S,file_root,date,$
;       /RFI,/patrol,/le_sig,/sum,/pex,/quantile_norm,save_file_name2=save_file_name,mask_save=mask_save

pro lofar_RFI,tmin,N_time,steps,fmin,fmax,patrol_value,le_sig_value,Beam,S,file_root,date,$
              RFI=RFI,patrol=patrol,le_sig=le_sig,sum=sum,PEX=PEX,$
              PS=PS,VERBOSE=VERBOSE,$
              allsave=allsave,full_clean_save=full_clean_save,$
              save_file_name2=save_file_name2,mask_save=mask_save,save_datagain=save_datagain,$
              avg_norm=avg_norm,q_norm=quantile_norm,$
              mask_bit=mask_bit,bit=bit,$
              Run_JUPITER=Run_JUPITER,REMOVE_BADF=REMOVE_BADF,$
              quantile=quantile,$
              Process_METHOD2=Process_METHOD2,$
              LOFAR=LOFAR,NENUFAR=NENUFAR,DYNSPECBF=DYNSPECBF,filename=filename

;start time
st= SYSTIME(1)
startmem_RFI =  MEMORY(/CURRENT)
paper_plot = 0

;Set Polarization
If S ne 0 then begin  
 quantile_norm = 0  ;Only use quantile for Stokes-I 
 avg_norm = 1       ;Use average background for other polarizations
 print,'Polarzation Selection: Normalize by Background'
endif

;---------------------------
;----------------------------
;        Check Flags
;----------------------------
;---------------------------

If keyword_set(VERBOSE) then begin
 print,'-------------------------------------'
 print,'---        Checking Flags         ---'
 print,'-------------------------------------'
endif 

If steps eq 1 then disp = 0  ;skip de-disp if steps = 1
If steps eq 1 then print, 'Steps equal 1, De-Dispersion not set'
If keyword_set(VERBOSE) then begin
  if keyword_set(Patrol) ne 0 then $
       print, 'Patrol: set' else $
       print, 'Patrol: not set'
  
  if keyword_set(sum) ne 0 then $
       print, 'Sum: set' else $
       print, 'Sum: not set'
  
  if keyword_set(le_sig) ne 0 then $
       print, 'Le Sig: set' else $
       print, 'Le Sig: not set'  
  
  if keyword_set(pex) ne 0 then $
       print, 'Pex: set' else $
       print, 'Pex: not set'    
  if keyword_set(Run_Jupiter)  then $
    print, 'Jupiter: set' else $
    print, 'Jupiter: not set'
  if keyword_set(bit) ne 0 then $
    print, 'Bit: set' else $
    print, 'Bit: not set'  
  If keyword_set(bit) then begin
    mask_bit = 1 
    print, 'Mask Bit Flag overwritten since bit = 1'
  Endif
  If not(keyword_set(mask_bit)) then begin
    bit = 0 
    print, 'Manipulate Bit Flag (bit) overwritten since mask_bit = 0'
  Endif
  if keyword_set(mask_bit) ne 0 then $
    print, 'mask_bit: set' else $
    print, 'mask_bit: not set'
  print, 'Post Process Method 2 (YES = 1; NO = 0)', Process_METHOD2  
  print,'---------------------------------------'
  print,'---------------------------------------'         
endif 
;------------------------
;-----------------------
;       Inputs
;-----------------------
;;-----------------------
quantile_string = +strtrim(string(cgnumber_formatter(quantile,decimals=2)),1)
If keyword_set(VERBOSE) then begin
  print,'Quantile: ', quantile
endif 

;--------------------------------------------------------
;;              Time
;;--------------------------------------------------------
;total_time = steps*N_time  ;max time for RFI cleaning

;*********************
;Create save_file_name 
;*********************
save_create,N_time,steps,patrol_value,le_sig_value,beam,RFI=RFI,patrol=patrol,le_sig=le_sig,sum=sum,pex=pex,run_JUPITER=run_JUPITER,$
  save_file_name=save_file_name

If keyword_set(VERBOSE) then begin
  print,'File name for Save: ', save_file_name
endif 

 ;Restore Cleaned Data and mask if not running RFI
If not(keyword_set(RFI)) then begin
  If keyword_set(VERBOSE) then begin
    print,'******************************'
    print,'Restoring Cleaned Data for RFI'
    print,'******************************'
    print, 'Restore name: ', 'dataclean_full'+save_file_name+'.sav'
  endif      
  restore,file='dataclean_full_'+save_file_name+'.sav
  ;restores: data_clean and p2_full
endif

    ;convert Beam and Pol values into strings 
    B_string = strtrim(string(beam),1)
    S_string2 = strtrim(string(S),1)        ;save filename 
      If S lt 4 then S_in = S   ; I, Q, U, abs(V)
      If S eq 4 then S_in = 3   ; V' (method 1)
      If S eq 5 then S_in = 3   ; V'
      If S eq 6 then begin     ; L = sqrt(Q^2 + U^2)  
         S_in = 1 
         Sin1 = 1              ;Q
         Sin2 = 2              ;V
      endif 
      If S lt 4 then S_string = strtrim(string(uint(S)),1)   ; I, Q, U, V   ;read file name
      If S eq 4 then S_string = strtrim(string(uint(3)),1)   ; V'
      If S eq 5 then S_string = strtrim(string(uint(3)),1)   ; V'
      If S eq 6 then begin     ; L = sqrt(Q^2 + U^2)  
         S_string = strtrim(string(uint(1)),1); S= 1 = Q
      endif 

    ;----------------------------------------------------------
    ;              Setup Important Info for Data
    ;----------------------------------------------------------
    ;file to read in headers and arrays
    If keyword_set(LOFAR) then begin 
      file = file_root+date+'/'+date+'_SAP000_B00'+B_string+'_S'+S_string+'_P000_bf'
      print, 'File: ',file
  
      READ_H5_DATA, file,0,1,fmin,fmax,Beam,S_in,x,t,xf,n_t,nf,Time=tf,/NSPECTRA,$
        NOF_SAMPLES=NOF_SAMPLES,START_MJD=START_MJD,MJD_sampling_time=MJD_sampling_time,$
           ANTENNA_SET=ANTENNA_SET,SAMPLING_TIME=SAMPLING_TIME,NOF_stations=NOF_stations,$
           channel_width=channel_width,REMOVE_BADF=REMOVE_BADF  
    endif 
    
    If keyword_set(NENUFAR) then begin
      Beam_in = Beam
      If S eq 0 then Stokes = 0 ;Stokes-I
      If S eq 1 then Stokes = 1 ;Stokes-Q
      If S eq 2 then Stokes = 2 ;Stokes-U
      If S eq 3 then Stokes = 3 ;Stokes-V
      If S eq 4 then Stokes = 4 ;Stokes-V Prime
      If S eq 0 then nstokes = 1 ELSE nstokes = 4
      read_nenufar, file_root+filename, temporary(data2),tt,ff,temporary(Beam_in),ndata, ntt,dt,nf,df,ns,Stokes,$
        nstokes=nstokes,tmin=0,tmax=1,/VERBOSE,CHAN=CHAN,BEAMLETS=BEAMLETS,/REMOVE_BADF,/CORRECT_BANDPASS,/DATA_BANDPASS
      sampling_time = dt
      channel_width = df
    Endif
    If keyword_set(DYNSPECBF) then begin
      If S eq 0 then POL = 1 ;Stokes-I
      If S eq 1 then POL = 2 ;Stokes-Q
      If S eq 2 then POL = 3 ;Stokes-U
      If S eq 3 then POL = 4 ;Stokes-V
      If S eq 4 then POL = 4 ;Stokes-V prime
      If S eq 5 then POL = 4 ;Stokes-V prime
      print,'File to read',file_root+filename
      read_dynspec_fits,file_root+filename+'_beam'+B_string,0,1,fmin,fmax,POL,temporary(data2),xt,xf,nt,nf,dt=dt,df=df
      sampling_time = dt
      channel_width = df
    Endif
         
    If keyword_set(PS) then begin   
        print,'----Making Plot----'
        set_plot,'PS'
        filename_ps = date+'_pol'+S_string2+'_'+save_file_name
        save_file_name = filename_ps
        print,'Filename Ps:', filename_ps
        device,filename=filename_ps+'_RFI.ps',/landscape
    endif 
   
   ;----------------------------------------------
   ;----------------------------------------------
   ;               Setup Time Array
   ;----------------------------------------------
   ;----------------------------------------------
   tmin_or = tmin
   If keyword_set(VERBOSE) then begin
   ; print, 'Max Time (s)' ,max(tf)
   endif 
  
  ;Setup Jupiter
  If keyword_set(run_JUPITER) then begin
     input_file_J = 'Input_Jupiter.dat'
     read_jupiter_input,input_file_J,tmin_Jup=tmin_J_main,tmax_Jup=tmax_J_main
     lofar_jupitersetup,input_file_J,data,xf,xt,S_in,tmin_J,tmax_J,data_new=data_new,$
       factor=factor,data_Jupiter=Jupiter_Norm,p2_Jupiter=p2_Jupiter,t_J=t_J,freq_J=freq_J,$
       filename_p2=filename_p2, tf_J=tf_J, nf_J=nf_J,fmin_Jup=fmin_Jup,fmax_Jup=fmax_Jup,target_beam_J=target_beam_J,$
       A_eff_L=A_eff_L, JA_eff_L=A_eff_L_Jup, Jupiter_background=Jupiter_background, scale=scale,LOFAR_RUN=LOFAR_RUN,$
       file_Jup=file_Jup,f_nearest=f_nearest,time_p=time_p,freq_p=freq_p,method2_Jup=method2_Jup,JBackV=Jupiter_background_V,JBackAll=Jupiter_Back_All,VSkyBeam=VSky_Prime
  endif ;endJupiter   

  ;---------------------------------
  ;          I Mask for all RFI
  ;Read I mask for all polarizations
  ;-------------------------------- 
  If S ge 1 then begin
    save_create,N_time,steps,patrol_value,le_sig_value,beam,RFI=RFI,patrol=patrol,le_sig=le_sig,sum=sum,pex=pex,run_JUPITER=run_JUPITER,save_file_name=save_file_name_I  ;All RFI  
    save_file_name3 = date+'_pol'+strtrim(string(0),1)+'_'+save_file_name_I
    If keyword_set(mask_bit) and keyword_set(bit)  then begin
     print,'Restore I Mask: ','mask_full_bit_steps'+save_file_name3+'.sav'
     restore,filename='mask_full_bit_steps_'+save_file_name3+'.sav'   
       ;xt,xf,xsize_step,p2_full_bit_p
      p2_full_bit_p_I = temporary(p2_full_bit_p)
      xsize_step_I    = temporary(xsize_step)
    endif
    
     If keyword_set(mask_bit) and not(keyword_set(bit)) then begin
      print,'Restore mask: ', 'mask_full_bit_'+save_file_name3+'.sav'
      restore,filename='mask_full_bit_'+save_file_name3+'.sav' ;flag
        ;xt,xf,xsize_step,p2_full_bit
      BITARRAY_TO_BYTEARRAY,p2_full_bit,xsize,p2_full
      p2_full_I =temporary(p2_full) 
     Endif
   
  endif 
   
  ;Setup Bit mask
  If keyword_set(mask_bit) and keyword_set(bit) then begin
      p2_full_bit_p =  ptrarr(steps)  ;pointer for mask
  Endif ;mask bit
 
  ;-------------------------------
  ;Setup Method 2 Post-Processing 
  ;-------------------------------
  If keyword_set(Process_METHOD2) then begin 
     freq_response  = dblarr(steps,nf)
     xt_rebin       = dblarr(steps)
     pf_full      = dblarr(steps,nf)  ;frequency mask for every step  
  endif
 
  ;-------------------------------------------------      
  ;--------------------------------------------------
  ;                 Time Loop
  ;-------------------------------------------------- 
  ;--------------------------------------------------
  for kl=0,steps-1 do begin  ;begin time loop
    tmax =  tmin + N_time  ;set tmax
    print,'-----------------------------------------------------------------------'
    print,'reading slice ',kl,' start time: ', tmin, ' end time: ', tmax
    print,'-----------------------------------------------------------------------'
    If kl eq 0 then VERB = 1 Else VERB = 0
      
    ;----------------------------------------------------------
    ;        Read Data for specific Beam and Polarization
    ;----------------------------------------------------------
    If keyword_set(LOFAR) then begin 
      file = file_root+date+'/'+date+'_SAP000_B00'+B_string+$
        '_S'+S_string+'_P000_bf'                                  ;file name
      If keyword_set(VERBOSE) then begin
        print,'File: ', file
      endif 
      READ_H5_DATA, file,tmin,tmax,fmin,fmax,Beam,S_in,data,time,xf,nt,nf,VERB=VERB,REMOVE_BADF=REMOVE_BADF      
    endif 
    
    If keyword_set(NENUFAR) then begin
      beam_in = beam
      If S eq 0 then Stokes = 0 ;Stokes-I
      If S eq 1 then Stokes = 1 ;Stokes-Q
      If S eq 2 then Stokes = 2 ;Stokes-U
      If S eq 3 then Stokes = 3 ;Stokes-V
      If S eq 4 then Stokes = 4 ;Stokes-V Prime
      If S eq 0 then nstokes = 1 ELSE nstokes = 4
      read_nenufar, file_root+filename,data,time,xf,temporary(beam_in),ndata,nt,dt,nf,df,ns,Stokes,nstokes=nstokes,tmin=tmin,tmax=tmax,fmin=fmin,fmax=fmax,VERBOSE=VERB,$
        REMOVE_BADF=REMOVE_BADF,/CORRECT_BANDPASS,/DATA_BANDPASS
      
      beam_in = beam
      read_nenufar, file_root+filename,temporary(d),temporary(tt),temporary(ff),temporary(beam_in),ndata,nt,dt,nf,df,ns,Stokes,nstokes=nstokes,tmin=tmin,tmax=tmax,fmin=fmin,fmax=fmax,VERBOSE=VERB,$
        REMOVE_BADF=REMOVE_BADF,/CORRECT_BANDPASS,/DATA_BANDPASS,RFI=1,mask=mask
      If keyword_set(VERBOSE) then begin
        print,'File: ' ,file_root+filename
      endif
    endif
    If keyword_set(DYNSPECBF) then begin
      If S eq 0 then POL = 1 ;Stokes-I
      If S eq 1 then POL = 2 ;Stokes-Q
      If S eq 2 then POL = 3 ;Stokes-U
      If S eq 3 then POL = 4 ;Stokes-V
      If S eq 4 then POL = 4 ;Stokes-V^2
      If S eq 5 then POL = 4 ;Stokes-V prime
      read_dynspec_fits,file_root+filename+'_beam'+B_string,tmin,tmax,fmin,fmax,POL,data,time,xf,nt,nf
    Endif
    
    ;*********************
    ;   Save Time Series
    ;*********************
    If kl eq 0 then begin
      xt = DINDGEN(steps*nt,INCREMENT=sampling_time, start = tmin)
    endif
    xt[nt*kl:(kl+1)*nt-1] = time
    
    If keyword_set(Process_METHOD2) then begin 
      reduce_array,time,[nt],temp_t
      xt_rebin[kl] = temporary(temp_t)
    endif   
    
        ;******Polarizations *******
        If S eq 3 then data = abs(data)  ;V = abs(V)
        If S eq 4 or S eq 5 then begin ;V'
          print,' Running Vprime'
          Print,'-------------------'
          Print,'Reading I Data'
          Print,'-------------------'
          If keyword_set(LOFAR) then begin
            S_I_string = strtrim(string(0),1)
            file_I = file_root+date+'/'+date+'_SAP000_B00'+B_string+$
              '_S'+S_I_string+'_P000_bf'
            print,'Tmin, Tmax, Fmin, Fmax: ',tmin, tmax,fmin,fmax
            READ_H5_DATA, file_I,tmin,tmax,fmin,fmax,Beam,0,data_I,t,f,nt,nf,REMOVE_BADF=REMOVE_BADF,VERB=VERB   ;read Stokes-I data
          endif 
          If keyword_set(NENUFAR) then begin
            beam_in = beam
            Stokes = 0
            nstokes_I = 1
            read_nenufar, file_root+filename,data_I,time,xf,temporary(beam_in),ndata,nt,dt,nf,df,ns,Stokes,nstokes=nstokes_I,tmin=tmin,tmax=tmax,fmin=fmin,fmax=fmax,VERBOSE=VERB,$
              REMOVE_BADF=REMOVE_BADF,/CORRECT_BANDPASS,/DATA_BANDPASS,RFI=1,mask=mask_I
              mask = mask*mask_I
          endif
          If kl eq 0 then STANDARD_PLOTS, data,p2,time,xf,'Raw V',xunit = 'secs',panel_ps=1
          If kl eq 0 then STANDARD_PLOTS, data_I,p2,time,xf,'Raw I',xunit = 'secs',panel_ps=1
          data = data/data_I   ;V/I
          ;If kl eq 0 then STANDARD_PLOTS, data,p2,time,xf,'V/I',xunit = 'secs',panel_ps=1
        endif
        If S eq 6 then begin ; L
          print,' Running L'
          data_Q = temporary(data)
           S_U_string = strtrim(string(2),1)
          file_U = file_root+date+'/'+date+'_SAP000_B00'+B_string+$
            '_S'+S_U_string+'_P000_bf'
          READ_H5_DATA, file_U,tmin,tmax,fmin,fmax,Beam,2,data_U,t,f,nt,nf,REMOVE_BADF=REMOVE_BADF 
          data = sqrt(data_Q^2.0d + data_U^2.0d)
        endif 
 
    If keyword_set(Process_METHOD2) and not(keyword_set(method2_Jup)) then begin
      print,'**********************************************************************'
      print,'*******Method 2 Processing: 10% quantile of data is used*********'
      print,'**********************************************************************'
      data_raw = data
    Endif
              
    ;*************************************
    ;         Add: Jupiter 
    ;*************************************
    If S eq 0 and keyword_set(run_JUPITER) then begin
      print,'********** Start Jupiter *********************'
      If kl eq 0 then tmin_J = tmin_J_main
      tmax_J =  tmin_J + N_time  ;set end time
      print,'Min and max',tmin_J,tmax_J
      
	;      If tmin le 6910. then begin          ;flag!!!!
	;        tmin_J = tmin_J_main
	;        tmax_J = tmin_J + N_time
	;        print,'Min and Max in skip loop', tmin_J,tmax_J
	;        print,'******Skipping Jupiter******'
	;        goto, SKIPJUPITER
	;      endif 
      
      If tmax_J ge tmax_J_main then begin
        print,'******Skipping Jupiter******'
        goto, SKIPJUPITER
      Endif
      lofar_jupiter,input_file_J,data,xf,xt,S,tmin_J,tmax_J,data_new=data_new,$
                  factor=factor,data_Jupiter=data_Jupiter,p2_Jupiter=p2_Jupiter,t_J=t_J,freq_J=freq_J,$
                  filename_p2=filename_p2, tf_J=tf_J, nf_J=nf_J,fmin_Jup=fmin_Jup,fmax_Jup=fmax_Jup,target_beam_J=target_beam_J,$
                  A_eff_L=A_eff_L, JA_eff_L=A_eff_L_Jup, Jupiter_background=Jupiter_background,scale=scale,LOFAR_RUN=LOFAR_RUN,$
                  file_Jup=file_Jup,f_nearest=f_nearest,time_p=time_p,freq_p=freq_p,$
                  UJupiter=data_Jupiter_U,QJupiter=data_Jupiter_Q      
        If kl eq 0 then begin
          STANDARD_PLOTS, data_Jupiter,data_Jupiter,t_J,freq_J,'Jupiter Norm',xunit = 'secs',panel_ps=1
          STANDARD_PLOTS, factor,factor,t_J,freq_J,'Jupiter Factor',panel_ps=1
        Endif
       
      data = temporary(data_new)
      SKIPJUPITER:
      print,'********** End Jupiter *********************'  
    Endif  ;end Jupiter
   
   If keyword_set(Process_METHOD2) and keyword_set(method2_Jup) then begin
    print,'**********************************************************************'
    print,'*******Method 2: 10% quantile of data WITH Jupiter is used*********'
    print,'**********************************************************************'
    data_raw = data
   Endif
    
    If kl eq 0 then begin
     If keyword_set(RFI) then begin
       If keyword_set(VERBOSE) then begin
         print, 'Setup data clean arrays'
       endif   
       
       If keyword_set(full_clean_save) then begin
         data_clean     = dblarr(long(steps*nt),nf)
       endif
       
       ;If keyword_set(save_datagain) then begin
        If S ne 5 then data_gain_save = dblarr(steps,nf)
        If S eq 5 or S eq 4 then begin
          data_mean_save = dblarr(steps,nf)
          data_mean2_save = dblarr(steps,nf)
          data_stdev_save = dblarr(steps,nf)
        Endif
      ; endif  
       
       If not(keyword_set(bit)) then begin
         p2_full        = bytarr(long(steps*nt),nf) +1b
       endif   
     Endif
    endif ;kl eq 0
               
    ;-------------------
    ;Create Mask
    ;-------------------
     p2=bytarr(nt,nf)+1b
     
     ;Update Mask from Nenufar
     If keyword_set(Nenufar) then begin
       If keyword_set(PEX) then begin 
         If not(keyword_set(PMIN_F)) then PMIN_F = 3
         If not(keyword_set(EXPF)) then EXPF = 5
         If not(keyword_set(PMIN_T)) then PMIN_T = 3
         If not(keyword_set(EXPT)) then EXPT = 5
         If not(keyword_set(PMIN_F2)) then PMIN_F2 = 5
         If not(keyword_set(EXPF2)) then EXPF2 = 12
         If not(keyword_set(PMIN_T2)) then PMIN_T2 = 5
         If not(keyword_set(EXPT2)) then EXPT2 = 12
       tex_frex, mask, [EXPF,PMIN_F], [EXPT,PMIN_T], /verbose
       tex_frex, mask, [EXPF2,PMIN_F2], [EXPT2,PMIN_T2], /verbose
      endif 
      
       w_35MhZ = Where(xf ge 35.3 and xf le 36.0,count)  ;Get rid of the
       If count gt 0 then p2[*,w_35MhZ] = 0b
       ;Combine Masks
       p2 = p2*mask
     endif 
     ;---------------------------------
     ; For polarization; read in I mask
     ;  Create Byte Mask from bit mask
     ;----------------------------------
      If S ne 0 and keyword_set(mask_bit) and keyword_set(bit) then begin
       pointer_to_array,p2_full_bit_p_I[kl],out_array=p2_bit_I ;convert step into bit array
       BITARRAY_TO_BYTEARRAY, p2_bit_I, xsize_step_I,p2_I      ;convert bit array to byte array
       If S gt 0 then p2 = p2*p2_I  ;update p2
      endif 

      If S ne 0 and keyword_set(mask_bit) and not(keyword_set(bit)) then begin
       p2_I = p2_full_I[nt*kl:(kl+1)*nt-1,*]
       p2 = p2*p2_I
      endif 

    ;---------------------------
    ;    Make Plot of Raw Data
    ;--------------------------
      If keyword_set(PS) and kl eq 0 then begin
         label   = 'Raw Data: S'+strtrim(string(uint(S_In)),1) 
         data = data
         If S eq 4 or S eq 5 then label='v=V/I'
         STANDARD_PLOTS, data,p2,time,xf,label,xunit = 'secs',panel_ps=1
         reduce_array,data,[nt,1],data_f2
      endif
      
     ;-----------------------------------------------
     ;         Normalize the Data: V Data
     ;-----------------------------------------------
     If S eq 5 or S eq 4 then begin      
       ;find data_mean of background
       MAKE_BACKGROUND,data*p2,'', data_mean,ss,nn
       LE_AUTO_S,data_mean,101,3.5,0,data_mean,pnet
       data_mean = smooth(data_mean,3,/EDGE_TRUNCATE)   
       If keyword_set(PS) and kl eq 0 then begin 
         cgplot,xf,data_mean,/xstyle,/ystyle,title='v=V/I mean ',xtitle='Freq (MHz)'
       endif 
       
       ;vprime  
       data = (data - rebin(reform(data_mean,1,nf),nt,nf))                        ;vprime = v - <v> 
       If kl eq 0 then STANDARD_PLOTS, data*p2,p2,time,xf,'vprime=(v-<v>)',xunit = 'secs',panel_ps=1
         
       ;Make I background and I data
 ;      MAKE_BACKGROUND,data_I*p2,'', data_mean_I,ss,nn
       data_gain = fltarr(nf) 
       for iii=0,nf-1 do data_gain[iii] = dyn_n(data_I(*,iii),quantile) 
       LE_AUTO_S,data_gain,101,patrol_value,0,data_gain,pnet 
       data_I = temporary(data_I)/rebin(reform(data_gain,1,nf),nt,nf)
       If keyword_set(PS) and kl eq 0 then begin
        STANDARD_PLOTS, data_I*p2_I,p2,time,xf,'Stokes-I Normalized',xunit = 'secs',panel_ps=1
       endif
       
       ;Vprime 
       data = data*data_I*p2    ;V(f) = vprime*I*mask
    endif 
      
    ;-----------------------------------------------
    ;         Normalize the Data: I Data
    ;-----------------------------------------------
    If S ne 4 then begin
       ;--------------------------
       ;Integrated Spectrum Curve; Background
       ;--------------------------
       If keyword_set(avg_norm) then begin
         ;data_avg = rebin(data,1,nf)
         MAKE_BACKGROUND,data,'',data_avg,ss,nn 
         plot,xf,data_avg,/xsty,xtit='Frequency (MHz)',ytit='Intensity',tit='Normalization Curve: Integrated Spectrum',$
          yrange=[min(rebin(data,1,nf)),2e14],/ystyle
       endif 
        
       ;-----------------------------------
       ;        10% Quantile
       ;-----------------------------------
       If keyword_set(quantile_norm) then begin
         data_gain = fltarr(nf) 
         for iii=0,nf-1 do data_gain[iii] = dyn_n(data(*,iii),quantile) 
        
          If keyword_set(PS) and kl eq 0 then begin
           plot,xf,data_gain,/xsty,xtit='Frequency (MHz)',ytit='Intensity',$
            tit= strtrim(string(cgnumber_formatter(quantile*100.,decimals=2)),1)+'% Quantile Normalization Curve',$
            yrange=[min(data_gain),max(data_gain)],/ystyle;,xrange=[50,70]
         endif
       endif 
    
       ;----------------------------------------
       ;Clean Normalization Curve of RFI Spikes
       ;----------------------------------------
       If keyword_set(avg_norm) then begin
         LE_AUTO_S,data_avg,101,patrol_value,0,data_avg,pnet
       endif 
        
       If keyword_set(quantile_norm) then begin
        LE_AUTO_S,data_gain,101,patrol_value,0,data_gain,pnet 
         If keyword_set(save_datagain) then begin 
          data_gain_save(kl,*) = data_gain
         endif    
       endif 
       
      ;----------
      ;Make Plot 
      ;---------- 
       If keyword_set(quantile_norm) then begin 
          If keyword_set(PS) and kl eq 0 then begin
            plot,xf,data_gain,/xsty,xtit='Frequency (MHz)',ytit='Time-integrated values',$
              title='Cleaned 10% Normalization Curve'+$
              strtrim(string(cgnumber_formatter(tmin,decimals=2)),1)+'_to_'+$
              strtrim(string(cgnumber_formatter(tmax,decimals=2)),1),$
             yrange=[min(data_gain),max(data_gain)],/ystyle
          endif 
       endif   
          
       If keyword_set(avg_norm) then begin
         If keyword_set(PS) and kl eq 0 then begin
           plot,xf,data_avg,/xsty,xtit='Frequency (MHz)',ytit='Time-integrated values',tit='Cleaned Average Normalization Curve',$
            /ystyle,yrange=[min(data_avg),max(data_avg)]
         endif
       endif
      
      ;*******************
      ;Normalize the Data 
      ;*******************
      If keyword_set(LOFAR) or keyword_set(NENUFAR) then begin 
        If keyword_set(avg_norm) then begin
          data = temporary(data)/rebin(reform(data_avg,1,nf),nt,nf)
        endif  
      
        If keyword_set(quantile_norm) then begin
          data = temporary(data)/rebin(reform(data_gain,1,nf),nt,nf) 
        endif 
     endif 
       ;Apply mask
       ;w =where(p2 eq 1)  ; good mask
       ;xm= median(data(w))
       ;w = where(p2 eq 0) ;bad mask
      ; data(w)= xm
      
      If keyword_set(PS) and kl eq 0 then begin
         label   = 'Normalized'
         STANDARD_PLOTS, data,p2,time,xf,label,xunit = 'secs',panel_ps=1
      endif      
   
      ;--------------------
      ;Normalize around 0
      ;--------------------
      data    = (data - median(data))
    
      If keyword_set(paper_plot) then data_norm=data
    endif ;S ne 4 and ne 5 
    
    If keyword_set(PS) and kl eq 0 then begin
        If S eq 5 or S eq 4 then label   = 'VPrime=vprime*Icor' else label = 'Zero-Mean Data'
        STANDARD_PLOTS, data,p2,time,xf,label,xunit = 'secs',panel_ps=1
    endif

    ;*************************************
    ;         Add: Jupiter 
    ;*************************************
    If S eq 4 and keyword_set(run_JUPITER) then begin
      print,'********** Start Jupiter *********************'
      If kl eq 0 then tmin_J = tmin_J_main
      tmax_J =  tmin_J + N_time  ;set end time
      print,'Min and max',tmin_J,tmax_J
            
      If tmax_J ge tmax_J_main then begin
        print,'******Skipping Jupiter******'
        goto, SKIPJUPITER2
      Endif
      lofar_jupiter,input_file_J,data,xf,xt,S,tmin_J,tmax_J,data_new=data_new,$
                  factor=factor,data_Jupiter=data_Jupiter,p2_Jupiter=p2_Jupiter,t_J=t_J,freq_J=freq_J,$
                  filename_p2=filename_p2, tf_J=tf_J, nf_J=nf_J,fmin_Jup=fmin_Jup,fmax_Jup=fmax_Jup,target_beam_J=target_beam_J,$
                  A_eff_L=A_eff_L, JA_eff_L=A_eff_L_Jup, Jupiter_background=Jupiter_background,scale=scale,LOFAR_RUN=LOFAR_RUN,$
                  file_Jup=file_Jup,f_nearest=f_nearest,time_p=time_p,freq_p=freq_p,$
                  UJupiter=data_Jupiter_U,QJupiter=data_Jupiter_Q,JBackV=Jupiter_background_V,JBackAll=Jupiter_Back_All;,VSkyBeam=VSky_Prime      
        If kl eq 0 then begin
          STANDARD_PLOTS, data_Jupiter,data_Jupiter,t_J,freq_J,'Jupiter Norm',xunit = 'secs',panel_ps=1
          STANDARD_PLOTS, factor,factor,t_J,freq_J,'Jupiter Factor',panel_ps=1
        Endif
       
      data = temporary(data_new);V3 
      ;data = temporary(data)/rebin(reform(data_gain,1,nf),nt,nf)

      ;-----------------
      ;I3Cor
      ;-------------
      ;S_I = 0 
      ;file_Jup = '/data/jake.turner/exoplanet/LC7_013/Jupiter/L568467/L568467_SAP000_B000_S0_P000_bf'
      ;lofar_jupiter,input_file_J,data_I,xf,xt,S_I,tmin_J,tmax_J,data_new=data_new,$
      ;            factor=factor,data_Jupiter=data_Jupiter,p2_Jupiter=p2_Jupiter,t_J=t_J,freq_J=freq_J,$
      ;            filename_p2=filename_p2, tf_J=tf_J, nf_J=nf_J,fmin_Jup=fmin_Jup,fmax_Jup=fmax_Jup,target_beam_J=target_beam_J,$
      ;;            A_eff_L=A_eff_L, JA_eff_L=A_eff_L_Jup, Jupiter_background=Jupiter_background,scale=scale,LOFAR_RUN=LOFAR_RUN,$
      ;            file_Jup=file_Jup,f_nearest=f_nearest,time_p=time_p,freq_p=freq_p,$
      ;            UJupiter=data_Jupiter_U,QJupiter=data_Jupiter_Q,JBackV=Jupiter_background_V,JBackAll=Jupiter_Back_All      
      ;;  If kl eq 0 then begin
      ;     STANDARD_PLOTS, data_Jupiter,data_Jupiter,t_J,freq_J,'Jupiter Norm (I)',xunit = 'secs',panel_ps=1
     ;     STANDARD_PLOTS, factor,factor,t_J,freq_J,'Jupiter Factor (I)',panel_ps=1
     ;   endif 
       ;MAKE_BACKGROUND,data_I*p2,'', data_mean_I,ss,nn
     ;    data_gain = fltarr(nf) 
      ; for iii=0,nf-1 do data_gain[iii] = dyn_n(data_new(*,iii),quantile) 
      ; LE_AUTO_S,data_gain,101,patrol_value,0,data_gain,pnet 
       
      ;   data = data*data_new/rebin(reform(data_gain,1,nf),nt,nf) ;Vnew = V3*I3 
       
 
      SKIPJUPITER2:
      print,'********** End Jupiter *********************'  
    Endif  ;end Jupiter
                  
    ;---------------------------
    ; Make Plot setup
    ;--------------------------
    nnv=nt/2000.        ; for visualization in ps file ;2000
    if nnv ne long(nnv) then nnv=long(nnv)+1 else nnv=long(nnv)
    nntv=long(nt/nnv)
    mmv=nf/1500.        ;1000
    if mmv ne long(mmv) then mmv=long(mmv)+1 else mmv=long(mmv)
    mmtv=long(nf/mmv)
    nwin=0 & xp=100 & yp=100
       
    If keyword_set(VERBOSE) then begin
      print,'-------------------------'
      print,'----Start RFI Mitigation----'
      print,'-------------------------'
    endif 
  
    ;********************************************************************
    ;********************************************************************
    ;                    RFI Mitigation 
    ;********************************************************************
    ;;********************************************************************
    ;Uses: RFI mitigate (Zarka) 
    ;       Inputs: x(t,f), p weigths, PATROL, le_sig runs le_auto_s
   
    ;*****************
    ;   Time Series
    ;*****************
     REDUCE_ARRAY, data, [1,nf], yavg_before,/dbl   ;create time series
     REDUCE_ARRAY, data, [nt,1], xavg_before,/dbl   ;create integrated spectrum
     If keyword_set(VERBOSE) then begin
      print,'Stdev of Integrated Spectra: ', stddev(xavg_before)
      print,'Stdev of Integrated Time Series: ', stddev(yavg_before)
     endif 
   
    IF keyword_set(RFI) then begin 
      If keyword_set(VERBOSE) then begin
       print,'*******************'
       print,'   Running RFI     '
       print,'*******************'
      endif   
        If kl eq 0 then begin
          PS_IN      = 1
          LAST_PS_IN = 0       
        Endif
        
        If kl eq 100 then begin
          PS_IN      = 1
          LAST_PS_IN = 0
        Endif
        
        If kl ne 0 and kl ne 100 then begin 
          PS_IN      = 0
          LAST_PS_IN = 1  
        endif 
        
        RFI_MITIGATE,data,p2,time=time,freq=xf,PATROL=patrol,LE_SIG=le_sig,$
          SIGMA_PATROL=patrol_value,LESIG_SIGMA=le_sig_value,SUM=SUM,PEX=PEX,xnet=xnet,PS=PS_IN,LAST_PS=LAST_PS_IN
     endif  ;end RFI
    
    ;If not running RFI- read in save file from last run
    If not(keyword_set(RFI)) then begin
      If keyword_set(VERBOSE) then begin
       print,'******************************'
       print,'Restoring Cleaned Data for RFI'
       print,'******************************'
      endif 
      data     = data_clean[nt*kl:(kl+1)*nt-1,*]
      p2       = p2_full[nt*kl:(kl+1)*nt-1,*]
    endif  
     
    ;***********************************************
    ;data is median where mask is located
    ;***********************************************
    data[where(p2 eq 0)]  =    mean(data[where(p2 eq 1)])

    ;If S eq 4 then begin
    ;  data = data*data_I*p2
    ;endif
    ;If S eq 4 then begin
    ;  STANDARD_PLOTS, data,p2,time,xf,'V=vprime*Icor',xunit = 'secs',panel_ps=1;,/correct_data
    ;  reduce_array,data,[nt,1],data_f
    ;  LE_AUTO_S,data_f,101,patrol_value,0,data_f,pnet
    ;  cgplot,xf,data_f,xtitle='Frequency',ytitle='Intensity' 
    ;endif
    

    ;------------------------------------------------------
    ;           Post-Processing Method2
    ;------------------------------------------------------
    If keyword_set(Process_METHOD2) then begin
      ;If S eq 0 then begin
       If kl eq 0 then begin
        lofar_process_method2,data_raw,p2,time,xf,S,quantile,SAMPLING_TIME,channel_width,patrol_value=patrol_value,freq_response=fres_step,p2_f=p2_f,$
          /PS,/VERBOSE,/RFI,F2=F2_step,STDEV=stdev_step,IDATA=data_I
       endif
      
       If kl ne 0 then begin
        lofar_process_method2,data_raw,p2,time,xf,S,quantile,SAMPLING_TIME,channel_width,patrol_value=patrol_value,freq_response=fres_step,p2_f=p2_f,$
        /RFI,F2=F2_step,IDATA=data_I,STDEV=stdev_step
       endif
       freq_response[kl,*]  = fres_step(*)
       pf_full[kl,*]        = p2_f(*)
       ;save step
      ;endif ;S =0 
      If S eq 4 then begin 
        data_mean_save(kl,*)  = fres_step(*)
        ;data_mean2_save(kl,*)  = F2_step(*)
       ;data_stdev_save(kl,*) = stdev_step(*)
      endif 
   endif   
            
    ;****************************
    ;    Save Cleaned data and mask
    ;****************************
    If keyword_set(RFI) then begin 
        If keyword_set(full_clean_save) then begin
          data_clean(kl*nt:((kl+1)*nt)-1,*) = data
        endif 
        
        ;Save Mask in bits and maniuplate in bits
        If keyword_set(mask_bit) and keyword_set(bit) then begin
          print,'Save mask in bits'
          BYTEARRAY_TO_BITARRAY, p2, xsize_step, p2_bit
          p2_full_bit_p[kl] = ptr_new(p2_bit)
        endif
        
        If not(keyword_set(bit)) then begin
          p2_full(kl*nt:(kl+1)*nt-1,*)      = p2   
        Endif
    endif 
    
    ;*************
    ; RFI Stats
    ;*************
    ww = where(p2 eq 0)  ;bad pixels 
    p0=n_elements(where(p2 eq 0))*100.d0/n_elements(p2)   ;percent of bad pixels 
    
    print,'# polluted channels = ',n_elements(ww),' / ',n_elements(p2)
    print,'RFI Mask:       ',p0,' % -> masked out'
    
    ;Plot showing Before & After (really good for papers)
    If keyword_set(paper_plot) then begin
        STANDARD_PLOTS, data_norm,p2*1.,time,xf,'Before RFI',/NO_ZOOM,/ONLY_PLOT,xunit = 'secs',panel_ps=1
        STANDARD_PLOTS, data,p2*1.,time,xf,'After RFI',/NO_ZOOM,/ONLY_PLOT,xunit = 'secs',panel_ps=1
    endif 

    ;****************
    ;update time loop
    ;****************
     tmin = tmin + N_time 
       
     If keyword_set(Run_JUPITER) then tmin_J = tmin_J + N_time
  endfor      ;end of time loop
   
 ;Save all the variables (useful for debugging but huge!)
 If keyword_set(allsave) then begin
    save,/ALL,file='all_'+save_file_name+'.sav' 
    print, 'Saved All Data: ', 'all_'+save_file_name+'.sav'
 endif  
 
 ;save the clean data (full; debugging)
 If keyword_set(full_clean_save) then begin
    save,xf,xt,data_clean,p2_full,file='dataclean_full_'+save_file_name+'.sav' 
 endif 
 
 ;----------------------------------------------------------------- 
 ;                        save mask 
 ;-----------------------------------------------------------------
 If keyword_set(mask_save) then begin
   ;--------------
   ;Mask in Bytes
   ;-------------
   If not(keyword_set(mask_bit)) and not(keyword_set(bit)) then begin   
    ;STANDARD_PLOTS, p2_full*1.,p2_full*1.,xt,xf,'Full Mask',/mask,xunit = 'secs'
    print,'************Save Byte mask**********'
    save,xt,xf,p2_full,file='mask_full_'+save_file_name+'.sav'
   endif
   
   ;-----------------------------------
   ;Save mask in bites for each step
   ;-----------------------------------
   If keyword_set(mask_bit) and keyword_set(bit)  then begin
    print,'************Save Bit mask (steps)**********'
    save,xt,xf,xsize_step,p2_full_bit_p,file='mask_full_bit_steps_'+save_file_name+'.sav'
   endif
   
   ;------------------------------
   ;Convert byte mask to bit mask 
   ;------------------------------
   If keyword_set(mask_bit) and not(keyword_set(bit)) then begin
     print,'************Save Bit mask (full)**********'
     print,'**Convert byte mask to bits***'
     ;lofar_bytebit,p2_full,nt,p2_out=p2_full_bit,/BYTE_TO_BIT
     BYTEARRAY_TO_BITARRAY, p2_full, xsize, p2_full_bit
     save,xt,xf,xsize,p2_full_bit,file='mask_full_bit_'+save_file_name+'.sav'
   endif
 endif ;save mask
 
 If not(keyword_set(mask_save)) then begin
   print, '************* WARNING *************************'
   print, 'You didnt save the mask: You will have problems....'
   print, '************* WARNING *************************'
 endif  

 If S eq 5 or S eq 4 then begin
  ;save,xt_rebin,xf,data_mean_save,data_mean2_save,data_stdev_save,pf_full,filename='DataMean_'+save_file_name+'.sav'
  save,xt_rebin,xf,data_mean_save,pf_full,filename='DataMean_'+save_file_name+'.sav'
 Endif 

 If keyword_set(save_datagain) then begin
   If S ne 5 and S ne 4 then save,data_gain_save,filename='DataGain_'+save_file_name+'.sav' 
 endif 
 
 If keyword_set(Process_METHOD2) then begin
    save,xt_rebin,xf,pf_full,freq_response,filename='Freq_Response_'+save_file_name+'.sav'
 endif  

 If keyword_set(PS) then begin
   device,/close
   cgfixps, filename_ps+'_RFI.ps'
   cgps2pdf,filename_ps+'_RFI.ps',/delete
   EXIT_PS
 endif

 ;update save_file_name
 save_file_name2 = save_file_name
 
 print,'SYSTIME entire RFI program=',SYSTIME(1) - st, ' sec'
 PRINT, 'Memory required for RFI: ', (MEMORY(/HIGHWATER))/1d9 ,' Gb'
 print,'The End of RFI Program'
 return
 
 ;free memory from any pointers
  HEAP_GC,/PTR,/VERBOSE
end
