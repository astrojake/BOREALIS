;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
;%%%%%%%%%%%%%%%%%%%%%%%%% Apply PCA Corrections %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
pro borealis_ApplyPCACorrections,tmin,N_time,steps,fmin,fmax,beam,S,file_root,date,save_file_name,$
  disp_values=disp_values,$
  PS=PS,VERBOSE=VERBOSE,$
  rebin_freq=rebin_freq, rebin_time=rebin_time,rebin_data=rebin_data,$
  TimeFreq_correction=TimeFreq_correction,Apply_Corrections=Apply_Corrections,skip_rebinloop=skip_rebinloop,$
  save_RFIMaskData=save_RFIMaskData,Smooth_Surface=Smooth_Surface,$
  mask_bit=mask_bit,bit=bit,$
  use_combinemask=use_combinemask,$
  Run_JUPITER=Run_JUPITER,REMOVE_BADF=REMOVE_BADF,N_smooth=N_smooth,Process_METHOD2=Process_METHOD2,$
  planet=planet,$
  LOFAR=LOFAR,NENUFAR=NENUFAR,filename=filename,N_SYSREM=N_SYSREM,data=data_rebin
  
  ;start time
  st_start= SYSTIME(1)

  If keyword_set(NENUFAR) then begin
    ;don't set bit and mask_bit because NENUFAR mask is in bytes
    ;Overide settings if differnet
    bit      = 0
    mask_bit = 0
  Endif

 ;Number of sysrem runs to do 
 If not(keyword_set(N_SYSREM)) then N_SYSREM = 3

  ;-----------------------
  ;Polarizations
  ;-----------------------
  ;convert Beam and Pol values into strings
  B_string = strtrim(string(beam),1)
  S_string2 = strtrim(string(S),1)        ;save filename
  If S lt 4 then S_in = S   ; I, Q, U, V
  If S eq 4 then S_in = 3   ; V^2
  If S eq 5 then S_in = 3   ; V'
  If S eq 6 then begin     ; L = sqrt(Q^2 + U^2)
    S_in = 1
    Sin1 = 1              ;Q
    Sin2 = 2              ;V
  endif
  If S lt 4 then S_string = strtrim(string(uint(S)),1)   ; I, Q, U, V   ;read file name
  If S eq 4 then S_string = strtrim(string(uint(3)),1)   ; V^2
  If S eq 5 then S_string = strtrim(string(uint(3)),1)   ; V'
  If S eq 6 then begin     ; L = sqrt(Q^2 + U^2)
    S_string = strtrim(string(uint(1)),1)
  endif

  If keyword_Set(disp_values) then begin
    RUN_DISP        = disp_values[0]
    size_disp       = disp_values[1]
    rebin_time_disp = disp_values[2]
    dm              = disp_values[3]
    If keyword_set(VERBOSE) then begin
      print, 'Size Disp', size_disp
      print, 'Disp Sample Time', rebin_time_disp
      print, 'dm', dm
    endif
  Endif


  ;***************************************************************************
  ;***************************************************************************
  If keyword_set(Apply_Corrections) then begin
    startmem_Apply = Memory(/current)

    If keyword_set(VERBOSE) then begin
      print,'----------------------'
      print,'  Apply Corrections   '
      print,'----------------------'
    endif

    ;------------------
    ;  Read Headers
    ;------------------
    If keyword_set(LOFAR) then begin
      
      file = file_root+date+'/'+date+'_SAP000_B00'+b_string+'_S'+S_string+'_P000_bf'
      If keyword_set(VERBOSE) then print, 'File:',file
      
      READ_H5_DATA, file,0,1,fmin,fmax,Beam,S_in,data,time,xf,nt,nf,Time=tf,/VERB,/NSPECTRA,$
        NOF_SAMPLES=NOF_SAMPLES,START_MJD=START_MJD,MJD_sampling_time=MJD_sampling_time,$
        ANTENNA_SET=ANTENNA_SET,SAMPLING_TIME=SAMPLING_TIME,NOF_stations=NOF_stations,$
        channel_width=channel_width
    endif ;end LOFAR 
    If keyword_set(NENUFAR) then begin
      Beam_in = Beam
      If S eq 0 then Stokes = 0 ;Stokes-I
      If S eq 1 then Stokes = 1 ;Stokes-Q
      If S eq 2 then Stokes = 2 ;Stokes-U
      If S eq 3 then Stokes = 3 ;Stokes-V
      If S eq 4 then Stokes = 4 ;Stokes-V Prime
      If S eq 0 then nstokes = 1 ELSE nstokes = 4
      temp       = filename+'*_'+strtrim(string(Beam_in),1)+'.spectra*'
      print, filename
      file_input = findfile(temp)
      print, 'file_input',file_input 
      READ_NU_FITS, file_input, command,param,variab, nt,dt,nf,df,ns,jd0,h0,fref,$
           x,xt,xf,beam_temp,ndata,corrt,corrf
         If S eq 0 then begin
           help, ndata
           help, x
           mask_nenufar = ndata[*,*]
           data_full = x[*,*,0] ; ;Stokes-I
         endif
         If S eq 3 then begin
           data_full = x[*,*,1] ; ;Stokes-V
           mask_nenufar = ndata[*,*]
         endif
         If S eq 4 then begin
           data_full = x[*,*,1] ;Stokes-V Prime
           mask_nenufar = ndata[*,*]
           data = data_full
         endif
         ;Restore I data and mask
         If S eq 3 or S eq 4 then begin
           save_file_name_I = date+'_'+planet+'_pol0_'+'beam'+strtrim(string(cgnumber_formatter(beam[0],decimals=0)),1)
           restore,filename='mask_full_'+save_file_name_I+'.sav'
            mask_nenufar_I = p2_full
            data_full_I = x[*,*,0]
            mask_nenufar = ndata[*,*]
            p2_combine = mask_nenufar_I*mask_nenufar
          ; p2_combine = mask_nenufar

         endif
      
      help, data3
      help, xt
      time = xt   ;fix everywhere
      help, nt_full 
      help, xf
      sampling_time = dt
      channel_width = df
      nof_samples = n_elements(nt)
      Start_MJD = jd0 - 2400000.5d 
      MJD_sampling_time = dt/86400.0d 
    Endif ;nenufar 
    ;----------------------------------------------
    ;----------------------------------------------
    ;               Setup Time Array
    ;----------------------------------------------
    ;----------------------------------------------
    MJD_t = dindgen(nt,start=START_MJD,INCREMENT=MJD_sampling_time)
    ;tmin= tmin_or
    
    ;-----------------------------------
    ;          Restore Mask
    ;-----------------------------------
    If keyword_set(VERBOSE) then print,'--Restore Mask----'
   save_file_name = 'test'
    filename_ps = save_file_name+'_processing'

;    ;-------------------
;    ;Full Mask in bytes
;    ;-------------------
;    If not(keyword_set(mask_bit)) then begin
;      If not(keyword_set(use_combinemask)) then begin
;        print,'Restore Mask: ', 'mask_full_'+save_file_name+'.sav'
;        restore,filename='mask_full_'+save_file_name+'.sav'
;        ;xt,xf,p2_full
;      endif
;      ;combined mask in bytes
;      If keyword_set(use_combinemask) then begin
;        save_file_name3 = strmid(save_file_name,0,strlen(save_file_name)-6); e.g. L429868_pol0_41sec_RFI_patrol5.5_pex_lesig3.5_sum
;        print,'Restore Mask: ','mask_combine_'+save_file_name3+'.sav'
;        restore,filename='mask_combine_'+save_file_name3+'.sav'
;        ;xt,xf,p2_full
;      Endif
;    endif
;
;    ;-----------------------------
;    ;Bit mask and convert to byte
;    ;-----------------------------
;    If keyword_set(mask_bit) and not(keyword_set(bit)) then begin
;      If not(keyword_set(use_combinemask)) then begin
;        print,'Restore Mask: ','mask_full_bit_'+save_file_name+'.sav'
;        restore,filename='mask_full_bit_'+save_file_name+'.sav' ;flag
;        ;xt,xf,xsize,p2_full_bit
;      endif
;      If keyword_set(use_combinemask) then begin
;        save_file_name3 = strmid(save_file_name,0,strlen(save_file_name)-6); e.g. L429868_pol0_41sec_RFI_patrol5.5_pex_lesig3.5_sum
;        print, 'Restore Mask: ','mask_combine_bit_'+save_file_name3+'.sav'
;        restore,filename='mask_combine_bit_'+save_file_name3+'.sav'
;      Endif
;    endif
;
;    ;-----------------------------
;    ;Bit mask and work in bits
;    ;-----------------------------
;    If keyword_set(mask_bit) and keyword_set(bit) then begin
;      If not(keyword_set(use_combinemask)) then begin
;        print,'Restore Mask: ','mask_full_bit_steps_'+save_file_name+'.sav'
;        restore,filename='mask_full_bit_steps_'+save_file_name+'.sav' ;flag
;        ;xt,xf,xsize_step,p2_full_bit_p
;      endif
;    Endif
;
;    nt = n_elements(xt)
;    nf = n_elements(xf)
;
;    ;---------------------------------------------------
;    ;Mask in bits but don't work in bits (not(bit))
;    ;---------------------------------------------------
;    If keyword_set(mask_bit) and not(keyword_set(bit)) then begin
;      If keyword_set(use_combinemask) then begin
;        BITARRAY_TO_BYTEARRAY,p2_full_all_bit,xsize,p2_full
;      Endif
;      If not(keyword_set(use_combinemask)) then begin
;        BITARRAY_TO_BYTEARRAY,p2_full_bit,xsize,p2_full
;      Endif
;    Endif

    If keyword_set(PS) and keyword_set(rebin_data) then begin
      set_plot,'PS'
      print,'Filename Ps:', filename_ps
      device,filename=filename_ps+'_rebin.ps',/landscape
    endif

    If keyword_set(run_JUPITER) then begin
      input_file_J = 'Input_Jupiter.dat'
      read_jupiter_input,input_file_J,tmin_Jup=tmin_J_main,tmax_Jup=tmax_J_main
      lofar_jupitersetup,input_file_J,data,xf,xt,S_in,tmin_J,tmax_J,data_new=data_new,$
        factor=factor,data_Jupiter=Jupiter_Norm,p2_Jupiter=p2_Jupiter,t_J=t_J,freq_J=freq_J,$
        filename_p2=filename_p2, tf_J=tf_J, nf_J=nf_J,fmin_Jup=fmin_Jup,fmax_Jup=fmax_Jup,target_beam_J=target_beam_J,$
        A_eff_L=A_eff_L, JA_eff_L=A_eff_L_Jup, Jupiter_background=Jupiter_background, scale=scale,LOFAR_RUN=LOFAR_RUN,$
        file_Jup=file_Jup,f_nearest=f_nearest,time_p=time_p,freq_p=freq_p,method2_Jup=method2_Jup,JBackV=Jupiter_background_V,JBackAll=Jupiter_Back_All
    Endif ;end Jupiter

    ;--------------------------------------------------
    ;                  Time Loop
    ;--------------------------------------------------
    ;for kl=0,steps-1 do begin  ;begin time loop-
    ;tmax =  tmin + N_time  ;set tmin

    ;  print,'-----------------------------------------------------------------------'
    ;  print,'reading slice ',kl,' start time: ', tmin, ' end time: ', tmax
    ;  print,'-----------------------------------------------------------------------'

      ;------------------------------
      ; Edit Time for selected input
      ;-----------------------------
    ;  u        = where(tf ge tmin and tf lt tmax)
    ;  mjd      = MJD_t(u)

      ;----------------------------------------------------------
      ;        Read Data for specific Beam and Polarization
      ;----------------------------------------------------------
      If keyword_set(NENUFAR) then begin
        print,'-----------------------------------------------------------------------'
        ;print,'Number of spectra in NENUFAR slice, nt: ', nt
        print,'-----------------------------------------------------------------------'
        If S ne 0 then begin
          data_I       = data_full_I
          p2_I         = mask_nenufar_I
        Endif
        
        data = data_full
        p2 = p2_combine
        
        If keyword_set(VERBOSE) then print,'File: ' ,file_input
      endif
     
      ;-------------------
      ;Create Mask
      ;-------------------
      ;p2=bytarr(nt,nf)+1b


;      If S eq 3 then begin
;        ;print,'Processing S=3, V with sign is used!!!'
;        print,'Processing S=3, abs(V) is used!!!'
;        data = abs(data)  ;V = |V|
;      endif ;s=3
;      If S eq 4 or S eq 5 then begin ;V'
;        print,' Running Vprime'
;        Print,'-------------------'
;        Print,'Reading I Data'
;        Print,'-------------------'
;        S_I_string = strtrim(string(0),1)
;        file_I = file_root+date+'/'+date+'_SAP000_B00'+B_string+$
;          '_S'+S_I_string+'_P000_bf'
;        READ_H5_DATA, file_I,tmin,tmax,fmin,fmax,Beam,0,data_I,t,f,nt,nf,REMOVE_BADF=REMOVE_BADF
;        If kl eq 0 then STANDARD_PLOTS, data,p2*1.,time,xf,'V Raw',xunit = 'secs',panel_ps=1
;        If kl eq 0 then STANDARD_PLOTS, data_I,p2*1.,time,xf,'I Raw',xunit = 'secs',panel_ps=1

        If S eq 4 then data = temporary(data)/data_I ; (V/I)
        If S eq 5 then data = temporary(data)/data_I  ;V/I
 ;       STANDARD_PLOTS, data,p2*1.,xt,xf,'V/I',xunit = 'secs',panel_ps=1

      ;endif ;s= 4 or 5
      If S eq 6 then begin ; L
        print,' Running L'
        data_U = temporary(data)
        S_Q_string = strtrim(string(2),1)
        file_I = file_root+date+'/'+date+'_SAP000_B00'+B_string+$
          '_S'+S_Q_string+'_P000_bf'
        READ_H5_DATA, file_I,tmin,tmax,fmin,fmax,Beam,2,data_Q,t,f,nt,nf,REMOVE_BADF=REMOVE_BADF
        data = sqrt(data_Q^2.0d + data_U^2.0d)
      endif ;L

      If keyword_set(VERBOSE) then begin
        print,'-------------------------'
        print,'----Done Reading --------'
        print,'-------------------------'
      endif

      ;*************************************
      ;         Add: Jupiter
      ;*************************************
      If S eq 0 and keyword_set(run_JUPITER) then begin
        print,'********** Start Jupiter *********************'
        If kl eq 0 then tmin_J = tmin_J_main
        tmax_J =  tmin_J + N_time  ;set end time
        print,'Min and max',tmin_J,tmax_J

        If tmax_J ge tmax_J_main then begin
          print,'******Skipping Jupiter******'
          goto, SKIPJUPITER
        Endif
        lofar_jupiter,input_file_J,data,xf,xt,S_in,tmin_J,tmax_J,data_new=data_new,$
          factor=factor,data_Jupiter=data_Jupiter,p2_Jupiter=p2_Jupiter,t_J=t_J,freq_J=freq_J,$
          filename_p2=filename_p2, tf_J=tf_J, nf_J=nf_J,fmin_Jup=fmin_Jup,fmax_Jup=fmax_Jup,target_beam_J=target_beam_J,$
          A_eff_L=A_eff_L, JA_eff_L=A_eff_L_Jup, Jupiter_background=Jupiter_background, scale=scale,LOFAR_RUN=LOFAR_RUN,$
          file_Jup=file_Jup,f_nearest=f_nearest,time_p=time_p,freq_p=freq_p
        data = temporary(data_new)  ;data updated with Jupiter
        SKIPJUPITER:
        print,'********** End Jupiter *********************'
      Endif  ;end Jupiter


;      ;-----------------
;      ;Restore Byte Mask
;      ;------------------
;      If not(keyword_set(bit)) then begin
;        p2 = p2_full[nt*kl:(kl+1)*nt-1,*]
;      endif

      ;---------------------------------
      ;Create Byte Mask from bit mask
      ;----------------------------------
      If keyword_set(mask_bit) and keyword_set(bit) then begin
        pointer_to_array,p2_full_bit_p[kl],out_array=p2_bit ;convert step into bit array
        BITARRAY_TO_BYTEARRAY, p2_bit, xsize_step,p2      ;convert bit array to byte array
      endif

   ;   If kl eq 0 then begin
        If keyword_set(PS) then begin
          label   = 'Data: Raw'
;          wbad = where(p2 eq 0)
;          wgood = where(p2 eq 1)
;          data[wbad] = median(data[wgood])
 ;         STANDARD_PLOTS, data,p2*1.,xt,xf,label,xunit = 'secs',panel_ps=1
        endif
     ; endif

      ;*************************************
     ; If kl eq 0 then begin
        ;---------------------------
        ;Set array for De-dispersion
        ;---------------------------
        If keyword_set(RUN_DISP) then begin
          If keyword_set(VERBOSE) then begin
            print,'***********************'
            print, 'Set arrays for de-dispersion'
            print,'**********************'
          endif
          If steps gt 1 then begin
            ;**************
            ;setup arrays
            ;**************

            ;Number of channels to combine to get the new channel width
            n_disp = long(size_disp/(channel_width/1e6))     ;channel_width in Hz
            mmv_1 = n_disp
            if mmv_1 ne long(mmv_1) then mmv_1=long(mmv_1)+1 else mmv_1=long(mmv_1)
            mmtv_1 = long(nf/mmv_1)

            ;New Sample Time
            N_reduce_time = long(rebin_time_disp/sampling_time) ;Number of bins to combine for time
            nnv_1 = N_reduce_time
            if nnv_1 ne long(nnv_1) then nnv_1=long(nnv_1)+1 else nnv_1=long(nnv_1)
            nntv_1=long(nt/nnv_1)

            x  =fltarr((steps-1)*nt,mmtv_1)
            w = x
            ;setup arrays
            y=fltarr(nt*2,nf)
            z = y
            p=bytarr(nt*2,nf)+1b
            q = p
          Endif ;steps gt 1

          ;De-disp doesn't make any sense
          If steps eq 1 then begin
            print, '****************************************'
            print,'You only have one step in de-dispersion!!'
            print, 'You might want one more'
            print, '****************************************'
          Endif
        endif ;end De-disp

        ;***************************************************************************
        ;                 Rebin array setup
        ;****************************************************************************
        If keyword_Set(rebin_data) then begin
          ;Number of channels to combine to get the new channel widtth
          Num_rebin_f = long((rebin_freq*1000.)/(channel_width))
          mmv_2 = Num_rebin_f
          mmtv_2 = long(nf/mmv_2)
          ;New Sample Time
          Num_rebin_t = long(rebin_time/sampling_time) ;Number of bins to combine for time
          nnv_2 = Num_rebin_t
          If Num_rebin_t eq 0 then Num_rebin_t = 1
          ;pointers
          data_rebin_p    = ptrarr(steps)
          p2_rebin_p      = ptrarr(steps)
          xt_rebin_p      = ptrarr(steps)
          mjd_rebin_p      = ptrarr(steps)
          UT_rebin_p            = ptrarr(steps)
          reduce_array,xf,[Num_rebin_f],xf_rebin,/dbl
        endif ;end rebin
      ;endif ;kl eq 0
      kl = 0
      ;-------------------------------------
      ;               Apply PCA 
      ;-------------------------------------
        If keyword_set(VERBOSE) then begin
          print,'Apply Time-Freq Correction at every freq'
        endif

      If not(keyword_set(gaussian_surface)) then begin
     ;     data = data*p2   ;apply mask 

          ;---------------
          ;T-F Correction
          ;---------------
          t_array = rebin(xt,nt,nf)                          ; (nt,nf)
          error_bar = 1.0d/sqrt(SAMPLING_TIME*channel_width)
          error = dblarr(nt,nf) + error_Bar
          
          If kl eq 0 then begin 
           reduce_array,data*p2,[nt,1],data_avg
           avg = rebin(reform(data_avg,1,nf),nt,nf)   ;just use the average of the first run
           ;plot,xt,avg,/xstyle,/ystyle
          endif 
          
          If S ne 0 then data_S =  data - avg    
          If S eq 0 then data_S = data/avg     
          data_S = temporary(data_S)*p2
         
          If kl eq 0 then STANDARD_PLOTS, data_S,p2,time,xf,'S: Pre-Sysrem',/NO_ZOOM 
          for i=0,N_SYSREM-1 do begin
            print,'V Systemaic: ',i
            sysrem,data_S,error,sys_err=sys_err;,/MAG
            data_S = data_S - sys_err
            STANDARD_PLOTS, sys_err,p2,time,xf,'SYSREM Surface: '+strtrim(string(i),1),/NO_ZOOM
            STANDARD_PLOTS, data_S,p2,time,xf,'S: Sysrem Number: '+strtrim(string(i),1),/NO_ZOOM
          endfor
          STANDARD_PLOTS, data_S,p2,time,xf,'Data After Sysrem',/NO_ZOOM
         
          
          data_S = data*p2

          If kl eq 0 then STANDARD_PLOTS, data_S,p2,time,xf,'S: Pre-Sysrem (no average)',/NO_ZOOM
          for i=0,N_SYSREM-1 do begin
            print,'V Systemaic: ',i
            sysrem,data_S,error,sys_err=sys_err,ITER=5;,/MAG
            data_S = data_S - sys_err
            STANDARD_PLOTS, sys_err,p2,time,xf,'SYSREM Surface (no average): '+strtrim(string(i),1),/NO_ZOOM,/panel_ps
            STANDARD_PLOTS, data_S,p2,time,xf,'S: Sysrem Number (no average): '+strtrim(string(i),1),/NO_ZOOM,/panel_ps
          endfor
          STANDARD_PLOTS, data_S,p2,time,xf,'Data After Sysrem (no average)',/NO_ZOOM,/panel_ps

          data  = data_S
;
;          data_S = rotate(data_S,1)
;          error  = rotate(error,1)
;          
;          for i=0,N_SYSREM-1 do begin
;            print,'V Systemaic (rotate): ',i
;
;            sysrem,data_S,error,sys_err=sys_err,ITER=5;,/MAG
;            data_S = data_S - sys_err
;     
;            data_S = data_S - sys_err
;            STANDARD_PLOTS, sys_err,rotate(p2,1),xf,time,'SYSREM Surface (rotate): '+strtrim(string(i),1),/NO_ZOOM
;            STANDARD_PLOTS, data,rotate(p2,1),xf,time,'S: Sysrem Number (rotate): '+strtrim(string(i),1),/NO_ZOOM
;          endfor
;
;          data = rotate(data,1)

          ;data = temporary(data_S)   
          
          ;If S eq 5 or S eq 4 then begin
            If S eq 5 and keyword_set(VERBOSE) then Print,'S = 5 (Vprime); subtract '
            If S eq 4 and keyword_set(VERBOSE) then  Print,'S = 4 (Vprime); subtract '
            
;            If S eq 4 then begin 
;               If kl eq 0 then begin
;                reduce_array,data_I*p2,[nt,1],data_avg_I
;                avg_I = rebin(reform(data_avg_I,1,nf),nt,nf)   ;just use the average of the first run
;              endif
              
;              data_S_I =  data_I/avg_I
;              for i=0,N_SYSREM-1 do begin
;                print,'I Systemaic: ',i
;                sysrem,data_S_I,error,sys_err=sys_err
;                data_S_I = data_S_I - sys_err
;                STANDARD_PLOTS, sys_err,p2,time,xf,'SYSREM Surface: '+strtrim(string(i),1),/NO_ZOOM 
;                STANDARD_PLOTS, data_S_I,p2,time,xf,'S: Sysrem Number; '+strtrim(string(i),1),/NO_ZOOM
;              endfor
              
           ;   data_S_I = data_S_I - 1.0 
              ;If kl eq 0 then STANDARD_PLOTS, data_I,p2,time,xf,'I (Sysrem); ',/NO_ZOOM

              ;data = data*p2              ;final Vhat
              ;If kl eq 0 then STANDARD_PLOTS, data_I,p2,time,xf,'Vprime (final); ',/NO_ZOOM

            data = data*p2   ;apply mask

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
                UJupiter=data_Jupiter_U,QJupiter=data_Jupiter_Q,JBackV=Jupiter_background_V,JBackAll=Jupiter_Back_All
              If kl eq 0 then begin
                STANDARD_PLOTS, data_Jupiter,data_Jupiter,t_J,freq_J,'Jupiter Norm',xunit = 'secs',panel_ps=1
                STANDARD_PLOTS, factor,factor,t_J,freq_J,'Jupiter Factor',panel_ps=1
              Endif

              data = temporary(data_new)
              SKIPJUPITER2:
              print,'********** End Jupiter *********************'
            Endif  ;end Jupiter

            If kl eq 0 then begin
              If keyword_set(PS) then begin
                label   = 'Data: Processing Vprime'
               ; STANDARD_PLOTS, data_prime,p2*1.,time,xf,label,/correct_data,xunit = 'secs',panel_ps=1
                label   = 'Data: Processing I '
              ;  STANDARD_PLOTS, data_I,p2*1.,time,xf,label,/correct_data,xunit = 'secs',panel_ps=1
              endif 
            endif;kl = 0
          ;endif ;S eq 4 or 5   
        endif  ;end not(keyword_set(gaussian_surface))

      ;------------------
      ;Apply SYSREM Mask
      ;------------------
     ; data = data - sys_err


      ;----------------------
      ;set RFI Channels
      ;----------------------
      If kl eq 0 then begin 
         wt = 0
         wn = 0  
      endif 
      ww=where(p2 eq 1)     ;good
      wbad=where(p2 eq 0)   ;bad
      wt = n_elements(p2) + wt
      wn = n_elements(wbad) + wn  ;add bad data to find total
      xm=median(data[ww])
      If S ne 5 and S ne 4 then data =  data   -  xm   ;subtract 1
      data(wbad) = 0        ;set masked data to mean

      If kl eq 0 then begin
        If keyword_set(PS) then begin
          label   = 'Data: Processing'
          STANDARD_PLOTS, data,p2*1.,time,xf,label,xunit = 'secs',panel_ps=1
        endif
      endif
      If keyword_set(rebin_data) then begin
        nf_rebin = n_elements(xf_rebin)
        nt_rebin = n_elements(xt_rebin)

        ;------------------------
        ;     Rebin in time/freq
        ;------------------------
        reduce_array,data*p2,[Num_rebin_t,Num_rebin_f],data_rebin2,/dbl
        reduce_array,p2*1.,[Num_rebin_t,Num_rebin_f],p2_rebin2,/dbl
        reduce_array,time,[Num_rebin_t],xt_rebin2,/dbl
        reduce_array,mjd,[Num_rebin_t],mjd_rebin2,/dbl
        reduce_array,UT,[Num_rebin_t],UT_rebin2,/dbl
        
        data_rebin_p[kl]  = ptr_new(data_rebin2)
        p2_rebin_p[kl]    = ptr_new(p2_rebin2)
        xt_rebin_p[kl]    = ptr_new(xt_rebin2)
        mjd_rebin_p[kl]   = ptr_new(mjd_rebin2)
        UT_rebin_p[kl]    = ptr_new(UT_rebin2)

        If keyword_set(VERBOSE) then begin
          label = 'Rebinned Data'
          SPDYNPS, data_rebin2/(p2_rebin2 + (p2_rebin2 eq 0)), min(xt_rebin2),max(xt_rebin2),min(xf_rebin),max(xf_rebin),$
            'Time (sec)','Frequency (MHz)',label,0,0,0,.05,.95,0,'.'
          SPDYNPS, p2_rebin2*1., min(xt_rebin2),max(xt_rebin2),min(xf_rebin),max(xf_rebin),$
            'Time (sec)','Frequency (MHz)',label,0,0,0,0,1,1,'.'
        endif
      Endif ;end rebin

      ;********************************************
      ;              De-dispersion
      ;********************************************
      ;Inputs:
      ; size_disp: size of disp in MHZ
      ; channel_width: width of frequency bins
      ; steps: number of steps to take in time
      ; nt: the number of time elements in one step
      ; num_bands: number of frequency bands
      ; n_disp: number of elements that correspond to new dispersion
      ; data: temp data of one slice
      ; kl: the time loop constant
      If keyword_set(RUN_DISP) then begin
        If keyword_set(VERBOSE) then begin
          print,'***********************'
          print, 'Running De-dispersion'
          print,'**********************'
        endif
        ;-------------------------------------------
        ;              De-dispersion
        ;------------------------------------------
        ;0 is special case since you don't de_disperse
        If kl eq 0 then begin
          y[nt:nt*2-1,*] = data;[0:nt-1,*]
          p[nt:nt*2-1,*] = p2;[0:nt-1,*]
        Endif

        If kl gt 0 then begin
          ;update y and p buffer (two buffers)
          y[0:nt-1,*]           =   y[nt:2*nt-1,*]      ; update the varabile from the last time
          y[nt:2*nt-1,*]        =   data[0:nt-1,*]         ;second block in temp new_data
          p[0:nt-1,*]           =   p[nt:2*nt-1,*]      ;second block temp p data set
          p[nt:2*nt-1,*]        =   p2[0:nt-1,*]
          for k=0,mmtv_1-1 do begin
            z(*,k*n_disp:(k*n_disp) + n_disp-1)= DE_DISP(y(*,k*n_disp:k*n_disp + n_disp-1),$
              xf(k*n_disp:k*n_disp + n_disp-1),dm,SAMPLING_TIME,xf[k*n_disp + n_disp-1])
          endfor

          for k=0,mmtv_1-1 do begin
            q(*,k*n_disp:k*n_disp  + n_disp-1)= DE_DISP(p(*,k*n_disp:k*n_disp+n_disp-1),$
              xf(k*n_disp:k*n_disp+n_disp-1),dm,SAMPLING_TIME,xf[k*n_disp + n_disp-1])
          endfor

          ;---------------------------------------------
          ;           Rebin in Frequency
          ;---------------------------------------------
          reduce_array,z[0:nt-1,*]*q[0:nt-1,*],[1,n_disp],x_temp
          x[(kl-1)*nt:(kl)*nt-1,*]= x_temp
          reduce_array,q[0:nt-1,*]*1.,[1,n_disp],w_temp
          w[(kl-1)*nt:(kl)*nt-1,*]=w_temp
        endif ;kl gt 0
      endif ;end de-disp

      ;****************
      ;update time loop
      ;****************
      ;tmin = tmin + N_time

      If keyword_set(Run_JUPITER) then tmin_J = tmin_J + N_time
   ; endfor ;endfor for time

    ;************************************
    ;        De-dispersed data
    ; ***********************************
    If keyword_set(RUN_DISP) then begin
      xt =  INDGEN((steps-1)*nt, /DOUBLE,increment=sampling_time,START=tmin_or) ;Create Time
      xf_disp=reverse(rebin(reverse(xf(0:mmtv_1*mmv_1-1)),mmtv_1,/sample))  ;rebin
      p = temporary(w)  ;mask
      dt = sampling_time
      save,x,p,xt,xf_disp,dt,dm,file=save_file_name+'_disp.sav'
      If keyword_set(VERBOSE) then begin
        print,'Just saved file: ', save_file_name+'_disp.sav'
      endif
    endif  ;end disp

    ;-----------------------------------------------
    ;Concatenate arrays for data_rebin and p2_Rebin
    ;-----------------------------------------------
    pointer_to_array,data_rebin_p,out_array=data_rebin
    pointer_to_array,p2_rebin_p,out_array=p2_rebin
    for j=0,steps-1 do begin ;n dates
      If j eq 0 then xt_rebin =  *xt_rebin_p[j]
      If j gt 0 then xt_rebin =  [ xt_rebin,*xt_rebin_p[j] ]
      If j eq 0 then mjd_rebin =  *mjd_rebin_p[j]
      If j gt 0 then mjd_rebin =  [ mjd_rebin,*mjd_rebin_p[j] ]
      If j eq 0 then UT_rebin =  *UT_rebin_p[j]
      If j gt 0 then UT_rebin =  [ UT_rebin,*UT_rebin_p[j] ]
    endfor

    ;********
    ;Plot
    ;********
    If keyword_set(PS) and keyword_set(rebin_data) then begin
      device,/close
      cgfixps, filename_ps+'_rebin.ps'
    endif

    If keyword_set(PS) then begin
      set_plot,'PS'
      print,'Filename Ps:', filename_ps
      device,filename=filename_ps+'_final.ps',/landscape
    endif

    ;Mask
    If keyword_set(PS) then begin
      If keyword_set(VERBOSE) then begin
        print, 'start mask'
      endif
      label   = 'Mask'
      STANDARD_PLOTS, p2_rebin*1.,p2_rebin*1.,xt_rebin,xf_rebin,label,/mask,xunit = 'secs',panel_ps=1
      If keyword_set(RUN_DISP) then begin
        label   = 'Mask (DISP)'
        STANDARD_PLOTS, p*1.,p*1.,xt,xf_disp,label,/mask,xunit = 'secs',panel_ps=1
      endif
    endif

    ;Data rebin
    If keyword_set(PS) and keyword_set(rebin_data) then begin
      label = 'Data: Rebin'
      STANDARD_PLOTS, data_rebin/(p2_rebin*1. + (p2_rebin*1. eq 0)),p2_rebin*1.,xt_rebin,xf_rebin,label,/histogram,xunit = 'secs',panel_ps=1
      STANDARD_PLOTS, data_rebin/(p2_rebin*1. + (p2_rebin*1. eq 0)),p2_rebin*1.,UT_rebin,xf_rebin,'Sky',/histogram,xunit = 'hour',xlabel='UT',/ONLY_PLOT,/NO_ZOOM,panel_ps=1
      If keyword_set(RUN_DISP) then begin
        label = 'Data: Rebin (DISP)'
        STANDARD_PLOTS, x/(p*1. + (p*1. eq 0)),p*1.,xt,xf_disp,label,/histogram,xunit = 'secs',panel_ps=1
      endif
    endif

    If keyword_set(PS) then begin
      device,/close
      cgfixps, filename_ps+'_final.ps'
      EXIT_PS
    endif

    ;*****
    ;Save
    ;*****
    save,xt_rebin,MJD_rebin,UT_rebin,xf_rebin,data_rebin,p2_rebin,filename='rebindata_'+save_file_name+'.sav'
    If keyword_set(VERBOSE) then begin
      print,'Just saved file: ', 'rebindata_'+save_file_name+'.sav'
    endif

    ;*************
    ; RFI Stats
    ;*************
    If keyword_set(mask_bit) then begin
      p0=wn*100.d0/wt    ;percent of bad pixels
      print,'# polluted channels = ',wn,' / ',wt
      print,'RFI Mask:       ',p0,' % -> masked out'
    endif

    ;mask byte
    If not(keyword_set(mask_bit)) then begin
      ww = where(p2_full eq 0)  ;bad pixels
      p0=n_elements(ww)*100.d0/n_elements(p2_full)   ;percent of bad pixels
      print,'# polluted channels = ',n_elements(ww),' / ',n_elements(p2_full)
      print,'RFI Mask:       ',p0,' % -> masked out'
    endif

    mem_apply = (MEMORY(/HIGHWATER))/1d9 ; Gb
  endif  ;End Apply Time-F Corrections

  If keyword_set(Apply_Corrections) then begin
    PRINT, 'Memory required for Applying Corrections: ', mem_apply ,' Gb'
  endif

  If keyword_set(Apply_Corrections) then begin
    PRINT, 'Time required for Applying Corrections: ', SYSTIME(1) - st_start ,' s'
  endif

  ;clear heap memory
  HEAP_GC,/PTR,/VERBOSE
end ;end Apply Corrections program

;;Processing file for the PCA analysis 
;pro borealis_processing_pca,tmin,N_time,steps,fmin,fmax,beam,S,file_root,date,save_file_name,$
;  disp_values=disp_values,$
;  PS=PS,VERBOSE=VERBOSE,$
;  rebin_freq=rebin_freq, rebin_time=rebin_time,rebin_data=rebin_data,$
;  TimeFreq_correction=TimeFreq_correction,Apply_Corrections=Apply_Corrections,skip_rebinloop=skip_rebinloop,$
;  save_RFIMaskData=save_RFIMaskData,Smooth_Surface=Smooth_Surface,$
;  mask_bit=mask_bit,bit=bit,$
;  use_combinemask=use_combinemask,$
;  Run_JUPITER=Run_JUPITER,REMOVE_BADF=REMOVE_BADF,N_smooth=N_smooth,Process_METHOD2=Process_METHOD2
;  
;  ;Save tmin
;  tmin_or = tmin
;  
;  ;************************************************
;  ;      Apply Time-Frequency Correction
;  ;************************************************
;
;  borealis_ApplyPCACorrections,tmin_or,N_time,steps,fmin,fmax,beam,S,file_root,date,save_file_name,$
;    disp_values=disp_values,$
;    PS=PS,VERBOSE=VERBOSE,$
;    rebin_freq=rebin_freq, rebin_time=rebin_time,rebin_data=rebin_data,$
;    TimeFreq_correction=TimeFreq_correction,Apply_Corrections=Apply_Corrections,skip_rebinloop=skip_rebinloop,$
;    save_RFIMaskData=save_RFIMaskData,Smooth_Surface=Smooth_Surface,$
;    mask_bit=mask_bit,bit=bit,$
;    use_combinemask=use_combinemask,$
;    Run_JUPITER=Run_JUPITER,REMOVE_BADF=REMOVE_BADF,N_smooth=N_smooth,Process_METHOD2=Process_METHOD2
;end   