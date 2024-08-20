;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
;%%%%%%%%%%%%%%%%%%%%%%%%% Apply PCA Corrections %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
pro borealis_applypcaCorrectionsv2,tmin,N_time,steps,fmin,fmax,beam,S,file_root,date,save_file_name,$
  disp_values=disp_values,$
  PS=PS,VERBOSE=VERBOSE,$
  rebin_freq=rebin_freq, rebin_time=rebin_time,rebin_data=rebin_data,$
  TimeFreq_correction=TimeFreq_correction,Apply_Corrections=Apply_Corrections,skip_rebinloop=skip_rebinloop,$
  save_RFIMaskData=save_RFIMaskData,Smooth_Surface=Smooth_Surface,$
  mask_bit=mask_bit,bit=bit,$
  use_combinemask=use_combinemask,$
  Run_JUPITER=Run_JUPITER,REMOVE_BADF=REMOVE_BADF,N_smooth=N_smooth,Process_METHOD2=Process_METHOD2,$
  planet=planet,$
  LOFAR=LOFAR,NENUFAR=NENUFAR,filename=filename,N_SYSREM=N_SYSREM,data=data
  
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
           data = data_full
         endif
         If S eq 3 then begin
           data_full = x[*,*,1] ; ;Stokes-V
           mask_nenufar = ndata[*,*]
           data = data_full
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
         endif   
            ;p2_combine = mask_nenufar_I*mask_nenufar
            ;p2_combine = mask_nenufar
      restore,filename='mask_full_'+save_file_name+'.sav'
      p2 = p2_full  ;restore mask
      
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
    
    
    ;-----------------------------------
    ;          Restore Mask
    ;-----------------------------------
    If keyword_set(VERBOSE) then print,'--Restore Mask----'
    filename_ps = save_file_name+'PCA_processing'

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
                
        If keyword_set(VERBOSE) then print,'File: ' ,file_input
      endif
     
      If S eq 4 then data = temporary(data)/data_I ; (V/I)
      If S eq 5 then data = temporary(data)/data_I  ;V/I
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
      kl = 0
      ;-------------------------------------
      ;               Apply PCA 
      ;-------------------------------------
        If keyword_set(VERBOSE) then begin
          print,'Apply PCA Time-Freq Correction at every freq'
        endif
        
          ;---------------
          ;T-F Correction
          ;---------------
          ;t_array = rebin(xt,nt,nf)                          ; (nt,nf)
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
         
          
;          If kl eq 0 then STANDARD_PLOTS, data_S,p2,time,xf,'S: Pre-Sysrem (no average)',/NO_ZOOM
;          for i=0,N_SYSREM-1 do begin
;            print,'V Systemaic: ',i
;            sysrem,data_S,error,sys_err=sys_err,ITER=5;,/MAG
;            data_S = data_S - sys_err
;            STANDARD_PLOTS, sys_err,p2,time,xf,'SYSREM Surface (no average): '+strtrim(string(i),1),/NO_ZOOM,/panel_ps
;            STANDARD_PLOTS, data_S,p2,time,xf,'S: Sysrem Number (no average): '+strtrim(string(i),1),/NO_ZOOM,/panel_ps
;          endfor
;          STANDARD_PLOTS, data_S,p2,time,xf,'Data After Sysrem (no average)',/NO_ZOOM,/panel_ps

          ;restore data 
          data  = data_S
          data = data*p2   ;apply mask
          
           If S eq 5 and keyword_set(VERBOSE) then Print,'S = 5 (Vprime); subtract '
           If S eq 4 and keyword_set(VERBOSE) then  Print,'S = 4 (Vprime); subtract '
           
      ;---------------------
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
      
    If keyword_set(PS) then begin
      ;set_plot,'PS'
      print,'Filename Ps:', filename_ps
      ;device,filename=filename_ps+'_final.ps',/landscape
    endif

    If keyword_set(PS) then begin
     ; device,/close
      ;cgfixps, filename_ps+'_final.ps'
      ;EXIT_PS
    endif

    ;*****
    ;Save
    ;*****
    save,data,filename='PCA_'+save_file_name+'.sav'
    If keyword_set(VERBOSE) then begin
      print,'Just saved file PCA: ', 'rebindata_'+save_file_name+'.sav'
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