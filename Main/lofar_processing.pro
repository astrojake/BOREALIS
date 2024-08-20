;***************************************************************
;      Processing Lofar exoplanet data
;**************************************************************
;***************************************************************
;Author: Jake Turner (UVA)
;        J.M. Griessmeier (CNRS) and Philippe Zarka (Obspm)
;
;Inputs: 
;      tmin     ;Starting time of observations (seconds)
;      N_time         ;Size of Steps in seconds
;      steps          ;number of steps to take in N_time
;      fmin     ;Starting Frequency
;      fmax       ;Ending Frequency
;      Beam           ;Beam (Beam 0 - Beam 4; 55cnc, Pulsar, Sky, Bright Source)
;      S              ;Polarization
;      file_root      ;root directory for raw data (eg. file_root = '/data/jake.turner/exoplanet/LC5_DDT_002/')
;      date           ;Dates in an array (eg. date = 'L429868'; folder must exist in root directory)
;      save_file_name ;save file name used for the sav files 
;
;Outputs
;     TimeFreqCorrect_array
;      
;Keywords
;      disp_values     ;De-disp values for the pulsar
;      disp            ;Run De-disp and fft for pulsar (values of pulsar are input as disp_values)
;      rebin           ;Rebin the data
;      rebin_freq      ;The requested freq sampling for rebin
;      rebin_time      ;The requested sampling time for rebin
;      time_correct    ;
;      spectral_correct
;      Apply_Corrections ; Apply the time-freq correction to the data
;      Smooth_Surface     ; Apply the 2-d surface fit (from the time-freq fit) 
;      mask_bit          ; Read the mask from the save file in bits (then it is converted to bytes)
;      bit               ; Manipulate the mask in bits 
;      use_combinemask      ; Read in combine mask instead of mask for each beam
;Uses:

pro lofar_timefreqcorrection,tmin,N_time,steps,fmin,fmax,beam,S,file_root,date,save_file_name,$
                     disp_values=disp_values,$
                     PS=PS,VERBOSE=VERBOSE,$
                     rebin_freq=rebin_freq,rebin_time=rebin_time,rebin_data=rebin_data,$
                     TimeFreq_correction=TimeFreq_correction,Apply_Corrections=Apply_Corrections,skip_rebinloop=skip_rebinloop,$
                     save_RFIMaskData=save_RFIMaskData,Smooth_Surface=Smooth_Surface,$
                     mask_bit=mask_bit,bit=bit,$
                     use_combinemask=use_combinemask,$
                     TimeFreqCorrect_array=TimeFreqCorrect_array,$
                     Run_JUPITER=Run_JUPITER,$
                     REMOVE_BADF=REMOVE_BADF,N_smooth=N_smooth,Process_METHOD2=Process_METHOD2,$
                     LOFAR=LOFAR,NENUFAR=NENUFAR,DYNSPECBF=DYNSPECBF,filename=filename,COMPLEX_NAME=COMPLEX_NAME,$
                     TO=TO,P=P
                     
  ;start time
  mem_start = (MEMORY(/HIGHWATER))/1d9 ; Gb
  st_start= SYSTIME(1)
  
  ;-----------------
  ;hardcoded values 
  ;-----------------    
  fac = 50.        ;rebin factor for saving raw data
  method2_again = 0   ;debuging
  If fac lt 11. then fac = 11
  If keyword_set(VERBOSE) then print, 'Rebin factor: ',fac
  
  ;----------------------
  ;Polarization Setup 
  ;----------------------
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

  ;----------------------------------------------------------
  ;              Setup Important Info for Data
  ;----------------------------------------------------------
  ;file to read in headers and arrays
  If keyword_set(LOFAR) then begin
   file = file_root+date+'/'+date+'_SAP000_B00'+b_string+'_S'+S_string+'_P000_bf'
   If keyword_set(VERBOSE) then begin
     print, 'File:',file
   endif
  endif 
 
  tmin_or = tmin
  filename_ps = save_file_name+'_processing'
 
  ;**********************************************
  ;       Find Time-Freq Correction
  ;**********************************************
     If keyword_set(TimeFreq_correction) then begin
       If keyword_set(PS) then begin
         set_plot,'PS'
         print,'Filename Ps:', filename_ps
         device,filename=filename_ps+'_TF.ps',/landscape
       endif  
      
       If keyword_set(VERBOSE) then begin
         print,'------------------------------------------------'
         print,'------------------------------------------------'
         print,'       Start Time-Frequency Correction Curve    '
         print,'------------------------------------------------'
         print,'------------------------------------------------'
       endif 
   
       ;-----------------------------------
       ;          Restore Mask
       ;-----------------------------------
       ;    Also read in Frequency range (xf)
       ;-----------------------------------   
       If keyword_set(VERBOSE) then begin
         print,'--Restore Mask----'
       endif
       
       ;--------------------------------
       ;Mask in bytes and work in bytes
       ;--------------------------------
       If not(keyword_set(mask_bit)) and not(keyword_set(bit)) then begin   
           If not(keyword_set(use_combinemask)) then begin
            print,'Restore Mask: ', 'mask_full_'+save_file_name+'.sav'
             restore,filename='mask_full_'+save_file_name+'.sav'
             ;xt,xf,p2_full   
           endif  
         If keyword_set(use_combinemask) then begin 
            save_file_name3 = strmid(save_file_name,0,strlen(save_file_name)-6); e.g. L429868_pol0_41sec_RFI_patrol5.5_pex_lesig3.5_sum
            print, 'Restore Mask:' ,'mask_combine_'+save_file_name3+'.sav'
            restore,filename='mask_combine_'+save_file_name3+'.sav'
             ;xt, xf, p2_full
         Endif
       endif ;not mask bit
        
      ;---------------------------------------------------
      ;Mask in bits but don't work in bits
      ;---------------------------------------------------
      If keyword_set(mask_bit) and not(keyword_set(bit)) then begin
        If not(keyword_set(use_combinemask)) then begin
          print,'Restore mask: ', 'mask_full_bit_'+save_file_name+'.sav'
          restore,filename='mask_full_bit_'+save_file_name+'.sav' ;flag
            ;xt,xf,xsize_step,p2_full_bit
          BITARRAY_TO_BYTEARRAY,p2_full_bit,xsize,p2_full
        Endif
         If keyword_set(use_combinemask) then begin
            save_file_name3 = strmid(save_file_name,0,strlen(save_file_name)-6); e.g. L429868_pol0_41sec_RFI_patrol5.5_pex_lesig3.5_sum
            print,'Restore mask: ','mask_combine_bit_'+save_file_name3+'.sav'
            restore,filename='mask_combine_bit_steps'+save_file_name3+'.sav'
              ;xt,xf,xsize,p2_full_bit
            BITARRAY_TO_BYTEARRAY,p2_full_bit,xsize,p2_full
        Endif
      Endif

      ;-------------------------------------
      ;Read I mask for all polarizations 
      ;-------------------------------------
      If S ge 1 then begin
        If keyword_set(COMPLEX_NAME) then begin 
          temp = strmid(save_file_name,12,strlen(save_file_name)-5-13) ; e.g. _pol0_41sec_RFI_patrol5.5_pex_lesig3.5_sum
          save_file_name3 = date+'_pol'+strtrim(string(0),1)+temp+'_beam'+b_String
        Endif        
        If not(keyword_set(COMPLEX_NAME)) then save_file_name3 = date+'_pol'+strtrim(string(0),1)+'_beam'+b_String
       If keyword_set(mask_bit) and keyword_set(bit)  then begin
         print,'Restore Mask: ','mask_full_bit_steps'+save_file_name3+'.sav'
         restore,filename='mask_full_bit_steps_'+save_file_name3+'.sav'   
          ;xt,xf,xsize_step,p2_full_bit_p
         p2_full_bit_p_I = temporary(p2_full_bit_p)
         xsize_step_I    = temporary(xsize_step)
       endif
      endif  
      
      ;-----------------------------
      ;Bit mask and work in bits
      ;-----------------------------
      If keyword_set(mask_bit) and keyword_set(bit) then begin
        If not(keyword_set(use_combinemask)) then begin
          print,'Restore Mask: ','mask_full_bit_steps_'+save_file_name+'.sav'
          restore,filename='mask_full_bit_steps_'+save_file_name+'.sav' ;flag
          ;xt,xf,xsize_step,p2_full_bit_p
        endif
      Endif
        
       ;----------------------------------------------
       ;               Setup Time Array/Headers
       ;---------------------------------------------- 
       If keyword_set(LOFAR) then begin
          READ_H5_DATA, file,0,1,fmin,fmax,Beam,S_in,x,t,xf,nt,nf,Time=tf,/NSPECTRA,$
         NOF_SAMPLES=NOF_SAMPLES,START_MJD=START_MJD,MJD_sampling_time=MJD_sampling_time,$
         ANTENNA_SET=ANTENNA_SET,SAMPLING_TIME=SAMPLING_TIME,NOF_stations=NOF_stations,$
         channel_width=channel_width
       endif 
       
       If keyword_set(NENUFAR) then begin
         Beam_in = Beam
         If S eq 0 then Stokes = 0 ;Stokes-I
         If S eq 1 then Stokes = 1 ;Stokes-Q
         If S eq 2 then Stokes = 2 ;Stokes-U
         If S eq 3 then Stokes = 3 ;Stokes-V
         If S eq 4 then Stokes = 4 ;Stokes-V Prime
         If S eq 0 then nstokes = 1 ELSE nstokes = 4
         read_nenufar, file_root+filename, temporary(data2),tt,ff,temporary(Beam_in),ndata, ntt,dt,nnf,df,nns,Stokes,nstokes=nstokes,$
           tmin=0,tmax=1,fmin=fmin,fmax=fmax,/VERBOSE,CHAN=CHAN,BEAMLETS=BEAMLETS,/REMOVE_BADF,/CORRECT_BANDPASS,/DATA_BANDPASS
           
         sampling_time = dt
         channel_width = df
       Endif
       
       nt = n_elements(xt)
       nf = n_elements(xf)
      ;------------------------------------------------------------------------------------------------
      ;------------------------------------------------------------------------------------------------
      ;                               Method 1: Line by line curve fitting
      ;-----------------------------------------------------------------------------------------------
      ;------------------------------------------------------------------------------------------------
      If not(keyword_set(Process_METHOD2)) or keyword_set(method2_again) then begin
          print,'**************************************'
          print,'************Method 1 Processing ******************'
          print,'**************************************'
          If keyword_set(run_JUPITER) then begin
            input_file_J = 'Input_Jupiter.dat'
            read_jupiter_input,input_file_J,tmin_Jup=tmin_J_main,tmax_Jup=tmax_J_main,/VERBOSE
            lofar_jupitersetup,input_file_J,data,xf,xt,S_in,tmin_J,tmax_J,data_new=data_new,$
              factor=factor,data_Jupiter=Jupiter_Norm,p2_Jupiter=p2_Jupiter,t_J=t_J,freq_J=freq_J,$
              filename_p2=filename_p2, tf_J=tf_J, nf_J=nf_J,fmin_Jup=fmin_Jup,fmax_Jup=fmax_Jup,target_beam_J=target_beam_J,$
              A_eff_L=A_eff_L, JA_eff_L=A_eff_L_Jup, Jupiter_background=Jupiter_background, scale=scale,LOFAR_RUN=LOFAR_RUN,$
              file_Jup=file_Jup,f_nearest=f_nearest,time_p=time_p,freq_p=freq_p
          Endif ;end Jupiter

        ;--------------------------------------------------
        ;                  Time Loop
        ;--------------------------------------------------
        If not(keyword_set(skip_rebinloop)) then begin
           for kl=0,steps-1 do begin  ;begin time loop
             tmax =  tmin + N_time  ;set end time
             
             If keyword_set(VERBOSE) then begin
               print,'-----------------------------------------------------------------------'
               print,'reading slice ',kl,' start time: ', tmin, ' end time: ', tmax
               print,'-----------------------------------------------------------------------'
             endif 
             If kl eq 0 then VERB = 1 Else VERB = 0
  
             If keyword_set(LOFAR) then begin
               READ_H5_DATA, file,tmin,tmax,fmin,fmax,Beam,S_in,data,time,xf,nt,nf,VERB=VERB,REMOVE_BADF=REMOVE_BADF
             Endif
             If keyword_set(NENUFAR) then begin
                beam_in = beam
                If S eq 0 then Stokes = 0 ;Stokes-I
                If S eq 1 then Stokes = 1 ;Stokes-Q
                If S eq 2 then Stokes = 2 ;Stokes-U
                If S eq 3 then Stokes = 3 ;Stokes-V
                If S eq 4 then Stokes = 4 ;Stokes-V Prime
                If S eq 0 then nstokes = 1 ELSE nstokes = 4
                read_nenufar, file_root+filename,data,time,xf,temporary(beam_in),ndata,nt,dt,nf,df,ns,Stokes,nstokes=nstokes,tmin=tmin,tmax=tmax,fmin=fmin,fmax=fmax,/REMOVE_BADF,VERBOSE=VERB,$
                  /CORRECT_BANDPASS,/DATA_BANDPASS
             endif
             If keyword_set(DYNSPECBF) then begin
               If S eq 0 then POL = 1 ;Stokes-I
               If S eq 1 then POL = 2 ;Stokes-Q
               If S eq 2 then POL = 3 ;Stokes-U
               If S eq 3 then POL = 4 ;Stokes-V
               If S eq 4 then POL = 4 ;Stokes-V^2
               If S eq 5 then POL = 4 ;Stokes-V prime
               read_dynspec_fits,file_root+filename+'_beam'+B_string,tmin,tmax,fmin,fmax,POL,data,time,xf,nt,nf
               stop
             Endif
             
             
             If keyword_set(VERBOSE) then begin
              print,'-------------------------'
              print,'----Done Reading --------'
              print,'-------------------------'
             endif 
             
             If S eq 4 or S eq 5 then begin ;V'
             	  print,' Running Vprime'
             	  Print,'-------------------'
             	  Print,'Reading I Data'
             	  Print,'-------------------'
             	  If keyword_set(NENUFAR) then begin 
               	  Beam_in = Beam
               	  Stokes = 0 ;Stokes-I
               	  nstokes = 1
               	  read_nenufar, file_root+filename,data_I,time,xf,temporary(Beam_in),ndata, nt,dt,nf,df,ns,Stokes,nstokes=nstokes,$
             	      tmin=tmin,tmax=tmax,fmax=fmax,fmin=fmin,/REMOVE_BADF,/CORRECT_BANDPASS,/DATA_BANDPASS
             	  endif 
             	  If keyword_set(LOFAR) then begin 
                	S_I_string = strtrim(string(0),1)
                	file_I = file_root+date+'/'+date+'_SAP000_B00'+B_string+$
                 	 '_S'+S_I_string+'_P000_bf'
                	READ_H5_DATA, file_I,tmin,tmax,fmin,fmax,Beam,0,data_I,t,f,nt,nf,REMOVE_BADF=REMOVE_BADF    
                endif 	 
               	If S eq 4 then data = temporary(data)/temporary(data_I) ; (V/I) 
               	If S eq 5 then data = temporary(data)/temporary(data_I)  ;V/I
             endif
             If S eq 6 then begin ; L
              print,' Running L'
              data_U = temporary(data)
              S_Q_string = strtrim(string(2),1)
              file_I = file_root+date+'/'+date+'_SAP000_B00'+B_string+$
                '_S'+S_Q_string+'_P000_bf'
              READ_H5_DATA, file_I,tmin,tmax,fmin,fmax,Beam,2,data_Q,t,f,nt,nf,REMOVE_BADF=REMOVE_BADF 
              data = sqrt(data_Q^2.0d + data_U^2.0d)
             end 
                              
             ;*************************************
             ;         Add: Jupiter 
             ;*************************************
             If keyword_set(run_JUPITER) then begin
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
                  A_eff_L=A_eff_L, JA_eff_L=A_eff_L_Jup, Jupiter_background=Jupiter_background,scale=scale,LOFAR_RUN=LOFAR_RUN,$
                  file_Jup=file_Jup,f_nearest=f_nearest,time_p=time_p,freq_p=freq_p
                data = data_new
                SKIPJUPITER:
                print,'********** End Jupiter *********************'  
             Endif  ;end Jupiter
             
             ;-----
             ;mask
             ;-----
             p2 = bytarr(nt,nf)+1b
            
             ;-----------------
             ;Restore Byte Mask
             ;------------------
             If not(keyword_set(bit)) then begin   
               p2 = p2_full[nt*kl:(kl+1)*nt-1,*]
             endif 
             
             ;---------------------------------
             ;Create Byte Mask from bit mask
             ;----------------------------------
             If keyword_set(mask_bit) and keyword_set(bit) then begin
               pointer_to_array,p2_full_bit_p[kl],out_array=p2_bit ;convert step into bit array
               BITARRAY_TO_BYTEARRAY, p2_bit, xsize_step,p2      ;convert bit array to byte array
             endif
             
             ;-------------------
             ;Method 2 Again
             ;---------------------
             If keyword_set(method2_again) then begin
                quantile     = 0.01
                print,'Quantile Method 2 Again',quantile
                patrol_value = 5.5
                RFI          = 1
                If kl eq 0 then begin
                  freq_response = dblarr(steps,nf) 
                  pf_full            = dblarr(steps,nf) 
                  xt_rebin       = dblarr(steps)
                endif  
                If kl eq 0 then VERBOSE_m = 1 Else VERBOSE_m = 0
                If kl eq 0 then  PS_m      = 1 Else Ps_m = 0 
           
                lofar_process_method2,data,p2,time,xf,quantile,SAMPLING_TIME,channel_width,$
                 patrol_value=patrol_value,freq_response=frequency_response,p2_f=p2_f,PS=PS_m,VERBOSE=VERBOSE_m,RFI=RFI
                
                freq_response(kl,*)= frequency_response
                pf_full(kl,*)           =p2_f
                
                ;If keyword_set(Process_METHOD2) then begin
                  reduce_array,time,[nt],temp_t
                 xt_rebin[kl] = temporary(temp_t)
                ;endif
                
                goto, skip_rebin 
             Endif
            
             ;rebin to save space 
             If kl eq 0 then begin
                 rebin_time2 = SAMPLING_TIME*fac   ;New Sample Time
                 num_kill = long(40./rebin_time2) ; remove last 40 seconds
                 num_rebin_t = fac
                 If Num_rebin_t eq 0 then Num_rebin_t = 1
                 data_rebin_p    = ptrarr(steps)  ;saves even more memory 
                 p2_rebin_p      = ptrarr(steps)
                 xt_rebin_p      = ptrarr(steps) 
             endif
             
             ;------------------------
             ;     Rebin in time
             ;------------------------ 
             reduce_array,data*p2,[Num_rebin_t,1],data_rebin2
             reduce_array,p2,[Num_rebin_t,1],p2_rebin2
             reduce_array,time,[Num_rebin_t],xt_rebin2
             data_rebin2 = data_rebin2/(p2_rebin2 + (p2_rebin2 eq 0))
             data_rebin_p[kl]  = ptr_new(data_rebin2)
             p2_rebin_p[kl]    = ptr_new(p2_rebin2)
             xt_rebin_p[kl]    = ptr_new(xt_rebin2)
             
             skip_rebin:
             
             ;****************
             ;update time loop
             ;****************
             tmin = tmin + N_time
             If keyword_set(Run_JUPITER) then tmin_J = tmin_J + N_time
          endfor ;endfor for time
          
          If keyword_set(method2_again) then goto,skip_all
            ;Concatenate arrays for data_rebin and p2_Rebin
            pointer_to_array,data_rebin_p,out_array=data_rebin
            pointer_to_array,p2_rebin_p,out_array=p2_rebin
            for j=0,steps-1 do begin 
              If j eq 0 then xt_rebin =  temporary(*xt_rebin_p[j])
              If j gt 0 then xt_rebin =  [ xt_rebin,temporary(*xt_rebin_p[j]) ]
            endfor ;end Concatenate arrays  
        
          If keyword_set(save_RFIMaskData) then begin
             print,'-----Save Raw Data--------'
             If S ne 3 then save,xf,xt_rebin,data_rebin,p2_rebin,filename  = save_file_name+'_dataraw_full.sav'
          endif 
        endif  ; not skiptime
            
        If keyword_set(skip_rebinloop) then begin
          If keyword_set(VERBOSE) then begin
            print,'Restore Data Raw'
          endif
          restore,filename  = save_file_name+'_dataraw_full.sav'
            ;xf,xt_rebin,data_rebin,p2_rebin
            
          If keyword_set(VERBOSE) then begin
             print, 'done'
          endif 
             
             rebin_time2 = sampling_time*fac
             num_kill = long(40./rebin_time2) ; remove last 40 seconds  
        endif     
            
        ;------------------------------------------
        ;         Time-Freq Correction
        ;------------------------------------------
          nt_rebin   = n_elements(xt_rebin)
          data_rebin = temporary(data_rebin[num_kill:nt_rebin-num_kill,*])  ; get rid of first and last 50 seconds for edge effects
          xt_rebin   = temporary(xt_rebin[num_kill:nt_rebin-num_kill])
          p2_rebin   = temporary(p2_rebin[num_kill:nt_rebin-num_kill,*])
          nt_rebin   = n_elements(xt_rebin)
          nf         = n_elements(xf)
          TimeFreqCorrect_array=dblarr(4,nf)
          nt_rebin   = n_elements(xt_rebin)
          !p.multi=[0,1,2]
          If keyword_set(VERBOSE) then begin
           print,'*****************************************'
           print,'Starting Time-Frequency Correction Curve'
           print,'*****************************************'
          endif
           
           skip = 0 
           ;start freq loop 
          for jjj=0,nf-1 do begin
             w  = where(p2_rebin(*,jjj) le 0.50,count_bad)  ;bad
             ww = where(p2_rebin(*,jjj) eq 1.,count_good)  ;good
             ratio_pixels = count_good/count_bad
             If ww[0] eq -1 or ratio_pixels lt 1 then begin 
               If keyword_set(VERBOSE) then begin
                Print,'Time all bad: ', jjj
               endif 
               If jjj eq 0 then begin
                fit_data = dblarr(3)
                 med = median(data_rebin(*,*))
                 TimeFreqCorrect_array(1,0)=med
                 TimeFreqCorrect_array(2,0)=0.
                 TimeFreqCorrect_array(3,0)=0.
                 fit_data[0] = TimeFreqCorrect_array(1,0)
                 fit_data[1] = TimeFreqCorrect_array(2,0)
                 fit_data[2] = TimeFreqCorrect_array(3,0)
               endif  
               If jjj gt 0 then begin
                fit_data[0] = TimeFreqCorrect_array(1,jjj-1) 
                fit_data[1] = TimeFreqCorrect_array(2,jjj-1) 
                fit_data[2] = TimeFreqCorrect_array(3,jjj-1) 
               endif 
               skip = 1 
               GoTo, skip   
             Endif 
             
             data_rebin[w,jjj] = median(data_rebin[ww,jjj])  ;set rfi = median 
             time              = xt_rebin(ww)                                    ;only good data
             tnet              = rebin(data_rebin(ww,jjj),n_elements(data_rebin(ww,jjj)),1)  
             tnet              = data_rebin(ww,jjj)
           
            mod_plot = jjj mod 1000
  
             win = 21. 
             If n_elements(tnet) le 40. then win=long(n_elements(tnet)/2.)
             LE_AUTO_S,tnet,win,3.5,0,tnet,pnet                 ;Additonal RFI
             tnet = smooth(tnet,win,/EDGE_TRUNCATE)
             ;fit data
             fit_data    = poly_fit(time,tnet,2)  ;fit
            
         skip:      
             ;Save fits
             TimeFreqCorrect_array(0,jjj) = xf[jjj]         ;frequency
             TimeFreqCorrect_array(1,jjj) = fit_data[0]    ;a
             TimeFreqCorrect_array(2,jjj) = fit_data[1]    ;b*x
             TimeFreqCorrect_array(3,jjj) = fit_data[2]    ;c*x^2    
            
           ;plot x000 
           mod_plot = jjj mod 1000
           
            If mod_plot eq 0 and skip eq 0 then begin
             If keyword_set(PS) then begin  
              If keyword_set(VERBOSE) then begin
               print, jjj
               print, 'Saved Plot'
              endif  
              cgplot,time,tnet,title='Time Correction Curve: Freq:'+strtrim(string(cgnumber_formatter(xf[jjj],decimals=4)),1)+' MHz',$
               xtitle='Time (s)',ytitle='Frequency-Integrated Values',$
               /xstyle, /ystyle,charsize=1;,yrange=[1.2e14,1.4e14]
              cgplot,time,(fit_data[0,0] + fit_data[0,1]*time + fit_data[0,2]*time^2.),color='red',/overplot
             endif  
            endif 
            
            skip = 0 
           
            If jjj eq nf-1 then begin
             If keyword_set(PS) then begin  
              If keyword_set(VERBOSE) then begin
               print, jjj
               print, 'Saved Plot'
              endif  
              cgplot,time,tnet,title='Time Correction Curve: Freq:'+strtrim(string(cgnumber_formatter(xf[jjj],decimals=4)),1)+' MHz',$
                xtitle='Time (s)',ytitle='Frequency-Integrated Values',/xstyle, /ystyle,charsize=1
              cgplot,time,(fit_data[0,0] + fit_data[0,1]*time + fit_data[0,2]*time^2.),color='red',/overplot
             endif 
            endif
          endfor ;end freq loop
          !p.multi=[0,1,1]
          
        skip_all:  
        If keyword_set(method2_again) then save,xt_rebin,xf,pf_full,freq_response,filename='Freq_Response_'+save_file_name+'.sav'
        end ;END METHOD 1 

     ;------------------------------------------------------------------------------------------------
     ;------------------------------------------------------------------------------------------
     ;                Method 2: Fitting Frequency Response (10% Quantile) over time
     ;------------------------------------------------------------------------------------------
     ;------------------------------------------------------------------------------------------------
        If keyword_set(Process_METHOD2) then begin
          print,'*************************************************'
          print,'***************Method 2**************************'
          print,'*************************************************'
                  
          If S ne 5 and S ne 4 then begin 
           restore,filename='Freq_Response_'+save_file_name+'.sav'
             ;--xt_rebin,xf,pf_full,freq_response(steps,nf)--
          endif 
          If S eq 5 or S eq 4 then begin
            print,'Restore: ', 'DataMean_'+save_file_name+'.sav'
            restore,filename='DataMean_'+save_file_name+'.sav'
              ;--xt_rebin,xf,data_mean_save,pf_full
            freq_response = data_mean_save
          endif 
           n_steps = n_elements(freq_response[*,0])
           nf      = n_elements(xf)
           TimeFreqCorrect_array=dblarr(4,nf)
           If S eq 5 or S eq 4 then begin 
            TimeFreqCorrect_array_STDEV=dblarr(4,nf)
            TimeFreqCorrect_array2=dblarr(4,nf)
           endif 
           skip = 0 
         
         ;------------------------------- 
         ;------------------------------
         ;Average mean across all time
         ;-------------------------------
         MeanV = dblarr(nf)
         for i=0,nf-1 do begin
            MeanV(i) = mean(freq_response[*,i])
         endfor 
         cgplot,xf,MeanV,title='V Mean',$
                  xtitle='Frequency (MHz)',ytitle='Intensity',$
                  /xstyle, /ystyle,charsize=1

	  ;--------------------------------------------
          !p.multi=[0,1,2]
          for i=0,nf-1 do begin
            ;fit data
            w           = where(pf_full[*,i] gt 0.90,count)
            yfit        = freq_response[w,i]
            If count ge 3 then begin 
             fit_data    = poly_fit(xt_rebin(w),yfit,2)  ;fit 
	    endif ;gt 3 blocks
            If count lt 3 then begin
              If i eq 0 then begin
                 fit_data=dblarr(3) 
                 fit_data[0] = yfit  
                 fit_data[1] = 0 
                 fit_data[2] = 0
              Endif
              If i gt 0 then begin
                fit_data[0] = TimeFreqCorrect_array(1,i-1)
                fit_data[1] = TimeFreqCorrect_array(2,i-1)
                fit_data[2] = TimeFreqCorrect_array(3,i-1)  
              endif
              skip = 1
            endif  ;lt 3 blocks
            
            ;Save fits
            TimeFreqCorrect_array(0,i)  = xf[i]         ;frequency
            TimeFreqCorrect_array(1,i) = fit_data[0]    ;a
            TimeFreqCorrect_array(2,i) = fit_data[1]    ;b*x
            TimeFreqCorrect_array(3,i) = fit_data[2]    ;c*x^2
 
            ;plot x000
            mod_plot = i mod 1000
            If mod_plot eq 0 and skip eq 0 then begin
              If keyword_set(PS) then begin
                If keyword_set(VERBOSE) then begin
                  print, i
                  print, 'Saved Plot'
                endif
                cgplot,xt_rebin,yfit,title='Time Correction Curve: Freq:'+strtrim(string(cgnumber_formatter(xf[i],decimals=4)),1)+' MHz',$
                  xtitle='Time (s)',ytitle='Frequency-Integrated Values',$
                  /xstyle, /ystyle,charsize=1;,yrange=[1.2e14,1.4e14]
                cgplot,xt_rebin,(fit_data[0] + fit_data[1]*xt_rebin + fit_data[2]*xt_rebin^2.),color='red',/overplot
              endif
            endif

            skip = 0

            If i eq nf-1 then begin
              If keyword_set(PS) then begin
                If keyword_set(VERBOSE) then begin
                  print, i
                  print, 'Saved Plot'
                endif
                cgplot,xt_rebin,yfit,title='Time Correction Curve: Freq:'+strtrim(string(cgnumber_formatter(xf[i],decimals=4)),1)+' MHz',$
                  xtitle='Time (s)',ytitle='Frequency-Integrated Values',/xstyle, /ystyle,charsize=1
                cgplot,xt_rebin,(fit_data[0] + fit_data[1]*xt_rebin + fit_data[2]*xt_rebin^2.),color='red',/overplot
              endif
            endif
          endfor ;loop freq
          
          time = xt_rebin
        endif; end METHOD 2
        
        ;************
        ;Surface
        ;************
        nt_new = n_elements(time)
        surface_fit = dblarr(nt_new,nf)
        
        If S eq 5 or S eq 4 then surface_fit_STDEV = dblarr(nt_new,nf)
        for j=0,nf-1 do begin
          surface_fit[*,j] = TimeFreqCorrect_array(1,j) + TimeFreqCorrect_array(2,j)*time +$
              TimeFreqCorrect_array(3,j)*time^2.
        endfor
        
        If keyword_set(Smooth_Surface) then begin
          surface_fit = smooth(surface_fit,[1,N_smooth])    
        endif
        
        label = '2d Surface From Time-Frequency Fit'
        If keyword_set(VERBOSE) then begin
          print,'Plot 2d-surface'
        endif 
         
        REDUCE_ARRAY, surface_fit, [1,nf], x_t
        REDUCE_ARRAY, surface_fit, [nt_new,1], x_f
     
        If keyword_set(PS) then begin    
          standard_plots,surface_fit,surface_fit,time,xf,label,/NO_ZOOM,xunit = 'secs',panel_ps=1
        endif 
        
        ;-------------------------------
        ;                RFI
        ;-------------------------------
        LE_AUTO_S,x_f,101,3.5,0,x_f,pnet
        wbad = where(pnet eq 0)
        wgood = where(pnet eq 1)
        w_n = n_elements(wbad)
        v = intarr(w_n)
        for i=0,w_n-1 do begin
         v = value_locate(wgood,wbad[i])
         TimeFreqCorrect_array(1,wbad[i]) = TimeFreqCorrect_array(1,wgood[v])
         TimeFreqCorrect_array(2,wbad[i]) = TimeFreqCorrect_array(2,wgood[v])
         TimeFreqCorrect_array(3,wbad[i]) = TimeFreqCorrect_array(3,wgood[v])
        endfor

        for j=0,nf-1 do begin
         surface_fit[*,j] = TimeFreqCorrect_array(1,j) + TimeFreqCorrect_array(2,j)*time[*] +$
           TimeFreqCorrect_array(3,j)*time[*]^2.
        endfor
       
        If keyword_set(PS) then begin    
          standard_plots,surface_fit,surface_fit,time,xf,label,/NO_ZOOM,xunit = 'secs',panel_ps=1
        endif 
       
        If keyword_set(VERBOSE) then begin
          print,'-----------------------------'
          print,'--End Time-Freq Correction---'
          print,'-----------------------------' 
        endif
        
        If keyword_set(PS) then begin
          device,/close
          cgfixps, filename_ps+'_TF.ps'
          cgps2pdf,filename_ps+'_TF.ps',/delete
        Endif
        
        ;Save Time Correction
        save,xt,TimeFreqCorrect_array,filename = save_file_name+'_TimeFreqCorrection.sav'
        If S eq 5 or S eq 4 then begin 
         save,xt,meanV,TimeFreqCorrect_array,filename = save_file_name+'_TimeFreqCorrection.sav'
	 ;save,xt,TimeFreqCorrect_array,TimeFreqCorrect_array2,TimeFreqCorrect_array_STDEV,STDEV_Final,filename = save_file_name+'_TimeFreqCorrection.sav'
        endif  
	If keyword_set(VERBOSE) then begin
          print,'Saved: ', save_file_name+'_TimeFreqCorrection.sav'
        endif
       
        If keyword_set(save_SurfaceFit) then begin
          save,xt_rebin,xf,surface_fit,filename=save_file_name+'_SurfaceFit.sav'
          print,'Saved: ', save_file_name+'_SurfaceFit.sav'
        endif
       
       mem_allTimeCorrect = (MEMORY(/HIGHWATER))/1d9
       PRINT, 'Memory required for all Time-Correction : ', mem_allTimeCorrect ,' Gb'
    
       ;clear heap memory
       HEAP_GC,/PTR,/VERBOSE
     endif ;end time-freq correction 
    
     ;*****************************************
     ;Read in Saved Time-Freq Correction Curve
     ;*****************************************
     If not(keyword_set(TimeFreq_correction)) then begin
         If keyword_set(VERBOSE) then begin
           print,'-------------------------'
           print,'Restoring Time-Freq Correction'
           print,'-------------------------'
         endif 
         
        If keyword_set(LOFAR) or keyword_set(NENUFAR) then begin  
         restore,filename = save_file_name+'_TimeFreqCorrection.sav'
            ;xt,TimeFreqCorrect_array     
         
         reduce_array,xt,[fac],xt_rebin
         xf = TimeFreqCorrect_array(0,*)          ;frequency
         nf = n_elements(xf) 
         ;************
         ;Surface
         ;************
         nt = n_elements(xt_rebin)
         surface_fit = dblarr(n_elements(xt_rebin),nf)
         for j=0,nf-1 do begin
           surface_fit[*,j] = TimeFreqCorrect_array(1,j) + TimeFreqCorrect_array(2,j)*xt_rebin[*] +$
             TimeFreqCorrect_array(3,j)*xt_rebin[*]^2.
         endfor
       
         If keyword_set(save_SurfaceFit) then begin
           save,xt_rebin,xf,surface_fit,filename=save_file_name+'_SurfaceFit.sav'
           print,'Saved: ', save_file_name+'_SurfaceFit.sav'
         endif
    
         If keyword_set(Smooth_Surface) then begin
           surface_fit = smooth(surface_fit,[1,11])
         endif
    
         label = '2d Surface From Time-Frequency Fit'
         ntmax=4500       & nfmax=1200
         nnt=ceil(nt/ntmax) & nnnt=long(nt/nnt)
         nnf=ceil(nf/nfmax) & nnnf=long(nf/nnf)
    
         If keyword_set(VERBOSE) then begin
           print,'Plot 2d-surface'
         endif
      
         If keyword_set(PS) then begin
           set_plot,'PS'
           print,'Filename Ps:', filename_ps
           device,filename=filename_ps+'_R_TF.ps',/landscape
         endif     
    
         REDUCE_ARRAY, surface_fit, [1,nf], x_t
         REDUCE_ARRAY, surface_fit, [nt,1], x_f
         If keyword_set(PS) then begin
           SPDYNPS, rebin(surface_fit(0:nnnt*nnt-1,0:nnnf*nnf-1),nnnt,nnnf), min(xt_rebin),max(xt_rebin),min(xf),max(xf),$
             'Time (sec)','Frequency (MHz)',label,0,0,0,.05,.95,0,'.'
           plot,xt_rebin,x_t,/xsty,xtit='Time (sec)',ytit='Frequency-integrated values',tit=label,$
             yrange=[min(x_t),max(x_t)],/ystyle
           plot,xf,x_f,/xsty,xtit='Frequency (MHz)',ytit='Time-integrated values',tit=label,$
             yrange=[min(x_f),max(x_f)],/ystyle
         endif
         
;         LE_AUTO_S,x_f,101,3.5,0,x_f,pnet
;         wbad = where(pnet eq 0)
;         wgood = where(pnet eq 1)
;         w_n = n_elements(wbad)
;        ; v = intarr(w_n)
;         for i=0,w_n-1 do begin
;           v = value_locate(wgood,wbad[i])
;           TimeFreqCorrect_array(1,wbad[i]) = TimeFreqCorrect_array(1,wgood[v])
;           TimeFreqCorrect_array(2,wbad[i]) = TimeFreqCorrect_array(2,wgood[v])
;           TimeFreqCorrect_array(3,wbad[i]) = TimeFreqCorrect_array(3,wgood[v])
;         endfor
;         
;         for j=0,nf-1 do begin
;           surface_fit[*,j] = TimeFreqCorrect_array(1,j) + TimeFreqCorrect_array(2,j)*xt_rebin[*] +$
;             TimeFreqCorrect_array(3,j)*xt_rebin[*]^2.
;         endfor
;    
;         If keyword_set(Smooth_Surface) then begin
;           surface_fit = smooth(surface_fit,[1,11])
;         endif
;         
;         REDUCE_ARRAY, surface_fit, [1,nf], x_t
;         REDUCE_ARRAY, surface_fit, [nt,1], x_f
;         If keyword_set(PS) then begin
;           SPDYNPS, rebin(surface_fit(0:nnnt*nnt-1,0:nnnf*nnf-1),nnnt,nnnf), min(xt_rebin),max(xt_rebin),min(xf),max(xf),$
;             'Time (sec)','Frequency (MHz)',label,0,0,0,.05,.95,0,'.'
;           plot,xt_rebin,x_t,/xsty,xtit='Time (sec)',ytit='Frequency-integrated values',tit=label,$
;             yrange=[min(x_t),max(x_t)],/ystyle
;           plot,xf,x_f,/xsty,xtit='Frequency (MHz)',ytit='Time-integrated values',tit=label,$
;             yrange=[min(x_f),max(x_f)],/ystyle
;         endif
         
         If keyword_set(PS) then begin
           device,/close 
           cgfixps, filename_ps+'_R_TF.ps'
           cgps2pdf,filename_ps+'_R_TF.ps',/delete
         endif  
         
         mem_allTimeCorrect = (MEMORY(/HIGHWATER))/1d9
         PRINT, 'Memory required for all Time-Correction: ', mem_allTimeCorrect ,' Gb'
      endif ;LOFAR and NENUFAR    
     endif  ;end reading time correction 
     
     print,'-----------------------------------------------------------------------'
     If keyword_set(TimeFreq_correction) then begin
       PRINT, 'Memory required for all Time-Correction: ', mem_allTimeCorrect ,' Gb'
     endif
     
     If keyword_set(Apply_Corrections) then begin
       PRINT, 'Time required for finding Time-Frequency Corrections: ', SYSTIME(1) - st_start ,' s'
     endif  
end  ;end Time-Freq Correction pro
 
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
;%%%%%%%%%%%%%%%%%%%%%%%%% Apply Corrections %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
pro lofar_ApplyCorrections,tmin,N_time,steps,fmin,fmax,beam,S,file_root,date,save_file_name,$
                     disp_values=disp_values,$
                     PS=PS,VERBOSE=VERBOSE,$
                     rebin_freq=rebin_freq, rebin_time=rebin_time,rebin_data=rebin_data,$
                     TimeFreq_correction=TimeFreq_correction,Apply_Corrections=Apply_Corrections,skip_rebinloop=skip_rebinloop,$
                     save_RFIMaskData=save_RFIMaskData,Smooth_Surface=Smooth_Surface,$
                     mask_bit=mask_bit,bit=bit,$
                     use_combinemask=use_combinemask,$
                     Run_JUPITER=Run_JUPITER,REMOVE_BADF=REMOVE_BADF,N_smooth=N_smooth,Process_METHOD2=Process_METHOD2,$
                     LOFAR=LOFAR,NENUFAR=NENUFAR,DYNSPECBF=DYNSPECBF,filename=filename,COMPLEX_NAME=COMPLEX_NAME 
      
      ;start time
      st_start= SYSTIME(1)

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
      
      tmin_or = tmin
      filename_ps = save_file_name+'_processing'
      
      ;----------------------------------------------------------
      ;              Setup Important Info for Data
      ;----------------------------------------------------------
      ;file to read in headers and arrays
      file = file_root+date+'/'+date+'_SAP000_B00'+b_string+'_S'+S_string+'_P000_bf'
      If keyword_set(VERBOSE) then begin
        print, 'File:',file
      endif
     
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
         READ_H5_DATA, file,0,1,fmin,fmax,Beam,S_in,data,time,xf,nt,nf,Time=tf,/VERB,/NSPECTRA,$
          NOF_SAMPLES=NOF_SAMPLES,START_MJD=START_MJD,MJD_sampling_time=MJD_sampling_time,$
          ANTENNA_SET=ANTENNA_SET,SAMPLING_TIME=SAMPLING_TIME,NOF_stations=NOF_stations,$
          channel_width=channel_width
       endif 
       If keyword_set(NENUFAR) then begin
         Beam_in = Beam
         If S eq 0 then Stokes = 0 ;Stokes-I
         If S eq 1 then Stokes = 1 ;Stokes-Q
         If S eq 2 then Stokes = 2 ;Stokes-U
         If S eq 3 then Stokes = 3 ;Stokes-V
         If S eq 4 then Stokes = 4 ;Stokes-V Prime
         If S eq 0 then nstokes = 1 ELSE nstokes = 4
         read_nenufar, file_root+filename, temporary(data2),tt,ff,temporary(Beam_in),ndata, ntt,dt,nf,df,ns,Stokes,nstokes=nstokes,$
           tmin=0,tmax=1,/VERBOSE,CHAN=CHAN,BEAMLETS=BEAMLETS,/REMOVE_BADF,NOF_SAMPLES=NOF_SAMPLES,ALLTIME=tf
         sampling_time = dt
         channel_width = df
       Endif
       If keyword_set(DYNSPECBF) then begin
         If S eq 0 then POL = 1 ;Stokes-I
         If S eq 1 then POL = 2 ;Stokes-Q
         If S eq 2 then POL = 3 ;Stokes-U
         If S eq 3 then POL = 4 ;Stokes-V
         If S eq 4 then POL = 4 ;Stokes-V^2
         If S eq 5 then POL = 4 ;Stokes-V prime
         read_dynspec_fits,file_root+filename+'_beam'+B_string,0,2000,fmin,fmax,POL,data,time,xf,nt,nf,dt=dt,df=df,tall=tf
         sampling_time = dt
         channel_width = df
       Endif
       ;----------------------------------------------
       ;----------------------------------------------
       ;               Setup Time Array
       ;----------------------------------------------
       ;----------------------------------------------
       IF keyword_set(LOFAR) then MJD_t = dindgen(nof_samples,start=START_MJD,INCREMENT=MJD_sampling_time)
       tmin= tmin_or
       
       ;---------------------------------------
       ;      Restore Time-Freq Correction
       ;---------------------------------------
       wn = 0. ;counting polluted channels
       wt = 0. ;counting good channels      
       If S eq 4 or S eq 5 then begin 
        If keyword_set(COMPLEX_NAME) then begin 
          temp = strmid(save_file_name,12,strlen(save_file_name)-5-13) ; e.g. _pol0_41sec_RFI_patrol5.5_pex_lesig3.5_sum
          save_file_name3 = date+'_pol'+strtrim(string(0),1)+temp+'_beam'+b_String
        Endif        
        If not(keyword_set(COMPLEX_NAME)) then save_file_name3 = date+'_pol'+strtrim(string(0),1)+'_beam'+b_String
        print,'Restore I Time Freq Correction: ', save_file_name3+'_TimeFreqCorrection.sav'
        restore,filename=save_file_name3+'_TimeFreqCorrection.sav'
            ;xt,TimeFreqCorrect_array
        TimeFreqCorrect_array_I =  TimeFreqCorrect_array
       endif 
       If keyword_set(LOFAR) or keyword_set(NENUFAR) then begin 
         restore,filename = save_file_name+'_TimeFreqCorrection.sav'
        ;xt,TimeFreqCorrect_array 
        ;xt,MeanV,TimeFreqCorrect_array if S = 4 or 5 
       endif 
       
       ;-----------------------------------
       ;          Restore Mask
       ;-----------------------------------
       If keyword_set(VERBOSE) then begin
         print,'--Restore Mask----'
       endif
       
      ;-------------------------------------
      ;Read I mask for all polarizations (should already be loaded up)
      ;-------------------------------------
      ;If S ge 1 then begin
      ;  temp = strmid(save_file_name,12,strlen(save_file_name)-5-13); e.g. _pol0_41sec_RFI_patrol5.5_pex_lesig3.5_sum
      ;  save_file_name3 = date+'_pol'+strtrim(string(0),1)+temp+'_beam'+b_String
        ;save_create,N_time,steps,patrol_value,le_sig_value,beam,RFI=RFI,patrol=patrol,le_sig=le_sig,sum=sum,$
                    ;pex=pex,run_JUPITER=run_JUPITER,save_file_name=save_file_name_I  ;All RFI  
        ;save_file_name3 = date+'_pol'+strtrim(string(0),1)+'_'+save_file_name_I
     ;  If keyword_set(mask_bit) and keyword_set(bit)  then begin
     ;    print,'Restore Mask: ','mask_full_bit_steps'+save_file_name3+'.sav'
     ;    restore,filename='mask_full_bit_steps_'+save_file_name3+'.sav'   
     ;     ;xt,xf,xsize_step,p2_full_bit_p
     ;    p2_full_bit_p_I = temporary(p2_full_bit_p)
     ;    xsize_step_I    = temporary(xsize_step)
      ; endif
      ;endif 

       ;-------------------
       ;Full Mask in bytes
       ;-------------------
       If not(keyword_set(mask_bit)) then begin
         If not(keyword_set(use_combinemask)) then begin
          print,'Restore Mask: ', 'mask_full_'+save_file_name+'.sav'
           restore,filename='mask_full_'+save_file_name+'.sav'
             ;xt,xf,p2_full
         endif
         ;combined mask in bytes
         If keyword_set(use_combinemask) then begin
           save_file_name3 = strmid(save_file_name,0,strlen(save_file_name)-6); e.g. L429868_pol0_41sec_RFI_patrol5.5_pex_lesig3.5_sum
           print,'Restore Mask: ','mask_combine_'+save_file_name3+'.sav'
           restore,filename='mask_combine_'+save_file_name3+'.sav'
             ;xt,xf,p2_full
         Endif
       endif
       
       ;-----------------------------
       ;Bit mask and convert to byte
       ;-----------------------------
       If keyword_set(mask_bit) and not(keyword_set(bit)) then begin
         If not(keyword_set(use_combinemask)) then begin
             print,'Restore Mask: ','mask_full_bit_'+save_file_name+'.sav'
             restore,filename='mask_full_bit_'+save_file_name+'.sav' ;flag
                ;xt,xf,xsize,p2_full_bit
         endif
         If keyword_set(use_combinemask) then begin
             save_file_name3 = strmid(save_file_name,0,strlen(save_file_name)-6); e.g. L429868_pol0_41sec_RFI_patrol5.5_pex_lesig3.5_sum
             print, 'Restore Mask: ','mask_combine_bit_'+save_file_name3+'.sav'
             restore,filename='mask_combine_bit_'+save_file_name3+'.sav'
         Endif
       endif
       
       ;-----------------------------
       ;Bit mask and work in bits 
       ;-----------------------------
       If keyword_set(mask_bit) and keyword_set(bit) then begin
         If not(keyword_set(use_combinemask)) then begin
           print,'Restore Mask: ','mask_full_bit_steps_'+save_file_name+'.sav'
           restore,filename='mask_full_bit_steps_'+save_file_name+'.sav' ;flag
           ;xt,xf,xsize_step,p2_full_bit_p
         endif  
       Endif
       
       nt = n_elements(xt)
       nf = n_elements(xf)
       
       ;---------------------------------------------------
       ;Mask in bits but don't work in bits (not(bit))
       ;---------------------------------------------------
       If keyword_set(mask_bit) and not(keyword_set(bit)) then begin
         If keyword_set(use_combinemask) then begin
           BITARRAY_TO_BYTEARRAY,p2_full_all_bit,xsize,p2_full
         Endif
         If not(keyword_set(use_combinemask)) then begin
           BITARRAY_TO_BYTEARRAY,p2_full_bit,xsize,p2_full
         Endif
       Endif
      
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
       for kl=0,steps-1 do begin  ;begin time loop-
         tmax =  tmin + N_time  ;set tmin
    
         print,'-----------------------------------------------------------------------'
         print,'reading slice ',kl,' start time: ', tmin, ' end time: ', tmax
         print,'-----------------------------------------------------------------------'
    
         ;------------------------------
         ; Edit Time for selected input
         ;-----------------------------
         u        = where(tf ge tmin and tf lt tmax)
         If keyword_SET(LOFAR) then mjd      = MJD_t(u)
        
         ;----------------------------------------------------------
         ;        Read Data for specific Beam and Polarization
         ;----------------------------------------------------------
         If keyword_set(LOFAR) then begin
           ;file name creation
           file = file_root+date+'/'+date+'_SAP000_B00'+B_string+$
             '_S'+S_string+'_P000_bf'
            READ_H5_DATA, file,tmin,tmax,fmin,fmax,Beam,S_in,data,time,xf,nt,nf,REMOVE_BADF=REMOVE_BADF,UT=UT;,VERB=VERB  
         endif 
         
         If keyword_set(NENUFAR) then begin 
           Beam_in = Beam
           If S eq 0 then Stokes = 0 ;Stokes-I
           If S eq 1 then Stokes = 1 ;Stokes-Q
           If S eq 2 then Stokes = 2 ;Stokes-U
           If S eq 3 then Stokes = 3 ;Stokes-V
           If S eq 4 then Stokes = 4 ;Stokes-V Prime
           If S eq 0 then nstokes = 1 ELSE nstokes = 4
           read_nenufar, file_root+filename,data,time,xf,temporary(Beam_in),ndata, nt,dt,nf,df,ns,Stokes,nstokes=nstokes,$
             tmin=tmin,tmax=tmax,fmax=fmax,fmin=fmin,/REMOVE_BADF,/CORRECT_BANDPASS,/DATA_BANDPASS
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
           
        If S eq 3 then begin 
          ;print,'Processing S=3, V with sign is used!!!'
          print,'Processing S=3, V is used!!!'
          data = data  ;V = |V|
        endif ;s=3
        If S eq 4 or S eq 5 then begin ;V'
          print,' Running Vprime'
          Print,'-------------------'
          Print,'Reading I Data'
          Print,'-------------------'
          If keyword_set(NENUFAR) then begin
            Beam_in = Beam
            Stokes = 0 ;Stokes-I
            nstokes = 1 
            read_nenufar, file_root+filename,data_I,time,xf,temporary(Beam_in),ndata, nt,dt,nf,df,ns,Stokes,nstokes=nstokes,$
              tmin=tmin,tmax=tmax,fmax=fmax,fmin=fmin,/REMOVE_BADF,/CORRECT_BANDPASS,/DATA_BANDPASS
          endif
          If keyword_set(LOFAR) then begin
            S_I_string = strtrim(string(0),1)
            file_I = file_root+date+'/'+date+'_SAP000_B00'+B_string+$
              '_S'+S_I_string+'_P000_bf'
            READ_H5_DATA, file_I,tmin,tmax,fmin,fmax,Beam,0,data_I,t,f,nt,nf,REMOVE_BADF=REMOVE_BADF     
          endif   
           If S eq 4 then data = temporary(data)/data_I ; (V/I) 
           If S eq 5 then data = temporary(data)/data_I  ;V/I
         endif ;s= 4 or 5
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
          
          ;-------------------
          ;Create Mask
          ;-------------------
          p2=bytarr(nt,nf)+1b 
           
          ;-----------------
          ;Restore Byte Mask
          ;------------------
           If not(keyword_set(bit)) then begin
             p2 = p2_full[nt*kl:(kl+1)*nt-1,*]
           endif

          ;---------------------------------
          ;Create Byte Mask from bit mask
          ;----------------------------------
           If keyword_set(mask_bit) and keyword_set(bit) then begin
              pointer_to_array,p2_full_bit_p[kl],out_array=p2_bit ;convert step into bit array
              BITARRAY_TO_BYTEARRAY, p2_bit, xsize_step,p2      ;convert bit array to byte array
           endif  
           
           If kl eq 0 then begin
             If keyword_set(PS) then begin
               label   = 'Data: Raw'
               wbad = where(p2 eq 0)
               wgood = where(p2 eq 1)
               data[wbad] = median(data[wgood]) 
               STANDARD_PLOTS, data,p2*1.,time,xf,label,xunit = 'secs',panel_ps=1
             endif
           endif
  
         ;*************************************
         If kl eq 0 then begin
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
             ;**** New Rebin Freq *******
             Num_rebin_f = long((rebin_freq*1000.)/(channel_width))   ;Number of channels to combine to get the new channel widtth
             mmv_2 = Num_rebin_f
             ;if mmv_2 ne long(mmv_2) then mmv_2=long(mmv_2)+1 else mmv_2=long(mmv_2)
             ;Num_rebin_f = mmv_2
             mmtv_2 = long(nf/mmv_2)
    
             ;**** New Sample Time ****'
             Num_rebin_t = long(rebin_time/sampling_time)             ;Number of bins to combine for time
             nnv_2 = Num_rebin_t
             If Num_rebin_t eq 0 then Num_rebin_t = 1
            ; if nnv_2 ne long(nnv_2) then nnv_2=long(nnv_2)+1 else nnv_2=long(nnv_2)
             ;Num_rebin_t = nnv_2
;             nntv_2=long(nt/nnv_2)
;             nntv_3 = long((nt*steps)/nnv_2)
               ;data(t,f); need this format for SPDYNPS

             ;pointers
             data_rebin_p    = ptrarr(steps)
             p2_rebin_p      = ptrarr(steps)
             xt_rebin_p      = ptrarr(steps)
             mjd_rebin_p     = ptrarr(steps)
             UT_rebin_p      = ptrarr(steps)
             
             
             reduce_array,xf,[Num_rebin_f],xf_rebin,/dbl
            
           endif ;end rebin
         endif ;kl eq 0
        
       ;Apply Time Correction Cir
       If not(keyword_set(Smooth_Surface)) then begin
         If keyword_set(VERBOSE) then begin
          print,'Apply Time-Freq Correction at every freq'
         endif  
        
         If not(keyword_set(gaussian_surface)) then begin
            data = data*p2

            ;---------------
            ;T-F Correction
            ;---------------
            t_array = rebin(time,nt,nf)                          ; (nt,nf)
            If keyword_set(LOFAR) or keyword_set(NENUFAR) then begin  
              If S ne 5 and S ne 4 then begin 
                 print,'S: (divide)'
                 TFC  = rebin(reform(TIMEFREQCORRECT_ARRAY),4,nf,nt)   ; (4,nf,nt)
  	             TFC  = transpose(TFC,[0,2,1])                         ; (4,nt,nf)
  	             surface = reform(TFC[1,*,*] + TFC[2,*,*]*t_array  + TFC[3,*,*]*t_array^(2.0d)) ;(nt,nf)
                 data = data/surface
                 data = data*p2
              endif 
              
         	    If S eq 5 or S eq 4 then begin
                     If S eq 5 and keyword_set(VERBOSE) then Print,'S = 5 (Vprime); subtract '
                     If S eq 4 and keyword_set(VERBOSE) then  Print,'S = 4 (Vprime); subtract '
                     ;v Mean Surface 
                     print,'Subtract by Surface Mean'
                     TFC_v  = rebin(reform(TIMEFREQCORRECT_ARRAY),4,nf,nt)   ; (4,nf,nt)
        	           TFC_v  = transpose(TFC_v,[0,2,1])                         ; (4,nt,nf)
        	           surface_v = reform(TFC_v[1,*,*] + TFC_v[2,*,*]*t_array  + TFC_v[3,*,*]*t_array^(2.0d)) ;(nt,nf)
                      
                     data_prime = data -surface_v  ;v - <v>time
                      
                     ;I 
                     TFC_I  = rebin(reform(TIMEFREQCORRECT_ARRAY_I),4,nf,nt)   ; (4,nf,nt)
        	           TFC_I  = transpose(TFC_I,[0,2,1])                         ; (4,nt,nf)
        	           surface_I = reform(TFC_I[1,*,*] + TFC_I[2,*,*]*t_array  + TFC_I[3,*,*]*t_array^(2.0d)) ;(nt,nf)
                     data_I     = data_I/surface_I ;- 1.0
                      
                     data = data_prime*data_I*p2              ;final Vhat 
        
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
                	  Endif  ;end S eq 4 for run Jupiter 
               
              		data = temporary(data_new)
              		SKIPJUPITER2:
              		print,'********** End Jupiter *********************'  
            	      Endif  ;end Jupiter
                  
                      If kl eq 0 then begin
                       If keyword_set(PS) then begin
                        label   = 'Data: Processing Vprime'
                        STANDARD_PLOTS, data_prime,p2*1.,time,xf,label,/correct_data,xunit = 'secs',panel_ps=1
                        label   = 'Data: Processing I '
                        STANDARD_PLOTS, (data_I-1)*p2,p2*1.,time,xf,label,/correct_data,xunit = 'secs',panel_ps=1
                       endif
                      endif;kl = 0 
        	    endif ;S eq 4 or 5
            endif ;end LOFAR and Nenufar
            
            ;If kl eq 0 then begin 
            ;  label = 'Surface'
            ;  STANDARD_PLOTS,surface,p2*1.,time,xf,label,xunit = 'secs',panel_ps=1
            ;endif ;end kl eq 0

            ;for ii=0,nf-1 do begin
            ; line           = (TimeFreqCorrect_array[1,ii] + TimeFreqCorrect_array[2,ii]*time + TimeFreqCorrect_array[3,ii]*time^2.)
            ; If S ne 5 and S ne 4 then begin 
            ;   print,'S: (divide)'
            ;   data(*,ii)     = data(*,ii)/reform(line,nt,1)
            ; endif   ;

             ;Pol S =4 or 5, v'   
            ;  If S eq 5 or S eq 4 then begin
           ;     If S eq 5 and ii eq 0 then Print,'S = 5 (Vprime); subtract '
            ;    If S eq 4 and ii eq 0 then  Print,'S = 4 (Vprime); subtract '
              ;  STDEV_line           = (TimeFreqCorrect_array_STDEV[1,ii] + TimeFreqCorrect_array_STDEV[2,ii]*time + $
               ;                        TimeFreqCorrect_array_STDEV[3,ii]*time^2.)
             ;   data_prime(*,ii)     = (data(*,ii)- reform(line,nt,1));/reform(STDEV_line,nt,1)
                ;line_I           = (TimeFreqCorrect_array_I[1,ii] + TimeFreqCorrect_array_I[2,ii]*time + TimeFreqCorrect_array_I[3,ii]*time^2.)
               ; data_I(*,ii)     = data_I(*,ii)/reform(line_I,nt,1) - 1.0 
               ; data(*,ii)    = data_prime(*,ii)*data_I(*,ii)*p2
            ;  endif 
            ;endfor ;end loop over nf 
            
          ; If kl eq 0 then begin
          ;  If keyword_set(PS) then begin
          ;    label   = 'Data: Processing'
          ;    STANDARD_PLOTS, data,p2*1.,time,xf,label,/correct_data,xunit = 'secs',panel_ps=1
          ;  endif
          ;endif

            ;surface_fit = dblarr(n_elements(time),nf)
            ;for j=0,nf-1 do begin
            ;  surface_fit[*,j] = TimeFreqCorrect_array(1,j) + TimeFreqCorrect_array(2,j)*time[*] +$
            ;    TimeFreqCorrect_array(3,j)*time[*]^2.
            ;endfor
            ;If kl eq 0 then begin 
            ;  label = 'Surface'
            ;  STANDARD_PLOTS,surface_fit,p2*1.,time,xf,label,xunit = 'secs',panel_ps=1
            ;endif ;end kl eq 0
        endif  ;end not(keyword_set(gaussian_surface))
       
;        ;********************
;        ;find guassian surface
;        ;*******************
;        If keyword_set(gaussian_surface) then begin
;          sigma_t = 11
;          sigma_f = 11
;          lofar_surface,data,p2,time,xf,sigma_t,sigma_f,V_prime=V_prime
;          data = data/V_prime
;        endif 
       Endif  
      
       If keyword_set(Smooth_Surface) then begin
         ;create surface for slice
         ;surface_fit = dblarr(n_elements(time),nf)
         ;for j=0,nf-1 do begin
         ;  surface_fit[*,j] = TimeFreqCorrect_array(1,j) + TimeFreqCorrect_array(2,j)*time[*] +$
         ;    TimeFreqCorrect_array(3,j)*time[*]^2.
         ;endfor
         surface_fit = smooth(surface,[1,N_smooth])
       
          If keyword_set(VERBOSE) then begin
           Print,'Apply Smoothed Surface'
          endif 
          data = data/surface_fit
       endif ;end apply surface 
       
        ;----------------------
        ;set RFI Channels
        ;----------------------
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
              STANDARD_PLOTS, data,p2*1.,time,xf,label,/correct_data,xunit = 'secs',panel_ps=1
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
           If keyword_set(LOFAR) then begin
            reduce_array,mjd,[Num_rebin_t],mjd_rebin2,/dbl
            reduce_array,UT,[Num_rebin_t],UT_rebin2,/dbl
           Endif
           
           ;data_rebin(kl*nntv_2:(kl+1.)*nntv_2-1.,*) = data_rebin2(0:nntv_2-1.,*)  ;update rebin array 
           ;p2_rebin(kl*nntv_2:(kl+1.)*nntv_2-1.,*)   = p2_rebin2(0:nntv_2-1.,*)*1.      ;update p2_rebin array
          
           data_rebin_p[kl]  = ptr_new(data_rebin2)
           p2_rebin_p[kl]    = ptr_new(p2_rebin2)
           xt_rebin_p[kl]    = ptr_new(xt_rebin2)
           If keyword_set(LOFAR) then begin
            mjd_rebin_p[kl]   = ptr_new(mjd_rebin2)
            UT_rebin_p[kl]    = ptr_new(UT_rebin2)
           endif
           
           If kl eq 0 then begin
             If keyword_set(PS) then begin
               label   = 'Data: Processing Rebinned'
               STANDARD_PLOTS, data_rebin2,p2*1.,xt_rebin2,xf_rebin,label,/correct_data,xunit = 'secs',panel_ps=1
             endif
           endif
           
           If keyword_set(PS) then begin
            label = 'Rebinned Data'
            !p.multi=[0,1,2]
            SPDYNPS, data_rebin2/(p2_rebin2 + (p2_rebin2 eq 0)), min(xt_rebin2),max(xt_rebin2),min(xf_rebin),max(xf_rebin),$
                    'Time (sec)','Frequency (MHz)',label,0,0,0,.05,.95,0,'.'      
            label = 'Rebinned Mask'
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
           ;-------------------------------------------
           ;as Kl goes forward,
           ;data chagnes because the starting and ending time change with loop.
           ;So data is a temp variable
    
           ;0 is special case since you don't de_disperse
           If kl eq 0 then begin
             y[nt:nt*2-1,*] = data;[0:nt-1,*]
             p[nt:nt*2-1,*] = p2;[0:nt-1,*]
           Endif
    
           If kl gt 0 then begin
             ;update y and p buffer (two buffers)
             y[0:nt-1,*]           =   y[nt:2*nt-1,*]      ; update the varabile from the last time
             ;y[nt:2*nt-1,*]        =   0.                                    ;reset
             y[nt:2*nt-1,*]        =   data[0:nt-1,*]         ;second block in temp new_data
    
             p[0:nt-1,*]           =   p[nt:2*nt-1,*]      ;second block temp p data set
             ;p[nt:2*nt-1,*]        =   1                                     ;reset
             p[nt:2*nt-1,*]        =   p2[0:nt-1,*]
    
             ;         if keyword_set(make_ps) then begin
             ;           SPDYNPS, data,0,slice,min(xf),max(xf),'#spectrum','MHz',' - slice'+strtrim(fix(kl),2),0,1,0,0.05,0.95,0,'.'
             ;        SPDYNPS, 1.-p2*1.,0,slice,min(xf),max(xf),'#spectrum','MHz',' - mask',0,0,0,0.05,0.95,0,'.'
             ;     endif
    
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
       tmin = tmin + N_time
       
       If keyword_set(Run_JUPITER) then tmin_J = tmin_J + N_time
     endfor ;endfor for time
     
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
        If keyword_set(LOFAR) then begin 
          If j eq 0 then mjd_rebin =  *mjd_rebin_p[j]
          If j gt 0 then mjd_rebin =  [ mjd_rebin,*mjd_rebin_p[j] ]
          If j eq 0 then UT_rebin =  *UT_rebin_p[j]
          If j gt 0 then UT_rebin =  [ UT_rebin,*UT_rebin_p[j] ]
        endif   
      endfor  
      
      ;********
      ;Plot
      ;********
       If keyword_set(PS) and keyword_set(rebin_data) then begin 
         device,/close
         cgfixps, filename_ps+'_rebin.ps'
         cgps2pdf,filename_ps+'_rebin.ps',/delete
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
        STANDARD_PLOTS, data_rebin/(p2_rebin*1. + (p2_rebin*1. eq 0)),p2_rebin*1.,xt_rebin/3600.,xf_rebin,label,/histogram,xunit = 'hour',/ONLY_PLOT,/NO_ZOOM,panel_ps=1
       If keyword_set(RUN_DISP) then begin
         label = 'Data: Rebin (DISP)'
         STANDARD_PLOTS, x/(p*1. + (p*1. eq 0)),p*1.,xt,xf_disp,label,/histogram,xunit = 'secs',panel_ps=1
       endif 
      endif 
     
      If keyword_set(PS) then begin 
        device,/close
        cgfixps, filename_ps+'_final.ps'
        cgps2pdf,filename_ps+'_final.ps',/delete
        EXIT_PS
      endif 
      
      ;*****
      ;Save
      ;*****
      
      ;phase setup 
      If keyword_set(LOFAR) then phase_rebin = (((MJD_REBIN - TO) mod P)/(P*1.d))
      If keyword_set(NENUFAR) then begin
        nt_new = n_elements(rebin_time)
        dt_days = rebin_time/86400.0d   ;rebin_time in seconds 
        JD_rebin = DINDGEN(nt_new,START = JD0, INCREMENT = dt_days)
        MJD_rebin = JD_rebin - 2400000.50d
        phase_rebin = (((MJD_REBIN - TO) mod P)/(P*1.d))
      endif

      save,xt_rebin,MJD_rebin,phase_rebin,UT_rebin,xf_rebin,data_rebin,p2_rebin,filename='rebindata_'+save_file_name+'.sav'
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

pro lofar_processing,tmin,N_time,steps,fmin,fmax,beam,S,file_root,date,save_file_name,$
  disp_values=disp_values,$
  PS=PS,VERBOSE=VERBOSE,$
  rebin_freq=rebin_freq, rebin_time=rebin_time,rebin_data=rebin_data,$
  TimeFreq_correction=TimeFreq_correction,Apply_Corrections=Apply_Corrections,skip_rebinloop=skip_rebinloop,$
  save_RFIMaskData=save_RFIMaskData,Smooth_Surface=Smooth_Surface,$
  mask_bit=mask_bit,bit=bit,$
  use_combinemask=use_combinemask,$
  Run_JUPITER=Run_JUPITER,REMOVE_BADF=REMOVE_BADF,N_smooth=N_smooth,Process_METHOD2=Process_METHOD2,$
  LOFAR=LOFAR,NENUFAR=NENUFAR,DYNSPECBF=DYNSPECBF,filename=filename,COMPLEX_NAME=COMPLEX_NAME,$
  TO=TO,P=P


  ;start time
  stp= SYSTIME(1)
  mem = (MEMORY(/HIGHWATER))/1d9 ; Gb
  
  print,'-------------------------------------'
  print,'---        Running Processing     ---'
  print,'-------------------------------------'
  print,'Date', date
  print,'Beam: ', beam
  print,'Start Frequency (MHz): ', fmin
  print,'End Frequency (MHz): ', fmax
  print,'Polarization', S
  print,'Start Time: ', tmin
  print,'Size of Steps (s): ', N_time
  print,'Steps to take: ', steps
  print,'Total Time (s): ', N_time*steps

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

  If steps eq 1 then disp_values[0] = 0     ;skip disp if steps = 1
  If keyword_set(VERBOSE) then begin
    If steps eq 1 then print, 'Steps equal 1, De-Dispersion not set'
    If disp_values[0] eq 1 then $
      print, 'De-Dispersion set' else $
      print, 'De-Dispersion not set'
      
    if TimeFreq_correction ne 0 then $
      print, 'Time-Freq Correct set' else $
      print, 'Time-Freq Correct not set'

    if keyword_set(Smooth_Surface) then $
      print, 'Smooth_Surface set' else $
      print, 'Smooth_Surface not set'
    
    if keyword_set(use_combinemask) then $
      print, 'use_combinemask set' else $
      print, 'use_combinemask not set'
    if keyword_set(Run_JUPITER) then $
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
         print,'---------------------------------------'
      print, 'Post Process Method 2 (YES = 1; NO = 0)', Process_METHOD2
  endif
   
   print,'------------------------------------------------------'
   print,'Rebin Frequency (KHz): ',rebin_freq
   print,'Rebin Time (s): ',rebin_time
   print,'------------------------------------------------------'
   
   ;Save tmin
   tmin_or = tmin
 
  ;************************************************
  ;       Time-Frequency Corrections
  ;************************************************
  lofar_timefreqcorrection,tmin,N_time,steps,fmin,fmax,beam,S,file_root,date,save_file_name,$
    disp_values=disp_values,$
    PS=PS,VERBOSE=VERBOSE,$
    rebin_freq=rebin_freq, rebin_time=rebin_time,rebin_data=rebin_data,$
    TimeFreq_correction=TimeFreq_correction,Apply_Corrections=Apply_Corrections,skip_rebinloop=skip_rebinloop,$
    save_RFIMaskData=save_RFIMaskData,Smooth_Surface=Smooth_Surface,$
    mask_bit=mask_bit,bit=bit,$
    use_combinemask=use_combinemask,$
    TimeFreqCorrect_array=TimeFreqCorrect_array,$
    Run_JUPITER=Run_JUPITER,REMOVE_BADF=REMOVE_BADF,N_smooth=N_smooth,Process_METHOD2=Process_METHOD2,$
    LOFAR=LOFAR,NENUFAR=NENUFAR,DYNSPECBF=DYNSPECBF,filename=filename,COMPLEX_NAME=COMPLEX_NAME
  
  ;************************************************
  ;      Apply Time-Frequency Correction 
  ;************************************************
  lofar_ApplyCorrections,tmin_or,N_time,steps,fmin,fmax,beam,S,file_root,date,save_file_name,$
    disp_values=disp_values,$
    PS=PS,VERBOSE=VERBOSE,$
    rebin_freq=rebin_freq, rebin_time=rebin_time,rebin_data=rebin_data,$
    TimeFreq_correction=TimeFreq_correction,Apply_Corrections=Apply_Corrections,skip_rebinloop=skip_rebinloop,$
    save_RFIMaskData=save_RFIMaskData,Smooth_Surface=Smooth_Surface,$
    mask_bit=mask_bit,bit=bit,$
    use_combinemask=use_combinemask,$
    Run_JUPITER=Run_JUPITER,REMOVE_BADF=REMOVE_BADF,N_smooth=N_smooth,Process_METHOD2=Process_METHOD2,$
    LOFAR=LOFAR,NENUFAR=NENUFAR,DYNSPECBF=DYNSPECBF,filename=filename,COMPLEX_NAME=COMPLEX_NAME,$
    TO=TO,P=P
   
  print, 'SYSTIME for Processing =',SYSTIME(1) - stp, ' sec'
  print,'-----------------------------------------------------------------------'
  print,'The End of Processing Program'
  return
end 
 
