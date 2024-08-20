;****************************************************************************************************
;                            LOFAR JUPITER: Program to add scaled Jupiter emission
;****************************************************************************************************
;Inputs:   input_file    ;Jupiter Input File (by default called Input_Jupiter.dat) 
;          data
;          xf
;          xt
;          S             ;polarization
;          tmin_J        ;Min time to read for Jupiter data slice
;          tmax_J        ;max time to read for Jupiter data slice
;
;Output    data_new      ;Data with scaled Jupiter emission added
;          factor        ;factor of Jupiter emission to add (optional)
;          data_Jupiter  ;Jupiter Data 
;              
pro lofar_jupitersetup,input_file,data,xf,xt,S,tmin_J,tmax_J,data_new=data_new,$
                  factor=factor,data_Jupiter=Jupiter_Norm,p2_Jupiter=p2_Jupiter,t_J=t_J,freq_J=freq_J,$
                  filename_p2=filename_p2, tf_J=tf_J, nf_J=nf_J,fmin_Jup=fmin_Jup,fmax_Jup=fmax_Jup,target_beam_J=target_beam_J,$
                  A_eff_L=A_eff_L, JA_eff_L_Jup=A_eff_L_Jup, Jupiter_background=Jupiter_background, scale=scale,LOFAR_RUN=LOFAR_RUN,file_Jup=file_Jup,$
                  f_nearest=f_nearest,time_p=time_p,freq_p=freq_p,method2_Jup=method2_Jup,JBackV=Jupiter_background_V,JBackAll=Jupiter_Back_All,VSkyBeam=VSky_Prime
  
 data_new = temporary(data)  ;create data_new

 ;-----------------------
 ;Read Jupiter Input File
 ;-----------------------
 read_jupiter_input,input_file,LOFAR_RUN=LOFAR_RUN,Nancay_RUN=Nancay_RUN,datafile_root=datafile_root,$
   filename=filename_J,target_beam=target_beam_J,RFI_filename=filename_p2,tmin_Jup=tmin_J_main,$
   tmax_Jup=tmax_J_main,fmin_Jup=fmin_Jup,fmax_Jup=fmax_Jup,backfile_root=backfile_root,$
   backfilename=backfilename,fmin_add=fmin_add,scale=scale,method2_Jup=method2_Jup,/VERBOSE

;If keyword_set(LOFAR_RUN) then begin
  file_Jup = datafile_root+filename_J
  READ_H5_DATA, file_Jup,0,1,fmin_Jup,fmax_Jup,target_beam_J,S,Jup_data,time,xf_Jup,nt_test,nf_J,time=Tf_j,$
    /NSPECTRA,/REMOVE_BADF,ANTENNA_SET=ANTENNA_SET,NOF_Stations=NOF_Stations,Freq=xf_Jup_all;,/VERB
    
  ;------------------
  ;  Effective Area
  ;------------------
  lofar_area,ANTENNA_SET,NOF_Stations,xf,A_eff_L=A_eff_L          ;LOFAR OFF Observations
  lofar_area,ANTENNA_SET,NOF_Stations,xf_Jup,A_eff_L=A_eff_L_Jup  ; LOFAR Jupiter Observations

  f_nearest = value_locate(xf,fmin_add)  ;nearest LOFAR frequency to input added frequency
  If f_nearest eq -1 then  f_nearest = value_locate(xf,fmin_Jup+df)
  print,'F_nearest',f_nearest

  ;-----------------------
  ;  Restore Background
  ;-----------------------
  If S eq 0 then begin 
   print,'Restore background file: ', backfile_root+backfilename
   restore,filename=backfile_root+backfilename
   ;restore, background
   v = where(xf_Jup_all ge fmin_Jup and xf_Jup_all le fmax_Jup)
   Jupiter_background  = temporary(background[v])
  endif 
  
  If S eq 3 then begin 
   print,'Restore background file: ', backfile_root+backfilename
   restore,filename=backfile_root+backfilename
   ;restore, background_I,background_V
   Jupiter_Back_All = background_I
   v = where(xf_Jup_all ge fmin_Jup and xf_Jup_all le fmax_Jup)
   Jupiter_background  = temporary(background_I[v])
   Jupiter_background_V  = temporary(background_V[v])
  endif

  ;--------------
  ;Read Mask Info  
  ;--------------
  file_p      = h5f_open(filename_p2)
  dataset_time = h5d_open(file_p,'Time')
  dataset_freq = h5d_open(file_p,'Frequency')
  time_p      = h5d_read(dataset_time)
  freq_p      = h5d_read(dataset_freq)
  
  h5d_close, dataset_time
  h5d_close, dataset_freq
  h5f_close, file_p                  
 ; endif  
end
