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
pro lofar_jupiter,input_file,data,xf,xt,S,tmin_J,tmax_J,data_new=data_new,$
                  factor=factor,data_Jupiter=Jupiter_Norm,p2_Jupiter=p2_Jupiter,t_J=t_J,freq_J=freq_J,$
                  filename_p2=filename_p2, tf_J=tf_J, nf_J=nf_J,fmin_Jup=fmin_Jup,fmax_Jup=fmax_Jup,target_beam_J=target_beam_J,$
                  A_eff_L=A_eff_L, JA_eff_L=A_eff_L_Jup, Jupiter_background=Jupiter_background, scale=scale,LOFAR_RUN=LOFAR_RUN,$
                  file_Jup=file_Jup,f_nearest=f_nearest,time_p=time_p,freq_p=freq_p,method2=method2,$
                  UJupiter=data_Jupiter_U,QJupiter=data_Jupiter_Q,JBackV=Jupiter_background_V,JBackAll=Jupiter_Back_All;,VSkyBeam=VSkyBeam
 
 smooth_Jupiter  = 0  ; don't smooth Jupiter
 LOFAR_RUN = 1        ;always run Lofar
  
 If S lt 4 then S_in = S   ; I, Q, U, abs(V)
 If S eq 4 then S_in = 3   ; V^2
 If S eq 5 then S_in = 3   ; V'
 If S eq 6 then begin     ; L = sqrt(Q^2 + U^2)  
   S_in = 1 
   Sin1 = 1              ;Q
   Sin2 = 2              ;V
 endif 
 print, 'S In',S_in 
 print, 'File_Jup',File_Jup
 data_new = temporary(data)  ;create data_new 
 If keyword_set(LOFAR_RUN) then begin  
   ;-----------------
   ;Read Jupiter
   ;-----------------
   READ_H5_DATA, file_Jup,tmin_J,tmax_J,fmin_Jup,fmax_Jup,target_beam_J,S_in,data_Jupiter,t_J,freq_J,nt_J,nf_J,$
     /REMOVE_BADF,Freq=center_freq,SAMPLING_TIME=SAMPLING_TIME,channel_width=channel_width
   print, '------S------- ' ,S
   If S eq 3 then begin
      data_Jupiter = abs(data_Jupiter)
   endif
   If S eq 4 then begin 
      file_JupI = strmid(file_Jup,0,73) ;/data/jake.turner/exoplanet/LC7_013/Jupiter/L568467/L568467_SAP000_B000_S
      file_JupI = file_JupI+'0_P000_bf' ;I
      READ_H5_DATA, file_JupI,tmin_J,tmax_J,fmin_Jup,fmax_Jup,target_beam_J,0,data_Jupiter_I,t_J,freq_J,nt_J,nf_J,$
       /REMOVE_BADF
 
     ; file_JupI = strmid(file_Jup,0,70)
     ; file_JupI_2 = file_JupI+'2_S0_P000_bf' ;I & beam 2 
     ; READ_H5_DATA, file_JupI_2,tmin_J,tmax_J,fmin_Jup,fmax_Jup,2,0,data_B2,t_J,freq_J,nt_J,nf_J,$
     ;  /REMOVE_BADF
     ; MAKE_BACKGROUND,data_B2,'', data_back_I,ss,nn       ;Background for I data (to Correct)
     data_Jupiter_V = data_Jupiter
     data_Jupiter   = data_Jupiter_V/data_Jupiter_I   ;V/I

     w = where(center_freq gt fmin_Jup and center_freq le fmax_Jup)
     ;Vskybeam = VskyBeam(w)
     Jupiter_background_V = Jupiter_background_V(w)
     Jupiter_background_I   =Jupiter_background(w)
   endif 
   If S eq 6 then begin
      data_Jupiter_Q = data_Jupiter
      ;file_Jup = /data/jake.turner/exoplanet/LC7_013/Jupiter/L568467/L568467_SAP000_B000_S0_P000_bf
      file_JupU = strmid(file_Jup,0,73) ;/data/jake.turner/exoplanet/LC7_013/Jupiter/L568467/L568467_SAP000_B000_S
      file_JupU = file_JupU+'2_P000_bf' ;U
      READ_H5_DATA, file_JupU,tmin_J,tmax_J,fmin_Jup,fmax_Jup,target_beam_J,2,data_Jupiter_U,t_J,freq_J,nt_J,nf_J,$
        /REMOVE_BADF
      print,'**Create L Pol**'
      data_Jupiter = sqrt(data_Jupiter_Q^(2.0d) +data_Jupiter_U^(2.0d)) ;L = sqrt(U^2 + Q^2) 
     
      ; STANDARD_PLOTS, data_Jupiter,data_Jupiter,t_J,freq_J,'Jupiter Norm L',xunit = 'secs'
      print,'Min & Max: ' ,min(data_Jupiter),max(data_Jupiter)
   endif

   ;--------------------
   ;    Read Mask 
   ;--------------------
   file_idJ_p2      = h5f_open(filename_p2)
   dataset_idJ_p2   = h5d_open(file_idJ_p2,'Mask')
   dataspace_idJ_p2 = H5D_GET_SPACE(dataset_idJ_p2)  
   u_J        = where(time_p ge tmin_J and time_p lt tmax_J)
   w_J        = where(freq_p ge fmin_Jup and freq_p le fmax_Jup) 
   nt_Jp = n_elements(u_J)
   nf_Jp = n_elements(w_J)
   startJ_p2 = [min(u_J),min(w_J)]
   countJ_p2 = [nt_Jp,nf_Jp]
   H5S_SELECT_HYPERSLAB, dataspace_idJ_p2, startJ_p2, countJ_p2, /RESET
   memory_space_idJ_p2 = H5S_CREATE_SIMPLE(countJ_p2)
   p2_Jupiter          = H5D_READ(dataset_idJ_p2, FILE_SPACE=dataspace_idJ_p2,MEMORY_SPACE=memory_space_idJ_p2)
   h5s_close, dataspace_idJ_p2
   H5S_CLOSE, memory_space_idJ_p2
   h5d_close, dataset_idJ_p2
   h5f_close, file_idJ_p2
   
   ;----------
   ;Apply Mask
   ;--------- 
   data_Jupiter = data_Jupiter*p2_Jupiter
   
   ;-------------------------------
   ;Smooth Jupiter Data ;Flag!!!!!!
   ;-------------------------------
   If keyword_set(smooth_Jupiter) then begin
     print,'Alert! Smooth Jupiter!!!'
    w_bad  = where(p2_Jupiter eq 0)
    w_good = where(p2_Jupiter eq 1)
    data_Jupiter[w_bad] = median(data_Jupiter[w_good])
    data_Jupiter = smooth(data_Jupiter,[953,1],/EDGE_TRUNCATE)  ;10 second smoothing with 0.01 msec resolution
    data_Jupiter = data_Jupiter*p2_Jupiter
   endif 

   ;---------------------
   ;Create Scaled Jupiter
   ;--------------------- 
   If S eq 0 then begin 
     factor       = data_Jupiter    ;create factor
     Jupiter_Norm = data_Jupiter    ;create Normalization of Jupiter
    for j=0,nt_J-1L do begin  
    for i=0,nf_J-1L do begin 
       If keyword_set(LOFAR_RUN) then begin
        ;A_fac       = A_eff_L[i+f_nearest]/A_eff_L_Jup[i]          ;area factor 
        ;freq_fac    = (xf[i+f_nearest]/freq_J[i])^2.55d            ;frequency factor
        Jupiter_Norm[j,i]  = data_Jupiter[j,i]/Jupiter_background[i]   ;Jupiter normalized to the background
        factor[j,i] = 1.0/(scale)*(Jupiter_Norm[j,i])*(30.0d/freq_J[i])^2.0d > 0.0d  ;30 MHz 
        ; factor[j,i] = 1.0/(scale)*temporary(A_fac)*temporary(freq_fac)*(Jupiter_Norm[j,i] - 1.0d) > 0.0d   ;complete factor for scaling
       Endif
      data_new[j,i+f_nearest] = data_new[j,i+f_nearest]*(1.0 +  factor[j,i])    ;add Jupiter to the data
    endfor
    endfor
  endif  ;S = 0 adding Jupiter
  
  If (S eq 4) or (S eq 5) or (S eq 6) then begin   ;Polarization 
     ;data_ICor = data_Jupiter_I/rebin(reform(Jupiter_background,1,nf_J),nt_J,nf_J)    ;I_{corrrected}      
    ; data_Vprime = (data_Jupiter -rebin(reform(Jupiter_background_V,1,nf_J),nt_J,nf_J))*data_Jupiter_I   ;Vprime =  [vJ(f) - <vS1(f)>]*IJupiter
     ;data_Jupiter=  data_Vprime - rebin(reform(VSkyBeam,1,nf_J),nt_J,nf_J)
     ;data_Jupiter = data_Jupiter*data_ICor*p2_Jupiter   ;V'J1 
 
     ;data_Jupiter ; v = V/I (Jupiter)
     ;Jupiter_background_V = Average sky background in V 
      data_Vprime = data_Jupiter -rebin(reform(Jupiter_background_V,1,nf_J),nt_J,nf_J) ; v1a - <v1b>
      Idiff       = (data_Jupiter_I - rebin(reform(Jupiter_background_I,1,nf_J),nt_J,nf_J))/rebin(reform(Jupiter_background_I,1,nf_J),nt_J,nf_J)
      data_Jupiter = data_Vprime*Idiff      

     alpha      = 1.0d/scale
     freq_fac   =dblarr(nf_J)
     back_fac   =dblarr(nf_J)      ;Ratio of Stokes-I backgrounds
     factor     = data_Jupiter    ;create factor
     for j=0,nt_J-1L do begin  
     for i=0,nf_J-1L do begin
       back_fac[i]                = Jupiter_Back_All[i+f_nearest]/Jupiter_Back_All[i]    ;Ratio of Stokes-I backgrounds
       freq_fac[i]                = (30.0d/freq_J[i])^2.0d				     ;Ratio of SEFDs 
       factor[j,i]                = alpha*freq_fac[i]*back_fac[i]*data_Jupiter[j,i]
       data_new[j,i+f_nearest]    = data_new[j,i+f_nearest] + factor[j,i]
     endfor
     endfor  
     Jupiter_Norm = data_Jupiter
  endif  ; pol 

Endif   ;lofar
 
return
end
