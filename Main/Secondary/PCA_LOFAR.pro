;PCA of data 
st= SYSTIME(1)
startmem =  MEMORY(/CURRENT)

tmin = 2.0d  ;spectra 
nspectra = 4000 ; 4000
fmin = 14.75
fmax = 62.41
;;;Beam = 2
S_in = 3
steps = 257   ;257 

;;read header 
READ_H5_DATA,'L570725_SAP000_B000_S0_P000_bf',0,1,fmin,fmax,0,0,$
  data_n,xt_s,xf_s,nt,nf,/VERB,/REMOVE_BADF,Time=xt,Freq=xf,/NSPECTRA,/READ_HEADER, SAMPLING_TIME=SAMPLING_TIME,channel_width=channel_width

SAMPLING_TIME  = 0.010485760d ;secs
channel_width  = 3051.7578    ;Hz


;N_time = sampling_time*(nspectra)
n_beam = 1
for j=n_beam,n_beam do begin   ;0, 1, 2
  Beam = j 
  Beam_string = strtrim(string(Beam),1)
  Pol_string = strtrim(string(S_in),1)
;  
  filename_ps='TEST_PCA_Beam'+Beam_string+'_S'+Pol_string
 ;cgPS_Open,filename_ps+'.ps'
;  
;  file = 'L570725_SAP000_B00'+Beam_string+'_S'+Pol_string+'_P000_bf'
;  data_S_p  = ptrarr(steps)
;  data_I_p  = ptrarr(steps)
;
  for i=0,steps-1 do begin 
   tmax =  tmin + N_time  ;set tmax
   print, 'step:', i,' Tmin:', tmin, ' Tmax: ',tmax
   READ_H5_DATA, file,tmin,tmax,fmin,fmax,Beam,S_in,data_S,xt,xf,nt,nf,/REMOVE_BADF
   READ_H5_DATA, 'L570725_SAP000_B00'+Beam_string+'_S0'+'_P000_bf',tmin,tmax,fmin,fmax,Beam,0,data_I,xt,xf,nt,nf,/REMOVE_BADF,$
    tMJD=MJD,UT=UT
    data_S_p[i]    = ptr_new(data_S)
    data_I_p[i]    = ptr_new(data_I)
    print, 'Tmin:', tmin, ' Tmax: ',tmax
    
    If i eq 0 then xt_all         = xt
    If i gt 0 then xt_all         = [xt_all,xt]
    If i eq 0 then UT_all        = UT
    If i gt 0 then UT_all        = [UT_all,UT]
    If i eq 0 then MJD_all        = MJD
    If i gt 0 then MJD_all        = [MJD_all,MJD]
    ;****************
    ;update time loop
    ;****************
    tmin = tmin + N_time
  endfor
;  pointer_to_array,temporary(data_S_p),out_array=data_S
;  pointer_to_array,temporary(data_I_p),out_array=data_I
; nf = n_elements(data_I[0,*])
; nt = n_elements(data_I[*,0])
;  If S_in ne 0 then data_S = temporary(data_S)/data_I  ;V, Q, U
;
;
;
;

;------- Read Data ---  
restore,filename=filename_ps+'.sav'

nt = n_elements(xt_all)
nf=  n_elements(xf) 


  error_bar = 1.0d/sqrt(SAMPLING_TIME*channel_width)
  error = dblarr(nt,nf) + error_Bar


  ;---------------------------------
  ;                V
  ;----------------------------------
  reduce_array,data_S*p,[nt,1],data_avg
  avg = rebin(reform(data_avg,1,nf),nt,nf)
  data_S =  data_S - avg
 ; STANDARD_PLOTS, data_S,p,xt,xf,'S: Pre-Sysrem',/ONLY_PLOT,/NO_ZOOM
  help, data_S
  help, nf,nt
  data_S = data_S*p
  for i=0,1 do begin 
    print,'I: ',i
    sysrem,data_S,error,sys_err=sys_err
    data_S = data_S - sys_err  
  endfor
 
  ;Rebin to save
  Num_rebin_t = long(1.0d/sampling_time)
  Num_rebin_f = long((45.0d*1000.)/(channel_width))
  reduce_array,data_S,[num_rebin_t, num_rebin_f], data_S
  reduce_array,xt_all,[num_rebin_t],xt_all
  reduce_array,UT_all,[num_rebin_t],UT_all
  reduce_array,mjd_all,[num_rebin_t],mjd_all
  reduce_array,xf,[num_rebin_f],xf

  save,data_S,xt_all,xf,UT_all,MJD_All,filename=filename_ps+'PCA_rebins.sav'

 
  HEAP_GC,/PTR,/VERBOSE 
endfor ;beams  
print,'SYSTIME entire PCA program=',SYSTIME(1) - st, ' sec'
PRINT, 'Memory required for PCA (1 beam): ', (MEMORY(/HIGHWATER))/1d9 ,' Gb'

end
