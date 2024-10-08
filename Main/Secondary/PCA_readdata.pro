;read all data for PCA
;
st= SYSTIME(1)

tmin = 2.0d  ;spectra
tmin_or = tmin
;tmax = 50
nspectra = 4000 ; 4000
fmin = 14.75
fmax = 62.41
;Beam = 2
S_in = 3
steps = 257   ;257

;read header
READ_H5_DATA,'L570725_SAP000_B000_S0_P000_bf',0,1,fmin,fmax,0,0,$
  data_S,xt_s,xf_s,nt,nf,/VERB,/REMOVE_BADF,Time=xt,Freq=xf,/NSPECTRA,/READ_HEADER, SAMPLING_TIME=SAMPLING_TIME,channel_width=channel_width

N_time = sampling_time*(nspectra)
n_beam = 2
for j=0,n_beam do begin   ;0, 1, 2
  Beam = j
  Beam_string = strtrim(string(Beam),1)
  Pol_string = strtrim(string(S_in),1)

  filename_ps='TEST_PCA_Beam'+Beam_string+'_S'+Pol_string
;  cgPS_Open,filename_ps+'.ps'

  file = 'L570725_SAP000_B00'+Beam_string+'_S'+Pol_string+'_P000_bf'
  data_S_p  = ptrarr(steps)
  data_I_p  = ptrarr(steps)
  tmin = tmin_or
  for i=0,steps-1 do begin
    tmax =  tmin + N_time  ;set tmax
    print, 'step:', i,' Tmin:', tmin, ' Tmax: ',tmax
    READ_H5_DATA, file,tmin,tmax,fmin,fmax,Beam,S_in,data_S,xt,xf,nt,nf,/REMOVE_BADF
    If S_in ne 0 then READ_H5_DATA, 'L570725_SAP000_B00'+Beam_string+'_S0'+'_P000_bf',tmin,tmax,fmin,fmax,Beam,0,data_I,xt,xf,nt,nf,$
                      /REMOVE_BADF,tMJD=MJD,UT=UT
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
  endfor ;steps 
  pointer_to_array,temporary(data_S_p),out_array=data_S
  pointer_to_array,temporary(data_I_p),out_array=data_I
  nf = n_elements(data_I[0,*])
  nt = n_elements(data_I[*,0])
  If S_in ne 0 then data_S = temporary(data_S)/data_I  ;V, Q, U
  p_V=bytarr(nt,nf)+1b
  p_I=bytarr(nt,nf)+1b
  p = p_V*p_I 
  
  If S_in eq 0 then begin 
    ;--------------------------------------------
    ;                     I
    ;--------------------------------------------
    data_gain = fltarr(nf)
    for iii=0,nf-1 do data_gain[iii] = dyn_n(data_I(*,iii),quantile)
    data_gain = rebin(reform(data_gain,1,nf),nt,nf)
    data_I = temporary(data_I)/temporary(data_gain)
    i_nan = where(~finite(data_I), /null)
    data_I[i_nan] = 0.
    p[i_nan] = 0
;    STANDARD_PLOTS, data_I,p,xt,xf,'I: Pre-Sysrem',/ONLY_PLOT,/NO_ZOOM
    data_I = temporary(data_I)*p
  endif
  
  If S_in ne 0 then begin 
    ;---------------------------------
    ;                V
    ;----------------------------------
    reduce_array,data_S*p,[nt,1],data_avg
    avg = rebin(reform(data_avg,1,nf),nt,nf)
    data_S =  temporary(data_S) - temporary(avg)
;    STANDARD_PLOTS, data_S,p,xt,xf,'S: Pre-Sysrem',/ONLY_PLOT,/NO_ZOOM
    help, data_S
    help, nf,nt
    data_S = temporary(data_S)*p
  endif 
  If S_in eq 0 then save,data_I,p,xt_all,xf,UT_all,MJD_All,filename=filename_ps+'.sav'
  If S_in ne 0 then save,data_S,p,xt_all,xf,UT_all,MJD_All,filename=filename_ps+'.sav'
end  ;beam  

end 
