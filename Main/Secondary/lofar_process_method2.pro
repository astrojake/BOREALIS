;This program is a sub-program for processing method2 with LOFAR
;
;What it does: This program calculates the frequency response of a slice of data after applying the mask. To be used in a loop. 
;              Not advised to use as stand-alone program. 
;              This program exists to allow for Method2 to be consistent between lofar_rfi and lofar_processing
;
;Where it is used: Added in lofar_rfi and lofar_processsing. 
;
;Inputs:
;       data_raw    ; The raw data
;       p2          ; The mask
;       time        ;time of slice
;       xf          ;freq of raw data
;       quantile 
;       sampling_time
;       channel_width
;       patrol_value
;
;Flags:
;       PS
;       verbose 
;
;Outputs:
;       freq_response    ;The frequency response over the time interval
;       p2_f             ;integrated spectrum mask
;       F2               ;;The frequency response over the time interval after S = 4
;Uses: 
;     reduce_array
;     dyn_n
;     Le_Auto_S
;     standard_plots
;***********************************************************************************************************************************
pro lofar_process_method2,data_raw,p2,time,xf,S,quantile,SAMPLING_TIME,channel_width,patrol_value=patrol_value,$
                          freq_response=freq_response,p2_f=p2_f,PS=PS,VERBOSE=VERBOSE,RFI=RFI,IDATA=IDATA,F2=F2,STDEV=stdev
 ;Flags 
 If keyword_set(VERBOSE) then print, 'Quantile for Method 2', quantile
 If not(keyword_set(patrol_value)) then patrol_value  = 4.5
 !p.multi=[0,1,2]
 nt = n_elements(time)
 nf = n_elements(xf)

 ;--------------------------
 ;Mask rebin (bad and good)
 ;-------------------------
 REDUCE_ARRAY, p2*1.0d, [nt,1],p2_f,/dbl   ;mask integrated spectrum
 REDUCE_ARRAY, p2*1.0d, [1,nf],p2_t,/dbl   ;mask time series
 wtime   = where(p2_t ge 0.50)        ;good time
 
 ;-----------
 ;Apply mask
 ;-----------
 data_raw = data_raw*p2

 ;------------------------------------------------------------
 ;factor to fix quantile (note sqrt(2) for two polarizations
 ;------------------------------------------------------------
 If S eq 0 then begin 
  dem = sqrt(2.0d)*sqrt(SAMPLING_TIME*channel_width*1.0d)
  num = gauss_cvf(quantile)*1.0d
  fac = 1.0d - num/dem  ;sqrt(2) comes from polairzation
  If keyword_set(verbose) then print,'Factor of quantile', fac
 endif 
 ;----------------------------------------
 ;Find Freq Response of good data (and fix)
 ;----------------------------------------
 freq_response = dblarr(nf)
 If S eq 0 then begin 
  for iii=0,nf-1 do begin
   freq_response[iii] = dyn_n(data_raw(wtime,iii),quantile)/fac
  endfor
 endif 

 ;---------------------
 ;Fix bad frequencies
 ;---------------------
 wbad = where(p2_f lt 0.50)
 wgood = where(p2_f gt 0.90)
 w_n = n_elements(wbad)
 for i=0,w_n-1 do begin
  v = value_locate(wgood,wbad[i])
  If v eq -1 then v = 0
  freq_response(wbad[i]) = freq_response(wgood[v])
 endfor

 ;---------
 ;More RFI
 ;----------
 If keyword_set(RFI) then begin
  freq_response_in = freq_response - smooth(freq_response,101,/EDGE_TRUNCATE)
  LE_AUTO_S,freq_response_in,31,3.0,0,freq_r,pnet1
  LE_AUTO_S,freq_r,31,3.0,0,freq_r,pnet
  pnet = pnet*pnet1 
   wbad = where(pnet lt 0.50)
   wgood = where(p2_f gt 0.90)
   w_n = n_elements(wbad)
   for i=0,w_n-1 do begin
    v = value_locate(wgood,wbad[i])
    If v eq -1 then v = 0
    freq_response(wbad[i]) = freq_response(wgood[v])
   endfor
 Endif ;end rfi 
 
 If S ne 0 then begin 
   MAKE_BACKGROUND,data_raw*p2,'', data_mean,ss,nn
   LE_AUTO_S,data_mean,101,3.5,0,data_mean,pnet
   data_mean = smooth(data_mean,3,/EDGE_TRUNCATE)  
   freq_response = data_mean 
 endif 

 If keyword_set(PS) then plot, xf, freq_response,xtitle='Frequency (MHz)',ytitle='Frequency Response',title='Cleaned Frequency Response',/ystyle,/xstyle

 ;Plot 
 If keyword_set(PS) then begin
  label   = 'Data (Normalized by method 2)'
  If S eq 0 then data = data_raw/rebin(reform(freq_response[*],1,nf),nt,nf)
  If S eq 4 then data = (data_raw-rebin(reform(freq_response[*],1,nf),nt,nf)) 
  w_good = where(p2 eq 1) 
  data_plot = data - mean(data[w_good])
  STANDARD_PLOTS,data_plot*p2,p2,time,xf,label,xunit = 'secs'
 endif
 
end
