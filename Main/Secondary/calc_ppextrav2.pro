;calculations for lofar post-process numbers
;input: x     (x ramp; threshold values for Q4 or time for Q1,/Flags determine)
;       y     (y ramp; ON-OFF for Q4 or ON-Beam for Q1) 
;       gauss (gaussian curve) 
;       
;Flags:
;       Q1          ; Run Q1 (need to input y2)
;       RunQ4       ; Run Q4 
;       y2          ; OFF curve for Q1  
;       gauss_file  ;Folder of LOFAR pipeline (e.g. /data/jake.turner/ where LOFAR_Bf_PipelineV3/ is located) 
;                     ;if not set then location is code on NANCEP2
;       MINT        ; Min threshold vallue used for false-postiive calc on Q4
;       MAXT        ; Max threshold value used for false-postiive calc on Q4
;         
;Output: 
;       out; numerical values for Q1 or Q4 
;Version; 1 Sept 2018
;         2 Oct 4 2018 Updated with probability of false detection     
pro calc_ppextrav2,x,y,gauss=gauss,out=out,RunQ1=Q1,Q1a=Q1a,Q1b=Q1b,RunQ4=Q4,Q4a=Q4a,Q4b=Q4b,Q4c=Q4c,Q4d=Q4d,Q4e=Q4e,Q4f=Q4f,y2=y2,file_g=file_g,MINT=MINT,MAXT=MAXT
  n = n_elements(x) 
  ny= n_elements(y) 
  If n ne ny then begin
   print,'Nx and Ny not equal!'
   stop 
  endif  

If keyword_set(Q4) then begin
  ;Out array
  nn = 13
  out = dblarr(nn) 
  
  ;-------------------------
  ;Integral of ON-OFF curve
  ;-------------------------
  out[0]    = INT_TABULATED(x,y)    
  
  ;-------------------------------------
  ;Ratio of ON-OFF Integral to Gaussian 
  ;-------------------------------------
  int_gauss1 = INT_TABULATED(x,gauss)              ;Integral of Gaussian
  int_gauss2 = INT_TABULATED(x,(gauss*2.0))
  int_gauss3 = INT_TABULATED(x,(gauss*3.0))
  out[1]     = float(out[0])/float(int_gauss1)     ;ratio
  out[2]     = float(out[0])/float(int_gauss2)
  out[3]     = float(out[0])/float(int_gauss3)

  ;--------------------------------
  ;Fraction of ON-OFF values >0 
  ;--------------------------------
  wabove = where(y gt 0,count_above)
  out[4] = (float(count_above)/float(ny))*100.0
 
  ;------------------------------------
  ;Number of values above 1,2, 3 sigma 
  ;------------------------------------
  w_1 = where(y gt gauss,count_above1)    
  w_2 = where(y gt gauss*2.0,count_above2)
  w_3 = where(y gt gauss*3.0,count_above3)
  out[5] = count_above1
  out[6] = count_above2
  out[7] = count_above3

  ;----------------------------------
  ;min/max threshold where ON-OFF >0
  ;----------------------------------
  out[8] = min(x[wabove])
  out[9] = max(x[wabove])

  ;--------------------------------------------------
  ;% of y values > 0 between min and max threshold
  ;--------------------------------------------------
  w_b = where(x[wabove] gt out[8] and x[wabove] lt  out[9],count_b_above) 
  out[10] = (float(count_b_above)/float(ny))*100. 
  
  ;------------------------
  ;% of values gt tau(max)
  ;-----------------------
  w_max = where(x gt out[9],count_allabove)        ;select all values above max threhold
  w_abovemax = where(y[w_max] lt 0,count_abovemax) ;
  out[11] = (float(count_abovemax)/float(n))*100.  ;y above max threshold and non-zero
  
  ;1 sigma
  test_detectionv2,y,x,gauss,PROP_FALSE=PROP_FALSE1,file_g=file_g,/RunQ4,Q4a=Q4a,Q4b=Q4b,Q4c=Q4c,Q4d=Q4d,Q4e=Q4e,Q4f=Q4f,/VERBOSE,MINT=MINT,MAXT=MAXT
  out[12] = PROP_FALSE1  ;Probability of Flase Detection
endif 

If keyword_set(Q1) then begin
  ;inputs: x = time, y  = ON-Beam, y2 = OFF-Beam
  If not(keyword_set(y)) then print,'OFF beam not input!! ERRROR!'
  out = fltarr(4) 
  diff = y - y2
  out[0] = INT_TABULATED(x,diff)               ;integral of ON-OFF
  wabove = where(y gt y2,count_above)        
  out[1] = (float(count_above)/float(n_elements(x)))*100. ;% of time where ON > OFF
  out[2] = stddev(diff)                     ;sddev of ON-OFF    
  test_detectionv2,diff,x,gauss,PROP_FALSE=PROP_FALSE1,file_g=file_g,/RunQ1,Q1a=Q1a,Q1b=Q1b,/VERBOSE
  out[3] = PROP_FALSE1
endif 
end
