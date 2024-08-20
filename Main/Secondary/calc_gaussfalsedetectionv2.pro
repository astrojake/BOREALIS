;Program: calc_GaussFalseDetection

;This program is edited to only find the maxium value to run millons to trillons of simulations. Not general 

;Gaussian for False Detection calculations:
stp          =  SYSTIME(1)        ;start time

  threshold = DINDGEN(45,start=1.5,   INCREMENT=0.1)    ;1.5 -- 5.9
  n_th = n_elements(threshold)
  N_new  = 1000
 
  folder = '/data/jaket.turner/exoplanet/Gauss/PP_Gauss'
  restore,filename='/Users/jaketurner/Documents/Exoplanet/Codes/LOFAR_Bf_PipelineV3/Main/Secondary/GAUSS/gauss_ref_FalsePostiveRate.sav';'gauss_ref_FalsePostiveRate_10000.sav'
  N_old = N
  res_Q1a_mean = dblarr(N)
  res_Q1b_mean = dblarr(N)
  for i =0,N-1 do begin
    res_Q1a_mean[i] = mean(r_Q1a(*,i))
    res_Q1b_mean[i] = mean(r_Q1b(*,i))
  endfor
 
   res_Q1a_max = max(res_Q1a_mean)
   res_Q1b_max = max(res_Q1b_mean)

 counta = 0
 countb = 0
  ;---------------------------------
  ;   Q1 False Postive Rates
  ;---------------------------------
  dt    = 1.0                    ;time res of rebinned data (1 sec)
  df    = 45.0*1e3               ;freq res of rebinned data (45 kHz or 45000 Hz)
  sigma = 1.0/sqrt(dt*df)        ;sigma on flux; 1/sqrt(b*t)
  nf    = 1080.                  ;The number of freq channels in our data at 45 kHz
  nt    = 1000.                  ;Max number of time steps (should always be less) about 4 hours

  ;Arrays
 ; res_Q1a = fltarr(nt,N)
 ; res_Q1b = fltarr(nf,N)

  for j=0,N_new-1 do begin
    IF (J MOD 100) EQ 0 THEN PRINT,J,' /',N_new
    xdyn = RANDOMN(seed, nt, nf)*sigma  + 1.0
    ydyn = RANDOMN(seed2, nt, nf)*sigma + 1.0

    REDUCE_ARRAY, xdyn, [1,nf], x_t,/dbl  ;time series
    REDUCE_ARRAY, ydyn, [1,nf], y_t,/dbl  ;time series
    REDUCE_ARRAY, xdyn, [nt,1], x_f,/dbl  ;int spectrum
    REDUCE_ARRAY, ydyn, [nt,1], y_f,/dbl  ;int spectrum

    res_Q1a_test= mean(( x_t - y_t)/(s_Q1a))  
    IF (J MOD 100) EQ 0 THEN PRINT, res_Q1a_test
    If res_Q1a_test gt res_Q1a_max then begin
      If counta eq 0 then begin
        N_Save_Q1a = j
        res_Q1a_save = res_Q1a_test
      Endif
      If counta gt 0 then begin
        N_Save_Q1a = [N_Save_Q1a,j]
        res_Q1a_save = [res_Q1a_save,res_Q1a_test]
      Endif
      counta = counta + 1
    Endif
    
    res_Q1b_test= mean(( x_f - y_f)/(s_Q1b > 1.e-5)) 
    IF (J MOD 100) EQ 0 THEN PRINT, res_Q1b_test
    If res_Q1b_test gt res_Q1b_max then begin
      If countb eq 0 then begin
        N_Save_Q1b = j
        res_Q1b_save = res_Q1b_test
      Endif
      If counta gt 0 then begin
        N_Save_Q1b = [N_Save_Q1b,j]
        res_Q1b_save = [res_Q1b_save,res_Q1b_test]
      Endif
      counta = counta + 1
    Endif
  endfor

  ;-------
  ;Time
  ;-------
;  s_Q1a =stddev(res_Q1a(*,*))  ;stdev
;  help, s_Q1a
;  r_Q1a=res_Q1a/rebin(reform(s_Q1a > num,1),nt,N)
;  help, r_Q1a
;
;  ;-------
;  ;Freq
;  ;-------
;  s_Q1b =stddev(res_Q1b(*,*))  ;stdev
;  help, s_Q1b
;  r_Q1b=res_Q1b/rebin(reform(s_Q1b > num,1),nt,N)
;  help, r_Q1b    
    
   save, N_Save_Q1b,res_Q1b_save, N_Save_Q1a,res_Q1a_save, file = 'Gauss_Save_N'+strtrim(string(cgnumber_formatter(N,decimals=0)),1)+'.sav'
   print,'SYSTIME for Simulation =',SYSTIME(1) - stp, ' sec'

end