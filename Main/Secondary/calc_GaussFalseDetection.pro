;Program Name: calc_GaussFalseDetection
;Purpose: Calculate Gaussian simulations for the false postive rate for the false postive rate

;Output: gauss_ref_FalsePostiveRate_'+N_string+'.sav'

;Notes: Threshold goes from 1 to 5.9 (like default pipeline)
;
;Version 1: Oct 16
;Version 2: Oct 31 ;includes Q1 false postive rates
;---------------------------------------------------------------------------------
stp          =  SYSTIME(1)        ;start time

;-------
;Steps
;-------
;N = 100.
;
;N = 500 000.0d
N = 1000000.0d
print,N
Nsim = 10000
N_string = STRTRIM(string(N), 1) 

;---------------------------------
;   Q4 False Postive Rates
;---------------------------------
  threshold = DINDGEN(50,start=1.0,   INCREMENT=0.1)    ;1.0 -- 5.9 threholds
  n_th = n_elements(threshold)
  ;array setup for resiudals and stdev 
  res_Q4a=fltarr(n_th,N)  & res_Q4b=fltarr(n_th,N)
  res_Q4c=fltarr(n_th,N)  & res_Q4d=fltarr(n_th,N)
  res_Q4e=fltarr(n_th,N)  & res_Q4f=fltarr(n_th,N)
  s_Q4a=fltarr(n_th)      & s_Q4b=fltarr(n_th) 
  s_Q4c=fltarr(n_th)      & s_Q4d=fltarr(n_th)
  s_Q4e=fltarr(n_th)      & s_Q4f=fltarr(n_th) 

  for j=0,N-1 do begin
    x=randomn(seed,Nsim) & y=randomn(seed2,Nsim) &   z=randomn(seed3,Nsim)  ;random arrays
     IF (J MOD 10000) EQ 0 THEN PRINT,J,' / N: ',N
    for i=0,n_th-1 do begin
      
      postprocessing_calc,x,y,threshold[i],slope=2,$  ;x vs y
        NxA=NxA,PNxA=Power_NxA,$
        Nx2A=Nx2A,PNx2A=Power_Nx2A,$
        Nx3Ae=NxA_3e,PNAx3f=Power_NxA_3f,$
        Nx3Be=NxB_3e,PNBx3f=Power_NxB_3f,$
        NxB=NxB,PNxB=Power_NxB,$
        Nx2B=Nx2B,PN2xB=Power_Nx2B
        
      postprocessing_calc,z,y,threshold[i],slope=2,$  ;z vs y
        NxA=NxA_2,PNxA=Power_NxA_2,$
        Nx2A=Nx2A_2,PNx2A=Power_Nx2A_2,$
        Nx3Ae=NxA_3e_2,PNAx3f=Power_NxA_3f_2,$
        Nx3Be=NxB_3e_2,PNBx3f=Power_NxB_3f_2,$
        NxB=NxB_2,PNxB=Power_NxB_2,$
        Nx2B=Nx2B_2,PN2xB=Power_Nx2B_2
      
        res_Q4a(i,j)=(NxA-NxA_2)                                                ;Q4a (ON-OFF)
        res_Q4b(i,j)=(Power_NxA-Power_NxA_2)                                    ;Q4b (ON-OFF)
        res_Q4c(i,j)=((NXA  - NxB) - (NXA_2  - NxB_2))                          ;Q4c (ON-OFF)
        res_Q4d(i,j)=((Power_NXA  - Power_NxB) - (Power_NXA_2  - Power_NxB_2))  ;Q4d (ON-OFF)
        res_Q4e(i,j)=(NxA_3e-NxA_3e_2)                                          ;Q4e (ON-OFF)
        res_Q4f(i,j)=(Power_NxA_3f-Power_NxA_3f_2)                              ;Q4f (ON - OFF)
      endfor ;thresholds
  endfor
 
  ;Stdev 
  for i=0,n_th-1 do begin 
    s_Q4a(i)=stddev(res_Q4a(i,*))
    s_Q4b(i)=stddev(res_Q4b(i,*))
    s_Q4c(i)=stddev(res_Q4c(i,*))
    s_Q4d(i)=stddev(res_Q4d(i,*))
    s_Q4e(i)=stddev(res_Q4e(i,*))
    s_Q4f(i)=stddev(res_Q4f(i,*))
  endfor  
  
  ;num = 1.0d/N
  num = 1.e-5
  print, 'Num', num
  ;r_Q4a=res_Q4a/rebin(reform(s_Q4a > (1.e-5),n_th,1),n_th,N)
  r_Q4a=res_Q4a/rebin(reform(s_Q4a > num,n_th,1),n_th,N)
  r_Q4b=res_Q4b/rebin(reform(s_Q4b > num,n_th,1),n_th,N)
  r_Q4c=res_Q4c/rebin(reform(s_Q4c > num,n_th,1),n_th,N)
  r_Q4d=res_Q4d/rebin(reform(s_Q4d > num,n_th,1),n_th,N)
  r_Q4e=res_Q4e/rebin(reform(s_Q4e > num,n_th,1),n_th,N)
  r_Q4f=res_Q4f/rebin(reform(s_Q4f > num,n_th,1),n_th,N)
  
  threshold_g = threshold

;---------------------------------
;   Q1 False Postive Rates
;---------------------------------
dt    = 1.0                    ;time res of rebinned data (1 sec)
df    = 45.0*1e3               ;freq res of rebinned data (45 kHz or 45000 Hz)
sigma = 1.0/sqrt(dt*df)        ;sigma on flux; 1/sqrt(b*t)
nf    = 1080.                  ;The number of freq channels in our data at 45 kHz
nt    = 1000.                  ;Max number of time steps (should always be less) about 4 hours

;Arrays
res_Q1a = fltarr(nt,N)
res_Q1b = fltarr(nf,N)

for j=0,N-1 do begin
   IF (J MOD 1000) EQ 0 THEN PRINT,J,' / N'
  xdyn = RANDOMN(seed, nt, nf)*sigma  + 1.0
  ydyn = RANDOMN(seed2, nt, nf)*sigma + 1.0

  REDUCE_ARRAY, xdyn, [1,nf], x_t,/dbl  ;time series
  REDUCE_ARRAY, ydyn, [1,nf], y_t,/dbl  ;time series
  REDUCE_ARRAY, xdyn, [nt,1], x_f,/dbl  ;int spectrum
  REDUCE_ARRAY, ydyn, [nt,1], y_f,/dbl  ;int spectrum 

  res_Q1a(*,j) = x_t - y_t 
  res_Q1b(*,j) = x_f - y_f 
endfor

  ;-------
  ;Time 
  ;-------
  s_Q1a =stddev(res_Q1a(*,*))  ;stdev 
  help, s_Q1a
  r_Q1a=res_Q1a/rebin(reform(s_Q1a > num,1),nt,N)
  help, r_Q1a

  ;-------
  ;Freq 
  ;-------
  s_Q1b =stddev(res_Q1b(*,*))  ;stdev 
  help, s_Q1b
  r_Q1b=res_Q1b/rebin(reform(s_Q1b > num,1),nt,N)
  help, r_Q1b

;----------------------------------------------------------------------------------------------
;sav file
;----------------------------------------------------------------------------------------------
 save,res_Q4a,res_Q4b,res_Q4c,res_Q4d,res_Q4e,res_Q4f,s_Q4a,s_Q4b,$
      s_Q4c,s_Q4d,s_Q4e,s_Q4f,r_Q4a,r_Q4b,r_Q4c,r_Q4d,r_Q4e,r_Q4f,r_Q1a,r_Q1b,s_Q1a,s_Q1b,N,res_Q1a,res_Q1b,threshold_g,file='gauss_ref_FalsePostiveRate_'+N_string+'.sav'

print,'Number of Gauss Sims: ', N
print,'SYSTIME for Processing =', SYSTIME(1) - stp, ' sec'

end
