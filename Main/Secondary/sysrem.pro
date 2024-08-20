;Name: sysrem
;PURPOSE: Remove systematic effects in large sets of light curves or wavelengths
;Citation: Tamuz, Mazeh, & Zucker, 2005, MNRAS, 356, 1466, "Correcting systematic effects in a large set of photometric light curves"
;         
;Input:
;      t[M]       = t[frames]                         ;time 
;      flux[N,M]  = Flux[stars or wavelengths,frames] ;flux
;      error[N,M]                                     ;error on flux
;
;Outputs:
;      Cflux[N,M]  ; The  corrected flux subtracting the systemaic error ; flux - sys_err
;
;Input Flags: 
;      iter        ;The number of iterations (default is 10)
;      MAG         ;Convert flux to mags 
;Output Flags:
;       A      
;       MEMSAV     ;Sav Memory (do not output A or C)
;https://github.com/stephtdouglas/PySysRem/blob/master/sysrem.py
;The Sys-Rem Detrending Algorithm:  Implementation and Testing; https://arxiv.org/pdf/astro-ph/0612418.pdf; Stopping Criteria
;https://arxiv.org/abs/1003.0427; on actual stars for flux
;
pro sysrem,flux,error,FixFLUX=CFLUX,ITER=ITER,A=A,C=C,sys_err=sys_err,MAG=MAG,MEMSAV=MEMSAV

If not(keyword_set(ITER)) then ITER = 10  ;numer of iterations

;Convert flux to magnitudes
If keyword_set(MAG) then begin 
  flux = -2.5*alog10(flux)     ;convert to mags 
  error = (2.5*1/flux)*error  ; sigma^2 = df/dx^2*error^2
Endif

;Elements
N  = n_elements(flux[*,0]) ;# of stars/light curves or wavelengths
M  = n_elements(flux[0,*]) ;# of frames

;Initialize c and a
c = fltarr(N)           
a = fltarr(M) + 1.0d    ;the starting value 

for k=0,ITER-1 do begin   ;iterations to be done 
  ;print, k
  ;---------
  ;Find C(i)
  ;---------
  num = 0.0d  &  dem = 0.0d   ;initalize
  for i=0,N-1 do begin
    err_sq  = total(error[i,*]^2.0d,/nan)
    num = (total(flux[i,*]*a[*],/nan))/err_sq 
    dem       = total(a[*]^2.0d,/nan)/err_sq  
    c[i] = num/dem                                    ;Equation 2
  endfor

  ;---------
  ;Find a(j)
  ;---------
  num_aj = 0.0d  &  dem_aj = 0.0d   ;initalize
  for j=0,M-1 do begin      ;check!!! 
    err_sq      = total(error[*,j]^2.0d,/nan)
    num_aj     = total(flux[*,j]*c[*],/nan)/err_sq
    dem_aj     = total(c[*]^2.0d,/nan)/err_sq
    a[j] = num_aj/dem_aj
  endfor  
endfor  ;end iterations (k) 

sys_err = fltarr(N,M)   ;c_{i}*a_{j}
for j=0,M-1 do begin
 for i=0,N-1 do begin
  sys_err[i,j]=c[i]*a[j]
 endfor
endfor 
 
 
 ;Remove the Systematic error
 CFLUX = flux - sys_err        ;r_{ij}(1) = r_{ij} - c_{i}*a_{j}
  
 If keyword_set(MAG) then begin
   CFLUX = 10^(-2.5*CFLUX)     ;convert to mags
 Endif
  
end  