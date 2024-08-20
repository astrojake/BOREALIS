;---------------------------------------------------------------------------------------
  pro make_background,tab,save_set,f,s,n,SUPMIN=SUPMIN,POSITIVE=POSITIVE,NSIG=NSIG
;---------------------------------------------------------------------------------------
; Calculation of Background and Sigma arrays of a 2-D distribution of 
; intensities, frequency by frequency. A save set can be put in 'save_set'

; tab =  (INPUT)  2-D array of intensities, function of (time,freq)
; save_set =  (INPUT)  name of save_set containing f & s
; f,s,n =  (OUTPUT)  arrays of background, fluctuations (1 sigma level), 
;                    and number of pixels used per frequency
; SUPMIN = only values >min(tab) are used
; POSITIVE = retains only tab values > 0
; NSIG = number of sigmas

  nf=n_elements(tab(0,*))
  f=fltarr(nf) & s=f & n=f
  if keyword_set(SUPMIN) then mintab=min(tab)
  for i=0,nf-1 do begin
    if keyword_set(SUPMIN) then begin
      test=where(tab(*,i) gt mintab)
      if test(0) ne -1 then begin
	background,tab(test,i),ff,ss,nn,POSITIVE=POSITIVE,NSIG=NSIG
	f(i)=ff & s(i)=ss & n(i)=nn
      endif
    endif else begin
      background,tab(*,i),ff,ss,nn,POSITIVE=POSITIVE,NSIG=NSIG
      f(i)=ff & s(i)=ss & n(i)=nn
    endelse
  endfor
  if save_set ne '' then save,f,s,n,filename=save_set
return
end
