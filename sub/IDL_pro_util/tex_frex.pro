
PRO tex_frex, para, tex, frex, verbose=verbose

; tex = [texp,tmin] = time   expansion of para for patches > tmin by texp (on each side)
; frex = [fexp,fmin] = frequency expansion of para for patches > fmin by fexp (on each side)

nf=n_elements(para[*,0]) & nt=n_elements(para[0,*])

  if n_elements(tex) eq 0   then tex=0
  if n_elements(tex) eq 1   then tex=[tex(0),0]
  tex=tex > 0
  if n_elements(frex) eq 0   then frex=0
  if n_elements(frex) eq 1   then frex=[frex(0),0]
  frex=frex > 0

  if tex(0) gt 0 then begin
    tmin=bytarr(tex(1)+1)+1b
    texp=bytarr(2*tex(0)+1)+1b
    for i=0,nf-1 do para(i,*)=para(i,*)*(1b-dilate(dilate(erode(1b-reform(para(i,*)),tmin),tmin),texp))
  endif 

  if frex(0) gt 0 then begin
    fmin=bytarr(frex(1)+1)+1b
    fexp=bytarr(2*frex(0)+1)+1b
    for i=0,nt-1 do para(*,i)=para(*,i)*(1b-dilate(dilate(erode(1b-reform(para(*,i)),fmin),fmin),fexp))
  endif

  if tex(0) gt 0 or frex(0) gt 0 then begin
    p0=n_elements(where(para eq 0))*100.d0/n_elements(para)
    if keyword_set(VERBOSE) then print,'para (frex/tex):       ',p0,' % -> masked out'
  endif

return
END