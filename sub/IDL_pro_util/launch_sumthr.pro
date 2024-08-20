;------------------------------------------------------------------------------
 pro launch_sumthr, x, p, Mf, Mt, BACK=BACK, THR0=THR0
;------------------------------------------------------------------------------
; SumThreshold of VarThreshold on data x(nt,nf) & previous flags p(nt,nf)
; flags rows then columns separately for the list of window sizes specified in
; Mf, Mt (preferably powers of 2)
; /BACK = pre-normalize data by background
; THR0 = constant defining the threshold series [default = 10.]

if not(keyword_set(thr0)) then thr0=10.
sx=size(x)
nt=sx[1] & nf=sx[2]
pf=p & pt=p

; processing frequencies
nm=n_elements(Mf)
thr=thr0/(1.5^(alog(Mf)/alog(2))) 
for i=0,nf-1 do begin
  y=reform(pf[*,i])
  if keyword_set(BACK) then begin
    BACKGROUND, reform(x[*,i]), mx, sx, ny
    xx=(reform(x[*,i])-mx)/sx
  endif else xx=reform(x[*,i])
  for j=0,nm-1 do sumthr, xx, Mf[j], thr[j], y
  pf[*,i]=y
endfor

; processing spectra
nm=n_elements(Mt)
thr=thr0/(1.5^(alog(Mt)/alog(2)))
for i=0,nt-1 do begin
  y=reform(pt[i,*])
  if keyword_set(BACK) then begin
    BACKGROUND, reform(x[i,*]), mx, sx, ny
    xx=(reform(x[i,*])-mx)/sx
  endif else xx=reform(x[i,*])
  for j=0,nm-1 do sumthr, xx, Mt[j], thr[j], y
  pt[i,*]=y
endfor

p=pf*pt

return
end
