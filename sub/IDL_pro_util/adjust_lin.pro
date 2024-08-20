; **********************
; *                    *
; *   ADJUST_LIN.PRO   *
; *                    *
; **********************

pro ADJUST_LIN, x,para, xnet, NPIX=NPIX

; Linear interpolation of bad pixels in a 1-D array. 
; Bad pixels at the beginning or end of the array
; are replaced by constant (first or last) values.

; (INPUT)
; x = 1-D array of original data
; para = weights identifying good (1) and bad (0) pixels in the original data
; NPIX = number of good pixels on either side of each bad pixel to interpolate
; (OUTPUT)
; xnet = cleaned data, with bad pixels replaced by interpolated values

  if not(keyword_set(NPIX)) then npix=10 else npix=fix(abs(npix))
  nx=n_elements(x)
  xx=lindgen(nx)
  xnet=x
  para2=para
  test=where(para2 ne 0b)
  if test(0) eq -1 then begin			; entirely bad arrays
	xnet=x-x & return			; are set to zero
  endif
  ntest=n_elements(test)
  xnet(0:test(0))=xnet(test(0))			& para2(0:test(0))=1b
  xnet(test(ntest-1):*)=xnet(test(ntest-1))	& para2(test(ntest-1):*)=1b
  test=where(para2 eq 0b)
  if test(0) eq -1 then return
  ntest=n_elements(test)

  for i=0,ntest-1 do begin
    j=test(i)
    if para2(j-1) eq 0b and NPIX gt 1 then xnet(j)=result(j-a1)  else $
    if para2(j-1) eq 0b and NPIX eq 1 then xnet(j)=xnet(j-1)     else begin
      avant=max(where(para2(0:j) eq 1b))
      apres=min(where(para2(j:*) eq 1b))+j
      if NPIX gt 1 then begin
	a1=max([0,avant-npix]) & a2=min([apres+npix,nx-1])
	aa=where(para(a1:a2) eq 1b)
	rien=linfit(xx(a1+aa),x(a1+aa))
	result=rien(0)+rien(1)*xx(a1:a2)
;	rien=poly_fit(xx(a1:a2),x(a1:a2),1,meas=1./(para2(a1:a2)*1.e6+1.e-6),yfit=result)
	xnet(j)=result(j-a1)
      endif else xnet(j)=(x(avant)+x(apres))/2.
    endelse
  endfor
return
end
