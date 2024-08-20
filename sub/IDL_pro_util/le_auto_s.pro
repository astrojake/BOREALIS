; *********************
; *                   *
; *   LE_AUTO_S.PRO   *
; *                   *
; *********************

pro LE_AUTO_S, x,largeur,nsigma,pm, xnet,para, logfile, $
	MEDIAN=MEDIAN, EDGE=EDGE, NPIX=NPIX, MAXITERATIONS=MAXITERATIONS, SAVSAMPLE=SAVSAMPLE, verbose=verbose

; Automated & iterative identification of spikes in 1-D array. Spikes are replaced by interpolated values. 
; Iteration is made until no peak >nsigma*sigmas remains.
; This program call subroutines: LIGHTNING, ADJUST_LIN, SMOOTHING.
; Difference with LE_AUTO consists in no correction for edge effects here.

; (INPUT)
; x = 1-D array of original data
; largeur = width of the window used for running average of original data as well as interpolation
; nsigma = number of standard deviations above which a value is considered as a spike
; pm = if>0 removes spikes > nsigma only, if<0 removes spikes < nsigma only, if=0 removes both (parameter passed to LIGHTNING)
; (OUTPUT)
; xnet = cleaned data, with spikes replaced by interpolated values
; para = weights identifying good (1) and bad (0) pixels in the original data
; logfile = optional log text file
; (KEYWORDS)
; MEDIAN, EDGE = passed to LIGHTNING
; NPIX = passed to ADJUST_LIN
; MAXITERATIONS = max number of iterations
; SAVSAMPLE = N --> save in sampleN.sav file examples for which iterations > 10
; verbose = comments printing

DIR_SEPARATOR,sep
width=largeur < n_elements(x)

  if not(keyword_set(MEDIAN)) then MEDIAN=0
  if not(keyword_set(EDGE)) then EDGE=0
  if not(keyword_set(NPIX)) then npix=10 else npix=fix(abs(npix))
  if not(keyword_set(MAXITERATIONS)) then maxiterations=20 else maxiterations=fix(maxiterations)
  nx=n_elements(x)
  x1=x & testpred=0
  para=bytarr(nx)+1b
  for i=1,maxiterations do begin
    LIGHTNING, x1,width,nsigma,pm, para1, MEDIAN=MEDIAN, EDGE=EDGE
    ADJUST_LIN, x1,para1, xnet, NPIX=NPIX
    para=para*para1
    test=n_elements(where(para eq 0b))
    if abs(test-testpred) le 1 then goto,fin
    x1=xnet
    testpred=test
  endfor
fin:
  if keyword_set(VERBOSE) then print,'Number of iterations = ',i,', ',long(total(1-para,/double)),' / ',long(n_elements(para))
  if n_elements(logfile) ne 0 and ((i gt 10 and keyword_set(SAVSAMPLE)) or keyword_set(VERBOSE)) then begin
	openw,v, logfile, /get_lun, /append
	if keyword_set(SAVSAMPLE) ne 0 then $
	  printf,v,'Number of iterations = ',i,', ',long(total(1-para,/double)),' / ',long(n_elements(para)),'   ',strtrim(fix(SAVSAMPLE),2)+'_S.sav' $
	  else printf,v,'Number of iterations = ',i,', ',long(total(1-para,/double)),' / ',long(n_elements(para))
	close, v & free_lun, v
  endif
  if keyword_set(SAVSAMPLE) ne 0 and (i gt 10 or keyword_set(VERBOSE)) then begin $
	save,x,largeur,width,nsigma,pm,xnet,para,i,MEDIAN,EDGE,NPIX,MAXITERATIONS,file='.'+sep+'samples'+sep+strtrim(fix(SAVSAMPLE),2)+'_S.sav'
	SAVSAMPLE=SAVSAMPLE+1
  endif

return
end
