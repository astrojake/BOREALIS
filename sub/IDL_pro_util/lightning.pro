; *********************
; *                   *
; *   LIGHTNING.PRO   *
; *                   *
; *********************

pro LIGHTNING, x,largeur,nsigma,pm, para, MEDIAN=MEDIAN, EDGE=EDGE

; Subtracts a running average of 1-D array and identifies spikes > nsigma
; standard deviations.

; (INPUT)
; x = 1-D array of original data
; largeur = width of the window used for running average
; nsigma = number of standard deviations above which a value is considered as a spike
; pm = if>0 removes spikes > nsigma only, if<0 removes spikes < nsigma only, if=0 removes both
; (OUTPUT)
; para = weights identifying good (1) and bad (0) pixels in the original data
; (KEYWORDS)
; MEDIAN, EDGE = passed to SMOOTHING

  if not(keyword_set(MEDIAN)) then MEDIAN=0
  if not(keyword_set(EDGE)) then EDGE=0
  SMOOTHING, x,largeur, sx, MEDIAN=MEDIAN, EDGE=EDGE
  rx=x-sx
  para=bytarr(n_elements(x))+1b
  sigma=stddev(rx,/double)
  if pm lt 0 then test=where(rx lt -nsigma*sigma)
  if pm eq 0 then test=where(abs(rx) gt nsigma*sigma)
  if pm gt 0 then test=where(rx gt nsigma*sigma)
  if test(0) ne -1 then para(test)=0b
return
end
