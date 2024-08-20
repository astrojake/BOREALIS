; *********************
; *                   *
; *   SMOOTHING.PRO   *
; *                   *
; *********************

pro SMOOTHING, x,w, r, MEDIAN=MEDIAN, EDGE=EDGE

; Calculates a running average over w points, with a window narrowing at
; the edges of the original array.

; (INPUT)
; x = 1-D array of original data
; w = width of the window used for running average
; (OUTPUT)
; r = result after smoothing
; (KEYWORDS)
; MEDIAN = median average (slower)
; EDGE = manual edge truncation in smooth (slow)

  nx=n_elements(x)
  if keyword_set(MEDIAN) then begin
    r=median(x,w,/even)
    for i=1,long(w/2) do r(i)=median(x(0:2*i),/even)
    for i=nx-1-long(w/2),nx-2 do r(i)=median(x(2*i-nx+1:*),/even)
  endif else begin
    if keyword_set(EDGE) then begin
      r=smooth(x,w)
      for i=1,long(w/2) do r(i)=total(x(0:2*i),/double)/(2*i+1.)
      for i=nx-1-long(w/2),nx-2 do r(i)=total(x(2*i-nx+1:*),/double)/(2*nx-2*i-1.)
    endif else r=smooth(x,w,/edge_trunc) 
  endelse
return
end
