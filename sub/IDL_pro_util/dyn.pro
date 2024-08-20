;--------------------------------------------------------------------
  pro DYN, tab,fracmin,fracmax, tabmin,tabmax, QUIET=QUIET
;--------------------------------------------------------------------
; Adjustment of dynamic range
; use LOG scale if large range of values in tab

; (INPUT)
; tab = array of values
; fracmin,fracmax = quantiles for which ...
; (OUTPUT)
; ... tabmin,tabmax values are required

  xmin=min(tab) & xmax=max(tab)
  if xmin eq xmax then begin
    if not(keyword_set(QUIET)) then print,'Null dynamic range'
    tabmin=xmin & tabmax=xmax & return
  endif
  dh=(xmax-xmin)/1000.
  h=histogram(float(tab),min=xmin,max=xmax,binsize=dh,loc=xh)
  nh=n_elements(h)
  for i=1L,nh-1L do h(i)=h(i)+h(i-1)
  th=h(nh-1)
  test=abs(h-fracmin*th)
  tabmin=xh(where(test eq min(test)))
  tabmin=tabmin(0)
  test=abs(h-fracmax*th)
  tabmax=xh(where(test eq min(test)))
  tabmax=tabmax(0)
return
end
