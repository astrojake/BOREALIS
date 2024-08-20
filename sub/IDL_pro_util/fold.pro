; -----------------------------------------------------------------------
pro FOLD, xfy,fy,range,xfyf,fyf,h
  ; -----------------------------------------------------------------------
  ; folding spectrum harmonics

  ;h=[0.5,1,2,3,4,5,6,7,8,9,10]
  if n_elements(h) eq 0 then h=[1,2,3,4,5,6]
  nh=n_elements(h)

  w=where(xfy ge range(0) and xfy le range(1))
  nw=n_elements(w)

  xfyf=xfy(w)
  fyf=fltarr(nw)

  for i=0,nw-1 do begin
    for j=0,nh-1 do begin
      ww=where(abs(xfy-xfy(w(i))/h(j)) eq min(abs(xfy-xfy(w(i))/h(j))))
      fyf(i)=fyf(i)+fy(ww(0))
    endfor
  endfor

  return
end
