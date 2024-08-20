;-------------------------------------
  function S_HMS, s, hms
;-------------------------------------
; converts sec to hms(3) = (hour,min,sec)
;-------------------------------------
  n=n_elements(s)
  hms=dblarr(3,n)
  for i=0,n-1 do begin
    hh=double(long(s(i)/3600.d0))
    mm=double(long((s(i)-hh*3600.d0)/60.d0))
    ss=s(i)-hh*3600.-mm*60.
    hms(*,i)=[hh,mm,ss]
  endfor
  hms=reform(hms)
return, hms
end
