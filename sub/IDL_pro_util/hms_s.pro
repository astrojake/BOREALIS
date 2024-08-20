;-------------------------------------
  function HMS_S, hms, s
;-------------------------------------
; converts hms(3) = (hour,min,sec) to sec
;-------------------------------------
  h=[3600.d0,60.d0,1.d0]
  sz=size(hms)
  if sz(0) eq 1 then s=total(transpose(hms)#h) else begin
    if sz(1) eq 3 then begin
      s=dblarr(sz(2))
      for i=0,sz(2)-1 do s(i)=total(transpose(reform(hms(*,i)))#h)
    endif else if sz(2) eq 3 then begin
      s=dblarr(sz(1))
      for i=0,sz(1)-1 do s(i)=total(transpose(reform(hms(i,*)))#h)
    endif else stop,' input format not understood'
  endelse
return, s
end
