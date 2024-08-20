function ANG_DIST, ra1,dec1, ra2,dec2, hour=hour,hms=hms,rdms=rdms,ddms=ddms

; ra1,dec1, ra2,dec2 = scalars or vectors (by default in °), ra direct, dec from equator
; /hour --> ra in hours
; /hms  --> ra in hms = ra(3,n)
; /rdms --> ra in °'" = ra(3,n)
; /ddms --> dec in °'" = dec(3,n)
; result in °

if keyword_set(hour) then ph1=ra1*15.d0*!dtor else if keyword_set(hms) then ph1=TENV(reform(ra1(0,*)),reform(ra1(1,*)),reform(ra1(2,*)))*15.d0*!dtor else if keyword_set(dms) then ph1=TENV(double(reform(ra1(0,*))),double(reform(ra1(1,*))),double(reform(ra1(2,*))))*!dtor else ph1=double(ra1)*!dtor
if keyword_set(hour) then ph2=ra2*15.d0*!dtor else if keyword_set(hms) then ph2=TENV(reform(ra2(0,*)),reform(ra2(1,*)),reform(ra2(2,*)))*15.d0*!dtor else if keyword_set(dms) then ph2=TENV(double(reform(ra2(0,*))),double(reform(ra2(1,*))),double(reform(ra2(2,*))))*!dtor else ph2=double(ra2)*!dtor
if keyword_set(dms) then l1=TENV(double(reform(dec1(0,*))),double(reform(dec1(1,*))),double(reform(dec1(2,*))))*!dtor else l1=double(dec1)*!dtor
if keyword_set(dms) then l2=TENV(double(reform(dec2(0,*))),double(reform(dec2(1,*))),double(reform(dec2(2,*))))*!dtor else l2=double(dec2)*!dtor

;x1=cos(ph1)*cos(l1)
;y1=sin(ph1)*cos(l1)
;z1=sin(l1)
;x2=cos(ph2)*cos(l2)
;y2=sin(ph2)*cos(l2)
;z2=sin(l2)
;ad=acos(x1*x2+y1*y2+z1*z2)*!radeg

ad=acos(sin(l1)*sin(l2)+cos(l1)*cos(l2)*cos((ph1-ph2)))*!radeg

return,ad
end
