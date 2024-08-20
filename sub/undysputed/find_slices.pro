pro FIND_SLICES, namefits, wbegin, wend, nslices

READ_NU_FITS, namefits, command,param,variab, nt,dt,nf,df,ns,jd0,h0,fref, data,time,freq,beam,ndata,corrt,corrf,nstokes=1,/quiet

w=[0,where(corrt(0:-2,3) ne corrt(1:-1,3) and corrt(0:-2,4) ne corrt(1:-1,4)),nt]
while max(w(1:*)-w) gt ceil(360./dt)+1 do begin
  wm=where(w(1:*)-w eq max(w(1:*)-w))+1 & wm=wm(0)
  w=[w(0:wm-1),w(wm)-360./dt,w(wm:*)]
endwhile
w=round(w)
nw=n_elements(w)

wind=[0]
for i=0,nw-1 do begin
  if w(i)-wind(-1) gt 10 then begin
    if w(i)-wind(-1) lt 1500 then wind=[wind,w(i)] else wind=[wind,w(i)-long((w(i)-w(i-1))/2),w(i)]
  endif
endfor
nslices=n_elements(wind)-1
wbegin=wind(0:nslices-1) & wend=wind(1:nslices)-1
end
