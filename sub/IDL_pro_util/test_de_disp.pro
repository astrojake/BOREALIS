pro TEST_DE_DISP, max=max

; /max = use max frequency of band for dedispersion (default = min freq)
; attention: use the same ramp of freq of 1st step of dedispersion for the 2nd step

x=bytarr(6000,5000)
x(3000:3020,*)=1
window,0,xs=500,ys=500
tvscl,rebin(x,500,500)
rien='' & read,rien

f=findgen(5000)*0.004+12
plot,f
rien='' & read,rien

y=de_disp(x*1.,f,-2.,0.02,32)
tvscl,rebin(x+y,500,500)
rien='' & read,rien

y=de_disp(x*1.,f,-2.,0.02,1e8)
tvscl,rebin(x+y,500,500)
rien='' & read,rien

z=x
if keyword_set(max) then for i=0,9 do z(*,i*500:(i+1)*500-1)=DE_DISP(y(*,i*500:(i+1)*500-1),f(i*500:(i+1)*500-1),2.,0.02,f((i+1)*500-1)) else $
for i=0,9 do z(*,i*500:(i+1)*500-1)=DE_DISP(y(*,i*500:(i+1)*500-1),f(i*500:(i+1)*500-1),2.,0.02,f(i*500))
tvscl,rebin(x+y+z,500,500)
rien='' & read,rien

u=x
fz=rebin(findgen(10)*2+12,5000,/sample)
if keyword_set(max) then fz=fz+2
plot,fz
rien='' & read,rien

u=DE_DISP(z,fz,2.,0.02,32)
tvscl,rebin(x+y+z+u,500,500)
rien='' & read,rien

u=DE_DISP(z,fz,2.,0.02,1e8)
tvscl,rebin(x+y+z+u,500,500)

end