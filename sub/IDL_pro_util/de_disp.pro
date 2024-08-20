function DE_DISP,x,xf,dm,dt,xfref

; x = data (time,frequency)
; xf = frequency ramp [MHz]
; dm = dispersion measure [pc.cm-3]
; dt = time resolution [sec]
; xfref = reference frequency for dispersion

k = 1.d0/2.410331d-4    ; = 4148.8 (Zakharenko et al., 2013)
if n_elements(xfref) eq 0 then xfref=min(xf)
delta=k*dm/(xf^2)    
deltaref=k*dm/(xfref^2)    
delta=delta-deltaref
dpix=long(delta/dt+0.5)
nf=n_elements(x(0,*))
y=x
for i=0,nf-1 do y(*,i)=shift(reform(y(*,i)),-dpix(i))
return,y
end
