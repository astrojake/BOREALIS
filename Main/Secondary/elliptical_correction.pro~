;NAME: elliptical_correction
;
;Inputs: 
;        x    ; Target data
;        y    ; Sky data 
;Outputs 
;        xout ; Fixed target data 
;        yout ; fixed sky data 
;
;Comments: based off sym_scatter5 from JMaG
;
;Version 2: Dec 5, 2017 (JMaG)
;
pro elliptical_correction,x,y,xout=xnew,yout=ynew

; --------------------------------------------
; calculate Eigenvalues and Eigenvectors
npts = N_Elements(x)
i11 = Total(y^2.0) / npts
i22 = Total(x^2.0) / npts 
i12 = -Total(x * y) / npts 
tensor = [ [i11, i12], [i12,i22] ]
evals = Eigenql(tensor, Eigenvectors=evecs)
evec = evecs[*,0]
b = ((ATAN(evec[1], evec[0]) * 180. / !dpi - 90.0) +180 ) mod 180
; has to be > 0 (we do not want the symmetrical solution at -135°)
print,'orientation (before correction)',b
ell_x = Sqrt(evals[0]) * 2.0
ell_y = Sqrt(evals[1]) * 2.0
print,'semimajor and semiminor axis',ell_x,ell_y,' ratio= ',ell_x/ell_y
; ell_x and ell_y are very different, meaning we are dealing with an ellipse
r=sqrt(x^2+y^2)
t=atan(y,x)
;plot,r,t,psym=4,/polar

; --------------------------------------------
; define fitted ellipse
phi=dindgen(3600)*0.1d
xfit=ell_x*cos((phi)/!radeg)*2
yfit=ell_y*sin((phi)/!radeg)*2
rfit=sqrt(xfit^2+yfit^2)
tfit=atan(yfit,xfit)+b/!radeg

;  recalculate to get angles ordered
xfit=rfit*cos(tfit)
yfit=rfit*sin(tfit)
tfit=atan(yfit,xfit)

; --------------------------------------------
; elliptical correction
for i=0,n_elements(x)-1 do begin
	scalefactor=rfit( where( (abs(0-tfit)) eq min(abs(0-tfit)) ) )
	r(i)=r(i)/(rfit( where( (abs(t(i)-tfit)) eq min(abs(t(i)-tfit)) ) ))*scalefactor
	x(i)=r(i)*cos(t(i))
	y(i)=r(i)*sin(t(i))
endfor

;cgplot,x,y,psym=4,color='red',/overplot;,xra=[-6,6],yra=[-6,6],/iso  ; corrected x and y make a symmetrical 2D gaussian
;oplot,[0,0],[-10,10],line=1
;oplot,[-10,10],[0,0],line=1

; --------------------------------------------
; calculate Eigenvalues and Eigenvectors
npts = N_Elements(x)
i11 = Total(y^2.) / npts
i22 = Total(x^2.) / npts 
i12 = -Total(x * y) / npts 
tensor = [ [i11, i12], [i12,i22] ]
evals = Eigenql(tensor, Eigenvectors=evecs)
evec = evecs[*,0]
b = ((ATAN(evec[1], evec[0]) * 180. / !dpi - 90.0) +180 ) mod 180
; has to be > 0 (we do not want the symmetrical solution at -135°)
print,'orientation (before correction)',b
ell_x = Sqrt(evals[0]) * 2.0
ell_y = Sqrt(evals[1]) * 2.0
print,'semimajor and semiminor axis',ell_x,ell_y,' ratio= ',ell_x/ell_y
; ell_x and ell_y are similar, meaning we are dealing with a circle
; the slope is unchanged, which is essential
  
  xnew = x
  ynew = y
end
