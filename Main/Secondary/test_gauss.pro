
gauss_steps = 5d5

N = 2  
y           = DINDGEN(N,start = 1000,increment=10000)

A = dblarr(N)
A_d = dblarr(N)
B = dblarr(N)
B_d = dblarr(N)
C = dblarr(N)
C_d = dblarr(N)
D = dblarr(N)
D_d = dblarr(N)
E = dblarr(N)
E_d = dblarr(N)
F = dblarr(N)
F_d = dblarr(N)

for i=0,N-1 do begin
  print, 'start: ', i
  gauss_lofar_postprocessing,gauss_steps,y(i),slope=2,threshold,$
    Q3a=Ao,dQ3a=A_do,Q3b=Bo,dQ3b=B_do,Q3c=Co,dQ3c=C_do,Q3d=dout,dQ3d=d_do,Q3e=eo,dQ3e=E_Do,Q3f=Fo,dQ3f=F_Do  
 
  A(i) = Ao     
  A_d(i) = A_do
  B(i) = Bo
  B_d(i) = B_do
  C(i) = Co
  C_d(i) = C_do
  D(i) = Dout
  D_d(i) = D_do
  E(i) = Eo
  E_d(i) = E_do
  F(i) = Fo
  F_d(i) = F_do
  print,'End: ', i
endfor

 Print,'A Ratio: ', A(N-1)/A(0)
 Print,'B Ratio: ', B(N-1)/B(0)
 print, 'C ratio (corrected): ', C(N-1)/C(0), sqrt(11)
 print, 'D ratio (corrected): ',D(N-1)/D(0), sqrt(y(N-1)/y(0))
 Print,'E Ratio: ', E(N-1)/E(0)
 Print,'F Ratio: ', F(N-1)/F(0)

end
