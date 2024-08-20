;This file creates the master Gaussian used for the LOFAR PP Analysis

;Number of Guassian steps: 1e6
;Orginal size of setup (y): 1e6
;Threshold range: The DEFAULT RANGE FOR LOFAR Pipeline V2: 1 to 5.9 sigma in increments of 0.1 sigma

;gauss_steps      = 1d4  ;gaussian steps
;N_Gauss          = 5d4  ;Size of setup
pro master_gauss,gauss_steps,N_Gauss,slope=slope

If not(keyword_set(slope)) then slope = 2.0d

;Threshold setup
N      = 50
threshold = dindgen(50,start=1,increment=0.1)
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
  print,i
  gauss_lofar_postprocessing,gauss_steps,N_Gauss,slope=slope,threshold(i),$
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
endfor  

save,Gauss_Steps,N_Gauss,threshold,A,A_d,B,B_d,C,C_d,D,D_d,E,E_d,F,F_d,filename='LOFAR_PP_Gaussian.sav'

end