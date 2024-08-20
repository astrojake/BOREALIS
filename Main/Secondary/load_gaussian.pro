;INPUTS: 
;        threshold_ORG   ; The threshold that the PP is using (used to check that the same thresholds are used)
;        file          ; Where the Gaussian Q values are located
;        N               ; Number of elements in Q2 or y (high pass filtered time series)
;OUTPUTS        
;        Q_Gv1               ; Gaussian Q values (a-f) to calculated
;
pro load_gaussian,threshold_ORG,N,Q=Q_Gv1,file=file

  n_thres = n_elements(threshold_ORG) ;number of thresholds to run 
  ;--------------------------
  ;Restore Gaussian Q values
  ;--------------------------
  restore,filename=file+'Gaussian_QA.sav'
  ;yfit_A,filename='Gaussian_QA.sav'
  restore,filename=file+'Gaussian_QA_D.sav'
  ;yfit_A_D,filename='Gaussian_QA_D.sav'
  restore,filename=file+'Gaussian_QB.sav'
  ;yfit_B,filename='Gaussian_QB.sav'
  restore,filename=file+'Gaussian_QB_D.sav'
  ;yfit_B_D,filename='Gaussian_QB_D.sav'
  restore,filename=file+'Gaussian_QC.sav'
  ;yfit_C,filename='Gaussian_QC.sav'
  restore,filename=file+'Gaussian_QC_D.sav'
  ;yfit_C_D,filename='Gaussian_QC_D.sav'
  restore,filename=file+'Gaussian_QD.sav'
  ;yfit_D,filename='Gaussian_QD.sav'
  restore,filename=file+'Gaussian_QD_D.sav'
  ;yfit_D_D,filename='Gaussian_QD_D.sav'
  restore,filename=file+'Gaussian_QE.sav'
  ;yfit_E,filename='Gaussian_QE.sav'
  restore,filename=file+'Gaussian_QE_D.sav'
  ;yfit_E_D,filename='Gaussian_QE_D.sav'
  restore,filename=file+'Gaussian_QF.sav'
  ;yfit_F,filename='Gaussian_QF.sav'
  restore,filename=file+'Gaussian_QF_D.sav'
  ;yfit_F_D,filename='Gaussian_QF_D.sav'


  for i=0,n_thres-1 do begin 
    If threshold_ORG[i] ne threshold[i] then begin 
      print,'We have a problem!!! Thresholds are not the same!'
      print,'Program will stop!! Did you change the default thresholds? Dont do that'
      print,'STOP!!!'
      stop
    Endif  
  endfor
  
  Q_Gv1 = dblarr(2,n_thres,6)    ;master Q Gaussian 
    ;e.g. q4_gauss[0,*,0] ;sum of 4a for Guass
    ;     q4_gauss[1,*,0] ;Difference 4a for Gauss
 
 ;--------- 
 ;Qa
 ;---------
 Qa         = yfit_A[0,*] + yfit_A[1,*]*N   ;a + bx, where x (N) is the number of elements
 Q_Gv1(0,*,0) = Qa
 Qa_d       = yfit_A_D[0,*] + yfit_A_D[1,*]*N + yfit_A_D[2,*]*N^(2.0d)   ;a + bx + cx^2, where x (N) is the number of elements
 Q_Gv1(1,*,0) = Qa_d   

 ;---------
 ;Qb
 ;---------
 Qb         = yfit_B[0,*] + yfit_B[1,*]*N   ;a + bx, where x (N) is the number of elements
 Q_Gv1(0,*,1) = Qb
 Qb_d       = yfit_B_D[0,*] + yfit_B_D[1,*]*N + yfit_B_D[2,*]*N^(2.0d)   ;a + bx + cx^2, where x (N) is the number of elements
 Q_Gv1(1,*,1) = QB_d

 ;---------
 ;Qe
 ;---------
 Qe         = yfit_E[0,*] + yfit_E[1,*]*N   ;a + bx, where x (N) is the number of elements
 Q_Gv1(0,*,4) = Qe
 Qe_d       = yfit_E_D[0,*] + yfit_E_D[1,*]*N + yfit_E_D[2,*]*N^(2.0d)   ;a + bx + cx^2, where x (N) is the number of elements
 Q_Gv1(1,*,4) = Qe_d

 ;---------
 ;Qf
 ;---------
 Qf         = yfit_F[0,*] + yfit_F[1,*]*N   ;a + bx, where x (N) is the number of elements
 Q_Gv1(0,*,5) = Qf
 Qf_d       = yfit_F_D[0,*] + yfit_F_D[1,*]*N + yfit_F_D[2,*]*N^(2.0d)   ;a + bx + cx^2, where x (N) is the number of elements
 Q_Gv1(1,*,5) = Qf_d

 ;---------
 ;Qc
 ;---------
 Qc        = yfit_C[0,*] + yfit_C[1,*]*N  + yfit_C[2,*]*N^(2.0d)  ;a + bx + cx^2, where x (N) is the number of elements
 Q_Gv1(0,*,2) = Qc
 Qc_d    = yfit_C_D[0,*] + yfit_C_D[1,*]*N + yfit_C_D[2,*]*N^(2.0d)   ;a + bx + cx^2, where x (N) is the number of elements
 Q_Gv1(1,*,2) = QC_d
 
 ;---------
 ;Qd
 ;---------
 Qd        = yfit_D[0,*] + yfit_D[1,*]*N  + yfit_D[2,*]*N^(2.0d)  ;a + bx + cx^2, where x (N) is the number of elements
 Q_Gv1(0,*,3) = Qd
 QD_d    = yfit_D_D[0,*] + yfit_D_D[1,*]*N + yfit_D_D[2,*]*N^(2.0d)   ;a + bx + cx^2, where x (N) is the number of elements
 Q_Gv1(1,*,3) = QD_d
 
end