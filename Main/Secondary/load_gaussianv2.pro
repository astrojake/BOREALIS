;Purpose: The purpose of this program is to find the guassian Q4 values for any given length N. 
;         The program uses the folder where Gaussian_Master.sav located in /LOFAR_Bf_PipelineV3/Main/Secondary/GAUSS/ for the master guassian curves
;          
;INPUTS: 
;        threshold_ORG ; The threshold that the PP is using (used to check that the same thresholds are used; should be 1 to 5.9)
;        file          ; The folder where Gaussian_Master.sav is located (default would be /LOFAR_Bf_PipelineV3/Main/Secondary/GAUSS/)
;        N             ; Number of elements in Q2 or y (high pass filtered time series)
;        
;SAV FILE INPUT: 
;        Gaussian_Master.sav 
;           ;3 arrays: gauss_fit --> Qn, gauss_fit_diff --> Qn diff plots , threshold --> thresholds used to create fits (1 to 5.9)
;           ;format (6,3,n_threshold); Qn (a-f),fit constants (abc where, a + bx + cx^2),threshold  
;           ;Notes: Gauss_fit and gauss_fit_diff are fits (either quad or linear) to many large N Gaussian simulations (has been verfied)
;           
;OUTPUTS        
;        Q_G               ; Gaussian Q values (a-f) to be calculated 
;                          ;  format: (2,n_threshold,6): 2 is for Q and Qdiff, 6 is for Qa-f  
;                          ;  e.g. q4_gauss[0,*,n] ;Q4n for Guassian
;                          ;       q4_gauss[1,*,n] ;Difference 4n for Gaussian
;Version : 
;         Version 1: Sept 18, 2018
;         Version 2: Nov 1, 2018 (consolidated the program and updated header) 
pro load_gaussianv2,threshold_ORG,N,Q=Q_G,file=file

  n_thres = n_elements(threshold_ORG) ;number of thresholds to run 
  ;--------------------------
  ;Restore Gaussian Q values
  ;--------------------------
  restore,filename=file+'Gaussian_Master.sav'
   ;threshold,gauss_fit,gauss_fit_diff

  for i=0,n_thres-1 do begin 
    If threshold_ORG[i] ne threshold[i] then begin 
      print,'We have a problem!!! Thresholds are not the same!'
      print,'Program will stop!! Did you change the default thresholds? Dont do that'
      print,'STOP!!!'
      stop
    Endif  
  endfor
  
  Q_G = dblarr(2,n_thres,6)    ;master Q Gaussian 
    ;e.g. q4_gauss[0,*,0] ;sum of 4a for Guass
    ;     q4_gauss[1,*,0] ;Difference 4a for Gauss
 
 ;---------------------------
 ;Find Q value for input N
 ;Gauss_fit[Qa-f,fit constant (abc),threshold]
 ;---------------------------
 for i=0,5 do begin          ;for loop over Qn (a-f)
    Q_G(0,*,i)  = gauss_fit[i,0,*]      + gauss_fit[i,1,*]*N          + gauss_fit[i,2,*]*N^2.0d             ;Q value
    Q_G(1,*,i)  = gauss_fit_diff[i,0,*] + gauss_fit_diff[i,1,*]*N     + gauss_fit_diff[i,2,*]*N^2.0d        ;Q diff  
 endfor 

end