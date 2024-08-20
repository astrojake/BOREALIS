;PROGRAM: test_detection
;
;INPUTS: 
;       Q_DIFF      :Q difference values (ON-OFF) to be tested of length threshold
;       THRESHOLD   ;threhold values to be tested
;       GQ_DIFF     ;Gaussian Q difference curve of reference (usually 1sigma)
;       MINT        ;Min threshold to test (e.g 2)
;       MAXT        ;Max threshold to test (e.g. 4.5)
;       file_g  ;Location of the Folder GAUSS folder in the LOFAR pipeline (e.g. /home/Codes/LOFAR_Bf_PipelineV3/Main/Secondary/GAUSS/) ;if not set then location is code on NANCEP2

;FLAGS: 
;     SAVE_GAUSS
;     READ_GAUSS
;     VERBOSE 
;     MINT        ; Min threshold vallue used for false-postiive calc on Q4
;     MAXT        ; Max threshold value used for false-postiive calc on Q4
;     
;OUTPUTS 
;    TEST_AVG    ; Average value above Gaussian ref curve
;    PROP_FALSE  ; Probability of Flase Detection 
;
;This program is used to test the detection of the Q4 values in post-processing code
;Used for Q4e and Q4f right now
;
;Based off Zarka's code (figs_4_7_Test.pro)
;Version 1: Sept 2018
;Version 2: Oct 2018 (updated where Gaussian is external)
pro test_detectionv2,Q_DIFF,THRESHOLD,GQ_DIFF,file_g = file_g,TEST_AVG=TEST_AVG,PROP_FALSE=PROP_FALSE,$
     RunQ4=Q4,Q4a=Q4a,Q4b=Q4b,Q4c=Q4c,Q4d=Q4d,Q4e=Q4e,Q4f=Q4f,$
     RunQ1=Q1,Q1a=Q1a,Q1b=Q1b,$
     VERBOSE=VERBOSE,MINT=MINT,MAXT=MAXT

print,'----Start Test Detection-----'

;---------------------
;Gaussian Simulations 
;---------------------
If not(keyword_set(file_g)) then file_g = '/data/jake.turner/LOFAR_Bf_PipelineV3/Main/Secondary/GAUSS/'  ;location on NANCEP2
restore,filename=file_g+'gauss_ref_FalsePostiveRate.sav'
;res_Q4a,res_Q4b,res_Q4c,res_Q4d,res_Q4e,res_Q4f,s_Q4a,s_Q4b,s_Q4c,s_Q4d,s_Q4e,$
;s_Q4f,r_Q4a,r_Q4b,r_Q4c,r_Q4d,r_Q4e,r_Q4f,N,threshold
;res_Q1a, Res_Q1b, s_Q1a, S_Q1b 

If keyword_set(Q4) then begin 
  ;flags
  If not(keyword_set(MINT)) then MINT = 1.0 ; min(THRESHOLD)
  If not(keyword_set(MAXT)) then MAXT = 5.9; max(THRESHOLD)
  If threshold_g[0] ne MINT then print,'ERROR! Min threshold not the same as in Guassian run!'
  print,'Min T: ', MINT
  print, 'Max T: ', MAXT
Endif

If keyword_set(Q1a) then r = R_Q1A
If keyword_set(Q1b) then r = R_Q1B
If keyword_set(Q4a) then r = R_Q4A
If keyword_set(Q4b) then r = R_Q4B
If keyword_set(Q4c) then r = R_Q4C
If keyword_set(Q4d) then r = R_Q4D
If keyword_set(Q4e) then r = R_Q4E
If keyword_set(Q4f) then r = R_Q4F

;-------------------------------
;Calculate probability for Q4
;------------------------------
If keyword_set(Q4) then begin
 TEST = Q_DIFF/GQ_DIFF   ; (ON-OFF)/1 sigma
 w=where(threshold ge MINT and threshold le MAXT)
 TEST_AVG=rebin(TEST(w),1)
 ww=where(threshold_g ge MINT and threshold_g le MAXT)
 nn=where(reform(rebin(r(ww,*),1,N)) ge TEST_AVG(0),count)
 NEL = n_elements(nn)
endif

;-------------------------------
;Calculate probability for Q1
;------------------------------
If keyword_set(Q1) then begin
 TEST = Q_DIFF/GQ_DIFF   ; (ON-OFF)/1 sigma
 TEST_AVG=rebin(TEST,1)
 nn=where(reform(rebin(r(*,*),1,N)) ge TEST_AVG(0),count)
 NEL = n_elements(nn)
endif


If keyword_set(VERBOSE) then print, 'Average Q Value above Gaussian ref curve: ', TEST_AVG
PROP_FALSE = (NEL*1.0d)/N
If keyword_set(VERBOSE) then print, 'Number of Gaussian Tests', N
If keyword_set(VERBOSE) then print, 'Number of Elements >= Avg Q Value: Method 1', NEL
If keyword_set(VERBOSE) then print, 'Number of Elements >= Avg Q Value: Method 2', count
If keyword_set(VERBOSE) then begin 
  If NEL ne count then print, 'Warming!!! Method 1 and Method 2 do not match that means we need to run more Gaussian simulations!'
Endif
If keyword_set(VERBOSE) then print, 'Probability of False Detection (Method1/N):', PROP_FALSE
If keyword_set(VERBOSE) then print, '% Probability of Real detection: ' ,100.0d - PROP_FALSE*100.0d

end
