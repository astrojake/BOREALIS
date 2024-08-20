;PROGRAM: test_detection
;
;INPUTS: 
;       Q_DIFF      :Q difference values (ON-OFF) to be tested of length threshold
;       THRESHOLD   ;threhold values to be tested
;       GQ_DIFF     ;Gaussian Q difference curve of reference (usually 1sigma)
;       MINT        ;Min threshold to test (e.g 2)
;       MAXT        ;Max threshold to test (e.g. 4.5)
;FLAGS: 
;     SAVE_GAUSS
;     READ_GAUSS
;     VERBOSE 
;
;OUTPUTS 
;    TEST_AVG    ; Average value above Gaussian ref curve
;    PROP_FALSE  ; Probability of Flase Detection 
;
;This program is used to test the detection of the Q4 values in post-processing code
;Used for Q4e and Q4f right now
;
;Based off Zarka's code (figs_4_7_Test.pro)
pro test_detection,Q_DIFF,THRESHOLD,GQ_DIFF,MINT=MINT,MAXT=MAXT,save_gauss=save_gauss,read_gauss=read_gauss,VERBOSE=VERBOSE,TEST_AVG=TEST_AVG,PROP_FALSE=PROP_FALSE

print,'----Start Test Detection-----'
;flags
If not(keyword_set(MINT)) then MINT = 1.5 ; min(THRESHOLD)
If not(keyword_set(MAXT)) then MAXT = 4.5; max(THRESHOLD)
print,'Min T ', MINT
print, 'Max T: ', MAXT 

;---------------------
;Gaussian Simulations 
;---------------------
If not(keyword_set(read_gauss)) then begin 
  n_th = n_elements(threshold)
  N = 10000.0d
  res=fltarr(n_th,N)
  for j=0,N-1 do begin
    x=randomn(seed,N)
    y=randomn(seed,N)
    z=randomn(seed,N)
    rx=fltarr(n_th) & rz=rx
    for i=0,n_th-1 do rx(i)=total(x(where(x ge threshold(i) and x ge 2.*y)))   ;Q4f
    for i=0,n_th-1 do rz(i)=total(z(where(z ge threshold(i) and z ge 2.*y)))   ;Q4f
    res(*,j)=(rx-rz)  ; /10000.*100.  ; %
  endfor
    s=fltarr(50) & for i=0,49 do s(i)=stddev(res(i,*))
    r=res/rebin(reform(s > (1.e-5),50,1),50,N)
    If keyword_set(save_gauss) then save,res,s,r,N,threshold,file='gauss_ref.sav'
endif 
If keyword_set(read_gauss) then restore,filename='gauss_ref.sav'

;----------------------
;Calculate probability
;----------------------
TEST = Q_DIFF/GQ_DIFF   ; (ON-OFF)/1 sigma

w=where(threshold ge MINT and threshold le MAXT)
TEST_AVG=rebin(TEST(w),1)
nn=where(reform(rebin(r(w,*),1,N)) ge TEST_AVG(0))
NEL = n_elements(nn)

print, 'Average Q Value: ', TEST_AVG
PROP_FALSE = (NEL*1.0d)/N
print, 'Probability of False Detection:', PROP_FALSE
print, '% Probability of Real detection: ' ,100.0d - PROP_FALSE*100.0d

end