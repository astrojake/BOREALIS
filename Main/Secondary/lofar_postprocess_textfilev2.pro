;Program to write/read output text file from postprocessing
;
;filename has to be a text file and created with the postprocessing code
;to run: .r lofar_postprocess_textfile.pro
;        lofar_postprocess_textfile,filename,/READ
;
;See OVERVIEW file for a description of the outputs 
;
;Versions
;    V1: Sept 2018 
;    V2; To include probability of false Detection, Oct 2018
pro lofar_postprocess_textfilev2,filename,Date=Date,S=S,fmin=fmin,fmax=fmax,BT=TBeam,BSky=SBeam,$
    	Q1aV=Q1aV,Q1bV=Q1bV,Q4aV=Q4aV,Q4bV=Q4bV,Q4cV=Q4CV,Q4dV=Q4dV,Q4eV=Q4eV,Q4fV=Q4fV,Detect=Detect,$
        WRITE=WRITE,READ=READ

If not(keyword_set(Date)) then date ='L570725'
If not(keyword_set(S)) then S = 0.0d
If not(Keyword_Set(fmin)) then fmin=0.0d
If not(Keyword_Set(fmax)) then fmax=0.0d
;If not(keyword_set(BT)) then tbeam= 0.0d
;If not(keyword_set(BSky)) then SBeam = 0.d 

If keyword_set(WRITE) then begin
 ;---------------
 ;Open
 ;--------------
 openw,lun,filename+'.txt',/get_lun ;flag
 
 ;------------------------------
 printf,lun,'***Output File*****'
 printf,lun,'*******************'
 printf,lun,filename,format='(A)'
 printf,lun,Date,format='(A-,T45,"Date")'
 printf,lun,S,format='(A-,T45,"Polarization")'
 printf,lun,fmin,format='(F-6.2,T45,"Min Frequency (MHz):")'
 printf,lun,fmax,format='(F-6.2,T45,"Max Frequency (MHz):")'
 printf,lun,TBeam,format='(I-1,T45,"Target Beam:")'
 printf,lun,SBeam,format='(I-1,T45,"Sky Beam:")'
 printf,lun,'********************'
 printf,lun,''
 printf,lun,'*********************************'
 printf,lun,'       Q1a: Time Series          '
 printf,lun,'*********************************'
 printf,lun,Q1aV[0],format='(F-,T45,"Integrated ON-OFF level (SEFD)")'
 printf,lun,Q1aV[1],format='(F-6.2,T45,"% where ON > OFF")'
 printf,lun,Q1aV[2],format='(F-,T45,"STDEV of ON-OFF (SEFD)")'
 printf,lun,Q1aV[3],format='(D11.8,T45,"Probability of False Detection")'
 printf,lun,'*********************************'
 printf,lun,'    Q1b: Integrated Spectrum     '
 printf,lun,'*********************************'
 printf,lun,Q1bV[0],format='(F-,T45,"Integrated ON-OFF level (SEFD)")'
 printf,lun,Q1bV[1],format='(F-6.2,T45,"% where ON > OFF")'
 printf,lun,Q1bV[2],format='(F-,T45,"STDEV of ON-OFF (SEFD)")'
 printf,lun,Q1bV[3],format='(D11.8,T45,"Probability of False Detection")'
 printf,lun,'*********************************'
 printf,lun,'    Q4a: Number of Peaks         '
 printf,lun,'*********************************'
 printf,lun,Q4aV[0],format='(F-8.2,T45,"Integral of ON-OFF curve")'
 printf,lun,Q4aV[1],format='(F-6.3,T45,"Ratio of ON-OFF Integral to 1 sigma Gaussian")'
 printf,lun,Q4aV[4],format='(F-6.2,T45,"% of ON-OFF values > 0:")'
 printf,lun,Q4aV[5],Q4aV[6],Q4aV[7],format='(I-6,I-6,I-6,T45,"Number of ON-OFF values above 1,2,3 sigma")'
 printf,lun,Q4aV[8],Q4aV[9],format='(F-4.2,2x,F-4.2,T45,"Min and Max Threshold value where ON-OFF > 0")'
 printf,lun,Q4aV[10],format='(-F6.2,T45,"% of ON-OFF values between min and max threshold & >0")'
 printf,lun,Q4aV[11],format='(F-6.2,T45,"% of ON-OFF values lt 0 and > max threshold")'
 printf,lun,Q4aV[12],format='(D11.8,T45,"Probability of False Detection")'
 printf,lun,'*********************************'
 printf,lun,'    Q4b: Power of Peaks          '
 printf,lun,'*********************************'
 printf,lun,Q4bV[0],format='(F-8.2,T45,"Integral of ON-OFF curve")'
 printf,lun,Q4bV[1],format='(F-6.3,T45,"Ratio of ON-OFF Integral to 1 sigma Gaussian")'
 printf,lun,Q4bV[4],format='(F-6.2,T45,"% of ON-OFF values > 0:")'
 printf,lun,Q4bV[5],Q4bV[6],Q4bV[7],format='(I-6,I-6,I-6,T45,"Number of ON-OFF values above 1,2,3 sigma")'
 printf,lun,Q4bV[8],Q4bV[9],format='(F-4.2,2x,F-4.2,T45,"Min and Max Threshold value where ON-OFF > 0")'
 printf,lun,Q4bV[10],format='(F-6.2,T45,"% of ON-OFF values between min and max threshold & >0")'
 printf,lun,Q4bV[11],format='(F-6.2,T45,"% of ON-OFF values lt 0 and > max threshold")'
 printf,lun,Q4bV[12],format='(D11.8,T45,"Probability of False Detection")'
 printf,lun,'*********************************'
 printf,lun,'    Q4c: Peak Asymmetry          '
 printf,lun,'*********************************'
 printf,lun,Q4cV[0],format='(F-8.2,T45,"Integral of ON-OFF curve")'
 printf,lun,Q4cV[1],format='(F-6.3,T45,"Ratio of ON-OFF Integral to 1 sigma Gaussian")'
 printf,lun,Q4cV[4],format='(F-6.2,T45,"% of ON-OFF values > 0:")'
 printf,lun,Q4cV[5],Q4cV[6],Q4cV[7],format='(I-6,I-6,I-6,T45,"Number of ON-OFF values above 1,2,3 sigma")'
 printf,lun,Q4cV[8],Q4cV[9],format='(F-4.2,2x,F-4.2,T45,"Min and Max Threshold value where ON-OFF > 0")'
 printf,lun,Q4cV[10],format='(F-6.2,T45,"% of ON-OFF values between min and max threshold & >0")'
 printf,lun,Q4cV[11],format='(F-6.2,T45,"% of ON-OFF values lt 0 and > max threshold")'
 printf,lun,Q4cV[12],format='(D11.8,T45,"Probability of False Detection")'
 printf,lun,'*********************************'
 printf,lun,'    Q4d: Power Asymmetry          '
 printf,lun,'*********************************'
 printf,lun,Q4dV[0],format='(F-8.2,T45,"Integral of ON-OFF curve")'
 printf,lun,Q4dV[1],format='(F-6.3,T45,"Ratio of ON-OFF Integral to 1 sigma Gaussian")'
 printf,lun,Q4dV[4],format='(F-6.2,T45,"% of ON-OFF values > 0:")'
 printf,lun,Q4dV[5],Q4dV[6],Q4dV[7],format='(I-6,I-6,I-6,T45,"Number of ON-OFF values above 1,2,3 sigma")'
 printf,lun,Q4dV[8],Q4dV[9],format='(F-4.2,2x,F4.2,T45,"Min and Max Threshold value where ON-OFF > 0")'
 printf,lun,Q4dV[10],format='(F-6.2,T45,"% of ON-OFF values between min and max threshold & >0")'
 printf,lun,Q4dV[11],format='(F-6.2,T45,"% of ON-OFF values lt 0 and > max threshold")'
 printf,lun,Q4dV[12],format='(D11.8,T45,"Probability of False Detection")'
 printf,lun,'*********************************'
 printf,lun,'    Q4e: Peak Offset             '
 printf,lun,'*********************************'
 printf,lun,Q4eV[0],format='(F-8.2,T45,"Integral of ON-OFF curve")'
 printf,lun,Q4eV[1],format='(F-6.3,T45,"Ratio of ON-OFF Integral to 1 sigma Gaussian")'
 printf,lun,Q4eV[4],format='(F-6.2,T45,"% of ON-OFF values > 0:")'
 printf,lun,Q4eV[5],Q4eV[6],Q4eV[7],format='(I-6,I-6,I-6,T45,"Number of ON-OFF values above 1,2,3 sigma")'
 printf,lun,Q4eV[8],Q4eV[9],format='(F-4.2,2x,F-4.2,T45,"Min and Max Threshold value where ON-OFF > 0")'
 printf,lun,Q4eV[10],format='(F-6.2,T45,"% of ON-OFF values between min and max threshold & >0")'
 printf,lun,Q4eV[11],format='(F-6.2,T45,"% of ON-OFF values lt 0 and > max threshold")'
 printf,lun,Q4eV[12],format='(D11.8,T45,"Probability of False Detection")'
 printf,lun,'*********************************'
 printf,lun,'    Q4f: Power Offset             '
 printf,lun,'*********************************'
 printf,lun,Q4fV[0],format='(F-8.2,T45,"Integral of ON-OFF curve")'
 printf,lun,Q4fV[1],format='(F-6.3,T45,"Ratio of ON-OFF Integral to 1 sigma Gaussian")'
 printf,lun,Q4fV[4],format='(F-6.2,T45,"% of ON-OFF values > 0:")'
 printf,lun,Q4fV[5],Q4fV[6],Q4fV[7],format='(I-6,I-6,I-6,T45,"Number of ON-OFF values above 1,2,3 sigma")'
 printf,lun,Q4fV[8],Q4fV[9],format='(F-4.2,2x,F-4.2,T45,"Min and Max Threshold value where ON-OFF > 0")'
 printf,lun,Q4fV[10],format='(F-6.2,T45,"% of ON-OFF values between min and max threshold & >0")'
 printf,lun,Q4fV[11],format='(F-6.2,T45,"% of ON-OFF values lt 0 and > max threshold")'
 printf,lun,Q4fV[12],format='(D11.8,T45,"Probability of False Detection")'
 If Q4fV[1] ge 0.85 then Detect = 1 else Detect = 0 
 printf,lun,Detect,format='(I-,T45,"Possible Detection (Yes=1 or NO=0)")'
 close,lun
 FREE_LUN, lun
endif ;write 

If keyword_set(READ) then begin
 ;---------------
 ;Open
 ;--------------
 ;openr,lun,filename+'_Output.txt',/get_lun
 openr,lun,filename+'.txt',/get_lun

 ;Intialize 
 filename=''  &Date = ''
 S = ''       &fmin =0.0
 fmax = 0.0   &Tbeam = 0
 Sbeam = 0    &skip='' &temp=''
 Q_test = 0.0d & Q_test1 = 0.0d & Q_test2 = 0.0d
 n = 13     ;number of elements
 Q1aV = fltarr(4) & Q1bV = fltarr(4) &Q4av = fltarr(n) & Q4bV = fltarr(n) & Q4cV = fltarr(n)
 Q4dV = fltarr(n) & Q4eV = fltarr(n) & Q4fV = fltarr(n) 

 ;------------------------------
 readf,lun,skip
 readf,lun,skip
 readf,lun,filename,format='(A-)'
 readf,lun,Date,skip,format='(A-8,T45,A)'
 readf,lun,S,skip,format='(A-3,T45,A)'
 readf,lun,fmin,skip,format='(F-6.2,T45,A)'
 readf,lun,fmax,skip,format='(F-6.2,T45,A)'
 readf,lun,TBeam,skip,format='(I-,T45,A)'
 readf,lun,SBeam,skip,format='(I-,T45,A)'
 readf,lun,skip
 readf,lun,skip
 readf,lun,skip;'*********************************'
 readf,lun,skip;'       Q1a: Time Series          '
 readf,lun,skip;'*********************************'
 readf,lun,Q_test,skip,format='(F,T45,A)'
 Q1aV[0] = Q_test
 readf,lun,Q_test,skip,format='(F5.1,T45,A)'
 Q1aV[1] = Q_test
 readf,lun,Q_test,skip,format='(F,T45,A)'
 Q1aV[2] = Q_test
 readf,lun,Q_test,skip,format='(D11.8,T45,A)'
 Q1aV[3] = Q_test
 readf,lun,skip;'*********************************'
 readf,lun,skip;'    Q1b: Integrated Spectrum     '
 readf,lun,skip;'*********************************'
 readf,lun,Q_test,skip,format='(F,T45,A)'
 Q1bV[0] = Q_test
 readf,lun,Q_test,skip,format='(F5.1,T45,A)'
 Q1bV[1] = Q_test
 readf,lun,Q_test,skip,format='(F,T45,A)'
 Q1bV[2] = Q_test
 readf,lun,Q_test,skip,format='(D11.8,T45,A)'
 Q1bV[3] = Q_test
 readf,lun,skip;'*********************************'
 readf,lun,skip;'    Q4a: Number of Peaks         '
 readf,lun,skip;'*********************************'
 readf,lun,Q_test,skip,format='(F-8.2,T45,A)'
 Q4aV[0] = Q_test
 readf,lun,Q_test,skip,format='(F-6.2,T45,A)'
 Q4aV[1] = Q_test
 Q4aV[2] = Q_test/2.0d
 Q4aV[3] = Q_test2/3.0d
 readf,lun,Q_test,skip,format='(F-6.2,T45,A)'
 Q4aV[4] = Q_test
 readf,lun,Q_test,Q_test1,Q_test2,skip,format='(I-6,I-6,I-6,T45,A)'
 Q4aV[5] = Q_test
 Q4aV[6] = Q_test1
 Q4aV[7] = Q_test2
 Q_Test1 = 0.0
 readf,lun,Q_test,Q_test1,skip,format='(F-4.2,2x,F-4.2,T45,A)'
 Q4aV[8] = Q_test
 Q4aV[9] = Q_test1
 readf,lun,Q_test,skip,format='(F-6.2,T45,A)'
 Q4aV[10] = Q_test
 readf,lun,Q_test,skip,format='(F-6.2,T45,A)'
 Q4aV[11] = Q_test
 readf,lun,Q_test,format='(D11.8,T45)'
 Q4aV[12] = Q_test
 readf,lun,skip;'*********************************'
 readf,lun,skip;'    Q4b: Power of Peaks          '
 readf,lun,skip;'*********************************'
 readf,lun,Q_test,skip,format='(F-8.2,T45,A)'
 Q4bV[0] = Q_test
 readf,lun,Q_test,skip,format='(F-6.2,T45,A)'
 Q4bV[1] = Q_test
 Q4bV[2] = Q_test/2.
 Q4bV[3] = Q_test/3.
 readf,lun,Q_test,skip,format='(F-6.2,T45,A)'
 Q4bV[4] = Q_test
 readf,lun,Q_test,Q_test1,Q_test2,skip,format='(I-6,I-6,I-6,T45,A)'
 Q4bV[5] = Q_test
 Q4bV[6] = Q_test1
 Q4bV[7] = Q_test2
 Q_Test1 = 0.0
 readf,lun,Q_test,Q_test1,skip,format='(F-4.2,2x,F-4.2,T45,A)'
 Q4bV[8] = Q_test
 Q4bV[9] = Q_test1
 readf,lun,Q_test,skip,format='(F-6.2,T45,A)'
 Q4bV[10] = Q_test
 readf,lun,Q_test,skip,format='(F-6.2,T45,A)'
 Q4bV[11] = Q_test
 readf,lun,Q_test,format='(D11.8,T45)'
 Q4bV[12] = Q_test
 readf,lun,skip;'*********************************'
 readf,lun,skip;'    Q4c: Peak Asymmetry          '
 readf,lun,skip;'*********************************'
 readf,lun,Q_test,skip,format='(F-8.2,T45,A)'
 Q4cV[0] = Q_test
 readf,lun,Q_test,skip,format='(F-6.2,T45,A)'
 Q4cV[1] = Q_test
 Q4cV[2] = Q_test/2.0d
 Q4cV[3] = Q_test/3.0d
 readf,lun,Q_test,skip,format='(F-6.2,T45,A)'
 Q4cV[4] = Q_test
 readf,lun,Q_test,Q_test1,Q_test2,skip,format='(I-6,I-6,I-6,T45,A)'
 Q4cV[5] = Q_test
 Q4cV[6] = Q_test1
 Q4cV[7] = Q_test2
 Q_Test1 = 0.0
 readf,lun,Q_test,Q_test1,skip,format='(F-4.2,2x,F-4.2,T45,A)'
 Q4cV[8] = Q_test
 Q4cV[9] = Q_test1
 readf,lun,Q_test,skip,format='(F-6.2,T45,A)'
 Q4cV[10] = Q_test
 readf,lun,Q_test,skip,format='(F-6.2,T45,A)'
 Q4cV[11] = Q_test
 readf,lun,Q_test,format='(D11.8,T45)'
 Q4cV[12] = Q_test
 readf,lun,skip;'*********************************'
 readf,lun,skip;'    Q4d: Power Asymmetry          '
 readf,lun,skip;'*********************************'
 readf,lun,Q_test,skip,format='(F-8.2,T45,A)'
 Q4dV[0] = Q_test
 readf,lun,Q_test,skip,format='(F-6.2,T45,A)'
 Q4dV[1] = Q_test
 Q4dV[2] = Q_test/2.0d
 Q4dV[3] = Q_test/2.0d
 readf,lun,Q_test,skip,format='(F-6.2,T45,A)'
 Q4dV[4] = Q_test
 readf,lun,Q_test,Q_test1,Q_test2,skip,format='(I-6,I-6,I-6,T45,A)'
 Q4dV[5] = Q_test
 Q4dV[6] = Q_test1
 Q4dV[7] = Q_test2
 Q_Test1 = 0.0
 readf,lun,Q_test,Q_test1,skip,format='(F-4.2,2x,F-4.2,T45,A)'
 Q4dV[8] = Q_test
 Q4dV[9] = Q_test1
 readf,lun,Q_test,skip,format='(F-6.2,T45,A)'
 Q4dV[10] = Q_test
 readf,lun,Q_test,skip,format='(F-6.2,T45,A)'
 Q4dV[11] = Q_test
 readf,lun,Q_test,format='(D11.8,T45)'
 Q4dV[12] = Q_test
 readf,lun,skip;'*********************************'
 readf,lun,skip;'    Q4e: Peak Offset             '
 readf,lun,skip;'*********************************'
 readf,lun,Q_test,skip,format='(F-8.2,T45,A)'
 Q4eV[0] = Q_test
 readf,lun,Q_test,skip,format='(F-6.2,T45,A)'
 Q4eV[1] = Q_test
 Q4eV[2] = Q_test/2.0d
 Q4eV[3] = Q_test/3.0d
 readf,lun,Q_test,skip,format='(F-6.2,T45,A)'
 Q4eV[4] = Q_test
 readf,lun,Q_test,Q_test1,Q_test2,skip,format='(I-6,I-6,I-6,T45,A)'
 Q4eV[5] = Q_test
 Q4eV[6] = Q_test1
 Q4eV[7] = Q_test2
 Q_Test1 = 0.0
 readf,lun,Q_test,Q_test1,skip,format='(F-4.2,2x,F-4.2,T45,A)'
 Q4eV[8] = Q_test
 Q4eV[9] = Q_test1
 readf,lun,Q_test,skip,format='(F-6.2,T45,A)'
 Q4eV[10] = Q_test
 readf,lun,Q_test,skip,format='(F-6.2,T45,A)'
 Q4eV[11] = Q_test
 readf,lun,Q_test,format='(D11.8,T45)'
 Q4eV[12] = Q_test
 readf,lun,skip;'*********************************'
 readf,lun,skip;'    Q4f: Power Offset             '
 readf,lun,skip;'*********************************'
  readf,lun,Q_test,skip,format='(F-8.2,T45,A)'
 Q4fV[0] = Q_test
 readf,lun,Q_test,skip,format='(F-6.2,T45,A)'
 Q4fV[1] = Q_test
 Q4fV[2] = Q_test/2.0d
 Q4fV[3] = Q_test/2.0d
 readf,lun,Q_test,skip,format='(F-6.2,T45,A)'
 Q4fV[4] = Q_test
 readf,lun,Q_test,Q_test1,Q_test2,skip,format='(I-6,I-6,I-6,T45,A)'
 Q4fV[5] = Q_test
 Q4fV[6] = Q_test1
 Q4fV[7] = Q_test2
 Q_Test1 = 0.0
 readf,lun,Q_test,Q_test1,skip,format='(F-4.2,2x,F-4.2,T45,A)'
 Q4fV[8] = Q_test
 Q4fV[9] = Q_test1
 readf,lun,Q_test,skip,format='(F-6.2,T45,A)'
 Q4fV[10] = Q_test
 readf,lun,Q_test,skip,format='(F-6.2,T45,A)'
 Q4fV[11] = Q_test
 readf,lun,Q_test,format='(D11.8,T45)'
 Q4fV[12] = Q_test
 readf,lun,Q_test,skip,format='(I-,T45,A)'
 Detect = Q_test
 close,lun
 FREE_LUN, lun
endif ;end read 

end
