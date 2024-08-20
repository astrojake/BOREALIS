;This program analyzes the PP text files for possible detections
;The input to this program is from plot_texfile 
;
;Input:  filename        ;filename containing a list of PP txt files
;        sfilename       ;Sky filename containing a list of PP txt files (same length as list)
;         
;Optional Input Flag:
;        PROB          ;Run probability 
;        SAVE          ;Save output
;
;Outputs:         
;What: finds the mean, median, stdev, max, min of each Qvalue [5]

pro analyze_pp,filename,sfilename=sfilename,PROB=PROB,Q1a_pp=Q1a_pp,$
    Q1b_pp=Q1b_pp,Q4a_pp=Q4a_pp,Q4b_pp=Q4b_pp,Q4c_pp=Q4c_pp,Q4d_pp=Q4d_pp,Q4e_pp=Q4e_pp,Q4f_pp=Q4f_pp,SAVE=SAVE

;Read in all text files and create plots and list possible detections
  plot_textfile,filename,Sfilename=Sfilename,Date=Date,Pol=S,fmin=fmin,fmax=fmax,TBeam=TBeam,Sbeam=Sbeam,$
    Q1aV=Q1a,Q1bV=Q1b,Q4aV=Q4a,Q4bV=Q4b,Q4cV=Q4c,Q4dV=Q4d,Q4eV=Q4e,Q4fV=Q4f,N=x,PROB=PROB
   ;e.g. Q1a [[N,4], Q4a [N,14] 

 ;init
 n_Q1 = 4 & n_Q4 = 13
 n_moments = 5
 Q1a_pp  = dblarr(n_moments,n_Q1)   & Q1b_pp =  dblarr(n_moments,n_Q1) 
 Q4a_pp  = dblarr(n_moments,n_Q4)   & Q4b_pp =  dblarr(n_moments,n_Q4)
 Q4c_pp  = dblarr(n_moments,n_Q4)   & Q4d_pp =  dblarr(n_moments,n_Q4)
 Q4e_pp  = dblarr(n_moments,n_Q4)   & Q4f_pp =  dblarr(n_moments,n_Q4)

 ;-----------------------------
 ;-------- Q1a and Q1b --------
 ;-----------------------------
 for i=0,n_Q1-1 do begin 
  find_moments,Q1a[*,i],data_out=data_out
  Q1a_pp[*,i] = data_out
  find_moments,Q1b[*,i],data_out=data_out
  Q1b_pp[*,i] = data_out
 endfor 
 
 ;-----------------------------
 ;--------  Q4a-Q4f    --------
 ;-----------------------------
 for i=0,n_Q4-1 do begin
   find_moments,Q4a[*,i],data_out=data_out
   Q4a_pp[*,i] = data_out
   find_moments,Q4b[*,i],data_out=data_out
   Q4b_pp[*,i] = data_out
   find_moments,Q4c[*,i],data_out=data_out
   Q4c_pp[*,i] = data_out
   find_moments,Q4d[*,i],data_out=data_out
   Q4d_pp[*,i] = data_out
   find_moments,Q4e[*,i],data_out=data_out
   Q4e_pp[*,i] = data_out
   find_moments,Q4f[*,i],data_out=data_out
   Q4f_pp[*,i] = data_out
 endfor
 
 save,Date,S,fmin,fmax,Tbeam,Sbeam,Q1a,Q1b,Q4a,Q4b,Q4c,Q4d,Q4e,Q4f,x,Q1a_pp,Q1b_pp,Q4a_pp,Q4b_pp,Q4c_pp,Q4d_pp,Q4e_pp,Q4f_pp,filename=filename+'.sav'

 
 OPENW, unit, filename+'_moments.txt',/GET_LUN
 ;Print Q1a
 print_array = strarr(1,5)
 print_array[0] = 'Mean' & print_array[1] = 'Median'
 print_array[2] = 'Max'  & print_array[3] = 'Min' & print_array[4] = 'STDEV'
 
 printf,unit,'**************************'
 printf,unit,'*** Q1a (time-series)  ***'
 printf,unit,'**************************'
  printf,unit,''
 print_head = strarr(n_Q1) 
 print_head[0] = 'Integral' & print_head[1] = 'ON-OFF >0' & print_head[2] = 'ON-OFF STDEV' 
 print_head[3] = 'False-Postive Probability'
 for i=0,n_Q1-1 do begin 
   print_all = strarr(2,5)
   print_all[*] = [print_array,string(transpose(cgNumber_Formatter(Q1a_pp[*,i],DECIMALS=2)))]
   printf,unit,'**********'
   printf,unit,print_head[i]
   printf,unit,'**********'
   printf,unit,print_all;,format='(A,2x,A)'
   printf,unit,''
 endfor 
 
 printf,unit,'**************************'
 printf,unit,'*** Q1b (time-series)  ***'
 printf,unit,'**************************'
 printf,unit,''
 print_head = strarr(n_Q1)
 print_head[0] = 'Integral' & print_head[1] = 'ON-OFF >0' & print_head[2] = 'ON-OFF STDEV'
 print_head[3] = 'False-Postive Probability'
 for i=0,n_Q1-1 do begin
   print_all = strarr(2,5)
   print_all[*] = [print_array,string(transpose(cgNumber_Formatter(Q1b_pp[*,i],DECIMALS=2)))]
   printf,unit,'**********'
   printf,unit,print_head[i]
   printf,unit,'**********'
   printf,unit,print_all;,format='(A,2x,A)'
   printf,unit,''
 endfor
 
 FREE_LUN, unit

end
