;plot_textfile.pro
;Plots of PP tex files
;
;Input:  list     ;filename containing a list of PP txt files  
;
;Opitonal Input:         
;                  sfilename; Sky filename (if the sky filename is input than the program will take the ratio of all Q values) 
;
;OUTPUT: 
;       N        ;Number of txt files (see filename+'_numbers.txt' for each N represents)
;       Date     ;Date [N]
;       Pol      ;Polarization [N]
;       fmin     ;fmin [N]
;       fmax     ;fmax [N]
;       TBeam    ;Target Beam [N] (ranges from 0-2)
;       Sbeam    ;Target Beam [N] (ranges from 0-2)
;       Q1aV     ;Q1a quantitative values [N,4]
;       Q1bV     ;Q1b quantitative values [N,4]
;       Q4aV     ;Q4a quantitative values [N,14]
;       Q4bV     ;Q4b quantitative values [N,14]
;       Q4cV     ;Q4c quantitative values [N,14]
;       Q4dV     ;Q4d quantitative values [N,14]
;       Q4eV     ;Q4e quantitative values [N,14]
;       Q4fV     ;Q4f quantitative values [N,14]
;
pro plot_textfile,filename,Sfilename=Sfilename,Date=Date,Pol=S,fmin=fmin,fmax=fmax,TBeam=TBeam,Sbeam=Sbeam,$
    Q1aV=Q1a,Q1bV=Q1b,Q4aV=Q4a,Q4bV=Q4b,Q4cV=Q4c,Q4dV=Q4d,Q4eV=Q4e,Q4fV=Q4f,N=x,PROB=PROB

Factor =  2.0d  ; factor that Q1a and Q1b are greater than

ratio = 0 
;If keyword_set(Sfilename) then ratio = 'Sky ratio' else ratio=''
;If keyword_set(Sfilename) then skyratio = 1                         ;If sfilename is set, take the skyratio 
;If keyword_set(skyratio) then ytitle = 'Sky Ratio'

n = FILE_LINES(filename)
file = STRARR(n)
 OPENR, unit, filename,/GET_LUN
 READF, unit, file
 FREE_LUN, unit

;intialize variables 
x = indgen(n)
n_values = 13   ;num values for Q4i 
Date=strarr(n) & S=strarr(n)	&fmin=dblarr(n)
fmax=dblarr(n) & Tbeam=dblarr(n)&SBeam=dblarr(n)
Q1a=dblarr(n,4) &Q1b=dblarr(n,4)	
Q4a=dblarr(n,n_values)
Q4b=dblarr(n,n_values)	&Q4c=dblarr(n,n_values)	&Q4d=dblarr(n,n_values)
Q4e=dblarr(n,n_values)	&Q4f=dblarr(n,n_values)

;---------------------------------------------------------------------
;sky variables (set all to 1 and then overide if sky ratio is set
;-------------- ------------------------------------------------------
Tbeam_s=dblarr(n)&     SBeam_s=dblarr(n)
Q1a_S=dblarr(n,4)+1.0d &Q1b_s=dblarr(n,4) +1.0d	
Q4a_s=dblarr(n,n_values)+1.0d & Q4b_s=dblarr(n,n_values)+1.0d  & Q4c_s=dblarr(n,n_values)+1.0d
Q4d_s=dblarr(n,n_values)+1.0d & Q4e_s=dblarr(n,n_values)+1.0d  & Q4f_s=dblarr(n,n_values)+1.0d


;-------------------------
;read all data for Target
;-------------------------
for i=0,n-1 do begin 
 print, file[i]
 lofar_postprocess_textfilev2,file[i],Date=Da,S=Sout,fmin=fmi,fmax=fma,BT=BT,BSky=BSky,$
    	Q1aV=Q1aV,Q1bV=Q1bV,Q4aV=Q4aV,Q4bV=Q4bV,Q4cV=Q4CV,Q4dV=Q4dV,Q4eV=Q4eV,Q4fV=Q4fV,$
        /READ
  date[i]  = da
  S[i]     = Sout
  fmin[i]  = fmi
  fmax[i]  = fma
  TBeam[i] = BT
  Sbeam[i] = Bsky
  Q1a[i,*] = Q1aV
  Q1b[i,*] = Q1bV
  Q4a[i,*] = Q4aV
  Q4b[i,*] = Q4bV
  Q4c[i,*] = Q4cV
  Q4d[i,*] = Q4dV
  Q4e[i,*] = Q4eV
  Q4f[i,*] = Q4fV
endfor 

filename = strmid(filename,0,STRPOS(filename, '.txt'))

OPENW, unit, filename+'_numbers.txt',/GET_LUN
print_array = strarr(7,n)
print_array[0,*] = x     ;# 
print_array[1,*] =date   ;Date
print_array[2,*] =S      ;Polarization 
print_array[3,*] =cgNumber_Formatter(fmin,decimals=1)  ;Fmin
print_array[4,*] =cgNumber_Formatter(fmax,decimals=1)   ;Fmax
print_array[5,*] =UINT(Tbeam)  ;Target Beam #
print_array[6,*] =UINT(Sbeam)  ;Sky Beam # 
printf,unit,print_array
FREE_LUN, unit

;output file name with detections 
If not(keyword_set(Sfilename)) then OPENW, unit, filename+'_Detection.txt',/GET_LUN
If keyword_set(Sfilename) then OPENW, unit, filename+'_Detection_DivSky.txt',/GET_LUN

;--------------
;   Read Sky
;-------------
If keyword_set(Sfilename) then begin
 file_s = STRARR(n)
 OPENR, unit2, Sfilename,/GET_LUN
 READF, unit2, file_s
 FREE_LUN, unit2

 for i=0,n-1 do begin 
  lofar_postprocess_textfilev2,file_s[i],Date=D,S=Sout,fmin=fmi,fmax=fma,BT=BT,BSky=BSky,$
    	Q1aV=Q1aV,Q1bV=Q1bV,Q4aV=Q4aV,Q4bV=Q4bV,Q4cV=Q4CV,Q4dV=Q4dV,Q4eV=Q4eV,Q4fV=Q4fV,$
        /READ
  S[i]     = Sout
  fmin[i]  = fmi
  fmax[i]  = fma
  TBeam_S[i] = BT
  Sbeam_s[i] = Bsky
  Q1a_s[i,*] = Q1aV
  Q1b_s[i,*] = Q1bV
  Q4a_s[i,*] = Q4aV
  Q4b_s[i,*] = Q4bV
  Q4c_s[i,*] = Q4cV
  Q4d_s[i,*] = Q4dV
  Q4e_s[i,*] = Q4eV
  Q4f_S[i,*] = Q4fV
 endfor 
endif ;end reading sky file

If not(keyword_set(Sfilename)) then filename=filename
If keyword_set(Sfilename) then filename = filename+'_DivSky'
cgPS_Open,filename+'.ps'

;--------------------------------------------------------
;--------------------------------------------------------
;                       Q1a (Time-Series)
;--------------------------------------------------------
;--------------------------------------------------------

sky = Q1a_s[*,*]
y = Q1a[*,*];/(Q1a_s[*,*] + (Q1a_s[*,*] eq 0))
q1_textplot,x,y,unit=unit,ytitle=ytitle,title=title,skyratio=skyratio,factor=factor,print_array=print_array,/Q1a,PROB=PROB
q1_textplot,x,y,unit=unit,ytitle=ytitle,title=title,skyratio=skyratio,factor=factor,print_array=print_array,/Q1a,PROB=PROB,/COMBINEHISTO,Psky=Q1a_s[*,*]

;--------------------------------------------------------
;--------------------------------------------------------
;                       Q1b (int spectrum)
;--------------------------------------------------------
;--------------------------------------------------------

y = Q1b[*,*];/(Q1b_s[*,*] + (Q1b_s[*,*] eq 0))
q1_textplot,x,y,unit=unit,ytitle=ytitle,title=title,skyratio=skyratio,factor=factor,print_array=print_array,/Q1b,PROB=PROB
q1_textplot,x,y,unit=unit,ytitle=ytitle,title=title,skyratio=skyratio,factor=factor,print_array=print_array,/Q1b,PROB=PROB,/COMBINEHISTO,Psky=Q1b_s[*,*]


;--------------------------------------------------------
;--------------------------------------------------------
;                       Q4 
;--------------------------------------------------------
;--------------------------------------------------------

;------------------
; Q4a
;------------------
q4_textplot,x,Q4a,'a',Q_sky=Q4a_s,skyratio=skyratio,print_array=print_array,unit=unit,PROB=PROB

;------------------
; Q4b
;------------------
q4_textplot,x,Q4b,Q_sky=Q4b_s,'b',skyratio=skyratio,print_array=print_array,unit=unit,PROB=PROB

;------------------
; Q4c
;------------------
q4_textplot,x,Q4c,Q_sky=Q4c_s,'c',skyratio=skyratio,print_array=print_array,unit=unit,PROB=PROB

;------------------
; Q4d
;------------------
q4_textplot,x,Q4d,Q_sky=Q4d_s,'d',skyratio=skyratio,print_array=print_array,unit=unit,PROB=PROB

;------------------
; Q4e
;------------------
q4_textplot,x,Q4e,Q_sky=Q4e_s,'e',skyratio=skyratio,print_array=print_array,unit=unit,PROB=PROB

;------------------
; Q4f
;------------------
q4_textplot,x,Q4f,Q_sky=Q4f_s,'f',skyratio=skyratio,print_array=print_array,unit=unit,PROB=PROB,/COMBINEHISTO,Psky=Q4f_S[*,*]


;q4_textplot,x,Q4f,Q_sky=Q4f_s,'f',skyratio=skyratio,print_array=print_array,unit=unit,PROB=PROB,/COMBINEHISTO


FREE_LUN, unit
cgPS_Close
cgfixps,filename+'.ps' 
ps_pdf
spawn,'rm '+filename+'.ps'

end

