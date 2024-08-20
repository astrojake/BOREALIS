;***************************************************************
;      Read Lofar data and make quicklook
;**************************************************************
;Uses: hdf5 programs
;***************************************************************
;Author: Jake Turner (UVA)
;        J.M. Griessmeier (CNRS) and Philippe Zarka (Obspm)
;        
;Things that you can change:
;   fmin
;   fmax
;   start_time
;   N_time
;   steps
;   file_root
;   date

st= SYSTIME(1)
j_name = 'time_'+strtrim(string(st),1)+'.dat'
journal,j_name

;*****************************************************
;                  Inputs
;*****************************************************
;Ep Andr 
fmin = 14.8d  ;MHz
fmax   = 62.3d   ;MHz

;55 Cnc 
;fmin = 26.07d
;fmax = 73.7d

;V830 
;fmin = 109.87
;fmax = 189.93
;fmin = 109.86328125d   ;MHz    ;**input**
;   fmax   = 189.993835d   ;MHz   ;**input**
;---------
;  Time 
;---------
start_time = 0.00
tmin = start_time

;N_time = 1800.d        ;L545209
;N_time = 1799.85       ;L554093
;N_time = 599.0d        ;L545211

;55 Cnc
;steps = 16               ; L429868
;steps  = 6
 
;Eps 
;steps  = 10.0d          ; L545209
;steps  = 1.0d                ; L545211

;V830
N_time = 400.0 
steps = 1.

total_t = N_time*steps  ;total time requested

;file_root = '/data/jake.turner/exoplanet/LC5_DDT_002/'
;date = 'L527649'

;file_root = '/data/jake.turner/exoplanet/LC6_010/'
;date = 'L545209'
;date = 'L545211'
;date = 'L545213'

;file_root = '/data/jake.turner/exoplanet/LC7_015/'
;date = 'L567727'
;date = 'L573685'
;
;tau
file_root = '/data/jake.turner/exoplanet/LC7_013/'
date = 'L568467'
;
;These strings are for Beam and Pol for loops
;strings for information
B_start = 0    ;Beam Start
B_end   = 0    ;Beam  End
S_start = 0   ;Pol Start
S_end   = 0   ;Pol End

plot_raw  = 1  ; 0=NO, 1= YES  ;plot the raw data
data_norm = 1 ; 0=NO, 1= YES   ;normalize the data by the 10% quantile and divide out slope 

;**************************************
;Don't need to change anything below here
;*********************************************
B_start_string = strtrim(string(b_start),1)
S_start_string = strtrim(string(S_start),1)
B_end_string = strtrim(string(b_end),1)
S_start_string = strtrim(string(S_start),1)
S_end_string = strtrim(string(S_end),1)

If S_start eq 4 then S_start_string = strtrim(string(3),1)
If S_end eq 4 then S_end_string = strtrim(string(3),1)
IF s_start eq 4 then S_start = 4
IF s_end eq 4 then S_end = 4

;-----------------------------
;       Create file name
;-----------------------------
for i=0L,0 do begin  ;loop for dates  
  ;----------------------------------------------------------
  ;        Start loops for Beam and Polarization
  ;----------------------------------------------------------
  ;loops  beams, and polarization
  for j=b_start,b_end do begin  ; loop for Beam: B000-
    for k=s_start,s_end do begin  ; loop for Polarization:S0-S3
      j_string = strtrim(string(j),1)
      k_string = strtrim(string(k),1)
      If k eq 4 then k_string = strtrim(string(3),1) ;V/I
      ;----------------------------------------------------------
      ;        Read Data for specific Beam and Polarization
      ;----------------------------------------------------------
      ;file name creation
      file = file_root+date[i]+'/'+date[i]+'_SAP000_B00'+j_string+$
        '_S'+k_string+'_P000_bf'
        
      
      print, file
      set_plot,'PS'
      device,filename=file_root+date[i]+'/'+date[i]+'_SAP000_B00'+j_string+$
        '_S'+k_string+'.ps'
      ;Time loop
      tmin = start_time
      for kkk=0L,steps-1 do begin
        
        tmax =  tmin + N_time  ;set end time 
        If k eq 0 then READ_H5_DATA, file,tmin,tmax,fmin,fmax,j,k,x,t,f,nt,nf;,nnt,nnf,dt;,/VERB;,/NSPECTRA
        If k eq 3 then READ_H5_DATA, file,tmin,tmax,fmin,fmax,j,k,x,t,f,nt,nf;,nnt,nnf,dt,POL=k_string;,/VERB;,/NSPECTRA
        
        
        If k eq 4 then READ_H5_DATA, file_root+date[i]+'/'+date[i]+'_SAP000_B00'+j_string+$
        '_S0_P000_bf',tmin,tmax,fmin,fmax,j_string,x_I,t,f,nt,nf,nnt,nnf,dt;,/VERB;,/NSPECTRA
        If k eq 4 then READ_H5_DATA, file,tmin,tmax,fmin,fmax,j_string,x_V,t,f,nt,nf,nnt,nnf,dt,POL=k_string;,/VERB;,/NSPECTRA
        
;        save,x,t,f,filename='quicklook_data.sav'
;        print, 'saved'
;        stop
        ;Reduce size
        nnv=nt/2000.        ; for visualization in ps file
        if nnv ne long(nnv) then nnv=long(nnv)+1 else nnv=long(nnv)
        nntv=long(nt/nnv)
        mmv = nf/2000.
        if mmv ne long(mmv) then mmv=long(mmv)+1 else mmv=long(mmv)
        mmtv=long(nf/mmv)

        nwin=0 & xp=100 & yp=100
  !p.multi=[0,1,1,1]      
  If plot_raw eq 1 and k ne 4 then begin
        SPDYNPS,rebin(x(0:nntv*nnv-1,0:mmtv*mmv-1),nntv,mmtv),min(t),max(t), $
          min(f),max(f), $
          'Time [s]','Frequency (MHz)','Raw: '+date+' Beam:'+j_string+' P:'+k_string, $
          0,0,0,.05,.95,0,'.'
  Endif
  
  help,x
  print, min(x)
 ; print,x[34370:34382,0]
  device,/close
  stop
  
  If plot_raw eq 1 and k eq 4 then begin
        SPDYNPS,rebin(x_V(0:nntv*nnv-1,0:mmtv*mmv-1),nntv,mmtv),min(t),max(t), $
          min(f),max(f), $
          'Time [s]','Frequency (MHz)','Raw: '+date+' Beam:'+j_string+' P:'+k_string, $
          0,0,0,.05,.95,0,'.'
  Endif
  
 If k eq 3 then begin ; V pol
   x = x*x  ;V^2 
   SPDYNPS,rebin(x(0:nntv*nnv-1,0:mmtv*mmv-1),nntv,mmtv),min(t),max(t), $
     min(f),max(f), $
     'Time [s]','Frequency (MHz)','Raw (V^2): '+date+' Beam:'+j_string+' P:'+k_string, $
     0,0,0,.05,.95,0,'.'
 Endif
 
 If k eq 4 then begin ; V pol
   x = x_V/x_I  ;V/I
   SPDYNPS,rebin(x(0:nntv*nnv-1,0:mmtv*mmv-1),nntv,mmtv),min(t),max(t), $
     min(f),max(f), $
     'Time [s]','Frequency (MHz)','Raw (V/I): '+date+' Beam:'+j_string+' P:'+k_string, $
     0,0,0,.05,.95,0,'.' 
   STANDARD_PLOTS, x,x,t,f,'V/I'
   
   y = (x - mean(x))/stddev(x)  ;V/I
   SPDYNPS,rebin(y(0:nntv*nnv-1,0:mmtv*mmv-1),nntv,mmtv),min(t),max(t), $
     min(f),max(f),$
     'Time [s]','Frequency (MHz)','Raw ( vprime): '+date+' Beam:'+j_string+' P:'+k_string, $
     0,0,0,.05,.95,0,'.'
   STANDARD_PLOTS, y,y,t,f,'v prime'  
     
 Endif
  
  If data_norm eq 1 then begin
        print,'-------------------------'
        print, 'Start Normalize Data'
        print,'-------------------------'
        ;-------------------
        ; Normalize the Data
        ;-------------------
        data_gain = fltarr(nf)
        upper_quantile = 0.1d  
        for iii=0L,nf-1 do data_gain[iii] = dyn_n(x(*,iii),upper_quantile)
       
         x = x/rebin(reform(data_gain,1,nf),nt,nf)  
         
         !p.multi=[0,1,2]
         SPDYNPS,rebin(x(0:nntv*nnv-1,0:mmtv*mmv-1),nntv,mmtv),min(t),max(t), $
           min(f),max(f), $
           '','Frequency (MHz)','Gain/Background: '+date+' Beam:'+j_string+' P:'+k_string, $
           0,1,0,.05,.95,0,'.'
         ;!P.Multi = 0
         ; exit_ps
         
         for jj=0L, nf-1 do begin
           y         = x(*,jj)
           lin_data  = linfit(t,y)
           x(*,jj) = x(*,jj)/(lin_data[1]*t + lin_data[0])
         endfor 
         
             
         SPDYNPS,rebin(x(0:nntv*nnv-1,0:mmtv*mmv-1),nntv,mmtv),min(t),max(t), $
               min(f),max(f), $
               'Time [s]','Frequency (MHz)','Gain/Background/slope: '+date+' Beam:'+j_string+' P:'+k_string, $
               0,1,0,.05,.95,0,'.'
   Endif       
     ;****************
     ;update time loop
     ;****************
     tmin = tmin + N_time     
    end ;time
   end  ;beam
  end   ;pol
end ; dates 
set_plot,'PS'
device,/close

print,'The End of Program: Goodbye'
print,'SYSTIME=',SYSTIME(1) - st, ' sec'
journal
end         
  