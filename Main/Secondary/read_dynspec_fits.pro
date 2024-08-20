;************************************************
;Purpose: Read beamformed data in fits format 
;************************************************
;
;Note: This program was created for dynspec from the LoTSS survey using DynspecMS by Cyril. 
;      However, it can be use for any dynspec in fits format as long as the headers in the correct format
;      
;;Author: Jake Turner (UVA)
;        J.M. Griessmeier (CNRS) and Philippe Zarka (Obspm)
;        
;Inputs: filename  ;Filename to read
;        ts_min      ;Min time to output data (sec)
;        ts_max      ;Max time to output data (sec)
;        fs_min      ;Min freq to output data (MHz)
;        fs_max      ;Max freq to output data (MHz)
;        POL       ;Polarzation (1=I, 2=Q, 3=U, 4=V)        
;        
;Outputs: data (data)
;         xt (time)
;         xf (freq) 
;         nt (# of spectra)
;         nf (# of channels)
;
;Flag 
;      min_freq               ; Min Freq of the full observation
;      max_freq               ; Max Freq of the full observation
;      dt                     ;time resolution (sec)              
;      df                     ;freq resolution (Hz)
;      ATmax                  ;lmax time (s)
;      tall                   ;all time (s)
;      MJD                    ;all time in MJD (s)
pro read_dynspec_fits,filename,ts_min,ts_max,fs_min,fs_max,POL,x,xt,xf,nt,nf,$
                      fmin=fmin_all,fmax=fmax_all,dt=dt,df=df,$
                      tmax=tmax,tall=xt_all,nnt=nnt,nnf=nnf,MJD=MJD_all,HEADER_ONLY=HEADER_ONLY,$
                      PRINT_HEAD=PRINT_HEAD

x = readfits(filename+'.fits',header)  ; (time, freq)
If POL eq 1 then x = x[*,*,0]  ; I
If POL eq 2 then x = x[*,*,1] ;  Q 
If POL eq 3 then x = x[*,*,2] ;  U
If POL eq 4 then x = x[*,*,3] ;  V 

;******* Headers *****************
  ObJName  = SXPAR(header,'Name')           ;Object Name
  date     = SXPAR(header,'OBSID')          ;Obs ID number
  dt       = SXPAR(header,'CDELT1')         ;time resolution (sec)
  df       = SXPAR(header,'CDELT2')         ;freq resolution (MHz) 
  fmin_all = SXPAR(header,'FRQ-MIN')/1d6        ;min Freq in full data (MHz)
  fmax_all = SXPAR(header,'FRQ-MAX')/1d6        ;max Freq in full data (MHz)
  nnt      = SXPAR(header,'NAXIS1')         ;# of time samples in the full dataset
  nnf      = SXPAR(header,'NAXIS2')         ;# of freq chanels in the full dataset
    
  ;------
  ; Freq
  ;-------
  xf_all = dindgen(nnf,start=fmin_all,increment=df) ;freq of full dataset ; MHZ
  
  ;----
  ;Time
  ;----
  xt_all              = dindgen(nnt,start-0.0d,INCREMENT=dt)   ;time in seconds of full dataset 
  tmax                = max(xt_all)
  UT_start            = SXPAR(header,'OBS-STAR');Start time ;YYYY-MM-DD HH:MM:SS.SS  (ISO standard)
  UT_end              = SXPAR(header,'OBS-STOP')  
  MJD_start           = DATE_CONV(UT_start,'MODIFIED')   ;MJD
  MJD_end             = DATE_CONV(UT_end,'MODIFIED')   ;MJD
  MJD_dt              = dt/86400.d
  MJD_all             = dindgen(nnt,start=MJD_START,INCREMENT=MJD_dt)

If not(keyword_set(HEADER_ONLY)) then begin
  ;***************************************
  ;******** Selection of data  ***********
  ;***************************************
  v        = where(xf_all ge fs_min and xf_all le fs_max)
  nf       = n_elements(v)
  xf       = xf_all(v)
  x        = x[*,v]
  
  u       = where(xt_all ge ts_min and xt_all lt ts_max)
  nt      = n_elements(u)
  xt       = xt_all(u)
  mjd     = mjd_all(u)
  x       = x[u,*] 
endif 

If keyword_set(PRINT_HEAD) then begin
  print,'date', date 
  print,'df (MHz)',df 
  print,'dt (s)',dt
  print,'Full fmin (MHz)',fmin_all
  print,'Full fmax (MHz)',fmax_all
  print,'# f',nnf
  print,'# t',nnt
  print,'Max t (s)', tmax
  print,'MJD start',MJD_start
  print,'MJD end',MJD_end
Endif

end