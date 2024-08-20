;***************************************************************
;      Read Lofar/Nenufar data and makes a quicklook
;**************************************************************
;Uses: hdf5 programs
;***************************************************************
;Author: Jake Turner (UVA)
;        J.M. Griessmeier (CNRS) and Philippe Zarka (Obspm)
;        
;Inputs:    planet      ;Planet Name
;           file_root   ;Root directory where the data is located (e.g. /databf/exoplanet/)
;           dates     
;           fmin
;           fmax,
;           tmin
;           N_time
;           steps
;           bmin  (starting beam)
;           bmax
;           Smin (starting Polarization)
;           Smax 
;         
;Required FLAGS 
;           LOFAR        ; Plots LOFAR Data 
;           NENUFAR      ; Plots NENUFAR Data         
;  
;Required NENUFAR Flags (only if NENUFAR is set)
;           NENUFAR_Name ; Filename of NENUFAR file
;          
;Keywords: 
;          PLOT_RAW       ; plot raw data (0=No, 1=YES)
;          DATA_NORM      ; Normalize data by 10% quantile (0=No, 1=YES)
;          DATA_SLOPE     ; Normalize and divide out long term slope (0=No, 1=YES)         
;          MASK           ; Apply MASK; input is the full mask array (same size as the full dataset)
;          FULL_QUICKLOOK ; plot all the data in one plot (0=NO, 1=YES) 
;          PAPERPLOTS     ; Removes header information (useful for test but not for papers) from plots (0=NO,1=YES)
;          time_UT        ; PLot the time axis in UT time
;          
;Note: Problem with numbers around 34379 (not zero) and 34380 (always zero); don't read full wavelength range (a HDF file problem)
;
pro lofar_quicklook,planet,file_root,date,fmin,fmax,tmin,N_time,steps,bmin,bmax,Smin,Smax,$
     LOFAR=LOFAR,NENUFAR=NENUFAR,DYNSPECBF=DYNSPECBF,filename=filename,$
     PLOT_RAW=PLOT_RAW,DATA_NORM=DATA_NORM,DATA_SLOPE=DATA_SLOPE,$
     MASK=MASK,FULL_QUICKLOOK=FULL_QUICKLOOK,PAPERPLOTS=PAPERPLOTS
             
st           = SYSTIME(1)
startmem_RFI =  MEMORY(/CURRENT)
 
Print, 'Filename',filename 
 
;Check Flags 
If keyword_set(LOFAR) then print,'***LOFAR Data ****'     
If keyword_set(NENUFAR) then begin
  print,'***NenuFAR Data****'
  If not(keyword_set(NENUFAR_Name)) then print,'ERROR; no filename set for NENUFAR data'
Endif   
If keyword_set(DYNSPECBF) then print,'***DYNSPEC Data ****'  
If not(keyword_set(FULL_QUICKLOOK)) then  FULL_QUICKLOOK= 1 ;plot all the data in one plot
If keyword_set(DYNSPECBF) then FULL_QUICKLOOK= 0
If steps lt 2 then FULL_QUICKLOOK=0.0  ;if only one stpe then no need for full quicklook 

If not(keyword_set(time_UT)) then time_UT  = 0         ; Plot time in UT 
time_sec = 1

If keyword_set(mask) then begin
   restore,filename=MASK
   BITARRAY_TO_BYTEARRAY, p2_full_bit, xsize, p2_full
Endif

If keyword_set(FULL_QUICKLOOK) then begin
  print, 'Full Quicklook'
  data_full_p  = ptrarr(steps)
  xt_rebin     =dblarr(steps)
Endif

;---------
;  Time 
;---------
start_time = tmin       ;save start_Time
total_t = N_time*steps  ;total time requested

;**************************************
;Don't need to change anything below here
;****************************************
  ;----------------------------------------------------------
  ;        Start loops for Beam and Polarization
  ;----------------------------------------------------------
  ;loops  beams, and polarization
  for jb=bmin,bmax do begin  ; loop for Beam: B000-
    for k=Smin,Smax do begin  ; loop for Polarization:S0-S3
      j_string = strtrim(string(uint(jb)),1)   ;Beam
      k_string2 = strtrim(string(uint(k)),1)   ;Polarization
      If k lt 4 then kin = k   ; I, Q, U, V
      If k eq 4 then kin = 3   ; V^2
      If k eq 5 then kin = 3   ; V'
      If k eq 6 then begin     ; L = sqrt(Q^2 + U^2)  
         kin1 = 1              ;Q
         kin  = kin1
         kin2 = 2              ;V
      endif 
      If k eq 0 then begin  ; I
        PolLabel = 'Stokes-I'       
      endif 
       If k eq 4 then begin  ; I
        PolLabel = 'Stokes-V'       
      endif 
      
      ;----------------------------------------------------------
      ;        Read Data for specific Beam and Polarization
      ;----------------------------------------------------------
      set_plot,'PS'
      
      If keyword_set(LOFAR) then filename_ps = date+'_SAP000_B00'+j_string+'_S'+k_string2+'_QuickLook'        ;file name creation
      If keyword_set(NENUFAR) then filename_ps = date+'_'+planet+'_QuickLook_Beam'+j_string+'_Pol'+k_string2        ;file name creation
      If keyword_set(DYNSPECBF) then filename_ps = date+'_QuickLook_Beam'+j_string+'_Pol'+k_string2       ;file name creation
      device,filename=filename_ps+'.ps',/landscape
       !p.multi=[0,1,2]

      If k lt 4 then k_string = strtrim(string(uint(k)),1)   ; I, Q, U, V
      If k eq 4 then k_string = strtrim(string(uint(3)),1)   ; V'
      If k eq 5 then k_string = strtrim(string(uint(3)),1)   ; V'
      If k eq 6 then begin     ; L = sqrt(Q^2 + U^2)  
         k_string = strtrim(string(uint(1)),1)
      endif 
      
      file = file_root+date+'/'+date+'_SAP000_B00'+j_string+'_S'+k_string+'_P000_bf'
      If keyword_set(LOFAR) then print, 'File for LOFAR: ', file ELSE print,'File: ',filename
      tmin = start_time  ;reset tmin
       
       ;-------------------------------------------------------------------
       ;                           Header info
       ;-------------------------------------------------------------------
       If keyword_set(LOFAR) then begin
        READ_H5_DATA, file,0,1,fmin,fmax,jb,kin,xx,tt,ff,ntt,nff,/NSPECTRA,min_freq=min_freq,$
          max_freq = max_freq,SAMPLING_TIME=SAMPLING_TIME,channel_width=channel_width,/READ_HEADER
          ;read h5 data has problems if you read the entire
          ;range of wavelengths for an array larger than 700 seconds
          If fmin le min_freq then fmin = min_freq + 0.01
          If fmax ge max_freq then fmax = max_freq - 0.01
       Endif
       If keyword_set(NENUFAR) then begin
         in_j = jb ;beam
         If k eq 0 then Stokes = 0 ;Stokes-I
         If k eq 1 then Stokes = 1 ;Stokes-Q
         If k eq 2 then Stokes = 2 ;Stokes-U
         If k eq 3 then Stokes = 3 ;Stokes-V
         If k eq 4 then Stokes = 4 ;Stokes-V Prime
         If k eq 0 then nstokes = 1 ELSE nstokes = 4
         temp       = filename+'*_'+strtrim(string(in_j),1)+'*'
         file_input = findfile(temp)
         READ_NU_FITS, file_input, command,param,variab, nt,dt,nf,df,ns,jd0,h0,fref,$
           temporary(data3),xt,xf,beam,ndata,corrt,corrf,/nodata,/quiet
        ; read_nenufar, file_root+filename, temporary(data2),tt,ff,temporary(in_j),ndata, ntt,dt,nff,df,ns,Stokes,nstokes=nstokes,$
        ;               tmin=0,tmax=1,fmin=9,fmax=70,/VERBOSE,CHAN=CHAN,BEAMLETS=BEAMLETS,/REMOVE_BADF,/CORRECT_BANDPASS,/DATA_BANDPASS
         sampling_time = dt 
         channel_width = df
       Endif
       If keyword_set(DYNSPECBF) then begin
         If k eq 0 then POL = 1 ;Stokes-I 
         If k eq 1 then POL = 2 ;Stokes-Q
         If k eq 2 then POL = 3 ;Stokes-U
         If k eq 3 then POL = 4 ;Stokes-V
         If k eq 4 then POL = 4 ;Stokes-V^2 
         If k eq 5 then POL = 4 ;Stokes-V prime
         print,'File to read',file_root+filename
         read_dynspec_fits,file_root+filename+'_beam'+j_string,0,1,fmin,fmax,POL,temporary(data2),xt,xf,nt,nf,dt=dt,df=df,ATmax=Tmax_all
         print,'max time',tmax_all
       Endif
       
       If keyword_set(mask) then N_time = sampling_time*4000. ;flag
       
        for kkk=0L,steps-1 do begin
          tmax =  tmin + N_time  ;set end time
          print,'-----------------------------------------------------------------------'
          print,'reading slice ',kkk,' start time: ', tmin, ' end time: ', tmax
          print,'-----------------------------------------------------------------------'
          
          If keyword_set(LOFAR) then begin 
           If k le 4 then READ_H5_DATA, file,tmin,tmax,fmin,fmax,jb,kin,x,t,f,nt,nf,UT=UT   ;I, U, Q, V, V^2
           If k eq 4 then begin  ;V'
              READ_H5_DATA, file_root+date+'/'+date+'_SAP000_B000'+'_S0_P000_bf',tmin,tmax,fmin,fmax,jb,0,x_I,t,f,nt,nf,UT=UT  ; I 
              READ_H5_DATA, file,tmin,tmax,fmin,fmax,jb,kin,x_V,t,f,nt,nf,UT=UT                                                ; V
           endif     
           If k eq 6 then begin ;L
              READ_H5_DATA, file,tmin,tmax,fmin,fmax,j,kin1,x_Q,t,f,nt,nf,UT=UT    ;Q
              READ_H5_DATA, file_root+date+'/'+date+'_SAP000_B000'+'_S2_P000_bf',tmin,tmax,fmin,fmax,jb,kin2,x_U,t,f,nt,nf,UT=UT  ; U 
             READ_H5_DATA, file_root+date+'/'+date+'_SAP000_B000'+'_S0_P000_bf',tmin,tmax,fmin,fmax,jb,0,x_I,t,f,nt,nf,UT=UT  ; I 
           endif 
          endif 
          
          If keyword_set(NENUFAR) then begin 
             beam_in = jb 
             If k eq 0 then Stokes = 0 ;Stokes-I
             If k eq 1 then Stokes = 1 ;Stokes-Q
             If k eq 2 then Stokes = 2 ;Stokes-U
             If k eq 3 then Stokes = 3 ;Stokes-V
             If k eq 4 then Stokes = 4 ;Stokes-V Prime
             If k eq 0 then nstokes = 1 ELSE nstokes = 4
             temp       = filename+'*_'+strtrim(string(beam_in),1)+'*'
             file_input = findfile(temp)
             READ_NU_FITS, file_input, command,param,variab, nt,dt,nf,df,ns,jd0,h0,fref,$
               x,t,f,beam,ndata,corrt,corrf,/quiet
            ; read_nenufar, file_root+filename, x,t,f,temporary(beam_in),ndata, nt,dt,nf,df,ns, Stokes,nstokes=nstokes,tmin=tmin,tmax=tmax,fmin=fmin,fmax=fmax,/REMOVE_BADF,$
            ;  /CORRECT_BANDPASS,/DATA_BANDPASS
          endif
          
          If keyword_set(DYNSPECBF) then begin
            If k eq 0 then POL = 1 ;Stokes-I 
            If k eq 1 then POL = 2 ;Stokes-Q
            If k eq 2 then POL = 3 ;Stokes-U
            If k eq 3 then POL = 4 ;Stokes-V
            If k eq 4 then POL = 4 ;Stokes-V^2 
            If k eq 5 then POL = 4 ;Stokes-V prime
            read_dynspec_fits,file_root+filename+'_beam'+j_string,tmin,tmax,fmin,fmax,POL,x,t,f,nt,nf
          Endif
                
           If kkk eq 0 then  freq_response = dblarr(steps,nf) 
           if keyword_set(Mask) then begin
             p2 = p2_full[kkk*nt:(kkk+1)*nt-1,0:nf-1]
           endif
         
          ;Reduce size
          nnv=nt/3000.d        ; for visualization in ps file
          if nnv ne long(nnv) then nnv=long(nnv)+1 else nnv=long(nnv)
          nntv=long(nt/nnv)
          mmv = nf/2000.d
          if mmv ne long(mmv) then mmv=long(mmv)+1 else mmv=long(mmv)
          mmtv=long(nf/mmv)
          nwin=0. & xp=100. & yp=100.
         
           If not(keyword_set(DATA_NORM)) then !p.multi=[0,1,2]
           If k eq 0 then k_string3 = 'I'
           If k eq 1 then k_string3 = 'Q'
           If k eq 2 then k_string3 = 'U'
           If k eq 3 then k_string3 = 'V'
           If k eq 4 then k_string3 = 'VPrime'
           If k eq 5 then k_string3 = 'Vprime'
           If k eq 6 then k_string3 = 'L'
           If not(keyword_set(PAPERPLOTS)) then label = planet+' Raw: '+date+' Beam:'+j_string+' P:'+k_string3 else label = planet
          
          ;Update X 
          If k le 3 then x = x        ;|V|             
          If k eq 4 then begin 
            y = x_V/X_I 
            MAKE_BACKGROUND,y,'', data_mean,ss,nn
            x = y - rebin(reform(data_mean,1,nf),nt,nf)                ;V^2
            ;x = x^2.0d
          endif 
          If k eq 5 then begin
            y = x_V/X_I 
            x = y - mean(y)                ;V'
          endif  
          If k eq 6 then begin
            x = sqrt(x_Q^2.0d + X_U^2.0d)  ;L
          endif 
          
          ;------------
          ;Plot Raw 
          ;------------
          If keyword_set(PLOT_RAW) then begin 
             If keyword_set(time_sec) then standard_plots,x,x*0.,t,f,label,/ONLY_PLOT,/NO_ZOOM,xunit='sec',panel_ps=1
             If keyword_set(time_UT) then standard_plots,x,x*0.,UT,f,label,/ONLY_PLOT,/NO_ZOOM,xunit='hour',xlabel='UT',panel_ps=1
          endif 
          If keyword_set(FULL_QUICKLOOK) then begin
             reduce_array,x,[10,1],x_r
             data_full_p[kkk]    = ptr_new(x_r)
            
             p2 = bytarr(nt,nf)+1b
             If keyword_set(LOFAR) then begin
               lofar_process_method2,x,p2,t,f,k,0.1,SAMPLING_TIME,channel_width,$
                 freq_response=frequ_response,/VERBOSE,RFI=1,/PS
                 freq_response(kkk,*)= frequ_response
             Endif
             If keyword_set(NENUFAR) then begin
               MAKE_BACKGROUND,x,'',data_avg,ss,nn 
               LE_AUTO_S,data_avg,101,4.5,0,frequ_response,pnet
               freq_response(kkk,*)= frequ_response   
             Endif
                          
             reduce_array,t,[10],t
             If keyword_set(time_UT) then reduce_array, UT,[10],UT
             If kkk eq 0 then t_full         = t
             If kkk gt 0 then t_full         = [t_full,t]
             If keyword_set(time_UT) then If kkk eq 0 then UT_full        = UT
             If keyword_set(time_UT) then If kkk gt 0 then UT_full        = [UT_full,UT]
           
             nt = n_elements(t) 
             reduce_array,t,[nt],temp_t
             xt_rebin[kkk] = temporary(temp_t)
          Endif ;full data

          If keyword_set(DATA_NORM) then begin
                print,'-------------------------'
                print, 'Start Normalize Data'
                print,'-------------------------'
                ;-------------------
                ; Normalize the Data
                ;-------------------
;                data_back = fltarr(nf)
;                upper_quantile = 0.1d  
;                make_background,x,'', data_back,ff,nn
;                LE_AUTO_S,data_back,101,4.5,0,data_back,pnet
;                data_back = smooth(data_back,5,/EDGE_TRUNCATE)
;                cgplot,f,data_back,title='Smoothed'
;             
;                If k ne 4 then x = x/rebin(reform(data_back,1,nf),nt,nf)  
;                If k eq 4 then x = x - rebin(reform(data_back,1,nf),nt,nf) 
  
                label = planet+' Normalized: '+date+' Beam:'+j_string+' P:'+k_string
                If not(keyword_set(PAPERPLOTS)) then label = planet+' Normalized: '+date+' Beam:'+j_string+' P:'+k_string3
                If keyword_set(PAPERPLOTS) then label = ''
                If keyword_set(time_sec) then standard_plots,x,x*0.,t,f,label,/ONLY_PLOT,/NO_ZOOM,xunit='sec',panel_ps=1,backgr=1
                If keyword_set(time_UT) then standard_plots,x,x*0.,UT,f,label,/ONLY_PLOT,/NO_ZOOM,xunit='hour',xlabel='UT',panel_ps=1,backgr=1
   
                 If keyword_set(MASK) then begin
                   x = x*p2
                 Endif

                 If keyword_Set(DATA_SLOPE) then begin 
                     !p.multi=[0,1,2]
                     reduce_array,x,[1,nf],y
                     lin_data  = linfit(t,y)
                     for jj=0L, nf-1 do begin
                       x(*,jj) = x(*,jj)/(lin_data[1]*t + lin_data[0])
                     endfor                      
                       
                    label = planet+' Normalized & De-Trended: '+date+' Beam:'+j_string+' P:'+k_string3  
                  	If not(keyword_set(PAPERPLOTS)) then label = label
                  	If keyword_set(PAPERPLOTS) then label = ''
                    If keyword_set(time_sec) then standard_plots,x,x*0.,t,f,label,/ONLY_PLOT,/NO_ZOOM,xunit='sec',panel_ps=1
              		  If keyword_set(time_UT) then standard_plots,x,x*0.,UT,f,label,/ONLY_PLOT,/NO_ZOOM,xunit='hour',xlabel='UT',panel_ps=1
                 endif ;slope
           Endif  ;data norm   
           
           ;****************
           ;update time loop
           ;****************
           tmin = tmin + N_time     
      end ;time
     
     If keyword_set(FULL_QUICKLOOK) then begin
       pointer_to_array,data_full_p,out_array=data_full
       
      ; If keyword_set(time_UT) then reduce_array,UT_full,[10],UT_full
       nt = n_elements(t_full)
       nf = n_elements(f) 
        
       n_steps = n_elements(freq_response[*,0])
       TimeFreqCorrect_array=dblarr(3,nf)

       for i=0,nf-1 do begin
        ;fit data
        yfit        = freq_response[*,i]
        fit_data    = poly_fit(xt_rebin,yfit,2)  ;fit
        TimeFreqCorrect_array(0,i) = fit_data[0]    ;a
        TimeFreqCorrect_array(1,i) = fit_data[1]    ;b*x
        TimeFreqCorrect_array(2,i) = fit_data[2]    ;c*x^2
       endfor          

       !P.Multi = 0
       !p.multi=[0,1,2]

       for ii=0,nf-1 do begin
          line           = TimeFreqCorrect_array[0,ii] + TimeFreqCorrect_array[1,ii]*t_full + TimeFreqCorrect_array[2,ii]*t_full^2.
          ;line           = TimeFreqCorrect_array[0,ii] + TimeFreqCorrect_array[1,ii]*t_full 
          If k ne 4 then data_full(*,ii)     = data_full(*,ii)/reform(line,nt,1)
          If k eq 4 then  data_full(*,ii)    = data_full(*,ii)-reform(line,nt,1)
       endfor 

       nnv=nt/3000.d        ; for visualization in ps file
       if nnv ne long(nnv) then nnv=long(nnv)+1 else nnv=long(nnv)
       nntv=long(nt/nnv)
       mmv = nf/2000.d
       if mmv ne long(mmv) then mmv=long(mmv)+1 else mmv=long(mmv)
       mmtv=long(nf/mmv)
        
       t_full = t_full/3600.0d ;hours
       xt_rebin = xt_rebin/3600.0d ;hours 
       label = planet+' Normalized (Full): '+date+' Beam:'+j_string+' P:'+k_string3
       If not(keyword_set(PAPERPLOTS)) then label = label
       If keyword_set(PAPERPLOTS) then label = planet+': '+Pollabel
	     If keyword_set(time_sec) then standard_plots,data_full,data_full*0.,t_full,f,label,/ONLY_PLOT,/NO_ZOOM,xunit='hours',panel_ps=1
	     If keyword_set(time_UT) then standard_plots,data_full,data_full*0.,UT_full,f,label,/ONLY_PLOT,/NO_ZOOM,xunit='hour',xlabel='UT',panel_ps=1
       HEAP_GC,/PTR,/VERBOSE 
     endif ;fullquickllook
     device,/close
     cgfixps,filename_ps+'.ps'
     cgps2pdf,filename_ps+'.ps',/delete
   end  ;beam
  end   ;pol
;set_plot,'PS'
print,'The End of Program: Goodbye'
PRINT, 'Memory required for Quicklook: ', (MEMORY(/HIGHWATER))/1d9 ,' Gb'
print,'SYSTIME=',SYSTIME(1) - st, ' sec'
end         
  
