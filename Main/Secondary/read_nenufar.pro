;Purpose: Read NenuFAR bf data for the pipeline 
;         - Removes bad frequencies
;
;USES: READ_NU_SPEC (from NenuFAR team and Phillippe)
;
pro read_nenufar, file, data,time,freq,beam,ndata, nt,dt,nf,df,ns, Stokes, tmin=tmin,tmax=tmax, fmin=fmin,fmax=fmax, $
      beams=beams, nchannels=nchannels, dm=dm, ntimes=ntimes, nstokes=nstokes,$
      fflat=fflat, fclean=fclean, tclean=tclean, info=info, help=help,$
      VERBOSE=VERBOSE,CHAN=CHAN,BEAMLETS=BEAMLETS,MINFREQ=MINFREQ,MAXFREQ=MAXFREQ,$
      NOF_SAMPLES=NOF_SAMPLES,ALLTIME=ALLTIME,$
      REMOVE_BADF=REMOVE_BADF,CORRECT_BANDPASS=CORRECT_BANDPASS,CEDRIC_BANDPASS=CEDRIC_BANDPASS,DATA_BANDPASS=DATA_BANDPASS,$
      RFI=RFI,mask=mask
 
 print,'Tmin', Tmin
 print, 'Tmax', Tmax
 
    ;------
    ;Flags
    ;------
   ; REMOVE_BADF      = 1   ;remove bad pixels (0=NO, 1=YES)
   ; CORRECT_BANDPASS = 1   ;Correct bandpass (0=NO, 1=YES)
   ; CEDRIC_BANDPASS  = 0   ;Use Cedric bandpass (0=NO, 1=YES)
   ; DATA_BANDPASS    = 1   ;Use empirical bandpass (0=NO, 1=YES)
   MINFREQ=1
   MAXFREQ=1
   If not(keyword_set(nstokes)) then nstokes = 1   
    If nstokes eq 1 then Stokes=0  ;nstokes=1 (I)
    If nstokes eq 2 then begin
      print,'You set nstokes =2, error! Set to Stokes-I'
      nstokes=1 ;    (I)
      Stokes=0
    Endif
    If nstokes eq 4 then begin
        print,'Nstokes=4, therefore read Stokes'
         If stokes eq 0 then istokes=0
         If stokes eq 1 then istokes=1
         If stokes eq 2 then istokes=2
         If stokes eq 3 then istokes=3
         If stokes eq 4 then istokes=3  ;Vprime 
          ;0=I, 1=U, Q=2, 3=V, 4=V'
    endif 
    file = file+'_'+strtrim(string(uint(beam)),1)+'.spectra
    
    ;Read Header info and full freq range
    ;nstokes=1 (I), 2 (Ix, Iy), 4 (full Stokes)
    READ_NU_SPEC, file,data,time,freq,beam,ndata, nt,dt,nf,df,ns,tmin=0,tmax=1,nstokes=nstokes,$
      CHAN=CHAN,BEAMLETS=BEAMLETS,MINFREQ=MINFREQ,MAXFREQ=MAXFREQ,VERBOSE=VERBOSE,NOF_SAMPLES=NOF_SAMPLES,ALLTIME=ALLTIME
    
;    pf = bytarr(nf)
;    pt = bytarr(nt)
  If not(keyword_set(RFI)) then begin
    READ_NU_SPEC2, file, data,time,freq,beam,ndata, nt,dt,nf,df,ns, tmin=tmin,tmax=tmax, $
      nchannels=CHAN, ntimes=ntimes, nstokes=nstokes
  endif 

   IF keyword_set(RFI) then begin
     print,'RFI in READ_NU_SPEC2'
     READ_NU_SPEC2, file, d,t,f,beam,ndata, nt,dt,nf,df,ns, tmin=tmin,tmax=tmax, $
       nchannels=CHAN, ntimes=ntimes, nstokes=nstokes,fclean=[5.5,51,0], tclean=[5.5,51,0]
     mask = bytarr(nt,nf) + 1b
     for i=0,nt-1 do begin 
       w_bad = where(d[i,*] eq 0,count)
       mask(i,w_bad) =  0
     endfor    
     STANDARD_PLOTS, d,p2,t,f,'TEst',xunit = 'secs',panel_ps=1
     STANDARD_PLOTS, mask,p2,t,f,'TEst',xunit = 'secs',panel_ps=1,/mask

   ENDIF
;   
    print,'N times', ntimes
    If nstokes eq 4 then begin
        data = data[*,*,istokes]   ;0=I, 1=U, Q=2, 3=V, 4=V'
    endif
    
    df = df*1d6  ; HZ 
   If keyword_set(VERBOSE) then begin
    print, 'Min Freq (MHz)', MINFREQ
    print, 'Max Freq (MHz)',MAXFREQ  
    print, 'Freq Res (Hz)', df
    print, 'Time Res (s)', dt 
   endif 
    
    If not(keyword_set(fmin)) then fmin = minfreq 
    If not(keyword_set(fmax)) then fmax = maxfreq
        
    ;-----------------------
    ;Correction of Bandpass
    ;-----------------------
    If CHAN eq 128 then begin
      If keyword_set(CORRECT_BANDPASS) then begin 
        print, '****** Correct BandPass *******'
        If keyword_set(CEDRIC_BANDPASS) then begin 
         readcol,'/cep/lofar/nenufar/nenufar-tf/LOFAR-Bandpass-Correction-nfft-128.dat',coeffs   ;Cedris bandpass
         coeffs = coeffs^2.0d 
        endif
        
        If keyword_set(DATA_BANDPASS) then begin
          print,' ***** Empirical bandpass correction applied ********'
          ;readcol,'/data/jake.turner/LOFAR_Bf_PipelineV4/Main/Secondary/LOFAR_Bandpass_data-128.dat',coeffs
          restore,filename='/data/jake.turner/LOFAR_Bf_PipelineV4/Main/Secondary/LOFAR_Bandpass_data-128.sav'
        Endif
        for i=0,beamlets-1 do begin
          data[*,i*(CHAN):(i+1)*(CHAN)-1] = data[*,i*(CHAN):(i+1)*(CHAN)-1]*transpose(rebin(coeffs,chan,nt))  
        endfor 
      endif ;correct bandpass   
    endif ;CHAN 128
    
    If CHAN eq 64 then begin
      If keyword_set(CORRECT_BANDPASS) then begin
        print, '****** Correct BandPass *******'
        If keyword_set(CEDRIC_BANDPASS) then begin
          readcol,'/cep/lofar/nenufar/nenufar-tf/LOFAR-Bandpass-Correction-nfft-064.dat',coeffs   ;Cedris bandpass
          ;restore,filename='/data/jake.turner/LOFAR_Bf_PipelineV4/Main/Secondary/LOFAR_Bandpass_data-064.dat'
          coeffs = coeffs^2.0d
        endif

        If keyword_set(DATA_BANDPASS) then begin
          print,' ***** Empirical bandpass correction applied ********'
          ;readcol,'/data/jake.turner/LOFAR_Bf_PipelineV4/Main/Secondary/LOFAR_Bandpass_data-64.dat',coeffs
          restore,filename='/data/jake.turner/LOFAR_Bf_PipelineV4/Main/Secondary/LOFAR_Bandpass_data-64.sav'
        Endif
        for i=0,beamlets-1 do begin
          data[*,i*(CHAN):(i+1)*(CHAN)-1] = data[*,i*(CHAN):(i+1)*(CHAN)-1]*transpose(rebin(coeffs,chan,nt))
        endfor
      endif ;correct bandpass
    endif ;CHAN 128
   
    If keyword_set(REMOVE_BADF) then remove_badf,data,freq,BEAMLETS*2.,CHAN/2.0,xFix=data,fFix=freq,XBAD=1
  
    ;update freqs 
    w = where(freq ge fmin and freq le fmax) 
    freq = freq(w) 
    nf = n_elements(freq) 
    data = data(*,w) 
    return 
    close, /all
end