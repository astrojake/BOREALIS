;---------------------------------------------------------
  pro CHECK_ES02_FITS, pathobs, files, nocorrt=nocorrt
;---------------------------------------------------------
; nocorrt before 2019/07/25

if strmid(pathobs,strlen(pathobs)-1,1) ne '/' then pathobs=pathobs+'/'
if n_elements(files) eq 0 then begin
  spawn,'ls -1 '+pathobs+'*.fits > list'
  files=['rien'] & buf=''
  on_ioerror,suite
  openr,v,'list',/get_lun
encore:
  readf,v,buf
  files=[files,DELPATH(buf)]
  goto,encore
suite:
  close,v & free_lun,v
  files=files[1:*]
endif

onoff=['ON','OFF1','OFF2','OFF3']

for ifile=0,n_elements(files)-1 do begin

  file=files(ifile)
  name=strmid(file,0,strlen(file)-13)
  set_ps,name+'.fits.ps',/portrait
  device,/color
  !p.charsize=1.3
  loadct,0
  ionoff=fix(strmid(file,strlen(file)-14,1))
  namefits=pathobs+file

  command=string(MRDFITS(namefits,0,h))
  print
  print,h,command

  param=MRDFITS(namefits,1,h)
  print
  print,h
  print,'version,versiondate: ',param.version,param.versiondate
  print,'tmin,tmax,julian:',param.tmin,param.tmax,param.julian
  print,'fmin,fmax,exactfreq:',param.fmin,param.fmax,param.exactfreq
  print,'beams:',param.beams
  print,'nchannels,ntimes,nobandpass:',param.nchannels,param.ntimes,param.nobandpass
  print,'ex_chan:',param.ex_chan
  print,'ex_beamlets:',param.ex_beamlets
  print,'fclean:',param.fclean
  print,'bclean:',param.bclean
  print,'tclean:',param.tclean
  print,'fflat,tflat:',param.fflat,' - ',param.tflat
  print,'dm:',param.dm
  print,'fcompress,tcompress,fill,round_times,round_freq:',param.fcompress,param.tcompress,param.fill,param.round_times,param.round_freq

  variab=MRDFITS(namefits,2,h)
  print
  print,h
  print,'fes,timestamp,blseqnum:',variab.fes,variab.timestamp,variab.blseqnum
  print,'fftlen,nfft2int,fftovl,apod,nffte:',variab.fftlen,variab.nfft2int,variab.fftovl,variab.apod,variab.nffte,variab.nbeamlets
  print,'nbeamlets,filesize,beamletsize,blocksize,nblocks:',variab.nbeamlets,variab.filesize,variab.beamletsize,variab.blocksize,variab.nblocks
  print,'jd0,dt:',variab.jd0,variab.dt,' ',variab.dtunit
  print,'df,fref:',variab.df,' ',variab.dfunit,variab.fref,' ',variab.frefunit
  nt=variab.nt & nf=variab.nf & ns=variab.ns
  print,'nt,nf,ns',nt,nf,ns
  print

  datasize=long64(4)*variab.nt*variab.nf*variab.ns
  datasizemax=long64(2)^31-1
  if datasize le datasizemax then begin
    data=MRDFITS(namefits,3,h)
    print,h
    k=4
  endif else begin
    data=fltarr(nt,nf,ns)
    for k=0,ns-1 do begin
      x=MRDFITS(namefits,3+k,h)
      print,h
      data(*,*,k)=x
    endfor
    k=k+3
  endelse
  help,data

  ndata=MRDFITS(namefits,k,h)
  print
  print,h
  help,ndata

  time=MRDFITS(namefits,k+1,h)
  print
  print,h
  help,time

  freq=MRDFITS(namefits,k+2,h)
  print
  print,h
  help,freq

  corrt=MRDFITS(namefits,k+4,h)
  print
  print,h
  help,corrt

  corrf=MRDFITS(namefits,k+5,h)
  print
  print,h
  help,corrf

  w=where(rebin(ndata,1,nf) eq 0, nw)
  if w(0) ne -1 then begin
    m=reform(median(data,dimension=2),nt,1,ns)
    for i=0,nw-1 do data(*,w(i),*)=m
  endif

  nplots=ceil(max(time)/3600.)		; number of plots (1 / hour)
  titre=file+', '+onoff(ionoff)+', proc'
  ;NU_PLOTS, data, freq, time, titre, ndata, /L, /FZ, TZ=[1,600]

  if not(keyword_set(nocorrt)) then begin
    loadct,13
    !p.multi=[0,1,4]
    s=size(corrt)
    plot,time/3600.,corrt(*,0),xtit='Time since file start (h)',ytit='Gain(t) correction : corrt(0)',tit=titre,/xsty,/ynoz
    if s(0) gt 1 then begin
      plot,time/3600.,corrt(*,1),xtit='Time since file start (h)',ytit='data freq.-integ. time profile & fit',tit=titre,/xsty,/ynoz
      oplot,time/3600.,corrt(*,2),color=250
      plot,time/3600.,corrt(*,3),xtit='Time since file start (h)',ytit='[a].log(t) + b',/xsty,/ynoz
      plot,time/3600.,corrt(*,4),xtit='Time since file start (h)',ytit='a.log(t) + [b]',/xsty,/ynoz
    endif
    w=where(time ge 1 and time le 600)
    plot,time(w)/60.,corrt(w,0),xtit='Time since file start (min)',ytit='Gain(t) correction : corrt(0)',tit=titre,/xsty,/ynoz
    if s(0) gt 1 then begin
      plot,time(w)/60.,corrt(w,1),xtit='Time since file start (min)',ytit='data freq.-integ. time profile & fit',tit=titre,/xsty,/ynoz
      oplot,time(w)/60.,corrt(w,2),color=250
    endif
    loadct,0
  endif

  !p.multi=[0,1,4]
  for iplot=0,nplots-1 do begin
    w=where(time ge iplot*3600. and time le (iplot+1)*3600., nw)
    if nw gt 10 then begin
      nw=long(nw/2)*2 & w=w(0:nw-1)
      SPDYNPS, rebin(10.*alog10(data(w,*,0)),nw/2,nf),min(time(w))/60,max(time(w))/60,min(freq),max(freq),'Time since file start (min)','Frequency (MHz)',titre+', Stokes I',0,1,-0.5,0.,0.98,0,'dB'
      SPDYNPS, rebin(data(w,*,1)/data(w,*,0),nw/2,nf),min(time(w))/60,max(time(w))/60,min(freq),max(freq),'Time since file start (min)','Frequency (MHz)',titre+', Stokes V/I',0,0,0,0.02,0.98,0,'V/I'
      SPDYNPS, rebin(data(w,*,2)/data(w,*,0),nw/2,nf),min(time(w))/60,max(time(w))/60,min(freq),max(freq),'Time since file start (min)','Frequency (MHz)',titre+', Stokes L/I',0,0,0,0.02,0.98,0,'L/I'
      SPDYNPS, rebin(ndata(w,*)*100,nw/2,nf),min(time(w))/60,max(time(w))/60,min(freq),max(freq),'Time since file start (min)','Frequency (MHz)',titre+', pixOK=1-flag',0,0,0,0.,1.,1,'%'
    endif
  endfor

  device,/close
  ps_pdf;,name+'.fits.ps',/rem

endfor

exit_ps
end
