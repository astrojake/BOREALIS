;---------------------------------------------------------------------------------------------------------------------
  pro READ_NU_FITS, file, command,param,variab, nt,dt,nf,df,ns,jd0,h0,fref, data,time,freq,beam,ndata,corrt,corrf, $
			nodata=nodata, nstokes=nstokes, quiet=quiet
;---------------------------------------------------------------------------------------------------------------------
; nstokes =0 [default] => read all, =1 => read I, =2 => read IV, =3 => read IVL

  if not(keyword_set(nstokes)) then nstokes=0

  command=string(MRDFITS(file,0,h))
  if not(keyword_set(quiet)) then print,command

  param=MRDFITS(file,1,h)
  if not(keyword_set(quiet)) then begin
    print
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
  endif

  variab=MRDFITS(file,2,h)
  nt=variab.nt & nf=variab.nf & ns=variab.ns
  dt=variab.dt & df=variab.df
  h0={fes:ulong64(0), timestamp:ulong64(0), blseqnum:ulong64(0), fftlen:0L, nfft2int:0L, fftovl:0L, apod:0L, nffte:0L, nbeamlets:0L}
  h0.fes=variab.fes+long64(2)^63 & h0.timestamp=variab.timestamp+long64(2)^63 & h0.blseqnum=variab.blseqnum+long64(2)^63
  h0.fftlen=variab.fftlen & h0.nfft2int=variab.nfft2int & h0.fftovl=variab.fftovl
  h0.apod=variab.apod & h0.nffte=variab.nffte & h0.nbeamlets=variab.nbeamlets
  jd0=variab.jd0 & fref=variab.fref
  if not(keyword_set(quiet)) then begin
    print
    print,'fes,timestamp,blseqnum:',h0.fes,h0.timestamp,h0.blseqnum
    print,'fftlen,nfft2int,fftovl,apod,nffte:',h0.fftlen,h0.nfft2int,h0.fftovl,h0.apod,h0.nffte
    print,'nbeamlets,filesize,beamletsize,blocksize,nblocks:',h0.nbeamlets,variab.filesize,variab.beamletsize,variab.blocksize,variab.nblocks
    print,'nt,nf,ns:',nt,nf,ns
    print,'dt,df:',dt,' ',variab.dtunit,df,' ',variab.dfunit
    print,'jd0,fref:',jd0,fref,' ',variab.frefunit
    print
  endif

  if not(keyword_set(nodata)) then begin
    datasize=long64(4)*nt*nf*ns
    datasizemax=long64(2)^31-1
    if datasize le datasizemax then begin
      data=MRDFITS(file,3,h)
      k=4
    endif else begin
      data=fltarr(nt,nf,ns)
      for k=0,ns-1 do begin
        x=MRDFITS(file,3+k,h)
        data(*,*,k)=x
      endfor
      k=k+3
    endelse
    if nstokes ne 0 then data=reform(data(*,*,0:nstokes-1))
    ndata=MRDFITS(file,k,h)
    time=MRDFITS(file,k+1,h)
    freq=MRDFITS(file,k+2,h)
    beam=MRDFITS(file,k+3,h)
    corrt=MRDFITS(file,k+4,h)
    corrf=MRDFITS(file,k+5,h)
    if not(keyword_set(quiet)) then help,file,command,param,variab,nt,dt,nf,df,ns,jd0,h0,fref,data,time,freq,beam,ndata,corrf,corrt
  endif else begin
    data=0. & time=0. & freq=0. & beam=0 & ndata=1 & corrt=1. & corrf=1.
    if not(keyword_set(quiet)) then help,file,command,param,variab,nt,dt,nf,df,ns,jd0,h0,fref
  endelse

end

