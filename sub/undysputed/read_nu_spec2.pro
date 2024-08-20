;--------------------------------------------------------------------------------------------------------------------------
  pro READ_NU_SPEC2, file, data,time,freq,beam,ndata, nt,dt,nf,df,ns, tmin=tmin,tmax=tmax, fmin=fmin,fmax=fmax, $
			beams=beams, nchannels=nchannels, ntimes=ntimes, nstokes=nstokes, dm=dm, $
			bandpass_corr=bandpass_corr, fflat=fflat, ex_chan=ex_chan, fclean=fclean, tclean=tclean, $
			data_compress=data_compress, round_times=round_times, beamlet_inc=beamlet_inc, info=info, help=help
;--------------------------------------------------------------------------------------------------------------------------
; READ_NU_SPEC2 reads NenuFAR UnDySPuTeD data, and performs on-the-fly integrations /t,f and pre-processing.
;
; [INPUT VARIABLE]
; file = filename TARGET_starting_date&time_n.spectra, with n=0,1,2,3 the output data lanes from UnDySPuTeD (with path /databf2/nenufar-tf/...)
;
; [INPUT KEYWORDS]
; tmin,tmax (in seconds since file start) [default = read all times]
; fmin,fmax (MHz) [default = read all frequencies]
; beams = list of beams to read from [default = read from all beams]
; nchannels = number of channels per beamlet after reading
; 	the recorded number of channels per beamlet is fftlen = 2^n with n=4 to 11, i.e. fftlen = 16 to 2048 channels
;	the number of channels per beamlet after reading should be 2^p with p < n, i.e. the data are rebinned on-the-fly
;	[default : p = n (2^n = fftlen)]
; ntimes = the rebin factor in time at reading (on-the-fly)
;	if the time resolution of recorded data is dt, read data will have a time
;	resolution ntimes*dt [default = 1 (full resolution)]
; nstokes = 1 (I), 2 (IV), 3 (IVL), 4 (IQUV) [default = 1]
; dm = dispersion measure (pc.cm-3) at which data within each integrated channel are de-dispersed, relative to the channel's center frequency [default = 0]
; /bandpass_corr = on-the-fly correction of the beamlet bandpass response (improves bandpass response but still imperfect)
; /fflat or fflat=1 = on-the-fly division of the dynamic spectrum by the average spectrum (if /fflat, then /bandpass_corr is useless)
;	fflat=2 = on-the-fly division of the dynamic spectrum by the median spectrum
; ex_chan = list of channel numbers to exclude from output dynamic spectrum (also when rebin on nchannels), e.g. [0,-1] to exclude beamlet edges, default = []
; fclean = [nsigmas,width,pm] for cleaning the data versus frequency (i.e. mitigate narrowband RFI)   within each beamlet, using subroutine LE_AUTO_
; tclean = [nsigmas,width,pm] for cleaning the data versus time      (i.e. mitigate broadband bursts) within each beamlet, using subroutine LE_AUTO_
;	nsigmas MUST be provided, pm or (width and pm) are optional and will be automatically computed if not provided
;	nsigmas is spikes detection threshold, width the sliding window size in kHz (fclean) or seconds (tclean), pm >1 / <1 / 0 = flag positive / negative spikes / both
; /data_compress = remove spectra and frequencies entirely masked by ex_chan, fclean, or tclean (implies irregular time or frequency ramp) [ default = none]
; round_times = rounds the number of spectra to the default exact multiple of round_times [default = none]
; beamlet_inc = beamlet increment for printing progress of reading/pre-processing [default = 10]
; /info = returns only file structure description on the screen
; /help = returns the syntax for launching READ_NU_SPEC2
;
; [OUTPUT VARIABLES]
; data(t,f,nstokes) = dynamic spectrum (linear values)
; time = time ramp since file start [seconds]
; freq = frequency ramp [MHz]
; beam = beam number affected to each channel
; ndata(t,f) = number of elementary data integrated per data element (after clean)
; nt,dt = number of time steps (dimension #1 of data) and time resolution [seconds] after reading
; nf,df = number of frequency steps (dimension #2 of data) and frequency resolution [MHz] after reading
; ns = number of Stokes (dimension #3 of data)


if keyword_set(help) then begin
  print
  print,'READ_NU_SPEC2, file, data,time,freq,beam,ndata, nt,dt,nf,df,ns, tmin=tmin,tmax=tmax, fmin=fmin,fmax=fmax, $'
  print,'		beams=beams, nchannels=nchannels, ntimes=ntimes, nstokes=nstokes, dm=dm, $'
  print,'		bandpass_corr=bandpass_corr, fflat=fflat, ex_chan=ex_chan, fclean=fclean, tclean=tclean, $'
  print,'		data_compress=data_compress, round_times=round_times, beamlet_inc=beamlet_inc, info=info, help=help'
  print
  return
endif

openr,u,file,/get_lun
f=fstat(u)

h0={fes:long64(0), timestamp:long64(0), blseqnum:long64(0), fftlen:0L, nfft2int:0L, fftovl:0L, apod:0L, nffte:0L, nbeamlets:0L}
readu,u,h0
print
print,'First Effective Sample         : ',h0.fes
print,'TIMESTAMP of first sample      : ',h0.timestamp
print,'BLOCKSEQNUMBER of first sample : ',h0.blseqnum
print,'fftlen = N Channels / Beamlet  : ',h0.fftlen
print,'nfft2int                       : ',h0.nfft2int
print,'FFToverlap 0/1                 : ',h0.fftovl
print,'apodisation                    : ',h0.apod 
print,'nffte spectra / block          : ',h0.nffte
print,'N Beamlets                     : ',h0.nbeamlets

filesize=f.size
beamletsize=4L*3+4L*2*h0.fftlen*h0.nffte*2
blocksize=8L*3+4L*6+beamletsize*h0.nbeamlets
nblocks=double(filesize)/blocksize
print,'filesize                      : ',filesize
print,'blocksize                     : ',blocksize
print,'beamletsize                   : ',beamletsize
print,'N blocks in file              : ',long(nblocks),'   +',nblocks-long(nblocks)

nblocks=long(nblocks)
df0=195312.5d0
df=df0/h0.fftlen
dt=1.d0/df*h0.nfft2int
nt=h0.nffte*nblocks
dtblock=dt*h0.nffte
time=dindgen(nt)*dt
nf=h0.nbeamlets*h0.fftlen
print,'N frequency channels          : ',nf
print,'df (Hz)                       : ',df
print,'N time steps                  : ',nt
print,'dt (s)                        : ',dt
print,'dt per block (s)              : ',dtblock
print,'Time (s) min,max              : ',min(time),max(time)

ilane=lonarr(h0.nbeamlets)
ibeam=lonarr(h0.nbeamlets)
ibeamlet=lonarr(h0.nbeamlets)
h1={ilane:0L,ibeam:0L,ibeamlet:0L}
for i=0L,h0.nbeamlets-1L do begin
  point_lun,u,8L*3+4L*6+i*beamletsize
  readu,u,h1
  ilane(i)=h1.ilane
  ibeam(i)=h1.ibeam
  ibeamlet(i)=h1.ibeamlet
endfor
df0MHz=df0/1.d6
freq=ibeamlet*df0MHz

slane=ilane(uniq(ilane,sort(ilane))) & nslane=n_elements(slane)
sbeam=ibeam(uniq(ibeam,sort(ibeam))) & nsbeam=n_elements(sbeam)

for i=0L,nslane-1 do begin
for j=0L,nsbeam-1 do begin
  w=where(ilane eq slane(i) and ibeam eq sbeam(j))
  if w(0) ne -1 then begin
    print,'Lane ',slane(i),'   Beam ',sbeam(j)
    print,'Beamlet min,max : ',min(ibeamlet(w)),max(ibeamlet(w))
    print,'Frequency (MHz) min,max : ',min(freq(w)),max(freq(w))+df0MHz
  endif
endfor
endfor

print
if keyword_set(info) then return

if not(keyword_set(tmin)) then tmin=min(time)
if not(keyword_set(tmax)) then tmax=max(time)
wtime=where(time ge tmin and time le tmax)
if wtime(0) eq -1 then stop,'No time step selected'

if not(keyword_set(fmin)) then fmin=min(freq)
if not(keyword_set(fmax)) then fmax=max(freq)
if not(keyword_set(beams)) then beams=[-1] else beams=[beams]
nbeams=n_elements(beams)
ok=bytarr(h0.nbeamlets)
for i=0L,h0.nbeamlets-1L do begin
  if freq(i)+df0MHz ge fmin and freq(i) le fmax then begin
    for j=0,nbeams-1 do if beams(j) eq -1 or beams(j) eq ibeam(i) then ok(i)=1b
  endif
endfor
wbeamlets=where(ok eq 1)
if wbeamlets(0) ne -1 then nwbeamlets=n_elements(wbeamlets) else stop,'No beamlet selected'

if not(keyword_set(nchannels)) then nchannels=h0.fftlen
if nchannels gt h0.fftlen then nchannels=h0.fftlen
if not(keyword_set(ntimes)) then ntimes=1
if not(keyword_set(nstokes)) then nstokes=1

if not(keyword_set(dm)) then dm=0.
fmin0=long(fmin/df0MHz)*df0MHz
textension=DT_DISP(fmin0,fmin0+df0MHz,dm)

blockmin=long((tmin-textension)/dtblock) > 0
blockmax=ceil((tmax+textension)/dtblock) < (nblocks-1)
blockmax=round(blockmax)
nblocksread=blockmax-blockmin+1
ntread=h0.nffte*nblocksread
time=time(h0.nffte*blockmin:h0.nffte*(blockmax+1)-1)
if keyword_set(bandpass_corr) then gain=rebin(reform(compute_PFB_corr(h0.fftlen),1,h0.fftlen,1),2,h0.fftlen,ntread)

beamlet_flag=bytarr(h0.fftlen)+1
if n_elements(ex_chan) gt 0 then for i=0,n_elements(ex_chan)-1 do beamlet_flag(ex_chan(i))=0
beamlet_flag=rebin(reform(beamlet_flag,h0.fftlen,1),h0.fftlen,ntread)

data =fltarr(nstokes,nchannels*nwbeamlets,long(ntread/ntimes))
ndata=intarr(nchannels*nwbeamlets,long(ntread/ntimes))

h2={ilane:0L,ibeam:0L,ibeamlet:0L,pol0:fltarr(2,h0.fftlen,h0.nffte),pol1:fltarr(2,h0.fftlen,h0.nffte)}
h3=replicate(h2,nblocksread)

st=systime(1)
if not(keyword_set(beamlet_inc)) then beamlet_inc=10

for i=0L,nwbeamlets-1L do begin
  if (i mod beamlet_inc) eq 0 then print,i,' / ',nwbeamlets,'   ',systime(1)-st,' sec'

  beamlet=fltarr(nstokes,h0.fftlen,ntread)
  for j=blockmin,blockmax do begin
    point_lun,u,j*long64(blocksize)+8L*3+4L*6+wbeamlets(i)*beamletsize
    readu,u,h2
    h3(j-blockmin)=h2
  endfor

  pol0=reform(h3.pol0,2,h0.fftlen,ntread)
  if keyword_set(bandpass_corr) then pol0=pol0*gain
  beamlet(0,0:h0.fftlen/2-1,*) = pol0(0,h0.fftlen/2:*,*)+pol0(1,h0.fftlen/2:*,*)		; I
  beamlet(0,h0.fftlen/2:*,*)   = pol0(0,0:h0.fftlen/2-1,*)+pol0(1,0:h0.fftlen/2-1,*)	; distribution of negative & positive frequencies
  if nstokes ge 2 then begin
    pol1=reform(h3.pol1,2,h0.fftlen,ntread)
    if keyword_set(bandpass_corr) then pol1=pol1*gain
    beamlet(1,0:h0.fftlen/2-1,*) = 2.*pol1(1,h0.fftlen/2:*,*)				; V
    beamlet(1,h0.fftlen/2:*,*)   = 2.*pol1(1,0:h0.fftlen/2-1,*)
  endif
  if nstokes eq 3 then begin								 	; L
    beamlet(2,0:h0.fftlen/2-1,*) = sqrt((pol0(0,h0.fftlen/2:*,*)-pol0(1,h0.fftlen/2:*,*))^2+(2.*pol1(0,h0.fftlen/2:*,*))^2)
    beamlet(2,h0.fftlen/2:*,*)   = sqrt((pol0(0,0:h0.fftlen/2-1,*)-pol0(1,0:h0.fftlen/2-1,*))^2+(2.*pol1(0,0:h0.fftlen/2-1,*))^2)
  endif
  if nstokes eq 4 then begin
    beamlet(2,0:h0.fftlen/2-1,*) = pol0(0,h0.fftlen/2:*,*)-pol0(1,h0.fftlen/2:*,*)		; Q
    beamlet(2,h0.fftlen/2:*,*)   = pol0(0,0:h0.fftlen/2-1,*)-pol0(1,0:h0.fftlen/2-1,*)	;
    beamlet(3,0:h0.fftlen/2-1,*) = 2.*pol1(0,h0.fftlen/2:*,*)				; U
    beamlet(3,h0.fftlen/2:*,*)   = 2.*pol1(0,0:h0.fftlen/2-1,*)
    beamlet = beamlet([0,2,3,1],*,*)
  endif

;  if nstokes eq 1 then begin
;    pol0=reform(h3.pol0,2,h0.fftlen,ntread)
;    if keyword_set(bandpass_corr) then pol0=pol0*gain
;    beamlet(0,0:h0.fftlen/2-1,*) = pol0(0,h0.fftlen/2:*,*)+pol0(1,h0.fftlen/2:*,*)		; distribution of negative
;    beamlet(0,h0.fftlen/2:*,*)   = pol0(0,0:h0.fftlen/2-1,*)+pol0(1,0:h0.fftlen/2-1,*)	; & positive frequencies
;  endif
;  if nstokes eq 2 then begin
;    pol0=reform(h3.pol0,2,h0.fftlen,ntread)
;    if keyword_set(bandpass_corr) then pol0=pol0*gain
;    beamlet(*,0:h0.fftlen/2-1,*) = pol0(*,h0.fftlen/2:*,*)
;    beamlet(*,h0.fftlen/2:*,*)   = pol0(*,0:h0.fftlen/2-1,*)
;  endif
;  if nstokes eq 4 then begin
;    pol0=reform(h3.pol0,2,h0.fftlen,ntread)
;    if keyword_set(bandpass_corr) then pol0=pol0*gain
;    beamlet(0:1,0:h0.fftlen/2-1,*) = pol0(*,h0.fftlen/2:*,*)
;    beamlet(0:1,h0.fftlen/2:*,*)   = pol0(*,0:h0.fftlen/2-1,*)
;    pol1=reform(h3.pol1,2,h0.fftlen,ntread)
;    if keyword_set(bandpass_corr) then pol1=pol1*gain
;    beamlet(2:3,0:h0.fftlen/2-1,*) = pol1(*,h0.fftlen/2:*,*)
;    beamlet(2:3,h0.fftlen/2:*,*)   = pol1(*,0:h0.fftlen/2-1,*)
;  endif

  if keyword_set(fflat) then begin
    if fflat eq 1 then beamlet=beamlet/rebin(rebin(beamlet,nstokes,h0.fftlen,1),nstokes,h0.fftlen,ntread)
    if fflat eq 2 then beamlet=beamlet/rebin(reform(median(beamlet,dimension=3),nstokes,h0.fftlen,1),nstokes,h0.fftlen,ntread)
  endif

  nbeamlet=bytarr(h0.fftlen,ntread)+1b
  if n_elements(ex_chan) gt 0 or keyword_set(fclean) or keyword_set(tclean) then begin
    if n_elements(ex_chan) gt 0 then nbeamlet=nbeamlet*beamlet_flag
    if keyword_set(fclean) then begin
      x=reform(rebin(beamlet(0,*,*),1,h0.fftlen,1))
      if n_elements(fclean) eq 3 then pm=fclean(2) else pm=0
      if n_elements(fclean) ge 2 then width=round(fclean(1)*1000./df) else width=round(50000./df)
      npixels=fix(h0.fftlen/8.+1) <10
      LE_AUTO,x,width,fclean(0),pm, xnet,para,EDGE=2,NPIX=npixels,MAXITERATIONS=5
      nbeamlet=nbeamlet*rebin(reform(para,h0.fftlen,1),h0.fftlen,ntread)
    endif
    if keyword_set(tclean) then begin
      x=reform(rebin(beamlet(0,*,*),1,1,ntread))
      if n_elements(tclean) eq 3 then pm=tclean(2) else pm=0
      if n_elements(tclean) ge 2 then width=round(tclean(1)/dt) else width=round(3./dt)
      LE_AUTO,x,width,tclean(0),pm, xnet,para,EDGE=2,MAXITERATIONS=5
      nbeamlet=nbeamlet*rebin(reform(para,1,ntread),h0.fftlen,ntread)
    endif
  endif
  for k=0,nstokes-1 do beamlet(k,*,*)=beamlet(k,*,*)*nbeamlet

  if dm ne 0. then begin			; de-dispersion within each beamlet
    xf=freq(wbeamlets(i))+dindgen(h0.fftlen)/h0.fftlen*df0MHz
    for k=0,nstokes-1 do beamlet(k,*,*)=transpose(DE_DISP(transpose(reform(beamlet(k,*,*))),xf,dm,dt,freq(wbeamlets(i))+df0MHz))
    nbeamlet=transpose(DE_DISP(transpose(nbeamlet),xf,dm,dt,freq(wbeamlets(i))+df0MHz))
  endif

 if nchannels ne h0.fftlen then begin		; spectral integration of each beamlet on nchannels
    beamlet=rebin(beamlet,nstokes,nchannels,ntread)
    nbeamlet=rebin(nbeamlet*h0.fftlen/nchannels,nchannels,ntread)
  endif

 if nchannels ne 1 and dm ne 0. then begin	; re-dispersion of beamlet channels
    xf=rebin(xf,nchannels)
    for k=0,nstokes-1 do beamlet(k,*,*)=transpose(DE_DISP(transpose(reform(beamlet(k,*,*))),xf,-dm,dt,freq(wbeamlets(i))+df0MHz))
    nbeamlet=transpose(DE_DISP(transpose(nbeamlet),xf,-dm,dt,freq(wbeamlets(i))+df0MHz))
  endif

  data(*,i*nchannels:(i+1)*nchannels-1,*)=rebin(beamlet(*,*,0:long(ntread/ntimes)*ntimes-1),nstokes,nchannels,long(ntread/ntimes))
  ndata(i*nchannels:(i+1)*nchannels-1,*)=rebin(nbeamlet(*,0:long(ntread/ntimes)*ntimes-1)*ntimes,nchannels,long(ntread/ntimes))
endfor

time=rebin(time(0:long(ntread/ntimes)*ntimes-1),long(ntread/ntimes))
wtime=where(time ge tmin and time le tmax)
time=time(wtime)
nt=n_elements(wtime)
dt=dt*ntimes

data=transpose(reform(data(*,*,wtime)))
ndata=transpose(reform(ndata(*,wtime)))

freq=reform(rebin(reform(freq(wbeamlets),1,nwbeamlets),nchannels,nwbeamlets)+rebin(reform(dindgen(nchannels)/nchannels,nchannels,1),nchannels,nwbeamlets)*df0MHz,nchannels*nwbeamlets)
beam=reform(rebin(reform(ibeam(wbeamlets),1,nwbeamlets),nchannels,nwbeamlets),nchannels*nwbeamlets)
nf=nchannels*nwbeamlets
df=df0MHz/nchannels
ns=nstokes

if keyword_set(data_compress) then begin
  w=where(reform(rebin(ndata,1,nf)) ne 0)
  data=data(*,w,*)
  ndata=ndata(*,w)
  freq=freq(w)
  beam=beam(w)
  nf=n_elements(freq)
  w=where(reform(rebin(ndata,nt,1)) ne 0)
  data=data(w,*,*)
  ndata=ndata(w,*)
  time=time(w)
  nt=n_elements(time)
endif

if keyword_set(round_times) then begin
  nt=long(nt/round_times)*round_times
  data=data(0:nt-1,*,*)
  ndata=ndata(0:nt-1,*)
  time=time(0:nt-1)
endif

close,u & free_lun,u
end
