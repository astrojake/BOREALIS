;--------------------------------------------------------------------------------------------------------------------------
  pro READ_NU_SPEC, file, data,time,freq,beam,ndata, nt,dt,nf,df,ns, tmin=tmin,tmax=tmax, fmin=fmin,fmax=fmax, $
			beams=beams, nchannels=nchannels, dm=dm, ntimes=ntimes, nstokes=nstokes, $
			fflat=fflat, fclean=fclean, tclean=tclean, info=info, help=help,$
			CHAN=CHAN,BEAMLETS=BEAMLETS,VERBOSE=VERBOSE,MINFREQ=MINFREQ,MAXFREQ=MAXFREQ,NOF_SAMPLES=NOF_SAMPLES,ALLTIME=ALLTIME
;--------------------------------------------------------------------------------------------------------------------------
; reads NenuFAR UnDySPuTeD data
; [INPUTS]
; from spectra file
; between tmin,tmax (s since file start) and fmin,fmax (MHz)  [default = all times & frequencies]
; within an optional list of beams
; integrating over nchannels=2^n per beamlet  [default = fftlen]
; after dedispersion at dm (pc.cm-3)  [default = 0, may be negative for dedispersing positive drifts]
; and over ntimes  [default = 1]
; nstokes=1 (I), 2 (Ix, Iy), 4 (full Stokes)
; if /info, only returns file structure description
; if /fflat then division by average spectrum
; fclean = Nsigmas for cleaning versus frequency within beamlets
; tclean = Nsigmas for cleaning versus time within beamlets
; [OUTPUT]
; data(t,f), time & freq ramps, beam # per channel
; ndata(t,f) = number of elementary data integrated per data element (after clean)
; nt,dt,nf,df,ns associated dimensions and resolutions
; [TBD]
; implement filling ndata, different fflat methods (background, median...), fclean, tclean
; weight by spectral response within beamlet depending on N Channels / Beamlet
; calculate and calibrate Stokes parameters from raw X & Y auto/cross-correlations

if keyword_set(help) then begin
  print
  print,'READ_NU_SPEC, file, data,time,freq,beam,ndata, nt,dt,nf,df,ns, tmin=tmin,tmax=tmax, fmin=fmin,fmax=fmax, $'
  print,'		beams=beams, nchannels=nchannels, dm=dm, ntimes=ntimes, nstokes=nstokes, $'
  print,'		fflat=fflat, fclean=fclean, tclean=tclean, info=info, help=help'
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
NOF_SAMPLES=nt
ALLTIME = time
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
MINFREQ = min(freq)
MAXFREQ = max(freq)

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

if not(keyword_set(nstokes)) then nstokes=1
if not(keyword_set(nchannels)) then nchannels=h0.fftlen
if nchannels gt h0.fftlen then nchannels=h0.fftlen
if not(keyword_set(ntimes)) then ntimes=1
if not(keyword_set(dm)) then dm=0.

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

if not(keyword_set(tmin)) then tmin=min(time)
if not(keyword_set(tmax)) then tmax=max(time)
wtime=where(time ge tmin and time le tmax)
if wtime(0) eq -1 then stop,'No time step selected'
fmin0=long(fmin/df0MHz)*df0MHz
textension=DT_DISP(fmin0,fmin0+df0MHz,dm)
blockmin=long((tmin-textension)/dtblock) > 0
blockmax=ceil((tmax+textension)/dtblock) < (nblocks-1)
blockmax=round(blockmax)
nblocksread=blockmax-blockmin+1
ntread=h0.nffte*nblocksread
time=time(h0.nffte*blockmin:h0.nffte*(blockmax+1)-1)

data=fltarr(nstokes,nchannels*nwbeamlets,long(ntread/ntimes))
; ndata=fltarr(nstokes,nchannels*nwbeamlets,long(ntread/ntimes))

h2={ilane:0L,ibeam:0L,ibeamlet:0L,pol0:fltarr(2,h0.fftlen,h0.nffte),pol1:fltarr(2,h0.fftlen,h0.nffte)}
h3=replicate(h2,nblocksread)

st=systime(1)
for i=0L,nwbeamlets-1L do begin
  if (i mod 10) eq 0 then print,i,' / ',nwbeamlets,'   ',systime(1)-st,' sec'
  beamlet=fltarr(nstokes,h0.fftlen,ntread)
; nbeamlet=beamlet+1			; used for filling ndata
  for j=blockmin,blockmax do begin
    point_lun,u,j*long64(blocksize)+8L*3+4L*6+wbeamlets(i)*beamletsize
    readu,u,h2
    h3(j-blockmin)=h2
  endfor
  if nstokes eq 1 then begin
    pol0=reform(h3.pol0,2,h0.fftlen,ntread)
;   beamlet(0,0:h0.fftlen/2+1,*) = pol0(0,h0.fftlen/2-2:*,*)+pol0(1,h0.fftlen/2-2:*,*)
;   beamlet(0,h0.fftlen/2-1:*,*) = pol0(0,0:h0.fftlen/2,*)+pol0(1,0:h0.fftlen/2,*)
    beamlet(0,0:h0.fftlen/2-1,*) = pol0(0,h0.fftlen/2:*,*)+pol0(1,h0.fftlen/2:*,*)		; seems to be the right
    beamlet(0,h0.fftlen/2:*,*)   = pol0(0,0:h0.fftlen/2-1,*)+pol0(1,0:h0.fftlen/2-1,*)	; distribution of negative & positive freq
;   beamlet(0,0:h0.fftlen/2,*)   = pol0(0,h0.fftlen/2-1:*,*)+pol0(1,h0.fftlen/2-1:*,*)
;   beamlet(0,h0.fftlen/2+1:*,*) = pol0(0,0:h0.fftlen/2-2,*)+pol0(0,0:h0.fftlen/2-2,*)
; old before inversion
;   beamlet=beamlet+h3.pol0(0,*,*,*)+h3.pol0(1,*,*,*)	; old
  endif
  if nstokes eq 2 then begin
    pol0=reform(h3.pol0,2,h0.fftlen,ntread)
    beamlet(*,0:h0.fftlen/2-1,*) = pol0(*,h0.fftlen/2:*,*)
    beamlet(*,h0.fftlen/2:*,*)   = pol0(*,0:h0.fftlen/2-1,*)
; old before inversion
;   beamlet=beamlet+h3.pol0
  endif
  if nstokes eq 4 then begin
    pol0=reform(h3.pol0,2,h0.fftlen,ntread)
    beamlet(0:1,0:h0.fftlen/2-1,*) = pol0(*,h0.fftlen/2:*,*)
    beamlet(0:1,h0.fftlen/2:*,*)   = pol0(*,0:h0.fftlen/2-1,*)
    pol1=reform(h3.pol1,2,h0.fftlen,ntread)
    beamlet(2:3,0:h0.fftlen/2-1,*) = pol1(*,h0.fftlen/2:*,*)
    beamlet(2:3,h0.fftlen/2:*,*)   = pol1(*,0:h0.fftlen/2-1,*)
; old before inversion
;   beamlet(0:1,*,*)=beamlet(0:1,*,*)+h3.pol0
;   beamlet(2:3,*,*)=beamlet(2:3,*,*)+h3.pol1
  endif
  if dm ne 0. then begin
    xf=freq(wbeamlets(i))+dindgen(h0.fftlen)/h0.fftlen*df0MHz
    for k=0,nstokes-1 do beamlet(k,*,*)=transpose(DE_DISP(transpose(reform(beamlet(k,*,*))),xf,dm,dt,freq(wbeamlets(i))+df0MHz))
  endif
  if keyword_set(fflat) then beamlet=beamlet/rebin(rebin(beamlet,nstokes,h0.fftlen,1),nstokes,h0.fftlen,ntread)
  if nchannels ne h0.fftlen then begin
    beamlet=rebin(beamlet,nstokes,nchannels,ntread)
;   nbeamlet=rebin(nbeamlet,nstokes,nchannels,ntread)*h0.fftlen/nchannels
  endif
  data(*,i*nchannels:(i+1)*nchannels-1,*)=rebin(beamlet(*,*,0:long(ntread/ntimes)*ntimes-1),nstokes,nchannels,long(ntread/ntimes))
; ndata(*,i*nchannels:(i+1)*nchannels-1,*)=rebin(nbeamlet(*,*,0:long(ntread/ntimes)*ntimes-1),nstokes,nchannels,long(ntread/ntimes))*ntimes
endfor

time=rebin(time(0:long(ntread/ntimes)*ntimes-1),long(ntread/ntimes))
wtime=where(time ge tmin and time le tmax)
time=time(wtime)
nt=n_elements(wtime)
dt=dt*ntimes

data=transpose(reform(data(*,*,wtime)))
; ndata=transpose(reform(ndata(*,*,wtime)))

freq=reform(rebin(reform(freq(wbeamlets),1,nwbeamlets),nchannels,nwbeamlets)+rebin(reform(dindgen(nchannels)/nchannels,nchannels,1),nchannels,nwbeamlets)*df0MHz,nchannels*nwbeamlets)
beam=reform(rebin(reform(ibeam(wbeamlets),1,nwbeamlets),nchannels,nwbeamlets),nchannels*nwbeamlets)
nf=nchannels*nwbeamlets
df=df0MHz/nchannels

BEAMLETS= h0.nbeamlets
CHAN=h0.fftlen
ns=nstokes
free_lun,u
end
