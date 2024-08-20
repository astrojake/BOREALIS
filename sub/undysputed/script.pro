!path= expand_path('+/cep/lofar/nenufar')+'/:'+!path

READ_NU_SPEC,/help

;--- Pulsar

file='/databf2/nenufar-tf/20190212_123000_20190212_123700_B1919+21_TRACKING_BHR/B1919+21_B1_TRACKING_BHR_20190212_123033_0.spectra'
READ_NU_SPEC,file,/info

READ_NU_SPEC, file, data,time,freq,beam,ndata,nt,dt,nf,df,ns, fmin=44,fmax=54,nstokes=2,/fflat
help
x=reform(data(*,*,0))
display_data,x>0.9<1.3,1800
xd=de_disp(x,freq,12.4405d0,dt,max(freq))
display_data,xd>0.9<1.3,1800
xdt=reform(rebin(xd,nt,1))
p=1.337261142d0
p/dt	;       255.06232109069825
plot,time,xdt-smooth(xdt,255),/xsty
 
READ_NU_SPEC, file, data,time,freq,beam,ndata,nt,dt,nf,df,ns, tmin=30,tmax=60,fmin=44,fmax=54,nstokes=2,/fflat
help
x=reform(data(*,*,0))
display_data,x>0.9<1.3,1800

READ_NU_SPEC, file, data,time,freq,beam,ndata,nt,dt,nf,df,ns, tmin=30,tmax=60,fmin=44,fmax=70,nstokes=2,/fflat
x=reform(data(*,*,0))
display_data,x>0.9<1.3,1800

READ_NU_SPEC, file, data,time,freq,beam,ndata,nt,dt,nf,df,ns,fmin=44,fmax=70,nstokes=2,/fflat
x=reform(data(*,*,0))
display_data,x>0.9<1.3,1800
x=reform(data(*,*,1))
display_data,x>0.9<1.3,1800

READ_NU_SPEC, file, data,time,freq,beam,ndata,nt,dt,nf,df,ns,/fflat
help,data
display_data,data>0.9<1.3,1800

READ_NU_SPEC, file, data,time,freq,beam,ndata,nt,dt,nf,df,ns,/fflat,nstokes=2
x=reform(data(*,*,0))
help,x
display_data,x>0.9<1.3,1800


;--- Jupiter

file='/databf2/nenufar-tf/20190308_022000_20190308_062100_JUPITER_TRACKING_BHR/JUPITER_TRACKING_BHR_20190308_022036_0.spectra'
READ_NU_SPEC,file,/info

READ_NU_SPEC, file, data,time,freq,beam,ndata,nt,dt,nf,df,ns, tmin=5000, tmax=5300, ntimes=4
help
display_data,10.*alog10(data),1800

; nstokes=2 => Faraday fringes
; full time resol => S-bursts (7000-9000 sur 14300 spectres de 20 ms)
; reduced frequency band
; flat ?


