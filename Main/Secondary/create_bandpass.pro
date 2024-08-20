;Find bandpass of Nenufar
;

set_plot,'PS'
filename_ps = '20190323_TauBootis_Bandpass'           ;file name creation
file = '/databf2/nenufar-tf/20190323_013900_20190323_045300_TAU_BOOTIS_TRACKING_BHR/TAU_BOOTIS_TRACKING_BHR_20190323_013935_1.spectra'
device,filename=filename_ps+'.ps',/landscape
!p.multi=[0,1,2]

read_nenufar, file, data,time,freq,beam,ndata, nt,dt,nf,df,ns, tmin=0,tmax=100,/VERBOSE,$
              CHAN=CHAN,BEAMLETS=BEAMLETS,MINFREQ=MINFREQ,MAXFREQ=MAXFREQ

reduce_array,data,[nt,1],data_f
LE_AUTO_S,data_f,101,4.5,0,data_f,pnet
cgplot,freq,data_f,/xstyle,/ystyle;,xrange=[36,46]

bandpass = dblarr(chan, beamlets)

for i=0,beamlets-2 do begin
  bandpass(*,i) = data_f[i*(CHAN):(i+1)*(CHAN)-1]
  bandpass(*,i) =  bandpass(*,i)/mean([bandpass(20:60,i),bandpass(70:100,i)])
endfor

readcol,'/cep/lofar/nenufar/nenufar-tf/LOFAR-Bandpass-Correction-nfft-128.dat',coeffs

avg = dblarr(chan)
for i=0,CHAN-1 do begin
  avg[i] = median(bandpass[i,70:beamlets-1])
endfor
cgplot,avg,/ystyle, /xstyle,xtitle='Channel',ytitle='Intensity',thick=2;,xrange=[min(avg),max(avg)]
cgplot,bandpass[*,70],/overplot,color='green'
cgplot,bandpass[*,100],/overplot,color='blue'
cgplot,bandpass[*,190],/overplot,color='Red'
;cgplot,coeffs^2.0d,/overplot,color='red'
al_legend,['Median',strtrim(string(cgnumber_formatter(freq[60*128])))+'MHz',$
  strtrim(string(cgnumber_formatter(freq[100*128])))+'MHz',strtrim(string(cgnumber_formatter(freq[128*190])))+'MHz'],/bottom,/center,textcolor=['Black','Green','Blue','Red']

ratio_avg = 1.0d/avg

cgplot,ratio_avg,/ystyle, /xstyle,xtitle='Channel',ytitle='Intensity',thick=2,yrange=[0.85,2.0]
cgplot,coeffs^2.0d,/overplot,color='red'
al_legend,['1/Median','Cedric Coeffs^2'],textcolor=['black','red'],/center,/top

cgplot,ratio_avg/(coeffs^2.0d),/ystyle, /xstyle,xtitle='Channel',ytitle='Ratio',title='Ratio between Cedric & Emperical Bandpass'

openw,1,'LOFAR_Bandpass_data-128.dat'
printf,1,transpose(ratio_avg)
close,1

set_plot,'PS'
device,/close
cgfixps,filename_ps+'.ps'
spawn,'ps2pdf '+filename_ps+'.ps'
spawn,'rm '+filename_ps+'.ps'




end