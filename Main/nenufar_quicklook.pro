;pro nenufar_quicklook,planet,file_root,date,fmin,fmax,tmin,N_time,steps,bmin,bmax,Smin,Smax,PLOT_RAW=PLOT_RAW,$
;  DATA_NORM=DATA_NORM,DATA_SLOPE=DATA_SLOPE,ALL_FREQ=ALL_FREQ,$
;  MASK=MASK,FULL_QUICKLOOK=FULL_QUICKLOOK,PAPERPLOTS=PAPERPLOTS

;;Setup Plot
set_plot,'PS'
;filename_ps = '20190307_Jupiter_QuickLook_ALL'           ;file name creation
;file = '/databf2/nenufar-tf/20190307_022000_20190307_062100_JUPITER_TRACKING_BHR/JUPITER_TRACKING_BHR_20190307_022033_0.spectra'

;Date 1 
filename_ps = '20190323_TauBootis_QuickLook_CedricCorrect_Beam0'           ;file name creation
file = '/databf2/nenufar-tf/20190323_013900_20190323_045300_TAU_BOOTIS_TRACKING_BHR/TAU_BOOTIS_TRACKING_BHR_20190323_013935_0.spectra'

;Date 2
filename_ps = '20190328_TauBootis_QuickLook_CedricCorrect_64_Beam2'
file = '/databf2/nenufar-tf/20190328_001800_20190328_043400_TAU_BOOTIS_TRACKING_BHR/TAU_BOOTIS_TRACKING_BHR_20190328_001833'
beam = 2
 
 
device,filename=filename_ps+'.ps',/landscape
!p.multi=[0,1,2]


;file = '/databf2/nenufar-tf/20190308_022000_20190308_062100_JUPITER_TRACKING_BHR/JUPITER_TRACKING_BHR_20190308_022036_0.spectra'
;n_steps = 1
;for i=0,n_steps-1 do begin
  ;print,'I', i
;  READ_NU_SPEC, file,data,time,freq,beam,ndata, nt,dt,nf,df,ns,tmin=500,tmax=510, fmin=5,fmax=70
;  If i eq 0 then standard_plots,data,data*0.,time,freq,'Raw NenuFAR Data',/NO_ZOOM,xunit='sec',panel_ps=1  
;  If i ne 0 then standard_plots,data,data*0.,time,freq,'Raw NenuFAR Data',/NO_ZOOM,xunit='sec',panel_ps=1,/ONLY_PLOT
;  close,/all
;endfor

;print,'Number of f Before:',nf 

;data_or = data
;freq_or = freq
;If keyword_set(REMOVE_BADF) then begin
;  channels_per_subband     = 16.0/2.0d
;  nof_subbands             = 151.0d*2.0
 ; remove_badf,data,freq,nof_subbands,channels_per_subband,xFix=x,fFix=freq
;endif

;w = where(freq gt 20 and freq lt 22)
;w_or = where(freq_or gt 20 and freq_or lt 22)
;print,'Number of f After:',nf 
;standard_plots,data_or(*,w_or),data_or(*,w_or)*0.,time,freq_or(w_or),'Raw NenuFAR Data ',/NO_ZOOM,xunit='sec',panel_ps=1
;standard_plots,x(*,w),x(*,w)*0.,time,freq(w),'Raw NenuFAR Data (Remove Bad F)',/NO_ZOOM,xunit='sec',panel_ps=1
n_step = 0
step_size = 30
for i=0,n_step do begin
  read_nenufar, file, data,time,freq,beam,ndata, nt,dt,nf,df,ns, tmin=i*step_size,tmax=(i+1)*step_size,fmin=28,fmax=30,/VERBOSE
  standard_plots,data,data*0.,time,freq,'NenuFAR: Raw',/NO_ZOOM,xunit='sec',panel_ps=1
endfor

i = 0 
file = '/databf2/nenufar-tf/20190328_001800_20190328_043400_TAU_BOOTIS_TRACKING_BHR/TAU_BOOTIS_TRACKING_BHR_20190328_001833'
beam = 2
for i=0,n_step do begin
read_nenufar, file, data,time,freq,beam,ndata, nt,dt,nf,df,ns, tmin=i*step_size,tmax=(i+1)*step_size,fmin=28,fmax=30,/VERBOSE,/REMOVE_BADF,$
  /CORRECT_BANDPASS,/Cedric_BANDPASS
standard_plots,data,data*0.,time,freq,'NenuFAR: Cedri-derived Bandpass',/NO_ZOOM,xunit='sec',panel_ps=1
endfor 

;n_steps = 26
;for i=0,n_steps-1 do begin
;  print, '*******Step *******', i
;  READ_NU_SPEC, file,data,time,freq,beam,ndata, nt,dt,nf,df,ns,tmin=i*120,tmax=(i+1)*120, fmin=10,fmax=70
;  If i eq 0 then standard_plots,data,data*0.,time,freq,'Raw NenuFAR Data',/NO_ZOOM,xunit='sec',panel_ps=1
;  If i ne 0 then standard_plots,data,data*0.,time,freq,'Raw NenuFAR Data',/NO_ZOOM,xunit='sec',panel_ps=1,/ONLY_PLOT
;  close,/all
;endfor



;READ_NU_SPEC, file,data,time,freq,beam,ndata, nt,dt,nf,df,ns,tmin=0,tmax=4, fmin=10,fmax=40,$
;  NCHAN=NCHAN,NBEAMLETS=NBEAMLETS,/VERBOSE 
;standard_plots,data,data*0.,time,freq,'Raw NenuFAR Data',/NO_ZOOM,xunit='sec',panel_ps=1
;label='Raw NenuFAR Data'
;REDUCE_ARRAY, data, [1,nf], x_t
;REDUCE_ARRAY, data, [nt,1], x_f


;cgplot,time,x_t,xtit='Time (sec)',ytit='Intensity !C (SEFD) ',tit='Time Series: '+label,$
;  yrange=[min(x_t),max(x_t)],/ystyle,/xsty
  
;for i=0,17 do begin
;  print, x_t[(167+42*i)-5:(167+42*i)+5]
;endfor


;cgplot,freq,x_f,xtit='Frequency (MHz)',ytit='Intensity !C (SEFD)',tit='Integrated Spectrum: '+label,$
;  yrange=[min(x_f),max(x_f)],/ystyle,/xstyle


set_plot,'PS'
device,/close
cgfixps,filename_ps+'.ps'
spawn,'ps2pdf '+filename_ps+'.ps'
spawn,'rm '+filename_ps+'.ps'


end