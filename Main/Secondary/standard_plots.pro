;----------------------------------
 pro STANDARD_PLOTS, x,pp,t,f,label,mask=mask,correct_data = correct_data,VERBOSE=VERBOSE,$
	ONLY_PLOT=ONLY_PLOT,histogram=histogram,NO_ZOOM=NO_ZOOM,skyplot=skyplot,xunit=xunit,$
	fracmin=fracmin,fracmax=fracmax,ytitle_bar=ytitle_bar,xlabel=xlabel,panel_ps=panel_ps,$
	backgr=backgr,PLOTMASK=PLOTMASK,DYNSPECMS=DYNSPECMS,multi_page=multi_page
;----------------------------------
; x(t,f)       ; data
; pp           ; mask(t,f) same size as data
; t            ; array of elements in time
; f            ; array of elements in frequency
; label        ; label for the plots
;
; Keywords:
;  Mask        : The x data is the mask and therefore will not create histogram for x
;  correct data: Display the data correctly with the mask using reduce_array 
;  ONLY_PLOT   : Only plot the dynamic spectrum and nothing else
;  correct_data ; correct the data with the mask 
;  NO_zoom      ; don't show zoom
;  sky_plot     ; plot sky below target 
;  xunit        ;unit of the the time axis
;  ytitle_bar    ; title of bar 
;  ;fracmin,fracmax =  (INPUT)  min and max threshold for optimization of the dynamic range of the plot
set_plot,'PS'
;If keyword_set(panel_ps) and not(keyword_set(multi_page)) then !p.multi=[0,1,2] else !p.multi=[0,1,1]
If keyword_set(panel_ps) then !p.multi=[0,1,2] else !p.multi=[0,1,1]

If not(keyword_set(backgr)) then backgr=0  ;no background subtract
If not(keyword_set(xunit)) then xunit = 'time unit'
If not(keyword_set(xlabel)) then xlabel='Time'
If not(keyword_set(fracmin)) then fracmin = 0.05
If not(keyword_set(fracmax)) then fracmax = 0.95
If not(keyword_set(ytitle_bar)) then ytitle_bar = 'Intensity (Arbitrary units)'
nt=n_elements(t) & nf=n_elements(f)
ntmax=3000.0d       & nfmax=3000.0d
ntmax=1500.0d       & nfmax=1500.0d 

nnt=ceil(nt/ntmax) & nnnt=long(nt/nnt)
if nnt ne long(nnt) then nnt=long(nnt)+1 else nnt=long(nnt)
nnf=ceil(nf/nfmax) & nnnf=long(nf/nnf)
if nnt ne long(nnf) then nnf=long(nnf)+1 else nnt=long(nnf)

nnv=nt/ntmax        ; for visualization in ps file
if nnv ne long(nnv) then nnv=long(nnv)+1 else nnv=long(nnv)
nntv=long(nt/nnv)
mmv = nf/nfmax
if mmv ne long(mmv) then mmv=long(mmv)+1 else mmv=long(mmv)
mmtv=long(nf/mmv)

If keyword_set(DYNSPECMS) then begin
  ntmax = nt 
  nfmax = nf
  print, 'nt',nt
  print, 'ntmax',ntmax
  print, 'nf',nf
  print,'nfmax',nfmax
Endif

If nt gt ntmax or nf gt nfmax then begin
   If not(keyword_set(mask)) then begin
      reduce_array,x,[nnv,mmv],x_r,/dbl
      SPDYNPS,x_r, min(t),max(t),min(f),max(f),xlabel+' ('+xunit+')','Frequency (MHz)',label,$
      0,backgr,0,fracmin,fracmax,0,'.',ytitle_bar=ytitle_bar
  
    If keyword_set(PLOTMASK) then begin
      reduce_array,pp,[nnv,mmv],pp_r,/dbl
      SPDYNPS, pp_r, min(t),max(t),min(f),max(f),xlabel+' ('+xunit+')','Frequency (MHz)',label,0,0,0,0,1,1,'.',ytitle_bar='Mask'
    Endif
  
    If keyword_set(skyplot) then begin
     print, 'Sky Plot'
     reduce_array,skyplot,[nnv,mmv],sky_r,/dbl
     SPDYNPS,sky_r, min(t),max(t),min(f),max(f),xlabel+' ('+xunit+')','Frequency (MHz)','Sky',$
       0,backgr,0,fracmin,fracmax,0,'.',ytitle_bar=ytitle_bar
    Endif  
  Endif
  If keyword_set(mask) then begin
    reduce_array,x,[nnv,mmv],x_r,/dbl
    SPDYNPS, x_r, min(t),max(t),min(f),max(f),xlabel+' ('+xunit+')','Frequency (MHz)',label,0,0,0,0,1,1,'.',ytitle_bar='Mask'
  Endif  
Endif

If nt le ntmax and nf le nfmax then begin
  If not(keyword_set(mask)) then begin
    SPDYNPS, x, min(t),max(t),min(f),max(f),xlabel+' ('+xunit+')','Frequency (MHz)',label,0,backgr,0,fracmin,fracmax,0,'.',ytitle_bar=ytitle_bar
    
    If keyword_set(PLOTMASK) then begin
      reduce_array,pp,[nnv,mmv],pp_r,/dbl
      SPDYNPS, pp_r, min(t),max(t),min(f),max(f),xlabel+' ('+xunit+')','Frequency (MHz)',label,0,0,0,0,1,1,'.',ytitle_bar='Mask'
    Endif
    
    If keyword_set(skyplot) then begin
      print, 'Sky Plot'
      SPDYNPS,skyplot, min(t),max(t),min(f),max(f),xlabel+' ('+xunit+')','Frequency (MHz)','',$
        0,backgr,0,fracmin,fracmax,0,'.',ytitle_bar=ytitle_bar
    Endif
  Endif

  If keyword_set(mask) then begin
    SPDYNPS, x, min(t),max(t),min(f),max(f),xlabel+' ('+xunit+')','Frequency (MHz)',label,0,0,0,0,1,1,'.',ytitle_bar='Mask'
  Endif
Endif

If not(keyword_set(NO_ZOOM)) then begin 
  if nt gt ntmax or nf gt nfmax then begin
    xx=x(((nt-ntmax)/2)>0:((nt+ntmax)/2)<(nt-1),((nf-nfmax)/2)>0:((nf+nfmax)/2)<(nf-1))
    tt=t(((nt-ntmax)/2)>0:((nt+ntmax)/2)<(nt-1))
    ff=f(((nf-nfmax)/2)>0:((nf+nfmax)/2)<(nf-1))
    If not(keyword_set(mask)) then begin
      SPDYNPS, xx, min(tt),max(tt),min(ff),max(ff),xlabel+' ('+xunit+')','Frequency (MHz)','Dynamic Spectrum '+label+' (zoom)',0,backgr,0,.05,.95,0,'.',ytitle_bar=ytitle_bar
    endif
    
    If keyword_set(mask) then begin
      SPDYNPS, xx, min(tt),max(tt),min(ff),max(ff),xlabel+' ('+xunit+')','Frequency (MHz)',label+' (zoom)',0,backgr,0,0.0,1.0,1,'.',ytitle_bar='Mask'
    Endif
    
  endif
endif 

If not(keyword_set(ONLY_PLOT)) then begin    
  If not(keyword_set(correct_data)) then begin
   REDUCE_ARRAY, x, [1,nf], x_t
   ;reduce_Array,x_t,[nnt],x_t
   ;reduce_array,t,[nnt],t
   If keyword_set(VERBOSE) then begin
     print,'Stdev of Time Series: ',stddev(x_t)
   endif 
   REDUCE_ARRAY, x, [nt,1], x_f
   If keyword_set(VERBOSE) then begin
     print,'Stdev of Time-integrated: ',stddev(x_f)
   endif 
   cgplot,t,x_t,xtit=xlabel+' ('+xunit+')',ytit='Intensity !C (SEFD) ',tit='Time Series: '+label,$
    yrange=[min(x_t),max(x_t)],/ystyle,/xsty    
   cgplot,f,x_f,xtit='Frequency (MHz)',ytit='Intensity !C (SEFD)',tit='Integrated Spectrum: '+label,$
    yrange=[min(x_f),max(x_f)],/ystyle,/xstyle
  Endif

 If keyword_set(correct_data) then begin
    REDUCE_ARRAY, x*pp, [1,nf], x_t
    REDUCE_ARRAY, pp*1., [1,nf], p_t
    x_t = x_t / (p_t + (p_t eq 0) )
   If keyword_set(VERBOSE) then begin
     print,'Stdev of Time Series: ',stddev(x_t)
   endif
   
    REDUCE_ARRAY, x*pp, [nt,1], x_f
    REDUCE_ARRAY, pp*1., [nt,1], p_f
    x_f = x_f / (p_f + (p_f eq 0) )
   If keyword_set(VERBOSE) then begin
     print,'Stdev of Time-integrated: ',stddev(x_f)
   endif 
       
    cgplot,t,x_t,/xsty,xtit=xlabel+' ('+xunit+')',ytit='Frequency-integrated values',tit=label,$
      yrange=[min(x_t),max(x_t)],/ystyle
    cgplot,f,x_f,/xsty,xtit='Frequency (MHz)',ytit='Time-integrated values',tit=label,$
      yrange=[min(x_f),max(x_f)],/ystyle
 Endif 

 If keyword_set(histogram) then begin
  ;only take good data for histogram
  wn=where(pp ne 0) 
  
  mx=mean(x(wn),/double) & sx=stddev(x(wn),/double)
  hx=histogram((x(wn)-mx)/sx,min=-5,max=10,nbins=301,loc=xh)
  xh=xh+0.025
  hxm=exp(-xh^2/2)
  hxm=(hxm/total(hxm,/double)*total(hx,/double))
  If not(keyword_set(mask)) then begin
      plot_io,xh,hx>1,psym=10,/xsty,xtit='Standard deviations',ytit='N',tit=label
      oplot,xh,hxm,thick=6.,line=2
      plot_io,xh,((hx>1)/(hxm>1)),psym=10,/xsty,xtit='Standard deviations',ytit='ratio'
      oplot,xh,xh^0,line=1
  endif 
 endif  
endif 

If keyword_set(panel_ps) then !p.multi=0
;If keyword_set(panel_ps) and not(keyword_set(multi_page)) then !p.multi=0
return
;!p.multi=[0,1,1]
end
