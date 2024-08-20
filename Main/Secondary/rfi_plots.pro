;Create plots for RFI
;Creates plots with input of data and mask 
;
;Input: 
;      data
;      p        
;      t
;      f 
;      label
;      
;
;Options
;      p_before   ;to make plot of difference of masks
;
pro rfi_plots,data,p,t,f,label=label,p_before=p_before,PS=PS,LAST_PS=LAST_PS
 
  If not(keyword_set(label)) then label = ''
  ;Set Mask
  data[where(p eq 0)]  =    median(data[where(p eq 1)])
  nf = n_elements(f)
  nt = n_elements(t) 
  
  ;Update data
  print, 'Median of Data :', median(data[where(p eq 1)])
  REDUCE_ARRAY, data*p, [1,nf], data_t,/dbl   ;create t series
  REDUCE_ARRAY, data*p, [nt,1], data_f,/dbl   ;create integrated spectrum
  print,'Stdev of Integrated Spectra: ', stddev(data_t)
  print,'Stdev of Integrated t Series: ', stddev(data_f)

  data  = data   -    mean(data[where(p eq 1)])
  If not(keyword_set(LAST_PS)) then begin
    If keyword_set(PS) then begin
      STANDARD_PLOTS, data,p*1.,t,f,'Data: '+label,VERBOSE=VERBOSE,panel_ps=1
    endif
  
    If keyword_set(PS) then begin
      STANDARD_PLOTS, p*1.,p*1.,t,f,'Mask: '+label,/mask,panel_ps=1
    endif
  
    If keyword_set(PS) then begin
      STANDARD_PLOTS, data,p*1.,t,f,'Data*Mask: '+label,/correct_data,/histogram,panel_ps=1
    endif
    
    If keyword_set(p_before) then begin
     If keyword_set(PS) then begin
      dmask = p_before*1. - p*1.
      STANDARD_PLOTS, dmask,p*1.,t,f,'Delta Mask: '+label,/mask,panel_ps=1
     endif
    endif
  endif ;LAST PS
  
  If keyword_set(LAST_PS) then begin
    !p.multi=[0,1,2]
    ntmax=2500       & nfmax=1200
    nt=n_elements(t) & nf=n_elements(f)
    nnt=ceil(nt/ntmax) & nnnt=long(nt/nnt)
    nnf=ceil(nf/nfmax) & nnnf=long(nf/nnf)
    label   = 'Data*Mask: After All RFI'
    SPDYNPS, rebin(data(0:nnnt*nnt-1,0:nnnf*nnf-1),nnnt,nnnf), min(t),max(t),min(f),max(f),'Time (sec)','Frequency (MHz)',label,0,0,0,.05,.95,0,'.'
    label = 'Mask'
    SPDYNPS, rebin(p(0:nnnt*nnt-1,0:nnnf*nnf-1)*1.,nnnt,nnnf), min(t),max(t),min(f),max(f),'Time (sec)','Frequency (MHz)',label,0,0,0,0.0,1.0,1,'.'
  endif

return  
      
end
