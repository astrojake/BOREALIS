;READ UTR-2 sav sets from Ian's pipeline
;
;Inputs: root  ;root folder (e.g.)
;        file  ;name of the file that contains the folders (e.g. foldernames.dat)
;pro read_utr2_sav,root,file

;Contents of interpointing save sets
;  mode - {0-Waveform, 1-Sum/Dif, 2-Correlation}
;  ndim - number of dimensions in data: n_elements(data(*,0,0))
;  ch {0,1,2} - what part of data was used (e.g. Sum or Dif or both Sum & Dif)
;  ip_size - (ndim,nf,nt)
;  new_df - new frequency resolution of reduced data
;  new_dt - new frequency resolution of reduced data
;  ff - reduction factor in frequency
;  ft - reduction factor in time
;  cat - RFI contamination category
;  red_data - reduced data array
;  red_freq_ramp - reduced frequency ramp
;  red_time_ramp - reduced time ramp
;  mask_original_bit - array, containing masks for each data dimension (e.g. Sum, Dif) after cleaning
;  red_mask_merged - reduced merged mask, not in bits
start_mem = MEMORY(/CURRENT) ;start memory
st= SYSTIME(1)               ;start time

filename_ps  = 'UTR2_test'
cgPS_Open,filename_ps+'.ps'
PS = 1 

;-------
;Inputs
;-------
root        ='/databf/lofar/EXO/utr-2/136-Exoplanet-TauBootes-P144/'
file        ='foldernames.dat'   ;filename of all folders (e.g. 260307_020315) 
VERBOSE = 1 
;-------------
;Analysis
;-------------

;Read folders 
readcol,file,folders,format='A'
date =STRMID(folders[0], 0, 6)
n_folders = n_elements(folders)

n_ip_all = dblarr(n_folders)
;Find the total number of interpointings 
 for j=0,n_folders-1  do begin ;loop over folders
  restore,filename=root+folders[j]+'/'+folders[j]+'.sav'
  ;ipnum-1     ;number of interpoints
  n_ip_all[j] = ipnum - 1
endfor  

ntot_ip = total(n_ip_all)
data_ON_p    = ptrarr(ntot_ip)
data_OFF_p   = ptrarr(ntot_ip)
mask_p       = ptrarr(ntot_ip) 
c_ip = 0 


for j=0,n_folders-1  do begin ;loop over folders
  Print,'*******   Folder  *****: ', folders[j]
  Print,'*******   Total IP  *****: ', n_ip_all[j]
 for i=0,n_ip_all[j]-1  do begin ;loop over ip 
  print,'******* Ip Number *****: ', i 
  If i lt 10 then restore,filename=root+folders[j]+'/ip=0'+strtrim(string(uint(i)),1)+'.sav'
  If i ge 10 then restore,filename=root+folders[j]+'ip='+strtrim(string(uint(i)),1)+'.sav' 
  If keyword_set(VERBOSE) and i eq 0 then begin
    print, 'File: ', root+folders[0]
    print, 'Interpointing: ',i
    print, 'Mode (0-Waveform, 1-Sum/Dif, 2-Correlation)', mode
    print, 'Dimensions of the org data',ip_size
    print, '# of freqs in reduced data', n_elements(red_freq_ramp)
    print, '# of time steps in reduced data', n_elements(red_time_ramp)
  Endif
  ;-------------------
  ;Read Sum/Diff Data 
  ;------------------
  If mode eq 1 then begin
    SUM = transpose(RED_DATA[0,*,*])   ; (f,t)  -- (t,f)
    DIFF = transpose(RED_DATA[1,*,*]) 
    MJD_t = temporary(red_time_ramp - 2400000.5d)
    If j eq 0 and i eq 0 then tmin = 0
    If j ne 0 or i ne 0 then begin
      ;diff_t = (min(MJD_t) - max_mjd_t)*86400.0d   ; diff in secs between ips
      ;tmin = tmax + diff_t
      tmin = tmax
    Endif
    Max_MJD_t = max(MJD_t) 
    nt = n_elements(mjd_t) 
    t = dindgen(nt,start=tmin,increment=new_dt) ;time array (s)
    tmax = max(t)                               ;max time (s) 
    f = red_freq_ramp 
    nf = n_elements(f)
    mask = transpose(red_mask_merged)
     
        
    If j eq 0 and keyword_set(PS) then begin
;       STANDARD_PLOTS, SUM,mask,t,f,'SUM: Raw',/correct_data
;       MAKE_BACKGROUND,sum*mask,'', data_mean,ss,nn
;       cgplot,f,data_mean,title='Background',xtitle='Freq',ytitle='flux'
;       STANDARD_PLOTS, mask,mask,t,f,'Mask',/ONLY_PLOT,/MASK
;       ;data_mean = smooth(data_mean,12,/EDGE_TRUNCATE)
;       SUM_norm = SUM/rebin(reform(data_mean,1,nf),nt,nf)
;       ;cgplot,f,data_mean,title='10% Quantile',xtitle='Freq',ytitle='flux'
;       STANDARD_PLOTS, SUM_norm,mask,t,f,'SUM Normalized',/ONLY_PLOT
    Endif
      
    ;Create Beams
    ON                = SUM - DIFF
    OFF               = DIFF
    
    If j eq 0 and keyword_set(PS) then begin
      STANDARD_PLOTS, SUM/(mask + (mask eq 0)),mask,t,f,'SUM: Raw',/ONLY_PLOT
      STANDARD_PLOTS, ON/(mask + (mask eq 0)),mask,t,f,'ON (Sum-Diff): Raw',/ONLY_PLOT
      STANDARD_PLOTS, OFF/(mask + (mask eq 0)),mask,t,f,'Diff (OFF): Raw',/ONLY_PLOT
    endif   
    
    FIX_Background = 1
    If keyword_set(FIX_BACKGROUND) then begin
      reduce_array,SUM,[nt,1],data_mean_SUM
      reduce_array,ON,[nt,1],data_mean_ON
      reduce_array,OFF,[nt,1],data_mean_OFF
      reduce_array, mask, [nt,1],mask_mean
      If j eq 0 then begin
        cgplot,f,data_mean_SUM,xtitle='Freq',ytitle='Flux',title='Mean Background'
        cgplot,f,data_mean_OFF,/overplot,color='red' ; OFF
        cgplot,f,data_mean_ON,/overplot,color='blue' ;ON 
        Al_legend,['SUM','DIFF','ON (SUM-DIFF)'],textcolor=['black','red','blue']
        ;cgplot,f,data_mean_OFF,xtitle='Freq',ytitle='Flux',title='Mean Background; OFF'
        ;cgplot,f,data_mean_ON,xtitle='Freq',ytitle='Flux',title='Mean Background; ON'
      Endif
      data_mean_SUM = data_mean_SUM/(mask_mean + (mask_mean eq 0))
      data_mean_OFF = data_mean_OFF/(mask_mean + (mask_mean eq 0))
      data_mean_ON = data_mean_ON/(mask_mean + (mask_mean eq 0)) 
      If j eq 0 then begin
        cgplot,f,data_mean_SUM,xtitle='Freq',ytitle='Flux',title='Mean Background/Mask'
        cgplot,f,data_mean_OFF,/overplot,color='red' ; OFF
        cgplot,f,data_mean_ON,/overplot,color='blue' ;ON
        Al_legend,['SUM','DIFF','ON (SUM-DIFF)'],textcolor=['black','red','blue']
      Endif
     ; If j eq 0 then cgplot,f,data_mean_OFF,xtitle='Freq',ytitle='Flux',title='Mean Background; OFF (mask)'
     ; If j eq 0 then cgplot,f,data_mean_ON,xtitle='Freq',ytitle='Flux',title='Mean Background; ON (mask)' 
    Endif
    
    
    SUM = SUM/rebin(reform(data_mean_SUM,1,nf),nt,nf)
    OFF = OFF/rebin(reform(data_mean_OFF,1,nf),nt,nf)
    ON = ON/rebin(reform(data_mean_ON,1,nf),nt,nf)
    
    If j eq 0 and keyword_set(PS) then begin
      STANDARD_PLOTS, SUM/(mask + (mask eq 0)),mask,t,f,'SUM Normalized',/ONLY_PLOT
      STANDARD_PLOTS, OFF/(mask + (mask eq 0)),mask,t,f,'DIFF Normalized',/ONLY_PLOT
      STANDARD_PLOTS, ON/(mask + (mask eq 0)),mask,t,f,'(SUM-DIFF) Normalized',/ONLY_PLOT
      STANDARD_PLOTS, (SUM-OFF)/(mask + (mask eq 0)),mask,t,f,'(SUM/Norm) -(DIFF/Norm): ON Normalized',/ONLY_PLOT
    Endif
    
    ;After Norma
    reduce_array,ON,[nt,1],data_mean_ON
    reduce_array,OFF,[nt,1],data_mean_OFF_After
    If j eq 0 and keyword_set(PS) then begin
      cgplot,f,data_mean_OFF_After/(mask_mean + (mask_mean eq 0)),xtitle='Freq',ytitle='Flux',title='Mean Background of DIFF After Normalization'
      cgplot,f,data_mean_ON/(mask_mean + (mask_mean eq 0)),xtitle='Freq',ytitle='Flux',title='Mean Background of ON Normalization'
    endif
    
    ;have the correct ON
    ON = SUM - OFF
    
    ;-----------
    ;rebin time
    ;-----------
    new_rebin_time = 1.0d   ; 1 second
    num_t = long(new_rebin_time/new_dt)
    reduce_array,ON, [num_t,1],ON
    reduce_array,OFF,[num_t,1],OFF
    reduce_array,mask,[num_t,1],mask
    reduce_array,t,[num_t],t
    reduce_array,MJD_t,[num_t],MJD_t
  
    data_ON_p[c_ip]   = ptr_new(temporary(ON))
    data_OFF_p[c_ip]  = ptr_new(temporary(OFF))
    mask_p[c_ip]      = ptr_new(temporary(mask))
     
    If j eq 0 and i eq 0 then begin
      MJD_all = [MJD_t]
      t_all   = [t]
    Endif
    If j gt 0 or i gt 0 then begin
      MJD_all = [MJD_all,MJD_t]
      t_all   = [t_all,t]
    Endif
    c_ip = c_ip + 1  ;count ips for the pointer
  endif
 endfor ;end loop over ip    
endfor

;-----------------------------------------------
;Concatenate arrays for ON, OFF, and mask
;-----------------------------------------------
pointer_to_array,data_ON_p,out_array=data_ON
pointer_to_array,data_OFF_p,out_array=data_OFF
pointer_to_array,mask_p,out_array=mask

If keyword_set(VERBOSE) then begin
  print,'*********************************************************'
  print,'**         After Combining All IP ***********************'
  print,'*********************************************************'
  print,'# of time elements in full dataset (t_all)',n_elements(t_all)
  print,'# of time steps in combined dataset', n_elements(data_ON(*,0))
  print,'# of freq channels in combined dataset', n_elements(data_OFF(0,*))
Endif

;Rebin Arrays

;Save arrays
xt_rebin   = t_all 
xf_rebin   = f 
MJD_rebin  = MJD_all 
data_rebin = data_ON 
p2_rebin   = mask
save,data_rebin,xt_rebin,xf_rebin,MJD_rebin,p2_rebin,filename='rebindata_'+date+'_beam0.sav'

data_rebin = data_OFF
save,data_rebin,xt_rebin,xf_rebin,MJD_rebin,p2_rebin,filename='rebindata_'+date+'_beam1.sav'

STANDARD_PLOTS, data_ON/(mask + (mask eq 0)),mask,t_all,f,'ON-All Normalized',/ONLY_PLOT
STANDARD_PLOTS, data_OFF/(mask + (mask eq 0)),mask,t_all,f,'OFF-All Normalized',/ONLY_PLOT


If keyword_Set(PS) then begin
  cgPS_Close
  cgfixps,filename_ps+'.ps'
  ;Create panels pdfs
  spawn,'ps2pdf -dPDFSETTINGS=/ebook '+filename_ps+'.ps'
  spawn, 'rm '+filename_ps+'.ps'
 endif 
  
 ;clear heap memory
 HEAP_GC,/PTR,/VERBOSE 
 PRINT, 'Memory required for UTR2: ', (MEMORY(/HIGHWATER))/1d9 ,' Gb'
 print,'SYSTIME=',SYSTIME(1) - st, ' sec'
end