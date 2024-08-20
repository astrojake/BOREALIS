;2d surface fit 
st= SYSTIME(1)
startmem_RFI =  MEMORY(/CURRENT)
print,'Start Memory Gb:' ,startmem_RFI/1d9
set_plot,'PS'
device,filename='surfacefit_Offringa.ps',/landscape

READ_H5_DATA, 'L429868_SAP000_B000_S0_P000_bf',0,4001,26,74,'0',x,t,f,nt,nf,nnt,nnf,dt,/NSPECTRA
STANDARD_PLOTS, x,x,t,f,'Data (RAW)',/ONLY_PLOT
;SPDYNPS, x, min(t),max(t),min(f),max(f),$
;  'Time (sec)','Frequency (MHz)','Data (Raw)',0,0,0,.05,.95,0,'.'

  x=reform(temporary(x),nt,64,244)
  x=x(*,1:62,*)
  x=reform(temporary(x),nt,62*244)
  f=reform(temporary(f),64,244)
        f=f(1:62,*)
  f=reform(temporary(f),62*244)
  nf = n_elements(f)
;restore,filename='/data/jake.turner/exoplanet/LC5_DDT_002/L429868/all_beams/mask_full_L429868_pol0_15770sec_RFI_patrol5.5_pex_lesig3.5_sum_beam0.sav'
;p2_full = p2_full[0:4000,*]
;save,p2_full,filename='p2_reduce.sav'
;restore, filename='p2_reduce.sav'
restore, filename='mask_full_L429868_pol0_41sec_RFI_patrol5.5_pex_lesig3.5_sum_beam0.sav'
print,'Done Restore'

xnew =x*p2_full
  
sigma_t = 10
sigma_f =10 
kernal = GAUSSIAN_FUNCTION([sigma_t,sigma_f])
V_1 = convol(xnew,kernal(*,0),/edge_mirror)
V_2 = convol(temporary(V_1),kernal(0,*),/edge_mirror)
V_b_1 = convol(p2_full*1.0,kernal(*,0),/edge_mirror)
v_b_2 = convol(temporary(v_b_1),kernal(0,*),/edge_mirror)
V_prime = temporary(V_2)/temporary(V_b_2)

STANDARD_PLOTS, V_prime,p2_full,t,f,'Surface'

new = x/V_prime
new = new - 1
new =new*p2_full
;STANDARD_PLOTS,new,p2,t,f,label,/correct_data
STANDARD_PLOTS, new,p2_full,t,f,'After Division',/correct_data


REDUCE_ARRAY, new, [1,nf], x_t
REDUCE_ARRAY, p2_full*1., [1,nf], p_t
x_t = x_t / (p_t + (p_t eq 0) )
REDUCE_ARRAY, new*p2_Full, [nt,1], x_f
REDUCE_ARRAY, p2_Full*1., [nt,1], p_f
x_f = x_f / (p_f + (p_f eq 0) )
LE_AUTO_S,x_f,101,4.5,0,x_f,pnet


p=rebin(reform(pnet,1,nf),nt,nf,/sample)
p2_full = p2_full*p
;new =new*p2_full
;;STANDARD_PLOTS,new,p2,t,f,label,/correct_data
;STANDARD_PLOTS, new,p2_full,t,f,'After Division',/correct_data

wbad = where(pnet eq 0)
wgood = where(pnet eq 1)
w_n = n_elements(wbad)
v = intarr(w_n)
for i=0,w_n-1 do begin
  v[i] = value_locate(wgood,wbad[i])
  V_prime(*,wbad[i]) = V_prime(*,v[i])
endfor

STANDARD_PLOTS, V_prime,p2_full,t,f,'Surface New'

new = x/V_prime
new = new - 1
new =new*p2_full
;STANDARD_PLOTS,new,p2,t,f,label,/correct_data
STANDARD_PLOTS, new,p2_full,t,f,'After Division: New Surface',/correct_data

mem_apply = (MEMORY(/HIGHWATER))/1d9
print,'SYSTIME entire program=',SYSTIME(1) - st, ' sec'
print,'Memory Gb', mem_apply
device,/close
cgfixps, 'surfacefit_Offringa.ps'
stop
;kernal = [1,2,3,2,1]/10.

;kernelSize = [20, 20]
;kernel = REPLICATE((1./(kernelSize[0]*kernelSize[1])),kernelSize[0], kernelSize[1])
;savgolFilter = SAVGOL(10, 10, 0, 4, /DOUBLE )


print, 'Start Convol'
surface = convol(xnew,kernal,invalid=0,/edge_mirror)
print,'End Convol'
;STANDARD_PLOTS, p2_full,p2_full,t,f,'Mask',/mask,/ONLY_PLOT
STANDARD_PLOTS, surface,p2_full,t,f,'Surface'

;REDUCE_ARRAY, surface, [nt,1], surface_f
;LE_AUTO_S,surface_f,101,3.5,0,x_f,pnet 
REDUCE_ARRAY, surface, [nt,1], surface_f
LE_AUTO_S,surface_f,101,4.5,0,surface_f,pnet 
p=rebin(reform(pnet,1,nf),nt,nf,/sample)
p2_full = p2_full*p


xnew =x*p2_full

sigma_t = 50
sigma_f = 200
kernal = GAUSSIAN_FUNCTION([sigma_t,sigma_f])

print, 'Start Convol'
surface = convol(xnew,kernal,invalid=0,/normalize,/edge_mirror)
print,'End Convol'
;STANDARD_PLOTS, p2_full,p2_full,t,f,'Mask',/mask,/ONLY_PLOT
STANDARD_PLOTS, surface,p2_full,t,f,'Surface new 0'

;SPDYNPS,p2_full, min(t),max(t),min(f),max(f),$
;  'Time (sec)','Frequency (MHz)','Mask',0,0,0,0,1,1,'.'
;SPDYNPS, result, min(t),max(t),min(f),max(f),$
;  'Time (sec)','Frequency (MHz)','Surface',0,0,0,.05,.95,0,'.'

new = x/surface
new = new - 1 
new =new*p2_full
;STANDARD_PLOTS,new,p2,t,f,label,/correct_data
 STANDARD_PLOTS, new,p2_full,t,f,'After Division',/correct_data
;  SPDYNPS, new*p2_full, min(t),max(t),min(f),max(f),$
;    'Time (sec)','Frequency (MHz)','Result',0,0,0,.05,.95,0,'.'
;STANDARD_PLOTS, p2_full(0:100,*),p2_full(0:100,*),t(0:100,*),f,'Mask (Zoom)',/mask,/ONLY_PLOT
;STANDARD_PLOTS, new(0:100,*),p2_full(0:100,*),t(0:100),f,'After Division (Zoom)',/ONLY_PLOT

  REDUCE_ARRAY, new, [1,nf], x_t
  REDUCE_ARRAY, p2_full*1., [1,nf], p_t
  x_t = x_t / (p_t + (p_t eq 0) )
  REDUCE_ARRAY, new*p2_Full, [nt,1], x_f
  REDUCE_ARRAY, p2_Full*1., [nt,1], p_f
  x_f = x_f / (p_f + (p_f eq 0) )
;  plot,t,x_t,/xsty,xtit='Time (sec)',ytit='Frequency-integrated values',tit=label,$
;    yrange=[min(x_t),max(x_t)],/ystyle
;  plot,f,x_f,/xsty,xtit='Frequency (MHz)',ytit='Time-integrated values',tit=label,$
;    yrange=[min(x_f),max(x_f)],/ystyle
   
   LE_AUTO_S,x_f,101,4.5,0,x_f,pnet 
   
;   plot,f,x_f,/xsty,xtit='Frequency (MHz)',ytit='Time-integrated values',tit=label,$
;     yrange=[min(x_f),max(x_f)],/ystyle
   
   p=rebin(reform(pnet,1,nf),nt,nf,/sample)
   p2_full = p2_full*p 
   new =new*p2_full
   ;STANDARD_PLOTS,new,p2,t,f,label,/correct_data
   STANDARD_PLOTS, new,p2_full,t,f,'After Division',/correct_data
   
   wbad = where(pnet eq 0)
   wgood = where(pnet eq 1)
   w_n = n_elements(wbad)
   v = intarr(w_n)
   for i=0,w_n-1 do begin
     v[i] = value_locate(wgood,wbad[i])
     surface(*,wbad[i]) = surface(*,v[i])
   endfor

   STANDARD_PLOTS, surface,p2_full,t,f,'Surface New'
   
   new = x/surface
   new = new - 1
   new =new*p2_full
   ;STANDARD_PLOTS,new,p2,t,f,label,/correct_data
   STANDARD_PLOTS, new,p2_full,t,f,'After Division: New Surface',/correct_data

device,/close
cgfixps, 'surfacefit_Offringa.ps'
mem_apply = (MEMORY(/HIGHWATER))/1d9
print,'SYSTIME entire program=',SYSTIME(1) - st, ' sec'
print,'Memory Gb', mem_apply
end