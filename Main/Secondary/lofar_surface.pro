;sigma_t = 5
;sigma_f =10

;From Offringa's thesis: Section 2.2
;Outputs: V_prime 

;restore,filename='mask_full_bit_L527649_pol0_10779sec_RFI_patrol5.5_pex_lesig3.5_sum_beam2.sav'

pro lofar_surface,x,p2,t,f,sigma_t,sigma_f,V_prime=V_prime
  p2_s = p2
  ;p2(*,*) = 1.0 ;flag 
  
  st= SYSTIME(1)
  mem_apply = (MEMORY(/HIGHWATER))/1d9 ; Gb
  print,'Memory Before',mem_apply
  
  nt = n_elements(t)
  nf = n_elements(f)
 
  x =x
  Wd = GAUSSIAN_FUNCTION([sigma_t,sigma_f])
  nt_wd = n_elements(wd(*,0))
  nf_wd = n_elements(wd(0,*)) 
  Uv = Wd(*,nt_wd/2) 
  Ut = Wd(nf_wd/2,*)
  
  V_1 = convol(x,Uv,/edge_mirror)              ;Wf x V 
  V_2 = convol(temporary(V_1),Ut,/edge_mirror)
  V_b_1 = convol(p2*1.0,Uv,/edge_mirror)
  v_b_2 = convol(temporary(v_b_1),Ut,/edge_mirror)
  V_prime = temporary(V_2)/temporary(V_b_2) 

  STANDARD_PLOTS, p2,p2,t,f,'Mask',/mask
  STANDARD_PLOTS, V_prime,p2,t,f,'Surface'
  
  ;x = (x*p2)/(p2 + (p2 eq 0))
  d = x/V_prime -1.
  d = d*p2_s
  
  STANDARD_PLOTS,d,p2_s,t,f,'Data/Surface',/correct_data
  
  If keyword_set(MORE) then begin
    ;xnew = x/V_prime
    ;xnew = xnew - 1
    ;xnew =xnew*p2
    xnew = d
    
    REDUCE_ARRAY, xnew, [1,nf], x_t
    REDUCE_ARRAY, p2*1., [1,nf], p_t
    x_t = x_t / (p_t + (p_t eq 0) )
    REDUCE_ARRAY, xnew*p2, [nt,1], x_f
    REDUCE_ARRAY, p2*1., [nt,1], p_f
    x_f = x_f / (p_f + (p_f eq 0) )
    
    ;RFI
    LE_AUTO_S,x_f,101,4.5,0,x_f,pnet
    p=rebin(reform(pnet,1,nf),nt,nf,/sample)
    p2 = p2*p
    
    wbad = where(pnet eq 0)
    wgood = where(pnet eq 1)
    w_n = n_elements(wbad)
    v = intarr(w_n)
    for i=0,w_n-1 do begin
      v[i] = value_locate(wgood,wbad[i])
      V_prime(*,wbad[i]) = V_prime(*,v[i])
    endfor
    
    STANDARD_PLOTS, V_prime,p2,t,f,'Surface Updated'
  endif 
  
  ;save,V_prime,t,f,filename='surafce.sav'
  device,/close
  mem_apply = (MEMORY(/HIGHWATER))/1d9 ; Gb
  print,'Memory After (Gb): ',mem_apply
  print,'SYSTIME lofar surface program=',SYSTIME(1) - st, ' sec'
end