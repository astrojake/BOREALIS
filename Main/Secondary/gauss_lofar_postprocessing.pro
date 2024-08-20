;Gauss similation for post-processing of lofar
;not general
pro gauss_lofar_postprocessing,gauss_steps,array_size,threshold,slope=slope,Q3a=Q3a,dQ3a=Q3a_diff,Q3b=Q3b,dQ3b=Q3b_diff,Q3c=Q3c,dQ3c=Q3c_diff,$
                               Q3d=Q3d,dQ3d=Q3d_diff,Q3e=Q3e,dQ3e=Q3e_diff,Q3f=Q3f,dQ3f=Q3f_diff,square=square
  ;start time
  st= SYSTIME(1)
  
  If not(keyword_set(slope)) then slope = 2.0 
  
  ;-----------------------------------------------------------
  ;                       Guassian
  ;-----------------------------------------------------------
  If not(keyword_set(gauss_steps)) then gauss_steps = 10000
  
  count_3a  = dblarr(gauss_steps)    & count_3a2 =   dblarr(gauss_steps)
  count_3c  = dblarr(gauss_steps)    & count_3c2 =   dblarr(gauss_steps)
  
  count_3e     = dblarr(gauss_steps)    
  count_3e_2   = dblarr(gauss_steps)    
  Power_3b  = dblarr(gauss_steps)   &Power_3b2  = dblarr(gauss_steps)
  Power_3f  = dblarr(gauss_steps)    &Power_3d  =  dblarr(gauss_steps)
  Power_3d2 = dblarr(gauss_steps)
   Power_3f_2  = dblarr(gauss_steps)   
   
  for i=0,gauss_steps-1 do begin
    gauss   = RANDOMN(seed,array_size)
    If keyword_set(square) then gauss = gauss^2.0d
    gauss2  = RANDOMN(seed,array_size)
    If keyword_set(square) then gauss2 = gauss2^2.0d
    gauss3  = RANDOMN(seed,array_size)
    If keyword_set(square) then gauss3 = gauss3^2.0d
    gauss4  = RANDOMN(seed,array_size)
    If keyword_set(square) then gauss4 = gauss4^2.0d
    
    index_Ag     = where(gauss ge threshold,count1)    ;number of peaks above threshold
    index_Ag_2   = where(gauss2 ge threshold,count1_2) ;# of peaks above threshold (gaussian 2)
    index_Bg     = where(gauss LE -threshold, count2)   ;number of peaks below threshold
    index_bg_2   = where(gauss2 LE -threshold, count2_2)   ;number of peaks below threshold (gaussian 2)
    
    index_Q3e    = WHERE(gauss ge threshold and gauss ge slope*gauss2, count3)  ;number of peaks above threshold and gt than off values by factor of  2
    index_Q3e_2  =  WHERE(gauss3 ge threshold and gauss3 ge slope*gauss4, count3_2) ;another 2 guassians
     
    count_3a(i)  = count1   ;3a
    count_3a2(i) = count1_2 ;3a (2) 
    
    count_3c(i)  = count2   ;3c
    count_3c2(i) = count2_2 ;3c(2)
    
    count_3e(i)  = count3 ;3e
    count_3e_2(i)  = count3_2 ;3e
    
    If count1 gt 0 then   Power_3b(i)  = total(gauss[index_Ag],/double)       ;Guassian Power of peaks above threshold
    If count1_2 gt 0 then Power_3b2(i) = total(gauss2[index_Ag_2],/double) ;Guassian Power of peaks above threshold (Gauss 2)
    If count2 gt 0 then Power_3d(i) =  total(abs(gauss[index_Bg]),/double )
    If count2_2 gt 0 then Power_3d2(i) = total(abs(gauss2[index_Bg_2]),/double )
    If count3 gt 0  then begin
      Power_3f(i) = total(gauss[index_Q3e],/double)      ;power of 3f (Guassian)
    Endif
    If count3_2 gt 0  then begin
      Power_3f_2(i) = total(gauss3[index_Q3e_2],/double)      ;power of 3f (Guassian)
    Endif
  endfor ;gauss steps = 10000

  Q3a              = mean(count_3a)                       ;3a
  Q3a_diff         = stddev(count_3a - count_3a2)         ;3a diff
  Q3b              = mean(Power_3b)                       ;3b
  Q3b_diff         = stddev(Power_3b - Power_3b2)         ;3b_diff
  Q3c              = stddev(count_3a - count_3c)  ;3c
  Q3c_diff         = stddev((count_3a - count_3c) - (count_3a2 - count_3c2))
  Q3d              = stddev(Power_3b - Power_3d)  ;3d
  Q3d_diff         = stddev((Power_3b - Power_3d) - (Power_3b2 - Power_3d2))
  Q3e              = mean(count_3e)                       ;3e
  Q3e_diff         = stddev(count_3e - count_3e_2)         ;3e Diff 
  Q3f              = mean(Power_3f)                       ;3f
  Q3f_diff         = stddev(Power_3f- Power_3f_2)                       ;3f
end
