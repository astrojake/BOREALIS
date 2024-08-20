;Guass similation for post-processing of lofar
;not general
pro gauss_lofar_postprocesssing,gauss_steps,array_size,threshold,Q3a=Q3a,Q3a_diff=Q3adiff,Q3b=Q3b,Q3bdiff=Q3b_diff,Q3c=Q3c,Q3d=Q3d,Q3e=Q3e,Q3f=Q3f
  ;start time
  st= SYSTIME(1)
  
  ;-----------------------------------------------------------
  ;                       Guassian
  ;-----------------------------------------------------------
  If not(keyword_set(gauss_steps)) then gauss_steps = 10000
  
  count_3a  = dblarr(gauss_steps)    & count_3a2 =   dblarr(gauss_steps)
  count_3c  = dblarr(gauss_steps)    &count_3e  = dblarr(gauss_steps)    
  Power_3b  = dblarr(gauss_steps)   &Power_3b2  = dblarr(gauss_steps)
  Power_3f  = dblarr(gauss_steps)    &Power_3d  =  dblarr(gauss_steps)
  
  for i=0,gauss_steps-1 do begin
    gauss   = RANDOMN(seed,array_size)
    gauss2  = RANDOMN(seed,array_size)
    index_Ag     = where(gauss ge threshold,count1)    ;number of peaks above threshold
    index_Ag_2   = where(gauss2 ge threshold,count1_2) ;# of peaks above threshold (gaussian 2)
    index_Bg     = where(gauss LE -threshold, count2)   ;number of peaks below threshold
    index_bg_2   = where(gauss2 LE -threshold, count2_2)   ;number of peaks below threshold (gaussian 2)
    index_Q3e    = WHERE(gauss ge threshold and gauss ge 2.0d*gauss2, count3)  ;number of peaks above threshold and gt than off values by factor of  2
     
    count_3a(i)  = count1   ;3a
    count_3a2(i) = count1_2 ;3a (2) 
    
    count_3c(i)  = count2 ;3c
    count_3e(i)  = count3 ;3e

    If count1 gt 0 then   Power_3b(i)  = total(gauss[index_Ag],/double)       ;Guassian Power of peaks above threshold
    If count1_2 gt 0 then Power_3b2(i) = total(gauss2[index_Ag_2],/double) ;Guassian Power of peaks above threshold (Gauss 2)
    
    If count2 gt 0 then begin
      Power_3d(i) =  total(abs(gauss[index_Bg]),/double )
    endif
    If count3 gt 0  then begin
      Power_3f(i) = total(gauss[index_Q3e],/double)      ;power of 3f (Guassian)
    Endif
  endfor ;gauss steps = 10000

  Q3a              = mean(count_3a)                       ;3a
  Q3a_diff         = stddev(count_3a - count_3a2)         ;3a diff
  Q3b              = mean(Power_3b)                       ;3b
  Q3b_diff         = stddev(Power_3b - Power_3b2)         ;3b_diff
  Q3c              = stddev(count_3a - count_3c)  ;3c
  Q3d              = stddev(Power_3b - Power_3d)  ;3d
  Q3e              = mean(count_3e)                       ;3e
  Q3f              = mean(Power_3f)                       ;3f
end