;This program is a subprogram for lofar_postprocessing
;
;Purpose: Calculate Q's for an input of two variables (y and y2) 
;         This program can be used in a loop or not
;
;Input;       y
;             y2        
;             threshold
;
pro postprocessing_calc,y,y2,threshold,slope=slope,$
                        NxA=NxA,PNxA=Power_NxA,$
                           Nx2A=Nx2A,PNx2A=Power_Nx2A,$
                           Nx3Ae=NxA_3e,PNAx3f=Power_NxA_3f,$
                           Nx3Be=NxB_3e,PNBx3f=Power_NxB_3f,$
                           NxB=NxB,PNxB=Power_NxB,$
                           Nx2B=Nx2B,PN2xB=Power_Nx2B
  
  ;create arrays
  NxA           = 0. & Nx2A         = 0.
  NxA_3e        = 0. & NxB_3e       = 0. 
  Power_NxA     = 0. & Power_Nx2A   = 0. 
  Power_NxA_3f  = 0. & Power_NxB_3f = 0.
  NxB           = 0. & Nx2B         = 0.  
  Power_NxB     = 0. & Power_Nx2B   = 0. 
  
  If not(keyword_set(slope)) then slope  = 2.0 
  If slope ne 2.0 then begin
    print, '----Slope Value Changed-----'
    print, 'Slope Value: ',slope
  Endif
  
  ;---------------
  ;Count peaks
  ;---------------
  index_NxA    = where(y GE threshold,count)                         ;above threshold ON
  NxA          = count                                               ;number of peaks for ON beam above theshold
  index_Nx2A   = where(y2 GE threshold,count2)                       ;above threshold OFF
  Nx2A         = count2                                              ;number of peaks for OFF beam above theshold
  index_Q3e_A  = WHERE(y ge threshold and y ge slope*y2, count_3e)    ;ON  (Q3e)
  NxA_3e       = count_3e                                            ;# of peaks above threshold for ON beam and gt than 2*OFF values
  index_Q3e_B  = WHERE(y2 ge threshold and y2 ge slope*y, count_3e2)  ;OFF flipping (Q3e)
  NxB_3e       = count_3e2                            ;Flipped! # of peaks above threshold for OFF beam & gt than 2*ON values

  ;---------------------
  ;Find count and Power
  ;---------------------
  If count gt 0 then begin
    Power_NxA   = total(y[index_NxA],/double)           ;power of peaks for on beam
  endif
  If count2 gt 0 then begin
    Power_Nx2A = total(y2[index_Nx2A],/double)          ;power of peaks for off beam
  Endif
  If count_3e gt 0 then begin
    Power_NxA_3f = total(y[index_Q3e_A],/double)        ;Power 3f (ON)
  endif
  If count_3e2 gt 0 then begin
    Power_NxB_3f = total(y2[index_Q3e_B],/double)       ;Power 3f (OFF)
  Endif

  ;-----------------------
  ;Asymmetry Parameter
  ;---------------------- 
  index_NxB   = where(y LE -threshold, count1b)     ;number of peaks for OM beam below theshold
  NxB          = count1b
  index_Nx2B  = where(y2 LE -threshold, count2b)    ;number of peaks for OFF beam below theshold
  Nx2B       = count2b
  
  If count1b gt 0 then begin
    Power_NxB = total(abs(y[index_NxB]),/double )   ;power of peaks for on beam below tau
  endif
  If count2 gt 0 then begin
    Power_Nx2B = total(abs(y2[index_Nx2B]),/double )   ;power of peaks for off beam below tau
  endif

return 
end