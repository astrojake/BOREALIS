;***************************************************************
;      Convert Byte or Bit arrays for lofar data 
;**************************************************************
;***************************************************************
;Author: Jake Turner (UVA)
;        J.M. Griessmeier (CNRS) and Philippe Zarka (Obspm)
;        
;Inputs: 
;        p2    ;  mask to be convereted
;        xsize ;  size of orginal mask
;        
;Outputs:       
;        p2_out:  converted mask. 2-D array
;        
;Keywords: BYTE_TO_BIT
;          BIT_TO_BIT
;
pro lofar_bytebit,p2,xsize,p2_out=p2_out,BIT_TO_BYTE=BIT_TO_BYTE,BYTE_TO_BIT=BYTE_TO_BIT
 
  nf =  n_elements(p2(0,*))

  If keyword_set(BIT_TO_BYTE) then begin
    BITARRAY_TO_BYTEARRAY, p2(*,0), xsize, p2_temp
    nt = n_elements(p2_temp)
  endif
  
  If keyword_set(BYTE_TO_BIT) then begin
    BYTEARRAY_TO_BITARRAY, p2(*,0), xsize, p2_temp
    nt = n_elements(p2_temp)
  endif

 p2_out = bytarr(nt,nf)+1b

  for i=0,nf-1 do begin
     If keyword_set(BIT_TO_BYTE) then begin
        BITARRAY_TO_BYTEARRAY, p2(*,i), xsize, p2_temp
        p2_out(*,i) = p2_temp
     endif   
     If keyword_set(BYTE_TO_BIT) then begin
        BYTEARRAY_TO_BITARRAY, p2(*,i), xsize, p2_temp 
        p2_out(*,i) = p2_temp
     endif   
  endfor

return
end