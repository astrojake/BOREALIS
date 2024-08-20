;Concatenate pointer array
;Input: Pointer (points to 1d or 2d arrays)
;
; Outputs: array (1d or 2d output) 
pro pointer_to_array,pointer,out_array=array 
 
  steps = n_elements(pointer)
  for j=0,steps-1 do begin ;n dates
    a                = temporary(*pointer[j])
    If j eq 0 then begin
      array  = temporary(a(*,*))
    endif

    If j gt 0 then begin
      temp = [array[*,*],temporary(a[*,*])]
      array = temporary(temp)
    Endif
  endfor ;end Concatenate arrays

end

