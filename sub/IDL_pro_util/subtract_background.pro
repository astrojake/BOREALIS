;--------------------------------------------------------------------
  pro SUBTRACT_BACKGROUND, f,s,seuil, tab
;--------------------------------------------------------------------
; Substraction of (background + seuil * sigmas) from a 2-D array of raw 
; data, and setting of resulting negative values to 0.

; f,s =  (INPUT)  arrays of background and fluctuations (1 sigma level)
; seuil =  (INPUT)  threshold value for the number of sigmas above background
; tab =  (INPUT/OUTPUT)  2-D array of raw data (time,freq) -> tab-(f+seuil*s)

  for i=0,n_elements(tab(*,0))-1 do tab(i,*)=(tab(i,*)-f-seuil*s) >0
return
end
