;--------------------------------------------------------------------
  function DYN_N, data, frac
;--------------------------------------------------------------------
; computes values of quantiles in a distribution

; (INPUT)
; data = array of values
; frac[n] = quantiles for which ...
; (OUTPUT)
; ... value[n] values are required

  data_sorted=data(sort(data))
  ndata=n_elements(data)
  n=n_elements(frac)
  value=double(frac)
  for i=0L,n-1 do value(i)=data_sorted((ndata-1.0d)*frac(i))
return,value
end
