;------------------------------------------------------------------------
  pro BYTEARRAY_TO_BITARRAY, xbyte, xsize, xbit
;------------------------------------------------------------------------
; xbyte = BYTE array 0 or 1 (or e.g. do x = x gt 0) to convert to bit array
; xsize = size of original byte array
  
  xsize=size(xbyte)
  ns=n_elements(xsize)
  nx=xsize(ns-1)
  nx8=ceil(nx/8.d0,/L64)*long64(8)
  if nx8 eq nx then xbit=reform(xbyte,8,nx8/8) else xbit=reform([reform(xbyte,nx),bytarr(nx8-nx)],8,nx8/8)
  xbit=reform(byte([128,64,32,16,8,4,2,1]#xbit))

return
end

