;------------------------------------------------------------------------
  pro BITARRAY_TO_BYTEARRAY, xbit, xsize, xbyte
;------------------------------------------------------------------------
; xbit = BIT array stored in byte format (per 8 bits)
; xsize = size of original byte array

  ns=n_elements(xsize)
  nx=xsize(ns-1)
  nx8=ceil(nx/8.d0,/L64)*long64(8)
  xbyte=reform(transpose(reform([(xbit and 128b)/128b,(xbit and 64b)/64b,(xbit and 32b)/32b,(xbit and 16b)/16b,(xbit and 8b)/8b,(xbit and 4b)/4b,(xbit and 2b)/2b,xbit and 1b],nx8/8,8)),nx8)
  ndim=xsize(0)
  if ndim eq 1 then xbyte=reform(xbyte(0:nx-1),xsize(1))
  if ndim eq 2 then xbyte=reform(xbyte(0:nx-1),xsize(1),xsize(2))
  if ndim eq 3 then xbyte=reform(xbyte(0:nx-1),xsize(1),xsize(2),xsize(3))
  if ndim eq 4 then xbyte=reform(xbyte(0:nx-1),xsize(1),xsize(2),xsize(3),xsize(4))

return
end

