;-----------------------------------
  pro REDUCE_ARRAY, x, n, y, dbl=dbl
;-----------------------------------
; x = data of any dimension up to 4
; n = reduction factor in each dimension
; y = reduced array (if reduced dimension is 1, it is removed by reform)

s=size(x)
ndim=s(0)
nn=lonarr(ndim)
for i=0,ndim-1 do nn(i)=floor(s(i+1)/n(i))
if keyword_set(dbl) then begin
  if ndim eq 1 then y=reform( rebin( double( x[0:nn(0)*n(0)-1] ) ,nn(0) ) )
  if ndim eq 2 then y=reform( rebin( double( x[0:nn(0)*n(0)-1,0:nn(1)*n(1)-1] ) ,nn(0),nn(1) ) )
  if ndim eq 3 then y=reform( rebin( double( x[0:nn(0)*n(0)-1,0:nn(1)*n(1)-1,0:nn(2)*n(2)-1] ) ,nn(0),nn(1),nn(2) ) )
  if ndim eq 4 then y=reform( rebin( double( x[0:nn(0)*n(0)-1,0:nn(1)*n(1)-1,0:nn(2)*n(2)-1,0:nn(3)*n(3)-1] ) ,nn(0),nn(1),nn(2),nn(3) ) )
endif else begin
  if ndim eq 1 then y=reform( rebin( float( x[0:nn(0)*n(0)-1] ) ,nn(0) ) )
  if ndim eq 2 then y=reform( rebin( float( x[0:nn(0)*n(0)-1,0:nn(1)*n(1)-1] ) ,nn(0),nn(1) ) )
  if ndim eq 3 then y=reform( rebin( float( x[0:nn(0)*n(0)-1,0:nn(1)*n(1)-1,0:nn(2)*n(2)-1] ) ,nn(0),nn(1),nn(2) ) )
  if ndim eq 4 then y=reform( rebin( float( x[0:nn(0)*n(0)-1,0:nn(1)*n(1)-1,0:nn(2)*n(2)-1,0:nn(3)*n(3)-1] ) ,nn(0),nn(1),nn(2),nn(3) ) )
endelse
return
end
