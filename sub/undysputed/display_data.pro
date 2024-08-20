pro DISPLAY_DATA,x,nx

if n_elements(nx) eq 0 then nx=1400

nt=n_elements(x(*,0))
nf=n_elements(x(0,*))

nnt=ceil(nt*1./nx)
nx=long(nt/nnt)

nnf=ceil(nf/900.)
ny=long(nf/nnf)

xx=rebin(x(0:nx*nnt-1,0:ny*nnf-1),nx*nnt,ny)
window,0,xs=nx,ys=ny

rien=''
for i=0,nx*nnt-1,nx do begin
  print,i,' -> ',i+nx
  tvscl,xx[i:i+nx-1,*]
  read,rien
endfor

return
end
