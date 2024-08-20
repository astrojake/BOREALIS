pro SET_PS, file, portrait=portrait, landscape=landscape

set_plot,'PS'
if keyword_set(landscape) then device, yoff=27., xsize=27., xoff=0.8, ysize=19.5, /landscape
if keyword_set(portrait) then device, yoff=0., ysize=27., xoff=1.2, xsize=19.5, /portrait

if n_elements(file) ne 0 then device, filename=file
device, bits=8
!p.font=0
!p.charsize=1.3

return
end