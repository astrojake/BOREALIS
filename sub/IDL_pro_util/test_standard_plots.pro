set_ps,'sp.ps',/landscape
x=randomn(seed,6000,4000)
t=findgen(6000)
f=findgen(4000)/3999*48+26

set_ps,'sp.ps',/landscape
t = findgen(3999)
f = findgen(61)
label='Test'
restore,filename='83sec_RFI_patrol5.5_pex_lesig3.5_beam1.sav'
STANDARD_PLOTS, x_new_disp,t,f,label
device,/close
ps_pdf,/rem

end

