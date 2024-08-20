function DT_DISP, f1,f2,dm

; f1,f2 = band or channel limits in MHz
; dm = dispersion measure [pc.cm-3]
; result is in seconds

dt=4150.*dm*(float(f1)^(-2.) - float(f2)^(-2.))
return,dt
end