;program for Jupiter Setup
;This program is not general
;
pro lofar_area,ANTENNA_SET,NOF_Stations,xf,A_eff_L=A_eff_L
  If ANTENNA_SET eq 'LBA_INNER' then begin
    print,'-----Running LBA_INNER----'
    ;------------------------
    ;Calibration info
    ;-----------------------
    lamda   =  [20.0,10.0,6.67,5.0,4.0]                  ;m
    A_in_meas =  [419.77,419.77,415.37,347.37,239.67]  ;m^2

    ;fit with two functions
    poly      = poly_fit(lamda[2:4],A_in_meas[2:4],2)
    poly2     = poly_fit(lamda[0:2],A_in_meas[0:2],1)

    ;Lofar Wavelengths
    x_wave_Lofar = 2.998e8/(xf*1.0d6)

    ;Lofar Area
    A_inner = (poly[0] + poly[1]*x_wave_Lofar + poly[2]*x_wave_Lofar^2.)
    w = where(x_wave_Lofar gt 6.67)
    A_inner(w) =  (poly2[0] + poly2[1]*x_wave_Lofar(w))
    A_eff_L = A_inner*NOF_stations
  Endif

  If ANTENNA_SET eq 'LBA_OUTER' then begin
    ;------------------------
    ;Calibration info
    ;-----------------------
    lamda   =  [20.,10.,6.67,5.0,4.0]                  ;m
    A_out_meas =  [1973.4,1343.5,693.61,398.18,256.00] ;m^2

    ;----------
    ;fit
    ;----------
    poly      = poly_fit(lamda,A_out_meas,2)

    ;Lofar Wavelengths
    x_wave_Lofar = 2.998e8/(xf*1.0d6)

    A_outer = (poly[0] + poly[1]*x_wave_Lofar + poly[2]*x_wave_Lofar^2. ) ;km
    A_eff_L = A_outer*NOF_stations ;km
  Endif
return  
end