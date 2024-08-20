;Beam distances
; calc the distance between beams 
;ON: List of ON ra and dec
;OFF; list of OFF ra and dec
;
; TNOLIST: Target is not a list
; ONOLIST; OFF is not a list
;calc_dis,TARGET_LIST,OFF_LIST,DISTANCE=DISTANCE,/OUT,/TNOLIST,Target_RA=Target_RA,Target_Dec=Target_Dec
pro calc_dis,TARGET_LIST,OFF_LIST,DISTANCE=DISTANCE,OUT=OUT,$
  TNOLIST=TNOLIST,Target_RA=Target_RA,Target_Dec=Target_Dec,$
  ONOLIST=ONOLIST,OFF_RA=OFF_RA,OFF_DEC=OFF_DEC

  If not(keyword_set(TNOLIST)) then begin
    readcol,TARGET_LIST,RA_string_t,Dec_string_t,format='(A,A)',DELIMITER=','
  Endif ELSE begin
     RA_string_t  = Target_RA
     Dec_string_t = Target_Dec
  ENDELSE
  
  If not(keyword_set(ONOLIST)) then begin
    readcol,OFF_LIST,RA_string,Dec_string,format='(A,A)',DELIMITER=','
  Endif ELSE begin
    RA_string  = OFF_RA
    Dec_string = OFF_Dec
  ENDELSE
  n= n_elements(RA_String)
  RA = dblarr(n)
  Dec = dblarr(n)
  DISTANCE = dblarr(n)
  RA_Target = ten(RA_string_t)
  DEC_Target =ten(Dec_string_t)

  number = indgen(n,start=1)
  ;Convert RA and Dec and find distance
  for i=0,n-1 do begin
    RA[i] = ten(RA_string[i])     ;(hours)
    DEC[i] = ten(Dec_string[i])   ;(deg)
    GCIRC, 1, RA_Target, DEC_Target, RA[i] , DEC[i], DIS   ;RA in hours, Dec in deg, DIS in arc seconds
    DISTANCE[i] = DIS/3600.0d                              ; degrees
  endfor

  If keyword_set(OUT) then begin
    print_array = dblarr(4,n)
    print_array[0,*] = number[SORT(DISTANCE)]
    print_array[1,*] = RA[SORT(DISTANCE)]
    print_array[2,*] = DEC[SORT(DISTANCE)]
    print_array[3,*] = DISTANCE[SORT(DISTANCE)]
    openw,1,'Distance.dat'
    printf,1,print_array
    close,1
  Endif


end