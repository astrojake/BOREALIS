; -------------------------------------------------
  Function AMJ_AJ, amj
; -------------------------------------------------
; date conversion AAAAMMJJ -> AAAAJJJ or AAMMJJ -> AAJJJ
;		  YYYYMMDD -> YYYYDDD or YYMMDD -> YYDDD
; data type = long integer or double precision, scalar or 1D array
; call : aj = AMJ_AJ(amj)
;	 yd = AMJ_AJ(ymd)

  mois=[0L,31,59,90,120,151,181,212,243,273,304,334,365]
  a=long(amj/10000)
  m=long((amj-a*10000L)/100)
  j=mois(m-1)
  test=float(a)/4.
  for i=0L, n_elements(a)-1L do $
    if test(i) eq float(fix(test(i))) and m(i) ge 3 then j(i)=j(i)+1
  aj=a*1000+j+(amj-a*10000-m*100)

return, aj
end
