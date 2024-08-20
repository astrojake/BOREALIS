; *******************
; *                 *
; *   DELPATH.PRO   *
; *                 *
; *******************

Function DELPATH, str

; removes the path in filename str

p = max([STRPOS(str,'/', /REVERSE_SEARCH),STRPOS(str,':', /REVERSE_SEARCH),STRPOS(str,'\', /REVERSE_SEARCH)])
return,strmid(str,p+1,strlen(str)-p)
end
