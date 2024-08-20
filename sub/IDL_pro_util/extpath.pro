; *******************
; *                 *
; *   EXTPATH.PRO   *
; *                 *
; *******************

Function EXTPATH, str

; retrieves the path of filename str

p = max([strpos(str,'/',/REVERSE_SEARCH),strpos(str,':',/REVERSE_SEARCH),strpos(str,'\',/REVERSE_SEARCH)])
return,strmid(str,0,p+1)
end
