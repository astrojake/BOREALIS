;Authors: Jake Turner (Cornell), Jean-Mathias G., Phillipe Zarka
;
;Purpose: Remove bad frequencies at the edges of subbands in LOFAR/NenuFAR data 
;
;Inputs:   x         ;data(t,f) 
;          f         ;Frequency array 
;          NSUBBAND  ;Number of Subbands 
;          NCHANNEL  ;Number of Channels
;
;Outputs   xFix     ;data with bad frequencies removed 
;          fFix     ;freq array with bad frequencies removed     
;          
;Flag:     xBAD     ;How many bad freq to remove around the subband edge (default = 1) 
;
;Version 1; March 21, 2019                           
pro remove_badf,x,f,NSUBAND,NCHANNEL,xFix=xFix,fFix=Ffix,XBAD=XBAD
  If not(keyword_set(XBAD)) then XBAD = 1           ;remove one freq element 
  nt = n_elements(x[*,0])
  
  x=reform(x,nt,NCHANNEL,NSUBAND )
  x =x(*,XBAD:NCHANNEL-(XBAD+1),*)
  xFix=reform(temporary(x),nt,(NCHANNEL-(XBAD*2))*NSUBAND )
  f=reform(temporary(f),NCHANNEL,NSUBAND )
  f=f(XBAD:NCHANNEL-(XBAD+1),*)
  fFix=reform(temporary(f),(NCHANNEL-(XBAD*2))*NSUBAND)
  return
end