;Input
;      x       ; x ramp 
;      Q       ; Q to be plotted
;      Q_l     ; Q label (e.g. Q4a) 
;flags:       
;      Q_sky   ; Q of the sky (note if not dividing by the sky then this value is 1)
;      skyratio; whether Q value is divided by the sky (If no, Q_Sky = 1)
;      unit    ; plot unit for the output file
;
;program is called in plot_textfile.pro
pro q4_textplot,x,Q,Q_l,Q_sky=Q_sky,skyratio=skyratio,print_array=print_array,unit=unit,PROB=PROB,COMBINEHISTO=COMBINEHISTO,Psky=sky
 
 Qlabel = 'Q4'+Q_l
 If keyword_set(COMBINEHISTO) then begin 
  !P.Multi = [0]
 Endif 
;-----Sky Ratio Setup-----------
If keyword_set(skyratio) then ratio='Sky Ratio' else ratio=''
If not(keyword_set(skyratio)) then begin 
   Q_sky = Q 
   Q_sky[*,*] = 1.0d
Endif 

;Setup
n_values = 12   ;num values for Q4i w/o PROB
If keyword_set(PROB) then n_values=13

;-------Master Array---------- 
y = Q_sky    
y[*,*] = Q[*,*]/(Q_sky[*,*] + (Q_sky[*,*] eq 0))

for i =0,n_values-1 do begin
  If i eq 0 then title  = 'Q4'+Q_l+': ON-OFF Integral'
  If i eq 1 then title  = 'Q4'+Q_l+': Integral Ratio to Gaussian (1 sig)'
  If i eq 4 then title  = 'Q4'+Q_l+': ON-OFF > 0'
  If i eq 5 then title  = 'Q4'+Q_l+': ON-OFF > 1 sig Gaussian'
  If i eq 10 then title  = 'Q4'+Q_l+': ON-OFF % between max and min'
  If i eq 11 then title  = 'Q4'+Q_l+': ON-OFF % > max'
  
  If i eq 0 then ytitle = 'Integral (SEFD)'
  If i eq 1 then ytitle = 'Ratio'
  If i eq 4 then ytitle = '%'
  If i eq 5 then ytitle = 'Number'
  If i eq 10 then ytitle = '%'
  If i eq 11 then ytitle = 'Probability'
  If keyword_set(skyratio) then ytitle = 'Sky Ratio'

  y[*,1] = Q[*,1]  ;No ratio for integral ratio
  If i eq 1 then ytitle = 'Ratio'
  
  ;skip these i values (don't plot)
  If i eq 2 or i eq 3 or i eq 6 or i eq 7 or i eq 8 or i eq 9 then begin 
    GOTO,skip
  Endif  
  
  ;plot
  If not(keyword_set(COMBINEHISTO)) then begin  
    cgplot,x,y[*,i],/xstyle, /ystyle,xrange=[min(x)-1,max(x)+1],yrange=[min(y[*,i])-max(y[*,i])*0.1,max(y[*,i])+max(y[*,i])*0.1],$
    xtitle='Run #',ytitle=ytitle,title=title,psym=16
  endif 
  
  If keyword_set(COMBINEHISTO) then begin 
;    !P.Multi = [0]
;    pos = gang_plot_pos(1,2,0)
;    pos2 = gang_plot_pos(1,2,1)
;;     If i eq 0 then  begin 
;      cgplot,x,y[*,i],/xstyle, /ystyle,xrange=[min(x)-1,max(x)+1],yrange=[min(y[*,i])-max(y[*,i])*0.1,max(y[*,i])+max(y[*,i])*0.1],$
;        xtitle='Run #',ytitle=ytitle,title=title,psym=16,position=pos.position      
;     cgHistoplot, y[*,i],NBINS=10,xtitle='Number',position=pos2.position,/rotate,$
;        yrange=[min(y[*,i])-max(y[*,i])*0.1,max(y[*,i])+max(y[*,i])*0.1],YTickformat='(A1)',thick=3,/FILL
;      AXIS, YAXIS = 1
;    If i eq 1 then cgHistoplot, y[*,1],title='Q4'+Q_l+': Ratio to Gaussian (1 sig)',NBINS=10,xtitle='Value',_EXTRA=gang_plot_pos(1,2,1)
;    If i eq 4 then cgHistoplot, y[*,4],title='Q4'+Q_l+': ON-OFF > 0',NBINS=10,xtitle='Value',_EXTRA=gang_plot_pos(1,2,1)
;    If i eq 5 then cgHistoplot, y[*,5],title='Q4'+Q_l+': ON-OFF Values > 1 sig',NBINS=10,xtitle='Value',_EXTRA=gang_plot_pos(1,2,1)
;    If i eq 10 then cgHistoplot, y[*,10],title='Q4'+Q_l+': % between max and min',NBINS=10,xtitle='Value',_EXTRA=gang_plot_pos(1,2,1)
;    If i eq 11 then cgHistoplot, y[*,11],title='Q4'+Q_l+': % > max',NBINS=10,xtitle='Value',_EXTRA=gang_plot_pos(1,2,1)
;   
;   Endif    
     !P.Multi = [0,1,1]
     If i eq 0 then title  = Qlabel+': N5'   ;ON-OFF Integral'
     If i eq 1 then title  = Qlabel+': N6'   ; Integral Ratio to Gaussian (1 sig)'
     If i eq 4 then title  = Qlabel+': N7'   ;ON-OFF > 0'
     If i eq 5 then title  = Qlabel+': N8'   ;ON-OFF > 1 sig Gaussian'
     If i eq 10 then title  = Qlabel+': N10' ;  ON-OFF % between max and min'
     If i eq 11 then title  = Qlabel+': N12'  ;  ON-OFF % > max'
     
     If keyword_set(sky) then Ymin = [min(y[*,i]),min(sky[*,i])] else Ymin = min(y[*,i])
     If keyword_set(sky) then  ymax = [max(y[*,i]),max(sky[*,i])] else ymax = max(y[*,i])
    
     cgplot,x,y[*,i],/xstyle, /ystyle,xrange=[min(x)-1,max(x)+1],yrange=[min(y[*,i])-min(y[*,i])*0.20,max(y[*,i])+max(y[*,i])*0.20],$
       xtitle='Run #',ytitle=ytitle,title=title,psym=16;,position=pos.position
     xyouts,66,0.9,'66'
     xyouts,129,1.4,'129'  
     cgplot,[45,45],[-1e6,1e6],/overplot,linestyle=2
      cgplot,[126,126],[-1e6,1e6],/overplot,linestyle=2
     If keyword_set(sky) then begin
       cgplot,x,sky[*,i],/overplot,psym=15,color='red'
       If i eq 0 then al_legend,['ON - OFF', 'OFF 1 - OFF 2 '],textcolor=['black','red'],psym=[16,15],colors=['black','red'],/top,/left,background_color='white'
       If i eq 1 then al_legend,['ON - OFF Ratio ', 'OFF 1 - OFF 2 Ratio'],textcolor=['black','red'],psym=[16,15],colors=['black','red'],/top,/left,background_color='white'
       If i eq 4 then al_legend,['ON - OFF > 0  ', 'OFF 1 - OFF 2 >0 '],textcolor=['black','red'],psym=[16,15],colors=['black','red'],/top,/left ,background_color='white'
       If i eq 5 then al_legend,['ON - OFF > 1'+cgGreek('sigma')+'Gauss', 'OFF 1 - OFF 2 > 1'+cgGreek('sigma')+'Gauss'],textcolor=['black','red'],$
        psym=[16,15],colors=['black','red'],/top,/left,background_color='white'
       
       ratio = y[*,i]/sky[*,i]
      
       cgplot,x,ratio,/xstyle, /ystyle,xrange=[min(x)-1,max(x)+1],$
         xtitle='Run #',ytitle='Ratio: ON/Sky',title=title,psym=16;,position=pos.position
      endif  
   !P.Multi = [0, 3, 2]
   endif 
  ;---------
  ;integral
  ;---------
  If i eq 0 and keyword_set(skyratio) then begin
    printf,unit,'**********Q4'+Q_l+' Detection (Integral Ratio > 1)********* '
    w = where(y[*,i] ge 1.00 and Q[*,0] ge 0,count) ;sky ratio ge 1 and value is postive
    If count ge 1 then printf,unit,print_array[*,w]
  endif
  
 ;---------------------------
 ; ON-OFF Ratio with 1 sig (i = 1), 2 sig (i=2), 3 sig (i=3)
 ;----------------------------
 If i eq 1 then begin  
   w = where(Q[*,i] ge 1.0,count)   ; value is postive
   If count ge 1 then begin
     printf,unit,'**********Q4'+Q_l+' Detection (Gaussian Ratio)********* '
     printf,unit,print_array[*,w]
   Endif   ;If count ge 1 then print,'Q[*,1]',Q[w,1]
   If keyword_set(skyratio) then begin
     printf,unit,'**********Q4'+Q_l+' Detection (Gaussian Ratio with Sky Ratio)********* '
     w = where(y[*,i] ge 1 and Q[*,i] ge 1,count) ; value is postive
     printf,unit,'1sigma'
     If count ge 1 then printf,unit,print_array[*,w]
   endif
 endif
 
 ;------------------------------------------
 ;--------  ON-OFF values > 0        ------- (i eq 4) 
 ;------------------------------------------
 If i eq 4 and keyword_set(skyratio) then begin
    printf,unit,'**********Q4'+Q_l+' Detection (On-Off>0)********* '
    w = where(y[*,i] ge 1 and Q[*,4] ge 0,count) ;sky ratio ge 1 and value is postive
    If count ge 1 then printf,unit,print_array[*,w]
 endif
 
 ;------------------------------------------
 ;   ON-OFF values above 1,2,3 sigma 
 ;   (1sig = i eq 5)(2sig = i eq 6)(3sig = i eq 7)
 ;------------------------------------------
 If keyword_set(skyratio) then begin
   printf,unit,'**********Q4'+Q_l+' Detection (Values above 1,2,3sigma)********* '
   w = where(y[*,5] ge 1 and Q[*,5] ge 1,count) ; value is postive
   If count ge 1 then printf,unit,print_array[*,w]
 endif
 
 ;---------------------------------------------------------------------
 ;---% of ON-OFF values between min and max threshold & >0
 ;   min (i eq 8) max (i eq 9) % i eq 10
 ;-------------------------------------------------------------------
 If i eq 10 and keyword_set(skyratio) then begin
   printf,unit,'**********Q4'+Q_l+' Detection (% between min and max)********* '
   w = where(y[*,10] ge 1 and Q[*,10] ge 0,count) ;sky ratio ge 1 and value is postive
   If count ge 1 then printf,unit,print_array[*,w]
 endif
 
 ;--------------------------------------------------
 ;      % of ON-OFF values lt 0 and > max threshold
 ;----------------------------------------------------
 If i eq 11 and keyword_set(skyratio) then begin
   printf,unit,'**********Q4'+Q_l+' Detection (%>max; Sky Ratio)********* '
   w = where(y[*,11] lt 1,count);sky ratio ge 1
   If count ge 1 then printf,unit,print_array[*,w]
 endif
  
 ;---------------------------------------------------
 ;       Probability of False Detection
 ;---------------------------------------------------
 If i eq 12 and keyword_set(PROB) then begin 
   If not(keyword_set(skyratio)) then begin
     w = where(y[*,12] lt 2e-3,count)       ;False Detectio lt 2e-3 (3.0 sigma)
     If count ge 1 then begin
       printf,unit,'**********Q4'+Q_l+' Prob of False Detection (3.0 sigma)********* '
       printf,unit,print_array[*,w]
     endif
   endif
   If keyword_set(skyratio) then begin
     w = where(Q[*,12] lt 2e-3 and y[*,12] le 1,count)       ;False Detectio lt 2e-3 (3.0 sigma)
     If count ge 1 then begin
       printf,unit,'**********Q4'+Q_l+' Prob of False Detection Ratio********* '
       printf,unit,print_array[*,w]
     endif
   endif
 endif ;eq 11
 SKIP:
 
endfor

;Histograms
If not(keyword_set(COMBINEHISTO)) then begin 
  cgHistoplot, y[*,0],title='Q4'+Q_l+': ON-OFF Integral',NBINS=10,xtitle='Value'
  cgHistoplot, y[*,1],title='Q4'+Q_l+': Ratio to Gaussian (1 sig)',NBINS=10,xtitle='Value'
  cgHistoplot, y[*,4],title='Q4'+Q_l+': ON-OFF > 0',NBINS=10,xtitle='Value'
  cgHistoplot, y[*,5],title='Q4'+Q_l+': ON-OFF Values > 1 sig',NBINS=10,xtitle='Value'
  cgHistoplot, y[*,10],title='Q4'+Q_l+': % between max and min',NBINS=10,xtitle='Value'
  cgHistoplot, y[*,11],title='Q4'+Q_l+': % > max',NBINS=10,xtitle='Value'
  If keyword_set(PROB) then cgHistoplot, y[*,12],title='Q4'+Q_l+'Prob of False',NBINS=10,xtitle='Value'
endif 
 
; ;------------------------------------------
; ;On-OFF Integral Ratio to Gaussian (1 sig)
; ;------------------------------------------
; y = Q[*,1]/(Q_Sky[*,1] +  (Q_Sky[*,1]  eq 0))
;; y2 = Q[*,2]/Q_sky[*,2]
;; y3 = Q[*,3]/Q_sky[*,3]
; cgplot,x,y,/xstyle, /ystyle,xrange=[min(x)-1,max(x)+1],yrange=[min(y)-max(y)*0.1,max(y)+max(y)*0.1], $
;        xtitle='Run #',ytitle='Ratio',title='Q4'+Q_l+': On-OFF Integral Ratio to Gaussian (1 sig)',psym=16
;  ;---print---
;   w = where(Q[*,1] ge 0.90,count)   ; value is postive 
;   If count ge 1 then begin 
;     printf,unit,'**********Q4'+Q_l+' Detection (Gaussian Ratio)********* '
;     printf,unit,print_array[*,w]
;   Endif   ;If count ge 1 then print,'Q[*,1]',Q[w,1]
;If keyword_set(skyratio) then begin
;   printf,unit,'**********Q4'+Q_l+' Detection (Gaussian Ratio with Sky Ratio)********* '
;   w = where(y1 ge 1 and Q[*,1] ge 1,count) ; value is postive 
;   printf,unit,'1sigma'
;   If count ge 1 then printf,unit,print_array[*,w]
;endif 
;
; ;------------------------------------------
; ;--------  ON-OFF values > 0        -------
; ;------------------------------------------
; y = Q[*,4]/(Q_sky[*,4] + (Q_sky[*,4] eq 0))
; cgplot,x,y,/xstyle, /ystyle,xrange=[min(x)-1,max(x)+1],yrange=[min(y)-max(y)*0.1,max(y)+max(y)*0.1], $
;  xtitle='Run #',ytitle='%'+' '+ratio,title='Q4'+Q_l+': ON-OFF values > 0',psym=16
;  ;---print---
;  If keyword_set(skyratio) then begin
;   printf,unit,'**********Q4'+Q_l+' Detection (On-Off>0)********* '
;   w = where(y ge 1 and Q[*,4] ge 0,count) ;sky ratio ge 1 and value is postive 
;   If count ge 1 then printf,unit,print_array[*,w]
;  endif 
;  
;------------------------------------------
;   ON-OFF values above 1,2,3 sigma
;------------------------------------------
; Y1 = Q[*,5]/(Q_SKY[*,5] + (Q_SKY[*,5] EQ 0))
; Y2 = Q[*,6]/(Q_SKY[*,6] + (Q_SKY[*,6] EQ 0))
; Y3 = Q[*,7]/(Q_SKY[*,7] + (Q_SKY[*,7] EQ 0))
;YS = [Y1,Y2,Y3]
; CGPLOT,X,Y1,/XSTYLE, /YSTYLE,XRANGE=[MIN(X)-1,MAX(X)+1],YRANGE=[MIN(YS)-MAX(YS)*0.1,MAX(YS)+MAX(YS)*0.1],$
;        XTITLE='RUN #',YTITLE='NUMBER'+' '+RATIO,TITLE='Q4'+Q_L+': ON-OFF VALUES > 1 SIGMA',PSYM=16,COLOR='BLACK'
;If keyword_set(skyratio) then begin
;   printf,unit,'**********Q4'+Q_l+' Detection (Values above 1,2,3sigma)********* '
;   w = where(y1 ge 1 and Q[*,5] ge 1,count) ; value is postive 
;   printf,unit,'1sigma'
;  If count ge 1 then  printf,unit,print_array[*,w]
;   w = where(y2 ge 1 and Q[*,6] ge 1,count) ; value is postive 
;   printf,unit,'2sigma'
;   If count ge 1 then printf,unit,print_array[*,w]
;   w = where(y3 ge 1 and Q[*,7] ge 1,count) ;value is postive 
;   printf,unit,'3sigma'
;  If count ge 1 then printf,unit,print_array[*,w]
;endif 

;---------------------------------------------------------------------
;---% of ON-OFF values between min and max threshold & >0        
;-------------------------------------------------------------------
; y = Q[*,10]/(Q_sky[*,10]  + (Q_sky[*,10] eq 0))
; cgplot,x,y,/xstyle, /ystyle,xrange=[min(x)-1,max(x)+1],yrange=[min(y)-max(y)*0.1,max(y)+max(y)*0.1],$
;        xtitle='Run #',ytitle='% '+ratio,title='Q4'+Q_l+': % between max and min',psym=16
; ;---print---
;If keyword_set(skyratio) then begin
;   printf,unit,'**********Q4'+Q_l+' Detection (% between min and max)********* '
;   w = where(y ge 1 and Q[*,10] ge 0,count) ;sky ratio ge 1 and value is postive 
;  If count ge 1 then printf,unit,print_array[*,w]
;endif 

;;--------------------------------------------------
;;      % of ON-OFF values lt 0 and > max threshold
;;----------------------------------------------------
; y = Q[*,11]/(Q_Sky[*,11] + (Q_sky[*,11] eq 0))
; finite = finite(y)
; y = y*finite
; cgplot,x,y,/xstyle, /ystyle,xrange=[min(x)-1,max(x)+1],yrange=[min(y)-max(y)*0.1,max(y)+max(y)*0.1],$
;   xtitle='Run #',ytitle='% '+ratio,title='Q4'+Q_l+': % > max',psym=16
;; ;---print---
;If keyword_set(skyratio) then begin
;   printf,unit,'**********Q4'+Q_l+' Detection (%>max; Sky Ratio)********* '
;   w = where(y lt 1,count);sky ratio ge 1 
;  If count ge 1 then printf,unit,print_array[*,w]
;endif

;;---------------------------------------------------
;;       PROBABILITY OF FALSE DETECTION
;;---------------------------------------------------
;IF KEYWORD_SET(PROB) THEN BEGIN
;  Y = Q[*,12]/(Q_SKY[*,12] + (Q_SKY[*,12] EQ 0))
;  CGPLOT,X,Y,/XSTYLE, /YSTYLE,XRANGE=[MIN(X)-1,MAX(X)+1],YRANGE=[MIN(Y)-MAX(Y)*0.1,MAX(Y)+MAX(Y)*0.1],$
;    XTITLE='RUN #',YTITLE='PROBABILITY '+RATIO,TITLE='Q4'+Q_L+': PROB OF FALSE DETECTION',PSYM=16,/YLOG
;  IF NOT(KEYWORD_SET(SKYRATIO)) THEN BEGIN
;    W = WHERE(Y LT 2E-3,COUNT)       ;FALSE DETECTIO LT 2E-3 (3.0 SIGMA)
;    IF COUNT GE 1 THEN BEGIN 
;      PRINTF,UNIT,'**********Q4'+Q_L+' PROB OF FALSE DETECTION (3.0 SIGMA)********* '
;      PRINTF,UNIT,PRINT_ARRAY[*,W]
;    ENDIF
;  ENDIF 
;  IF KEYWORD_SET(SKYRATIO) THEN BEGIN
;     W = WHERE(Q[*,12] LT 2E-3 AND Y LE 1,COUNT)       ;FALSE DETECTIO LT 2E-3 (3.0 SIGMA)
;     IF COUNT GE 1 THEN BEGIN 
;      PRINTF,UNIT,'**********Q4'+Q_L+' PROB OF FALSE DETECTION RATIO********* '
;      PRINTF,UNIT,PRINT_ARRAY[*,W]
;     ENDIF  
;  ENDIF
;ENDIF 

end