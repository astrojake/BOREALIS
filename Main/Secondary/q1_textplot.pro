;This program plots Q1a and Q1b
;
;Input x       ; x ramp
;      y      ; Q to be plotted  (Q[*,4]) ---->> 4 numbers
;
;program is called in plot_textfile.pro
pro q1_textplot,x,y,unit=unit,ytitle=ytitle,title=title,skyratio=skyratio,factor=factor,print_array=print_array,Q1a=Q1a,Q1b=Q1b,PROB=PROB,COMBINEHISTO=COMBINEHISTO,Psky=sky

If not(keyword_set(Psky)) then begin 
  !P.Multi = [0, 3, 2]
  If keyword_set(PROB) then !P.Multi = [0, 4, 2]
Endif
If keyword_set(Q1a) then Qlabel = 'Q1a'
If keyword_set(Q1b) then Qlabel = 'Q1b'
If not(keyword_set(factor)) then factor =2 

n_numbers = 4  ;number of quant #s 
If not(keyword_set(PROB)) then n_numbers=3

for i =0,n_numbers-1 do begin 
  
  If i eq 0 then title  = Qlabel+': ON-OFF Integral'
  If i eq 1 then title  = Qlabel+': ON-OFF >0'
  If i eq 2 then title  = Qlabel+': STDDEV On-OFF'
  If i eq 3 then title  = Qlabel+': ON-OFF Integral'
  If i eq 4 then title  = Qlabel+': PROB of False Detection'
  If not(keyword_set(skyratio)) then begin 
    If i eq 0 then ytitle = 'Integral (SEFD)'
    If i eq 1 then ytitle = '%'
    If i eq 2 then ytitle = 'STDEV (SEFD)'
    If i eq 3 then ytitle = 'Probability'
  Endif  
  
  If not(keyword_set(COMBINEHISTO)) then begin
    print,'Not COMBINEHISTO'
    If keyword_set(sky) then Ymin = [min(y[*,i]),min(sky[*,i])] else Ymin = min(y[*,i])
    If keyword_set(sky) then  ymax = [max(y[*,i]),max(sky[*,i])] else ymax = max(y[*,i])

    cgplot,x,y[*,i],/xstyle, /ystyle,xrange=[min(x)-1,max(x)+1],$
      yrange=[Ymin-Ymin*0.15,ymax+ymax*0.15],$
      xtitle='Run #',ytitle=ytitle,title=title,psym=16
  endif 
  If keyword_set(COMBINEHISTO) then begin
    print,'********COMBINEHISTO********'

     !P.Multi = [0,1,1]
;    pos = gang_plot_pos(1,2,0)
;    pos2 = gang_plot_pos(1,2,1)
    ;     If i eq 0 then  begin
    If i eq 0 then title  = Qlabel+': N1'
    If i eq 1 then title  = Qlabel+': N2'
    If i eq 2 then title  = Qlabel+': N3'
    
   If keyword_set(sky) then Ymin = [min(y[*,i]),min(sky[*,i])] else Ymin = min(y[*,i])
   If keyword_set(sky) then  ymax = [max(y[*,i]),max(sky[*,i])] else ymax = max(y[*,i])
   
    cgplot,x,y[*,i],/xstyle, /ystyle,xrange=[min(x)-1,max(x)+1],yrange=[Ymin-Ymin*0.15,ymax+ymax*0.15],$
      xtitle='Run #',ytitle=ytitle,title=title,psym=16;,position=pos.position     
    If keyword_set(sky) then begin
      print, i
      cgplot,x,sky[*,i],/overplot,psym=15,color='red'
      cgplot,[45,45],[-1e6,1e6],/overplot,linestyle=2
      cgplot,[126,126],[-1e6,1e6],/overplot,linestyle=2
      If i eq 0 then al_legend,['ON - OFF', 'OFF 1 - OFF 2 '],textcolor=['black','red'],psym=[16,15],colors=['black','red'],/top,/left,background_color='white'
      If i eq 1 then al_legend,['ON - OFF > 0 ', 'OFF 1 - OFF 2 >0'],textcolor=['black','red'],psym=[16,15],colors=['black','red'],/top,/left,background_color='white'
      If i eq 2 then al_legend,['ON - OFF', 'OFF 1 - OFF 2 '],textcolor=['black','red'],psym=[16,15],colors=['black','red'],/top,/left,background_color='white'

;      cgHistoplot, y[*,i],NBINS=10,xtitle='Number',position=pos2.position,/rotate,$
;        YTickformat='(A1)',thick=3,/FILL,binsize=binsize,yrange=[min(y[*,i])-min(y[*,i])*0.1,max(y[*,i])+max(y[*,i])*0.1] 
      ;cgHistoplot, sky[*,i],/oplot,POLYCOLOR='royal blue',/rotate,BINSIZE=binsize,yrange=[min(y[*,i])-min(y[*,i])*0.1,max(y[*,i])+max(y[*,i])*0.1]
      ;DELVAR,binsize
     !P.Multi = [0, 3, 2]
    Endif      ;yrange=[min(y[*,i])-min(y[*,i])*0.1,max(y[*,i])+max(y[*,i])*0.1],
      ;  AXIS, YAXIS = 2
    ;    If i eq 1 then cgHistoplot, y[*,1],title='Q4'+Q_l+': Ratio to Gaussian (1 sig)',NBINS=10,xtitle='Value',_EXTRA=gang_plot_pos(1,2,1)
    ;    If i eq 4 then cgHistoplot, y[*,4],title='Q4'+Q_l+': ON-OFF > 0',NBINS=10,xtitle='Value',_EXTRA=gang_plot_pos(1,2,1)
    ;    If i eq 5 then cgHistoplot, y[*,5],title='Q4'+Q_l+': ON-OFF Values > 1 sig',NBINS=10,xtitle='Value',_EXTRA=gang_plot_pos(1,2,1)
    ;    If i eq 10 then cgHistoplot, y[*,10],title='Q4'+Q_l+': % between max and min',NBINS=10,xtitle='Value',_EXTRA=gang_plot_pos(1,2,1)
    ;    If i eq 11 then cgHistoplot, y[*,11],title='Q4'+Q_l+': % > max',NBINS=10,xtitle='Value',_EXTRA=gang_plot_pos(1,2,1)
    ;   Endif
  endif
  
  ;-----------------------
  ;---- Integral----------
  ;-----------------------
  If i eq 0 and keyword_set(skyratio) then begin
    printf,unit,'**********'+Qlabel+' Detection (integral Ratio > 1)********* '
    w = where(y[*,0] ge Factor,count) ;sky ratio ge 1 and value is postive
    If count ge 1 then printf,unit,print_array[*,w]
  endif
  
  ;-----------------------
  ;---(On-Off>0)-----------
  ;-----------------------
  If i eq 1 then begin 
    If not(keyword_set(skyratio)) then begin
      w = where(y[*,1] ge 80,count )     ;% where On-Off>0 gt 80%
      printf,unit,'**********'+Qlabel+' Detection (On-Off>0)********* '
      If count ge 1 then printf,unit,print_array[*,w]
    Endif
    If keyword_set(skyratio) then begin
      printf,unit,'**********'+Qlabel+' Detection (On-Off>0 Ratio)********* '
      w = where(y[*,1] ge Factor,count ) ;sky ratio ge 1 and value is postive
      If count ge 1 then printf,unit,print_array[*,w]
    endif
  endif
  
  ;-----------------------
  ;---STDEV-----------
  ;-----------------------
  If i eq 2 and keyword_set(skyratio) then begin
    printf,unit,'**********'+Qlabel+' Detection (STDDEV Ratio)********* '
    w = where(y[*,2] ge Factor ge 0,count) ;sky ratio ge 1 and value is postive
    If count ge 1 then printf,unit,print_array[*,w]
  endif
  
  ;-----------------------
  ;---PROBability of False Detection-----------
  ;-----------------------
  If i eq 3 then begin 
    If not(keyword_set(skyratio)) then begin
      w = where(y[*,3] lt 2e-3,count)       ;False Detectio lt 2e-3 (3.0 sigma)
      If count ge 1 then begin
        printf,unit,'**********'+Qlabel+': PROB of False Detection (3.0 sigma)********* '
        printf,unit,print_array[*,w]
      endif
    endif
    If keyword_set(skyratio) then begin
      w = where(y[*,3] le 1,count)       ;False Detectio lt 2e-3 (3.0 sigma)
      If count ge 1 then begin
        printf,unit,'**********'+Qlabel+' PROB of False Detection Ratio********* '
        printf,unit,print_array[*,w]
      endif
    endif
  Endif  
endfor

If not(keyword_set(COMBINEHISTO)) then begin
  cgHistoplot, y[*,0],title=Qlabel+' (Integral)',NBINS=10,xtitle='Value'
  cgHistoplot, y[*,1],title=Qlabel+' (On-Off>0)',NBINS=10,xtitle='Value'
  cgHistoplot, y[*,2],title=Qlabel+' (STDDEV)',NBINS=10,xtitle='Value'
  If keyword_set(PROB) then cgHistoplot, y[*,3],title=Qlabel+' (False PROB)',NBINS=10,xtitle='Value'
endif 


end