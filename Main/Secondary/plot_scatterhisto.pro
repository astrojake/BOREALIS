pro plot_scatterHisto,x,y,xlabel=xlabel,ylabel=ylabel
  page_height = 27.94
  page_width = 21.59
  
  xmin = -10 
  xmax = 10 
  ymin = -10
  ymax = 10
  
  yhist = cgHistogram(y,binsize=0.1,NBINS=NBINS)
  ;print,yhist
  ;---------------------
  ;Q2 Histogram
  ;---------------------
  ; define position of plot
  ; bottom left corner (in cm)
  plot_left = 3
  plot_bottom = 12.4

  ; define plot size (in cm)
  xsize = 8
  ysize = 8
  ; create rotated histogram
  cghistoplot, x, $
    charsize = 1.2, $
    ;xtitle = "N", $
    ytitle = "", $
    xrange= [xmin,xmax], $
    yrange= [0., max(yhist)/6.], $
    thick = 5, $
    charthick = 3, $
    xstyle = 9, $
    xtickname = replicate(' ', 60), $
    bin = 0.1, $
    datacolorname = "grn5",$
    /normal, $
   ;  /outline , $
    ; /log,$
    position = [plot_left /page_width, plot_bottom/page_height, (plot_left + xsize)/page_width, (plot_bottom + ysize)/page_height]

  ; second histogram rotated  ____________________________________________________
  ; define position of plot
  ; bottom left corner (in cm)
  plot_left = 11.4
  plot_bottom = 4

  ; define plot size (in cm)
  xsize = 8
  ysize = 8


  ; create rotated histogram
  cghistoplot, y, $
    charsize = 1.2, $
    xtitle = "", $
    ;ytitle = "data", $
    thick = 5, $
    charthick = 3, $
    yrange= [ymin, ymax], $
    xrange= [0., max(yhist)/6.], $
    ystyle = 9, $
    ytickname = replicate(' ', 60), $
    bin = 0.1, $
    rotate=1, $
    /noerase, $
    ;/outline , $
    datacolorname = "grn5",$
    /normal, $
   ; /log,$
    position = [plot_left /page_width, plot_bottom/page_height, (plot_left + xsize)/page_width, (plot_bottom + ysize)/page_height]

  ; second histogram rotated  ____________________________________________________
  ; main plot  ____________________________________________________
  ; define position of plot
  ; bottom left corner (in cm)
  plot_left = 3
  plot_bottom = 4

  ; define plot size (in cm)
  xsize = 8
  ysize = 8

  cgplot, x, y , $
    xcharsize=0.8  , ycharsize=0.8   ,$
    xrange= [xmin, xmax], $
    yrange= [ymin, ymax], $
    charthick = 3, $
    ytitle =ylabel,$
    xtitle =xlabel,$
    xthick = 4, $
    ythick = 4, $
    /noerase, $
    xstyle =1,ystyle =1 , /nodata, $
    /normal, $
    position = [plot_left /page_width, plot_bottom/page_height, (plot_left + xsize)/page_width, (plot_bottom + ysize)/page_height]


  cgplot, x, y , $
    psym = symcat(16), /overplot,/noerase, $
    color = "blu7" ,symsize= 0.5
  ; main plot  ____________________________________________________


end