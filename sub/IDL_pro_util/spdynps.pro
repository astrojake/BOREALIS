; ******************
; *                *
; *   SPDYNPS.PRO  *
; *                *
; ******************

;--------------------------------------------------------------------
  pro SPDYNPS, image,xmin,xmax,ymin,ymax,name_x,name_y,name_plot, $
       grid,backgr,thr,fracmin,fracmax,color,bar_titre, $
	log=log, zeromin=zeromin, iso=iso, screen=screen, newwin=newwin,ytitle_bar=ytitle_bar
;--------------------------------------------------------------------

; PostScript / Screen plot of dynamic spectrum.

; image =  (INPUT)  dynamic spectrum
; xmin,xmax,ymin,ymax =  (INPUT)  ranges of x & y coordinates
; name_x,name_y,name_plot =  (INPUT)  x & y captions + title of plot
; grid =  (INPUT)  =1 to superimpose a grid to the plot, =0 otherwise
; backgr =  (INPUT)  =1 to substract a background to the data before plotting, =0 otherwise
; thr =  (INPUT)  =number of sigmas -in addition to backgr- to substract from data
; fracmin,fracmax =  (INPUT)  min and max threshold for optimization of the dynamic range of the plot
; color =  (INPUT)  grey scale of plot: =0 -> black=intense, =1 -> white=intense
; bar_titre =  (INPUT)  title of dynamic bar (unit = dB, except if '.')
; ytitle_bar = (INPUT)  ytitle of dynamic bar 
; Last parameters can be omitted (beginning anywhere after 'image')

; Default values = 0,n_spectra,0,n_frequencies,' ',' ',' ',0,0,0,.05,.95,0,' '
  if n_elements(xmin) eq 0 or n_elements(xmax) eq 0 then begin
    xmin=0 & xmax=n_elements(image(*,0))
  endif
  if n_elements(ymin) eq 0 or n_elements(ymax) eq 0 then begin
    ymin=0 & ymax=n_elements(image(0,*))
  endif
  if n_elements(name_x) eq 0 then name_x=' '
  if n_elements(name_y) eq 0 then name_y=' '
  if n_elements(name_plot) eq 0 then name_plot=' '
  if n_elements(grid) eq 0 then grid=0
  if n_elements(backgr) eq 0 then backgr=0
  if n_elements(thr) eq 0 then thr=0
  if n_elements(fracmin) eq 0 or n_elements(fracmax) eq 0 then begin
    fracmin=0.05 & fracmax=0.95
  endif
  if  n_elements(color) eq 0 then color=0
  if n_elements(bar_titre) eq 0 then bar_titre=' '
  if n_elements(ytitle_bar) eq 0 then ytitle_bar=' '

  LOADCT, 0
  tab=image

; Background
  if backgr eq 1 then begin
    MAKE_BACKGROUND, tab, '', f,s,n
    for i=0,n_elements(tab(*,0))-1 do tab(i,*)=(tab(i,*)-f-thr*s) >0
  endif

; Dynamic range
  if (fracmin ne 0.) or (fracmax ne 1.) then begin
    tabmin=DYN_N(tab,fracmin) & tabmax=DYN_N(tab,fracmax)
    if (tabmax eq -1) then tabmax=max(tab)
  endif else begin
    tabmin=min(tab)
    tabmax=max(tab)
  endelse
  if keyword_set(zeromin) then tabmin=0.
  if tabmin eq tabmax then begin
    print,'tabmin,max=',tabmin,tabmax
    RETURN
  endif
  tab=bytscl(tab,min=tabmin,max=tabmax)

  if color eq 0 then tab=255-tab  		; 0->noir=intense, 1->blanc=intense
  tl=-0.015
  if grid eq 1 then tl=1.0

  nx=!p.multi(1) & ny=!p.multi(2)
  if nx eq 0 then nx=1
  if ny eq 0 then ny=1

  x=!p.multi(0)
  if x ne 0 then x=nx*ny-x
  i=x mod nx
  j=ny-(x-i)/nx-1

  if keyword_set(screen) then begin
    n1=n_elements(image(*,0)) & n2=n_elements(image(0,*))
    c=1.4 & f=c/2+0.5 & c=c+c*(nx gt 2 or ny gt 2)
    xs1=n1+83*f & ys1=n2+68*f
    if keyword_set(newwin) then window,0,xs=nx*xs1,ys=ny*ys1
    if not keyword_set(log) then $
      plot, [xmin,xmax], [ymin,ymax], /nodata, xra=[xmin,xmax], xstyle=13, yra=[ymin,ymax], ystyle=13,$
            title=name_plot, ticklen=tl, charsize=c, /device, pos=[i*xs1+63*f,j*ys1+45*f,i*xs1+63*f+n1,j*ys1+45*f+n2], iso=iso $
    else plot_io, [xmin,xmax], [ymin,ymax], /nodata, xra=[xmin,xmax], xstyle=13, yra=[ymin,ymax], ystyle=13, title=name_plot, ticklen=tl, charsize=c, /device, pos=[i*xs1+63*f,j*ys1+45*f,i*xs1+63*f+n1,j*ys1+45*f+n2], iso=iso
    tv, tab, !x.window(0), !y.window(0), xsize=!x.window(1)-!x.window(0), ysize=!y.window(1)-!y.window(0),/normal
    axis, xaxis=0, xra=[xmin,xmax], xstyle=1, xticklen=tl, xtitle=name_x, charsize = c
    axis, xaxis=1, xra=[xmin,xmax], xstyle=1, xticklen=tl, xtickname=replicate(' ',9)
    axis, yaxis=0, yra=[ymin,ymax], ystyle=1, yticklen=tl, ytitle=name_y, charsize = c
    axis, yaxis=1, yra=[ymin,ymax], ystyle=1, yticklen=tl, ytickname=replicate(' ',9)
  endif

  if not keyword_set(screen) then begin
    c=1.2 & f=c/2+0.5 & c=c+0.6*(nx gt 2 or ny gt 2)
    ;c=1.5 & f=c/2+0.5 & c=c+0.6*(nx gt 2 or ny gt 2)

    bar_tick = [' ','dB',' ']
    if bar_titre eq '.' then begin
      bar_titre=' ' & bar_tick = [' ',' ',' ']
    endif
    cgplot, [0,1], [0,1], /nodata, xst=1, yst=12, pos=[(i+0.96)/nx-0.06,(j+0.1)/ny+0.075,(i+1.)/nx-0.06,(j+0.95)/ny-0.035], $
      /norm, xticklen=0, xticks=2, xtickv=[0,0.5,1], xtickname=bar_tick, xtitle=bar_titre, charsize=1.7,charthick=3,thick=3
    cgaxis, /yaxis, yra=[tabmin,tabmax], ticklen=-0.02, /ysty, charsize=0.9*1.45,ycharsize=0.9*1.2,ytitle=ytitle_bar,charthick=3
    bcb = bytscl(replicate(1,10)#bindgen(256))
    if color eq 0 then bcb=255-bcb
    tv, bcb, !x.window(0), !y.window(0), xsize=!x.window(1)-!x.window(0), ysize=!y.window(1)-!y.window(0), /normal
    ;tv, bcb,/normal
    oplot, [1,0,0,1], [1,1,0,0]

    if not keyword_set(log) then $
      cgplot, [xmin,xmax], [ymin,ymax], /nodata, /noerase, xra=[xmin,xmax], xstyle=13, yra=[ymin,ymax], ystyle=13, title=name_plot,$
           ticklen=tl, charsize=1.5, /norm, pos=[i*1./nx+0.06,j*1./ny+0.075,(i+0.92)/nx-0.06,(j+1.)/ny-0.035], iso=iso,charthick=3 $
           ;ticklen=tl, charsize=c, /norm, iso=iso $
    else plot_io, [xmin,xmax], [ymin,ymax], /nodata, /noerase, xra=[xmin,xmax], xstyle=13, yra=[ymin,ymax], ystyle=13, title=name_plot,$
         ticklen=tl, charsize=1.5, /norm, pos=[i*1./nx+0.06,j*1./ny+0.075,(i+0.92)/nx-0.06,(j+1.)/ny-0.035], iso=iso
         ; ticklen=tl, charsize=c, /norm, iso=iso
    tv, tab, !x.window(0), !y.window(0), xsize=!x.window(1)-!x.window(0), ysize=!y.window(1)-!y.window(0),/normal
    ;tv, tab,/normal
    cgaxis, xaxis=0, xra=[xmin,xmax], xstyle=1, xticklen=tl, xtitle=name_x, charsize = 1.5,charthick=3  ;xaxis 
    cgaxis, xaxis=1, xra=[xmin,xmax], xstyle=1, xticklen=tl, xtickname=replicate(' ',9),charthick=3
    cgaxis, yaxis=0, yra=[ymin,ymax], ystyle=1, yticklen=tl, ytitle=name_y, charsize = 1.5,charthick=3
    cgaxis, yaxis=1, yra=[ymin,ymax], ystyle=1, yticklen=tl, ytickname=replicate(' ',9),charthick=3
  endif

return
end


;--------------------------------------------------------------------
  pro HELP_SPDYNPS
;--------------------------------------------------------------------
  print,'SPDYNPS, image, xmin,xmax,ymin,ymax, name_x,name_y,name_plot, $'
  print,'         grid,backgr,thr,fracmin,fracmax,color,bar_titre, $'
  print,'         [/log], [/iso], [/zeromin], [/screen], [/newwin]'
return
end
