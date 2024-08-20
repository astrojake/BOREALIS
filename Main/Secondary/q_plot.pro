pro q_plot,ON,OFF,tau
print, tau
cgScatter2D,ON,OFF,title='Q2 Scatter Plot: '+cgGreek('tau')+'=2'+cgGreek('sigma'),xtitle='ON (sigma)',ytitle='OFF (sigma)',fit=0,$
            /iso,xrange=[-10,10],yrange=[-10,10];,/ystyle,/xstyle
            
cgplot,[0,0],[-1e6,1e6],/overplot,color='black'
cgplot,[-1e6,1e6],[0,0],/overplot
;cgplot,[tau,tau],[-100,100],/overplot,color='red',thick=5
;cgplot,[-100,100],[tau,tau],/overplot,color='red',thick=5

;Fill On beam 
slope = 2.0
xfill = [tau,tau,!X.CRange[1],!X.CRange[1]]
yfill = [tau/slope,!Y.CRange[0],!Y.CRange[0],!X.CRange[1]/slope]
cgcolorfill,xfill,yfill,color='Gray'

;fill Off Beam
cgcolorfill,yfill,xfill,color='Yellow';,color='blue'
cgplotS,xfill,yfill,color='black'
cgplotS,yfill,xfill,color='black'

cgScatter2D,ON,OFF,/overplot,fit=0
cgplot,[0,0],[-1e6,1e6],/overplot,color='black'
cgplot,[-1e6,1e6],[0,0],/overplot
xyouts,!X.CRange[1]/3,!Y.CRange[0]/2,'ON BEAM',charsize=2
xyouts,!X.CRange[0]/1.5,!Y.CRange[1]/1.5,'OFF BEAM',charsize=2

end

