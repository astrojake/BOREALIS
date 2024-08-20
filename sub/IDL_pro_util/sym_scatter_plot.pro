; creates an elongated scatter plot

set_plot,'ps'
device,filename='Scatter_TEST.ps'
angle=45.

rien=''
x=randomn(seed,1000)
y=randomn(seed,1000)
plot,x,y,psym=4,/iso ,title='symmetrical 2D gaussian' ; symmetrical 2D gaussian
oplot,[0,0],[-10,10],line=1
oplot,[-10,10],[0,0],line=1
;read, rien

x=x*2
plot,x,y,psym=4,/iso,title='2D gaussian elongated x'  ; 2D gaussian elongated along x axis
oplot,[0,0],[-10,10],line=1
oplot,[-10,10],[0,0],line=1
;read, rien

r=sqrt(x^2+y^2)
t=atan(y,x)
x=r*cos(t+angle*!dtor)
y=r*sin(t+angle*!dtor)
plot,x,y,psym=4,/iso,title='2D gaussian elongated along an axis at angle deg'  ; 2D gaussian elongated along an axis at angle deg.
res=poly_fit(x,y,1,yfit=yfit)
oplot,x,yfit,color=200
oplot,[0,0],[-10,10],line=1
oplot,[-10,10],[0,0],line=1
;read, rien

; and corrects it

t=(atan(y,x)*!radeg+180) mod 180
BACKGROUND,t,b,s,n,nsig=2
print,b,s,n ; b is the average slope of the cloud, 1 sigma error  s, computed on n points retained
r=sqrt(x^2+y^2)
t=atan(y,x)
x=r*cos(t-b*!dtor)
y=r*sin(t-b*!dtor)
x=x/stddev(x)
y=y/stddev(y)
plot,x,y,psym=4,/iso,title='Corrected x,y to make symmetrical 2D gaussian'  ; corrected x and y make a symmetrical 2D gaussian
oplot,[0,0],[-10,10],line=1
oplot,[-10,10],[0,0],line=1
device,/close
end
