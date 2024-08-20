function SIGN, x
y=x
w=where(x ne 0)
if w(0) ne -1 then y(w)=x(w)/abs(x(w))
y=fix(temporary(y))
return,y
end