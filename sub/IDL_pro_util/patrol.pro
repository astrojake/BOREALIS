pro PATROL, data, p, s_sk, m_sk, lc, ls

; v. 14 Apr 2015 

; data(nf,nt) - dynamic spectrum
; p - bad pixel map from previous stages
; l - level above which everything is considered RFI

s=size(data)
nf=s[1] & nt=s[2]
w=where(total(p,2) gt 0)
op=dblarr(nf)-1.d0
op[w]=s_sk[w]/m_sk[w]
Background, op[w], mt, st, nused
testarr=where(abs(op-mt) gt st*lc)
if(testarr[0] gt (-1)) then p[testarr,*]=0

;op=total(data,1)/nf
op=dblarr(nt)
;for i=0, nt-1 do op[i]= total(data[where(p[*,i] eq 1),i])/n_elements(where(p[*,i] eq 1))
op=total(data*p,1)/total(p,1)
w=where(finite(op) ne 1)
if w[0] gt -1 then op[w]=median(op)
Background, op, mt, st, nused
testarr=where(abs(op-mt) gt st*ls)
if(testarr[0] gt (-1)) then p[*,testarr]=0

return
end
