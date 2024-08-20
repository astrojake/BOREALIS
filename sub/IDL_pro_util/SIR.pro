;-----------------------------------------
  pro SIR, para, agg
;-----------------------------------------
; Offringa's Scale-Invariant Rank operator to expand flag table

; para - input 1D array of flags (0=bad , 1=good)
; agg - aggressiveness parameter (0=less, 1=more)

N=n_elements(para)
psi=agg-para

M=[0.,psi]
if N gt 1 then for i=2L,N do M(i)=M(i)+M(i-1)

P=fltarr(N)
for i=1L,N-1 do if M(P(i-1)) gt M(i) then P(i)=i else P(i)=P(i-1)

Q=fltarr(N)+N
for i=N-2L,0,-1 do if M(Q(i+1)) lt M(i+1) then Q(i)=i+1 else Q(i)=Q(i+1)

para=M(Q) lt M(P)

return
end