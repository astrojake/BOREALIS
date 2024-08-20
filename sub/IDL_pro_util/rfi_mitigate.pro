; -----------------------------------------------------------------------
pro RFI_MITIGATE, x, p, time=time,freq=freq,PATROL=PATROL, LE_SIG=LE_SIG, SUM=SUM, SIR_AGG=SIR_AGG, PEX=PEX,$
                  QUIET=QUIET,xnet=xnet,PS=PS,LAST_PS=LAST_PS,$
                  SIGMA_PATROL=SIGMA_PATROL,LESIG_SIGMA=LESIG_SIGMA,NFREQ_PATROL=NFREQ_PATROL,NTIME_PATROL=NTIME_PATROL,$
                  NFREQ_LESIG=NFREQ_LESIG,NTIME_LESIG=NTIME_LESIG,$
                  MF=MF,MT=MT,THRESHOLD=THRESHOLD,PMIN_F=PMIN_F,EXPF=EXPF,PMIN_T=PMIN_T,EXPT=EXPT,$
                  PMIN_F2=PMIN_F2,EXPF2=EXPF2,PMIN_T2=PMIN_T2,EXPT2=EXPT2
  ; -----------------------------------------------------------------------
  ; x(t,f) 
  ; p=0,1 (weights)
 
  ;Flags
  If not(keyword_set(SIGMA_PATROL)) then SIGMA_PATROL = 4.5
  If not(keyword_set(LESIG_SIGMA)) then LESIG_SIGMA = 3.5
  IF NOT(KEYWORD_SET(NFREQ_PATROL)) THEN NFREQ_PATROL = 81
  IF NOT(KEYWORD_SET(NTIME_PATROL)) THEN NTIME_PATROL = 101
  IF NOT(KEYWORD_SET(NFREQ_LESIG)) THEN NFREQ_LESIG = 81
  IF NOT(KEYWORD_SET(NTIME_LESIG)) THEN NTIME_LESIG =  101
;  If not(keyword_set(NFREQ_PATROL)) then NFREQ_PATROL = 51
;  If not(keyword_set(NTIME_PATROL)) then NTIME_PATROL = 51
;  If not(keyword_set(NFREQ_LESIG)) then NFREQ_LESIG = 51
;  If not(keyword_set(NTIME_LESIG)) then NTIME_LESIG =  51
  If not(keyword_set(MF)) then MF = [2,4,6,8,10,16,32,128,256,512]
  If not(keyword_set(MT)) then MT =  [2,4,8,10,16,32,64,128,256]
  If not(keyword_set(THRESHOLD)) then THRESHOLD = 3.5
  If not(keyword_set(BACK)) then BACK = 0 
  
  If not(keyword_set(PMIN_F)) then PMIN_F = 3 
  If not(keyword_set(EXPF)) then EXPF = 5
  If not(keyword_set(PMIN_T)) then PMIN_T = 3
  If not(keyword_set(EXPT)) then EXPT = 5
  
  If not(keyword_set(PMIN_F2)) then PMIN_F2 = 5
  If not(keyword_set(EXPF2)) then EXPF2 = 12
  If not(keyword_set(PMIN_T2)) then PMIN_T2 = 5
  If not(keyword_set(EXPT2)) then EXPT2 = 12
  
  p0 = p  ;input p
  nt=n_elements(x(*,0)) & nf=n_elements(x(0,*))
  wf=lindgen(nf) & wt=lindgen(nt)
   
  if keyword_set(PATROL) then begin
    if not(keyword_set(QUIET)) then print,'patrol ...'
    xx=reform(rebin(x,1,nf))
    print, 'Patrol, nfreq:', NFREQ_PATROL 
    LE_AUTO_S,xx,NFREQ_PATROL,SIGMA_PATROL,0,xnet,pnet         ;f
    p=rebin(reform(pnet,1,nf),nt,nf,/sample)
    p=p*p0
    pnet = reform(rebin(p,1,nf,/sample))
    wf=where(pnet eq 1)
    
    xx=reform(rebin(x(*,wf),nt,1))    
    print, 'Patrol, ntime:', NTIME_PATROL
    print, 'Patrol, Sigma:', SIGMA_PATROL
    LE_AUTO,xx,NTIME_PATROL,SIGMA_PATROL,0,xnet,pnet         ;t
    p=p*rebin(reform(pnet,nt,1),nt,nf,/sample)
    pnet = reform(rebin(p,nt,1,/sample))
    wt=where(pnet eq 1)
          
    ;Plot
    label = 'After Patrol'
    p_before= p 
    x_before =x 
    x2 = x
    ;LOADCT, 0
    rfi_plots,x2,p,time,freq,label=label,PS=PS
    ;x_after = x
   ; stop
  endif ;end patrol

  if keyword_set(LE_SIG) then begin
    if not(keyword_set(QUIET)) then print,'le_auto ...'
    print, 'Le Sig, nfreq:', NFREQ_LESIG
    for i=0,n_elements(wt)-1 do begin
      LE_AUTO_S,x(wt(i),*),NFREQ_LESIG,LESIG_SIGMA,0,xnet,pnet  ;f  ;flag
      p(wt(i),*)=p(wt(i),*)*pnet
    endfor
    
    print, 'Le Sig, ntime:', NTIME_LESIG
    print, 'LE_SIG Sigma:' ,LESIG_SIGMA
    
    for i=0,n_elements(wf)-1 do begin
      LE_AUTO,x(*,wf(i)),NTIME_LESIG,LESIG_SIGMA,0,xnet,pnet  ;t
      p(*,wf(i))=p(*,wf(i))*pnet
    endfor
    
    ;Plot   
    label = 'After Le SIG'
    x2 = x
    ;LOADCT, 0
    rfi_plots,x2,p,time,freq,label=label,PS=PS,p_before=p_before
    p_before = p
  endif

  if keyword_set(SUM) then begin
    if not(keyword_set(QUIET)) then print,'sumthreshold ...'
    xx=x(wt,*) & xx=xx(*,wf)
    pp=p(wt,*) & pp=pp(*,wf)
    
    print,'Sum, mf', MF 
    print, 'Sum, mt', MT
    print,'Threshold',THRESHOLD
    print,'PP t size: ',n_elements(pp[*,0]),'PP f size: ',n_elements(pp[0,*]) 

    If n_elements(pp[*,0]) gt 10 and n_elements(pp[0,*]) gt 10 then begin 
       launch_sumthr, xx, pp,MF,MT, THR0=THRESHOLD ,BACK=BACK
       for i=0,n_elements(wt)-1 do p(wt(i),wf)=pp(i,*)
    endif 

    ;Plot
    label = 'After SUM'
    x2 = x
    ;LOADCT, 0
    rfi_plots,x2,p,time,freq,label=label,PS=PS,p_before=p_before
     p_before = p
  endif

  if keyword_set(SIR_AGG) then begin
    if not(keyword_set(QUIET)) then print,'sir ...'
    for i=0,n_elements(wt)-1 do begin
      pp=reform(p(wt(i),*))
      SIR, pp, SIR_AGG
      p(wt(i),*)=pp
    endfor
    for i=0,n_elements(wf)-1 do begin
      pp=reform(p(*,wf(i)))
      SIR, pp, SIR_AGG
      p(*,wf(i))=pp
    endfor
  endif

  if keyword_set(PEX) then begin
    if not(keyword_set(QUIET)) then print,'tfex ...'
    
    tex_frex, p, [EXPF,PMIN_F], [EXPT,PMIN_T], /verbose
    tex_frex, p, [EXPF2,PMIN_F2], [EXPT2,PMIN_T2], /verbose

    ;Plot
    label = 'After PEX'
    x2 = x
    ;LOADCT, 0
    rfi_plots,x2,p,time,freq,label=label,PS=PS
  endif
  
  If keyword_set(LAST_PS) then begin 
    x2 = x
    rfi_plots,x2,p,time,freq,label=label,/LAST_PS
  endif 
  
  return
end
