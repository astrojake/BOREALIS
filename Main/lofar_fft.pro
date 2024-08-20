;;**************************************************************************
;;**************************************************************************
;;                   Lofar FFT on pulsar
;;**************************************************************************
;;**************************************************************************
;Author: Jake Turner (UVA) 
;        J.M. Griessmeier (CNRS) and Philippe Zarka (Obspm)
;        
;Inputs: Pulsarnanme
;        dm
;        Period
;        Pmin
;        fmin
;        fmax 
;        ddm 
;        save_file_name 
; 
;Hard coded: 15 and 7.5 min intervals
;      
;Uses:   Zarka programs     
;         
pro lofar_FFT,pulsarname,dm,Period,pmin,fmin,fmax,ddm,save_file_name

  ;Inputs
  ;      save_file_name
  ;      x
  ;      p
  ;      xt
  ;      xf
  ;      pmin
  ;      fmin
  ;      fmax
  ;      ddm

  ;------------------------
  ;        Pulsar info
  ;------------------------
  ;pulsarname = 'B0823+26'

  ;save_file_name = 'Zarka_test'
  restore,file=save_file_name+'_disp.sav'
  print,'file='+save_file_name+'_disp.sav'

  print,'-------------------------'
  print,'----Start FFT -----------'
  print,'-------------------------'
  ;update inputs
  ;inputs
;    x = x_new_disp
;    p = p_new_disp
;    xt = xt_new_disp
 ; dt = 0.0104857
  xf = xf_disp
 ; xt = xt[0:332000-1]

  xt_org = xt
  x_org  = x
  w1 = where(finite(x,/Nan)) 
  w2 = where(finite(x,/INFINITY)) 
  x[w1] = 0.
  x[w2] = 0.  
  p_org  = p
  p[w1] = 0. 
  p[w2] = 0.
  ;x      = x
    
  nt=n_elements(xt) & nf=n_elements(xf)

  nnv=nt/1000.       ; for visualization in ps file ;2000
  if nnv ne long(nnv) then nnv=long(nnv)+1 else nnv=long(nnv)
  nntv=long(nt/nnv)
  mmv=nf/1000.        ;1000
  if mmv ne long(mmv) then mmv=long(mmv)+1 else mmv=long(mmv)
  mmtv=long(nf/mmv)


  set_plot,'PS'
  device,filename=save_file_name+'_fft_pulsar.ps',/color
;******************************************************
;******************************************************
;If max(xt) lt 960. then begin    ; lt 16 mins just run one fft 
  !p.multi=[0,1,1]
  SPDYNPS,rebin(x(0:nntv*nnv-1,0:mmtv*mmv-1),nntv,mmtv),min(xt),max(xt), $
    min(xf),max(xf),'#spectrum','MHz',pulsarname+' - first slice',0,0,0,0.05,0.95,1,''

  h = [1,2,3,4]
  !p.multi=[0,5,4]
  xm=x/(p+(p eq 0))
 ; xxt=rebin(xt,nt)/3600. & xxtm=rebin(xt,1)/3600.
;  for i=0,nf-1 do begin ;flag
;    plot,xxt,rebin(xm(*,i),nt),xtit='Hour',tit=strtrim(xf(i)-0.390625,2)+'+/- 0.4 MHz',/ynoz,/xsty ; 128*3.0517578 kHz
;    oplot,[0,0]+xxtm(0),[0,max(xm(*,i))*10.],line=1
;  endfor

  !p.multi=[0,1,2]
  xmin=fltarr(nf) & xmax=xmin
  for i=0,nf-1 do begin
    xm(*,i)=median(xm(*,i))
    xmin(i)=min(xm(*,i))
    xmax(i)=max(xm(*,i))
  endfor
  plot_io,xf,xm(0,*)>1,/xsty,xtit='MHz',ytit='median spectrum',psym=-4,tit=pulsarname
 oplot,xf,xmin,line=2
  oplot,xf,xmax,line=2
  for i=0,n_elements(fmin)-1 do oplot,[0,0]+fmin(i),[1,10.*max(xm)],line=1
  for i=0,n_elements(fmax)-1 do oplot,[0,0]+fmax(i),[1,10.*max(xm)],line=1

  pp=reform(rebin(p,1,nf))
  plot,xf,pp,/xsty,xtit='MHz',ytit='fraction of good pixels',psym=-4,tit=pulsarname
  for i=0,n_elements(pmin)-1 do oplot,[0,100],[0,0]+pmin(i),line=1
  nn=0
  ntot=n_elements(pmin)*n_elements(fmin)*n_elements(fmax)*n_elements(ddm)

  for i=0,n_elements(pmin)-1 do begin
    for j=0,n_elements(fmin)-1 do begin
      for k=0,n_elements(fmax)-1 do begin
        for l=0,n_elements(ddm)-1  do begin
          nn=nn+1
          print,nn,' / ',ntot
          titre='>'+strtrim(fix(pmin(i)*100),2)+'%, '+strtrim(fix(fmin(j)),2)+'-'+strtrim(fix(fmax(k)),2)+' MHz'  ; , DM '+strtrim(ddm(l),2)
          xd=DE_DISP(x,xf,dm+ddm(l),dt,max(xf))
          pd=DE_DISP(p,xf,dm+ddm(l),dt,max(xf))
          w=where(pp ge pmin(i) and xf ge fmin(j) and xf le fmax(k))
          print,'LIN - FFT'
          xdi=reform(rebin(xd(*,w),nt,1)) & pdi=reform(rebin(pd(*,w),nt,1))
          y=xdi/(pdi+(pdi eq 0))
          nnt=(1.d0/2.410331d-4)*(dm+ddm(l))*((fmin(j))^(-2)-(fmax(k))^(-2))/dt
          nnt=ceil(nnt/1000.)*1000
          y=y(0:nt-nnt-1)
          fy=abs(fft(y,-1))^2
          ;fy=abs(fft(y,-1))
          xfy=dindgen(nt-nnt)/((nt-nnt)*dt)
          fy=fy(1:*) & xfy=1.d0/xfy(1:*)
          !p.multi=[0,2,2]
          plot,xfy,fy,xra=[0.1,1.2]*period,/xsty,xtit='Period (sec)',ytit='FFT pow (lin)',tit=titre
          oplot,[0,0]+period,[0,max(fy)*10.],line=1
          plot,xfy,fy,xra=[0.8,1.2]*period,/xsty,xtit='Period (sec)',ytit='FFT pow (lin)',tit=titre
          oplot,[0,0]+period,[0,max(fy)*10.],line=1
          FOLD,xfy,fy,[0.9,1.1]*period,xfyf,fyf,h
          ;FOLD,xfy,fy,[0,15728],xfyf,fyf,h
          plot,xfyf,fyf,/xsty,xtit='Period (sec)',ytit='FFT pow (lin) Folded Harmonics',tit=titre
          oplot,[0,0]+period,[0,max(fyf)*10.],line=1
          plot,xfyf,fyf,/xsty,xtit='Period (sec)',ytit='FFT Power (Folded Harmonics)',tit=titre,xra=[0.95,1.05]*period,/ylog
          oplot,[0,0]+period,[0,max(fyf)*10.],line=1
        
          !p.multi=[0,1,1]

          plot,xfy,fy,xra=[0.9,1.1]*period,/xsty,xtit='Period (sec)',ytit='FFT Power',tit=''
          oplot,[0,0]+period,[0,max(fy)*10.],line=1,color='red'

        ;  save,y,xfy,filename='data_BeforeFFT_new.sav'
        ;  
          ;---------------------
          ;     Find SN
          ;---------------------
         uu = where(xfyf gt 0.95d*period and xfyf lt 1.05d*period)
         ; uu = where(xfyf gt 0.80d*period and xfyf lt 1.20d*period)
          fyf_new  = fyf[uu]
          xfyf_new = xfyf[uu]

          max_fft = max(fyf_new, n_max,/nan)
          yfit   = GAUSSFIT(xfyf_new,fyf_new,CHISQ=chi)

          nyy = where(xfyf_new lt 0.99d*period)
          ny = where(xfyf_new gt 1.01*period)
          nnn = [nyy,ny]
          stdev_fft     = STDDEV(fyf_new[nnn])
          mean_baseline = mean(fyf_new[nnn]) 
          max_fft_from_base = max_fft - mean_baseline 
           
          SN      = max_fft/stdev_fft
          SN_base = max_fft_from_base/stdev_fft 
          
          print, '**************************************'
          print, '**************************************'
          print, 'Max FFT Value: ',max_fft
          print, 'Max FFT Value from baseline',mean_baseline
          print, 'Max FFT (s): ', xfyf_new[n_max]
          print, 'Chi-squared of gauss fit', chi
          print, 'STDEV of baseline', stdev_fft
          print, 'S/N', SN
          print, 'S/N from baseline',SN_base

          print,fmin[j], fmax[k], pmin[i], SN, SN_base
         
         ;plot FFT with SN
          !p.multi=[0,1,1]
          cgplot,xfyf_new,fyf_new,/xsty,/ysty,xtit='Period (sec)',ytit='FFT Power (Folded Harmonics)',tit='',charsize=0.9*1.2,/ylog,thick=2;,/xlog;,xra=[0.98,1.05]
          cgplot,[0,0]+period,[1e-14,1e-6],line=1,/overplot,color='red',thick=1.5
          ;oplot,[0,0]+period+5.3e-5,[0,max(fyf_new)*10.],line=1
          ;oplot,xfyf_new,yfit,line=1;,color='red'
          al_legend,[$
            ;'Max FFt: '+strtrim(string(cgnumber_formatter(mean_baseline,decimals=1)),1),$
            ;'STDEV: '+strtrim(string(cgnumber_formatter(Stdev_fft,decimals=1)),1),$
            'SNR: '+strtrim(string(cgnumber_formatter(SN_base,decimals=0)),1),$
            'Total Time: '+strtrim(string(cgnumber_formatter((max(xt) - min(xt)),decimals=3)),1)+' secs'],/top,/right,charsize=0.9*1.2
        
        ;--------Find SNR but not folded----------
        ;---------------------
        ;     Find SN
        ;---------------------
        uu = where(xfyf gt 0.95d*period and xfyf lt 1.05d*period)
        ; uu = where(xfyf gt 0.80d*period and xfyf lt 1.20d*period)
        fyf_new  = fy[uu]
        xfyf_new = xfy[uu]

        max_fft = max(fyf_new, n_max,/nan)
        yfit   = GAUSSFIT(xfyf_new,fyf_new,CHISQ=chi)

        nyy = where(xfyf_new lt 0.99d*period)
        ny = where(xfyf_new gt 1.01*period)
        nnn = [nyy,ny]
        stdev_fft     = STDDEV(fyf_new[nnn])
        mean_baseline = mean(fyf_new[nnn])
        max_fft_from_base = max_fft - mean_baseline

        SN      = max_fft/stdev_fft
        SN_base = max_fft_from_base/stdev_fft

        print, '**************************************'
        print, '**************************************'
        print, 'Not folded (lin)'
        print, 'Max FFT Value: ',max_fft
        print, 'Max FFT Value from baseline',mean_baseline
        print, 'Max FFT (s): ', xfyf_new[n_max]
        print, 'Chi-squared of gauss fit', chi
        print, 'STDEV of baseline', stdev_fft
        print, 'S/N', SN
        print, 'S/N from baseline',SN_base

        print,fmin[j], fmax[k], pmin[i], SN, SN_base

        ;plot FFT with SN
        !p.multi=[0,1,1]
        cgplot,xfyf_new,fyf_new,/xsty,/ysty,xtit='Period (sec)',ytit='FFT Power',tit='',charsize=0.9*1.2,/ylog,thick=2;,/xlog;,xra=[0.98,1.05]
        cgplot,[0,0]+period,[1e-14,1e-6],line=1,/overplot,color='red',thick=1.5
        ;oplot,[0,0]+period+5.3e-5,[0,max(fyf_new)*10.],line=1
        ;oplot,xfyf_new,yfit,line=1;,color='red'
        al_legend,[$
          ;'Max FFt: '+strtrim(string(cgnumber_formatter(mean_baseline,decimals=1)),1),$
          ;'STDEV: '+strtrim(string(cgnumber_formatter(Stdev_fft,decimals=1)),1),$
          'SNR: '+strtrim(string(cgnumber_formatter(SN_base,decimals=0)),1),$
          'Total Time: '+strtrim(string(cgnumber_formatter((max(xt) - min(xt)),decimals=3)),1)+' secs'],/top,/right,charsize=0.9*1.2

        endfor
      endfor
    endfor
  endfor
;  device,/close
;  EXIT_PS
;endif   ;  lt 16 mins
  
;******************************************************
;******************************************************
If (max(xt)-min(xt)) gt 2000. then begin    ; gt 16 mins just run one fft  ;flag

  block_time  = 1680.           ; 28 mins 
  block_time  = 60.            ;15 mins
  
  blocks      = long(max(xt_org)/(block_time))  ;how many blocks of block_time mins
  print,'Max T', max(xt_org)
  print, blocks
  small_steps = (blocks*2) - 1      ;how many steps of 15 mins 
  ;small_steps = 2 ;flag!!
  
  SN_block = dblarr(small_steps)
  reduce_array,xt_org,small_steps,xt_reduce

  min_time = min(xt_org)  
  max_time = min(xt_org) +  block_time
;**************************************  
 for iii=0,small_steps-1 do begin
;**************************************  
  
  ww = where(xt_org ge min_time and xt_org le max_time)  ; find time 
  
  ;If iii eq 2 then begin
  ;  ww = where(xt_org ge min(xt_org) and xt_org le max(xt_org))  ; find time
  ;Endif
  
  xt = xt_org[ww]
  x = x_org[ww,*]
  p = p_org[ww,*]
  ;x = x*p 
  
  nt=n_elements(xt) & nf=n_elements(xf)

  nnv=nt/1000.       ; for visualization in ps file ;2000
  if nnv ne long(nnv) then nnv=long(nnv)+1 else nnv=long(nnv)
  nntv=long(nt/nnv)
  mmv=nf/1000.        ;1000
  if mmv ne long(mmv) then mmv=long(mmv)+1 else mmv=long(mmv)
  mmtv=long(nf/mmv)
  
  !p.multi=[0,1,1]
;  SPDYNPS,rebin(x(0:nntv*nnv-1,0:mmtv*mmv-1),nntv,mmtv),min(xt),max(xt), $
;    min(xf),max(xf),'#spectrum','MHz',pulsarname+' - first slice',0,0,0,0.05,0.95,1,''

  h = [1,2,3,4]
  !p.multi=[0,5,4]
  xm=x/(p+(p eq 0))
  ; xxt=rebin(xt,nt)/3600. & xxtm=rebin(xt,1)/3600.
  ;  for i=0,nf-1 do begin ;flag
  ;    plot,xxt,rebin(xm(*,i),nt),xtit='Hour',tit=strtrim(xf(i)-0.390625,2)+'+/- 0.4 MHz',/ynoz,/xsty ; 128*3.0517578 kHz
  ;    oplot,[0,0]+xxtm(0),[0,max(xm(*,i))*10.],line=1
  ;  endfor

  !p.multi=[0,1,2]
  xmin=fltarr(nf) & xmax=xmin
  for i=0,nf-1 do begin
    xm(*,i)=median(xm(*,i))
    xmin(i)=min(xm(*,i))
    xmax(i)=max(xm(*,i))
  endfor
  ;plot_io,xf,xm(0,*)>1,/xsty,xtit='MHz',ytit='median spectrum',psym=-4,tit=pulsarname
  ;oplot,xf,xmin,line=2
  ;oplot,xf,xmax,line=2
  for i=0,n_elements(fmin)-1 do oplot,[0,0]+fmin(i),[1,10.*max(xm)],line=1
  for i=0,n_elements(fmax)-1 do oplot,[0,0]+fmax(i),[1,10.*max(xm)],line=1

  pp=reform(rebin(p,1,nf))
  ;plot,xf,pp,/xsty,xtit='MHz',ytit='fraction of good pixels',psym=-4,tit=pulsarname
  for i=0,n_elements(pmin)-1 do oplot,[0,100],[0,0]+pmin(i),line=1
  nn=0
  ntot=n_elements(pmin)*n_elements(fmin)*n_elements(fmax)*n_elements(ddm)

  for i=0,n_elements(pmin)-1 do begin
    for j=0,n_elements(fmin)-1 do begin
      for k=0,n_elements(fmax)-1 do begin
        for l=0,n_elements(ddm)-1  do begin
          nn=nn+1
          print,nn,' / ',ntot
          titre='>'+strtrim(fix(pmin(i)*100),2)+'%, '+strtrim(fix(fmin(j)),2)+'-'+strtrim(fix(fmax(k)),2)+' MHz'  ; , DM '+strtrim(ddm(l),2) 
          xd=DE_DISP(x,xf,dm+ddm(l),dt,max(xf))
          pd=DE_DISP(p,xf,dm+ddm(l),dt,max(xf))
          w=where(pp ge pmin(i) and xf ge fmin(j) and xf le fmax(k))
          print,'LIN - FFT'
          xdi=reform(rebin(xd(*,w),nt,1)) & pdi=reform(rebin(pd(*,w),nt,1))
          y=xdi/(pdi+(pdi eq 0))
          nnt=(1.d0/2.410331d-4)*(dm+ddm(l))*((fmin(j))^(-2)-(fmax(k))^(-2))/dt
          nnt=ceil(nnt/1000.)*1000
          y=y(0:nt-nnt-1)
          fy=abs(fft(y,-1))^2
          xfy=dindgen(nt-nnt)/((nt-nnt)*dt)
          fy=fy(1:*) & xfy=1.d0/xfy(1:*)
          !p.multi=[0,2,2]
          plot,xfy,fy,xra=[0.1,1.15],/xsty,xtit='Period (sec)',ytit='FFT pow (lin)',tit=titre
          oplot,[0,0]+period,[0,max(fy)*10.],line=1
          plot,xfy,fy,xra=[0.8,1.2]*period,/xsty,xtit='Period (sec)',ytit='FFT pow (lin)',tit=titre
          oplot,[0,0]+period,[0,max(fy)*10.],line=1
          FOLD,xfy,fy,[0.9,1.1]*period,xfyf,fyf,h
          plot,xfyf,fyf,/xsty,xtit='Period (sec)',ytit='FFT pow (lin) Folded Harmonics',tit=titre
          oplot,[0,0]+period,[0,max(fyf)*10.],line=1
          plot,xfyf,fyf,/xsty,xtit='Period (sec)',ytit='FFT pow (lin) Folded Harmonics',tit=titre,xra=[0.95,1.05]*period
          oplot,[0,0]+period,[0,max(fyf)*10.],line=1

          ;---------------------
          ;     Find SN
          ;---------------------
          uu = where(xfyf gt 0.95d*period and xfyf lt 1.05d*period)
          ; uu = where(xfyf gt 0.80d*period and xfyf lt 1.20d*period)
          fyf_new  = fyf[uu]
          xfyf_new = xfyf[uu]

          max_fft = max(fyf_new, n_max,/nan)
          ;yfit   = GAUSSFIT(xfyf_new,fyf_new,CHISQ=chi)

          nyy = where(xfyf_new lt 0.99d*period)
          ny = where(xfyf_new gt 1.01*period)
          nnn = [nyy,ny]
          stdev_fft     = STDDEV(fyf_new[nnn])
          mean_baseline = mean(fyf_new[nnn])
          max_fft_from_base = max_fft - mean_baseline

          SN      = max_fft/stdev_fft
          SN_base = max_fft_from_base/stdev_fft

          print, '**************************************'
          print, '**************************************'
          print, 'Max FFT Value: ',max_fft
          print, 'Max FFT Value from baseline',mean_baseline
          print, 'Max FFT (s): ', xfyf_new[n_max]
          print, 'Chi-squared of gauss fit', chi
          print, 'STDEV of baseline', stdev_fft
          print, 'S/N', SN
          print, 'S/N from baseline',SN_base

          print,fmin[j], fmax[k], pmin[i], SN, SN_base

          ;plot FFT with SN
          !p.multi=[0,1,1]
          plot,xfyf_new,fyf_new,/xsty,xtit='Period (sec)',ytit='FFT pow (flat lin) Folded Harmonics',tit=titre;,xra=[0.98,1.05]*period
          oplot,[0,0]+period,[0,max(fyf_new)*10.],line=1
          oplot,[0,0]+period+5.3e-5,[0,max(fyf_new)*10.],line=1
          ;oplot,xfyf_new,yfit,line=1;,color='red'
          al_legend,['Max FFt: '+strtrim(string(cgnumber_formatter(mean_baseline,decimals=1)),1),$
            'STDEV: '+strtrim(string(cgnumber_formatter(Stdev_fft,decimals=1)),1),$
            'SN: '+strtrim(string(cgnumber_formatter(SN_base,decimals=1)),1),$
            'Total Time: '+strtrim(string(cgnumber_formatter((max(xt) - min(xt)),decimals=3)),1),$
            'Min Time: '+strtrim(string(cgnumber_formatter(min(xt),decimals=3)),1),$
            'Max Time: '+strtrim(string(cgnumber_formatter((max(xt)),decimals=3)),1)],/top,/right
        endfor
      endfor
    endfor
  endfor
 
  min_time = min_time + block_time/2.
  max_time = max_time + block_time/2.
  SN_block[iii] = SN_base  ;save SN
endfor  ;end 7.5 mins blocks 
save,xt_reduce,SN_block,filename='FFT_time.sav'
nntv = long(n_elements(xt_reduce)/n_elements(SN_Block))
reduce_array,xt_reduce,nntv,xt_block

cgplot,xt_block,SN_block,/xstyle,/ystyle,title='SNR of Pulsar',xtitle='Time Blocks (15 mins)',ytitle='FFT SNR',charsize='1.1'
endif 



device,/close
cgfixps, save_file_name+'_fft_pulsar.ps'
cgps2pdf,save_file_name+'_fft_pulsar.ps',/delete

EXIT_PS
 
end
