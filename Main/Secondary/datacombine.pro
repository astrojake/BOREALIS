;***************************************************************
;      Combine Lofar exoplanet data for lofar_postprocessing 
;**************************************************************
;***************************************************************
;Author: Jake Turner (UVA)
;        J.M. Griessmeier (CNRS) and Philippe Zarka (Obspm)
;
;Inputs:  file_root      
;         save_file_name ;save file name used for the sav files (e.g. pol0_RFI_patrol5.5_pex_lesig3.5_sum)
;         dates          ;['L429868','L433872','L441630']

;file_root  = '/data/jake.turner/exoplanet/LC5_DDT_002/'
out_ps = 'combine'+save_file_name

set_plot,'PS'
device,filename=out_ps+'.ps',/landscape

  restore,filename=file_root+dates[0]+'/all_beams/rebindata_'+dates[0]+$
                             dates_filenames[0]+'_beam'+$
                             strtrim(string(cgnumber_formatter(beam,decimals=0)),1)+'.sav'
  data_rebin0_1 = data_rebin
  mjd_rebin0_1 = mjd_rebin
  p2_rebin0_1 = p2_rebin
  xt_rebin0_1 = xt_rebin

  restore,filename=file_root+dates[1]+'/all_beams/rebindata_'+dates[1]+$
                             dates_filenames[1]+'_beam'+$
                             strtrim(string(cgnumber_formatter(beam,decimals=0)),1)+'.sav'
  data_rebin0_2 = data_rebin
  mjd_rebin0_2 = mjd_rebin
  p2_rebin0_2 = p2_rebin
  xt_rebin0_2 = xt_rebin + max(xt_rebin0_1) + 1500. ;gap of 5 mins
  
  restore,filename=file_root+dates[2]+'/all_beams/rebindata_'+dates[2]+$
                             dates_filenames[2]+'_beam'+$
                             strtrim(string(cgnumber_formatter(beam,decimals=0)),1)+'.sav'
  data_rebin0_3 = data_rebin
  mjd_rebin0_3 = mjd_rebin
  p2_rebin0_3 = p2_rebin
  xt_rebin0_3 = xt_rebin + max(xt_rebin0_2) + 1500. ;gap of 5 mins

    ;combine
    mjd_rebin0       = [mjd_rebin0_1,mjd_rebin0_2,mjd_rebin0_3]
    xt_rebin0        = [xt_rebin0_1,xt_rebin0_2,xt_rebin0_3]
    xt = xt_rebin0
    nt = n_elements(xt_rebin0)
    nf = n_elements(xf_rebin)
    data_rebin0 = dblarr(nt,nf) 
    p2_rebin0   = dblarr(nt,nf)
    
    for i=0,nf-1 do begin
      data_rebin0(*,i) = [data_rebin0_1[*,i], data_rebin0_2[*,i],data_rebin0_3[*,i]]
      p2_rebin0(*,i) = [p2_rebin0_1[*,i], p2_rebin0_2[*,i],p2_rebin0_3[*,i] ]
    endfor
    
    data_rebin = data_rebin0
    p2_rebin   =  p2_rebin0
    mjd_rebin = mjd_rebin0
    xt_rebin = xt_rebin0
    save,data_rebin,p2_rebin,mjd_rebin,xt_rebin,xf_rebin,filename=out_ps+'_beam0.sav'
    
    ;If keyword_set(PS) and kl eq 0 then begin
      label   = 'Data Combine: Beam 0'
      STANDARD_PLOTS, data_rebin0,p2_rebin0,xt_rebin0,xf_rebin,label
     
      
    ;endif
    
     Num_rebin_t = long(120)
     P  = 0.7365417d         ; (days) exoplanets.eu
     T0  = 2455962.0697d     ;exoplanets.eu for 55 Cnc e (BJD)
     T0  = 2455962.063374324d;JD (from http://astroutils.astronomy.ohio-state.edu/time/bjd2utc.html)
     T0  = T0 - 2400000.500d ;(MJD)
    REDUCE_ARRAY, mjd_rebin0_2,[Num_rebin_t],mjd_rebin0_2_q1,/dbl  ;rebin time
    REDUCE_ARRAY, mjd_rebin0_1,[Num_rebin_t],mjd_rebin0_1_q1,/dbl  ;rebin time
    mjd_rebin0_q1  =[mjd_rebin0_2_q1,mjd_rebin0_1_q1]
;    phase1 = (((mjd_rebin0 - T0) mod P)/(P*1.d)) ;mod 1
;
;    sortp = sort(phase1)
;    mjd_rebin0 = mjd_rebin0[sortp]
;    phase1 = phase1[sortp]
;    xt_rebin0     = xt_rebin0[sortp]
;        data_rebin0 = data_rebin0[sortp,*]
;        p2_rebin0 = p2_rebin0[sortp,*]
        
;REDUCE_ARRAY, mjd_rebin0,[Num_rebin_t],MJD_rebin0_q1,/dbl
    REDUCE_ARRAY, xt_rebin0,[Num_rebin_t],xt_rebin0,/dbl    ;rebin time
    REDUCE_ARRAY, data_rebin0, [Num_rebin_t,nf], x_q1,/dbl              ;rebin data for time series
    REDUCE_ARRAY, p2_rebin0*1., [Num_rebin_t,nf], px_q1,/dbl         ;rebin mask
    x_q1 = x_q1/(px_q1 + (px_q1 eq 0))                         ;apply mask correctly
    
    xt_rebin0 = xt_rebin0/86400.
    
   ; MJD_rebin0_q1 = MJD_rebin0_q1[sort(MJD_rebin0_q1)]
    
   label='3 Dates'
    ytitle = 'Frequency-Integrated Values'
    P  = 0.7365417d         ; (days) exoplanets.eu
    T0  = 2455962.0697d     ;exoplanets.eu for 55 Cnc e (BJD)
    T0  = 2455962.063374324d;JD (from http://astroutils.astronomy.ohio-state.edu/time/bjd2utc.html)
    T0  = T0 - 2400000.500d ;(MJD)
    
    phase1 = (((MJD_rebin0_q1 - T0) mod P)/(P*1.d)) ;mod 1
    ;phase2 = (((t2 - T0) mod P)/(P*1.d)) ;mod 1
    
    cgplot,phase1-0.5,x_q1,title=label,xtitle='Phase',$
      ytitle=ytitle,charsize='1.1',/xstyle,/ystyle,psym=46,yrange=[-1.8e-4,1.9e-4]
    
    diagnostic_plots,x_q1,x_q1,MJD_rebin0_q1,xt_rebin0,P,T0,label,ytitle,t2=MJD_rebin2_q1,/TARGET_ONLY,$
      ; ymax1=0.00025,ymin1=-0.00047,$; 1 hour
      ; ymax2=1.2e-4,ymin2=-2.35e-4;
      ymax1=max(x_q1)*2.0,ymin1=min(x_q1)*2.0,$; 1 hour
      ymax2=max(x_q1)*2.0,ymin2=min(x_q1)*2.0;
      device,/close
end