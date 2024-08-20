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
;
;Note: This program is general with the exception of where the files are located
;       eg. file_root+/date/ or /data/jake.turner/exoplanet/LC5_DDT_002/L527649 
;
pro BOREALIS_combine,file_root,save_file_name,dates,dates_filenames,beam,LIST=LIST,TO=TO,P=P;,save_file_name2=save_file_name2


out_ps = 'combine_'+save_file_name+'beam'+strtrim(string(fix(beam)),1)
set_plot,'PS'
device,filename=out_ps+'.ps',/landscape
  
  n    = n_elements(dates)
  LIST = dindgen(2,n) 
  data_rebin_p    = ptrarr(n)   ;pointer
  mjd_rebin_p     = ptrarr(n)
  xt_rebin_p      = ptrarr(n)
  p2_rebin_p      = ptrarr(n)
  UT_rebin_p      = ptrarr(n)
  phase_p         = ptrarr(n)

  ;*********************************************************
  ;                         Combine dates
  ;*********************************************************
  for i=0,n-1 do begin 
    restore_filename = file_root+'/'+dates[i]+'/rebindata_'+$
      dates_filenames[i]+'beam'+$
      strtrim(string(fix(beam)),1)+'.sav'

    print,'Restore file: ',restore_filename
    restore,filename=restore_filename

    phase = (((MJD_REBIN - TO) mod P)/(P*1.d))

   ;setup pointers 
    data_rebin_p[i]  = ptr_new(data_rebin)
    mjd_rebin_p[i]   = ptr_new(mjd_rebin)
    ut_rebin_p[i]   = ptr_new(UT_rebin)
    If i eq 0 then xt_rebin_p[i]    = ptr_new(xt_rebin)
    If i gt 0 then xt_rebin_p[i]    = ptr_new(xt_rebin + max( *xt_rebin_p[i-1] ) + 500. )
    p2_rebin_p[i]    = ptr_new(p2_rebin)
    phase_p[i]    = ptr_new(phase)

    LIST[0,i]  = 0 
    LIST[1,i]  = LIST[0,i] + n_elements(xt_rebin[i,*]) 
  
    ;Concatenate arrays 
    If i eq 0 then begin
      xt_rebin_all =  *xt_rebin_p[i]
      mjd_rebin_all =  *mjd_rebin_p[i]
      UT_rebin_all =  *UT_rebin_p[i]
      phase_all =  *phase_p[i]
    Endif
    If i gt 0 then begin
      xt_rebin_all  =  [ xt_rebin_all,*xt_rebin_p[i] ]
      mjd_rebin_all =  [ mjd_rebin_all,*mjd_rebin_p[i] ]
      UT_rebin_all  =  [ UT_rebin_all,*UT_rebin_p[i] ]
      phase_all     =  [ phase_all,*phase_p[i] ]
    endif
    
    If i eq 0 then LIST[0,i]  = 0
    If i gt 0 then LIST[0,i]  = LIST[1,i-1] + 1L 
    LIST[1,i]  = LIST[0,i] + n_elements(*xt_rebin_p[i]) - 1L
  endfor 
  
  nt = n_elements(xt_rebin_all)
  nf = n_elements(xf_rebin)
  data_rebin = dblarr(nt,nf)
  p2_rebin   = dblarr(nt,nf)
  na = intarr(n)
  
  start = 0L &enda  = 0L
  ;Concatenate arrays for data_rebin and p2_Rebin
  for j=0,n-1 do begin ;n dates
    a                = *data_rebin_p[j]
    b                = *p2_rebin_p[j]
    na[j] = n_elements(a[*,0])
    If j eq 0 then data_rebin(0:na[j]*1L-1L,*)  = a(*,*)
    If j eq 0 then p2_rebin(0:na[j]*1L-1L,*)    = b(*,*)
    If j gt 0 then begin
      start = start + na[j-1]*1L
      enda   = start + na[j]*1L
      data_rebin(start:enda-1L,*)  = a(*,*)
      p2_rebin(start:enda-1L,*)    = b(*,*)
    Endif
  endfor
  
  ;*********************
  ;        Save
  ;*********************
  xt_rebin = xt_rebin_all
  mjd_rebin = mjd_rebin_all
  save,data_rebin,p2_rebin,mjd_rebin,xt_rebin,xf_rebin,filename='rebindata_'+out_ps+'.sav'

  label   = 'Data Combine: Beam '+strtrim(string(cgnumber_formatter(beam,decimals=0)),1)+'.sav'
  STANDARD_PLOTS, data_rebin,p2_rebin,xt_rebin,xf_rebin,label
  print,'********Done Combining Beam '+strtrim(string(cgnumber_formatter(beam,decimals=0)),1)+'*********'


;  ;Phase sort
;  sort_index =sort(phase)
;  phase     = phase[sort_index]
;  mjd_rebin = mjd_rebin[sort_index]
;  UT_rebin = UT_rebin[sort_index]
;  data_rebin = data_rebin[sort_index,*]
;  p2_rebin = p2_rebin[sort_index,*]
;
;  ;save phase 
;  save,phase,data_rebin,p2_rebin,mjd_rebin,UT_rebin,xt_rebin,xf_rebin,filename='rebindata_'+out_ps+'.sav'

 
  cgfixps,out_ps+'.ps'
      device,/close
    
    ;clear heap memory
      ptr_free,data_rebin_p
      ptr_free,mjd_rebin_p
      ptr_free,xt_rebin_p
      ptr_free,p2_rebin_p
      ptr_free,UT_rebin_p
      ptr_free,phase_p

end