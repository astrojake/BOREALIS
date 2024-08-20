;***************************************************************
;      Combine masks of Lofar exoplanet data for lofar_processing etc...
;**************************************************************
;***************************************************************
;Author: Jake Turner (UVA)
;        J.M. Griessmeier (CNRS) and Philippe Zarka (Obspm)
;
;Inputs:  
;         save_file_name ;save file name used for the sav files (e.g. L429868_pol0_15770sec_RFI_patrol5.5_pex_lesig3.5_sum)
; 
;Flags    /bit_mask       ;Input bit masks for beams (save lots of virtual memory)
;         /PS
;         
;Assumes:  4 beams
pro lofar_combinemaskold,save_file_name,bit_mask=bit_mask,save_bitmask=save_bitmask,verbose=verbose,PS=PS

  If keyword_set(verbose) then begin
    if n_elements(bit_mask) ne 0 then $
      print, 'Mask Bit: set' else $
      print, 'Mask Bit: not set'
  endif 
  
  If keyword_set(PS) then begin
    print,'----Making Plot----'
    set_plot,'PS'
    device,filename=save_file_name+'_maskcombine.ps',/landscape
  endif
  
;bit_mask = 0  ;0 = N, 1 = YES
;save_file_name='L429868_pol0_15770sec_RFI_patrol5.5_pex_lesig3.5_sum'

If not(keyword_Set(bit_mask)) then begin
  print,'Bit mask not set'
  ;---------------------------------
  ;          Restore Beam 0 
  ;---------------------------------
  restore, filename='mask_full_'+save_file_name+'_beam0.sav'
  ;xt,xf,p2_full restored
  p2_full_0 = temporary(p2_full)
    ww0 = where(p2_full_0 eq 0)  ;bad pixels
    p0=n_elements(ww0)*100.d0/n_elements(p2_full_0)   ;percent of bad pixels
    print,'# polluted channels for Beam 0 = ',n_elements(ww0),' / ',n_elements(p2_full_0)
    print,'RFI Mask for Beam 0:       ',p0,' % -> masked out'
  
  nt = n_elements(xt) 
  If keyword_Set(save_bitmask) then begin
    ;lofar_bytebit,p2_full_0,nt,p2_out=p2_full_bit_0,/BYTE_TO_BIT  
    BYTEARRAY_TO_BITARRAY,p2_full_0,xsize,p2_full_bit_0
    save,xt,xf,xsize,p2_full_bit_0,filename='mask_full_bit_'+save_file_name+'_beam0.sav'
  endif 
  
  If keyword_set(PS) then begin
    reduce_array,p2_full_0,[1000,1],p2_plot
    reduce_array,xt,[1000],xt_new
    STANDARD_PLOTS, p2_plot*1.0,p2_plot*1.0,xt_new,xf,'Mask: Beam 0',/mask,/ONLY_PLOT
  endif 
 
    ;---------------------------------
    ;          Restore Beam 1
    ;---------------------------------  
  restore, filename='mask_full_'+save_file_name+'_beam1.sav'
  ;xt,xf,p2_full restored
  p2_full_1 = temporary(p2_full)
    ww1 = where(p2_full_1 eq 0)  ;bad pixels
    p0=n_elements(ww1)*100.d0/n_elements(p2_full_1)   ;percent of bad pixels
    print,'# polluted channels for Beam 1 = ',n_elements(ww1),' / ',n_elements(p2_full_1)
    print,'RFI Mask for Beam 1:       ',p0,' % -> masked out'
    
     If keyword_Set(save_bitmask) then begin
       lofar_bytebit,p2_full_1,nt,p2_out=p2_full_bit_1,/BYTE_TO_BIT   
       save,xt,xf,p2_full_bit_1,filename='mask_full_bit_'+save_file_name+'_beam1.sav'
     endif 
   
 If keyword_set(PS) then begin
    reduce_array,p2_full_1,[1000,1],p2_plot
    reduce_array,xt,[1000],xt_new
    STANDARD_PLOTS, p2_plot*1.0,p2_plot*1.0,xt_new,xf,'Mask: Beam 1',/mask,/ONLY_PLOT
  endif   
  ;---------------------------------
  ;          Restore Beam 2
  ;---------------------------------
  restore, filename='mask_full_'+save_file_name+'_beam2.sav'
  ;xt,xf,p2_full restored
  p2_full_2 = temporary(p2_full)
    ww2 = where(p2_full_2 eq 0)  ;bad pixels
    p0=n_elements(ww2)*100.d0/n_elements(p2_full_2)   ;percent of bad pixels
    print,'# polluted channels for Beam 2 = ',n_elements(ww2),' / ',n_elements(p2_full_2)
    print,'RFI Mask for Beam 2:       ',p0,' % -> masked out'
   
   If keyword_Set(save_bitmask) then begin
     lofar_bytebit,p2_full_2,nt,p2_out=p2_full_bit_2,/BYTE_TO_BIT      
     save,xt,xf,p2_full_bit_2,filename='mask_full_bit_'+save_file_name+'_beam2.sav'
   endif
   
 If keyword_set(PS) then begin
    reduce_array,p2_full_2,[1000,1],p2_plot
    reduce_array,xt,[1000],xt_new
    STANDARD_PLOTS, p2_plot*1.0,p2_plot*1.0,xt_new,xf,'Mask: Beam 2',/mask,/ONLY_PLOT
  endif 
    ;---------------------------------
    ;          Restore Beam 3
    ;---------------------------------
  restore, filename='mask_full_'+save_file_name+'_beam3.sav'
  ;xt,xf,p2_full restored
  p2_full_3 = temporary(p2_full)
    ww3 = where(p2_full_3 eq 0)  ;bad pixels
    p0=n_elements(ww3)*100.d0/n_elements(p2_full_3)   ;percent of bad pixels
    print,'# polluted channels for Beam 3 = ',n_elements(ww3),' / ',n_elements(p2_full_3)
    print,'RFI Mask for Beam 3:      ',p0,' % -> masked out'
   
   If keyword_Set(save_bitmask) then begin
    lofar_bytebit,p2_full_3,nt,p2_out=p2_full_bit_3,/BYTE_TO_BIT  
    save,xt,xf,p2_full_bit_3,filename='mask_full_bit_'+save_file_name+'_beam3.sav'
   endif  
   
 If keyword_set(PS) then begin
    reduce_array,p2_full_3,[1000,1],p2_plot
    reduce_array,xt,[1000],xt_new
    STANDARD_PLOTS, p2_plot*1.0,p2_plot*1.0,xt_new,xf,'Mask: Beam 3',/mask,/ONLY_PLOT
  endif 
   
  ;---------------------------------
  ;          Combine 
  ;---------------------------------
  p2_full_all = temporary(p2_full_0)*temporary(p2_full_1)*temporary(p2_full_2)*temporary(p2_full_3)
  
    ww = where(p2_full_all eq 0)  ;bad pixels
    p0=n_elements(ww)*100.d0/n_elements(p2_full_all)   ;percent of bad pixels
    print,'# polluted channels for all Beams = ',n_elements(ww),' / ',n_elements(p2_full_all)
    print,'RFI Mask for all beams:       ',p0,' % -> masked out'
    
    If keyword_set(PS) then begin
    reduce_array,p2_full_all,[1000,1],p2_plot
    reduce_array,xt,[1000],xt_new
    STANDARD_PLOTS, p2_plot*1.0,p2_plot*1.0,xt_new,xf,'Mask: All',/mask,/ONLY_PLOT
  endif 
    
  If keyword_Set(save_bitmask) then begin  
    lofar_bytebit,p2_full_all,nt,p2_out=p2_full_all_bit,/BYTE_TO_BIT  
    save,xt,xf,p2_full_all_bit,filename='mask_combine_bit_'+save_file_name+'_beam0.sav'
  endif 
  save,xt,xf,p2_full_all,filename='mask_combine_'+save_file_name+'.sav'
Endif

;************************************************************
;                    Mask Bit
;************************************************************
If keyword_Set(bit_mask) then begin
  ;------------------
  ;  Restore Beam 0 
  ;------------------
  restore, filename='mask_full_bit_'+save_file_name+'_beam0.sav'
    ;xt,xf,p2_full_bit restored
  nt = n_elements(xt)
  p2_full_bit_0 = temporary(p2_full_bit)
  lofar_bytebit,p2_full_bit_0,nt,p2_out=p2_full_0,/BIT_TO_BYTE
    
  ;------------------
  ;  Restore Beam 1
  ;------------------
  restore, filename='mask_full_bit_'+save_file_name+'_beam1.sav'
  ;xt,xf,p2_full restored
  p2_full_bit_1 = temporary(p2_full_bit)
  lofar_bytebit,p2_full_bit_1,nt,p2_out=p2_full_1,/BIT_TO_BYTE
  
  
  ;------------------
  ;  Restore Beam 2
  ;------------------
  restore, filename='mask_full_bit_'+save_file_name+'_beam2.sav'
  ;xt,xf,p2_full restored
  p2_full_bit_2 = temporary(p2_full_bit)
  lofar_bytebit,p2_full_bit_2,nt,p2_out=p2_full_2,/BIT_TO_BYTE
  
  ;------------------
  ;  Restore Beam 3
  ;------------------
  restore, filename='mask_full_bit_'+save_file_name+'_beam3.sav'
  ;xt,xf,p2_full restored
  p2_full_bit_3 = temporary(p2_full_bit)
  lofar_bytebit,p2_full_bit_3,nt,p2_out=p2_full_3,/BIT_TO_BYTE
  
  ;---------------------------------
  ;          Combine 
  ;---------------------------------
  p2_full_all = temporary(p2_full_0)*temporary(p2_full_1)*temporary(p2_full_2)*temporary(p2_full_3)
    ww = where(p2_full_all eq 0)  ;bad pixels
    p0=n_elements(ww)*100.d0/n_elements(p2_full_all)   ;percent of bad pixels
    print,'# polluted channels for all Beams = ',n_elements(ww),' / ',n_elements(p2_full_all)
    print,'RFI Mask for all beams:       ',p0,' % -> masked out'
   
  If keyword_set(PS) then begin
   nt = n_elements(xt)
   Num_rebin_t = 1000. ;reduction factor 
   num_time  = 4000.  ;size of steps 
   steps     = floor(nt/num_time) ;number of steps
   p2_plot = dblarr(floor(nt/Num_rebin_t),1)    ;size of rebined dataa
   
   nnv_2 = Num_rebin_t
   if nnv_2 ne long(nnv_2) then nnv_2=long(nnv_2)+1 else nnv_2=long(nnv_2)
   Num_rebin_t = nnv_2
   nntv_2=long(num_time/nnv_2)
   nntv_3 = long((num_time*steps)/nnv_2)
   
   
   nt_new = floor(nt/Num_rebin_t)    ;ultra memory saver
   nntv_3 = long((num_time*steps)/nnv_2) 
   
   for i=0,steps-1 do begin
     reduce_array,p2_full_all(i*num_time:(i+1)*num_time-1,*),[Num_rebin_t,1],p2_temp
     p2_plot(i*nntv_2:(i+1)*nntv_2-1,*) = p2_temp(0:nntv_2-1,*)  ;update rebin array
     reduce_array,xt(i*num_time:(i+1)*num_time-1,*),[Num_rebin_t],xt_r
     xt_new(i*nntv_2:(i+1)*nntv_2-1,*) = xt_r(0:nntv_2-1,*)
   endfor
    
   ;reduce_array,p2_full_all,[time_reduce,1],p2_plot
   
   STANDARD_PLOTS, p2_plot*1.0,p2_plot*1.0,xt_new,xf,'Mask: All',/mask,/ONLY_PLOT
  endif ;end plot
    
  lofar_bytebit,p2_full_all,nt,p2_out=p2_full_all_bit,/BYTE_TO_BIT   
  ;save,xt,xf,p2_full_all_bit,filename='mask_combine_bit_'+save_file_name+'_beam0.sav'
  save,xt,xf,p2_full_all,filename='mask_combine_frombit_'+save_file_name+'.sav'
Endif

If keyword_set(PS) then begin
  device,/close
  cgfixps, save_file_name+'_maskcombine.ps'
  EXIT_PS
endif

end