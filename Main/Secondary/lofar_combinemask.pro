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
;Flags    /mask_bit       ;Input mask is in bits for beams (save lots of virtual memory)
;         /save_bitmask   ;If mask not in bits, save in bits
;         /PS
;         /verbose
;         /bitstep        ;mask is saved as a bit pointer (will save memory)
;         
;Last Update: July 21, 2017 
;
;Updates to Perform: - Use bit pointer to save memory
;                    - create arbitary number of masks to combine 
;         
;RAM Usage: 110.16529 Gb for byte masks
;Assumes:  Hardcoded to allow up to four beams to be combined 

pro lofar_combinemask,save_file_name,mask_bit=mask_bit,save_bitmask=save_bitmask,verbose=verbose,PS=PS,BEAMS=BEAMS,bitstep=bitstep
  st= SYSTIME(1)
  startmem_RFI =  MEMORY(/CURRENT)

  If keyword_set(verbose) then begin
    if n_elements(mask_bit) ne 0 then $
      print, 'Mask Bit: set' else $
      print, 'Mask Bit: not set'
    print,'Beams to run: ', BEAMS
  endif 
  
  If keyword_set(PS) then begin
    print,'----Making Plot----'
    set_plot,'PS'
    device,filename=save_file_name+'_maskcombine.ps',/landscape
  endif
  
If not(keyword_Set(mask_bit)) then begin
  ;---------------------------------
  ;          Restore Beam 0
  ;---------------------------------
  run_0 = where(BEAMS eq 0,/null)
  If run_0 ne !Null then begin 
    restore, filename='mask_full_'+save_file_name+'_beam0.sav'
    ;xt,xf,p2_full restored
    p2_full_0 = temporary(p2_full)
      ww0 = where(p2_full_0 eq 0)  ;bad pixels
      p0=n_elements(ww0)*100.d0/n_elements(p2_full_0)   ;percent of bad pixels
      print,'# polluted channels for Beam 0 = ',n_elements(ww0),' / ',n_elements(p2_full_0)
      print,'RFI Mask for Beam 0:       ',p0,' % -> masked out'
    
    nt = n_elements(xt) 
    If keyword_Set(save_bitmask) then begin
      BYTEARRAY_TO_BITARRAY,p2_full_0,xsize,p2_full_bit_0 
      save,xt,xf,xsize,p2_full_bit_0,filename='mask_full_bit_'+save_file_name+'_beam0.sav'
    endif 
    
    If keyword_set(PS) then begin
      reduce_array,p2_full_0,[1000,1],p2_plot
      reduce_array,xt,[1000],xt_new
      STANDARD_PLOTS, p2_plot*1.0,p2_plot*1.0,xt_new,xf,'Mask: Beam 0',/mask,/ONLY_PLOT
    endif 
     
    p2_full_all = temporary(p2_full_0) 
  endif
  ;---------------------------------
  ;          Restore Beam 1
  ;---------------------------------
  run_1 = where(BEAMS eq 1,/null)
  If run_1 ne !Null then begin 
      restore, filename='mask_full_'+save_file_name+'_beam1.sav'
      ;xt,xf,p2_full restored
      p2_full_1 = temporary(p2_full)
        ww1 = where(p2_full_1 eq 0)  ;bad pixels
        p0=n_elements(ww1)*100.d0/n_elements(p2_full_1)   ;percent of bad pixels
        print,'# polluted channels for Beam 1 = ',n_elements(ww1),' / ',n_elements(p2_full_1)
        print,'RFI Mask for Beam 1:       ',p0,' % -> masked out'
        
         If keyword_Set(save_bitmask) then begin
           BYTEARRAY_TO_BITARRAY,p2_full_1,xsize,p2_full_bit_1  
           save,xt,xf,xsize_p2_full_bit_1,filename='mask_full_bit_'+save_file_name+'_beam1.sav'
         endif 
       
     If keyword_set(PS) then begin
        reduce_array,p2_full_1,[1000,1],p2_plot
        reduce_array,xt,[1000],xt_new
        STANDARD_PLOTS, p2_plot*1.0,p2_plot*1.0,xt_new,xf,'Mask: Beam 1',/mask,/ONLY_PLOT
      endif   
  endif 
  
  If run_0 eq !Null and run_1 ne !NUll then p2_full_all = temporary(p2_full_1)
  If run_0 ne !Null and run_1 ne !NUll then p2_full_all = temporary(p2_full_1)*p2_full_all
  
  ;---------------------------------
  ;          Restore Beam 2
  ;---------------------------------
  run_2 = where(BEAMS eq 2,/null)
  If run_2 ne !Null then begin
        restore, filename='mask_full_'+save_file_name+'_beam2.sav'
        ;xt,xf,p2_full restored
        p2_full_2 = temporary(p2_full)
        ww2 = where(p2_full_2 eq 0)  ;bad pixels
        p0=n_elements(ww2)*100.d0/n_elements(p2_full_2)   ;percent of bad pixels
        print,'# polluted channels for Beam 2 = ',n_elements(ww2),' / ',n_elements(p2_full_2)
        print,'RFI Mask for Beam 2:       ',p0,' % -> masked out'
       
      If keyword_Set(save_bitmask) then begin
         BYTEARRAY_TO_BITARRAY,p2_full_2,xsize,p2_full_bit_2        
         save,xt,xf,xsize,p2_full_bit_2,filename='mask_full_bit_'+save_file_name+'_beam2.sav'
       endif
       
       If keyword_set(PS) then begin
        reduce_array,p2_full_2,[1000,1],p2_plot
        reduce_array,xt,[1000],xt_new
        STANDARD_PLOTS, p2_plot*1.0,p2_plot*1.0,xt_new,xf,'Mask: Beam 2',/mask,/ONLY_PLOT
       endif 
  endif 
  
   ;0= no, 1 = no, 2 = yes
  If run_0 eq !Null and run_1 eq !NUll and run_2 ne !NULL then p2_full_all = temporary(p2_full_2)
  If run_0 ne !Null or run_1 ne !NUll then begin
    If run_2 ne !Null then p2_full_all = p2_full_all*temporary(p2_full_2)
  Endif 
    ;---------------------------------
    ;          Restore Beam 3
    ;---------------------------------
    run_3 = where(BEAMS eq 3,/null)
    If run_3 ne !Null then begin
      restore, filename='mask_full_'+save_file_name+'_beam3.sav'
        ;xt,xf,p2_full restored
      p2_full_3 = temporary(p2_full)
      ww3 = where(p2_full_3 eq 0)  ;bad pixels
      p0=n_elements(ww3)*100.d0/n_elements(p2_full_3)   ;percent of bad pixels
      print,'# polluted channels for Beam 3 = ',n_elements(ww3),' / ',n_elements(p2_full_3)
      print,'RFI Mask for Beam 3:      ',p0,' % -> masked out'
     
      If keyword_Set(save_bitmask) then begin
         BYTEARRAY_TO_BITARRAY,p2_full_3,xsize,p2_full_bit_3 
         save,xt,xf,xsize,p2_full_bit_3,filename='mask_full_bit_'+save_file_name+'_beam3.sav'
      endif  
     
      If keyword_set(PS) then begin
       reduce_array,p2_full_3,[1000,1],p2_plot
       reduce_array,xt,[1000],xt_new
       STANDARD_PLOTS, p2_plot*1.0,p2_plot*1.0,xt_new,xf,'Mask: Beam 3',/mask,/ONLY_PLOT
      endif 
   endif 
   
     If run_0 eq !Null and run_1 eq !NUll and run_2 eq !NULL and run_3 ne !NULL then p2_full_all = temporary(p2_full_3)
     If run_0 ne !Null or run_1 ne !NUll or run_2 ne !NUll then begin
       If run_3 ne !Null then p2_full_all = p2_full_all*temporary(p2_full_3)
     Endif
     
  ;---------------------------------
  ;          Combine 
  ;--------------------------------- 
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
     BYTEARRAY_TO_BITARRAY, p2_full_all, xsize, p2_full_all_bit
    save,xt,xf,xsize,p2_full_all_bit,filename='mask_combine_bit_'+save_file_name+'.sav'
  endif 
  print,save_file_name
  save,xt,xf,p2_full,filename='mask_combine_'+save_file_name+'.sav'
Endif

;************************************************************
;                    Mask Bit
;************************************************************
If keyword_Set(mask_bit) then begin
  run_0 = where(BEAMS eq 0,/null)
  If run_0 ne !Null then begin
    print,'**Run Beam 0**'
    ;------------------
    ;  Restore Beam 0 
    ;------------------
    restore, filename='mask_full_bit_'+save_file_name+'_beam0.sav'
     ;xt,xf,xsize,p2_full_bit restored
    BITARRAY_TO_BYTEARRAY, p2_full_bit, xsize, p2_full_0
    p2_full_all = temporary(p2_full_0) 
    
    ww = where(p2_full_all eq 0)  ;bad pixels
    n_ww = n_elements(ww)
    n_p2_all = n_elements(p2_full_all)
    p0=n_ww*100.d0/n_p2_all   ;percent of bad pixels
    print,'# polluted channels for 0 Beam = ',n_ww,' / ',n_p2_all
    print,'RFI Mask for all beams:       ',p0,' % -> masked out'  
    PRINT, 'Memory after Beam 0', (MEMORY(/HIGHWATER))/1d9 ,' Gb'
  endif 
    
  ;------------------
  ;  Restore Beam 1
  ;------------------
  run_1 = where(BEAMS eq 1,/null)
  If run_1 ne !Null then begin
    print,'**Run Beam 1**'
    restore, filename='mask_full_bit_'+save_file_name+'_beam1.sav'
      ;xt,xf,xsize,p2_full restored
    p2_full_bit_1 = temporary(p2_full_bit)
    BITARRAY_TO_BYTEARRAY, p2_full_bit_1, xsize, p2_full_1
    
    If run_0 eq !Null and run_1 ne !NUll then p2_full_all = temporary(p2_full_1)
    If run_0 ne !Null and run_1 ne !NUll then p2_full_all = temporary(p2_full_1)*temporary(p2_full_all)
   
    ww = where(p2_full_all eq 0)  ;bad pixels
    n_ww = n_elements(ww)
    n_p2_all = n_elements(p2_full_all)
    p0=n_ww*100.d0/n_p2_all   ;percent of bad pixels
    print,'# polluted channels for 0 and 1 Beams = ',n_ww,' / ',n_p2_all
    print,'RFI Mask for all beams:       ',p0,' % -> masked out'
    PRINT, 'Memory after Beam 1', (MEMORY(/HIGHWATER))/1d9 ,' Gb'
  endif
  
  ;------------------
  ;  Restore Beam 2
  ;------------------
  run_2 = where(BEAMS eq 2,/null)
  If run_2 ne !Null then begin
     print,'**Run Beam 2**'
     restore, filename='mask_full_bit_'+save_file_name+'_beam2.sav'
      ;xt,xf,xsize,p2_full restored
    p2_full_bit_2 = temporary(p2_full_bit)
    BITARRAY_TO_BYTEARRAY, p2_full_bit_2, xsize, p2_full_2
    
    If run_0 eq !Null and run_1 eq !NUll and run_2 ne !NULL then p2_full_all = temporary(p2_full_2)
    If run_0 ne !Null or run_1 ne !NUll then begin
      If run_2 ne !Null then p2_full_all = temporary(p2_full_all)*temporary(p2_full_2)
    Endif
   
    ww = where(p2_full_all eq 0)  ;bad pixels
    n_ww = n_elements(ww)
    n_p2_all = n_elements(p2_full_all)
    p0=n_ww*100.d0/n_p2_all   ;percent of bad pixels
    print,'# polluted channels for 0,1,2 Beams = ',n_ww,' / ',n_p2_all
    print,'RFI Mask for all beams:       ',p0,' % -> masked out'
    PRINT, 'Memory after Beam 2', (MEMORY(/HIGHWATER))/1d9 ,' Gb'
  endif 
  
  ;------------------
  ;  Restore Beam 3
  ;------------------
  run_3 = where(BEAMS eq 3,/null)
  If run_3 ne !Null then begin
    print,'**Run Beam 3**'
    restore, filename='mask_full_bit_'+save_file_name+'_beam3.sav'
     ;xt,xf,xsize,p2_full_bit restored
    p2_full_bit_3 = temporary(p2_full_bit)
    BITARRAY_TO_BYTEARRAY, p2_full_bit_3, xsize, p2_full_3
    
    If run_0 eq !Null and run_1 eq !NUll and run_2 eq !NULL and run_3 ne !NULL then p2_full_all = temporary(p2_full_3)
    If run_0 ne !Null or run_1 ne !NUll or run_2 ne !NUll then begin
      If run_3 ne !Null then p2_full_all = temporary(p2_full_all)*temporary(p2_full_3)
    Endif
   
  endif 
  
  ;---------------------------------
  ;          Combine 
  ;---------------------------------
    ww = where(p2_full_all eq 0)  ;bad pixels
    n_ww = n_elements(ww)
    n_p2_all = n_elements(p2_full_all) 
    p0=n_ww*100.d0/n_p2_all   ;percent of bad pixels
    print,'# polluted channels for all Beams = ',n_ww,' / ',n_p2_all
    print,'RFI Mask for all beams:       ',p0,' % -> masked out'

    PRINT, 'Memory after combination', (MEMORY(/HIGHWATER))/1d9 ,' Gb'
   BYTEARRAY_TO_BITARRAY, temporary(p2_full_all), xsize, p2_full_bit    
  save,xt,xf,xsize,p2_full_bit,filename='mask_combine_bit_'+save_file_name+'.sav'
Endif

If keyword_set(PS) then begin
  device,/close
  EXIT_PS
  ;cgfixps, save_file_name+'_maskcombine.ps'
endif
print,'SYSTIME=',SYSTIME(1) - st, ' sec'
PRINT, 'Memory required for combining masks: ', (MEMORY(/HIGHWATER))/1d9 ,' Gb'

end