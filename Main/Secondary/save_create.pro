;Create save file name for lofar pipeline
;
;Not general 
;
;Includes RUN_JUPiTER

pro save_create,N_time,steps,patrol_value,le_sig_value,Beam,RFI=RFI,patrol=patrol,le_sig=le_sig,sum=sum,pex=pex,run_JUPITER=run_JUPITER,$
                save_file_name=save_file_name,NOBEAM=NOBEAM,LIST=LIST,COMPLEX_NAME=COMPLEX_NAME
               
   ;max time
   total_time = steps*N_time 
  
   If not(keyword_set(COMPLEX_NAME)) then begin
     save_file_name = 'beam'+strtrim(string(cgnumber_formatter(beam,decimals=0)),1)
   Endif
  
  ;*********************
  ;Filenames
  ;*********************
  If keyword_set(COMPLEX_NAME) then begin
    ;Patrol
    If keyword_set(patrol) and not(keyword_set(pex)) and not(keyword_set(le_sig)) and not(keyword_set(sum)) then begin
      save_file_name =  strtrim(string(fix(total_time)),1)+'sec_RFI_patrol'+strtrim(string(cgnumber_formatter(patrol_value,decimals=1)),1)+'_beam'$
        +strtrim(string(cgnumber_formatter(beam,decimals=0)),1)
    Endif
  
    ;Le sig
    If not(keyword_set(patrol)) and not(keyword_set(pex)) and keyword_set(le_sig) and not(keyword_set(sum)) then begin
      save_file_name =  strtrim(string(fix(total_time)),1)+'sec_RFI_lesig'+strtrim(string(cgnumber_formatter(le_sig_value,decimals=1)),1)+'_beam'$
        +strtrim(string(cgnumber_formatter(beam,decimals=0)),1)
    Endif
  
    ;sum
    If not(keyword_set(patrol)) and not(keyword_set(pex)) and not(keyword_set(le_sig)) and keyword_set(sum) then begin
      save_file_name =  strtrim(string(fix(total_time)),1)+'sec_RFI_sum_beam'$
        +strtrim(string(cgnumber_formatter(beam,decimals=0)),1)
    Endif
  
    ;Patrol and Pex
    If keyword_set(patrol) and keyword_set(pex) and not(keyword_set(le_sig)) and not(keyword_set(sum)) then begin
      save_file_name =  strtrim(string(fix(total_time)),1)+'sec_RFI_patrol'+strtrim(string(cgnumber_formatter(patrol_value,decimals=1)),1)+'_pex_beam'$
        +strtrim(string(cgnumber_formatter(beam,decimals=0)),1)
    Endif
  
    ;Patrol, Pex, Le Sig
    If keyword_set(patrol) and keyword_set(pex) and keyword_set(le_sig) and not(keyword_set(sum)) then begin
      save_file_name =  strtrim(string(fix(total_time)),1)+'sec_RFI_patrol'+strtrim(string(cgnumber_formatter(patrol_value,decimals=1)),1)+'_pex_lesig'+$
        strtrim(string(cgnumber_formatter(le_sig_value,decimals=1)),1)+'_beam'$
        +strtrim(string(cgnumber_formatter(beam,decimals=0)),1)
    Endif
  
    ;Pex, Le Sig
    If not(keyword_set(patrol)) and keyword_set(pex) and keyword_set(le_sig) and not(keyword_set(sum)) then begin
      save_file_name =  strtrim(string(fix(total_time)),1)+'sec_pex_lesig'+$
        strtrim(string(cgnumber_formatter(le_sig_value,decimals=1)),1)+'_beam'$
        +strtrim(string(cgnumber_formatter(beam,decimals=0)),1)
    Endif
  
    ;Patrol, Pex, Le sig, sum
    If keyword_set(patrol) and keyword_set(pex) and keyword_set(le_sig) and keyword_set(sum) then begin
      save_file_name =  strtrim(string(fix(total_time)),1)+'sec_RFI_patrol'+strtrim(string(cgnumber_formatter(patrol_value,decimals=1)),1)+'_pex_lesig'+$
        strtrim(string(cgnumber_formatter(le_sig_value,decimals=1)),1)+'_sum_beam'$
        +strtrim(string(cgnumber_formatter(beam,decimals=0)),1)
    Endif
  
    ;patrol, Pex, Sum
    If keyword_set(patrol) and keyword_set(pex) and not(keyword_set(le_sig)) and keyword_set(sum) then begin
      save_file_name =  strtrim(string(fix(total_time)),1)+'sec_RFI_patrol'+strtrim(string(cgnumber_formatter(patrol_value,decimals=1)),1)+'_pex_sum_beam'$
        +strtrim(string(cgnumber_formatter(beam,decimals=0)),1)
    Endif
  
    ;patrol and le sig
    If keyword_set(patrol) and not(keyword_set(pex)) and keyword_set(le_sig) and not(keyword_set(sum)) then begin
      save_file_name =  strtrim(string(fix(total_time)),1)+'sec_RFI_patrol'+strtrim(string(cgnumber_formatter(patrol_value,decimals=1)),1)+'_lesig'+$
        strtrim(string(cgnumber_formatter(le_sig_value,decimals=1)),1)+'_beam'$
        +strtrim(string(cgnumber_formatter(beam,decimals=0)),1)
    Endif
  
    ;Patrol, le sig, sum
    If keyword_set(patrol) and not(keyword_set(pex)) and keyword_set(le_sig) and keyword_set(sum) then begin
      save_file_name =  strtrim(string(fix(total_time)),1)+'sec_RFI_patrol'+strtrim(string(cgnumber_formatter(patrol_value,decimals=1)),1)+'_lesig'+$
        strtrim(string(cgnumber_formatter(le_sig_value,decimals=1)),1)+'_sum_beam'$
        +strtrim(string(cgnumber_formatter(beam,decimals=0)),1)
    Endif
  
    ;patrol and sum
    If keyword_set(patrol) and not(keyword_set(pex)) and not(keyword_set(le_sig)) and keyword_set(sum) then begin
      save_file_name =  strtrim(string(fix(total_time)),1)+'sec_RFI_patrol'+strtrim(string(cgnumber_formatter(patrol_value,decimals=1)),1)+'_sum_beam'$
        +strtrim(string(cgnumber_formatter(beam,decimals=0)),1)
    Endif
  endif 
  
  If keyword_set(run_JUPITER) then begin
    input_file_J = 'Input_Jupiter.dat'
    read_jupiter_input,input_file_J,scale=scale,fmin_add=fmin_add
    l = STRLEN(save_file_name)
    save_file_name = STRMID(save_file_name, 0, l-6) ;get rid of beam  
    save_file_name = save_file_name+'Jupiter_reduce'+strtrim(string(ulong(scale)),1)+'_'$
      +strtrim(string(fix(fmin_add)),1)+'MHz_beam'$
      +strtrim(string(cgnumber_formatter(beam,decimals=0)),1)
  Endif

  If keyword_set(NOBEAM) then begin
    l = STRLEN(save_file_name)
    save_file_name = STRMID(save_file_name, 0, l-6) ;get rid of beam
  Endif
  
  If keyword_set(LIST) then begin
    l = STRLEN(save_file_name)
    l2 = STRPOS(save_file_name, 'RFI')
    save_file_name = STRMID(save_file_name, l2, l-6) ;no timing and no beam
  Endif
end