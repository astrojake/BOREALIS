;Create filename list for the combined dates program
pro create_filenames,file_root,dates,dates_steps,N_time,patrol_value,le_sig_value,S,$
                     RFI=RFI,patrol=patrol,le_sig=le_sig,sum=sum,pex=pex,dates_filenames=dates_filenames,$
                     save_file_name=save_file_name,COMPLEX_NAME=COMPLEX_NAME
                        
n               = n_elements(dates)
DATES_FILENAMES = strarr(n)

for i=0,n-1 do begin
  
  save_create,N_time,dates_steps[i],patrol_value,le_sig_value,0,RFI=RFI,patrol=patrol,le_sig=le_sig,sum=sum,pex=pex,$
    save_file_name=save_file_name,/NOBEAM,COMPLEX_NAME=COMPLEX_NAME
    
  dates_filenames[i] = dates[i]+'_pol'+strtrim(string(S),1)+'_'+save_file_name
endfor

;find save_file_name for combine file
l = STRLEN(save_file_name)
l2 = STRPOS(save_file_name, 'RFI')
save_file_name = STRMID(save_file_name, l2, l)

end 