
PRO ps_to_pdf, file, local_path
; convert 1 .ps file to a .pdf file
; v. 13 Mar 2015

DIR_SEPARATOR, sep

shortname = strmid(DELPATH(file),0,strlen(DELPATH(file))-4)
fsplit=strsplit(file,sep)
addfolder=strmid(file,fsplit[n_elements(fsplit)-2],(fsplit[n_elements(fsplit)-1]-fsplit[n_elements(fsplit)-2]))

command='ps2pdf '+local_path+addfolder+shortname+sep+shortname+'.ps'
spawn, command

command='del '+local_path+addfolder+shortname+sep+shortname+'.ps'
spawn, command

return
END