;;Batch run of LOFAR_BF_PIPELINEV3
;Inputs: filename   (this is the filename containing the PP of interest 
;
;To run 
; .r LOFAR_BF_PipelineV3
; .r batch_PPrun.pro
; filename = 'Target_TauBoo_dir.dat'
pro batch_PPrun,filename
readcol,filename,dir,format='(A)'

n_dir = n_elements(dir)

for i=0,n_dir-1 do begin 
 cd,dir[i]  
 ;spawn,'rm log*'
 ;RESOLVE_ROUTINE, 'LOFAR_BF_PipelineV3.pro' 
 LOFAR_BF_PipelineV3,'Inputv4.dat'
endfor

end