;Creat panels for post-processing setup
;Input:   filename
;Output:  filename_Output.pdf
pro panel_pdf,filename,Q3=Q3,REMOVE=REMOVE

;------------
;Reduce size
;------------
;sizes = prepress, screen 
;string='gs -sDEVICE=pdfwrite -dPDFSETTINGS=/prepress -dNOPAUSE -dQUIET -dBATCH -sOutputFile='+filename+'_Output.pdf'+' '+filename+'_Out.pdf'
string1='gs -sDEVICE=pdfwrite -dCompatibilityLevel=1.4 -dPDFSETTINGS=/screen -dNOPAUSE -dQUIET -dBATCH -sOutputFile='+filename+'_Output.pdf'+' '+filename+'.pdf'
spawn,string1 ;compress pdf 

;update filename
filename1 = +filename+'_Output.pdf'

;----------------
;Rearrage pages:
;-----------------
;page_string = '1,3,11,2,4,12,7,5,9,8,6,10,13,17,21,14,18,22,15,19,23,16,20,24'
If not(keyword_set(Q3)) then page_string =  '1,3,11,2,4,12,9,5,7,10,6,8,13,17,21,14,18,22,15,19,23,16,20,24'
If keyword_set(Q3)      then page_string =  '1,3,11,2,4,12,9,5,7,10,6,8,13,17,21,14,18,22,15,19,23,16,20,24,25,29,33,26,30,34,27,31,35,28,32,36'
string2='pdfjam '+filename1+' "'+page_string+'" -o '+filename1
spawn,string2

;-----------
;crop pages 
;-----------
string3='pdfcrop '+filename1+' '+filename1
spawn,string3 

;--------------
;create panels 
;--------------
string4 = 'pdfnup --nup 3x4 --no-landscape '+filename1+' -o '+filename1
spawn, string4 

;-------------
;crop page 
;-------------
spawn,string3 

;remove pdf
string5 = 'rm '+filename+'.pdf'
If keyword_set(REMOVE) then spawn,string5 

end
