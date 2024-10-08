pro latex_pdf,filename,MinF,MaxF

 ;create latex and pdf 
  date = strmid(filename,0,7)
  tex_output,filename,date,MinF,MaxF 
  run_string = 'pdflatex '+filename+'_Out.tex'
  spawn,run_string ;create pdf from tex
  spawn,'rm '+filename+'*.aux' ;remove tex files 
  spawn,'rm '+filename+'*.dvi' 
  spawn,'rm '+filename+'*.log' 
     
  out_string='gs -sDEVICE=pdfwrite -dPDFSETTINGS=/prepress -dNOPAUSE -dQUIET -dBATCH -sOutputFile='+filename+'_Output.pdf'+' '+filename+'_Out.pdf'
  spawn,out_string ;compress pdf 
  string_rm = 'rm '+filename+'_Out.pdf'  
  ;spawn,string_rm
end
