;all inputs are strings
;eg. filename='L570725_pol0_10779sec_RFI_patrol5.5_pex_lesig3.5_sumJupiter_reduce1000_50MHz14.7to62.4MHz_thresT0.10_thresF0.10_postprocessing'
pro tex_output,filename,date,MinF,MaxF,Q3=Q3

openw,lun,filename+'_Out.tex',/get_lun

printf,lun,'\documentclass{article}'
printf,lun,'\usepackage{graphicx}'
printf,lun,'\usepackage[letterpaper,margin=2in]{geometry}'
printf,lun,'\newcommand\filename{{'+filename+'}}'
printf,lun,'\setlength{\tabcolsep}{1.6em}'
printf,lun,'\begin{document}'
printf,lun,'\begin{figure}'
printf,lun,'\centering'
printf,lun,'\begin{tabular}{lll}'
printf,lun,'\textbf{a.)} & \textbf{b.)} &\textbf{c.)}\\'
printf,lun,'\includegraphics[width=0.39\textwidth,trim=4cm 1cm 1cm 2.2cm,page=1]{\filename.pdf} &'
printf,lun,'\includegraphics[width=0.39\textwidth,trim=4cm 1cm 1cm 2.2cm,page=3]{\filename.pdf} &'
printf,lun,'\includegraphics[width=0.39\textwidth,trim=4cm 1cm 1cm 2.2cm,page=11]{\filename.pdf} \\'
printf,lun,'\textbf{d.)} & \textbf{e.)} & \textbf{f.)}\\'
printf,lun,'\includegraphics[width=0.39\textwidth,trim=4cm 1cm 1cm 2.2cm,page=2]{\filename.pdf} &'
printf,lun,'\includegraphics[width=0.39\textwidth,trim=4cm 1cm 1cm 2.2cm,page=4]{\filename.pdf} &
printf,lun,'\includegraphics[width=0.39\textwidth,trim=4cm 1cm 1cm 2.2cm,page=12]{\filename.pdf} \\'
printf,lun,'\textbf{g.)} & \textbf{h.)} & \textbf{i.)}\\'
printf,lun,'\includegraphics[width=0.39\textwidth,trim=4cm 1cm 1cm 2.2cm,page=7]{\filename.pdf} &'
printf,lun,'\includegraphics[width=0.39\textwidth,trim=4cm 1cm 1cm 2.2cm,page=5]{\filename.pdf} &'
printf,lun,'\includegraphics[width=0.39\textwidth,trim=4cm 1cm 1cm 2.2cm,page=9]{\filename.pdf} \\'
printf,lun,'\textbf{j.)} & \textbf{k.)} & \textbf{l.)}\\'
printf,lun,'\includegraphics[width=0.39\textwidth,trim=4cm 1cm 1cm 2.2cm,page=8]{\filename.pdf} &'
printf,lun,'\includegraphics[width=0.39\textwidth,trim=4cm 1cm 1cm 2.2cm,page=6]{\filename.pdf} &'
printf,lun,'\includegraphics[width=0.39\textwidth,trim=4cm 1cm 1cm 2.2cm,page=10]{\filename.pdf} \\'
printf,lun,'\end{tabular}'   
printf,lun,'\caption{'+date+' from '+Minf+' to '+Maxf+'. \textbf{a.)} ON-beam dynamic spectrum. \textbf{b.)} Q1a (Time-series). \textbf{c.)} Q2 vs time. \textbf{d.)} OFF-beam dyanmic spectrum. \textbf{e.)} Q1a (ON-OFF). \textbf{f.)} Q2 vs tiem (On-OFF). \textbf{g.)} On-beam high-pass filtered dynamic spectrum. \textbf{h.)} Q1b (integrated-spectrum). \textbf{i.)} Q2 (Before El Correction). \textbf{j.)} OFF-beam high-pass filtered dynamic spectrum. \textbf{k.)} Q1b (ON-OFF). \textbf{l.)} Q2 (After El Correction).}'
printf,lun,'\end{figure}'

;If keyword_set(Q3) then begin
;---------------- Q4-------------------------------------------------------------------
 printf,lun,'\begin{figure}'
 printf,lun,'\begin{center}'
 printf,lun,'\begin{tabular}{lll}'
 printf,lun,'\textbf{a.)} & \textbf{b.)} &\textbf{c.)}\\'
 printf,lun,'\includegraphics[width=0.39\textwidth,trim=4cm 1cm 1cm 2.2cm,page=13]{\filename.pdf} &'
 printf,lun,'\includegraphics[width=0.39\textwidth,trim=4cm 1cm 1cm 2.2cm,page=17]{\filename.pdf} &'
 printf,lun,'\includegraphics[width=0.39\textwidth,trim=4cm 1cm 1cm 2.2cm,page=21]{\filename.pdf} \\'
printf,lun,'\textbf{d.)} & \textbf{e.)} & \textbf{f.)}\\'
 printf,lun,'\includegraphics[width=0.39\textwidth,trim=4cm 1cm 1cm 2.2cm,page=14]{\filename.pdf} &'
 printf,lun,'\includegraphics[width=0.39\textwidth,trim=4cm 1cm 1cm 2.2cm,page=18]{\filename.pdf} &
 printf,lun,'\includegraphics[width=0.39\textwidth,trim=4cm 1cm 1cm 2.2cm,page=22]{\filename.pdf} \\'
 printf,lun,'\textbf{g.)} & \textbf{h.)} & \textbf{i.)}\\'
 printf,lun,'\includegraphics[width=0.39\textwidth,trim=4cm 1cm 1cm 2.2cm,page=15]{\filename.pdf} &'
 printf,lun,'\includegraphics[width=0.39\textwidth,trim=4cm 1cm 1cm 2.2cm,page=19]{\filename.pdf} &'
 printf,lun,'\includegraphics[width=0.39\textwidth,trim=4cm 1cm 1cm 2.2cm,page=23]{\filename.pdf} \\'
 printf,lun,'\textbf{j.)} & \textbf{k.)} & \textbf{l.)}\\'
 printf,lun,'\includegraphics[width=0.39\textwidth,trim=4cm 1cm 1cm 2.2cm,page=16]{\filename.pdf} &'
 printf,lun,'\includegraphics[width=0.39\textwidth,trim=4cm 1cm 1cm 2.2cm,page=20]{\filename.pdf} &'
 printf,lun,'\includegraphics[width=0.39\textwidth,trim=4cm 1cm 1cm 2.2cm,page=24]{\filename.pdf} \\'
 printf,lun,'\end{tabular}'   
 printf,lun,'\caption{'+date+' from '+Minf+' to '+Maxf+'. \textbf{a.)} Q4a (number of peaks) \textbf{b.)} Q4c (peak asymmetry) \textbf{c.)} Q4e (peak offset) \textbf{d.)} Q4a (ON-OFF)  \textbf{e.)} Q4c (ON-OFF) \textbf{f.)} Q4e (ON-OFFF)  \textbf{g.)} Q4b (power of peaks) \textbf{h.)} Q4d (power asymmetry)  \textbf{i.)} Q4f (power offset) \textbf{j.)} Q4b (ON-OFF) \textbf{k.)} Q4d (ON-OFF) \textbf{l.)} Q4f (ON-OFF) }'
 printf,lun,'\end{center}'
 printf,lun,'\end{figure}'
;endif 

printf,lun,'\end{document}'

FREE_LUN, lun ;close text file
end

