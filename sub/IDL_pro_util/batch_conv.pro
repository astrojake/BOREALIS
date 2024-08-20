; ------------------------------------------------
  pro BATCH_CONV, conv
; ------------------------------------------------

  path=extpath(conv)
  conv_out=conv+'_'
  openr, 1, conv
  openw, 2, conv_out
  on_ioerror,suite
  buf=''

encore:
  err=0
  readf, 1, buf
  f=delpath(buf)
;  printf,2, 'ps2pdfwr '+f
   printf,2, 'ps2pdf -dPDFSETTINGS=/ebook '+f
; printf,2, 'rm '+f
  goto,encore

suite:
  close,1 & close,2
return
end
