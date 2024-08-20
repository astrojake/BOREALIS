pro EXIT_PS
  !p.multi=0
  case strlowcase(!VERSION.OS_FAMILY) of
	'x': set_plot,'X'
	'unix': set_plot,'X'
	'windows': set_plot, 'WIN'
  endcase
return
end
