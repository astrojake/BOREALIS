pro DIR_SEPARATOR, sep
  case strlowcase(!VERSION.OS_FAMILY) of
	'unix': sep='/'
	'windows': sep='\'
  endcase
return
end
