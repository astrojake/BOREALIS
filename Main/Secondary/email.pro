;Send Email
pro email,email_address,string,send=send 

 If keyword_set(send) then begin
  binTime = BIN_DATE(systime(/ut))  ;
  timestamp_string = TIMESTAMP(YEAR = binTime[0], MONTH = binTime[1],$
    DAY = binTime[2], HOUR = binTime[3], MINUTE = binTime[4], SECOND =binTime[5])  ;2012-09-04T11:25:15Z
  string2 = 'mail -s "'+string+'at '+timestamp_string+'"'
  spawn, string2+' '+email_address+'< /dev/null'
 Endif
end

