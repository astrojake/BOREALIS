;Find the moments of the data
;
pro find_moments,data,data_out=data_out
  data_out = dblarr(5)
  data_out[0]  = mean(data,/nan,/double)
  data_out[1]  = median(data,/double)
  data_out[2]  = max(data,/nan)
  data_out[3]  = min(data,/nan)
  data_out[4] =  stddev(data,/nan,/double)
end