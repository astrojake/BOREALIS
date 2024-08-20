;*****************************************************
;************* Read H5 DATA Jupiter Nancay ********************
;****************************************************
;          Name of dataset is Data 
;         Data in format (frequency, time) 
;;-------------------------------------
 pro READ_H5_DATA_Jupiter,filename,x=x,freq=freq,time=time,HEADER=HEADER;,tmin,tmax,fmin,fmax

  file_id             = h5f_open(filename+'.h5')
  dataset_id          = h5d_open(file_id,'Data')
  dataspace_id        = H5D_GET_SPACE(dataset_id)
  dimensions          = H5S_GET_SIMPLE_EXTENT_DIMS(dataspace_id)
  
  If keyword_set(HEADER) then begin 
    dataset_id2         = h5d_open(file_id,'Freq')
      dataset_id3         = h5d_open(file_id,'Time')
      dataspace_id2       = H5D_GET_SPACE(dataset_id2)
      dataspace_id3       = H5D_GET_SPACE(dataset_id3)
      freq                = h5d_read(dataset_id2)
      Time                = h5d_read(dataset_id3)
      
       H5D_CLOSE, dataset_id2
       H5S_CLOSE, dataspace_id2
  endif
  
  
  ;time                = h5a_read(h5a_open_name(dataset_id,'TIME') )
  ;freq                = h5a_read(h5a_open_name(dataset_id,'FREQ') ) 
  nf = dimensions[0]
  nt = dimensions[1]
  
;  u = where(time ge tmin and time lt tmax) 
;  t       = time(u)
;  v = where(center_freq ge fmin and center_freq le fmax)
;
;start   = [min(v),min(u)]
;count    = [n_elements(u),n_elements(v)]
;
;  h5s_select_hyperslab,dataspace_id,start,count,/RESET
;  memory_space_id   = h5s_create_simple(count)
  x                 = h5d_read(dataset_id)
  x                 = transpose(temporary(x))
  
  
 print, min(freq)

  
  help,x
  
 ; H5S_CLOSE, memory_space_id
  H5S_CLOSE, dataspace_id

  H5D_CLOSE, dataset_id
  H5F_CLOSE, file_id
end