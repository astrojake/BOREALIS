;***************************************************************
;      Read Lofar data
;      This read program is only for polarization (QUV)
;**************************************************************      
;Uses: hdf5 programs
;***************************************************************     
;Author: Jake Turner (UVA)
;Things that you can change: 
;   start_freq
;   end_freq 
;   start_time 
;   end_time
;   file_root
;   date

;initial inputs
st= SYSTIME(1)
j_name = 'time_'+strtrim(string(st),1)+'.dat'
journal,j_name

; check if all parameters are reasonable
;Start Freq: ;End Freq: 
start_freq = 26.07d  ;MHz
end_freq   = 73.7d   ;MHz

;start time and end time
start_time = 100.  ; works good

;N_time = 1919.0d  ;amount of time for each plot
N_time = 10.d 
steps = 1.0d    ;number of steps to take of N_time 
;steps = 15.0d ; 1919 all data

total_t = N_time*steps  ;total time requested

  ;------------
  ;File Format
  ;-----------
  ;L248553_SAP000_B000_S0_P000_bf.h5'
  ;L is the date 
  ;B000 - B003: 4 beams: 4 coherent (55cnc, Pulsar, empty, source)
  ;S0-S3:  4 polizations: IQUV
  
  ;root directory for raw data
  file_root = '/data/jake.turner/exoplanet/LC5_DDT_002/'
 ; file_root = '/databf/lofar/EXO/LC5_DDT_002/'
  ;file_root = '/databf/lofar/EXO/LC2_018/'
  ;---------
  ;Read Dates
  ;----------
  date = ''
  date = 'L429868'  ;on /data
  
  
  ;These strings are for Beam and Pol for loops 
  ;strings for information
  B_start = 2    ;Beam Start 
  B_end   = 2    ;Beam  End
  
  ;You always have to start with I because you need it for the QUV 
  
  ;polarization to run
  ;QU       S_run = 1; L = sqrt(Q^2 + U^2)
  ;V        S_run = 2;
  ;QU and V S_run = 3; 
 
  S_run   = 1    ;polarization Start
  ;S_run   = 2 
 
  ;**************************************
  ;Don't need to change anything below here
  ;*********************************************
  B_start_string = strtrim(string(b_start),1)
;  S_start_string = strtrim(string(S_start),1)
  B_end_string = strtrim(string(b_end),1)
 ; S_start_string = strtrim(string(S_start),1)
;  S_end_string = strtrim(string(S_end),1)
  
  S_run_string = strtrim(string(s_run),1)
  
  ;-----------------------------
  ;       Create file name 
  ;-----------------------------
for i=0L,0 do begin  ;loop for dates
  
  ;----------------------------------------------------------
  ;              Setup Important Info for Data
  ;----------------------------------------------------------
  ;file name creation
  ; file = file_root+date[i]+'/'+date[i]+'_SAP000_B00'+j_string+$
  ;       '_S'+k_string+'_P000_bf.h5'
  
  ;file to read in headers and arrays
  file = file_root+date[i]+'/'+date[i]+'_SAP000_B00'+b_start_string+'_S0'+'_P000_bf.h5'
  ;file = '/data/jake.turner/raw/L248553_link/'+date[i]+'_SAP000_B000_S0_P000_bf.h5'
  print, file
  ;--------------------------------------------------
  ;--------------------------------------------------
  ;         READ Files, groups, and data bound
  ;--------------------------------------------------
  ; Check file, group, and dataset names with: hdfview
  ;--------------------------------------------------
  ;pay attention to hdf, h5g, and h5d
  file_id     = h5f_open(file)
  sub_file_id = h5g_open(file_id,'SUB_ARRAY_POINTING_000')
  group_id    = h5g_open(file_id,'SUB_ARRAY_POINTING_000/BEAM_00'+b_start_string)
  dataset_id  = h5d_open(group_id,'STOKES_0')
  
  dataspace_id = H5D_GET_SPACE(dataset_id)
  dimensions   = H5S_GET_SIMPLE_EXTENT_DIMS(dataspace_id)
  
  
  ;--------------------------------------------------
  ;--------------------------------------------------
  ;             Get Useful Information
  ;--------------------------------------------------
  ;--------------------------------------------------
  
  ;-------------------------
  ;Bandwidth and Freq Range
  ;-------------------------
  Bandwidth      =  h5a_read( h5a_open_name(file_id,'BANDWIDTH' ))  ; units in MHz
  Bandwidth_unit =  h5a_read( h5a_open_name(file_id,'BANDWIDTH_UNIT' )) ; should be MHz
  Min_Freq       = h5a_read( h5a_open_name(file_id,'OBSERVATION_FREQUENCY_MIN' ))
  Max_Freq       = h5a_read( h5a_open_name(file_id,'OBSERVATION_FREQUENCY_MAX' ))
  Freq_Unit      = h5a_read( h5a_open_name(file_id,'OBSERVATION_FREQUENCY_UNIT' )) ;MHz
  Channel_Width  = h5a_read( h5a_open_name(group_id,'CHANNEL_WIDTH') ) ; Hz
  Channel_Width_unit  = h5a_read( h5a_open_name(group_id,'CHANNEL_WIDTH_UNIT') ) ; Hz
  
  ;--------
  ;  Time
  ;-------
  START_MJD          =  h5a_read( h5a_open_name(sub_file_id,'EXPTIME_START_MJD' ))
  END_MJD            =  h5a_read( h5a_open_name(sub_file_id,'EXPTIME_END_MJD' ))
  TOT_INT_TIME       = h5a_read( h5a_open_name(sub_file_id,'TOTAL_INTEGRATION_TIME'));units in secs
  SAMPLING_TIME      = h5a_read( h5a_open_name(group_id,'SAMPLING_TIME') )
  NOF_SAMPLES        = h5a_read(h5a_open_name(dataset_id,'NOF_SAMPLES') )
  
  ;-------
  ; Bands
  ;-------
  NOF_SUBBANDS       = h5a_read(h5a_open_name(dataset_id,'NOF_SUBBANDS') )
  ;NOF_SUBBANDS_file   =reform(NOF_SUBBANDS)
   print,strtrim(NOF_SUBBANDS,2)+' subbands in this file'
  
  ;----------------------------------------------
  ;----------------------------------------------
  ;               Setup Frequency Array 
  ;----------------------------------------------
  ;----------------------------------------------
  ;Master Frequency Array
  
  ;-------------------------------
  ;New group id for Frequcy array
  ;------------------------------
  group_id_1    = h5g_open(file_id,'SUB_ARRAY_POINTING_000/BEAM_00'+b_start_string+'/COORDINATES/COORDINATE_1')
  
  ;---------------
  ;Find Frequency
  ;---------------
  center_freq = h5a_read( h5a_open_name(group_id_1,'AXIS_VALUES_WORLD' )) ;Hz
  center_freq = center_freq/1.0d6 ;MHz
  
  ;close groups!
  H5S_CLOSE, dataspace_id
  H5D_CLOSE, dataset_id
  h5g_close, group_id
  h5g_close, sub_file_id
  H5F_CLOSE, file_id
  ;------------------------------
  ; Edit Frequency for user input
  ;-----------------------------
  w        = where(CENTER_FREQ ge start_freq and CENTER_FREQ le end_freq)
  num_bands = n_elements(w) 
  
  ;----------------------------------------------
  ;----------------------------------------------
  ;               Setup Time Array
  ;----------------------------------------------
  ;----------------------------------------------

  
  ;-----------
  ;Find Time
  ;-----------
   time = dblarr(nof_samples)
  for ii=1d,nof_samples-1 do begin
    time[ii] = sampling_time*ii
  endfor
  ii = 1.d ; reset ii

  print, 'Max Time (s)' ,max(time)
  
    ;--------------------------------
    ;Update Frequency with user input
    ;--------------------------------
    center_freq_user = center_freq(w)
    print, 'Min User Frequency (MHz): ', min(center_freq_user)
    print, 'Max User Frequency (MHz): ' , max(center_freq_user)

    ;----------------------------------------------------------
    ;        Start loops for Beam and Polarization
    ;----------------------------------------------------------
    ;loops  beams, and polarization
for j=b_start,b_end do begin  ; loop for Beam: B000-B0003
      ;for k=0,0 do begin  ; loop for Polarization:S0-S3

   ;QU = 1; L = sqrt(Q^2 + U^2)
   ;V  = 2
   ;QU and V = 3
   If S_run eq 1 then begin
    k_string = 'QU' 
   Endif
   
   If S_run eq 2 then begin
     k_string = 'V'
   Endif
   
   If S_run eq 3 then begin
     k_string = 'QUV'
   Endif
   

        j_string = strtrim(string(j),1)
     ;   k_string = strtrim(string(k),1)
        ;----------------------------------------------------------
        ;        Read Data for specific Beam and Polarization
        ;----------------------------------------------------------
        ;file name creation
        set_plot,'PS'
        device,filename=file_root+date[i]+'/'+date[i]+'_SAP000_B00'+j_string+$
          '_S'+k_string+'.ps'

    ;Time loop 
 for kkk=0L,steps-1 do begin
  
              print, 'Beam:' ,j_string
            ;  print, 'Pol:', k_string
              print, 'Time step:',kkk
        
              end_time =  start_time + N_time
              print, 'start time', start_time
              print, 'end_Time: ', end_time  
          
          ;------------------------------
          ; Edit Time for user input
          ;-----------------------------
          u        = where(time ge start_time and time  le end_time)
          num_time = n_elements(u) 
          
            ;--------------------------------
            ;Update Frequency with user input
            ;--------------------------------
            time_user = time(u)
            print, 'Min User Time (s): ', min(time_user)
            print, 'Max User Time (s): ' , max(time_user)
          
         
              ;--------------------------------------------------
              ;--------------------------------------------------
              ;         READ Files, groups, and data bound
              ;--------------------------------------------------
              ; Check file, group, and dataset names with: hdfview
              ;--------------------------------------------------
              ;pay attention to hdf, h5g, and h5d
              
              ;we always need to read in I data for normalization
              file_I = file_root+date[i]+'/'+date[i]+'_SAP000_B00'+j_string+$
                '_S0'+'_P000_bf.h5'
              print, file
              
              
              file_id     = h5f_open(file_I)
              group_id    = h5g_open(file_id,'SUB_ARRAY_POINTING_000/BEAM_00'+j_string)
              dataset_id  = h5d_open(group_id,'STOKES_0')
              
              ; Open up the dataspace
              dataspace_id = H5D_GET_SPACE(dataset_id)
              dimensions   = H5S_GET_SIMPLE_EXTENT_DIMS(dataspace_id)
              print,'Orignal Size of Data:', dimensions[0],dimensions[1]
              
              ;-----------------------------------------------
              ;          Now choose our area of interest
              ; What elements of the data array do you want? 
              ;-----------------------------------------------
              start = [min(w),min(u)]      ;w and u are the array inputs
              count = [num_bands,num_time]
              
              ;Select hyperslab
              ; Be sure to use /RESET to turn off all other selected elements.
              H5S_SELECT_HYPERSLAB, dataspace_id, start, count, /RESET
              
              ; Create a simple dataspace to hold the result. Would have problems if
              ; not
              memory_space_id = H5S_CREATE_SIMPLE(count)
              
              ;--------------
              ;Read the data
              ;--------------
              data_I   = H5D_READ(dataset_id, FILE_SPACE=dataspace_id, $
                       MEMORY_SPACE=memory_space_id)
              ;data(f,t) 
              ;data   = H5D_READ(dataset_id)
              size_data = size(data_I)
              print,'Size of Data Now:', size_data[1], size_data[2] 
              print,'-------------------------'
              print,'----Done Reading File----'
              print,'-------------------------'
              
              ;-------------------
              ;Transpose the data
              ;-------------------
              data_I = transpose(temporary(data_I))
               ;data_new(t,f); need this format for SPDYNPS
               
              H5S_CLOSE, dataspace_id
              H5D_CLOSE, dataset_id
              h5g_close, group_id
              H5F_CLOSE, file_id
              
   ;***************************************************************
   ;                  Read Polarization Data
   ;***************************************************************
          If S_run eq 1 then begin
              k_string = 'QU'
          Endif

          If S_run eq 2 then begin
            k_string = 'V'
          Endif

          If S_run eq 3 then begin
            k_string = 'QUV'
        Endif
    
    ;Linear Polarization need to read both Q and U and combine    
    If S_run eq 1 then begin
        ;***************************************
        ;         Read Q 
        ;***************************************         
              ;we always need to read in I data for normalization
              file_Q = file_root+date[i]+'/'+date[i]+'_SAP000_B00'+j_string+$
                '_S1'+'_P000_bf.h5'
              print, file


              file_id     = h5f_open(file_Q)
              group_id    = h5g_open(file_id,'SUB_ARRAY_POINTING_000/BEAM_00'+j_string)
              dataset_id  = h5d_open(group_id,'STOKES_1')   ;!!!! Stokes 1

              ; Open up the dataspace
              dataspace_id = H5D_GET_SPACE(dataset_id)
              dimensions   = H5S_GET_SIMPLE_EXTENT_DIMS(dataspace_id)
              print,'Orignal Size of Data:', dimensions[0],dimensions[1]

              ;-----------------------------------------------
              ;          Now choose our area of interest
              ; What elements of the data array do you want?
              ;-----------------------------------------------
              start = [min(w),min(u)]      ;w and u are the array inputs
              count = [num_bands,num_time]

              ;Select hyperslab
              ; Be sure to use /RESET to turn off all other selected elements.
              H5S_SELECT_HYPERSLAB, dataspace_id, start, count, /RESET

              ; Create a simple dataspace to hold the result. Would have problems if
              ; not
              memory_space_id = H5S_CREATE_SIMPLE(count)

              ;--------------
              ;Read the data
              ;--------------
              data_Q   = H5D_READ(dataset_id, FILE_SPACE=dataspace_id, $
                MEMORY_SPACE=memory_space_id)
              ;data(f,t)
              ;data   = H5D_READ(dataset_id)
              size_data = size(data_Q)
              print,'Size of Data Now:', size_data[1], size_data[2]
              print,'-------------------------'
              print,'----Done Reading File----'
              print,'-------------------------'

              ;-------------------
              ;Transpose the data
              ;-------------------
              data_Q = transpose(temporary(data_Q))
              ;data_new(t,f); need this format for SPDYNPS

              H5S_CLOSE, dataspace_id
              H5D_CLOSE, dataset_id
              h5g_close, group_id
              H5F_CLOSE, file_id
              
        ;***************************************
       ;         Read U
        ;***************************************
              ;we always need to read in I data for normalization
              file_U = file_root+date[i]+'/'+date[i]+'_SAP000_B00'+j_string+$
                '_S2'+'_P000_bf.h5'
              print, file


              file_id     = h5f_open(file_U)
              group_id    = h5g_open(file_id,'SUB_ARRAY_POINTING_000/BEAM_00'+j_string)
              dataset_id  = h5d_open(group_id,'STOKES_2')  ;!!!! Stokes 2

              ; Open up the dataspace
              dataspace_id = H5D_GET_SPACE(dataset_id)
              dimensions   = H5S_GET_SIMPLE_EXTENT_DIMS(dataspace_id)
              print,'Orignal Size of Data:', dimensions[0],dimensions[1]

              ;-----------------------------------------------
              ;          Now choose our area of interest
              ; What elements of the data array do you want?
              ;-----------------------------------------------
              start = [min(w),min(u)]      ;w and u are the array inputs
              count = [num_bands,num_time]

              ;Select hyperslab
              ; Be sure to use /RESET to turn off all other selected elements.
              H5S_SELECT_HYPERSLAB, dataspace_id, start, count, /RESET

              ; Create a simple dataspace to hold the result. Would have problems if
              ; not
              memory_space_id = H5S_CREATE_SIMPLE(count)

              ;--------------
              ;Read the data
              ;--------------
              data_U   = H5D_READ(dataset_id, FILE_SPACE=dataspace_id, $
                MEMORY_SPACE=memory_space_id)
              ;data(f,t)
              ;data   = H5D_READ(dataset_id)
              size_data = size(data_U)
              print,'Size of Data Now:', size_data[1], size_data[2]
              print,'-------------------------'
              print,'----Done Reading File----'
              print,'-------------------------'

              ;-------------------
              ;Transpose the data
              ;-------------------
              data_U = transpose(temporary(data_U))
              ;data_new(t,f); need this format for SPDYNPS

              H5S_CLOSE, dataspace_id
              H5D_CLOSE, dataset_id
              h5g_close, group_id
              H5F_CLOSE, file_id
                
    ;********************************
    ; Combine Data together 
    ;********************************   
    Data_Lin = sqrt(data_Q^2.0d + data_U^2.0d)
         
   ;****************************************************************
   ;                        Normalize the Data
   ;****************************************************************       
             
   print,'-------------------------'
   print, 'Start Normalize Data'
   print,'-------------------------'
              ;-------------------
              ; Normalize the Data
              ;-------------------
              nf=n_elements(data_lin(0,*))
              nt=n_elements(data_lin(*,0))

  ;normalize
  data_new = data_lin/data_I          
              
              print, ''
              print,'-------------------------'
              print, 'Normalize Data Done'
              print,'-------------------------'
   
              ; ------------------------------------------------------------------------
              ;  Preparing screen output
              ; ------------------------------------------------------------------------
              nf=n_elements(data_new(0,*))   
              nt=n_elements(data_new(*,0))
              
              nv=nt/1000.       ; for visualization on screen
              if nv ne long(nv) then nv=long(nv)+1 else nv=long(nv)
              ntv=long(nt/nv)
              mv=nf/500.
              if mv ne long(mv) then mv=long(mv)+1 else mv=long(mv)
              mtv=long(nf/mv)
              
              nnv=nt/2000.        ; for visualization in ps file
              if nnv ne long(nnv) then nnv=long(nnv)+1 else nnv=long(nnv)
              nntv=long(nt/nnv)
              ;mmv=nf/1000.
              mmv = nf/2000.
              if mmv ne long(mmv) then mmv=long(mmv)+1 else mmv=long(mmv)
              mmtv=long(nf/mmv)
              
              nwin=0 & xp=100 & yp=100
              
              
              ;---------------------------
              ; Make Plot
              ;--------------------------
              print,'----Making Plot----'
             ;set_plot,'PS'
                !p.multi=[0,1,2,1]
               SPDYNPS,rebin(data_lin(0:nntv*nnv-1,0:mmtv*mmv-1),nntv,mmtv),min(time_user),max(time_user), $
                      min(center_freq_user),max(center_freq_user), $
                      '','Frequency (MHz)','Raw: '+date+' Beam:'+j_string+' Linear Polarization', $
                      0,0,0,.05,.95,0,'.'
            
                    SPDYNPS,rebin(data_I(0:nntv*nnv-1,0:mmtv*mmv-1),nntv,mmtv),min(time_user),max(time_user), $
                      min(center_freq_user),max(center_freq_user), $
                      'Time [s]','Frequency (MHz)','Raw: '+date+' Beam:'+j_string+' I Polarization', $
                      0,0,0,.05,.95,0,'.'
            
                    xavg1 = total(data_lin,1)
                    yavg1 = total(data_lin,2)
                        xavg2= total(data_I,1) 

                    !p.multi=[0,1,2]
                    plot,CONGRID(center_freq_user,1000),CONGRID(xavg1/mean(xavg1),1000),/xsty,/ynoz,xtit='Frequency',ytit='Intensity',tit='Integrated spectrum of Lin Polarization: Before RFI MITIGATE';,$
                    yrange=[0.985,1.0]
                    plot,CONGRID(center_freq_user,1000),CONGRID(xavg2/mean(xavg2),1000),/xsty,/ynoz,xtit='Frequency',ytit='Intensity',tit='Integrated spectrum of I Polarization: Before RFI MITIGATE';,$
                    yrange=[0.985,1.0]
                   ; plot,CONGRID(time_user,5000),CONGRID(yavg1/mean(yavg1),5000),/xsty,/ynoz,xtit='Time',ytit='Intensity',tit='Integrated Time Series: After RFI MITIGATE'   
            
                SPDYNPS,rebin(data_new(0:nntv*nnv-1,0:mmtv*mmv-1),nntv,mmtv),min(time_user),max(time_user), $
                      min(center_freq_user),max(center_freq_user), $
                      'Time [s]','Frequency (MHz)','Normalized: '+date+' Beam:'+j_string+'  Linear Polarization', $
                      0,0,0,.05,.95,0,'.'
            ;!P.Multi = 0
             ; exit_ps 
               
        
;             SPDYNPS,rebin(data_new_slope(0:nntv*nnv-1,0:mmtv*mmv-1),nntv,mmtv),min(time_user),max(time_user), $
;               min(center_freq_user),max(center_freq_user), $
;               'Time [s]','Frequency (MHz)','Gain/Background/slope: '+date+' Beam:'+j_string+' P:'+k_string, $
;               0,1,0,.05,.95,0,'.'
               
             xavg_before = total(data_new,1)
            yavg_before = total(data_new,2)
            
              !p.multi=[0,1,2]
             plot,CONGRID(center_freq_user,1000),CONGRID(xavg_before/mean(xavg_before),1000),/xsty,/ynoz,xtit='Frequency',ytit='Intensity',tit='Integrated spectrum: Before RFI MITIGATE';,$
               yrange=[0.985,1.0]
              plot,CONGRID(time_user,5000),CONGRID(yavg_before/mean(yavg_before),5000),/xsty,/ynoz,xtit='Time',ytit='Intensity',tit='Integrated Time Series: After RFI MITIGATE'
              ;,$
             ;  yrange =[0.95,1.1]
      
  endif
  
  
 

  ;Linear Polarization need to read both Q and U and combine
  If S_run eq 2 then begin
    ;***************************************
    ;         Read V
    ;***************************************
    ;we always need to read in I data for normalization
    file_V = file_root+date[i]+'/'+date[i]+'_SAP000_B00'+j_string+$
      '_S3'+'_P000_bf.h5'
    print, file


    file_id     = h5f_open(file_V)
    group_id    = h5g_open(file_id,'SUB_ARRAY_POINTING_000/BEAM_00'+j_string)
    dataset_id  = h5d_open(group_id,'STOKES_3')   ;!!!! Stokes 3

    ; Open up the dataspace
    dataspace_id = H5D_GET_SPACE(dataset_id)
    dimensions   = H5S_GET_SIMPLE_EXTENT_DIMS(dataspace_id)
    print,'Orignal Size of Data:', dimensions[0],dimensions[1]

    ;-----------------------------------------------
    ;          Now choose our area of interest
    ; What elements of the data array do you want?
    ;-----------------------------------------------
    start = [min(w),min(u)]      ;w and u are the array inputs
    count = [num_bands,num_time]

    ;Select hyperslab
    ; Be sure to use /RESET to turn off all other selected elements.
    H5S_SELECT_HYPERSLAB, dataspace_id, start, count, /RESET

    ; Create a simple dataspace to hold the result. Would have problems if
    ; not
    memory_space_id = H5S_CREATE_SIMPLE(count)

    ;--------------
    ;Read the data
    ;--------------
    data_V   = H5D_READ(dataset_id, FILE_SPACE=dataspace_id, $
      MEMORY_SPACE=memory_space_id)
    ;data(f,t)
    ;data   = H5D_READ(dataset_id)
    size_data = size(data_V)
    print,'Size of Data Now:', size_data[1], size_data[2]
    print,'-------------------------'
    print,'----Done Reading File----'
    print,'-------------------------'

    ;-------------------
    ;Transpose the data
    ;-------------------
    data_V = transpose(temporary(data_V))
    ;data_new(t,f); need this format for SPDYNPS

    H5S_CLOSE, dataspace_id
    H5D_CLOSE, dataset_id
    h5g_close, group_id
    H5F_CLOSE, file_id
  
    ;****************************************************************
    ;                        Normalize the Data
    ;****************************************************************

    print,'-------------------------'
    print, 'Start Normalize Data'
    print,'-------------------------'
    ;-------------------
    ; Normalize the Data
    ;-------------------
   ; nf=n_elements(data_lin(0,*))
  ;  nt=n_elements(data_lin(*,0))

    ;normalize
    data_new = data_V/data_I
    
    
    print, ''
    print,'-------------------------'
    print, 'Normalize Data Done'
    print,'-------------------------'

    ; ------------------------------------------------------------------------
    ;  Preparing screen output
    ; ------------------------------------------------------------------------
    nf=n_elements(data_new(0,*))
    nt=n_elements(data_new(*,0))

    nv=nt/1000.       ; for visualization on screen
    if nv ne long(nv) then nv=long(nv)+1 else nv=long(nv)
    ntv=long(nt/nv)
    mv=nf/500.
    if mv ne long(mv) then mv=long(mv)+1 else mv=long(mv)
    mtv=long(nf/mv)

    nnv=nt/2000.        ; for visualization in ps file
    if nnv ne long(nnv) then nnv=long(nnv)+1 else nnv=long(nnv)
    nntv=long(nt/nnv)
    ;mmv=nf/1000.
    mmv = nf/2000.
    if mmv ne long(mmv) then mmv=long(mmv)+1 else mmv=long(mmv)
    mmtv=long(nf/mmv)

    nwin=0 & xp=100 & yp=100


    ;---------------------------
    ; Make Plot
    ;--------------------------
    print,'----Making Plot----'
    ;set_plot,'PS'
    !p.multi=[0,1,2,1]
;    SPDYNPS,rebin(data_lin(0:nntv*nnv-1,0:mmtv*mmv-1),nntv,mmtv),min(time_user),max(time_user), $
;      min(center_freq_user),max(center_freq_user), $
;      '','Frequency (MHz)','Raw: '+date+' Beam:'+j_string+' Linear Polarization', $
;      0,0,0,.05,.95,0,'.'

    SPDYNPS,rebin(data_V(0:nntv*nnv-1,0:mmtv*mmv-1),nntv,mmtv),min(time_user),max(time_user), $
      min(center_freq_user),max(center_freq_user), $
      'Time [s]','Frequency (MHz)','Raw: '+date+' Beam:'+j_string+' V Polarization', $
      0,0,0,.05,.95,0,'.'

    xavg1 = total(data_lin,1)
    yavg1 = total(data_lin,2)
    xavg2= total(data_I,1)

    !p.multi=[0,1,2]
    plot,CONGRID(center_freq_user,1000),CONGRID(xavg1/mean(xavg1),1000),/xsty,/ynoz,xtit='Frequency',ytit='Intensity',tit='Integrated spectrum of V Polarization: Before RFI MITIGATE';,$
    yrange=[0.985,1.0]
    plot,CONGRID(center_freq_user,1000),CONGRID(xavg2/mean(xavg2),1000),/xsty,/ynoz,xtit='Frequency',ytit='Intensity',tit='Integrated spectrum of I Polarization: Before RFI MITIGATE';,$
    yrange=[0.985,1.0]
    ; plot,CONGRID(time_user,5000),CONGRID(yavg1/mean(yavg1),5000),/xsty,/ynoz,xtit='Time',ytit='Intensity',tit='Integrated Time Series: After RFI MITIGATE'

    SPDYNPS,rebin(data_new(0:nntv*nnv-1,0:mmtv*mmv-1),nntv,mmtv),min(time_user),max(time_user), $
      min(center_freq_user),max(center_freq_user), $
      'Time [s]','Frequency (MHz)','Normalized: '+date+' Beam:'+j_string+'  V Polarization', $
      0,0,0,.05,.95,0,'.'
    ;!P.Multi = 0
    ; exit_ps


    ;             SPDYNPS,rebin(data_new_slope(0:nntv*nnv-1,0:mmtv*mmv-1),nntv,mmtv),min(time_user),max(time_user), $
    ;               min(center_freq_user),max(center_freq_user), $
    ;               'Time [s]','Frequency (MHz)','Gain/Background/slope: '+date+' Beam:'+j_string+' P:'+k_string, $
    ;               0,1,0,.05,.95,0,'.'

    xavg_before = total(data_new,1)
    yavg_before = total(data_new,2)

    !p.multi=[0,1,2]
    plot,CONGRID(center_freq_user,1000),CONGRID(xavg_before/mean(xavg_before),1000),/xsty,/ynoz,xtit='Frequency',ytit='Intensity',tit='Integrated spectrum: Before RFI MITIGATE';,$
    yrange=[0.985,1.0]
    plot,CONGRID(time_user,5000),CONGRID(yavg_before/mean(yavg_before),5000),/xsty,/ynoz,xtit='Time',ytit='Intensity',tit='Integrated Time Series: After RFI MITIGATE'
    ;,$
    ;  yrange =[0.95,1.1]
  
   endif   
   
   start_time = start_time + N_time  ; update start_time
   print, 'New Start time:' ,start_time     
          endfor ;time loop
  ;endfor ; end for pol 
  device,/close


  
  
  endfor ; beam
  
endfor  ;end data loop
set_plot,'PS'
device,/close

print,'The End of Program: Goodbye'
print,'SYSTIME=',SYSTIME(1) - st, ' sec'
journal
end
