;*****************************************************
;************* Read LOFAR H5 DATA ********************
;*****************************************************
;
;Author: Jake Turner (UVA)
;        J.M. Griessmeier (CNRS) and Philippe Zarka (Obspm)
;        
;Inputs: filename (eg. filename='/databf/lofar/EXO/LC5_DDT_002/L429868/L429868_SAP000_B001_S0_P000_bf')
;        tmin, tmax or (if /NSPECTRA) in number of spectra
;        fmin, fmax in MHz
;        b; Beam
;        S; Pol  polarization IQUV (I =0, Q = 1, U = 2, V = 3)
;
;Outputs: x (data)
;         t (time)
;         f (freq) 
;         nt (# of spectra)
;         nf (# of channels)
;
;OPTIONAl Outputs:
;         nnt    ;# of spectra in the full dataset
;         nnf    ;# of freq chanels in the full dataset
;         dt     ;sampling_time
;         Time   (Time array of the entire observation)
;Flags: 
;      /SAV          ;save the requested data in a sav file
;      /VERB         ;print info about the file
;      /CREATE       ;create the symbolic link to the data
;      /REMOVE       ;remove symbolic link 
;      /NSPECTRA     ;Use number of spectra for tmin and tmax
;      /REMOVE_BADF  ;remove 0 and 64 element of each subband)
;
;HEADER INFO Flags:
;      /channels_per_subband   ; Channels per subband 
;      /NOF_SUBBANDS	       ; # of subbands
;      /SAMPLING_TIME          ; SAMPLING time 
;      /ANTENNA_SET	       ; what was the ANTENNA setup
;      /NOF_SAMPLES	       ; # of time samples 
;      /START_MJD	       ; Starting time in MJD
;      /MJD_sampling_time      ; Sampling time in MJD (days)
;      /NOF_stations           ; # of stations used
;      /channel_width          ; Width of the freq channel
;      /min_freq               ; Min Freq of the observation
;      /max_freq               ; Max Freq of the observation
;      /READ_HEADER            ; ONLY Read the header of the observations
;      /tMJD                   ; The full MJD time of the entire observation
;      /UT                     ; UT time of the selected time range 
;      /AllUT                  ; UT time of the entire observation
;
pro READ_H5_DATA, filename,tmin,tmax,fmin,fmax,b,S,x,t,f,nt,nf,$
        nnt=nnt,nnf=nnf,dt=dt,Time=Time,Freq=center_freq,$
	SAV=SAV,VERB=VERB,CREATE=CREATE,REMOVE=REMOVE,NSPECTRA=NSPECTRA,$
	REMOVE_BADF=REMOVE_BADF,$
	channels_per_subband=channels_per_subband,NOF_SUBBANDS=NOF_SUBBANDS,SAMPLING_TIME=SAMPLING_TIME,$
	ANTENNA_SET=ANTENNA_SET,NOF_SAMPLES=NOF_SAMPLES,START_MJD=START_MJD,MJD_sampling_time=MJD_sampling_time,$
	NOF_stations=NOF_stations,channel_width=channel_width,min_freq=min_freq,max_freq = max_freq,$
	READ_HEADER=READ_HEADER,tMJD=MJD,AMJD=MJD_All,UT=UT,AllUT=UT_all
; e.g. filename='/databf/lofar/EXO/LC5_DDT_002/L429868/L429868_SAP000_B001_S0_P000_bf'
; tmin,tmax in sec or (if /NSPECTRA) in number of spectra
; fmin,fmax in MHz
; b ; beam number (0-4); as a string
; /CREATE, /REMOVE = symbolic file link

if keyword_set(CREATE) then begin
  spawn,'ln -s '+filename+'.h5'
  spawn,'ln -s '+filename+'.raw'
endif

;file=DELPATH(filename)
file		= filename
file_id		= h5f_open(file+'.h5')                                                  
sub_file_id     = h5g_open(file_id,'SUB_ARRAY_POINTING_000')                    
group_id        = h5g_open(file_id,'SUB_ARRAY_POINTING_000/BEAM_00'+strtrim(string(uint(B)),1))             
dataset_id      = h5d_open(group_id,'STOKES_'+strtrim(string(uint(S)),1))                                   
dataspace_id    = H5d_get_space(dataset_id)
dimensions      = H5s_get_simple_extent_dims(dataspace_id)
bandwidth       = h5a_read(h5a_open_name(file_id,'BANDWIDTH'))
bandwidth_unit  = h5a_read(h5a_open_name(file_id,'BANDWIDTH_UNIT'))			; MHz
min_freq 	= h5a_read(h5a_open_name(file_id,'OBSERVATION_FREQUENCY_MIN'))
max_freq	= h5a_read(h5a_open_name(file_id,'OBSERVATION_FREQUENCY_MAX'))
freq_unit	= h5a_read(h5a_open_name(file_id,'OBSERVATION_FREQUENCY_UNIT'))	; MHz
NOF_stations    = h5a_read( h5a_open_name(file_id,'OBSERVATION_NOF_STATIONS' ))
ANTENNA_SET     = h5a_read( h5a_open_name(file_id,'ANTENNA_SET' ))
channel_width   = h5a_read(h5a_open_name(group_id,'CHANNEL_WIDTH'))	
channel_width_unit= h5a_read(h5a_open_name(group_id,'CHANNEL_WIDTH_UNIT'))		; Hz
sampling_time   = h5a_read(h5a_open_name(group_id,'SAMPLING_TIME'))         ;sec
START_MJD       = h5a_read(h5a_open_name(sub_file_id,'EXPTIME_START_MJD'))
end_mjd	        = h5a_read(h5a_open_name(sub_file_id,'EXPTIME_END_MJD'))
start_UT	= h5a_read(h5a_open_name(sub_file_id,'EXPTIME_START_UTC')) ;(e.g. 2016-09-08T00:00:00.000000000Z)
tot_int_time	= h5a_read(h5a_open_name(sub_file_id,'TOTAL_INTEGRATION_TIME'))	; sec
NOF_SUBBANDS    = h5a_read(h5a_open_name(dataset_id,'NOF_SUBBANDS') )
channels_per_subband = h5a_read(h5a_open_name(group_id,'CHANNELS_PER_SUBBAND'))
NOF_SAMPLES     =  h5a_read(h5a_open_name(dataset_id,'NOF_SAMPLES') )
nof_subbands	= h5a_read(h5a_open_name(dataset_id,'NOF_SUBBANDS'))
group_id_1	= h5g_open(file_id,'SUB_ARRAY_POINTING_000/BEAM_00'+strtrim(string(uint(B)),1)+'/COORDINATES/COORDINATE_1')
center_freq     = h5a_read(h5a_open_name(group_id_1,'AXIS_VALUES_WORLD'))		; Hz
center_freq	= center_freq/1.0d6							; MHz
nf		= n_elements(center_freq)
samples		= dindgen(NOF_SAMPLES)
time		= samples*SAMPLING_TIME   ;entire time of observation
;---------
;MJD time 
;---------
MJD_sampling_time   = SAMPLING_TIME/86400.d
MJD_all                 = dindgen(NOF_SAMPLES,start=START_MJD,INCREMENT=MJD_sampling_time)
;--------
;UT time
JD_all = MJD_all+2400000.5d
caldat,JD_all,Month,Day,Year,Hour,Min,Sec
UT_all = Hour + Min/60.d  + Sec/(3600.0d)   ;UT time in hours (3.5 hours = 3 hours 30 minutes)

If not(keyword_set(READ_HEADER)) then begin 
  if keyword_set(NSPECTRA) then u = where(samples ge tmin and samples lt tmax) else u = where(time ge tmin and time lt tmax)
  nt			= n_elements(u) 
  t			= time(u)
  mjd 			= mjd_all(u)
  UT 			= UT_all(u) 

  If keyword_set(REMOVE_BADF) then begin
    v = where(center_freq ge min_freq and center_freq le max_freq)  ;If removing, first read full frequency range
  Endif
  If not(keyword_set(REMOVE_BADF)) then v = where(center_freq ge fmin and center_freq le fmax)
  nf			= n_elements(v)
  f			  = center_freq(v)
  
  start	  = [min(v),min(u)]
  count    = [nf,nt]
  
  h5s_select_hyperslab,dataspace_id,start,count,/RESET
  memory_space_id		= h5s_create_simple(count)
  x			            = h5d_read(dataset_id,FILE_SPACE=dataspace_id,MEMORY_SPACE=memory_space_id)
  x			            = transpose(temporary(x))
  
  H5S_CLOSE, memory_space_id
  if keyword_set(REMOVE) then begin
    spawn,'rm ./'+file+'.h5'
    spawn,'rm ./'+file+'.raw'
  endif
  
  If keyword_set(REMOVE_BADF) then begin
    x_bad = 1
    x=reform(temporary(x),nt,channels_per_subband,NOF_SUBBANDS )
    ;x=x(*,1:channels_per_subband-2,*)
    x =x(*,x_bad:channels_per_subband-(x_bad+1),*)
    ;x=reform(temporary(x),nt,(channels_per_subband-2)*NOF_SUBBANDS )
    x=reform(temporary(x),nt,(channels_per_subband-(x_bad*2))*NOF_SUBBANDS )
    f=reform(temporary(f),channels_per_subband,NOF_SUBBANDS )
    ;f=f(1:channels_per_subband-2,*)
    f=f(x_bad:channels_per_subband-(x_bad+1),*)
    ;f=reform(temporary(f),(channels_per_subband-2)*NOF_SUBBANDS )
    f=reform(temporary(f),(channels_per_subband-(x_bad*2))*NOF_SUBBANDS)
    nf=n_elements(f)
    center_freq = f

    v = where(f ge fmin and f le fmax)
    nf      = n_elements(v)
    f       = f(v)
    x       = x[*,v]
  endif
endif  ;end read data
 h5s_close, dataspace_id
 h5d_close, dataset_id
 h5g_close, group_id
 h5g_close, sub_file_id
 h5f_close, file_id

  if keyword_set(VERB) then begin
    print,filename
    print,'Original size of data:',dimensions[0],dimensions[1]
    print,'Bandwidth:',bandwidth,' ',bandwidth_unit
    print,'Min,max Observation frequency:',min_freq,max_freq,' ',freq_unit
    print,'Channel width:',channel_width,' ',channel_width_unit
    print,'Start,stop MJD:',start_mjd,end_mjd
    print,'Min,max Time:',min(time),max(time),' sec'
    print,'Total integration time, Sampling time:',tot_int_time,sampling_time,' sec'
    print,'#samples, #subbands, #frequencies:',nof_samples,nof_subbands,nf
    
   If not(keyword_set(READ_HEADER)) then begin 
    print,'Min,max User center frequency:',min(f),max(f),' MHz'
    print,'Min,max User Time:',min(t),max(t),' sec'
    If keyword_set(REMOVE_BADF) then print,'Number of channels at subband edges that are bad,', x_bad 
   endif  
  endif

if keyword_set(VERB) then help,x,f,nf,t,nt
if keyword_set(SAV)  then save,x,f,nf,t,nt,filename=file+'.sav'

nnf=dimensions[0] & nnt=dimensions[1] & dt=sampling_time

return
end
