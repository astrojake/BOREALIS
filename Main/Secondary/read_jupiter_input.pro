;NAME: 
;    read_input
;    
;INPUT: filename with inputs (Input.dat) ;example file is supplied with the code 
;
;OUTPUT: 
;
; MODIFICATION HISTORY: v1: creation
;                       v2: Update for NenuFAR (Jun 13, 2023), works with input file v2 

pro read_jupiter_input,input_file,LOFAR_RUN = LOFAR_RUN,Nancay_RUN = Nancay_RUN, NenuFAR_RUN=NenuFAR,$
                        datafile_root=datafile_root,filename=filename,RFI_filename=RFI_filename,$
                        target_beam=target_beam,$
                        tmin_Jup=tmin_Jup,tmax_Jup=tmax_Jup,fmin_Jup=fmin_Jup,fmax_Jup=fmax_Jup,$
                        backfile_root=backfile_root,backfilename=backfilename,$
                        fmin_add=fmin_add,scale=scale,method2_Jup=method2_Jup,VERBOSE=VERBOSE
                        
                              
;**********************
;Declare Variables (create loop)
;**********************
;line = strarr(290,1)
line = ''

;*************************************
;Open File
;;*************************************
openr,lun,input_file,/get_lun

;-----------------------
;initalize variables 
;-----------------------
LOFAR_RUN = 0     & Nancay_RUN = 0      & NenuFAR_RUN = 0 
datafile_root=''  & filename=''         & target_beam=0    
tmin_Jup=0.0      & tmax_Jup=0.0        & fmin_Jup=0.0
fmax_Jup=0.0      & back_root=''        & RFI_filename=''
backfilename=''   & backfile_root = ''  & fmin_add=0.0      
& scale=0.0       method2_Jup = 1
;**********************
;Read Variables
;**********************
Readf, lun, line ;**************** ; line1
Readf, lun, line ;**************** ; line2
Readf, lun, line ;**************** ; line3
Readf, lun, line ;**************** ; line4
Readf, lun, line ;**************** ; line5
Readf, lun, line ;**************** ; line6
Readf, lun, line ;**************** ; line7
Readf, lun, line ;**************** ; line8
Readf, lun, LOFAR_RUN,format='(T57,I)'; line9
Readf, lun, NANCAY_RUN,format='(T57,I)' ;**************** ; line10
Readf, lun, line ;**************** ; line
Readf, lun, line ;**************** ; line
Readf, lun, line ;**************** ; line
Readf, lun, line ;**************** ; line
Readf, lun, datafile_root,format='(T57,A)'
Readf, lun, filename,format='(T57,A)'
readf, lun, target_beam,format='(T57,I)'
readf, lun, RFI_filename,format='(T57,A)'
readf, lun, tmin_Jup,format='(T57,F)'
readf,lun, tmax_Jup, format='(T57,F)'
readf, lun, fmin_Jup,format='(T57,F)'
readf,lun, fmax_Jup, format='(T57,F)'
Readf, lun, line ;**************** ; line
Readf, lun, line ;**************** ; line
Readf, lun, line ;**************** ; line
Readf, lun, line ;**************** ; line
Readf, lun, backfile_root,format='(T57,A)'
Readf, lun, backfilename,format='(T57,A)'
Readf, lun, line ;**************** ; line
Readf, lun, line ;**************** ; line
Readf, lun, line ;**************** ; line
Readf, lun, line ;**************** ; line
Readf, lun, fmin_add, format='(T57,F)'
Readf, lun, scale,format='(T57,F)'
Readf, lun, line ;**************** ; line
Readf, lun, line ;**************** ; line
Readf, lun, line ;**************** ; line
Readf, lun, line ;**************** ; line
readf, lun, method2_Jup,format='(T57,I)'
Free_Lun, lun

If keyword_set(VERBOSE) then begin
  print, '********************************************************'
  print, '********** Jupiter Input Parameters ********************'
  print, '********************************************************'
  print, 'Input Jupiter file name ',input_file
  print, 'Run LOFAR Jupiter Observations (0=NO, 1=YES): ',LOFAR_RUN
  print, 'Run Nancay Jupiter Observations (0=NO, 1=YES): ',Nancay_RUN
  print, 'Root file of the data :',datafile_root
  print, 'Name of H5 file to run:',filename
  print, 'Target Beam',target_beam
  print, 'RFI File for Target',RFI_filename
  print, 'start time of Jupiter data',tmin_Jup
  print, 'end time of Jupiter data',tmax_Jup
  print,'Min Frequency (MHZ)' ,fmin_Jup
  print,'Max Frequency (MHZ) ',fmax_Jup
  print,'root file of background',backfile_root
  print,'Name of file:' , backfilename
  print,'Min frequency in LOFAR data to add Jupiter data',fmin_add
  print,'Reduction Factor', scale
  print,'Default Method 2 Processing (1=Default),',method2_Jup
  print, '********************************************************'
endif   
end
