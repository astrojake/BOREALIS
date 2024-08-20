;NAME: 
;    read_input
;    
;INPUT: filename with inputs (Input.dat) ;example file is supplied with the code 
;
;OUTPUT: 
;
; MODIFICATION HISTORY: Version 3 (June 14, 2017; Added Process_METHOD2) 
;                       Version 4 (March 22, 2019; Added NENUFAR and DYNSPECMS, updated all and removed unnecessary or redunant keywords) 
;                       Version 5 (March 22, 2021; Updated various inputs for NENUFAR and DynspecMS)
pro read_inputv5,input_file,planet=planet,period=period,TO=TO,date=date,fileroot=file_root,$
               RQUICKLOOK=run_QUICKLOOK,RRFI=Run_RFI,RProcessing=Run_processing,RFFT=run_FFT,RDISP=run_DISP,RPostProcessing=run_postprocessing,$
               PLOT_RAW=PLOT_RAW,DATA_NORM=DATA_NORM,DATA_SLOPE=DATA_SLOPE,fminql=fmin_ql,fmaxql=fmax_ql,qltmin=tmin_ql,N_timeql=N_time_ql,$
               qlsteps=steps_ql,bminql=bmin_ql,bmaxql=bmax_ql,Sminql=Smin_ql,Smaxql=Smax_ql,$
               start_beam=start_beam,end_beam=end_beam,$
               Polarization=S,tmin=tmin,nspectra=nspectra,steps=steps,$
               PS=PS,VERBOSE=VERBOSE,mask_bit=mask_bit,bit=bit,send_email=send_email,email_address=email_address,$
               sigma_patrol=patrol_value,Si_le_sig=le_sig_value,patrol=patrol,LE_SIG=LE_SIG,SUM=SUM,PEX=PEX,$
               TF_correction=TimeFreq_correction,Apply_Corrections=Apply_Corrections,$
               Drebin_data=rebin_data,Frebin_freq=rebin_freq,Trebin_time=rebin_time,$
               beam_DISP=beam_DISP,pulsarname=pulsarname,dm=dm,size_disp=size_disp,$
               beam_FFT=beam_FFT,PFFT=Period_FFT,pmin=pmin,fmaxp=fmaxp,fminp=fminp,$
               Target_Beam = Target_Beam, Sky_Beam=Sky_Beam,Q1=Q1,Q2=Q2,Q3=Q3,Q4=Q4,$
               REMOVE_BADF=REMOVE_BADF,$
               SELECT_Freq,fminimum=fmin,fmaximum=fmax,$
               quantile_norm,quantile=quantile,mask_save=mask_save,$
               NFREQ_PATROL=NFREQ_PATROL,NTIME_PATROL=NTIME_PATROL,NFREQ_LESIG=NFREQ_LESIG,NTIME_LESIG=NTIME_LESIG,$
               MF=MF,MT=MT,THRESHOLD=THRESHOLD,FPMIN=PMIN_F,ExpansionF1=EXPF,TPMIN=PMIN_T,ExpansionT1=EXPT,F2PMIN=PMIN_F2,EXPF2=EXPF2,T2PMIN=PMIN_T2,EXPT2=EXPT2,$
               RFI=RFI,full_clean_save=full_clean_save,save_datagain=save_datagain,$
               RUN_DISP_TIME=RUN_DISP_TIME, rebin_time_disp=rebin_time_disp,$
               ddm=ddm,$
               create_combine_mask=create_combine_mask,combine_mask=combine_mask,NUM_MASKS=NUM_MASKS,COMBO_BEAMS=COMBO_BEAMS,$
               Smooth_Surface=Smooth_Surface,N_smooth=N_smooth,Process_METHOD2=Process_METHOD2,$
               skip_rebinloop=skip_rebinloop,save_RFIMaskData=save_RFIMaskData,$
               TI=TI,Default_MASK=Default_MASK,Run_Mask=Run_Mask,p_threshold=threshold_mask,Default_FreqRange=Default_FreqRange,Run_FreqRange=Run_FreqRange,$
               freq_minextra1=freq_minextra1,freq_maxextra2=freq_maxextra2,Default_Time=Default_Time,QTmin=QTmin,QTmax=QTmax,Default_Threshold=Default_Threshold,Exthreshold_extra=threshold_extra,$
               max_threshold=max_threshold,pthreshold=threshold_p,tau_Qplot=tau_Qplot,slope=slope,Qwin=Qwin,ElC=ElCorrect,save_gauss=save_gauss,gauss_steps=gauss_steps,file_g=file_g,$
               root_gauss =root_gauss,save_Q=save_Q,OnBeam_Label=OnBeam_Label,OffBeam_Label=OffBeam_Label,ExB1=Extra_Beam1,ExB2=Extra_Beam2,$
               Run_Rebintime2=Run_Rebintime2,rebin_time2=rebin_time2,Run_RebinFreq2=Run_RebinFreq2,rebin_freq2=rebin_freq2,$
               file_root_combine=file_root_combine,filenname_combinedates=filenname_combinedates,$
               create_dates_combine=create_dates_combine,Combine_NBeams=Combine_NBeams,dates_combine,Combine_Beams=Combine_Beams,$
               RUN_JUPITER=RUN_JUPITER,input_file_J=input_file_J,$
               PREMOVE=REMOVE,PAPERPLOTS=PAPERPLOTS,MINT=MINT,MAXT=MAXT,$
               LOFAR=LOFAR,NENUFAR=NENUFAR,DYNSPECMS=DYNSPECMS,filename=filename,$
               COMPLEX_NAME=COMPLEX_NAME,PCA=PCA,Full_RFI=Full_RFI,use_combinemask=use_combinemask
             
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
planet            = ' '     & period            = 0.0d 
TO                = 0.0d    & date              =''
file_root         =''       & email_address     = ''
start_beam        = 0       & end_beam          = 0 
S                 = 0       & tmin              = 0.0d 
nspectra          = 0.0d      &steps              = 0.0d
run_QUICKLOOK     = 0       &run_RFI            = 0 
run_processing    = 0       &run_FFT            = 0
run_DISP          = 0       &run_postprocessing = 0 
patrol_value      = 0.0d    &le_sig_value       = 0.0d
Patrol            = 0       &LE_SIG            = 0
SUM               = 0      & PEX               = 0
TimeFreq_correction = 0     &Apply_Corrections  = 0 
rebin_freq        =0.       &rebin_time         = 0.
rebin_freq2       =0.
beam_DISP         = 0.       &pulsarname         =''
dm                = 0.0d    &size_disp          = 0.0d  
Period_FFT        = 0.0d    &pmin 		          = 0.0d 
fmaxp             = 0.0d    &fminp              = 0.0d 
PLOT_RAW          = 1.       &DATA_NORM          = 1. 
DATA_SLOPE        = 1.       &fmin_ql            = 0.0d
 fmax_ql          = 0.0d    &tmin_ql            = 0.0d
 N_time_ql        = 0.0d    &steps_ql           = 1.
 bmin_ql          = 0       & bmax_ql           = 0 
Smin_ql           = 0       &Smax_ql            = 0 
REMOVE_BADF       = 1.       &fmin               = 0.0d
fmax              = 0.0d    &quantile           = 0.0d
mask_save         = 0.0     &MF                 = ''
Mt                = ''      & Process_METHOD2 = 1
rebin_time_disp   =0.       & Target_Beam = 1
beam1             = 0       &Sky_Beam = 1
beam2             = 0       & gauss_steps = 10000.
file_root_combine = ''      & file_g = ''
filenname_combinedates =''  & root_gauss = ''
create_dates_combine = 0.   & OnBeam_label = 'On Beam'
Combine_NBeams = 0.         & OffBeam_label = 'Off Beam'
rebin_time_disp =0.         &freq_minextra1 = 1.0 & freq_maxextra2 = 1.0
mask1 = 0 & mask2 = 0 & mask3 = 0 & mask4 = 0
RFI  = 0  & LOFAR = 0 & NENUFAR = 0 & DYNSPECMS = 0 
input_file_J = ''  & filename = '' &REMOVE = 0 & PAPERPLOTS = 0 & MINT  = 0.0 & MAXT = 0.0 
COMPLEX_NAME = 0   & PCA= 0  & MT_file ='' & Mf_file = ''  &root_code =''

;**********************
;Read Variables
;**********************
i =0 ;reset (counting lines)
a = ''
Readf, lun, line ;****************************************************************
Readf, lun, line ;****************************************************************
Readf, lun, line ;Input File for BOREALISv4 (LOFAR, NenuFAR, or DynSpecMS)
Readf, lun, line ;******************(Version 6 of Input file): March 26, 2021)****
Readf, lun, line ;****************************************************************
Readf, lun, line ;-------------------------------------------------------------
Readf, lun, line ;     Target Info
Readf, lun, line ;-------------------------------------------------------------
Readf, lun, planet,format='(T57,A)'         ;Target Name (e.g. TauBootis)
Readf, lun, period,format='(T57,D)'         ;Target Period (days)
Readf, lun, TO,    format='(T57,D)'         ;Reference Time (Mid-Transit/periastron) (MJD)
Readf, lun, date,  format='(T57,A)'         ;Date (e.g. LOFAR format: L429868, NENUFAR: UT date):
Readf, lun, file_root,format='(T57,A)'      ;Root Directory of Data (e.g. /databf/exoplanet/LC5/):
Readf, lun, line ;                          ;
Readf, lun, line ;-------------------------------------------------------------
Readf, lun, line ;   Telescope Information (which data to run?)
Readf, lun, line ;-------------------------------------------------------------
Readf, lun, LOFAR,format='(T57,A)'          ;RUN LOFAR (0=NO, 1=YES)
Readf, lun, NENUFAR,format='(T57,A)'        ;RUN NENUFAR (0=NO, 1=YES)
Readf, lun, DYNSPECMS,format='(T57,A)'      ;RUN DYNSPECMS (0=NO, 1=YES)
Readf, lun, line ;
Readf, lun, line ;-------------------------------------------------------------
Readf, lun, line ;File to Run (required for NENUFAR and DYNSPECMS)
Readf, lun, line ;-------------------------------------------------------------
Readf, lun, filename,format='(T57,A)'       ;Filename to run (see README for format) 
Readf, lun, line ;
Readf, lun, line ;-------------------------------------------------------------
Readf, lun, line ;                 What do you want to Run?
Readf, lun, line ;-------------------------------------------------------------
Readf, lun, run_QUICKLOOK,format='(T57,F)'        ;Run QUICKLOOK (0=NO, 1=YES): 
Readf, lun, run_RFI,format='(T57,F)'              ;Run RFI (0=NO, 1=YES):
Readf, lun, run_processing,format='(T57,F)'       ;RUN Processing (0=NO, 1=YES):
Readf, lun, run_FFT,format='(T57,F)'              ;Run FFT (0=NO, 1=YES):
Readf, lun, run_DISP,format='(T57,F)'             ;RUN DE-DISP (0=NO, 1=YES):
Readf, lun, run_postprocessing,format='(T57,F)'   ;Run Post-processing (0=NO, 1=YES):
Readf, lun, line ;  
Readf, lun, line ;------------------------------------------------------------------------------------
Readf, lun, line ;     QuickLook Settings
Readf, lun, line ;------------------------------------------------------------------------------------
Readf, lun, PLOT_RAW,format='(T57,F)'             ;Plot Raw Data (0=No, 1 = YES):
Readf, lun, DATA_NORM ,format='(T57,F)'           ;Plot Normalized Data (0=No, 1 = YES):
Readf, lun, DATA_SLOPE,format='(T57,F)'           ;Plot Data Normalized by a slope (0=No, 1 = YES):
Readf, lun,fmin_ql,format='(T57,F)'               ;Min Frequency (MHz):
Readf, lun, fmax_ql,format='(T57,F)'              ;Max Frequency (MHz):
Readf, lun, N_time_ql,format='(T57,F)'            ;Amount of time per step (secs):
Readf, lun, N_Quick_full,format='(T57,F)'         ;Run full file (0=NO, 1= YES)
Readf, lun, tmin_ql,format='(T57,F)'              ;Min Time to run (secs, Not used if running full file):
Readf, lun, steps_ql,format='(T57,F)'             ;Number of steps (Not used if running full file):
Readf, lun, bmin_ql,format='(T57,I)'              ;Starting Beam:
Readf, lun, bmax_ql,format='(T57,I)'              ;Ending Beam:
Readf, lun, Smin_ql,format='(T57,I)'              ;Starting Polarization (I=0,Q=1,U=2,|V|=3,V^2=4,V'=5):
Readf, lun, Smax_ql,format='(T57,I)'              ;Ending Polarization   (L=6): 
Readf, lun, line 
Readf, lun, line ;----------------------------------------------------------------------------------------------
Readf, lun, line ;All Settings Below are for the Main Code (RFI,processing,post-processing)
Readf, lun, line ;----------------------------------------------------------------------------------------------
Readf, lun, line 
Readf, lun, line ;----------------------------------------------------------------------------------------------
Readf, lun, line ;         Beams to run in loop
Readf, lun, line ;----------------------------------------------------------------------------------------------
Readf, lun, start_beam,format='(T57,I)';Starting Beam:
Readf, lun, end_beam,format='(T57,I)'  ;Ending Beam
Readf, lun, line 
Readf, lun, line ;----------------------------------------------------------------------------------------------
Readf, lun, line ;        Polarization (see README for more details)
Readf, lun, line ;----------------------------------------------------------------------------------------------
Readf, lun, S,format='(T57,I)'  ;Polarization (I=0,Q=1,U=2,V=3,V_prime=4,L=5):
Readf, lun, line ;
Readf, lun, line ;----------------------------------------------------------------
Readf, lun, line ;               Time setup
Readf, lun, line ;----------------------------------------------------------------
Readf, lun, tmin,format='(T57,F)'     ;Min Time to run (secs):
Readf, lun, nspectra,format='(T57,F)' ;Number of spectra (e.g. 4000):
Readf, lun, Full_RFI,format='(T57,F)' ;Run full file (0=NO, 1=YES)
Readf, lun, steps,format='(T57,F)'    ;Number of steps (Ignored if running full file):  
Readf, lun, line 
Readf, lun, line ;----------------------------------------------------------------
Readf, lun, line ;              General Settings
Readf, lun, line ;----------------------------------------------------------------
Readf, lun, PS,format='(T57,F)'              ;Create PS Files (0=NO, 1=YES):
Readf, lun, VERBOSE,format='(T57,F)'         ;Verbose (0=No, 1=YES):
Readf, lun, mask_bit,format='(T57,F)'        ;Read/Save Mask in Bits (0=NO, 1=YES)
Readf, lun, bit,format='(T57,F)'             ;Manipulate the mask in Bits (0=NO, 1=YES):
Readf, lun, send_email,format='(T57,F)'      ;Send an email about code progress (0=NO, 1=YES):
Readf, lun, email_address,format='(T57,A)'   ;Email Address:
Readf, lun, line 
Readf, lun, line ;----------------------------------------------------------------
Readf, lun, line ;        RFI Settings
Readf, lun, line ;----------------------------------------------------------------
Readf, lun, patrol_value,format='(T57,F)' ;Sigma Value for Patrol (Default: 5.5):
Readf, lun, le_sig_value,format='(T57,F)' ;Sigma Value for Le_SIF (Default: 5.5):
Readf, lun, patrol,format='(T57,F)'       ;Run Patrol (0=NO, 1=YES):
Readf, lun, le_sig,format='(T57,F)'       ;Run LESIG (0=NO, 1=YES):
Readf, lun, sum,format='(T57,F)'          ;Run SUM (0=NO, 1=YES):
Readf, lun, pex,format='(T57,F)'          ;Run PEX (0=NO, 1=YES):
Readf, lun, line ;
Readf, lun, line ;----------------------------------------------------------------
Readf, lun, line ;Processing Settings
Readf, lun, line ;----------------------------------------------------------------
Readf, lun, TimeFreq_correction,format='(T57,F)' ;Find the Time-Frequency Correction (0=NO, 1=YES):
Readf, lun, Apply_Corrections,format='(T57,F)'   ;Apply the Time-Frequency Correction (0=NO, 1=YES):
readf, lun, Process_METHOD2,format='(T57,I)'     ;Run Method 2 (0=NO, 1 = YES) (Default: YES):
Readf, lun, line ;
Readf, lun, line ;----------------------------------------------------------------
Readf, lun, line ;           Rebin Settings
Readf, lun, line ;----------------------------------------------------------------
Readf, lun, rebin_data,format='(T57,F)' ;Rebin the data (0=NO, 1=YES):
Readf, lun, rebin_freq,format='(T57,F)' ;Rebin Frequency (kHz) (Default: 45.0):
Readf, lun, rebin_time,format='(T57,F)' ;Rebin Time (secs) (Default: 1.0): 
Readf, lun, line ;
Readf, lun, line ;----------------------------------------------------------------
Readf, lun, line ; De-dispersion Settings
Readf, lun, line ;----------------------------------------------------------------
Readf, lun, beam_DISP,format='(T57,F)' ;What beam to run de-dispersion (any beam):
Readf, lun, pulsarname,format='(T57,A)';Name of de-disperson beam (e.g. B0823+26): 
Readf, lun, dm,format='(T57,F)'        ;dm (pc^3 cm):
Readf, lun, size_disp,format='(T57,F)' ;Rebin Frequency for De-dispersion (MHz):
Readf, lun, line
Readf, lun, line ;----------------------------------------------------------------
Readf, lun, line ;    FFT Settings 
Readf, lun, line ;----------------------------------------------------------------
Readf, lun, beam_FFT,format='(T57,F)'   ;What beam to run FFT (pulsar beam):
Readf, lun, Period_FFT,format='(T57,F)' ;Period (usually Pulsar) for FFT (secs):
Readf, lun, pmin,format='(T57,F)'       ;Min Mask Threshold for FFT:
Readf, lun, fminp,format='(T57,F)'      ;Min Frequency for FTT (MHz):
Readf, lun, fmaxp,format='(T57,F)'      ;Max Frequency for FFT (MHz):  
Readf, lun, line  
Readf, lun, line ;----------------------------------------------------------------
Readf, lun, line ; Post-Processing Settings
Readf, lun, line ; Note:Post-processing is very complicated. Make sure you read README
Readf, lun, line ;----------------------------------------------------------------
Readf, lun, Target_Beam,format='(T57,I)'  ;Target Beam to Run (Default: 0):
Readf, lun, Sky_Beam,format='(T57,I)'     ;Sky Beam to Run (Default: 2):
Readf, lun, Q1,format='(T57,F)'           ;Run Q1 (0=No, 1 = YES):
Readf, lun, Q2,format='(T57,F)'           ;Run Q2 (0=No, 1 = YES):
Readf, lun, Q3,format='(T57,F)'           ;Run Q3 (0=No, 1 = YES):
Readf, lun, Q4,format='(T57,F)'           ;Run Q4 (0=No, 1 = YES):
Readf, lun, line
Readf, lun, line;---------------------------------------------------------
Readf, lun, line;---------------------------------------------------------
Readf, lun, line;           Caution!!!! Caution!!!! Caution!!!! Caution!!!!Caution!!!!
Readf, lun, line;           Everything below does not need to be changed in normal operation     
Readf, lun, line;         (Advanced users can change settings but read the README if interested)    
Readf, lun, line;---------------------------------------------------------
Readf, lun, line;---------------------------------------------------------
Readf, lun, line
Readf, lun, line;---------------------------------------------------------
Readf, lun, line;Advanced Setting for LOFAR DATA
Readf, lun, line;---------------------------------------------------------
Readf, lun, line;Note: It was discovered that the first and last channel
Readf, lun, line;      in each subband was causing problems. Therefore by
Readf, lun, line;      default we eliminate them during the analysis
Readf, lun, line;---------------------------------------------------------
Readf, lun, REMOVE_BADF,format='(T57,F)'  ;Remove Bad Frequencies (0=NO, 1=YES):
Readf, lun, line
Readf, lun, line;-------------------------------------------------------------
Readf, lun, line;  Advanced Settings for Frequency Selection for RFI/Processing
Readf, lun, line;-------------------------------------------------------------
Readf, lun, line;Note: Changing the Frequency Settings for RFI/processing will
Readf, lun, line;      cause problems because of the elimination of the bad channels
Readf, lun, line;      at the start and end of each subband
Readf, lun, line;Fmin and Fmax by default are taken from the header
Readf, lun, line;-------------------------------------------------------------
Readf, lun, SELECT_Freq,format='(T57,F)';Change Frequency Settings (0=NO, 1=YES) (Default: NO):
Readf, lun, fmin,format='(T57,F)'       ;Min Frequency (MHz):
Readf, lun, fmax,format='(T57,F)'       ;Max Frequency (MHz):
Readf, lun, line
Readf, lun, line;-------------------------------------------------------------
Readf, lun, line;                 RFI Advanced Settings
Readf, lun, line;-------------------------------------------------------------
Readf, lun, quantile_norm,format='(T57,F)'  ;
Readf, lun, quantile,format='(T57,F)'       ;
Readf, lun, mask_save,format='(T57,F)'      ;
Readf, lun, line
Readf, lun, line;-------------------------------------------------------------
Readf, lun, line;              RFI Advanced Settings 2
Readf, lun, line;            (Parameters for PATROL,LE_SIG,SUM,PEX)
Readf, lun, line;-------------------------------------------------------------
Readf, lun, NFREQ_PATROL,format='(T57,F)'      ;
Readf, lun, NTIME_PATROL,format='(T57,F)'
Readf, lun, NFREQ_LESIG,format='(T57,F)'
Readf, lun, NTIME_LESIG,format='(T57,F)'
Readf, lun, RUN_SUM_dir,format='(T57,F)'           ;Run Sliding Windows files in directory (0=NO, 1=YES)   
Readf, lun, MF_file,format='(T57,A)'
Readf, lun, MT_file,format='(T57,A)'
Readf, lun, THRESHOLD,format='(T57,F)'
Readf, lun, PMIN_F,format='(T57,F)'
Readf, lun, EXPF,format='(T57,F)'
Readf, lun, PMIN_T,format='(T57,F)'
Readf, lun, EXPT,format='(T57,F)'
Readf, lun, PMIN_F2,format='(T57,F)'
Readf, lun, EXPF2,format='(T57,F)'
Readf, lun, PMIN_T2,format='(T57,F)'
Readf, lun, EXPT2,format='(T57,F)'        
Readf, lun, line
Readf, lun, line;----------------------------------------------------------------------------
Readf, lun, line;            RFI Debug Tools
Readf, lun, line;----------------------------------------------------------------------------
Readf, lun, RFI,format='(T81,I)'
Readf, lun, full_clean_save,format='(T81,F)' ;
Readf, lun, save_datagain,format='(T81,F)'   ;
Readf, lun, COMPLEX_NAME,format='(T81,I)'   ;Complex name for outputs includes RFI settings (0=No, 1 = YES) (Default: NO):   0
Readf, lun, line;
Readf, lun, line;----------------------------------------------------------------------------
Readf, lun, line;          De-dispersion Advanced Settings (After De-dispersion)
Readf, lun, line;----------------------------------------------------------------------------
Readf, lun, RUN_DISP_TIME,format='(T81,F)'
Readf, lun, rebin_time_disp,format='(T81,F)' ;
Readf, lun, line
Readf, lun, line;----------------------------------------------------------------------------
Readf, lun, line;FFT advanced Setting
Readf, lun, line;----------------------------------------------------------------------------
Readf, lun, ddm,format='(T81,F)' ;Change in dm (pc^3 cm) (Default: 0.0):
Readf, lun, line
Readf, lun, line;----------------------------------------------------------------------------
Readf, lun, line;     Combined Mask: Processing Settings       
Readf, lun, line;----------------------------------------------------------------------------
Readf, lun, create_combine_mask,format='(T81,F)' 
Readf, lun, use_combinemask,format='(T81,F)'
Readf, lun, NUM_MASKS,format='(T81,F)'   
 If NUM_MASKS eq 0 then Readf, lun, mask1,format='(T81,F)'    ;What masks do you want to combine together (e.g. 0,1):
 If NUM_MASKS eq 1 then Readf, lun, mask1,format='(T81,F)'    ;What masks do you want to combine together (e.g. 0,1):
 If NUM_MASKS eq 2 then begin 
   Readf, lun, mask1,mask2,format='(T81,F,F)'                 ;What masks do you want to combine together (e.g. 0,1):
   COMBO_BEAMS = [mask1,mask2]
 endif 
 If NUM_MASKS eq 3 then begin 
   Readf, lun, mask1,mask2,mask3,format='(T81,F,F,F)'         ;What masks do you want to combine together (e.g. 0,1):
   COMBO_BEAMS = [mask1,mask2,mask3]
 endif 
 If NUM_MASKS eq 4 then begin
   Readf, lun, mask1,mask2,mask3,mask4,format='(T81,F,F,F,F)'  ;What masks do you want to combine together (e.g. 0,1):
   COMBO_BEAMS = [mask1,mask2,mask3,mask4]
 endif 
Readf, lun, line
Readf, lun, line;----------------------------------------------------------------------------
Readf, lun, line;              Additional Advanced Processing Settings 
Readf, lun, line;----------------------------------------------------------------------------
Readf, lun, Smooth_Surface,format='(T81,F)'   ;Smooth the 2-d Time-Frequency Surface (0=No, 1 = YES) (Default: NO):
Readf, lun, N_smooth,format='(T81,F)'         ; Number of spectra to smooth surface over (Default: 11):
Readf, lun, PCA,format='(T81,I)'              ;Run PCA analysis (-- under development) (0=No, 1 = YES) (Default: NO):     0
Readf, lun, line
Readf, lun, line;----------------------------------------------------------------------------
Readf, lun, line;                     Debugging Tools for Processsing
Readf, lun, line;----------------------------------------------------------------------------
Readf, lun, skip_rebinloop,format='(T81,F)'     ;Skip Rebin Step in Time-Freq Correction (0=No, 1=YES) (Default: NO):            0
Readf, lun, save_RFIMaskData,format='(T81,F)'   ;Save Rebinned Data in Time-Freq Correction (0=No, 1=YES) (Default: NO):         0
Readf, lun, line;
Readf, lun, line;----------------------------------------------------------------------------
Readf, lun, line;                   Post-Processing Advanced Settings
Readf, lun, line; Note:Post-processing is very complicated
Readf, lun, line;       See Turner et al. 2019, 2021 A&A and the README for more details
Readf, lun, line;----------------------------------------------------------------------------
Readf, lun, Run_Rebintime2,format='(T113,F)'
Readf, lun, rebin_time2,format='(T113,F)'
Readf, lun, Run_RebinFreq2,format='(T113,F)'
Readf, lun, rebin_freq2,format='(T113,F)'
Readf, lun, TI,format='(T113,F)'
Readf, lun, Default_MASK,format='(T113,F)'
Readf, lun, Run_Mask,format='(T113,F)'
Readf, lun, threshold_mask,format='(T113,F)'
Readf, lun, Default_FreqRange,format='(T113,F)'
Readf, lun, Run_FreqRange,format='(T113,F)'
Readf, lun, freq_minextra1,freq_maxextra2,format='(T113,F,F)'
Readf, lun, Default_Time,format='(T113,F)'
Readf, lun, QTmin,QTmax,format='(T113,F,F)'
Readf, lun, Default_Threshold,format='(T113,F)'
Readf, lun, threshold_extra,format='(T113,F)' 
Readf, lun, max_threshold,format='(T113,F)'
Readf, lun,threshold_p,format='(T113,F)'
Readf, lun,tau_Qplot,format='(T113,F)'
Readf, lun,slope,format='(T113,F)'
Readf, lun,Qwin,format='(T113,F)'
Readf, lun,ElCorrect,format='(T113,F)'
Readf, lun,NOElCorrect,format='(T113,F)'
file_loc_all = FILE_WHICH('BOREALISv4.pro')
size = STRLEN(file_loc_all) 
root_code = STRMID(file_loc_all,size,size-14,/REVERSE_OFFSET)
print, 'Root Code', Root_code
;Readf, lun,root_code ,format='(T113,A)'                           ;Folder name for location of BOREALISv4 (e.g. /home/jake.turner/codes/BOREALISv4/)
root_gauss = root_code+'/Secondary/GAUSS/'
Readf, lun,save_Q,format='(T113,F)'
Readf, lun,OnBeam_Label,format='(T113,A)'
Readf, lun,OffBeam_Label,format='(T113,A)'
Readf, lun,Extra_Beam1,format='(T113,F)'
Readf, lun,Extra_Beam2,format='(T113,F)'
Readf, lun,REMOVE,format='(T113,F)'
Readf, lun,PAPERPLOTS,format='(T113,F)'
Readf, lun,MINT,format='(T113,F)'
Readf, lun,MAXT,format='(T113,F)'
Readf, lun, line
Readf, lun, line;----------------------------------------------------------------------------
Readf, lun, line; Post-Processing Advanced Settings for Combining Dates
Readf, lun, line; (Note: Combining Dates requires some special steps, See README) -- under development
Readf, lun, line;----------------------------------------------------------------------------
Readf, lun, file_root_combine,format='(T73,A)'       ; 
Readf, lun, filenname_combinedates,format='(T73,A)'  ;
Readf, lun, create_dates_combine,format='(T73,F)'    ; 
Readf, lun, Combine_NBeams,format='(T73,F)'          ;
If Combine_NBeams eq 0 then Readf, lun, beam1,format='(T73,F)'
 If Combine_NBeams eq 1 then Readf, lun, beam1,format='(T73,F)'
 If Combine_NBeams eq 2 then begin
  Readf, lun, beam1,beam2,format='(T73,F,F)'
  Combine_Beams = [beam1,beam2]
 endif
 If Combine_NBeams eq 3 then begin
  Readf, lun, beam1,beam2,beam3,format='(T73,F,F,F)'
  Combine_Beams = [beam1,beam2,beam3]
 endif
 If Combine_NBeams eq 4 then begin
  Readf, lun, beam1,beam2,beam3,beam4,format='(T73,F,F,F,F)' 
  Combine_Beams = [beam1,beam2,beam3,beam4]
 endif
Readf, lun,dates_combine,format='(T73,F)';line248
Readf, lun, line;
Readf, lun, line;----------------------------------------------------------------------------
Readf, lun, line;----------------------------------------------------------------------------
Readf, lun, line;Caution!!!! Caution!!!! Caution!!!! Caution!!!!
Readf, lun, line;   Everything below here are very advanced features to the code
Readf, lun, line;    Normal users will never deal with these parameters
Readf, lun, line;----------------------------------------------------------------------------
Readf, lun, line;----------------------------------------------------------------------------
Readf, lun, line
Readf, lun, line;----------------------------------------------------------------------------
Readf, lun, line;Very Advanced Settings for Jupiter Emission Tests 
Readf, lun, line;  (see Jupiter_Input.dat for more inputs)
Readf, lun, line;----------------------------------------------------------------------------
Readf, lun,RUN_JUPITER,format='(T57,F)'
Readf, lun,input_file_J,format='(T57,A)'
 Free_Lun, lun

;----------------------------------------------------
;                   NENUFAR
;If selected to run full file for Quicklook or RFI 
;----------------------------------------------------
If keyword_set(NENUFAR) then begin
  temp       = filename+'*_'+strtrim(string(start_beam),1)+'.spectra*'
  file_input = findfile(temp)
  print, file_input
  If keyword_set(RUN_QUICKLOOK) then beam_in = bmin_ql ELSE beam_in = start_beam
  If S eq 0 then Stokes = 0 ;Stokes-I
  If S eq 1 then Stokes = 1 ;Stokes-Q
  If S eq 2 then Stokes = 2 ;Stokes-U
  If S eq 3 then Stokes = 3 ;Stokes-V
  If S eq 4 then Stokes = 3 ;Stokes-V^2
  If S eq 5 then Stokes = 3 ;Stokes-V prime
  If S eq 0 then nstokes = 1 ELSE nstokes = 4
  READ_NU_FITS, file_input, command,param,variab, nt,dt,nf,df,ns,jd0,h0,fref,$
    temporary(data3),xt,xf,beam,ndata,corrt,corrf,/nodata,/quiet
  sampling_time = dt

  ;----------
  ;QuickLook 
  ;----------
  If N_Quick_full eq 1 then begin
    tmin_ql = 0
    min_t =param.tmin
    max_t =param.tmax
    steps_ql = FLOOR((max_t-tmin_ql)/N_time_ql)

    ;change the amount of time per step to not lose any data
    N_time_ql = (max_t-tmin_ql)/steps_ql
  Endif
  
  ;-----
  ;RFI
  ;-----
  If Full_RFI eq 1 then begin
    namefits = file_input
    FIND_SLICES, namefits, wbegin, wend, nslices
    ;max_t =param.tmax
    ;nt_new = FLOOR(((max_t - tmin)/dt))
    ;steps  = round(nt_new/nspectra)
    ;nspectra = round(nt_new/steps)
    
    steps = nslices
    print,'*******************************'
    Print,'Running full file'
    Print,'Number of steps', steps
    ;Print,'Number of spectra', nspectra
    print,'*******************************'
  Endif
  ;    read_nenufar, file_root+filename, temporary(data2),tt,ff,beam_in,ndata, ntt,dt,nff,df,ns,Stokes,tmin=0,tmax=1,fmin=9,fmax=70,nstokes=nstokes,/VERBOSE,MINFREQ=min_freq,MAXFREQ=max_freq
Endif  ; end NENUFAR

;-----------------------------------------------
;By Default; fmin and fmax are taken from header
;-----------------------------------------------
If not(keyword_set(SELECT_Freq)) then begin
  ;B_string = strtrim(string(start_beam),1)
   If S eq 4 then S_in = 3   ; V^2
  If S eq 5 then S_in = 3   ; V'
  If S eq 6 then S_in = 1   ; L = sqrt(Q^2 + U^2) 
  If S lt 4 then S_in = S   ; S = S (I, U, Q, V)
  If keyword_set(LOFAR) then begin
    file = file_root+date+'/'+date+'_SAP000_B00'+strtrim(string(uint(start_beam)),1)+'_S'+strtrim(string(uint(S_in)),1)+'_P000_bf'
    READ_H5_DATA, file,0,1,fmin,fmax,start_Beam,S_in,x,t,xf,nt,nf,/NSPECTRA,min_freq=min_freq,max_freq = max_freq,/READ_HEADER
  Endif
  If keyword_set(NENUFAR) then begin
    temp       = filename+'*_'+strtrim(string(start_beam),1)+'.spectra*'
    file_input = findfile(temp)
    print, file_input
    If keyword_set(RUN_QUICKLOOK) then beam_in = bmin_ql ELSE beam_in = start_beam
    If S eq 0 then Stokes = 0 ;Stokes-I
    If S eq 1 then Stokes = 1 ;Stokes-Q
    If S eq 2 then Stokes = 2 ;Stokes-U
    If S eq 3 then Stokes = 3 ;Stokes-V
    If S eq 4 then Stokes = 3 ;Stokes-V^2
    If S eq 5 then Stokes = 3 ;Stokes-V prime
    If S eq 0 then nstokes = 1 ELSE nstokes = 4
    READ_NU_FITS, file_input, command,param,variab, nt,dt,nf,df,ns,jd0,h0,fref,$
      temporary(data3),xt,xf,beam,ndata,corrt,corrf,/nodata,/quiet ;tmax=tmax
    sampling_time = dt
    min_freq =param.FMIN
    max_freq = param.FMAX
  Endif
  If keyword_Set(DYNSPECMS) then begin
    print, S
    If S eq 0 then POL = 1 ;Stokes-I
    If S eq 1 then POL = 2 ;Stokes-Q
    If S eq 2 then POL = 3 ;Stokes-U
    If S eq 3 then POL = 4 ;Stokes-V
    If S eq 4 then POL = 4 ;Stokes-V^2
    If S eq 5 then POL = 4 ;Stokes-V prime
    print,'filename for DYNSPECMS: ',filename
    read_dynspec_fits,filename+'_beam'+strtrim(string(uint(start_beam)),1),0,1,$
    fmin,fmax,POL,x,xt,xf,nt,nf,fmin=min_freq,fmax=max_freq,/PRINT_HEAD,nnt=nnt
  Endif
  fmin = min_freq
  fmax = max_freq 
Endif

;Full RFI for DYNSPECMS
If keyword_Set(DYNSPECMS) then begin
  If S eq 0 then POL = 1 ;Stokes-I
  If S eq 1 then POL = 2 ;Stokes-Q
  If S eq 2 then POL = 3 ;Stokes-U
  If S eq 3 then POL = 4 ;Stokes-V
  If S eq 4 then POL = 4 ;Stokes-V^2
  If S eq 5 then POL = 4 ;Stokes-V prime
  read_dynspec_fits,filename+'_beam'+strtrim(string(uint(start_beam)),1),0,1,$
    fmin,fmax,POL,x,xt,xf,nt,nf,fmin=min_freq,fmax=max_freq,/HEADER_ONLY,tmax=tmax,nnt=nnt
  If N_Quick_full eq 1 then begin
    N_time_ql = (tmax)/(steps_ql)
  Endif 
  
 If Full_RFI eq 1 then begin
  nspectra = nnt
  print,'*******************************'
  Print,'Running full file'
  Print,'Number of steps', steps
  Print,'Number of spectra', nspectra
  print,'*******************************'
 Endif
endif

;-----------------------------
;Setup RFI MF and MT for SUM
;-----------------------------
If keyword_set(RUN_SUM_dir) then begin
  CD, CURRENT=c
  dir_SUM = c
READCOL,dir_SUM+'/'+MF_file,MF
READCOL,dir_SUM+'/'+TF_file,MT
Endif
If not(keyword_set(RUN_SUM_dir)) then begin
  size = STRLEN(root_code)
  dir = STRMID(root_code,size,size-5,/REVERSE_OFFSET)
  dir_SUM = dir+'sub/'
  READCOL,dir_SUM+'mf.dat',MF
  READCOL,dir_SUM+'mt.dat',MT
Endif

;filename = file_root+filename



;*********************
;Create save_file_name
;*********************
;save_create,N_time,steps,patrol_value,le_sig_value,beam,RFI=RFI,patrol=patrol,le_sig=le_sig,sum=sum,pex=pex,$
;  save_file_name=save_file_name
;binTime = BIN_DATE(systime(/ut))  ;
;timestamp_string = TIMESTAMP(YEAR = binTime[0], MONTH = binTime[1],$
;  DAY = binTime[2], HOUR = binTime[3], MINUTE = binTime[4], SECOND =binTime[5])
;save,/ALL,filename='Inputs_'+timestamp_string+'.sav'

;Readf, lun,NOElCorrect

;Update elliptical correction flags 
If ElCorrect eq 0 and NOElCorrect eq 1 then ElCorrect=3  ;Run non-elliptical correction but no elliptical 
If ElCorrect eq 1 and NOELCorrect eq 0 then ElCorrect=0  ;Run elliptical and not non-ell
If ElCorrect eq 1 and NOELCorrect eq 1 then ElCorrect=1  ;Run both
return

end
