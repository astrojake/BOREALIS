;NAME: 
;    read_input
;    
;INPUT: filename with inputs (Input.dat) ;example file is supplied with the code 
;
;OUTPUT: 
;
; MODIFICATION HISTORY: Version 3 (June 14, 2017; Added Process_METHOD2) 
;                       Version 4 (March 22, 2019; Added NENUFAR and DynspecBF, updated all and removed unnessarcy or redunant keywords) 
pro read_input,input_file,planet=planet,period=period,T0=T0,date=date,fileroot=file_root,$
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
               LOFAR=LOFAR,NENUFAR=NENUFAR,DYNSPECBF=DYNSPECBF,filename=filename
             
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
T0                = 0.0d    & date              =''
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
mask1 = 0
mask2 = 0
mask3 = 0
mask4 = 0
RFI  = 0 
LOFAR = 0 & NENUFAR = 0 & DYNSPECBF = 0 
input_file_J = ''  & filename = ''
REMOVE = 0 & PAPERPLOTS = 0 & MINT  = 0.0 & MAXT = 0.0 

;**********************
;Read Variables
;**********************
i =0 ;reset (counting lines)
a = ''
Readf, lun, line ;**************** ; line1
Readf, lun, line ;**************** ; line2
Readf, lun, line ;**************** ; line3
Readf, lun, line ;**************** ; line4
Readf, lun, line ;**************** ; line5
Readf, lun, line ;**************** ; line6
Readf, lun, line ;**************** ; line7
Readf, lun, line ;**************** ; line8
Readf, lun, planet,format='(T57,A)'; line9
Readf, lun, period,format='(T57,D)';line10
Readf, lun, T0,    format='(T57,D)';line11
Readf, lun, date,  format='(T57,A)';line12
Readf, lun, file_root,format='(T57,A)';line13
Readf, lun, line ;**************** ;line14
Readf, lun, line ;**************** ;line15
Readf, lun, line ;****************;line16
Readf, lun, line ;****************;line17
Readf, lun, LOFAR,format='(T57,A)';***************;
Readf, lun, NENUFAR,format='(T57,A)' ;***************;
Readf, lun, DYNSPECBF,format='(T57,A)' ;***************;
Readf, lun, line ;**************** ;line
Readf, lun, line ;**************** ;line
Readf, lun, line ;****************;line
Readf, lun, line ;****************;line
Readf, lun, filename,format='(T57,A)' ;***************;
Readf, lun, line ;**************** ;line
Readf, lun, line ;**************** ;line
Readf, lun, line ;****************;line
Readf, lun, line ;****************        ;line
Readf, lun, run_QUICKLOOK,format='(T57,F)' ;line18
Readf, lun, run_RFI,format='(T57,F)'       ;line19
Readf, lun, run_processing,format='(T57,F)';line20
Readf, lun, run_FFT,format='(T57,F)'       ;line21
Readf, lun, run_DISP,format='(T57,F)'      ;line22
Readf, lun, run_postprocessing,format='(T57,F)';line23
Readf, lun, line ;**************** ;line24  
Readf, lun, line ;**************** ;line25
Readf, lun, line ;****************;line26
Readf, lun, line ;****************;line27
Readf, lun, PLOT_RAW,format='(T57,F)' ;line28
Readf, lun, DATA_NORM ,format='(T57,F)';line29
Readf, lun, DATA_SLOPE,format='(T57,F)';line30
Readf, lun,fmin_ql,format='(T57,F)'    ;line31
Readf, lun, fmax_ql,format='(T57,F)'   ;line32
Readf, lun, tmin_ql,format='(T57,F)'   ;line32
Readf, lun, N_time_ql,format='(T57,F)' ;line34
Readf, lun, steps_ql,format='(T57,F)'  ;line35
Readf, lun, bmin_ql,format='(T57,I)'   ;line36
Readf, lun, bmax_ql,format='(T57,I)'   ;lin37
Readf, lun, Smin_ql,format='(T57,I)'   ;line38
Readf, lun, Smax_ql,format='(T57,I)'   ;line39
Readf, lun, line ;****************     ;line40
Readf, lun, line ;***************      ;line41
Readf, lun, line ;***************      ;line42 
Readf, lun, line ;***************      ;line43
Readf, lun, line;***************       ;line44
Readf, lun, line;***************       ;line45
Readf, lun, line;***************       ;line46
Readf, lun, line;***************       ;line47
Readf, lun, start_beam,format='(T57,I)';line48
Readf, lun, end_beam,format='(T57,I)'  ;line49
Readf, lun, line ;****************   ;line50
Readf, lun, line ;***************    ;line51
Readf, lun, line ;***************    ;line52
Readf, lun, line;***************    ;line53
Readf, lun, S,format='(T57,I)'  ;line54
Readf, lun, line ;****************   ;line55
Readf, lun, line ;***************    ;line56
Readf, lun, line ;***************    ;line57
Readf, lun, line ;***************    ;line58
Readf, lun, tmin,format='(T57,F)'    ;line59
Readf, lun, nspectra,format='(T57,F)';line60
Readf, lun, steps,format='(T57,F)';line61
Readf, lun, line ;**************** 62
Readf, lun, line ;***************  63
Readf, lun, line ;***************  64
Readf, lun, line ;***************  65
Readf, lun, PS,format='(T57,F)'    ;line66
Readf, lun, VERBOSE,format='(T57,F)';line67
Readf, lun, mask_bit,format='(T57,F)' ;line68
Readf, lun, bit,format='(T57,F)'      ;line69
Readf, lun, send_email,format='(T57,F)';line70
Readf, lun, email_address,format='(T57,A)';line71
Readf, lun, line ;****************      ;line72
Readf, lun, line ;***************       ;line73
Readf, lun, line ;***************       ;line74
Readf, lun, line ;***************       ;line75
Readf, lun, patrol_value,format='(T57,F)' ;line76
Readf, lun, le_sig_value,format='(T57,F)' ;line77
Readf, lun, patrol,format='(T57,F)'       ;line78
Readf, lun, le_sig,format='(T57,F)'       ;line79
Readf, lun, sum,format='(T57,F)'          ;line80
Readf, lun, pex,format='(T57,F)'          ;line81
Readf, lun, line ;****************
Readf, lun, line ;***************
Readf, lun, line ;***************
Readf, lun, line ;***************
Readf, lun, TimeFreq_correction,format='(T57,F)' ;line86
Readf, lun, Apply_Corrections,format='(T57,F)'   ;line87
readf,lun,Process_METHOD2,format='(T57,I)'
Readf, lun, line ;****************
Readf, lun, line ;***************
Readf, lun, line ;***************
Readf, lun, line ;***************
Readf, lun, rebin_data,format='(T57,F)' ;line92
Readf, lun, rebin_freq,format='(T57,F)' ;line93
Readf, lun, rebin_time,format='(T57,F)' ;line94
Readf, lun, line ;****************
Readf, lun, line ;***************
Readf, lun, line ;***************
Readf, lun, line ;***************
Readf, lun, beam_DISP,format='(T57,F)' ;line99
Readf, lun, pulsarname,format='(T57,A)';line100
Readf, lun, dm,format='(T57,F)'        ;line101
Readf, lun, size_disp,format='(T57,F)' ;line102
Readf, lun, line
Readf, lun, line ;****************
Readf, lun, line ;***************
Readf, lun, line ;***************
Readf, lun, beam_FFT,format='(T57,F)'  ;line107
Readf, lun, Period_FFT,format='(T57,F)'
Readf, lun, pmin,format='(T57,F)'
Readf, lun, fminp,format='(T57,F)'
Readf, lun, fmaxp,format='(T57,F)'    ;line111
Readf, lun, line;****************
Readf, lun, line ;***************
Readf, lun, line ;***************
Readf, lun, line ;***************
Readf, lun, line;***************
Readf, lun, Target_Beam,format='(T57,I)';line117
Readf, lun, Sky_Beam,format='(T57,I)'
Readf, lun, Q1,format='(T57,F)'
Readf, lun, Q2,format='(T57,F)'
Readf, lun, Q3,format='(T57,F)'
Readf, lun, Q4,format='(T57,F)';line122
Readf, lun, line;***************
Readf, lun, line;***************
Readf, lun, line;***************
Readf, lun, line;***************
Readf, lun, line;***************
Readf, lun, line;***************
Readf, lun, line;***************
Readf, lun, line;***************;line130
Readf, lun, line;***************
Readf, lun, line;***************
Readf, lun, line;***************
Readf, lun, line;***************
Readf, lun, line;***************
Readf, lun, line;***************
Readf, lun, line;***************
Readf, lun, line;***************
Readf, lun, REMOVE_BADF,format='(T57,F)';line139
Readf, lun, line;***************;line140
Readf, lun, line;***************
Readf, lun, line;***************
Readf, lun, line;***************
Readf, lun, line;***************
Readf, lun, line;***************
Readf, lun, line;***************
Readf, lun, line;***************
Readf, lun, line;***************
Readf, lun, SELECT_Freq,format='(T57,F)';line149
Readf, lun, fmin,format='(T57,F)'
Readf, lun, fmax,format='(T57,F)'       ;line151
Readf, lun, line;***************
Readf, lun, line;***************
Readf, lun, line;***************
Readf, lun, line;***************
Readf, lun, quantile_norm,format='(T57,F)' ;;line156
Readf, lun, quantile,format='(T57,F)'
Readf, lun, mask_save,format='(T57,F)'    ;line158
Readf, lun, line;***************
Readf, lun, line;***************
Readf, lun, line;***************
Readf, lun, line;***************
Readf, lun, line;***************
Readf, lun, NFREQ_PATROL,format='(T57,F)'      ;line164
Readf, lun, NTIME_PATROL,format='(T57,F)'
Readf, lun, NFREQ_LESIG,format='(T57,F)'
Readf, lun, NTIME_LESIG,format='(T57,F)'
Readf, lun, MF,format='(T57,A)'
Readf, lun, MT,format='(T57,A)'
Readf, lun, THRESHOLD,format='(T57,F)'
Readf, lun, PMIN_F,format='(T57,F)'
Readf, lun, EXPF,format='(T57,F)'
Readf, lun, PMIN_T,format='(T57,F)'
Readf, lun, EXPT,format='(T57,F)'
Readf, lun, PMIN_F2,format='(T57,F)'
Readf, lun, EXPF2,format='(T57,F)'
Readf, lun, PMIN_T2,format='(T57,F)'
Readf, lun, EXPT2,format='(T57,F)'     ;line178
Readf, lun, line;***************
Readf, lun, line;***************
Readf, lun, line;***************
Readf, lun, line;***************
Readf, lun, RFI,format='(T81,I)'
Readf, lun, full_clean_save,format='(T81,F)' ;line183
Readf, lun, save_datagain,format='(T81,F)'   ;line184
Readf, lun, line;***************
Readf, lun, line;***************
Readf, lun, line;***************
Readf, lun, line;***************
Readf, lun, RUN_DISP_TIME,format='(T81,F)'
Readf, lun, rebin_time_disp,format='(T81,F)' ;line190
Readf, lun, line;***************
Readf, lun, line;***************
Readf, lun, line;***************
Readf, lun, line;***************
Readf, lun, ddm,format='(T81,F)' ;line195
Readf, lun, line;***************
Readf, lun, line;***************
Readf, lun, line;***************
Readf, lun, line;***************
Readf, lun, create_combine_mask,format='(T81,F)' ;line200
Readf, lun, use_combinemask,format='(T81,F)'
Readf, lun, NUM_MASKS,format='(T81,F)'   
 If NUM_MASKS eq 0 then Readf, lun, mask1,format='(T81,F)'
 If NUM_MASKS eq 1 then Readf, lun, mask1,format='(T81,F)'
 If NUM_MASKS eq 2 then begin 
   Readf, lun, mask1,mask2,format='(T81,F,F)'
   COMBO_BEAMS = [mask1,mask2]
 endif 
 If NUM_MASKS eq 3 then begin 
   Readf, lun, mask1,mask2,mask3,format='(T81,F,F,F)'
   COMBO_BEAMS = [mask1,mask2,mask3]
 endif 
 If NUM_MASKS eq 4 then begin
   Readf, lun, mask1,mask2,mask3,mask4,format='(T81,F,F,F,F)'  ;line203
   COMBO_BEAMS = [mask1,mask2,mask3,mask4]
 endif 
Readf, lun, line;***************
Readf, lun, line;***************
Readf, lun, line;***************
Readf, lun, line;***************
Readf, lun, Smooth_Surface,format='(T81,F)'   ;line208
Readf, lun, N_smooth,format='(T81,F)'
Readf, lun, line;***************
Readf, lun, line;***************
Readf, lun, line;***************
Readf, lun, line;***************
Readf, lun,skip_rebinloop,format='(T81,F)' ;line214
Readf, lun, save_RFIMaskData,format='(T81,F)';line215
Readf, lun, line;***************
Readf, lun, line;***************
Readf, lun, line;***************
Readf, lun, line;***************
Readf, lun, line;***************;line220
Readf, lun, line;***************;line221
Readf, lun, Run_Rebintime2,format='(T113,F)'
Readf, lun, rebin_time2,format='(T113,F)'
Readf, lun, Run_RebinFreq2,format='(T113,F)'
Readf, lun, rebin_freq2,format='(T113,F)';line237
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
Readf, lun, threshold_extra,format='(T113,F)' ;line232
Readf, lun, max_threshold,format='(T113,F)'
Readf, lun,threshold_p,format='(T113,F)'
Readf, lun,tau_Qplot,format='(T113,F)'
Readf, lun,slope,format='(T113,F)'
Readf, lun,Qwin,format='(T113,F)'
Readf, lun,ElCorrect,format='(T113,F)'
Readf, lun,NOElCorrect,format='(T113,F)'
Readf, lun,root_gauss ,format='(T113,A)'
Readf, lun,save_Q,format='(T113,F)'
Readf, lun,OnBeam_Label,format='(T113,A)'
Readf, lun,OffBeam_Label,format='(T113,A)'
Readf, lun,Extra_Beam1,format='(T113,F)'
Readf, lun,Extra_Beam2,format='(T113,F)'
Readf, lun,REMOVE,format='(T113,F)'
Readf, lun,PAPERPLOTS,format='(T113,F)'
Readf, lun,MINT,format='(T113,F)'
Readf, lun,MAXT,format='(T113,F)'
Readf, lun, line;***************
Readf, lun, line;***************
Readf, lun, line;***************line240
Readf, lun, line;***************
Readf, lun, line;***************line242
Readf, lun, file_root_combine,format='(T73,A)'
Readf, lun, filenname_combinedates,format='(T73,A)'
Readf, lun, create_dates_combine,format='(T73,F)'
Readf, lun, Combine_NBeams,format='(T73,F)'
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
Readf, lun, line;***************
Readf, lun, line;***************
Readf, lun, line;***************
Readf, lun, line;***************
Readf, lun, line;***************
Readf, lun, line;***************
Readf, lun, line;***************
Readf, lun, line;***************
Readf, lun, line;***************
Readf, lun, line;***************
Readf, lun, line;***************
Readf, lun, line;***************
Readf, lun, line;***************
Readf, lun,RUN_JUPITER,format='(T57,F)'
Readf, lun,input_file_J,format='(T57,A)'
 Free_Lun, lun

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
    If keyword_set(RUN_QUICKLOOK) then beam_in = bmin_ql ELSE beam_in = start_beam
    If S eq 0 then Stokes = 0 ;Stokes-I
    If S eq 1 then Stokes = 1 ;Stokes-Q
    If S eq 2 then Stokes = 2 ;Stokes-U
    If S eq 3 then Stokes = 3 ;Stokes-V
    If S eq 4 then Stokes = 3 ;Stokes-V^2
    If S eq 5 then Stokes = 3 ;Stokes-V prime
    If S eq 0 then nstokes = 1 ELSE nstokes = 4
    read_nenufar, file_root+filename, temporary(data2),tt,ff,beam_in,ndata, ntt,dt,nff,df,ns,Stokes,tmin=0,tmax=1,fmin=9,fmax=70,nstokes=nstokes,/VERBOSE,MINFREQ=min_freq,MAXFREQ=max_freq
  Endif
  If keyword_Set(DYNSPECBF) then begin
    If S eq 0 then POL = 1 ;Stokes-I
    If S eq 1 then POL = 2 ;Stokes-Q
    If S eq 2 then POL = 3 ;Stokes-U
    If S eq 3 then POL = 4 ;Stokes-V
    If S eq 4 then POL = 4 ;Stokes-V^2
    If S eq 5 then POL = 4 ;Stokes-V prime
    read_dynspec_fits,file_root+filename+'_beam'+strtrim(string(uint(start_beam)),1),0,1,fmin,fmax,POL,x,xt,xf,nt,nf,min_freq=min_freq,max_freq=max_freq
  Endif
  fmin = min_freq
  fmax = max_freq 
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
If ElCorrect eq 0 and NOElCorrect eq 1 then ElCorrect=3  ;Run non-elliptical correction but no elliptical 
If ElCorrect eq 1 and NOELCorrect eq 0 then ElCorrect=0  ;Run elliptical and not non-ell
If ElCorrect eq 1 and NOELCorrect eq 1 then ElCorrect=1  ;Run both
return

end
