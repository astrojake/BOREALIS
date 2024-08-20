;***************************************************************
;      Post-Processing Lofar diagnostic plots
;**************************************************************
;***************************************************************
;;Author: Jake Turner (UVA)
;         J.M. Griessmeier (CNRS) and Philippe Zarka (Obspm)
;        
;Inputs:
;        x               ; Target Observable to plot (On Beam)
;        x2              ; Off beam of observable to plot (Sky Beam)
;        xt              ; Time in Hour 
;        t               ; Time series of Observable (MJD) [as an 2d array for both target and off beam)
;        P               ; Period of planet (Hour) 
;        T0              ; Time of 0 phase for the exoplanet from first transit (MJD)
;        label           ; Name of Observable to include in plot 
;        ytitle           ;ytitle for the graphs (input ytitle)
;        
;Keywords:  
;       t2               ;Time in MJD 
;       ExBeam1           ;Include ExBeam1 in plots with On and Off Beam (t [MJD],xt [Hour; sorted for different dates],x3)
;       ExBeam2           ;Include ExBeam2 Source in plots (t [MJD],xt [Hour; sorted for different dates],x4)
;       /TARGET_ONLY     ;Only Plot Target 
;       Rebin_time       ;Smooth the data using the orginal rebin_time 
;       
;Output: 
;       Plot 1: time variation of x
;       Plot 2: smoothed (boxcar-filtered) time variation of x with 6-second sliding window 
;       Plot 3: Variation of x vs. orbital phase of the planet 
;       Plot 4: Smoothed (boxcar filter) variation of x vs. orbital phase of with short (6 sec) and long (12 secs) sliding windows 
;
;History
;   On OSUC: LIST and old Qs 3/21/17: (Update: Jupiter, Vprime Pol) --- All LIST and Qs programmming is there. 
;   
;   5/23/2017 - All list and Qs deleted until needed
;      
;--------------------------------------------
pro diagnostic_plots,x,x2,t,xt,P,T0,label,ytitle,t2=t2,$
    ExBeam1=ExBeam1,ExBeam2=ExBeam2,TARGET_ONLY=TARGET_ONLY,$
    LIST=LIST,single_plot=single_plot,single_phase=single_phase,$
    Error=Error,ONBEAM_LABEL=ONBEAM_LABEL,OFFBEAM_LABEL=OFFBEAM_LABEL,$
    Q_Gauss=Q_Gauss,NEG_GAUSS=NEG_GAUSS,diff=diff,line=line,UT=UT,min_scale=min_scale,max_scale=max_scale,ALL_Gauss=ALL_Gauss
 ;--------------------   
 ;Set Keywords
 ;---------------------
 dim_exbeam1 = size(ExBeam1,/N_DIMENSIONS)
 If dim_exbeam1 eq 0 then ExBeam1 = 0 Else ExBeam1=ExBeam1
 dim_exbeam2 = size(ExBeam2,/N_DIMENSIONS)
 If dim_exbeam2 eq 0 then ExBeam2 = 0 Else ExBeam2=ExBeam2
 If not(keyword_set(ONBEAM_LABEL)) then ONBEAM_LABEL = 'Target'
 If not(keyword_set(OFFBEAM_LABEL)) then OFFBEAM_LABEL = 'Sky'
 If not(keyword_set(error)) then error = 0.
 If not(keyword_set(max_scale)) then max_scale = 1.10 
 If not(keyword_set(min_scale)) then min_scale = 1.10 
 set_plot,'PS'
 
 
 ;------------------------------------------------
 ;       Set Max and Min Values for Plots 
 ;------------------------------------------------
  char_size = 2.0
  leg_charsize = 2.0
  If not(Keyword_set(TARGET_ONLY)) then begin
   ymin = (min([x,x2])-error[0]*min_scale)  & ymax = (max([x,x2])+error[0])*max_scale
  Endif
  
  If not(Keyword_set(TARGET_ONLY)) and keyword_set(diff) then begin
    ymin = (min([x,x2,x-x2])-error[0]*min_scale)  & ymax = (max([x,x2,x-x2])+error[0])*max_scale
  Endif
  
  If not(Keyword_set(TARGET_ONLY)) and keyword_set(Q_Gauss) and not(keyword_set(NEG_GAUSS)) then begin
    ymin = (min([x,x2])-error[0]*min_scale)  & ymax = max([(max([x,x2])+error[0]),Q_Gauss])*max_scale
  Endif
  If not(Keyword_set(TARGET_ONLY)) and keyword_set(Q_Gauss) and keyword_set(NEG_GAUSS) then begin
    ymin = min([(min([x,x2])-error[0]),-Q_Gauss])*min_scale  & ymax = max([max([x,x2])+error[0],Q_Gauss])*max_scale
  Endif
  If Keyword_set(TARGET_ONLY) then begin
    ymin = (min([x])-error[0]*min_scale)  & ymax = (max(x)+error[0])*max_scale
  Endif
  If Keyword_set(TARGET_ONLY) and keyword_set(Q_Gauss) and keyword_set(NEG_GAUSS) then begin
    ymin = min([min(x)-error[0],-Q_Gauss])*min_scale  & ymax = max([max(x)+error[0],Q_Gauss])*max_scale
  Endif
  
  If keyword_set(ExBeam1) and not(Keyword_set(TARGET_ONLY)) and not(Keyword_set(ExBeam2)) then begin
    ymin = (min([x,x2,ExBeam1[1,*]])-error)*min_scale  & ymax = (max([x,x2,ExBeam1[1,*]])+error)*max_scale
  endif 
  If keyword_set(ExBeam2) and not(Keyword_set(TARGET_ONLY)) and not(Keyword_set(ExBeam1)) then begin
    ymin = (min([x,x2,ExBeam2[1,*]])-error)*min_scale  & ymax = (max([x,x2,ExBeam2[1,*]])+error)*max_scale
  endif
  If keyword_set(ExBeam2) and Keyword_set(ExBeam1) and not(Keyword_set(TARGET_ONLY)) then begin
    ymin = (min([x,x2,reform(ExBeam1[1,*]),reform(ExBeam2[1,*])])-error)*min_scale  & ymax = (max([x,x2,reform(ExBeam1[1,*]),reform(ExBeam2[1,*])])+error)*max_scale
  endif
  
 ;---------------------------------------------------------------------------------
 ;                           Single Plot; time
 ;;--------------------------------------------------------------------------------
 If keyword_Set(single_plot) then begin
;   !p.multi=[1,1,1]
   ;**********************
   ;         Plot 1
   ;**********************
   If keyword_set(line) then psym_plot_target = -14 Else psym_plot_target = 46
   If keyword_set(line) then psym_plot_sky = -15 Else psym_plot_sky = 15
   If not(keyword_set(LIST)) then begin
     ;On Beam Plot
      If not(keyword_set(UT)) then begin 
       cgplot,xt-min(xt),x,err_yhigh=error,err_ylow=error,title=label,xtitle='Time (Hour)',ytitle=ytitle,$
       charsize=char_size,/xstyle,/ystyle,psym=psym_plot_target,yrange=[ymin,ymax],symsize=0.9,thick=2
      endif 
      If keyword_set(UT) then begin 
        cgplot,UT,x,err_yhigh=error,err_ylow=error,title=label,xtitle='UT (hour)',ytitle=ytitle,$
       charsize=char_size,/xstyle,/ystyle,psym=psym_plot_target,yrange=[ymin,ymax],symsize=0.9,thick=2
      endif 

     ;Off Beam
     If not(keyword_set(TARGET_ONLY)) then begin
       If not(keyword_set(UT)) then cgplot,xt-min(xt),x2,err_yhigh=error,err_ylow=error,/overplot,color='red',psym=psym_plot_sky,symsize=0.9
       If keyword_set(UT) then cgplot,UT,x2,err_yhigh=error,err_ylow=error,/overplot,color='red',psym=psym_plot_sky,symsize=0.9
     Endif
     
     ;Difference of ON and OFF 
     If not(keyword_set(TARGET_ONLY)) and keyword_set(diff) then begin
       If not(keyword_set(UT)) then cgplot,xt-min(xt),x-x2,err_yhigh=error*sqrt(2.),err_ylow=error*sqrt(2.),/overplot,color='Green',psym=-46
       If keyword_set(UT) then cgplot,UT,x-x2,err_yhigh=error*sqrt(2.),err_ylow=error*sqrt(2.),/overplot,color='Green',psym=-46
     Endif
     
     ;Plot ExBeam1
     If keyword_set(ExBeam1) then begin
       cgplot,ExBeam1[2,*]-min(ExBeam1[2,*]),ExBeam1[1,*],err_yhigh=error,err_ylow=error,/overplot,color='blue',psym=-16,symsize=0.9
     Endif

     ;Plot ExBeam2
     If keyword_set(ExBeam2) then begin
       cgplot,ExBeam2[2,*]-min(ExBeam2[2,*]),ExBeam2[1,*],err_yhigh=error,err_ylow=error,/overplot,color='Orange',psym=-17,symsize=0.9
     Endif
     
     If keyword_set(Q_Gauss) then begin
      n = n_elements(xt)
      Q_G = dblarr(n)
      Q_G(*) = Q_Gauss ;all the same
      If not(keyword_set(UT)) then cgplot,xt-min(xt),Q_G,/overplot,LineStyle=2
      If keyword_set(UT) then cgplot,UT,Q_G,/overplot,LineStyle=2 
      If keyword_set(ALL_Gauss) then cgplot,xt-min(xt),Q_G*2.0,/overplot,LineStyle=2
      If keyword_set(ALL_Gauss) then cgplot,xt-min(xt),Q_G*3.0,/overplot,LineStyle=2
      If keyword_set(NEG_GAUSS) then begin
        cgplot,xt-min(xt),-Q_G,/overplot,LineStyle=2
      ;  cgplot,xt-min(xt),-Q_G*2.0,/overplot,LineStyle=2
      ;  cgplot,xt-min(xt),-Q_G*3.0,/overplot,LineStyle=2
      Endif
     Endif
     
     If not(keyword_set(UT)) then cgplot,[0,max(xt)],[0,0],LineStyle=0,/overplot
     If keyword_set(UT) then cgplot,[0,max(UT)],[0,0],LineStyle=0,/overplot
   Endif ;no list
   
   ;----------------------------------------
   ;            Legend
   ;----------------------------------------
   ;normal: Target vs Sky
   If not(keyword_set(ExBeam1)) and not(keyword_set(ExBeam2)) and not(keyword_set(TARGET_ONLY)) and not(keyword_set(diff))  then begin
     al_legend,[OnBeam_Label,OffBeam_Label],colors=['Black','red'],psym=[46,15],/right,/top,charsize=leg_charsize,symsize=2
   Endif
   
   ;Target vs Sky and difference 
   If not(keyword_set(ExBeam1)) and not(keyword_set(ExBeam2)) and not(keyword_set(TARGET_ONLY)) and keyword_set(diff) then begin
     al_legend,[OnBeam_Label,OffBeam_Label,'On - Off'],colors=['Black','red','Green'],psym=[46,46,45],/right,/top,charsize=leg_charsize,symsize=2
   Endif
   
   ;Sky, ExBeam1, ExBeam2
   If keyword_set(ExBeam1) and keyword_set(ExBeam2) and not(keyword_set(TARGET_ONLY)) then begin
     al_legend,[OnBeam_Label,OffBeam_Label,'B0823+26 (pulsar)','Bright Source'],psym=[14,15,16,17],$
       colors=['Black','red','blue','orange'],colors=['Black','red','blue','orange'],/left,/top,charsize=leg_charsize,symsize=1
   Endif

   ;Sky, ExBeam1
   If keyword_set(ExBeam1) and not(keyword_set(ExBeam2)) and not(keyword_set(TARGET_ONLY)) then begin
     al_legend,[OnBeam_Label,OffBeam_Label,'Extra Beam 1'],psym=[46,46,46,46],$
       colors=['Black','red','blue'],/right,/top,charsize=leg_charsize,symsize=2
   endif

   ;Sky, ExBeam2
   If keyword_set(ExBeam2) and not(keyword_set(ExBeam1)) and not(keyword_set(TARGET_ONLY)) then begin
     al_legend,[OnBeam_Label,OffBeam_Label,'Extra Beam 2'],psym=[46,46,46,46],$
       textcolor=['Black','red','green'],/right,/top,charsize=leg_charsize,symsize=2
   endif
 Endif ;end single plot
 
 ;---------------------------------------------------------------------------------
 ;                           Single Plot; Phase
 ;;--------------------------------------------------------------------------------
 If keyword_Set(single_phase) then begin
 ;  !p.multi=[1,1,1]
   If not(keyword_set(LIST)) then begin
     ;target plot
     phase1 = (((t - T0) mod P)/(P*1.d))
     s1 = sort(phase1)
     cgplot,phase1[s1],x[s1],title=label,xtitle='Phase ',ytitle=ytitle,$
       charsize=char_size,psym=46,/xstyle,/ystyle,yrange=[ymin,ymax]

     If not(keyword_set(TARGET_ONLY)) then begin
       cgplot,phase1[s1],x2[s1],/overplot,color='red',psym=46
     endif

     If keyword_set(ExBeam1) then begin
       cgplot,phase1[s1],ExBeam1[1,s1],/overplot,color='blue',psym=46
     Endif

     If keyword_set(ExBeam2) then begin
      ; phase4 = (((ExBeam2[0,*] - T0) mod P)/(P*1.d)) ;mod 1
      ; s4 = sort(phase4)
       cgplot,phase1[s1],ExBeam2[1,s1],/overplot,color='green',psym=46
     Endif

     ;--------------------------------
     ;               Legend
     ;--------------------------------
     ;normal: Target vs Sky
     If not(keyword_set(ExBeam1)) and not(keyword_set(ExBeam2)) and not(keyword_set(TARGET_ONLY)) then begin
       al_legend,[OnBeam_Label,OffBeam_Label],textcolor=['Black','red'],psym=[46,46],/right,/top,charsize=leg_charsize,symsize=2
     Endif
     ;Sky, ExBeam1, ExBeam2
     If keyword_set(ExBeam1) and keyword_set(ExBeam2) and not(keyword_set(TARGET_ONLY)) then begin
       al_legend,[OnBeam_Label,OffBeam_Label,'Extra Beam 1','Extra Beam 2'],psym=[46,47,48,49],$
         textcolor=['Black','red','blue','green'],/right,/top,charsize=leg_charsize,symsize=2
     Endif
     ;Sky, ExBeam1
     If keyword_set(ExBeam1) and not(keyword_set(ExBeam2)) and not(keyword_set(TARGET_ONLY)) then begin
       al_legend,[OnBeam_Label,OffBeam_Label,'Extra Beam 1'],psym=[46,47,48],$
         textcolor=['Black','red','blue'],/right,/top,charsize=leg_charsize,symsize=2
     endif
     ;Sky, ExBeam2
     If keyword_set(ExBeam2) and not(keyword_set(ExBeam1)) and not(keyword_set(TARGET_ONLY)) then begin
       al_legend,[OnBeam_Label,OffBeam_Label,'Extra Beam 2'],psym=[46,47,49],$
         textcolor=['Black','red','green'],/right,/top,charsize=leg_charsize,symsize=2
     endif
   endif ;if no list
 Endif ;end single phase
      
 ;-------------------------------------------------------------------------------------------
 ;                                        Multi-plot
 ;-------------------------------------------------------------------------------------------
 If keyword_Set(multi_plot) then begin
 !p.multi=[0,2,2]
     
 ;**********************
 ;         Plot 1 
 ;**********************
 If not(keyword_set(LIST)) then begin 
   ;On Beam Plot
   cgplot,xt-min(xt),x,err_yhigh=error,err_ylow=error,title=label,xtitle='Time (Hour)',ytitle=ytitle,$
     charsize=char_size,/xstyle,/ystyle,psym=-46,yrange=[ymin,ymax]
      
   ;Off Beam
   If not(keyword_set(TARGET_ONLY)) then begin
     cgplot,xt-min(xt),x2,err_yhigh=error,err_ylow=error,/overplot,color='red',psym=-47
   Endif
     
   ;Plot ExBeam1
   If keyword_set(ExBeam1) then begin
      cgplot,ExBeam1[2,*]-min(ExBeam1[2,*]),ExBeam1[1,*],/overplot,color='blue',psym=-48 
   Endif
      
   ;Plot ExBeam2
   If keyword_set(ExBeam2) then begin
      cgplot,ExBeam2[2,*]-min(ExBeam2[2,*]),ExBeam2[1,*],/overplot,color='green',psym=-49
   Endif 
 Endif ;no list 
   
  ;normal: Target vs Sky
  If not(keyword_set(ExBeam1)) and not(keyword_set(ExBeam2)) and not(keyword_set(TARGET_ONLY)) then begin
      al_legend,[OnBeam_Label,OffBeam_Label],textcolor=['Black','red'],psym=[46,47],/right,/top,charsize=leg_charsize,symsize=2
  Endif
  
  ;Sky, ExBeam1, ExBeam2
  If keyword_set(ExBeam1) and keyword_set(ExBeam2) and not(keyword_set(TARGET_ONLY)) then begin
    al_legend,[OnBeam_Label,OffBeam_Label,'Extra Beam 1','Extra Beam 2'],psym=[46,47,48,49],$
      textcolor=['Black','red','blue','green'],/right,/top,charsize=leg_charsize,symsize=2
  Endif
  
  ;Sky, ExBeam1
  If keyword_set(ExBeam1) and not(keyword_set(ExBeam2)) and not(keyword_set(TARGET_ONLY)) then begin
    al_legend,[OnBeam_Label,OffBeam_Label,'Extra Beam 1'],psym=[46,47,48],$
      textcolor=['Black','red','blue'],/right,/top,charsize=leg_charsize,symsize=2
  endif
  
  ;Sky, ExBeam2
  If keyword_set(ExBeam2) and not(keyword_set(ExBeam1)) and not(keyword_set(TARGET_ONLY)) then begin
    al_legend,[OnBeam_Label,OffBeam_Label,'Extra Beam 2'],psym=[46,47,49],$
      textcolor=['Black','red','green'],/right,/top,charsize=leg_charsize,symsize=2
  endif
     
 ;**********************
 ;         Plot 2
 ;**********************
 If not(keyword_set(LIST)) then begin
   ; 12sec time intevrals
   If keyword_set(rebin_time) then win = long(12/rebin_time)
   If not(keyword_set(rebin_time)) then win = 2.
     x_plot2  = smooth(x,win,/EDGE_TRUNCATE)
     x2_plot2 =  smooth(x2,win,/EDGE_TRUNCATE) 
   
     cgplot,xt-min(xt),x_plot2,title=label,xtitle='Time (Hour)',$
       ytitle=ytitle,charsize=char_size,/xstyle,/ystyle,yrange=[ymin,ymax]
   
   If not(keyword_set(TARGET_ONLY)) then begin
     cgplot,xt-min(xt),x2_plot2,/overplot,color='red'
   endif 
   
   If keyword_set(ExBeam1) then begin
    ExBeam1_plot2 = smooth(reform(ExBeam1[1,*]),win,/EDGE_TRUNCATE)
     cgplot,ExBeam1[2,*]-min(ExBeam1[2,*]),ExBeam1_plot2,/overplot,color='blue';,psym=1
   Endif
   
   If keyword_set(ExBeam2) then begin
    ExBeam2_plot2 = smooth(reform(ExBeam2[1,*]),win,/EDGE_TRUNCATE)
     cgplot,ExBeam2[2,*]-min(ExBeam2[2,*]),ExBeam2_plot2,/overplot,color='green';,psym=1
   Endif
     
   ;--------------------------------
   ;               Legend 
   ;--------------------------------
   ;normal: Target vs Sky
   If not(keyword_set(ExBeam1)) and not(keyword_set(ExBeam2)) and not(keyword_set(TARGET_ONLY)) then begin
     al_legend,[OnBeam_Label,OffBeam_Label],textcolor=['Black','red'],psym=[46,47],/right,/top,charsize=leg_charsize,symsize=2
   Endif
   ;Sky, ExBeam1, ExBeam2
   If keyword_set(ExBeam1) and keyword_set(ExBeam2) and not(keyword_set(TARGET_ONLY)) then begin
     al_legend,[OnBeam_Label,OffBeam_Label,'Extra Beam 1','Extra Beam 2'],psym=[46,47,48,49],$
       textcolor=['Black','red','blue','green'],/right,/top,charsize=leg_charsize,symsize=2
   Endif
   ;Sky, ExBeam1
   If keyword_set(ExBeam1) and not(keyword_set(ExBeam2)) and not(keyword_set(TARGET_ONLY)) then begin
     al_legend,[OnBeam_Label,OffBeam_Label,'Extra Beam 1'],psym=[46,47,48],$
       textcolor=['Black','red','blue'],/right,/top,charsize=leg_charsize,symsize=2
   endif
   ;Sky, ExBeam2
   If keyword_set(ExBeam2) and not(keyword_set(ExBeam1)) and not(keyword_set(TARGET_ONLY)) then begin
     al_legend,[OnBeam_Label,OffBeam_Label,'Extra Beam 2'],psym=[46,47,49],$
       textcolor=['Black','red','green'],/right,/top,charsize=leg_charsize,symsize=2
   endif
 endif ;end no list 

 ;**********************
 ;         Plot 3
 ;**********************
 If not(keyword_set(LIST)) then begin
   ;target plot
   phase1 = (((t - T0) mod P)/(P*1.d))
   s1 = sort(phase1)
   cgplot,phase1[s1],x[s1],title=label,xtitle='Phase ',ytitle=ytitle,$
     charsize=char_size,psym=46,/xstyle,/ystyle,yrange=[ymin,ymax]
   
   If not(keyword_set(TARGET_ONLY)) then begin
    phase2 = (((t2 - T0) mod P)/(P*1.d)) 
    s2 = sort(phase2)
    cgplot,phase2[s2],x2[s2],/overplot,color='red',psym=46
   endif 
   
   If keyword_set(ExBeam1) then begin
     phase3       = (((ExBeam1[0,*] - T0) mod P)/(P*1.d)); mod 1
     s3 = sort(phase3)
     cgplot,phase3[s3],ExBeam1[1,s3],/overplot,color='blue',psym=46
   Endif
  
   If keyword_set(ExBeam2) then begin
    phase4 = (((ExBeam2[0,*] - T0) mod P)/(P*1.d)) ;mod 1
    s4 = sort(phase4) 
    cgplot,phase4[s4],ExBeam2[1,s4],/overplot,color='green',psym=46
   Endif
   
   ;--------------------------------
   ;               Legend
   ;--------------------------------
   ;normal: Target vs Sky
   If not(keyword_set(ExBeam1)) and not(keyword_set(ExBeam2)) and not(keyword_set(TARGET_ONLY)) then begin
     al_legend,[OnBeam_Label,OffBeam_Label],textcolor=['Black','red'],psym=[46,47],/right,/top,charsize=leg_charsize,symsize=2
   Endif
   ;Sky, ExBeam1, ExBeam2
   If keyword_set(ExBeam1) and keyword_set(ExBeam2) and not(keyword_set(TARGET_ONLY)) then begin
     al_legend,[OnBeam_Label,OffBeam_Label,'Extra Beam 1','Extra Beam 2'],psym=[46,47,48,49],$
       textcolor=['Black','red','blue','green'],/right,/top,charsize=leg_charsize,symsize=2
   Endif
   ;Sky, ExBeam1
   If keyword_set(ExBeam1) and not(keyword_set(ExBeam2)) and not(keyword_set(TARGET_ONLY)) then begin
     al_legend,[OnBeam_Label,OffBeam_Label,'Extra Beam 1'],psym=[46,47,48],$
       textcolor=['Black','red','blue'],/right,/top,charsize=leg_charsize,symsize=2
   endif
   ;Sky, ExBeam2
   If keyword_set(ExBeam2) and not(keyword_set(ExBeam1)) and not(keyword_set(TARGET_ONLY)) then begin
     al_legend,[OnBeam_Label,OffBeam_Label,'Extra Beam 2'],psym=[46,47,49],$
       textcolor=['Black','red','green'],/right,/top,charsize=leg_charsize,symsize=2
   endif
  endif ;if no list
 
 ;**********************
 ;         Plot 4
 ;**********************
 If keyword_set(rebin_time) then win = long(12/rebin_time) ;12 time intervals
 If not(keyword_set(rebin_time)) then win = 2.

 If not(keyword_set(LIST)) then begin
     x_plot4 = smooth(x,win,/EDGE_WRAP)
     x2_plot4 = smooth(x2,win,/EDGE_WRAP)

   cgplot,phase1,x_plot4,title=label,xtitle='Phase',$
       ytitle=ytitle,charsize=char_size,/xstyle,/ystyle,yrange=[ymin,ymax],psym=46
  
   If not(keyword_set(TARGET_ONLY)) then begin
     cgplot,phase1,x2_plot4,/overplot,color='red',psym=47
   endif 
   
   If keyword_set(ExBeam1) then begin
     ExBeam1_plot4 = smooth(reform(ExBeam1[1,*]),win,/EDGE_WRAP)
     cgplot,phase3,ExBeam1_plot4,/overplot,color='blue',psym=48;,psym=1
   Endif
  
   If keyword_set(ExBeam2) then begin
     ExBeam2_plot4 = smooth(reform(ExBeam2[1,*]),win,/EDGE_WRAP)
     cgplot,phase4,ExBeam2_plot4,/overplot,color='green',psym=49;,psym=1
   Endif
   
   ;--------------------------------
   ;               Legend
   ;--------------------------------
   ;normal: Target vs Sky
   If not(keyword_set(ExBeam1)) and not(keyword_set(ExBeam2)) and not(keyword_set(TARGET_ONLY)) then begin
     al_legend,[OnBeam_Label,OffBeam_Label],textcolor=['Black','red'],psym=[46,47],/right,/top,charsize=leg_charsize,symsize=2
   Endif
   ;Sky, ExBeam1, ExBeam2
   If keyword_set(ExBeam1) and keyword_set(ExBeam2) and not(keyword_set(TARGET_ONLY)) then begin
     al_legend,[OnBeam_Label,OffBeam_Label,'Extra Beam 1','Extra Beam 2'],psym=[46,47,48,49],$
       textcolor=['Black','red','blue','green'],/right,/top,charsize=leg_charsize,symsize=2
   Endif
   ;Sky, ExBeam1
   If keyword_set(ExBeam1) and not(keyword_set(ExBeam2)) and not(keyword_set(TARGET_ONLY)) then begin
     al_legend,[OnBeam_Label,OffBeam_Label,'Extra Beam 1'],psym=[46,47,48],$
       textcolor=['Black','red','blue'],/right,/top,charsize=leg_charsize,symsize=2
   endif
   ;Sky, ExBeam2
   If keyword_set(ExBeam2) and not(keyword_set(ExBeam1)) and not(keyword_set(TARGET_ONLY)) then begin
     al_legend,[OnBeam_Label,OffBeam_Label,'Extra Beam 2'],psym=[46,47,49],$
       textcolor=['Black','red','green'],/right,/top,charsize=leg_charsize,symsize=2
   endif
 endif ;end no LIST  
 endif ; multiplot
return 
end
