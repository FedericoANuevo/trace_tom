; This script plots, for each field-line of CR-2099, the Ne(r) and Tm(r)
; profiles and their fits for AIA, Mk4, and LASCO-C2, respectively. 


pro call_load_structure
  
  common data, N_fl, Npt_max, Npt_v, x_A, y_A, z_A, rad_A, lat_A, lon_A,$
     Ne_aia_A, Tm_aia_A, WT_aia_A, ldem_flag_aia_A, index_aia_A, index_sampling_aia_A,$
     Ne_euvia_A, Tm_euvia_A,  WT_euvia_A, ldem_flag_euvia_A, index_euvia_A, index_sampling_euvia_A,$
     Ne_euvib_A, Tm_euvib_A, WT_euvib_A, ldem_flag_euvib_A, index_euvib_A, index_sampling_euvib_A,$
     Ne_eit_A, Tm_eit_A, WT_eit_A, ldem_flag_eit_A, index_eit_A, index_sampling_eit_A,$
     Ne_mk4_A, index_mk4_A, index_sampling_mk4_A,$
     Ne_kcor_A, index_kcor_A, index_sampling_kcor_A,$
     Ne_c2_A, index_c2_A, index_sampling_c2_A,$
     rad_fit_aia_A, Ne_fit_aia_A, Tm_fit_aia_A, fitflag_aia_A,scN_fit_aia_A,scT_fit_aia_A,$
     rad_fit_euvia_A, Ne_fit_euvia_A, Tm_fit_euvia_A, fitflag_euvia_A,scN_fit_euvia_A,scT_fit_euvia_A,$
     rad_fit_euvib_A, Ne_fit_euvib_A, Tm_fit_euvib_A, fitflag_euvib_A,scN_fit_euvib_A,scT_fit_euvib_A,$
     rad_fit_eit_A, Ne_fit_eit_A, Tm_fit_eit_A, fitflag_eit_A,scN_fit_eit_A,scT_fit_eit_A,$
     rad_fit_c2_A, Ne_fit_c2_A, fitflag_c2_A,scN_fit_c2_A,$
     rad_fit_mk4_A, Ne_fit_mk4_A, fitflag_mk4_A,scN_fit_mk4_A,$
     N0_fit_aia_A,lN_fit_aia_A,T0_fit_aia_A,dTdr_fit_aia_A,$
     N0_fit_euvia_A,lN_fit_euvia_A,T0_fit_euvia_A,dTdr_fit_euvia_A,$
     N0_fit_euvib_A,lN_fit_euvib_A,T0_fit_euvib_A,dTdr_fit_euvib_A,$
     N0_fit_eit_A,lN_fit_eit_A,T0_fit_eit_A,dTdr_fit_eit_A,$
     N1_fit_aia_A,N2_fit_aia_A,p1_fit_aia_A,p2_fit_aia_A,$
     N1_fit_euvia_A,N2_fit_euvia_A,p1_fit_euvia_A,p2_fit_euvia_A,$
     N1_fit_euvib_A,N2_fit_euvib_A,p1_fit_euvib_A,p2_fit_euvib_A,$
     N1_fit_eit_A,N2_fit_eit_A,p1_fit_eit_A,p2_fit_eit_A,$
     N0_fit_mk4_A,lN_fit_mk4_A,$
     N1_fit_mk4_A,N2_fit_mk4_A,p1_fit_mk4_A,p2_fit_mk4_A,$
     N1_fit_c2_A,N2_fit_c2_A,p1_fit_c2_A,p2_fit_c2_A,$
     lN_fit_c2_A,$
     fit_F_Ne_aia,fit_F_Ne_mk4,fit_F_Ne_c2,$
     fit_F_Ne_euvia,fit_F_Ne_euvib,fit_F_eit_c2

; Directory, File-name and Window
; dir  = '/data1/DATA/fieldlines_judit/CR2099/map1/'
; file = 'list.map1.txt-tracing-structure-merge_aia_mk4_lascoc2.sav' & win=0
  dir  = '/data1/DATA/fieldlines_judit/CR2099/map7/'
  file = 'list.map7.txt-tracing-structure-merge_aia_mk4_lascoc2.sav' & win=0

; Load the traced-data-structure (stored in common data )
  load_traced_data_structure,dir=dir,structure_filename=file,/aia,/mk4,/lascoc2

; visualize each field-line profiles for AIA data 
; goto,skip_AIA
  print,'AIA profiles'
  print,'Basic Statistics:'
  i_fit = where(fitflag_aia_A eq +1)
  print,'field-lines with fit:',n_elements(i_fit), ' of ',N_fl
  i_goodfit = where( scN_fit_aia_A gt 0. and  scN_fit_aia_A lt 0.2)
  print,'field-lines with GOOD fit in Ne:',n_elements(i_goodfit), ' of ',N_fl
  i_goodfit = where( scT_fit_aia_A gt 0. and  scT_fit_aia_A lt 0.2)
  print,'field-lines with GOOD fit in Tm:',n_elements(i_goodfit), ' of ',N_fl
  print
  !p.multi=[0,1,2]
  for ifl = 0,N_fl-1 do begin
     print,'field-line #',ifl,'/',N_fl-1
     print,'init coord.',rad_A(ifl,0),lat_A(ifl,0),lon_A(ifl,0)
     print,'final coord.',rad_A(ifl,Npt_v(ifl)-1),lat_A(ifl,Npt_v(ifl)-1),lon_A(ifl,Npt_v(ifl)-1)
;    Plot geometry of the field-line
     window,win+2,xs=700,ys=1400
     plot,rad_A(ifl,0:Npt_v(ifl)-1),lat_A(ifl,0:Npt_v(ifl)-1),$
          xtitle='rad [Rsun]',ytitle='lat [DEG]',charsize=2
     plot,rad_A(ifl,0:Npt_v(ifl)-1),lat_A(ifl,0:Npt_v(ifl)-1),$
          xr=[1.02,1.25],xstyle=1,xtitle='rad [Rsun]',ytitle='lat [DEG]',charsize=2
;    Ne(r) and Tm(r) profiles  and their fits
     tmp = reform(index_sampling_aia_A(ifl,*))
     ind_samp_aia = where(tmp eq 1)
     if ind_samp_aia[0] ne -1 then begin
        radsamp = reform(rad_A(ifl,ind_samp_aia)) ; Rsun
        Nesamp  = reform(Ne_aia_A(ifl,ind_samp_aia))
        Tmsamp  = reform(Tm_aia_A(ifl,ind_samp_aia))
        WTsamp  = reform(WT_aia_A(ifl,ind_samp_aia))
        LFsamp  = reform(LDEM_FLAG_aia_A(ifl,ind_samp_aia))
        i1      = where( LFsamp eq 0.)
        i2      = where( LFsamp ne 0.)
        if fitflag_aia_A(ifl) eq  +1 then begin
           radfit = rad_fit_aia_A
           Nefit  = Ne_fit_aia_A(ifl,*)
           Tmfit  = Tm_fit_aia_A(ifl,*)
        endif
;       Plot Ne(r) profile and overplot fit      
        !p.multi=[0,1,2]
        window,win,xs=700,ys=1400
        plot,radsamp,Nesamp/1.e8,psym=4,xr=[1.02,1.25],xstyle=1,$
             xtitle='r [R!dSUN!N]', ytitle= 'N!de!n [10!u8!n cm!u-3!n]',$
             charsize=2.,yr=[0.,4.],ystyle=1,title='AIA',/nodata
        if i1[0] ne -1 then oplot,radsamp(i1),Nesamp(i1)/1.e8,psym=4,symsize=2
        if i2[0] ne -1 then oplot,radsamp(i2),Nesamp(i2)/1.e8,psym=5,symsize=2
        if fitflag_aia_A(ifl) eq  +1  then $
           oplot,radfit,Nefit/1.e8
;       Plot Tm(r) profile and overplot fit      
        plot,radsamp,Tmsamp/1.e6,psym=4,xr=[1.02,1.25],xstyle=1,$
             xtitle='r [R!dSUN!N]', ytitle= 'T!dm!n [MK]',$
             charsize=2.,yr=[0.5,4.],ystyle=1,/nodata
        if i1[0] ne -1 then oplot,radsamp(i1),Tmsamp(i1)/1.e6,psym=4,symsize=2
        if i2[0] ne -1 then oplot,radsamp(i2),Tmsamp(i2)/1.e6,psym=5,symsize=2
        errplot,radsamp,Tmsamp/1.e6-WTsamp/1.e6,Tmsamp/1.e6+WTsamp/1.e6
        if fitflag_aia_A(ifl) eq  +1  then $
           oplot,radfit,Tmfit/1.e6
;       Print quality fit indicator and density scale height  
        if fitflag_aia_A(ifl) eq +1 then begin
           print,'Ne fit quality:', scN_fit_aia_A(ifl)
           print,'Tm fit quality:', scT_fit_aia_A(ifl)
           print,'density scale height',lN_fit_aia_A(ifl)
        endif else begin
           PRINT,'NO FIT'
        endelse
     endif else begin
        print,'NO DATA'
     endelse
     STOP
  endfor
  !p.multi=0
  skip_AIA:

; Visualize each field-line profiles for Mk4 data 
  goto,skip_MK4
  print,'Mk4 profiles'
  print,'Basic Statistics:'
  i_fit = where(fitflag_mk4_A eq +1)
  print,'field-lines with fit:',n_elements(i_fit), ' of ',N_fl
  i_goodfit = where( scN_fit_mk4_A gt 0. and  scN_fit_mk4_A lt 0.2)
  print,'field-lines with GOOD fit in Ne:',n_elements(i_goodfit), ' of ',N_fl
  print
  for ifl = 0,N_fl-1 do begin
     print,'field-line #',ifl,'/',N_fl-1
     print,'init coord.',rad_A(ifl,0),lat_A(ifl,0),lon_A(ifl,0)
     print,'final coord.',rad_A(ifl,Npt_v(ifl)-1),lat_A(ifl,Npt_v(ifl)-1),lon_A(ifl,Npt_v(ifl)-1)
;    Plot field-line geometry
     window,win+2
     plot,rad_A(ifl,0:Npt_v(ifl)-1),lat_A(ifl,0:Npt_v(ifl)-1),$
          xr=[1.15,1.5],xstyle=1,xtitle='rad [Rsun]',ytitle='lat [DEG]',charsize=2
;    Ne(r) profile  and its fit
     tmp = reform(index_sampling_mk4_A(ifl,*))
     ind_samp_mk4 = where(tmp eq 1)
     if ind_samp_mk4[0] ne -1 then begin
        radsamp = reform(rad_A(ifl,ind_samp_mk4))      ; Rsun
        Nesamp  = reform(Ne_mk4_A(ifl,ind_samp_mk4))   ; cm-3
        if fitflag_mk4_A(ifl) eq  +1 then begin
           radfit = rad_fit_mk4_A
           Nefit  = Ne_fit_mk4_A(ifl,*)
        endif
;       Plot Ne(r) profile and overplot fit          
        window,win
        plot,radsamp,Nesamp/1.e8,psym=4,$
             xtitle='r [R!dSUN!N]', ytitle= 'N!de!n [10!u8!n cm!u-3!n]',$
             charsize=2.,xr=[1.15,1.5],xstyle=1,title='Mk4'
        if fitflag_mk4_A(ifl) eq  +1  then $
           oplot,radfit,Nefit/1.e8
;       Print quality fit indicator and density scale height  
        if fitflag_mk4_A(ifl) eq +1 then begin
           print,'Ne fit quality:',scN_fit_mk4_A(ifl)
           print,'density scale height',lN_fit_mk4_A(ifl)
        endif else begin
           PRINT,'NO FIT'
        endelse
     endif else begin
        print,'NO DATA'
     endelse
     STOP
  endfor
  skip_MK4:
    
; Visualize each field-line profiles for LASCO-C2 data 
    goto,skip_C2
    print,'LASCO-C2 profiles'
    print,'Basic Statistics:'
    i_fit = where(fitflag_c2_A eq +1)
    print,'field-lines with fit:',n_elements(i_fit), ' of ',N_fl
    i_goodfit = where( scN_fit_c2_A gt 0. and  scN_fit_c2_A lt 0.2)
    print,'field-lines with GOOD fit in Ne:',n_elements(i_goodfit), ' of ',N_fl
    print
    for ifl = 0,N_fl-1 do begin
       print,'field-line #',ifl,'/',N_fl-1
       print,'init coord.',rad_A(ifl,0),lat_A(ifl,0),lon_A(ifl,0)
       print,'final coord.',rad_A(ifl,Npt_v(ifl)-1),lat_A(ifl,Npt_v(ifl)-1),lon_A(ifl,Npt_v(ifl)-1)
;      Plot field-line geometry
       window,win+2
       plot,rad_A(ifl,0:Npt_v(ifl)-1),lat_A(ifl,0:Npt_v(ifl)-1),$
            xr=[2.5,6.0],xstyle=1,xtitle='rad [Rsun]',ytitle='lat [DEG]',charsize=2
;      Ne(r) profile and its fit      
       tmp = reform(index_sampling_c2_A(ifl,*))
       ind_samp_c2 = where(tmp eq 1)
       if ind_samp_c2[0] ne -1 then begin
          radsamp = reform(rad_A(ifl,ind_samp_c2)) ; Rsun
          Nesamp  = reform(Ne_c2_A(ifl,ind_samp_c2))
          if fitflag_c2_A(ifl) eq  +1 then begin
             radfit = rad_fit_c2_A
             Nefit  = Ne_fit_c2_A(ifl,*)
          endif
;         Plot Ne(r) profile and overplot fit          
          window,win           
          plot,radsamp,Nesamp/1.e5,psym=4,$
               xtitle='r [R!dSUN!N]', ytitle= 'N!de!n [10!u5!n cm!u-3!n]',$
               charsize=2.,xr=[2.5,6.0],xstyle=1,title='LASCO-C2'
          if fitflag_c2_A(ifl) eq  +1  then $
             oplot,radfit,Nefit/1.e5
;         Print quality fit indicator and density scale height          
          IF fitflag_c2_A(ifl) eq +1 then begin
             print,'Ne fit quality:',scN_fit_c2_A(ifl)
             print,'density scale height',lN_fit_c2_A(ifl)
          endif else begin
             PRINT,'NO FIT'
          endelse
       endif else begin
          print,'NO DATA'
       ENDELSE
       STOP
    endfor
    skip_C2:


    ;STOP
    return
 end
