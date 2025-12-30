;;
; Calling sequence examples:
;
; line_analysis, cr_number = 2082, map_number = 1, /euvia, /lascoc2
; line_analysis, cr_number = 2082, map_number = 7, /euvia, /lascoc2
; line_analysis, cr_number = 2099, map_number = 1, /aia  , /lascoc2, /mk4
; line_analysis, cr_number = 2099, map_number = 7, /aia  , /lascoc2, /mk4
;;

pro line_analysis, rel_sqrt_chisqr_crit=rel_sqrt_chisqr_crit,$
                   map_number=map_number, cr_number=cr_number, $
                   aia = aia, euvia = euvia, euvib = euvib, eit = eit,$
                   mk4 = mk4, kcor = kcor, lascoc2 = lascoc2,$
                   npower = npower, dpower = dpower, linfit = linfit

; 0) Set a default value for rel_sqrt_chisqr_crit
  if not keyword_set(rel_sqrt_chisqr_crit) then rel_sqrt_chisqr_crit = 0.2
  
; 1) Set filename of the structure corresponding to selected CR_NUM and MAP_NUM.
;    Also define groups of field lines according to terminal longitude.
  dir = './'
  if cr_number eq 2082 then begin
     if map_number eq 1 then begin
       ;structure_filename = 'CR2082_AWSoM-map1-tracing-structure-merge_euvia_lascoc2.sav'
        dir                = '/data1/DATA/trace_tom_files/CR2082/field_lines_geometry_map1/'
        structure_filename = 'list.map1.new.txt-tracing-structure-merge_euvia_lascoc2_sampled.sav'
        CritTermLon = [0.,100.,180.,270.,360.]
     endif  
     if map_number eq 7 then begin
       ;structure_filename = 'CR2082_AWSoM-map7-tracing-structure-merge_euvia_lascoc2.sav'
        dir                = '/data1/DATA/trace_tom_files/CR2082/field_lines_geometry_map7/'
        structure_filename = 'list.map7.new.txt-tracing-structure-merge_euvia_lascoc2_sampled.sav'
        CritTermLon = [0.,100.,180.,270.,360.]
     endif
  endif
  if cr_number eq 2099 then begin
     if map_number eq 1 then begin
       ;structure_filename = 'CR2099_AWSoM-map1_tracing-structure-merge_aia_mk4_lascoc2.sav'
        dir                = '/data1/DATA/trace_tom_files/CR2099/field_lines_geometry_map1/'
        structure_filename = 'list.map1.new.txt-tracing-structure-merge_aia_mk4_lascoc2_sampled.sav'
        CritTermLon = [0.,100.,180.,270.,360.]
     endif  
     if map_number eq 7 then begin
       ;structure_filename = 'CR2099_AWSoM-map7_tracing-structure-merge_aia_mk4_lascoc2.sav'
        dir                = '/data1/DATA/trace_tom_files/CR2099/field_lines_geometry_map7/'
        structure_filename = 'list.map7.new.txt-tracing-structure-merge_aia_mk4_lascoc2_sampled.sav'
        CritTermLon = [0.,100.,180.,270.,310.,360.]
     endif
  endif
  if NOT keyword_set(structure_filename) then begin
     print, 'Invalid CR and/or MAP number, try again.'
     return
  endif

; Restore the structure with the traced-data  
  restore,dir+structure_filename

; Number of field-lines  
  N_fl = *trace_data.N_FL
;;
; Determine the Lat and Lon of the footppoint and terminalpoint of each field line.
; Group field lines according to the Carrington Longitude of their terminalpoint.
;
; 1D Arrays: Footpoint and Terminalpoint Lon and Lat for each field line:
  Footpoint_Lon = *trace_data.FOOTPOINT_Lon
  Footpoint_Lat = *trace_data.FOOTPOINT_Lat
  Termpoint_Lon = *trace_data.TERMPOINT_Lon
  Termpoint_Lat = *trace_data.TERMPOINT_Lat

; Tag groups of field lines by means of user-defined ranges of terminal Longitudes.
  Ngroups     = n_elements(CritTermLon)-1
; Tag each field line (ifl) with the group number (ig) to which it belongs: 
  line_groupID  = intarr(N_FL)
  for ig=0,Ngroups-1 do begin
     ifl_A = where(Termpoint_Lon gt CritTermLon(ig) AND Termpoint_Lon lt CritTermLon(ig+1))
     line_groupID(ifl_A) = ig
  endfor
;;

; Lat/Lon plots of FootPoint and TermPoint
  ps1,'./'+structure_filename+'_connectivity-map.eps'
  np=1000
  !p.multi=[0,1,2]
  loadct,0
  !p.color=0
  !p.background=255
  csz=1
  plot,footpoint_lon,footpoint_lat,xr=[0,360],yr=[-90,+90],xstyle=1,ystyle=1,/nodata,charsize=csz,title=strmid(structure_filename,0,17)+'  r = 1.0 Rs',ytitle='Carrington Latitude [deg]'
  ctbl = 12                     ; 16-LEVEL
  colors = 16 + 190 * (indgen(Ngroups))/float(Ngroups-1)
  loadct,ctbl
  for ifl=0,N_FL-1 do oplot,[Footpoint_Lon(ifl)],[Footpoint_Lat(ifl)],psym=4,color=colors(line_groupID(ifl)),th=2
  loadct,0
  plot,footpoint_lon,footpoint_lat,xr=[0,360],yr=[-90,+90],xstyle=1,ystyle=1,/nodata,charsize=csz,xtitle='Carrington Longitude [deg]',title='r = 23.68 Rs',ytitle='Carrington Latitude [deg]'
  loadct,ctbl
  for ifl=0,N_FL-1 do oplot,[Termpoint_Lon(ifl)],[Termpoint_Lat(ifl)],psym=4,color=colors(line_groupID(ifl)),th=2
  loadct,0
  !p.multi=0
  ps2


;;
; Compute and plot average trends <Ne(r)> and <Te(r)> for available
; instruments and for each group of field lines:

; Set up graphical stuff
  nx   = ngroups                ; number of horizontal panels
  ny   = 2                      ; number of vertical   panels
  csz  = 0.75                   ; charsize for plots
  ctbl = 40                     ; color table to use for individual field lines

  if keyword_set(aia) then begin
; -------------------------------------------
     fitflag_AIA_A = *trace_data.fitflag_AIA
     rad_fit_aia_A = *trace_data.rad_fit_aia
      Ne_fit_aia_A = *trace_data. Ne_fit_aia
     scN_fit_aia_A = *trace_data.scN_fit_aia
      Tm_fit_aia_A = *trace_data. Tm_fit_aia
     scT_fit_aia_A = *trace_data.scT_fit_aia
; -------------------------------------------
   ; AIA Ne 
     Ne_fit_aia_groupavg = fltarr(Ngroups,n_elements(rad_fit_aia_A))
;---- Tag field lines for which the fit is positive at all heights ----
     tag_pos = fltarr(N_fl)
     for ifl=0,N_fl-1 do begin
        if min(Ne_fit_aia_A(ifl,*)) gt 0. then tag_pos(ifl)=+1
     endfor
;----------------------------------------------------------------------
     scN_crit = rel_sqrt_chisqr_crit
     !p.multi=[0,nx,ny]
     ps1,'./'+structure_filename+'_AIA-DEMT_Ne-profiles.eps'
     for ig=0,Ngroups-1 do begin
        ifl = where(line_groupID eq ig AND fitflag_AIA_A eq +1 AND tag_pos eq +1 AND scN_fit_aia_A le scN_crit)
   if ifl(0) eq -1 then begin
      plot,rad_fit_aia_A,Ne_fit_aia_groupavg(ig,*),charsize=csz,xtitle='r [Rsun]',title='AIA-DEMT Ne(r) [cm!U-3!N]. Group #'+strmid(string(ig+1),7,1),th=4, ystyle=2, xr=[1.02,1.25], xstyle=1,/nodata
      goto,skip_group_aia_Ne
   endif
   Ne_fit_aia_groupavg(ig,*) = total( Ne_fit_aia_A(ifl,*) , 1 ) / float(n_elements(ifl))
   loadct,0
   plot,rad_fit_aia_A,Ne_fit_aia_groupavg(ig,*),charsize=csz,xtitle='r [Rsun]',title='AIA-DEMT Ne(r) [cm!U-3!N]. Group #'+strmid(string(ig+1),7,1),th=4, ystyle=2, xr=[1.02,1.25], xstyle=1,/nodata
   loadct,ctbl
   color_index_step = fix(256./n_elements(ifl))
   for index=0,n_elements(ifl)-1 do begin
 ; ---------------------------------     
      Nsamp   = (*trace_data.Npt_aia)(ifl(index))
      rad_aia = (*trace_data.rad_aia)(ifl(index),0:Nsamp-1)
       Ne_aia = (*trace_data. Ne_aia)(ifl(index),0:Nsamp-1)
      oplot,rad_fit_aia_A                 ,Ne_fit_aia_A(ifl(index),*)       ,color=(index)*color_index_step
      oplot,rad_aia,Ne_aia,color=(index)*color_index_step,psym=4
 ; ---------------------------------     
   endfor
   loadct,0
   oplot,rad_fit_aia_A,Ne_fit_aia_groupavg(ig,*),th=4
   skip_group_aia_Ne:
endfor
  
for ig=0,Ngroups-1 do begin
   ifl = where(line_groupID eq ig AND fitflag_AIA_A eq +1 AND tag_pos eq +1 AND scN_fit_AIA_A le scN_crit)
   if ifl(0) eq -1 then begin
      plot,CritTermLon(ig:ig+1),[1.02,1.25],charsize=csz,xtitle='Lon [deg]',title='Rad [Rsun]. Group #'+strmid(string(ig+1),7,1),th=4, ystyle=1,xstyle=1,/nodata
      goto,skip_group_aia_geo
   endif
   loadct,0
   plot,CritTermLon(ig:ig+1),[1.02,1.25],charsize=csz,xtitle='Lon [deg]',title='Rad [Rsun]. Group #'+strmid(string(ig+1),7,1),th=4, ystyle=1,xstyle=1,/nodata
   loadct,ctbl
   color_index_step = fix(256./n_elements(ifl))
   for index=0,n_elements(ifl)-1 do begin
  ; ---------------------------------     
      Nsamp =   (*trace_data.Npt_aia)(ifl(index))
      rad_aia = (*trace_data.rad_aia)(ifl(index),0:Nsamp-1)
      lon_aia = (*trace_data.lon_aia)(ifl(index),0:Nsamp-1)
      oplot,lon_aia,rad_aia,color=(index)*color_index_step
   ; ---------------------------------        
   endfor
   loadct,0
   skip_group_aia_geo:
endfor
ps2


; AIA Te
Tm_fit_aia_groupavg = fltarr(Ngroups,n_elements(rad_fit_aia_A))
;---- Tag field lines for which the fit is positive at all heights ----
   tag_pos = fltarr(N_fl)
   for ifl=0,N_fl-1 do begin
      if min(Tm_fit_aia_A(ifl,*)) gt 0. then tag_pos(ifl)=+1
   endfor
;---- Tag field lines for which the fit is positive at all heights ----
   tag_pos = fltarr(N_fl)
   for ifl=0,N_fl-1 do begin
      if min(Ne_fit_aia_A(ifl,*)) gt 0. then tag_pos(ifl)=+1
   endfor
;----------------------------------------------------------------------
scT_crit = rel_sqrt_chisqr_crit
!p.multi=[0,nx,ny]
ps1,'./'+structure_filename+'_AIA-DEMT_Te-profiles.eps'
for ig=0,Ngroups-1 do begin
   ifl = where(line_groupID eq ig AND fitflag_AIA_A eq +1 AND tag_pos eq +1 AND scN_fit_aia_A le scT_crit)
   if ifl(0) eq -1 then begin
      plot,CritTermLon(ig:ig+1),[1.02,1.25],charsize=csz,xtitle='Lon [deg]',title='Rad [Rsun]. Group #'+strmid(string(ig+1),7,1),th=4, ystyle=1,xstyle=1,/nodata
      goto,skip_group_aia_Te
   endif
   Tm_fit_aia_groupavg(ig,*) = total( Tm_fit_aia_A(ifl,*) , 1 ) / float(n_elements(ifl))
   loadct,0
   plot,rad_fit_aia_A,Tm_fit_aia_groupavg(ig,*),charsize=csz,xtitle='r [Rsun]',title='AIA-DEMT Te(r) [MK]. Group #'+strmid(string(ig+1),7,1),th=4, yr=[0.5e6,3.0e6], ystyle=1, xr=[1.02,1.25], xstyle=1,/nodata
   loadct,ctbl
   color_index_step = fix(256./n_elements(ifl))
   for index=0,n_elements(ifl)-1 do begin
      Nsamp   = (*trace_data.Npt_aia)(ifl(index))
      rad_aia = (*trace_data.rad_aia)(ifl(index),0:Nsamp-1)
      Tm_aia  = (*trace_data.Tm_aia) (ifl(index),0:Nsamp-1)
      oplot,rad_fit_aia_A                 ,Tm_fit_aia_A(ifl(index),*)       ,color=(index)*color_index_step
      oplot,rad_aia,Tm_aia,color=(index)*color_index_step,psym=4
   endfor
   loadct,0
   oplot,rad_fit_aia_A,Tm_fit_aia_groupavg(ig,*),th=4
   skip_group_aia_Te:
endfor
for ig=0,Ngroups-1 do begin
   ifl = where(line_groupID eq ig AND fitflag_AIA_A eq +1 AND tag_pos eq +1 AND scN_fit_AIA_A le scN_crit)
   if ifl(0) eq -1 then begin
      plot,CritTermLon(ig:ig+1),[1.02,1.25],charsize=csz,xtitle='Lon [deg]',title='Rad [Rsun]. Group #'+strmid(string(ig+1),7,1),th=4, ystyle=1,xstyle=1,/nodata
      goto,skip_group_aia_geoT
   endif
   loadct,0
   plot,CritTermLon(ig:ig+1),[1.02,1.25],charsize=csz,xtitle='Lon [deg]',title='Rad [Rsun]. Group #'+strmid(string(ig+1),7,1),th=4, ystyle=1,xstyle=1,/nodata
   loadct,ctbl
   color_index_step = fix(256./n_elements(ifl))
   for index=0,n_elements(ifl)-1 do begin
      Nsamp = (*trace_data.Npt_aia)(ifl(index))
      rad_aia = (*trace_data.rad_aia)(ifl(index),0:Nsamp-1)
      lon_aia = (*trace_data.lon_aia)(ifl(index),0:Nsamp-1)
      oplot,lon_aia,rad_aia,color=(index)*color_index_step
   endfor
   loadct,0
   skip_group_aia_geoT:
endfor
ps2
endif
 

if keyword_set(euvia) then begin
; EUVI-A Ne
Ne_fit_euvia_groupavg = fltarr(Ngroups,n_elements(rad_fit_euvia_A))
;---- Tag field lines for which the fit is positive at all heights ----
   tag_pos = fltarr(N_fl)
   for ifl=0,N_fl-1 do begin
      if min(Ne_fit_euvia_A(ifl,*)) gt 0. then tag_pos(ifl)=+1
   endfor
;----------------------------------------------------------------------
scN_crit = rel_sqrt_chisqr_crit
!p.multi=[0,nx,ny]
ps1,'./'+structure_filename+'_EUVIA-DEMT_Ne-profiles.eps'
for ig=0,Ngroups-1 do begin
   ifl = where(line_groupID eq ig AND fitflag_EUVIA_A eq +1 AND tag_pos eq +1 AND scN_fit_euvia_A le scN_crit)
   if ifl(0) eq -1 then begin
      plot,rad_fit_euvia_A,Ne_fit_euvia_groupavg(ig,*),charsize=csz,xtitle='r [Rsun]',title='EUVIA-DEMT Ne(r) [cm!U-3!N]. Group #'+strmid(string(ig+1),7,1),th=4, ystyle=2, xr=[1.02,1.25], xstyle=1,/nodata
      goto,skip_group_euvia_Ne
   endif
   Ne_fit_euvia_groupavg(ig,*) = total( Ne_fit_euvia_A(ifl,*) , 1 ) / float(n_elements(ifl))
   loadct,0
   plot,rad_fit_euvia_A,Ne_fit_euvia_groupavg(ig,*),charsize=csz,xtitle='r [Rsun]',title='EUVIA-DEMT Ne(r) [cm!U-3!N]. Group #'+strmid(string(ig+1),7,1),th=4, ystyle=2, xr=[1.02,1.25], xstyle=1,/nodata
   loadct,ctbl
   color_index_step = fix(256./n_elements(ifl))
   for index=0,n_elements(ifl)-1 do begin
      tmp = reform(index_sampling_euvia_A(ifl(index),*))
      ind_samp_euvia = where(tmp eq 1)
      oplot,rad_fit_euvia_A                 ,Ne_fit_euvia_A(ifl(index),*)       ,color=(index)*color_index_step
      oplot,rad_A(ifl(index),ind_samp_euvia),Ne_euvia_A(ifl(index),ind_samp_euvia),color=(index)*color_index_step,psym=4
   endfor
   loadct,0
   oplot,rad_fit_euvia_A,Ne_fit_euvia_groupavg(ig,*),th=4
   skip_group_euvia_Ne:
endfor
for ig=0,Ngroups-1 do begin
   ifl = where(line_groupID eq ig AND fitflag_EUVIA_A eq +1 AND tag_pos eq +1 AND scN_fit_EUVIA_A le scN_crit)
   if ifl(0) eq -1 then begin
      plot,CritTermLon(ig:ig+1),[1.02,1.25],charsize=csz,xtitle='Lon [deg]',title='Rad [Rsun]. Group #'+strmid(string(ig+1),7,1),th=4, ystyle=1,xstyle=1,/nodata
      goto,skip_group_euvia_geo
   endif
   loadct,0
   plot,CritTermLon(ig:ig+1),[1.02,1.25],charsize=csz,xtitle='Lon [deg]',title='Rad [Rsun]. Group #'+strmid(string(ig+1),7,1),th=4, ystyle=1,xstyle=1,/nodata
   loadct,ctbl
   color_index_step = fix(256./n_elements(ifl))
   for index=0,n_elements(ifl)-1 do begin
      tmp = reform(index_sampling_euvia_A(ifl(index),*))
      ind_samp_euvia = where(tmp eq 1)
      oplot,lon_A(ifl(index),ind_samp_euvia),rad_A(ifl(index),ind_samp_euvia),color=(index)*color_index_step
   endfor
   loadct,0
   skip_group_euvia_geo:
endfor
ps2
; EUVI-A Te
Tm_fit_euvia_groupavg = fltarr(Ngroups,n_elements(rad_fit_euvia_A))
;---- Tag field lines for which the fit is positive at all heights ----
   tag_pos = fltarr(N_fl)
   for ifl=0,N_fl-1 do begin
      if min(Tm_fit_euvia_A(ifl,*)) gt 0. then tag_pos(ifl)=+1
   endfor
;---- Tag field lines for which the fit is positive at all heights ----
   tag_pos = fltarr(N_fl)
   for ifl=0,N_fl-1 do begin
      if min(Ne_fit_euvia_A(ifl,*)) gt 0. then tag_pos(ifl)=+1
   endfor
;----------------------------------------------------------------------
scT_crit = rel_sqrt_chisqr_crit
!p.multi=[0,nx,ny]
ps1,'./'+structure_filename+'_EUVIA-DEMT_Te-profiles.eps'
for ig=0,Ngroups-1 do begin
   ifl = where(line_groupID eq ig AND fitflag_EUVIA_A eq +1 AND tag_pos eq +1 AND scN_fit_euvia_A le scT_crit)
   if ifl(0) eq -1 then begin
      plot,CritTermLon(ig:ig+1),[1.02,1.25],charsize=csz,xtitle='Lon [deg]',title='Rad [Rsun]. Group #'+strmid(string(ig+1),7,1),th=4, ystyle=1,xstyle=1,/nodata
      goto,skip_group_euvia_Te
   endif
   Tm_fit_euvia_groupavg(ig,*) = total( Tm_fit_euvia_A(ifl,*) , 1 ) / float(n_elements(ifl))
   loadct,0
   plot,rad_fit_euvia_A,Tm_fit_euvia_groupavg(ig,*),charsize=csz,xtitle='r [Rsun]',title='EUVIA-DEMT Te(r) [MK]. Group #'+strmid(string(ig+1),7,1),th=4, yr=[0.5e6,3.0e6], ystyle=1, xr=[1.02,1.25], xstyle=1,/nodata
   loadct,ctbl
   color_index_step = fix(256./n_elements(ifl))
   for index=0,n_elements(ifl)-1 do begin
      tmp = reform(index_sampling_euvia_A(ifl(index),*))
      ind_samp_euvia = where(tmp eq 1)
      oplot,rad_fit_euvia_A                 ,Tm_fit_euvia_A(ifl(index),*)       ,color=(index)*color_index_step
      oplot,rad_A(ifl(index),ind_samp_euvia),Tm_euvia_A(ifl(index),ind_samp_euvia),color=(index)*color_index_step,psym=4
   endfor
   loadct,0
   oplot,rad_fit_euvia_A,Tm_fit_euvia_groupavg(ig,*),th=4
   skip_group_euvia_Te:
endfor
for ig=0,Ngroups-1 do begin
   ifl = where(line_groupID eq ig AND fitflag_EUVIA_A eq +1 AND tag_pos eq +1 AND scN_fit_EUVIA_A le scN_crit)
   if ifl(0) eq -1 then begin
      plot,CritTermLon(ig:ig+1),[1.02,1.25],charsize=csz,xtitle='Lon [deg]',title='Rad [Rsun]. Group #'+strmid(string(ig+1),7,1),th=4, ystyle=1,xstyle=1,/nodata
      goto,skip_group_euvia_geoT
   endif
   loadct,0
   plot,CritTermLon(ig:ig+1),[1.02,1.25],charsize=csz,xtitle='Lon [deg]',title='Rad [Rsun]. Group #'+strmid(string(ig+1),7,1),th=4, ystyle=1,xstyle=1,/nodata
   loadct,ctbl
   color_index_step = fix(256./n_elements(ifl))
   for index=0,n_elements(ifl)-1 do begin
      tmp = reform(index_sampling_euvia_A(ifl(index),*))
      ind_samp_euvia = where(tmp eq 1)
      oplot,lon_A(ifl(index),ind_samp_euvia),rad_A(ifl(index),ind_samp_euvia),color=(index)*color_index_step
   endfor
   loadct,0
   skip_group_euvia_geoT:
endfor
ps2
endif

if keyword_set(mk4) then begin
; -----------------------
     fitflag_mk4_A = *trace_data.fitflag_mk4
     rad_fit_mk4_A = *trace_data.rad_fit_mk4
     Ne_fit_mk4_A  = *trace_data.Ne_fit_mk4
     scN_fit_mk4_A = *trace_data.scN_fit_mk4
; -----------------------        
; Mk4 Ne
Ne_fit_mk4_groupavg = fltarr(Ngroups,n_elements(rad_fit_mk4_A))
;---- Tag field lines for which the fit is positive at all heights ----
   tag_pos = fltarr(N_fl)
   for ifl=0,N_fl-1 do begin
      if min(Ne_fit_mk4_A(ifl,*)) gt 0. then tag_pos(ifl)=+1
   endfor
;----------------------------------------------------------------------
scN_crit = rel_sqrt_chisqr_crit
!p.multi=[0,nx,ny]
ps1,'./'+structure_filename+'_Mk4-SRT_Ne-profiles.eps'
for ig=0,Ngroups-1 do begin
   ifl = where(line_groupID eq ig AND fitflag_mk4_A eq +1 AND tag_pos eq +1 AND scN_fit_mk4_A le scN_crit)
   if ifl(0) eq -1 then begin
      plot,rad_fit_mk4_A,Ne_fit_mk4_groupavg(ig,*),charsize=csz,xtitle='r [Rsun]',title='Mk4-SRT Ne(r) [cm!U-3!N]. Group #'+strmid(string(ig+1),7,1),th=4, ystyle=2, xr=[1.15,1.5], xstyle=1,/nodata
      goto,skip_group_mk4_Ne
   endif
   Ne_fit_mk4_groupavg(ig,*) = total( Ne_fit_mk4_A(ifl,*) , 1 ) / float(n_elements(ifl))
   loadct,0
   plot,rad_fit_mk4_A,Ne_fit_mk4_groupavg(ig,*),charsize=csz,xtitle='r [Rsun]',title='Mk4-SRT Ne(r) [cm!U-3!N]. Group #'+strmid(string(ig+1),7,1),th=4, ystyle=2, xr=[1.15,1.5], xstyle=1,/nodata
   loadct,ctbl
   color_index_step = fix(256./n_elements(ifl))
   for index=0,n_elements(ifl)-1 do begin
 ; ---------------------------------     
      Nsamp   = (*trace_data.Npt_mk4)(ifl(index))
      rad_mk4 = (*trace_data.rad_mk4)(ifl(index),0:Nsamp-1)
      Ne_mk4  = (*trace_data.Ne_mk4) (ifl(index),0:Nsamp-1)
      oplot,rad_fit_mk4_A                 ,Ne_fit_mk4_A(ifl(index),*),color=(index)*color_index_step
      oplot,rad_mk4,Ne_mk4,color=(index)*color_index_step,psym=4
 ; ---------------------------------     
   endfor
   loadct,0
   oplot,rad_fit_mk4_A,Ne_fit_mk4_groupavg(ig,*),th=4
   skip_group_mk4_Ne:
endfor
for ig=0,Ngroups-1 do begin
   ifl = where(line_groupID eq ig AND fitflag_mk4_A eq +1 AND tag_pos eq +1 AND scN_fit_mk4_A le ScN_crit)
   if ifl(0) eq -1 then begin
      plot,CrittermLon(ig:ig+1),[1.15,1.5],charsize=csz,xtitle='Lon [deg]',title='Rad [Rsun]. Group #'+strmid(string(ig+1),7,1),th=4, ystyle=1,xstyle=1,/nodata
      goto,skip_group_mk4_geo
   endif
   loadct,0
   plot,CrittermLon(ig:ig+1),[1.15,1.5],charsize=csz,xtitle='Lon [deg]',title='Rad [Rsun]. Group #'+strmid(string(ig+1),7,1),th=4, ystyle=1,xstyle=1,/nodata
   loadct,ctbl
   color_index_step = fix(256./n_elements(ifl))
   for index=0,n_elements(ifl)-1 do begin
  ; ---------------------------------     
      Nsamp =   (*trace_data.Npt_mk4)(ifl(index))
      rad_mk4 = (*trace_data.rad_mk4)(ifl(index),0:Nsamp-1)
      lon_mk4 = (*trace_data.lon_mk4)(ifl(index),0:Nsamp-1)
      oplot,lon_mk4,rad_mk4,color=(index)*color_index_step
   ; ---------------------------------        
   endfor
   loadct,0
   skip_group_mk4_geo:
endfor
ps2
endif


if keyword_set(lascoc2) then begin
; -----------------------
   fitflag_c2_A = *trace_data.fitflag_c2
   rad_fit_c2_A = *trace_data.rad_fit_c2
   Ne_fit_c2_A  = *trace_data.Ne_fit_c2
   scN_fit_c2_A = *trace_data.scN_fit_c2
; -----------------------        
; C2 Ne
Ne_fit_c2_groupavg = fltarr(Ngroups,n_elements(rad_fit_c2_A))
;---- Tag field lines for which the fit is positive at all heights ----
   tag_pos = fltarr(N_fl)
   for ifl=0,N_fl-1 do begin
      if min(Ne_fit_c2_A(ifl,*)) gt 0. then tag_pos(ifl)=+1
   endfor
;----------------------------------------------------------------------
scN_crit = rel_sqrt_chisqr_crit
!p.multi=[0,nx,ny]
ps1,'./'+structure_filename+'_C2-SRT_Ne-profiles.eps'
for ig=0,Ngroups-1 do begin
   ifl = where(line_groupID eq ig AND fitflag_C2_A eq +1 AND tag_pos eq +1 AND scN_fit_c2_A le scN_crit)
   if ifl(0) eq -1 then begin
         plot,rad_fit_c2_A,Ne_fit_c2_groupavg(ig,*),charsize=csz,xtitle='r [Rsun]',title='C2-SRT Ne(r) [cm!U-3!N]. Group #'+strmid(string(ig+1),7,1),th=4, ystyle=2, xr=[2.5,6.0], xstyle=1,/nodata
      goto,skip_group_c2_Ne
   endif
   Ne_fit_c2_groupavg(ig,*) = total( Ne_fit_c2_A(ifl,*) , 1 ) / float(n_elements(ifl))
   loadct,0
   plot,rad_fit_c2_A,Ne_fit_c2_groupavg(ig,*),charsize=csz,xtitle='r [Rsun]',title='C2-SRT Ne(r) [cm!U-3!N]. Group #'+strmid(string(ig+1),7,1),th=4, ystyle=2, xr=[2.5,6.0], xstyle=1,/nodata
   loadct,ctbl
   color_index_step = fix(256./n_elements(ifl))
   for index=0,n_elements(ifl)-1 do begin
      Nsamp   = (*trace_data.Npt_c2)(ifl(index))
      rad_c2 = (*trace_data.rad_c2)(ifl(index),0:Nsamp-1)
      Ne_c2  = (*trace_data.Ne_c2) (ifl(index),0:Nsamp-1)
      oplot,rad_fit_c2_A                 ,Ne_fit_c2_A(ifl(index),*),color=(index)*color_index_step
      oplot,rad_c2,Ne_c2,color=(index)*color_index_step,psym=4
 
   endfor
   loadct,0
   oplot,rad_fit_c2_A,Ne_fit_c2_groupavg(ig,*),th=4
   skip_group_c2_Ne:
endfor
for ig=0,Ngroups-1 do begin
   ifl = where(line_groupID eq ig AND fitflag_c2_A eq +1 AND tag_pos eq +1 AND scN_fit_c2_A le scN_crit)
   if ifl(0) eq -1 then begin
      plot,CritTermLon(ig:ig+1),[2.5,6.0],charsize=csz,xtitle='Lon [deg]',title='Rad [Rsun]. Group #'+strmid(string(ig+1),7,1),th=4, ystyle=1,xstyle=1,/nodata
      goto,skip_group_c2_geo
   endif
   loadct,0
   plot,CritTermLon(ig:ig+1),[2.5,6.0],charsize=csz,xtitle='Lon [deg]',title='Rad [Rsun]. Group #'+strmid(string(ig+1),7,1),th=4, ystyle=1,xstyle=1,/nodata
   loadct,ctbl
   color_index_step = fix(256./n_elements(ifl))
   for index=0,n_elements(ifl)-1 do begin
      Nsamp   = (*trace_data.Npt_c2)(ifl(index))
      rad_c2 = (*trace_data.rad_c2)(ifl(index),0:Nsamp-1)
      lon_c2 = (*trace_data.lon_c2)(ifl(index),0:Nsamp-1)
      oplot,lon_c2,rad_c2,color=(index)*color_index_step

   endfor
   loadct,0
   skip_group_c2_geo:
endfor
ps2
endif

for i_fl =0, N_fl-1 do begin


   if i_fl lt     10                     then ifl_str = '00' + strmid(string(long(i_fl)),11,1)
   if i_fl ge     10 and i_fl lt     100 then ifl_str = '0'  + strmid(string(long(i_fl)),10,2)
   if i_fl ge    100 and i_fl lt    1000 then ifl_str = ''   + strmid(string(long(i_fl)), 9,3)

 
   
   if keyword_set(aia) then begin
      Nsamp   = (*trace_data.Npt_aia)(i_fl)
      rad_aia = reform((*trace_data.rad_aia)(i_fl,0:Nsamp-1))
      Ne_aia  = reform((*trace_data.Ne_aia) (i_fl,0:Nsamp-1))
      rad_concat = rad_aia
      NE_concat  = Ne_aia
      inst_concat = 'aia' + strarr(n_elements(rad_aia))
      Inst_list   = ['aia']
   endif
   if keyword_set(mk4) then begin
      Nsamp   = (*trace_data.Npt_mk4)(i_fl)
      rad_mk4 = reform((*trace_data.rad_mk4)(i_fl,0:Nsamp-1))
      Ne_mk4  = reform((*trace_data.Ne_mk4) (i_fl,0:Nsamp-1))
      rad_concat = [rad_concat,rad_mk4]
      Ne_concat  = [Ne_concat,Ne_mk4]
      inst_concat = [inst_concat , 'mk4' + strarr(n_elements(rad_mk4))]
      if not keyword_set(Inst_list) then Inst_list=['mk4'] else Inst_list=[Inst_list,'mk4']
   endif
   if keyword_set(lascoc2) then begin
      Nsamp   = (*trace_data.Npt_c2)(i_fl)
      rad_c2  = reform((*trace_data.rad_c2)(i_fl,0:Nsamp-1))
      Ne_c2   = reform((*trace_data.Ne_c2) (i_fl,0:Nsamp-1))
      rad_concat = [rad_concat,rad_c2]
      Ne_concat  = [Ne_concat,Ne_c2]
      inst_concat = [inst_concat , 'c2' + strarr(n_elements(rad_c2))]
      if not keyword_set(Inst_list) then Inst_list=['mk4'] else Inst_list=[Inst_list,'c2']
   endif

   ps1,'line_'+ifl_str+'.eps'
   loadct,0   
   plot,rad_concat,Ne_concat,/ylog,/nodata,title='Field line #'+string(i_fl),xtitle='rad [Rsun]',ytitle='log!d10!N(Ne [cm!U-3!N] )'
   loadct,39
   color = fltarr(n_elements(Ne_concat))
   index = where(inst_concat eq 'aia')
   xpos = 0.2
   ypos = 0.9

   if index(0) ne -1 then begin
      color = 100
      oplot,rad_concat(index),Ne_concat(index),psym=4,color=color
      xyouts,[xpos],[ypos],['AIA'],color=color,/normal
   endif
   
   index = where(inst_concat eq 'mk4')
   if index(0) ne -1 then begin
      ypos = ypos - 0.03
      color = 200
      oplot,rad_concat(index),Ne_concat(index),psym=4,color=color
      xyouts,[xpos],[ypos],['Mk4'],color=color,/normal
   endif

   index = where(inst_concat eq 'c2' )
   if index(0) ne -1 then begin
      ypos = ypos - 0.03
      color = 250
      oplot,rad_concat(index),Ne_concat(index),psym=4,color=color
      xyouts,[xpos],[ypos],['C2'],color=color,/normal
   endif

   if index(0) ne -1 then oplot,rad_concat(index),Ne_concat(index),psym=4,color=250

   if keyword_set(npower) then $
   power_concat_fit, Inst_list, rad_concat, Ne_concat, inst_concat, A, chisq, /weighted
   if keyword_set(dpower) then $
   double_power_concat_fit, rad_concat, Ne_concat, inst_concat, A, chisq, /weighted
   if keyword_set(linfit) then $
      linear_fit, rad_concat, alog(Ne_concat), A, r2, chisqr = chisqr, /linfit_idl
   
   
   loadct,0
   Nrad_fit = 100
   rad_concat_fit = min(rad_concat) + (max(rad_concat)-min(rad_concat)) * findgen(Nrad_fit)/float(Nrad_fit-1)

   ypos = ypos -0.03
   xyouts,[xpos],[ypos],['Sqrt(Chisq)/Mean(Ne)='+strmid(string(sqrt(chisq)/mean(Ne_concat)),6,5)],/normal

  if keyword_set(npower) then begin
     if n_elements(Inst_list) eq 1 then $
        oplot, rad_concat_fit, A[0] / rad_concat_fit^A[1] ,th=2
     if n_elements(Inst_list) eq 2 then $
        oplot, rad_concat_fit, A[0] / rad_concat_fit^A[1] + A[2] / rad_concat_fit^A[3],th=2
     if n_elements(Inst_list) eq 3 then $
        oplot, rad_concat_fit, A[0] / rad_concat_fit^A[1] + A[2] / rad_concat_fit^A[3] + A[4] / rad_concat_fit^A[5], th=2
     ypos   = ypos -0.03
     xyouts, [xpos],[ypos],['Coeffs [cm!U-3!N]: '+string(A[0])+string(A[2])+string(A[4])],/normal
     ypos   = ypos -0.03
     xyouts, [xpos],[ypos],['Powers: '+string(A[1])+string(A[3])+string(A[5]) ],/normal
  endif

  if keyword_set(dpower) then begin
     oplot, rad_concat_fit, A[0] / rad_concat_fit^A[1] + A[2] / rad_concat_fit^A[3],th=2
     ypos   = ypos -0.03
     xyouts, [xpos],[ypos],['Coeffs [cm!U-3!N]: '+string(A[0])+string(A[2])],/normal
     ypos   = ypos -0.03
     xyouts, [xpos],[ypos],['Powers: '+string(A[1])+string(A[3]) ],/normal
  endif
   ps2
endfor


  return
end
 
PRO ps1,archivo
set_plot,'ps'
device,filename=archivo,bits_per_pixel=8,/color,/encapsulated
return
end

PRO ps2
device,/close
set_plot,'x'
!p.multi=0
return
end
