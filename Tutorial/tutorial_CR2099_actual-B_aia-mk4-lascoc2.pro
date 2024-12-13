pro mini_tutorial
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
  
; 1) Declare the DIR where the structure is located, and the filename.
  dir = './'
  structure_filename = 'CR2099_AWSoM-map1_tracing-structure-merge_aia_mk4_lascoc2.sav'
 ;structure_filename = 'CR2099_AWSoM-map7_tracing-structure-merge_aia_mk4_lascoc2.sav'

; 2) Load structure into memory and extract all available arrays from it.
  load_traced_data_structure, dir=dir, structure_filename=structure_filename, trace_data=trace_data, /aia, /mk4, /lascoc2

 goto,plots
 
print, 'Press SPACE BAR to continue.'
pause

; 3) See the full contents of the structure.
print
print, 'Bellow is the list of the full contents of the structure, named "trace_data", loaded into memory by the command:'
print
print, 'load_traced_data_structure, dir=dir, structure_filename=structure_filename, trace_data=trace_data, /aia, /mk4, /lascoc2'
print
print, 'For this particular structure only AIA, MK4 and C2 results exist.'
print
print,'-------------------------------'
print,'help, trace_data, /str'
print
print, 'Press SPACE BAR to continue.'
pause
print
help, trace_data, /str
print,'-------------------------------'
print

print, 'Press SPACE BAR to continue.'
pause

; 4) See the arrays extracted from the structure.
print
print, 'Besides loading into memory the structure above, the command "load_traced_data_structure...'
print, 'also creates in memory several variables and arrays with the results stored in the structure,'
print, 'so they are convenientely ready to be used.'
print
print, 'Let us see the output of the command "help" for all variables and arrays that MAY BE created from the structure:'
print
print, 'Press SPACE BAR to continue.'
pause
print
print,'-------------------------------'
help,N_fl, Npt_max, Npt_v, x_A, y_A, z_A, rad_A, lat_A, lon_A,$
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
print,'-------------------------------'
print

; 5) Some explanations follow.
print,'In this case only arrays wih filed line geometry, or pertaining AIA, Mk4 or C2, are defined.'
print
print,'In this case there are "N_fl=4" field lines that were traced between rmin=1 Rs and rmax=10 Rs. For each field line a maximum of "Npt_max=10100" points is allowed. This number will be optimized later on, depending on what we find to be a typical number of points for Judit field lines. For now we just set it as very large number. These are synthetic radial field lines with a uniform radial step of 0.0015 Rs, so that these four lines contain "Npt_v = [6000, 6000, 6000, 6000]" points (because (10-1)/0.0015=6000.'
print
print,'All variables named "name_A" are 2D arrays of dimensions "(N_fl,Npt_max)". For example: "rad_A(ifl,*)" contains the radial points of field line "ifl (=0,..,3 in this example)". Similarly, "Ne_aia_A(ifl,*)" contains the AIA-based DEMT result for Ne at points "rad_A(ifl,*)". Etc.'
print
print, 'Press SPACE BAR to continue.'
pause
print
print,'There is a KEY array, provided for each instrument, that deserves a separate explanation: "index_sampling_INSTRUMENT_A", with INSTRUMENT being in this case: "aia, mk4 or c2". For field line ifl, "index_sampling_INSTRUMENT_A(ifl,*)" contains values 0 and 1: the 1 values indicate a single point within a given tomographic cell grid that samples that cell: it is the mean point within the cell, and its value is 1 ONLY if there is a result for that INSTRUMENT in that cell. If its value is 0, that cell is off the reconstructed volume for that instrument, or the cell is within the volume reconstructed for that instrument but that specific cell has not been reconstructed (being either a ZDA for both VL-tomography and DEMT, or a location where the LDEM does not predict accurately all FBEs, in the case of DEMT only).' 
print
print, 'Press SPACE BAR to continue.'
pause

plots:

;;
; Determine the Lat and Lon of the footppoint and terminalpoint of each field line.
; Group field lines according to the Lon of their terminalpoint
;
; 1D Arrays: radial index corresponding to Rmin and Rmax for each field line:
irmin=intarr(N_FL)
irmax=intarr(N_FL)
for i=0,N_FL-1 do irmax(i)=where(    rad_A(i,*)  eq max(    rad_A(i,*)) )
for i=0,N_FL-1 do irmin(i)=where(abs(rad_A(i,*)) eq min(abs(rad_A(i,*))))
; 1D Arrays: Footpoint and Terminalpoint Lon and Lat for each field line:
Footpoint_Lon = fltarr(N_FL)
Footpoint_Lat = fltarr(N_FL)
 Terminal_Lon = fltarr(N_FL)
 Terminal_Lat = fltarr(N_FL)
for ifl = 0,N_FL-1 do Footpoint_Lon(ifl) = lon_A(ifl,irmin[ifl])
for ifl = 0,N_FL-1 do Footpoint_Lat(ifl) = lat_A(ifl,irmin[ifl])
for ifl = 0,N_FL-1 do  Terminal_Lon(ifl) = lon_A(ifl,irmax[ifl])
for ifl = 0,N_FL-1 do  Terminal_Lat(ifl) = lat_A(ifl,irmax[ifl])
;
; Tag groups of field lines by means of user-defined ranges of terminal Longitudes.
if structure_filename eq 'CR2099_AWSoM-map1_tracing-structure-merge_aia_mk4_lascoc2.sav' then CritTermLon = [0.,100.,180.,270.,360.]
if structure_filename eq 'CR2099_AWSoM-map7_tracing-structure-merge_aia_mk4_lascoc2.sav' then CritTermLon = [0.,100.,180.,270.,310.,360.]
;
Ngroups     = n_elements(CritTermLon)-1
; Tag each field line (ifl) with the group number (ig) to which it belongs: 
line_groupID  = intarr(N_FL)
for ig=0,Ngroups-1 do begin
   ifl_A = where(Terminal_Lon gt CritTermLon(ig) AND Terminal_Lon lt CritTermLon(ig+1))
   line_groupID(ifl_A) = ig
endfor
;;

; Set up graph options.
; color table for fieldlines

; Lat/Lon plots of FootPoint and TerminalPoint
;window,0,xs=1000,ys=1000
ps1,'./'+structure_filename+'_connectivity-map.eps'
np=1000
!p.multi=[0,1,2]
loadct,0
!p.color=0
!p.background=255
csz=1
plot,lon_A,lat_A,xr=[0,360],yr=[-90,+90],xstyle=1,ystyle=1,/nodata,charsize=csz,title=strmid(structure_filename,0,17)+'  r = 1.0 Rs',ytitle='Carrington Latitude [deg]'
ctbl = 12 ; 16-LEVEL
colors = 16 + 190 * (indgen(Ngroups))/float(Ngroups-1)
loadct,ctbl
for ifl=0,N_FL-1 do oplot,[Footpoint_Lon(ifl)],[Footpoint_Lat(ifl)],psym=4,color=colors(line_groupID(ifl)),th=2
loadct,0
plot,lon_A,lat_A,xr=[0,360],yr=[-90,+90],xstyle=1,ystyle=1,/nodata,charsize=csz,xtitle='Carrington Longitude [deg]',title='r = 23.68 Rs',ytitle='Carrington Latitude [deg]'
loadct,ctbl
for ifl=0,N_FL-1 do oplot,[Terminal_Lon(ifl)],[Terminal_Lat(ifl)],psym=4,color=colors(line_groupID(ifl)),th=2
loadct,0
!p.multi=0
ps2

;;
; Compute average trends <Ne(r)> and <Te(r)> for each group of field lines:
;if ngroups eq 4 then nx=2
;if ngroups eq 5 then nx=3
nx = ngroups
ny = 2
csz  = 0.75
cltb = 40

; AIA Ne
Ne_fit_aia_groupavg = fltarr(Ngroups,n_elements(rad_fit_aia_A))
;---- Tag field lines for which the fit is positive at all heights ----
   tag_pos = fltarr(N_fl)
   for ifl=0,N_fl-1 do begin
      if min(Ne_fit_aia_A(ifl,*)) gt 0. then tag_pos(ifl)=+1
   endfor
;----------------------------------------------------------------------
scN_crit = 0.5
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
   loadct,cltb
   color_index_step = fix(256./n_elements(ifl))
   for index=0,n_elements(ifl)-1 do begin
      tmp = reform(index_sampling_aia_A(ifl(index),*))
      ind_samp_aia = where(tmp eq 1)
      oplot,rad_fit_aia_A                 ,Ne_fit_aia_A(ifl(index),*)       ,color=(index)*color_index_step
      oplot,rad_A(ifl(index),ind_samp_aia),Ne_aia_A(ifl(index),ind_samp_aia),color=(index)*color_index_step,psym=4
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
      tmp = reform(index_sampling_aia_A(ifl(index),*))
      ind_samp_aia = where(tmp eq 1)
      oplot,lon_A(ifl(index),ind_samp_aia),rad_A(ifl(index),ind_samp_aia),color=(index)*color_index_step
   endfor
   loadct,0
   skip_group_aia_geo:
endfor
ps2

; Mk4 Ne
Ne_fit_mk4_groupavg = fltarr(Ngroups,n_elements(rad_fit_mk4_A))
;---- Tag field lines for which the fit is positive at all heights ----
   tag_pos = fltarr(N_fl)
   for ifl=0,N_fl-1 do begin
      if min(Ne_fit_mk4_A(ifl,*)) gt 0. then tag_pos(ifl)=+1
   endfor
;----------------------------------------------------------------------
scN_crit = 0.5
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
      tmp = reform(index_sampling_mk4_A(ifl(index),*))
      ind_samp_mk4 = where(tmp eq 1)
      oplot,rad_fit_mk4_A                 ,Ne_fit_mk4_A(ifl(index),*)       ,color=(index)*color_index_step
      oplot,rad_A(ifl(index),ind_samp_mk4),Ne_mk4_A(ifl(index),ind_samp_mk4),color=(index)*color_index_step,psym=4
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
      tmp = reform(index_sampling_mk4_A(ifl(index),*))
      ind_samp_mk4 = where(tmp eq 1)
      oplot,lon_A(ifl(index),ind_samp_mk4),rad_A(ifl(index),ind_samp_mk4),color=(index)*color_index_step
   endfor
   loadct,0
   skip_group_mk4_geo:
endfor
ps2

; C2 Ne
Ne_fit_c2_groupavg = fltarr(Ngroups,n_elements(rad_fit_c2_A))
;---- Tag field lines for which the fit is positive at all heights ----
   tag_pos = fltarr(N_fl)
   for ifl=0,N_fl-1 do begin
      if min(Ne_fit_c2_A(ifl,*)) gt 0. then tag_pos(ifl)=+1
   endfor
;----------------------------------------------------------------------
scN_crit = 0.5
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
   loadct,cltb
   color_index_step = fix(256./n_elements(ifl))
   for index=0,n_elements(ifl)-1 do begin
      tmp = reform(index_sampling_c2_A(ifl(index),*))
      ind_samp_c2 = where(tmp eq 1)
      oplot,rad_fit_c2_A                 ,Ne_fit_c2_A(ifl(index),*)      ,color=(index)*color_index_step
      oplot,rad_A(ifl(index),ind_samp_c2),Ne_c2_A(ifl(index),ind_samp_c2),color=(index)*color_index_step,psym=4
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
      tmp = reform(index_sampling_c2_A(ifl(index),*))
      ind_samp_c2 = where(tmp eq 1)
      oplot,lon_A(ifl(index),ind_samp_c2),rad_A(ifl(index),ind_samp_c2),color=(index)*color_index_step
   endfor
   loadct,0
   skip_group_c2_geo:
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
scT_crit = 0.5
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
   loadct,cltb
   color_index_step = fix(256./n_elements(ifl))
   for index=0,n_elements(ifl)-1 do begin
      tmp = reform(index_sampling_aia_A(ifl(index),*))
      ind_samp_aia = where(tmp eq 1)
      oplot,rad_fit_aia_A                 ,Tm_fit_aia_A(ifl(index),*)       ,color=(index)*color_index_step
      oplot,rad_A(ifl(index),ind_samp_aia),Tm_aia_A(ifl(index),ind_samp_aia),color=(index)*color_index_step,psym=4
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
      tmp = reform(index_sampling_aia_A(ifl(index),*))
      ind_samp_aia = where(tmp eq 1)
      oplot,lon_A(ifl(index),ind_samp_aia),rad_A(ifl(index),ind_samp_aia),color=(index)*color_index_step
   endfor
   loadct,0
   skip_group_aia_geoT:
endfor
ps2


  STOP
  return
end
