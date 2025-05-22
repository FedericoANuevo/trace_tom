pro Ne3_analysis, load=load, LonLimits=LonLimits, LatLimits=LatLimits, plot_filename_suffix=plot_filename_suffix,$
                  aia=aia, kcor=kcor, ucomp=ucomp

    common data, N_fl, Npt_max, Npt_v, x_A, y_A, z_A, rad_A, lat_A, lon_A,$
     Ne_aia_A, Tm_aia_A, WT_aia_A, ldem_flag_aia_A, index_aia_A, index_sampling_aia_A,$
     Ne_euvia_A, Tm_euvia_A,  WT_euvia_A, ldem_flag_euvia_A, index_euvia_A, index_sampling_euvia_A,$
     Ne_euvib_A, Tm_euvib_A, WT_euvib_A, ldem_flag_euvib_A, index_euvib_A, index_sampling_euvib_A,$
     Ne_eit_A, Tm_eit_A, WT_eit_A, ldem_flag_eit_A, index_eit_A, index_sampling_eit_A,$
     Ne_mk4_A, index_mk4_A, index_sampling_mk4_A,$
     Ne_kcor_A, index_kcor_A, index_sampling_kcor_A,$
     Ne_ucomp_A, index_ucomp_A, index_sampling_ucomp_A,$
     Ne_c2_A, index_c2_A, index_sampling_c2_A,$
     rad_fit_aia_A, Ne_fit_aia_A, Tm_fit_aia_A, fitflag_aia_A,scN_fit_aia_A,scT_fit_aia_A,$
     rad_fit_euvia_A, Ne_fit_euvia_A, Tm_fit_euvia_A, fitflag_euvia_A,scN_fit_euvia_A,scT_fit_euvia_A,$
     rad_fit_euvib_A, Ne_fit_euvib_A, Tm_fit_euvib_A, fitflag_euvib_A,scN_fit_euvib_A,scT_fit_euvib_A,$
     rad_fit_eit_A, Ne_fit_eit_A, Tm_fit_eit_A, fitflag_eit_A,scN_fit_eit_A,scT_fit_eit_A,$
     rad_fit_c2_A, Ne_fit_c2_A, fitflag_c2_A,scN_fit_c2_A,$
     rad_fit_mk4_A, Ne_fit_mk4_A, fitflag_mk4_A,scN_fit_mk4_A,$
     rad_fit_kcor_A, Ne_fit_kcor_A, fitflag_kcor_A,scN_fit_kcor_A,$
     rad_fit_ucomp_A, Ne_fit_ucomp_A, fitflag_ucomp_A,scN_fit_ucomp_A,$
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
     fit_F_Ne_euvia,fit_F_Ne_euvib,fit_F_eit_c2,$
     tomgrid_aia_hdr_A,tomgrid_aia_A,fitgrid_aia_hdr_A,fitgrid_aia_A,$
     tomgrid_euvia_hdr_A,tomgrid_euvia_A,fitgrid_euvia_hdr_A,fitgrid_euvia_A,$
     tomgrid_euvib_hdr_A,tomgrid_euvib_A,fitgrid_euvib_hdr_A,fitgrid_euvib_A,$
     tomgrid_mk4_hdr_A,tomgrid_mk4_A,fitgrid_mk4_hdr_A,fitgrid_mk4_A,$
     tomgrid_kcor_hdr_A,tomgrid_kcor_A,fitgrid_kcor_hdr_A,fitgrid_kcor_A,$
     tomgrid_ucomp_hdr_A,tomgrid_ucomp_A,fitgrid_ucomp_hdr_A,fitgrid_ucomp_A,$
     tomgrid_c2_hdr_A,tomgrid_c2_A,fitgrid_c2_hdr_A,fitgrid_c2_A,$
     leg_label_A,$
     Footpoint_Rad_A, Footpoint_Lon_A, Footpoint_Lat_A
  
    dir = '/data1/DATA/trace_tom_files/CR2254/field_lines_geometry_aunifgrid_1.15Rs_5x5deg/'
    structure_filename='fdips_field_150x180x360_mrmqs220221t2004c2254_000.ubdat_fline-filenames_list.txt-tracing-structure-merge_aia_kcor_ucomp.sav'

    if keyword_set(load) then $
    load_traced_data_structure, dir=dir, structure_filename=structure_filename, trace_data=trace_data,/aia,/ucomp,/kcor,/opcl

; Define plot filename suffix
  if not keyword_set(plot_filename_suffix) then plot_filename_suffix='footpoints-map'

; Define box of footpoint Lat/Lot to analyze
  if not keyword_set(LonLimits) then LonLimits = [  0., 360.]
  if not keyword_set(LatLimits) then LatLimits = [-90.,+ 90.]

; Tag field lines with footpoints within the BOX.
    tag_box_A = fltarr(N_fl) - 678.
    ifl_A = where( (Footpoint_Lon_A ge LonLimits(0) AND Footpoint_Lon_A le LonLimits(1)) AND $
                   (Footpoint_Lat_A ge LatLimits(0) AND Footpoint_Lat_A le LatLimits(1)) )
    tag_box_A(ifl_A) = +1.

; Tag field lines for which the fit is positive at all heights below
; R_max, as defined below.
; Also,
; Tag field lines for which the fit is positive at least for Nsamp
; points below R_crit, and for Nsamp values in range [R_crit,R_max],
; as defined next.
Nsamp = 2
; Set R_crit and R_max using the "bottleneck" instrument, achieved by organizing
; these conditionals in increasing "bottleneckless".
;
if keyword_set(kcor ) then begin
   r_crit = median(rad_fit_kcor_A )
   r_max  = max   (rad_fit_kcor_A )
endif
if keyword_set(aia  ) then begin
   r_crit = median(rad_fit_aia_A  )
   r_max  = max   (rad_fit_aia_A  )
endif
if keyword_set(ucomp) then begin
   r_crit = median(rad_fit_ucomp_A)
   r_max  = max   (rad_fit_ucomp_A)
endif
;
tag_pos_aia_A          = fltarr(N_fl) - 678.
tag_fullrange_aia_A    = fltarr(N_fl) - 678.
tag_pos_kcor_A         = fltarr(N_fl) - 678.
tag_fullrange_kcor_A   = fltarr(N_fl) - 678.
tag_pos_ucomp_A        = fltarr(N_fl) - 678.
tag_fullrange_ucomp_A  = fltarr(N_fl) - 678.
for ifl=0,N_fl-1 do begin
   if keyword_set(aia) then begin
      ifitpos_aia = where(reform(Ne_fit_aia_A  (ifl,*)) gt 0. and rad_fit_aia_A le R_max)
      ifitrad_aia = where(rad_fit_aia_A le R_max) 
      if ifitpos_aia(0) ne -1 then begin
         if n_elements(ifitpos_aia) eq n_elements(ifitrad_aia) then tag_pos_aia_A(ifl) = +1.
         if n_elements(where(rad_fit_aia_A(ifitpos_aia) lt r_crit))                                         ge Nsamp AND  $
            n_elements(where(rad_fit_aia_A(ifitpos_aia) gt r_crit AND rad_fit_aia_A(ifitpos_aia) lt r_max)) ge Nsamp then tag_fullrange_aia_A(ifl) = +1.
      endif
   endif else begin
      tag_pos_aia_A(ifl) = +1
      tag_fullrange_aia_A(ifl) = +1
   endelse
   if keyword_set(kcor) then begin
      ifitpos_kcor = where(reform(Ne_fit_kcor_A (ifl,*)) gt 0. and rad_fit_kcor_A le R_max)
      ifitrad_kcor = where(rad_fit_kcor_A le R_max) 
      if ifitpos_kcor(0) ne -1 then begin
         if n_elements(ifitpos_kcor) eq n_elements(ifitrad_kcor) then tag_pos_kcor_A(ifl) = +1.
         if n_elements(where(rad_fit_kcor_A(ifitpos_kcor) lt r_crit))                                           ge Nsamp AND  $
            n_elements(where(rad_fit_kcor_A(ifitpos_kcor) gt r_crit AND rad_fit_kcor_A(ifitpos_kcor) lt r_max)) ge Nsamp then tag_fullrange_kcor_A(ifl) = +1.
      endif
   endif else begin
      tag_pos_kcor_A(ifl) = +1
      tag_fullrange_kcor_A(ifl) = +1
   endelse
   if keyword_set(ucomp) then begin
      ifitpos_ucomp = where(reform(Ne_fit_ucomp_A(ifl,*)) gt 0. and rad_fit_ucomp_A le R_max)
      ifitrad_ucomp = where(rad_fit_ucomp_A le R_max) 
      if ifitpos_ucomp(0) ne -1 then begin
         if n_elements(ifitpos_ucomp) eq n_elements(ifitrad_ucomp) then tag_pos_ucomp_A(ifl) = +1.
         if n_elements(where(rad_fit_ucomp_A(ifitpos_ucomp) lt r_crit))                                             ge Nsamp AND  $
            n_elements(where(rad_fit_ucomp_A(ifitpos_ucomp) gt r_crit AND rad_fit_ucomp_A(ifitpos_ucomp) lt r_max)) ge Nsamp then tag_fullrange_ucomp_A(ifl) = +1.
      endif
   endif else begin
      tag_pos_ucomp_A(ifl) = +1
      tag_fullrange_ucomp_A(ifl) = +1
   endelse
endfor

;-----------------PLOTS SECTION--------------------------------------------------------------
; Define a few color codes.
blue  = 100
red   = 200
green =  16

; Lat/Lon plots of FootPoints
ps1,'./'+structure_filename+'_'+plot_filename_suffix+'.eps'
np=1000
!p.multi=[0,1,2]
loadct,0
!p.color=0
!p.background=255

csz=1
plot,lon_A,lat_A,xr=[0,360],yr=[-90,+90],xstyle=1,ystyle=1,/nodata,charsize=csz,$
     title='Location of footpoints',ytitle='Carrington Latitude [deg]',font=0
oplot,Footpoint_Lon_A,Footpoint_Lat_A,psym=4

loadct,12
; Highlight in blue field lines with footpoints within BOX.
ifl_A  = where(tag_box_A eq +1)
oplot,Footpoint_Lon_A(ifl_A),Footpoint_Lat_A(ifl_A),psym=4,th=2,color=blue

; Set a maximum threshold for scN
scN_crit = 0.2

; Highlight in red Box field lines for which tag_pos_INSTRUMENT_A = +1 and there
; is a fit with scN < scN_crit for all instruments.
ifl_aia_A   = indgen(N_fl)
ifl_kcor_A  = indgen(N_fl)
ifl_ucomp_A = indgen(N_fl)
if keyword_set(aia)   then ifl_aia_A   = where(tag_box_A eq +1 AND tag_pos_aia_A   eq +1 AND fitflag_aia_A   eq +1 AND scN_fit_aia_A   le scN_crit)
if keyword_set(kcor)  then ifl_kcor_A  = where(tag_box_A eq +1 AND tag_pos_kcor_A  eq +1 AND fitflag_kcor_A  eq +1 AND scN_fit_kcor_A  le scN_crit)
if keyword_set(ucomp) then ifl_ucomp_A = where(tag_box_A eq +1 AND tag_pos_ucomp_A eq +1 AND fitflag_ucomp_A eq +1 AND scN_fit_ucomp_A le scN_crit)
; Make an index array with the common elements of all ifl_INSTRUMENT_A
                           ifl_A = indgen(N_fl)
if keyword_set(aia)   then ifl_A = intersect(ifl_A,ifl_aia_A  )
if keyword_set(kcor)  then ifl_A = intersect(ifl_A,ifl_kcor_A )
if keyword_set(ucomp) then ifl_A = intersect(ifl_A,ifl_ucomp_A)

oplot,Footpoint_Lon_A(ifl_A),Footpoint_Lat_A(ifl_A),psym=4,th=2,color=red
; Compute the average Ne(r) of each instrument for lines indexed ifl_A
if keyword_set(aia) then begin
   Ne_fit_aia_Avg = fltarr(n_elements(rad_fit_aia_A)) - 678.
   Ne_fit_aia_Avg = total( Ne_fit_aia_A(ifl_A,*) , 1 ) / float(n_elements(ifl_A))
endif
if keyword_set(kcor) then begin
   Ne_fit_kcor_Avg = fltarr(n_elements(rad_fit_kcor_A)) - 678.
   Ne_fit_kcor_Avg = total( Ne_fit_kcor_A(ifl_A,*) , 1 ) / float(n_elements(ifl_A))
endif
if keyword_set(ucomp) then begin
   Ne_fit_ucomp_Avg = fltarr(n_elements(rad_fit_ucomp_A)) - 678.
   Ne_fit_ucomp_Avg = total( Ne_fit_ucomp_A(ifl_A,*) , 1 ) / float(n_elements(ifl_A))
endif
skip_tag_pos:

goto,skip_tag_fullrange
; Highlight in red Box field lines for which tag_fullrange_A = +1 and there is a fit with scN < scN_crit
ifl_A = where(tag_box_A eq +1 AND tag_fullrange_A eq +1 AND fitflag_AIA_A eq +1 AND scN_fit_aia_A le scN_crit)
oplot,Footpoint_Lon_A(ifl_A),Footpoint_Lat_A(ifl_A),psym=4,th=2,color=red
; Compute and plot the average Ne(r) of the RED field lines
Ne_fit_aia_Avg = fltarr(n_elements(rad_fit_aia_A)) - 678.
index_selected_A = intarr(N_fl)
for ir=0,n_elements(rad_fit_aia_A)-1 do begin
   Nptsavg = 0
   for ifl=0,N_fl-1 do begin
      if tag_box_A[ifl] eq +1 AND tag_fullrange_A[ifl] eq +1 AND fitflag_AIA_A[ifl] eq +1 AND scN_fit_aia_A[ifl] le scN_crit then begin
         index_selected_A(ifl) = +1
         if Ne_fit_aia_A(ifl,ir) gt 0. then begin
            Nptsavg = Nptsavg + 1
            Ne_fit_aia_Avg(ir) = Ne_fit_aia_Avg(ir) + Ne_fit_aia_A(ifl,ir)
         endif
      endif      
   endfor ; ifl
   Ne_fit_aia_Avg(ir) = Ne_fit_aia_Avg(ir) / float(Nptsavg)
endfor ; ir
ifl_A = where(index_selected_A eq +1)
oplot,Footpoint_Lon_A(ifl_A),Footpoint_Lat_A(ifl_A),psym=4,th=2,color=red
skip_tag_fullrange:

; Plot the average Ne(r) for all selected instruments.

xrange = [1.095,1.195]
yrange = [0.63 ,1.33 ]
unit           = 1.e8
unit_power_str =   '8'
plot,rad_fit_kcor_A,Ne_fit_kcor_A(0,*)/unit,charsize=csz,font=0,$
     title='Average N!De!N(r) along red-colored field lines',$
     ytitle = 'Ne(r) [x 10!U'+unit_power_str+'!N cm!U-3!N]',yr=yrange,ystyle=1,$
     xtitle = 'r [Rsun]'                                   ,xr=xrange,xstyle=1,$
     /nodata

loadct,12
if keyword_set(aia)   then oplot,rad_fit_aia_A  ,Ne_fit_aia_Avg  /unit,color=blue ,th=2
if keyword_set(kcor)  then oplot,rad_fit_kcor_A ,Ne_fit_kcor_Avg /unit,color=red  ,th=2
if keyword_set(ucomp) then oplot,rad_fit_ucomp_A,Ne_fit_ucomp_Avg/unit,color=green,th=2
loadct,0
!p.multi=0
ps2

stop
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
