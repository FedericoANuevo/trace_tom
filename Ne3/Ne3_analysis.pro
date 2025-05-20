pro Ne3_analysis, load=load

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
     fit_F_Ne_euvia,fit_F_Ne_euvib,fit_F_eit_c2,$
     tomgrid_aia_hdr_A,tomgrid_aia_A,fitgrid_aia_hdr_A,fitgrid_aia_A,$
     tomgrid_euvia_hdr_A,tomgrid_euvia_A,fitgrid_euvia_hdr_A,fitgrid_euvia_A,$
     tomgrid_euvib_hdr_A,tomgrid_euvib_A,fitgrid_euvib_hdr_A,fitgrid_euvib_A,$
     tomgrid_mk4_hdr_A,tomgrid_mk4_A,fitgrid_mk4_hdr_A,fitgrid_mk4_A,$
     tomgrid_kcor_hdr_A,tomgrid_kcor_A,fitgrid_kcor_hdr_A,fitgrid_kcor_A,$
     tomgrid_ucomp_hdr_A,tomgrid_ucomp_A,fitgrid_ucomp_hdr_A,fitgrid_ucomp_A,$
     tomgrid_c2_hdr_A,tomgrid_c2_A,fitgrid_c2_hdr_A,fitgrid_c2_A,$
     leg_label_A

    dir = '/data1/DATA/trace_tom_files/CR2254/field_lines_geometry_aunifgrid_1.15Rs_5x5deg/'
    structure_filename='fdips_field_150x180x360_mrmqs220221t2004c2254_000.ubdat_fline-filenames_list.txt-tracing-structure-merge_aia_kcor_ucomp.sav'

    if keyword_set(load) then $
    load_traced_data_structure, dir=dir, structure_filename=structure_filename, trace_data=trace_data,/aia,/ucomp,/kcor,/opcl

; Define boxes of footpoint Lat/Lot to analyze
    Nbox = 2
    LonLimits = fltarr(Nbox,2)
    LatLimits = fltarr(Nbox,2)
   ;Box-0: 
    LonLimits(0,*) = [  0.,  90.]
    LatLimits(0,*) = [-30.,+ 30.]
   ;Box-1: 
    LonLimits(1,*) = [ 100., 150.]
    LatLimits(1,*) = [+ 60.,+ 80.]

; Determine the Lat and Lon of the footppoint of each field line.
; 1D Arrays: radial index corresponding to Rmin and Rmax for each field line:
    irmin=intarr(N_fl)
    for i=0,N_fl-1 do irmin(i)=where(abs(rad_A(i,*)) eq min(abs(rad_A(i,*))))
; 1D Arrays: Footpoint Lon and Lat for each field line:
    Footpoint_Lon = fltarr(N_fl)
    Footpoint_Lat = fltarr(N_fl)
    for ifl = 0,N_fl-1 do Footpoint_Lon(ifl) = lon_A(ifl,irmin[ifl])
    for ifl = 0,N_fl-1 do Footpoint_Lat(ifl) = lat_A(ifl,irmin[ifl])

    
; Tag each field line (ifl) with the BOX number (ibox) to which it belongs.
; If tag is -1 the line footpoint is not within any BOX.
    line_boxID  = intarr(N_fl) - 1.
    for ibox=0,Nbox-1 do begin
       ifl_A = where( (Footpoint_Lon ge LonLimits(ibox,0) AND Footpoint_Lon le LonLimits(ibox,1)) AND $
                      (Footpoint_Lat ge LatLimits(ibox,0) AND Footpoint_Lat le LatLimits(ibox,1)) )
       line_boxID(ifl_A) = ibox
    endfor

; Lat/Lon plots of FootPoints
ps1,'./'+structure_filename+'_footpoints-map.eps'
np=1000
!p.multi=0;[0,1,2]
loadct,0
!p.color=0
!p.background=255
csz=1
plot,lon_A,lat_A,xr=[0,360],yr=[-90,+90],xstyle=1,ystyle=1,/nodata,charsize=csz,title='Location of footpoints',ytitle='Carrington Latitude [deg]'
oplot,Footpoint_Lon,Footpoint_Lat,psym=4
loadct,12
green =  16
blue  = 100
red   = 200
ibox = 0
ifl  = where(line_boxID eq ibox)
oplot,Footpoint_Lon(ifl),Footpoint_Lat(ifl),psym=4,th=2,color=red
ibox = 1
ifl  = where(line_boxID eq ibox)
oplot,Footpoint_Lon(ifl),Footpoint_Lat(ifl),psym=4,th=2,color=blue
loadct,0
!p.multi=0
ps2

;---- Tag field lines for which the fit is positive at all heights ----
   tag_pos = fltarr(N_fl)
   for ifl=0,N_fl-1 do begin
      if min(Ne_fit_aia_A(ifl,*)) gt 0. then tag_pos(ifl)=+1
   endfor
;----------------------------------------------------------------------
scN_crit = 0.2
Ne_fit_aia_BoxAvg = fltarr(Nbox,n_elements(rad_fit_aia_A))
for ibox=0,Nbox-1 do begin
   ifl = where(line_BoxID eq ibox AND fitflag_AIA_A eq +1 AND tag_pos eq +1 AND scN_fit_aia_A le scN_crit)
   if ifl(0) ne -1 then begin
      Ne_fit_aia_BoxAvg(ibox,*) = total( Ne_fit_aia_A(ifl,*) , 1 ) / float(n_elements(ifl))
   endif
   if ifl(0) eq -1 then stop
endfor

for ibox=0,Nbox-1 do begin
   window,ibox
   plot,rad_fit_aia_A,Ne_fit_aia_BoxAvg(ibox,*)/1.e8,xr=[1.09,1.21],yr=[0.6,1.4],xstyle=1,ystyle=1
endfor

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
