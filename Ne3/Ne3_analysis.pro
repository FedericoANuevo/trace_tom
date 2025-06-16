pro Ne3_analysis, load=load, LonLimits=LonLimits, LatLimits=LatLimits, $
                  plot_filename_suffix=plot_filename_suffix,$
                  aia=aia, kcor=kcor, ucomp=ucomp,$
                  plotaia=plotaia, plotkcor=plotkcor, plotucomp=plotucomp,$
                  open=open,closed=closed,connect=connect,$
                  positparam=positparam, ilstep=ilstep, constep=constep,$
                  r_max=r_max, only_loops=only_loops

  common data, N_fl, Npt_max, Npt_v, x_A, y_A, z_A, s_A, Br_A, Bth_A, Bph_A, B_A, rad_A, lat_A, lon_A,$
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
     N1_fit_kcor_A,N2_fit_kcor_A,p1_fit_kcor_A,p2_fit_kcor_A,$
     N1_fit_ucomp_A,N2_fit_ucomp_A,p1_fit_ucomp_A,p2_fit_ucomp_A,$
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
     leg_label_A,leg_footbfield_A,leg_length_A,$
     Footpoint_Rad_A, Footpoint_Lon_A, Footpoint_Lat_A

;===============================================================================================
; Select the project to analyze:
  PROJECT_NAME                   = 'CR2254'

; Select structure to read:
 ;structure_filename='fdips_field_150x180x360_mrmqs220221t2004c2254_000.ubdat_fline-filenames_list.txt-tracing-structure-merge_aia_kcor_ucomp.sav'
 ;structure_filename='fdips_field_150X180X360_mrmqs220831t1302c2261_000.ubdat_fline-filenames_list.txt-tracing-structure-merge_aia_kcor_ucomp.sav'
 ;structure_filename='fdips_field_150X180X360_hmi.Synoptic_Mr.2254.ubdat_fline-filenames_list.txt-tracing-structure-merge_aia_kcor_ucomp.sav'
  structure_filename='fdips_field_150X180X360_hmi.Synoptic_Mr_polfil.2254.ubdat_fline-filenames_list.txt-tracing-structure-merge_aia_kcor_ucomp.sav'

; Select dir where the structure is located (labeled after the selection of starting points) 
 ;field_line_geometry_suffix_dir = '_aunifgrid_multirad_5x5deg_HMI/'
 ;field_line_geometry_suffix_dir = '_aunifgrid_multirad_5x5deg_HMI-PolFil/'
 ;field_line_geometry_suffix_dir = '_aunifgrid_2.50Rs_2x2deg_HMI/'
  field_line_geometry_suffix_dir = '_aunifgrid_2.50Rs_2x2deg_HMI-PolFil/'
;===============================================================================================
  
; Load structure if so requested:
  dir = '/data1/DATA/trace_tom_files/'+PROJECT_NAME+'/field_lines_geometry'+field_line_geometry_suffix_dir
  if keyword_set(load) then $
     load_traced_data_structure, dir=dir, structure_filename=structure_filename, trace_data=trace_data,/aia,/ucomp,/kcor,/trace_Bs

; Define plot filename suffix
  if not keyword_set(plot_filename_suffix) then plot_filename_suffix='footpoints-map'

; Default index line step and connect step, for plots only, this does
; NOT affect average trend computations
  if not keyword_set( ilstep) then  ilstep=1
  if not keyword_set(constep) then constep=1
  
; Define box of footpoint Lat/Lot to analyze
  if not keyword_set(LonLimits) then LonLimits = [  0., 360.]
  if not keyword_set(LatLimits) then LatLimits = [-90.,+ 90.]

; Define a polarity tag
  tag_polarity_A = reform(leg_footBfield_A(*,0)/abs(leg_footBfield_A(*,0)))
  if ( where(finite(tag_polarity_A) eq 0) )(0) ne -1 then STOP  ; assumes foot-Br is not ZERO
  
; Tag field lines with footpoints within the BOX., either open or closed.
    tag_box_A = fltarr(N_fl)
    index_box = where( (Footpoint_Lon_A ge LonLimits(0) AND Footpoint_Lon_A le LonLimits(1)) AND $
                       (Footpoint_Lat_A ge LatLimits(0) AND Footpoint_Lat_A le LatLimits(1)) )
    tag_box_A(index_box) = +1.

; Now, UNTAG CLOSED field lines whose other leg's footpoint is NOT wihin the BOX,
;      as well as CLOSED field lines that do NOT comply with opposite polarity.
 if keyword_Set(only_loops) then begin
    ifl=0
    while ifl le N_fl-2 do begin
       if leg_label_A(ifl) eq 0. then begin
          ifl=ifl+1
       endif else begin
          if leg_label_A(ifl) ne leg_label_A(ifl+1) then STOP ; !this should never happen.
          if (tag_box_A(ifl) eq +1 AND tag_box_A(ifl+1) eq  0) then tag_box_A(ifl  )=0
          if (tag_box_A(ifl) eq  0 AND tag_box_A(ifl+1) eq +1) then tag_box_A(ifl+1)=0
          if (tag_polarity_A(ifl) eq tag_polarity_A(ifl+1))    then tag_box_A(ifl:ifl+1)=0 
          ifl=ifl+2
       endelse
    endwhile
 endif
 
; Create index arrays for closed and open field lines.
  ifl_closed_A = where(leg_label_A ne 0)
  ifl_open_A   = where(leg_label_A eq 0)
 
; Create index arrays for positive/negative footpoint Brad lines
  ifl_pos_A = where(tag_polarity_A eq +1.)
  ifl_neg_A = where(tag_polarity_A eq -1.)

; Tag field lines for which the fit is positive at all heights below
; R_max, as defined below.
; Also,
; Tag field lines for which the fit is positive at least for Nsamp
; points below R_crit, and for Nsamp values in range [R_crit,R_max],
; as defined next.
Nsamp = 2
; Set two parameters, named R_crit and R_max, using the "bottleneck" instrument,
; achieved by organizing these conditionals in increasing "bottleneckless".
;
if NOT keyword_set(r_max) then begin
   stop
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
endif
;
 tag_pos_aia_A          = fltarr(N_fl) - 678.
;tag_fullrange_aia_A    = fltarr(N_fl) - 678.
 tag_pos_kcor_A         = fltarr(N_fl) - 678.
;tag_fullrange_kcor_A   = fltarr(N_fl) - 678.
 tag_pos_ucomp_A        = fltarr(N_fl) - 678.
;tag_fullrange_ucomp_A  = fltarr(N_fl) - 678.
for ifl=0,N_fl-1 do begin
   if keyword_set(aia) then begin
      ifitpos_aia = where(reform(Ne_fit_aia_A  (ifl,*)) gt 0. and rad_fit_aia_A le R_max)
      ifitrad_aia = where(rad_fit_aia_A le R_max) 
      if ifitpos_aia(0) ne -1 then begin
         if n_elements(ifitpos_aia) eq n_elements(ifitrad_aia) then tag_pos_aia_A(ifl) = +1.
;         if n_elements(where(rad_fit_aia_A(ifitpos_aia) lt r_crit))                                         ge Nsamp AND  $
;            n_elements(where(rad_fit_aia_A(ifitpos_aia) gt r_crit AND rad_fit_aia_A(ifitpos_aia) lt r_max)) ge Nsamp then tag_fullrange_aia_A(ifl) = +1.
      endif
   endif else begin
      tag_pos_aia_A(ifl) = +1
;      tag_fullrange_aia_A(ifl) = +1
   endelse
   if keyword_set(kcor) then begin
      ifitpos_kcor = where(reform(Ne_fit_kcor_A (ifl,*)) gt 0. and rad_fit_kcor_A le R_max)
      ifitrad_kcor = where(rad_fit_kcor_A le R_max) 
      if ifitpos_kcor(0) ne -1 then begin
         if n_elements(ifitpos_kcor) eq n_elements(ifitrad_kcor) then tag_pos_kcor_A(ifl) = +1.
;         if n_elements(where(rad_fit_kcor_A(ifitpos_kcor) lt r_crit))                                           ge Nsamp AND  $
;            n_elements(where(rad_fit_kcor_A(ifitpos_kcor) gt r_crit AND rad_fit_kcor_A(ifitpos_kcor) lt r_max)) ge Nsamp then tag_fullrange_kcor_A(ifl) = +1.
      endif
   endif else begin
      tag_pos_kcor_A(ifl) = +1
;      tag_fullrange_kcor_A(ifl) = +1
   endelse
   if keyword_set(ucomp) then begin
      ifitpos_ucomp = where(reform(Ne_fit_ucomp_A(ifl,*)) gt 0. and rad_fit_ucomp_A le R_max)
      ifitrad_ucomp = where(rad_fit_ucomp_A le R_max) 
      if ifitpos_ucomp(0) ne -1 then begin
         if n_elements(ifitpos_ucomp) eq n_elements(ifitrad_ucomp) then tag_pos_ucomp_A(ifl) = +1.
;         if n_elements(where(rad_fit_ucomp_A(ifitpos_ucomp) lt r_crit))                                             ge Nsamp AND  $
;            n_elements(where(rad_fit_ucomp_A(ifitpos_ucomp) gt r_crit AND rad_fit_ucomp_A(ifitpos_ucomp) lt r_max)) ge Nsamp then tag_fullrange_ucomp_A(ifl) = +1.
      endif
   endif else begin
      tag_pos_ucomp_A(ifl) = +1
;      tag_fullrange_ucomp_A(ifl) = +1
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
                            opcl_str='open/closed'
if keyword_set(open)   then opcl_str='open '
if keyword_set(closed) then opcl_str='closed '
plot,lon_A,lat_A,xr=[0,360],yr=[-90,+90],xstyle=1,ystyle=1,/nodata,charsize=csz,$
     title='Location of '+opcl_str+'footpoints',ytitle='Carrington Latitude [deg]',font=0
if not keyword_set(open) and not keyword_set(closed) then oplot,Footpoint_Lon_A              ,Footpoint_Lat_A              ,psym=4
if     keyword_set(open)                             then oplot,Footpoint_Lon_A(ifl_open_A)  ,Footpoint_Lat_A(ifl_open_A)  ,psym=4
if                               keyword_set(closed) then oplot,Footpoint_Lon_A(ifl_closed_A),Footpoint_Lat_A(ifl_closed_A),psym=4

loadct,12

; Set a maximum threshold for scN
scN_crit = 0.25

; Independently for each instrument, tag field lines for which their DPL Ne-fit has all parameters positive
tag_posfit_aia_A = fltarr(N_Fl)
index            = where(N1_fit_aia_A gt 0. AND N2_fit_aia_A gt 0. AND p1_fit_aia_A gt 0. AND p2_fit_aia_A gt 0.)
tag_posfit_aia_A(index) = +1
;
tag_posfit_kcor_A = fltarr(N_Fl)
index             = where(N1_fit_kcor_A gt 0. AND N2_fit_kcor_A gt 0. AND p1_fit_kcor_A gt 0. AND p2_fit_kcor_A gt 0.)
tag_posfit_kcor_A(index) = +1
;
tag_posfit_ucomp_A = fltarr(N_Fl)
index              = where(N1_fit_ucomp_A gt 0. AND N2_fit_ucomp_A gt 0. AND p1_fit_ucomp_A gt 0. AND p2_fit_ucomp_A gt 0.)
tag_posfit_ucomp_A(index) = +1

; Independently for each instrument, index lines that:
; 1) are in the box, 2) have a fit, 3) have Ne_fit>0 up to R_max, and 4) have a low chisq fit to Ne
ifl_aia_A   = where(tag_box_A eq +1 AND tag_pos_aia_A   eq +1 AND fitflag_aia_A   eq +1 AND scN_fit_aia_A   le scN_crit)
ifl_kcor_A  = where(tag_box_A eq +1 AND tag_pos_kcor_A  eq +1 AND fitflag_kcor_A  eq +1 AND scN_fit_kcor_A  le scN_crit)
ifl_ucomp_A = where(tag_box_A eq +1 AND tag_pos_ucomp_A eq +1 AND fitflag_ucomp_A eq +1 AND scN_fit_ucomp_A le scN_crit)
; Also independently for each instrument,
; filter OUT lines for which NOT all their DPL Ne-fit parameters are positive
if keyword_set(positparam) then begin
   indpositparam_aia   = where(tag_posfit_aia_A   eq +1)  &  ifl_aia_A   = intersect(ifl_aia_A  ,indpositparam_aia  ) 
   indpositparam_kcor  = where(tag_posfit_kcor_A  eq +1)  &  ifl_kcor_A  = intersect(ifl_kcor_A ,indpositparam_kcor ) 
   indpositparam_ucomp = where(tag_posfit_ucomp_A eq +1)  &  ifl_ucomp_A = intersect(ifl_ucomp_A,indpositparam_ucomp) 
endif

; Make an index array with the common elements of all ifl_INSTRUMENT_A
                            ifl_A = indgen(N_fl) ; start index with ALL lines
if keyword_set(aia)    then ifl_A = intersect(ifl_A,ifl_aia_A  )
if keyword_set(kcor)   then ifl_A = intersect(ifl_A,ifl_kcor_A )
if keyword_set(ucomp)  then ifl_A = intersect(ifl_A,ifl_ucomp_A)
; Intersect with CLOSED or OPEN, if so requested
if keyword_set(closed) then ifl_A = intersect(ifl_A,ifl_closed_A)
if keyword_set(open)   then ifl_A = intersect(ifl_A,ifl_open_A  )

; Color-highlight all footpoints indicated by ifl_A accordind to their polarity
indxpos_A = intersect(ifl_A,ifl_pos_A)
if n_elements(indxpos_A) gt 1 then $
oplot,Footpoint_Lon_A(indxpos_A),Footpoint_Lat_A(indxpos_A),psym=4,th=2,color=green
indxneg_A = intersect(ifl_A,ifl_neg_A)
if n_elements(indxneg_A) gt 1 then $
oplot,Footpoint_Lon_A(indxneg_A),Footpoint_Lat_A(indxneg_A),psym=4,th=2,color=red

; Add to plot the connectivity of closed loops if requested
if keyword_set(connect) AND keyword_set(closed) then begin
   ifl=0
   closed_loop_count = -1
   while ifl le N_fl-2 do begin
      if leg_label_A(ifl) ne 0 then begin
         closed_loop_count = closed_loop_count + 1
         if (tag_box_A(ifl) eq +1) AND (closed_loop_count mod constep eq 0) then $
            oplot,reform(Footpoint_Lon_A(ifl:ifl+1)),reform(Footpoint_Lat_A(ifl:ifl+1)),color=blue
         ifl=ifl+2
      endif else begin
         ifl=ifl+1
      endelse
   endwhile
endif

; Compute the average Ne(r) of each instrument for lines indexed ifl_A
if keyword_set(aia) then begin
   Ne_fit_aia_Avg   = total( Ne_fit_aia_A(ifl_A,*)   , 1 ) / float(n_elements(ifl_A)) 
endif
if keyword_set(kcor) then begin
   Ne_fit_kcor_Avg  = total( Ne_fit_kcor_A(ifl_A,*)  , 1 ) / float(n_elements(ifl_A))
endif
if keyword_set(ucomp) then begin
   Ne_fit_ucomp_Avg = total( Ne_fit_ucomp_A(ifl_A,*) , 1 ) / float(n_elements(ifl_A))
endif
skip_tag_pos:

; Plot the average Ne(r) for all selected instruments.
loadct,0
unit           = 1.e8 ; cm-3
unit_power_str =   '8'
xrange = [1.095,r_max]
yrflag = -1
if keyword_set(aia) then begin
   r    = rad_fit_aia_A
   f    = Ne_fit_aia_Avg/unit
   miny = min(f(where(f gt 0. and r le max(xrange))))
   maxy = max(f(where(f gt 0. and r le max(xrange))))
   yrange = [miny,maxy]
   yrflag = +1
endif
if keyword_set(kcor) then begin
   r    = rad_fit_kcor_A
   f    = Ne_fit_kcor_Avg/unit
   miny = min(f(where(f gt 0. and r le max(xrange))))
   maxy = max(f(where(f gt 0. and r le max(xrange))))
   if yrflag eq -1. then yrange = [miny,maxy]
   if yrflag ne -1. then yrange = [min([miny,yrange(0)]),max([maxy,max(yrange(1))])]
   yrflag = +1
endif
if keyword_set(ucomp) then begin
   r    = rad_fit_ucomp_A
   f    = Ne_fit_ucomp_Avg/unit
   miny = min(f(where(f gt 0. and r le max(xrange))))
   maxy = max(f(where(f gt 0. and r le max(xrange))))
   if yrflag eq -1. then yrange = [miny,maxy]
   if yrflag ne -1. then yrange = [min([miny,yrange(0)]),max([maxy,max(yrange(1))])]
   yrflag = +1
endif
plot,rad_fit_kcor_A,Ne_fit_kcor_A(0,*)/unit,charsize=csz,font=0,$
     title  = '<N!De!N(r)>   Solid: AIA; Dashed: KCOR; Dot-Dashed: UCoMP',$
     ytitle = 'Ne(r) [x 10!U'+unit_power_str+'!N cm!U-3!N]',yr=yrange,ystyle=1,$
     xtitle = 'r [Rsun]'                                   ,xr=xrange,xstyle=1,$
     /nodata

if keyword_set(plotaia) or keyword_set(plotkcor) or keyword_set(plotucomp) then begin
loadct,39
Nifl = n_elements(ifl_A)
col = 255 * findgen(Nifl)/float(Nifl-1)
for ifl=0,Nifl-1,ilstep do begin
   il = (ifl_A(ifl))(0)
   if keyword_set(plotaia) then begin
      tmp      = reform(index_sampling_aia_A(il,*))
      ind_samp = where(tmp eq 1)
      oplot,rad_fit_aia_A     ,Ne_fit_aia_A(il,*)   /unit,color=col(ifl)
      oplot,rad_A(il,ind_samp),Ne_aia_A(il,ind_samp)/unit,color=col(ifl),psym=4
   endif
   if keyword_set(plotkcor) then begin
      tmp      = reform(index_sampling_kcor_A(il,*))
      ind_samp = where(tmp eq 1)
      oplot,rad_fit_kcor_A    ,Ne_fit_kcor_A(il,*)   /unit,color=col(ifl)
      oplot,rad_A(il,ind_samp),Ne_kcor_A(il,ind_samp)/unit,color=col(ifl),psym=4
   endif
   if keyword_set(plotucomp) then begin
      tmp      = reform(index_sampling_ucomp_A(il,*))
      ind_samp = where(tmp eq 1)
      oplot,rad_fit_ucomp_A   ,Ne_fit_ucomp_A(il,*)   /unit,color=col(ifl)
      oplot,rad_A(il,ind_samp),Ne_ucomp_A(il,ind_samp)/unit,color=col(ifl),psym=4
   endif
endfor
endif

loadct,0
color=0
if keyword_set(aia)   then oplot,rad_fit_aia_A  ,Ne_fit_aia_Avg  /unit,color=color,th=8
if keyword_set(kcor)  then oplot,rad_fit_kcor_A ,Ne_fit_kcor_Avg /unit,color=color,th=8,linestyle=2
if keyword_set(ucomp) then oplot,rad_fit_ucomp_A,Ne_fit_ucomp_Avg/unit,color=color,th=8,linestyle=3
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


