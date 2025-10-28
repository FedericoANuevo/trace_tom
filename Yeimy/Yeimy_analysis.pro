pro Yeimy_analysis, load=load, LonLimits=LonLimits, LatLimits=LatLimits,$
                    plot_filename_suffix=plot_filename_suffix,$
                    aia=aia,lascoc2=lascoc2,$
                    open=open,closed=closed,connect=connect,$
                    not_Bfield=not_Bfield

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
  
 
  if keyword_set(not_Bfield) then begin
     dir = '/data1/DATA/trace_tom_files/April24/field_lines_geometry_yeimy/'
     structure_filename='list_yeimy-fl.txt-tracing-structure-merge_aia_lascoc2.sav'
  endif else begin
   ; dir = '/data1/DATA/trace_tom_files/April24/field_lines_geometry/'
     dir = '/data1/DATA/trace_tom_files/April24/field_lines_geometry_1x1deg_multirad/'
     structure_filename='fdips_field_150X180X360_mrbqs240414t1304c2283_230_shift.ubdat_fline-filenames_list.txt-tracing-structure-merge_aia.sav'
  endelse
  if keyword_set(load) then begin
     if not keyword_set(not_Bfield) then load_traced_data_structure, dir=dir, structure_filename=structure_filename, trace_data=trace_data,/aia,/trace_Bs
     if     keyword_set(not_Bfield) then load_traced_data_structure, dir=dir, structure_filename=structure_filename, trace_data=trace_data,/aia,/lascoc2
  endif
; Define plot filename suffix
  if not keyword_set(plot_filename_suffix) then plot_filename_suffix='footpoints-map'

; Define box of footpoint Lat/Lot to analyze
  if not keyword_set(LonLimits) then LonLimits = [  0., 360.]
  if not keyword_set(LatLimits) then LatLimits = [-90.,+ 90.]

  if     keyword_set(not_Bfield) then       leg_label_A = fltarr(N_fl)
  if not keyword_set(not_Bfield) then begin
; Define a polarity tag
     tag_polarity_A = reform(leg_footBfield_A(*,0)/abs(leg_footBfield_A(*,0)))
     if ( where(finite(tag_polarity_A) eq 0) )(0) ne -1 then STOP ; assumes foot-Br is not ZERO
  endif else begin
     tag_polarity_A = leg_label_A ; Este es un truco para el caso de las f-lines provistas por Yeimy
  endelse
  
; Tag field lines with footpoints within the BOX., either open or closed.
  tag_box_A = fltarr(N_fl)
  ifl_A = where( (Footpoint_Lon_A ge LonLimits(0) AND Footpoint_Lon_A le LonLimits(1)) AND $
                 (Footpoint_Lat_A ge LatLimits(0) AND Footpoint_Lat_A le LatLimits(1)) )
  tag_box_A(ifl_A) = +1.

; Now, UNTAG CLOSED field lines whose other leg's footpoint is NOT wihin the BOX,
;      as well as CLOSED field lines that do NOT comply with opposite polarity.
   ifl=0
    while ifl le N_fl-2 do begin
       if leg_label_A(ifl) eq 0. then begin
          ifl=ifl+1
       endif else begin
          if leg_label_A(ifl) ne leg_label_A(ifl+1) then STOP ; !this should never happen.
          if (tag_box_A(ifl) eq +1 AND tag_box_A(ifl+1) eq  0) then tag_box_A(ifl  )=0
          if (tag_box_A(ifl) eq  0 AND tag_box_A(ifl+1) eq +1) then tag_box_A(ifl+1)=0
          if not keyword_set(not_bfield) then begin
             if (tag_polarity_A(ifl) eq tag_polarity_A(ifl+1))    then tag_box_A(ifl:ifl+1)=0
          endif
          ifl=ifl+2
       endelse
    endwhile
    
; Create index arrays for closed and open field lines.
  ifl_closed_A = where(leg_label_A ne 0)
  ifl_open_A   = where(leg_label_A eq 0)

; Create index arrays for positive/negative footpoint Brad lines
  ifl_pos_A = where(tag_polarity_A eq +1.)
  ifl_neg_A = where(tag_polarity_A eq -1.)
  ifl_zer_A = where(tag_polarity_A eq  0.)

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

if keyword_set(aia  ) then begin
   r_crit = median(rad_fit_aia_A  )
   r_max  = max   (rad_fit_aia_A  )
endif

if keyword_set(lascoc2) then begin
   r_crit = median(rad_fit_c2_A  )
   r_max  = max   (rad_fit_c2_A  )
endif


;
tag_pos_aia_A          = fltarr(N_fl) - 678.
tag_fullrange_aia_A    = fltarr(N_fl) - 678.
tag_pos_c2_A           = fltarr(N_fl) - 678.
tag_fullrange_c2_A     = fltarr(N_fl) - 678.

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
   if keyword_set(lascoc2) then begin
      ifitpos_c2 = where(reform(Ne_fit_c2_A  (ifl,*)) gt 0. and rad_fit_c2_A le R_max)
      ifitrad_c2 = where(rad_fit_c2_A le R_max) 
      if ifitpos_c2(0) ne -1 then begin
         if n_elements(ifitpos_c2) eq n_elements(ifitrad_c2) then tag_pos_c2_A(ifl) = +1.
         if n_elements(where(rad_fit_c2_A(ifitpos_c2) lt r_crit))                                        ge Nsamp AND  $
            n_elements(where(rad_fit_c2_A(ifitpos_c2) gt r_crit AND rad_fit_c2_A(ifitpos_c2) lt r_max))  ge Nsamp then tag_fullrange_c2_A(ifl) = +1.
      endif
   endif else begin
      tag_pos_c2_A(ifl) = +1
      tag_fullrange_c2_A(ifl) = +1
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
if keyword_set(aia)     then !p.multi=[0,1,3]
if keyword_set(lascoc2) then !p.multi=[0,1,2]
loadct,0
!p.color=0
!p.background=255

csz=1
                            opcl_str='open/closed'
if keyword_set(open)   then opcl_str='open '
if keyword_set(closed) then opcl_str='closed '
plot,lon_A,lat_A,xr=[0,360],yr=[-90,+90],xstyle=1,ystyle=1,/nodata,charsize=csz,$
     title='Location of '+opcl_str+'footpoints',ytitle='Carrington Latitude [deg]',font=0
oplot,Footpoint_Lon_A,Footpoint_Lat_A,psym=4

loadct,12

; Set a maximum threshold for scN and scT
scN_crit = 0.25
scT_crit = 0.25
; Index lines of each instrument that have a good fit to Ne and are in
; the box
if keyword_set(aia)    then ifl_aia_A  = where(tag_box_A eq +1 AND tag_pos_aia_A   eq +1 AND fitflag_aia_A   eq +1 AND scN_fit_aia_A   le scN_crit and scT_fit_aia_A le scT_crit)
if keyword_set(lascoc2)then ifl_c2_A   = where(tag_box_A eq +1 AND tag_pos_c2_A   eq +1 AND fitflag_c2_A   eq +1 AND scN_fit_c2_A   le scN_crit)

; Make an index array with the common elements of all ifl_INSTRUMENT_A
ifl_A = indgen(N_fl)
if keyword_set(aia)       then ifl_A = intersect(ifl_A,ifl_aia_A  )
if keyword_set(lascoc2)   then ifl_A = intersect(ifl_A,ifl_c2_A  )
                           
; Intersect with CLOSED or OPEN, if so requested
if keyword_set(closed) then ifl_A = intersect(ifl_A,ifl_closed_A)
if keyword_set(open)   then ifl_A = intersect(ifl_A,ifl_open_A  )

; Plot connectivity is requested
if keyword_set(connect) AND keyword_set(closed) then begin
   ifl_A_orig = ifl_A 
   ifl=0
   while ifl le N_fl-2 do begin
      if leg_label_A(ifl) ne 0 then begin
         if tag_box_A(ifl) eq +1 then oplot,reform(Footpoint_Lon_A(ifl:ifl+1)),reform(Footpoint_Lat_A(ifl:ifl+1)),color=blue
         ifl=ifl+2
      endif else begin
         ifl=ifl+1
      endelse
   endwhile
endif

; Now, color-highlight the footpoints indicated by ifl_A
indxpos_A = intersect(ifl_A,ifl_pos_A)
if n_elements(indxpos_A) gt 1 then $
   oplot,Footpoint_Lon_A(indxpos_A),Footpoint_Lat_A(indxpos_A),psym=4,th=2,color=green
indxneg_A = intersect(ifl_A,ifl_neg_A)
if n_elements(indxneg_A) gt 1 then $
   oplot,Footpoint_Lon_A(indxneg_A),Footpoint_Lat_A(indxneg_A),psym=4,th=2,color=red
indxzer_A = intersect(ifl_A,ifl_zer_A)
if n_elements(indxzer_A) gt 1 then $
   oplot,Footpoint_Lon_A(indxzer_A),Footpoint_Lat_A(indxzer_A),psym=4,th=2,color=blue

if NOT keyword_set(closed) and NOT keyword_set(open) then begin
   loadct,0
   oplot,Footpoint_Lon_A(ifl_A),Footpoint_Lat_A(ifl_A),psym=4,th=2
endif

; Compute the average Ne(r) and Tm(r) of each instrument for lines indexed ifl_A
if keyword_set(aia) then begin
   Ne_fit_aia_Avg = total( Ne_fit_aia_A(ifl_A,*) , 1 ) / float(n_elements(ifl_A))
   Tm_fit_aia_Avg = total( Tm_fit_aia_A(ifl_A,*) , 1 ) / float(n_elements(ifl_A))
endif
if keyword_set(lascoc2) then begin
   Ne_fit_c2_Avg = total( Ne_fit_c2_A(ifl_A,*) , 1 ) / float(n_elements(ifl_A))
endif


skip_tag_pos:

; Plot the average Ne(r) for all selected instruments.
if keyword_set(aia) then begin
   loadct,0
   unit           = 1.e8        ; cm-3
   unit_power_str =   '8'
   xrange = [1.02,1.25]
   yrflag = -1
   r    = rad_fit_aia_A
   f    = Ne_fit_aia_Avg/unit
   miny = min(f(where(f gt 0. and r le max(xrange))))
   maxy = max(f(where(f gt 0. and r le max(xrange))))
   yrange = [miny*0.5,maxy*2.0]
   yrflag = +1
   plot,rad_fit_AiA_A,Ne_fit_AIA_A(0,*)/unit,charsize=csz,font=0,$
        title='<N!De!N(r)>   Solid: AIA',$
        ytitle = 'Ne(r) [x 10!U'+unit_power_str+'!N cm!U-3!N]',yr=yrange,ystyle=1,$
        xtitle = 'r [Rsun]'                                   ,xr=xrange,xstyle=1,$
        /nodata  
  ;loadct,12
   color=0
endif
if keyword_set(lascoc2) then begin
   loadct,0
   unit           = 1.e5        ; cm-3
   unit_power_str =   '5'
   xrange = [2.5,6.0]
   yrflag = -1
   r    = rad_fit_c2_A
   f    = Ne_fit_c2_Avg/unit
   miny = min(f(where(f gt 0. and r le max(xrange))))
   maxy = max(f(where(f gt 0. and r le max(xrange))))
   yrange = [miny,maxy]
   yrflag = +1
   plot,rad_fit_c2_A,Ne_fit_c2_A(0,*)/unit,charsize=csz,font=0,$
        title='<N!De!N(r)>   Solid: C2',$
        ytitle = 'Ne(r) [x 10!U'+unit_power_str+'!N cm!U-3!N]',yr=yrange,ystyle=1,$
        xtitle = 'r [Rsun]'                                   ,xr=xrange,xstyle=1,$
        /nodata  
  ;loadct,12
   color=0
endif

loadct,39
if keyword_set(lascoc2) then begin
   Nifl = n_elements(ifl_A)
   col = 255 * findgen(Nifl)/float(Nifl-1)
   ilstep = 1
   for ifl=0,Nifl-1,ilstep do begin
      il = (ifl_A(ifl))(0)
      tmp      = reform(index_sampling_c2_A(il,*))
      ind_samp = where(tmp eq 1)
      oplot,rad_fit_c2_A     ,Ne_fit_c2_A(il,*)   /unit,color=col(ifl)
      oplot,rad_A(il,ind_samp),Ne_c2_A(il,ind_samp)/unit,color=col(ifl),psym=4
   endfor
endif
if keyword_set(aia) then begin
   Nifl = n_elements(ifl_A)
   col = 255 * findgen(Nifl)/float(Nifl-1)
   ilstep = 1
   for ifl=0,Nifl-1,ilstep do begin
      il = (ifl_A(ifl))(0)
      tmp      = reform(index_sampling_aia_A(il,*))
      ind_samp = where(tmp eq 1)
      oplot,rad_fit_aia_A     ,Ne_fit_aia_A(il,*)   /unit,color=col(ifl)
      oplot,rad_A(il,ind_samp),Ne_aia_A(il,ind_samp)/unit,color=col(ifl),psym=4
   endfor
endif




loadct,0
color = 0
if keyword_set(aia)     and n_elements(ifl_A) gt 1  then oplot,rad_fit_aia_A  ,Ne_fit_aia_Avg  /unit,color=color,th=8
if keyword_set(lascoc2) and n_elements(ifl_A) gt 1  then oplot,rad_fit_c2_A  ,Ne_fit_c2_Avg    /unit,color=color,th=8
loadct,0

if keyword_set(aia) then begin
yrange = [min(Tm_fit_aia_Avg)/1.e6*0.5,max(Tm_fit_aia_Avg)/1.e6*2.]
plot,rad_fit_AiA_A,Tm_fit_AIA_A(0,*)/unit,charsize=csz,font=0,$
     title='<T!Dm!N(r)>   Solid: AIA',$
     ytitle = 'Tm(r) [MK]',yr=yrange,ystyle=1,$
     xtitle = 'r [Rsun]'                                   ,xr=xrange,xstyle=1,$
     /nodata
;loadct,12
color=0
endif

loadct,39
if keyword_set(aia) then begin
   Nifl = n_elements(ifl_A)
   col = 255 * findgen(Nifl)/float(Nifl-1)
   ilstep = 1
   for ifl=0,Nifl-1,ilstep do begin
      il = (ifl_A(ifl))(0)
      tmp      = reform(index_sampling_aia_A(il,*))
      ind_samp = where(tmp eq 1)
      oplot,rad_fit_aia_A     ,Tm_fit_aia_A(il,*)   /1.e6,color=col(ifl)
      oplot,rad_A(il,ind_samp),Tm_aia_A(il,ind_samp)/1.e6,color=col(ifl),psym=4
   endfor
endif

loadct,0
color =0
if keyword_set(aia) and  n_elements(ifl_A) gt 1 then oplot,rad_fit_aia_A  ,Tm_fit_aia_Avg/1.e6,color=color,th=6
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


