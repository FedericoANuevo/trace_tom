
pro Yeimy_check,aia=aia,lascoc2=lascoc2
  
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
  
 

  dir = '/data1/DATA/trace_tom_files/April24/field_lines_geometry_yeimy/'
  structure_filename='list_yeimy-fl.txt-tracing-structure-merge_aia_lascoc2.sav'
  load_traced_data_structure, dir=dir, structure_filename=structure_filename, trace_data=trace_data,/aia,/lascoc2


  rad_concat = [1.25]
  Ne_concat  = [1.e8]
  
  ps1,'check_fig5_Yeimy.eps'
  loadct,0
  plot,rad_concat-1,Ne_concat,/xlog,/ylog,/nodata,title='AIA and C2 Ne',xtitle='rad [Rsun]-1',ytitle='Ne [cm!U-3!N]',xr=[1.e-2,1.e2],yr=[1.e1,1.e9]
  loadct,39
  
  if keyword_set(lascoc2) then begin
     xpos = 1.0 - 0.2
     ypos = 0.8- 0.05
     ilstep = 1
     color = 100
     xyouts,[xpos],[ypos],['C2'],color=color,/normal,charsize = 2.,charthick = 2
     for ifl=0,N_fl-1,ilstep do begin
        tmp      = reform(index_sampling_c2_A(ifl,*))
        ind_samp = where(tmp eq 1)
        oplot,rad_A(ifl,ind_samp)-1,Ne_c2_A(ifl,ind_samp),psym=4,color=color
     endfor
  endif
  
  if keyword_set(aia) then begin
     ilstep = 1
     color = 200
     xpos = 1.0 -0.2
     ypos = 0.8
     xyouts,[xpos],[ypos],['AIA'],color=color,/normal,charsize = 2.,charthick = 2
     for ifl=0,N_fl-1,ilstep do begin 
        tmp      = reform(index_sampling_aia_A(ifl,*))
        ind_samp = where(tmp eq 1)
        oplot,rad_A(ifl,ind_samp)-1,Ne_aia_A(ifl,ind_samp),psym=4,color = color
     endfor
  endif


  ps2
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


