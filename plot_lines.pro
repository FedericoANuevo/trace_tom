pro plot_lines, dir=dir, structure_filename=structure_filename, $
                aia = aia, euvia = euvia, euvib = euvib, eit = eit, $
                mk4 = mk4, kcor = kcor, lascoc2 = lascoc2

  common data, N_fl, Npt_max, x_A, y_A, z_A, rad_A, lat_A, lon_A,$
     Ne_aia_A, Tm_aia_A, index_aia_A, index_sampling_aia_A,$
     Ne_euvia_A, Tm_euvia_A, index_euvia_A, index_sampling_euvia_A,$
     Ne_euvib_A, Tm_euvib_A, index_euvib_A, index_sampling_euvib_A,$
     Ne_eit_A, Tm_eit_A, index_eit_A, index_sampling_eit_A,$
     Ne_mk4_A, index_mk4_A, index_sampling_mk4_A,$
     Ne_kcor_A, index_kcor_A, index_sampling_kcor_A,$
     Ne_lascoc2_A, index_lascoc2_A, index_sampling_lascoc2_A

  load_traced_data_structure, dir=dir, structure_filename=structure_filename, trace_data=trace_data,$
                              /aia
 
; Overplot AIA Ne(r) for all lines, and their averaged profile.
  if keyword_set(aia) then begin
     ps1, dir+structure_filename+'_AIA-Ne.eps',0 
     tmp = INDEX_SAMPLING_AIA_A(0,*)
     ind = where(tmp eq 1)
     plot,[0,1],[0,1],/nodata,$
          xr = [1,1.25],xstyle=1,$    
          yr = [0,max(Ne_AIA_A)],ystyle=1,$
          charsize=1,font=0,$
          xtitle = 'r [R!Dsun!N]', title='AIA-DEMT Ne [cm!U-3!N]'
     for ifl=0,N_fl-1 do begin
        ind_samp = reform( INDEX_SAMPLING_AIA_A(ifl,*) )
        ind_samp = where(ind_samp eq 1)
        oplot,(rad_A(ifl,*))(ind_samp),(Ne_AIA_A(ifl,*))(ind_samp)
        interpol_fl,xv=(rad_A(ifl,*))(ind_samp),yv=(Ne_AIA_A(ifl,*))(ind_samp),xi=xi,yi=yi,/aia
        oplot,xi,yi
        if ifl eq 0 then yi_avg =          yi/float(N_fl)
        if ifl gt 0 then yi_avg = yi_avg + yi/float(N_fl)
     endfor
     loadct,12
     oplot,xi,yi_avg,th=8,color=100
     loadct,0
     ps2
  endif
  stop
  return
end
