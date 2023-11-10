pro plot_lines, dir=dir, structure_filename=structure_filename, $
                aia = aia, euvia = euvia, euvib = euvib, eit = eit, $
                mk4 = mk4, kcor = kcor, lascoc2 = lascoc2, $
                log = log

  common data, N_fl, Npt_max, x_A, y_A, z_A, rad_A, lat_A, lon_A,$
     Ne_aia_A, Tm_aia_A, index_aia_A, index_sampling_aia_A,$
     Ne_euvia_A, Tm_euvia_A, index_euvia_A, index_sampling_euvia_A,$
     Ne_euvib_A, Tm_euvib_A, index_euvib_A, index_sampling_euvib_A,$
     Ne_eit_A, Tm_eit_A, index_eit_A, index_sampling_eit_A,$
     Ne_mk4_A, index_mk4_A, index_sampling_mk4_A,$
     Ne_kcor_A, index_kcor_A, index_sampling_kcor_A,$
     Ne_c2_A, index_c2_A, index_sampling_c2_A

  load_traced_data_structure, dir=dir, structure_filename=structure_filename, trace_data=trace_data,$
                              aia = aia, euvia = euvia, euvib = euvib, eit = eit, $
                              mk4 = mk4, kcor = kcor, lascoc2 = lascoc2

; Set optimal-height-range and optimal-Ne-range for requested data.
     ranger = 0
     rangeN = 0
     if keyword_set(aia)     then begin
        ranger = [1.0,1.25]
        ipos_aia = where(Ne_aia_A gt 0.)
        rangeN_aia = [min(Ne_aia_A(ipos_aia)),max(Ne_aia_A(ipos_aia))]
        rangeN     = rangeN_aia 
     endif
     if keyword_set(mk4)     then begin
        ranger = [1.15,1.5]
        ipos_mk4 = where(Ne_mk4_A gt 0.)
        rangeN_mk4 = [min(Ne_mk4_A(ipos_mk4)),max(Ne_mk4_A(ipos_mk4))]
        if NOT keyword_set(rangeN) then rangeN = rangeN_mk4 else $
        if     keyword_set(rangeN) then rangeN = [min([rangeN(0),rangeN_mk4(0)]),max([rangeN(1),rangeN_mk4(1)])]
     endif
     if keyword_set(lascoc2) then begin
        if NOT keyword_set(ranger) then ranger = [2.5,6.0] else $
        if     keyword_set(ranger) then ranger = [1.0,6.0] else stop
        ipos_c2 = where(Ne_c2_A gt 0.)
        rangeN_c2 = [min(Ne_c2_A(ipos_c2)),max(Ne_c2_A(ipos_c2))]
        if NOT keyword_set(rangeN) then rangeN = rangeN_c2 else $
        if     keyword_set(rangeN) then rangeN = [min([rangeN(0),rangeN_c2(0)]),max([rangeN(1),rangeN_c2(1)])]
     endif

; Ne(r) for all lines and instruments, and their averaged profile.
     ps1, dir+structure_filename+'_Ne(r).eps',0
     if NOT keyword_set(log) then $
     plot,[0,1],[0,1],/nodata ,xr = ranger,xstyle=1,yr = rangeN,ystyle=1,charsize=1,font=0,$
          xtitle = 'r [R!Dsun!N]', title='Ne [cm!U-3!N]'
     if     keyword_set(log) then $
     plot,[0,1],[0,1],/nodata ,xr = ranger,xstyle=1,yr = rangeN,ystyle=1,charsize=1,font=0,$
          xtitle = 'r [R!Dsun!N]', title='Ne [cm!U-3!N]',/ylog
     
     for ifl=0,N_fl-1 do begin
        if keyword_set(aia) then begin
           tmp = reform(index_sampling_aia_A(ifl,*)) & ind_samp = where(tmp eq 1)
           oplot,(rad_A(ifl,*))(ind_samp),(Ne_AIA_A(ifl,*))(ind_samp)
           interpol_fl,xv=(rad_A(ifl,*))(ind_samp),yv=(Ne_AIA_A(ifl,*))(ind_samp),xi=xi_aia,yi=yi,/aia
           oplot,xi_aia,yi
           if ifl eq 0 then yi_aia_avg =              yi/float(N_fl)
           if ifl gt 0 then yi_aia_avg = yi_aia_avg + yi/float(N_fl)
        endif
        if keyword_set(mk4) then begin
           tmp = reform(index_sampling_mk4_A(ifl,*)) & ind_samp = where(tmp eq 1)
           oplot,(rad_A(ifl,*))(ind_samp),(Ne_mk4_A(ifl,*))(ind_samp)
           interpol_fl,xv=(rad_A(ifl,*))(ind_samp),yv=(Ne_mk4_A(ifl,*))(ind_samp),xi=xi_mk4,yi=yi,/mk4
           oplot,xi_mk4,yi
           if ifl eq 0 then yi_mk4_avg =              yi/float(N_fl)
           if ifl gt 0 then yi_mk4_avg = yi_mk4_avg + yi/float(N_fl)
        endif
        if keyword_set(lascoc2) then begin
           tmp = reform(index_sampling_c2_A(ifl,*)) & ind_samp = where(tmp eq 1)
           oplot,(rad_A(ifl,*))(ind_samp),(Ne_c2_A(ifl,*))(ind_samp)
           interpol_fl,xv=(rad_A(ifl,*))(ind_samp),yv=(Ne_c2_A(ifl,*))(ind_samp),xi=xi_c2,yi=yi,/lascoc2
           oplot,xi_c2,yi
           if ifl eq 0 then yi_c2_avg =              yi/float(N_fl)
           if ifl gt 0 then yi_c2_avg = yi_c2_avg + yi/float(N_fl)
        endif
     endfor
     loadct,12
     if keyword_set(aia)     then oplot,xi_aia,yi_aia_avg,th=8,color=100
     if keyword_set(mk4)     then oplot,xi_mk4,yi_mk4_avg,th=8,color=100
     if keyword_set(lascoc2) then oplot,xi_c2 ,yi_c2_avg ,th=8,color=100
     loadct,0
     ps2
     
  stop
  return
end
