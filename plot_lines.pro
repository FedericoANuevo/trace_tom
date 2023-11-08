pro plot_lines, dir=dir, structure_filename=structure_filename, $
                aia = aia, euvia = euvia, euvib = euvib, eit = eit, $
                mk4 = mk4, kcor = kcor, lascoc2 = lascoc2

  common data, N_fl, Npt_max

  load_traced_data_structure, dir=dir, structure_filename=structure_filename, trace_data=trace_data  
 
; Overplot AIA Ne(r) for all lines, and their averaged profile.
  if keyword_set(aia) then begin
     ps1, dir+structure_filename+'_AIA-Ne.eps',0 
     tmp = reform( (*trace_data.INDEX_SAMPLING_AIA)(0,*) )
     ind = where(tmp eq 1)
     plot,[0,1],[0,1],/nodata,$
          xr = [1,1.25],xstyle=1,$    
          yr = [0,max((*trace_data.Ne_AIA))],ystyle=1,$
          charsize=1,font=0,$
          xtitle = 'r [R!Dsun!N]', title='AIA-DEMT Ne [cm!U-3!N]'
     for ifl=0,N_fl-1 do begin
        ind_samp = reform( (*trace_data.INDEX_SAMPLING_AIA)(ifl,*) )
        ind_samp = where(ind_samp eq 1)
        oplot,((*trace_data.rad)(ifl,*))(ind_samp),((*trace_data.Ne_AIA)(ifl,*))(ind_samp)
        interpol_fl,xv=((*trace_data.rad)(ifl,*))(ind_samp),yv=((*trace_data.Ne_AIA)(ifl,*))(ind_samp),xi=xi,yi=yi,/aia
        oplot,xi,yi
        if ifl eq 0 then yi_avg = yi
        if ifl gt 0 then yi_avg = yi_avg + yi
     endfor
     yi_avg = yi_avg / float(N_fl)
     loadct,12
     oplot,xi,yi_avg,th=8,color=100
     loadct,0
     ps2
  endif
  
  stop
  return
end
