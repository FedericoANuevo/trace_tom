;
; PURPOSE: This routine loads into memory a specified traced-data
; structure, and convenientely puts its contents back into arrays.
;
; INPUTS:
; dir and structure_filename: STRINGS; dir where the SAV file is
; located and its filename.
;
; OUTPUTS:
; trace_data: STRUCTURE; containing the tracing of tomographic
; products along field lines, as well as their geometry.
; several arrays and variables: all listed in COMMON BLOCK 'DATA'.
;
; FLAGS:
; /aia,.../lascoc2, set them to define in memory arrays with results
; from every instrument.
;
; HISTORY: V1.0, AMV, November 2023, IAFE.
;
pro load_traced_data_structure, dir=dir, structure_filename=structure_filename, trace_data=trace_data, $
                                aia = aia, euvia = euvia, euvib = euvib, eit = eit, $
                                mk4 = mk4, kcor = kcor, lascoc2 = lascoc2, $
                                fits = fits
  
  common data, N_fl, Npt_max, Npt_v, x_A, y_A, z_A, rad_A, lat_A, lon_A,$
     Ne_aia_A, Tm_aia_A, index_aia_A, index_sampling_aia_A,$
     Ne_euvia_A, Tm_euvia_A, index_euvia_A, index_sampling_euvia_A,$
     Ne_euvib_A, Tm_euvib_A, index_euvib_A, index_sampling_euvib_A,$
     Ne_eit_A, Tm_eit_A, index_eit_A, index_sampling_eit_A,$
     Ne_mk4_A, index_mk4_A, index_sampling_mk4_A,$
     Ne_kcor_A, index_kcor_A, index_sampling_kcor_A,$
     Ne_c2_A, index_c2_A, index_sampling_c2_A,$
     r_fit_aia_A, Ne_fit_aia_A, Te_fit_aia_A

  restore, filename = dir + structure_filename

  N_fl    = *trace_data.N_fl
  Npt_max = *trace_data.Npt_max
  Npt_v   = *trace_data.Npt_v
  x_A     = *trace_data.x
  y_A     = *trace_data.y
  z_A     = *trace_data.z
  rad_A   = *trace_data.rad
  lat_A   = *trace_data.lat
  lon_A   = *trace_data.lon
  if keyword_set(aia) then begin
                Ne_aia_A = *trace_data.Ne_aia 
                Tm_aia_A = *trace_data.Tm_aia
             index_aia_A = *trace_data.index_aia
    index_sampling_aia_A = *trace_data.index_sampling_aia
     if keyword_set(fits) then begin
             r_fit_aia_A = *trace_data.r_fit_aia             
            Ne_fit_aia_A = *trace_data.Ne_fit_aia
            Te_fit_aia_A = *trace_data.Tm_fit_aia
     endif
  endif
  if keyword_set(euvia) then begin
                Ne_euvia_A = *trace_data.Ne_euvia 
                Tm_euvia_A = *trace_data.Tm_euvia
             index_euvia_A = *trace_data.index_euvia
    index_sampling_euvia_A = *trace_data.index_sampling_euvia
  endif
  if keyword_set(euvib) then begin
                Ne_euvib_A = *trace_data.Ne_euvib 
                Tm_euvib_A = *trace_data.Tm_euvib
             index_euvib_A = *trace_data.index_euvib
    index_sampling_euvib_A = *trace_data.index_sampling_euvib
  endif
  if keyword_set(eit) then begin
                Ne_eit_A = *trace_data.Ne_eit 
                Tm_eit_A = *trace_data.Tm_eit
             index_eit_A = *trace_data.index_eit
    index_sampling_eit_A = *trace_data.index_sampling_eit
  endif
  if keyword_set(mk4) then begin
                Ne_mk4_A = *trace_data.Ne_mk4 
             index_mk4_A = *trace_data.index_mk4
    index_sampling_mk4_A = *trace_data.index_sampling_mk4
  endif
  if keyword_set(kcor) then begin
                Ne_kcor_A = *trace_data.Ne_kcor 
             index_kcor_A = *trace_data.index_kcor
    index_sampling_kcor_A = *trace_data.index_sampling_kcor
  endif
  if keyword_set(lascoc2) then begin
                Ne_c2_A = *trace_data.Ne_c2 
             index_c2_A = *trace_data.index_c2
    index_sampling_c2_A = *trace_data.index_sampling_c2
  endif
  return
end
