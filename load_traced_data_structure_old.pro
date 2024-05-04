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
;          V1.1, AMV, January  2023, IAFE. Added fitted results.
;          V1.2, AMV, January  2023, IAFE. Added fits' parameters.
;
pro load_traced_data_structure, dir=dir, structure_filename=structure_filename, trace_data=trace_data, $
                                aia = aia, euvia = euvia, euvib = euvib, eit = eit, $
                                mk4 = mk4, kcor = kcor, lascoc2 = lascoc2
  
  common data, N_fl, Npt_max, Npt_v, x_A, y_A, z_A, rad_A, lat_A, lon_A,$
     Ne_aia_A, Tm_aia_A, index_aia_A, index_sampling_aia_A,$
     Ne_euvia_A, Tm_euvia_A, index_euvia_A, index_sampling_euvia_A,$
     Ne_euvib_A, Tm_euvib_A, index_euvib_A, index_sampling_euvib_A,$
     Ne_eit_A, Tm_eit_A, index_eit_A, index_sampling_eit_A,$
     Ne_mk4_A, index_mk4_A, index_sampling_mk4_A,$
     Ne_kcor_A, index_kcor_A, index_sampling_kcor_A,$
     Ne_c2_A, index_c2_A, index_sampling_c2_A,$
     rad_fit_aia_A, Ne_fit_aia_A, Tm_fit_aia_A, fitflag_aia_A,scN_fit_aia_A,scT_fit_aia_A,$
     rad_fit_c2_A, Ne_fit_c2_A, fitflag_c2_A,scN_fit_c2_A,$
     rad_fit_mk4_A, Ne_fit_mk4_A, fitflag_mk4_A,scN_fit_mk4_A,$
     N0_fit_aia_A,lN_fit_aia_A,T0_fit_aia_A,dTdr_fit_aia_A,$
     N1_fit_aia_A,N2_fit_aia_A,p1_fit_aia_A,p2_fit_aia_A,$
     N0_fit_mk4_A,lN_fit_mk4_A,$
     N1_fit_mk4_A,N2_fit_mk4_A,p1_fit_mk4_A,p2_fit_mk4_A,$
     N1_fit_c2_A,N2_fit_c2_A,p1_fit_c2_A,p2_fit_c2_A,$
     lN_fit_c2_A,$
     fit_F_Ne_aia,fit_F_Ne_mk4,fit_F_Ne_c2

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
           rad_fit_aia_A = *trace_data. rad_fit_aia             
            Ne_fit_aia_A = *trace_data.  Ne_fit_aia
            Tm_fit_aia_A = *trace_data.  Tm_fit_aia
           scN_fit_aia_A = *trace_data. scN_fit_aia
           scT_fit_aia_A = *trace_data. scT_fit_aia
           fitflag_aia_A = *trace_data. fitflag_aia
            T0_fit_aia_A = *trace_data.  T0_fit_aia
          dTdr_fit_aia_A = *trace_data.dTdr_fit_aia
          fit_F_Ne_aia   = *trace_data.fit_F_Ne_aia
            lN_fit_aia_A = *trace_data.  lN_fit_aia
          if fit_F_Ne_aia eq 'IHS' then $
            N0_fit_aia_A = *trace_data.  N0_fit_aia
          if fit_F_Ne_aia eq 'DPL' then begin
             N1_fit_aia_A = *trace_data. N1_fit_aia
             N2_fit_aia_A = *trace_data. N2_fit_aia
             p1_fit_aia_A = *trace_data. p1_fit_aia
             p2_fit_aia_A = *trace_data. p2_fit_aia
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
           rad_fit_mk4_A = *trace_data. rad_fit_mk4             
            Ne_fit_mk4_A = *trace_data.  Ne_fit_mk4
           scN_fit_mk4_A = *trace_data. scN_fit_mk4
           fitflag_mk4_A = *trace_data. fitflag_mk4
          fit_F_Ne_mk4   = *trace_data.fit_F_Ne_mk4
            lN_fit_mk4_A = *trace_data.  lN_fit_mk4 
          if fit_F_Ne_mk4 eq 'IHS' then $
             N0_fit_mk4_A = *trace_data. N0_fit_mk4
          if fit_F_Ne_mk4 eq 'SPL' then begin
             N1_fit_mk4_A = *trace_data. N1_fit_mk4
             p1_fit_mk4_A = *trace_data. p1_fit_mk4
          endif
          if fit_F_Ne_mk4 eq 'DPL' then begin
             N1_fit_mk4_A = *trace_data. N1_fit_mk4
             N2_fit_mk4_A = *trace_data. N2_fit_mk4
             p1_fit_mk4_A = *trace_data. p1_fit_mk4
             p2_fit_mk4_A = *trace_data. p2_fit_mk4
          endif
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
           rad_fit_c2_A = *trace_data. rad_fit_c2             
            Ne_fit_c2_A = *trace_data.  Ne_fit_c2
           scN_fit_c2_A = *trace_data. scN_fit_c2
           fitflag_c2_A = *trace_data. fitflag_c2
          fit_F_Ne_c2   = *trace_data.fit_F_Ne_c2
            lN_fit_c2_A = *trace_data.  lN_fit_c2 
          if fit_F_Ne_c2 eq 'SPL' then begin
             N1_fit_c2_A = *trace_data. N1_fit_c2
             p1_fit_c2_A = *trace_data. p1_fit_c2
          endif
          if fit_F_Ne_c2 eq 'DPL' then begin
             N1_fit_c2_A = *trace_data. N1_fit_c2
             N2_fit_c2_A = *trace_data. N2_fit_c2
             p1_fit_c2_A = *trace_data. p1_fit_c2
             p2_fit_c2_A = *trace_data. p2_fit_c2
          endif
  endif
  return
end
