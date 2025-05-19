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
; /aia,.../lascoc2, set one flag for each specific instrument whose
; results is to be loaded into memory. The available instruments for a
; particular SAV file are explicitly indicated as suffixes of its
; filename.
;
; HISTORY: V1.0, AMV, November 2023, IAFE.
;          V1.1, AMV, January  2024, IAFE. Added fitted results.
;          V1.2, AMV, January  2024, IAFE. Added fits' parameters.
;          V1.3, FAN, May      2024, ClaSP. Added euvia, euvib, eit.
;          V1.4, AMV, May      2025, CLaSP. Expanded to ucomp, 
;                                    opcl, tom and fit grids.
;
pro load_traced_data_structure, dir=dir, structure_filename=structure_filename, trace_data=trace_data, $
                                aia=aia, euvia=euvia, euvib=euvib, eit=eit, $
                                mk4=mk4, kcor=kcor, lascoc2=lascoc2, ucomp=ucomp, $
                                opcl=opcl
  
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

  if keyword_set(opcl) then leg_label_A = *trace_data.leg_label
  
  if keyword_set(aia) then begin
                Ne_aia_A = *trace_data.Ne_aia 
                Tm_aia_A = *trace_data.Tm_aia
                WT_aia_A = *trace_data.WT_aia
         ldem_flag_aia_A = *trace_data.ldem_flag_aia
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
          if keyword_set(opcl) then begin
             tomgrid_aia_hdr_A = *trace_data. tomgrid_aia_hdr
             tomgrid_aia_A     = *trace_data. tomgrid_aia
             fitgrid_aia_hdr_A = *trace_data. fitgrid_aia_hdr
             fitgrid_aia_A     = *trace_data. fitgrid_aia
          endif          
  endif

  if keyword_set(euvia) then begin
                Ne_euvia_A = *trace_data.Ne_euvia 
                Tm_euvia_A = *trace_data.Tm_euvia
                WT_euvia_A = *trace_data.WT_euvia
         ldem_flag_euvia_A = *trace_data.ldem_flag_euvia
             index_euvia_A = *trace_data.index_euvia
    index_sampling_euvia_A = *trace_data.index_sampling_euvia
           rad_fit_euvia_A = *trace_data. rad_fit_euvia             
            Ne_fit_euvia_A = *trace_data.  Ne_fit_euvia
            Tm_fit_euvia_A = *trace_data.  Tm_fit_euvia
           scN_fit_euvia_A = *trace_data. scN_fit_euvia
           scT_fit_euvia_A = *trace_data. scT_fit_euvia
           fitflag_euvia_A = *trace_data. fitflag_euvia
            T0_fit_euvia_A = *trace_data.  T0_fit_euvia
          dTdr_fit_euvia_A = *trace_data.dTdr_fit_euvia
          fit_F_Ne_euvia   = *trace_data.fit_F_Ne_euvia
            lN_fit_euvia_A = *trace_data.  lN_fit_euvia
          if fit_F_Ne_euvia eq 'IHS' then $
            N0_fit_euvia_A = *trace_data.  N0_fit_euvia
          if fit_F_Ne_euvia eq 'DPL' then begin
             N1_fit_euvia_A = *trace_data. N1_fit_euvia
             N2_fit_euvia_A = *trace_data. N2_fit_euvia
             p1_fit_euvia_A = *trace_data. p1_fit_euvia
             p2_fit_euvia_A = *trace_data. p2_fit_euvia
          endif
          if keyword_set(opcl) then begin
             tomgrid_euvia_hdr_A = *trace_data. tomgrid_euvia_hdr
             tomgrid_euvia_A     = *trace_data. tomgrid_euvia
             fitgrid_euvia_hdr_A = *trace_data. fitgrid_euvia_hdr
             fitgrid_euvia_A     = *trace_data. fitgrid_euvia
          endif
       endif

  if keyword_set(euvib) then begin
                Ne_euvib_A = *trace_data.Ne_euvib 
                Tm_euvib_A = *trace_data.Tm_euvib
                WT_euvib_A = *trace_data.WT_euvib
         ldem_flag_euvib_A = *trace_data.ldem_flag_euvib
             index_euvib_A = *trace_data.index_euvib
    index_sampling_euvib_A = *trace_data.index_sampling_euvib
           rad_fit_euvib_A = *trace_data. rad_fit_euvib
            Ne_fit_euvib_A = *trace_data.  Ne_fit_euvib
            Tm_fit_euvib_A = *trace_data.  Tm_fit_euvib
           scN_fit_euvib_A = *trace_data. scN_fit_euvib
           scT_fit_euvib_A = *trace_data. scT_fit_euvib
           fitflag_euvib_A = *trace_data. fitflag_euvib
            T0_fit_euvib_A = *trace_data.  T0_fit_euvib
          dTdr_fit_euvib_A = *trace_data.dTdr_fit_euvib
          fit_F_Ne_euvib   = *trace_data.fit_F_Ne_euvib
            lN_fit_euvib_A = *trace_data.  lN_fit_euvib
          if fit_F_Ne_euvib eq 'IHS' then $
            N0_fit_euvib_A = *trace_data.  N0_fit_euvib
          if fit_F_Ne_euvia eq 'DPL' then begin
             N1_fit_euvib_A = *trace_data. N1_fit_euvib
             N2_fit_euvib_A = *trace_data. N2_fit_euvi
             p1_fit_euvib_A = *trace_data. p1_fit_euvib
             p2_fit_euvib_A = *trace_data. p2_fit_euvib
          endif
          if keyword_set(opcl) then begin
             tomgrid_euvib_hdr_A = *trace_data. tomgrid_euvib_hdr
             tomgrid_euvib_A     = *trace_data. tomgrid_euvib
             fitgrid_euvib_hdr_A = *trace_data. fitgrid_euvib_hdr
             fitgrid_euvib_A     = *trace_data. fitgrid_euvib
          endif
       endif
  

  if keyword_set(eit) then begin
                Ne_eit_A = *trace_data.Ne_eit 
                Tm_eit_A = *trace_data.Tm_eit
                WT_eit_A = *trace_data.WT_eit
         ldem_flag_eit_A = *trace_data.ldem_flag_eit
             index_eit_A = *trace_data.index_eit
    index_sampling_eit_A = *trace_data.index_sampling_eit
           rad_fit_eit_A = *trace_data. rad_fit_eit             
            Ne_fit_eit_A = *trace_data.  Ne_fit_eit
            Tm_fit_eit_A = *trace_data.  Tm_fit_eit
           scN_fit_eit_A = *trace_data. scN_fit_eit
           scT_fit_eit_A = *trace_data. scT_fit_eit
           fitflag_eit_A = *trace_data. fitflag_eit
            T0_fit_eit_A = *trace_data.  T0_fit_eit
          dTdr_fit_eit_A = *trace_data.dTdr_fit_eit
          fit_F_Ne_eit   = *trace_data.fit_F_Ne_eit
            lN_fit_eit_A = *trace_data.  lN_fit_eit
          if fit_F_Ne_eit eq 'IHS' then $
            N0_fit_eit_A = *trace_data.  N0_fit_eit
          if fit_F_Ne_eit eq 'DPL' then begin
             N1_fit_eit_A = *trace_data. N1_fit_eit
             N2_fit_eit_A = *trace_data. N2_fit_eit
             p1_fit_eit_A = *trace_data. p1_fit_eit
             p2_fit_eit_A = *trace_data. p2_fit_eit
          endif
          if keyword_set(opcl) then begin
             tomgrid_eit_hdr_A = *trace_data. tomgrid_eit_hdr
             tomgrid_eit_A     = *trace_data. tomgrid_eit
             fitgrid_eit_hdr_A = *trace_data. fitgrid_eit_hdr
             fitgrid_eit_A     = *trace_data. fitgrid_eit
          endif
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
          if keyword_set(opcl) then begin
             tomgrid_mk4_hdr_A = *trace_data. tomgrid_mk4_hdr
             tomgrid_mk4_A     = *trace_data. tomgrid_mk4
             fitgrid_mk4_hdr_A = *trace_data. fitgrid_mk4_hdr
             fitgrid_mk4_A     = *trace_data. fitgrid_mk4
          endif
       endif
  

  if keyword_set(kcor) then begin
                Ne_kcor_A = *trace_data.Ne_kcor 
             index_kcor_A = *trace_data.index_kcor
    index_sampling_kcor_A = *trace_data.index_sampling_kcor
           rad_fit_kcor_A = *trace_data. rad_fit_kcor             
            Ne_fit_kcor_A = *trace_data.  Ne_fit_kcor
           scN_fit_kcor_A = *trace_data. scN_fit_kcor
           fitflag_kcor_A = *trace_data. fitflag_kcor
          fit_F_Ne_kcor   = *trace_data.fit_F_Ne_kcor
            lN_fit_kcor_A = *trace_data.  lN_fit_kcor 
          if fit_F_Ne_kcor eq 'DPL' then begin
             N1_fit_kcor_A = *trace_data. N1_fit_kcor
             N2_fit_kcor_A = *trace_data. N2_fit_kcor
             p1_fit_kcor_A = *trace_data. p1_fit_kcor
             p2_fit_kcor_A = *trace_data. p2_fit_kcor
          endif
          if keyword_set(opcl) then begin
             tomgrid_kcor_hdr_A = *trace_data. tomgrid_kcor_hdr
             tomgrid_kcor_A     = *trace_data. tomgrid_kcor
             fitgrid_kcor_hdr_A = *trace_data. fitgrid_kcor_hdr
             fitgrid_kcor_A     = *trace_data. fitgrid_kcor
          endif
       endif

  if keyword_set(ucomp) then begin
                Ne_ucomp_A = *trace_data.Ne_ucomp 
             index_ucomp_A = *trace_data.index_ucomp
    index_sampling_ucomp_A = *trace_data.index_sampling_ucomp
           rad_fit_ucomp_A = *trace_data. rad_fit_ucomp             
            Ne_fit_ucomp_A = *trace_data.  Ne_fit_ucomp
           scN_fit_ucomp_A = *trace_data. scN_fit_ucomp
           fitflag_ucomp_A = *trace_data. fitflag_ucomp
          fit_F_Ne_ucomp   = *trace_data.fit_F_Ne_ucomp
            lN_fit_ucomp_A = *trace_data.  lN_fit_ucomp 
          if fit_F_Ne_ucomp eq 'DPL' then begin
             N1_fit_ucomp_A = *trace_data. N1_fit_ucomp
             N2_fit_ucomp_A = *trace_data. N2_fit_ucomp
             p1_fit_ucomp_A = *trace_data. p1_fit_ucomp
             p2_fit_ucomp_A = *trace_data. p2_fit_ucomp
          endif
          if keyword_set(opcl) then begin
             tomgrid_ucomp_hdr_A = *trace_data. tomgrid_ucomp_hdr
             tomgrid_ucomp_A     = *trace_data. tomgrid_ucomp
             fitgrid_ucomp_hdr_A = *trace_data. fitgrid_ucomp_hdr
             fitgrid_ucomp_A     = *trace_data. fitgrid_ucomp
          endif
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
          if keyword_set(opcl) then begin
             tomgrid_c2_hdr_A = *trace_data. tomgrid_c2_hdr
             tomgrid_c2_A     = *trace_data. tomgrid_c2
             fitgrid_c2_hdr_A = *trace_data. fitgrid_c2_hdr
             fitgrid_c2_A     = *trace_data. fitgrid_c2
          endif
  endif

    return    
 end

  

