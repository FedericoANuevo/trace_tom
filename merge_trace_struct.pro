;
; PURPOSE: This code merges the tracing of tomographic products based
; on data from different instruments along given AWSoM/FDIPS-PFSS fieldlines, and
; store all the information in a pointer structure. 
;
; It also performs, if requested (/fits) analytical fits to the
; various tomographic products along individual field lines, including
; the fits as part of the structure.
;
;
; INPUTS:
; fir_fl and fl_list: STRINGS. directory where field lines are
; located, and filename of list of names of field lines filenames.
;
; FLAGS: one per allowed instrument, use those for which you want to
; merge results, provided in any order.
;
; HISTORY: V1.0 FAN, IAFE, July 2025.

pro merge_trace_struct, fl_dir=fl_dir, fl_list=fl_list, $
                        aia=aia, euvia=euvia, euvib=euvib, eit=eit, $
                        mk4=mk4, kcor=kcor, ucomp=ucomp, lascoc2=lascoc2, $
                        structure_filename=structure_filename,$
                        trace_Bs=trace_Bs

  common datastructure, trace_data
  
  if not keyword_set(fl_dir) or not keyword_set(fl_list) then STOP

; Set up filename for output structure:
  if not keyword_set(structure_filename) then structure_filename = fl_list
  structure_filename = structure_filename + '-tracing-structure-merge'

; Read the list with the field-lines  
  N_fl     = 0L
  filename = ''
  openr,1,fl_dir+fl_list
  readf,1,N_fl

; Maximum number of points along each fieldline,
; set based on experience. These values may need
; to be changed. Also, expand if new instruments
; are added.
; NOTE: for an instrument that varies its FoV a lot, lik Metis, we
; will need a formula such as: Npt_max = 1.5 * Nrad.
  Npt_max_aia     =  40
  Npt_max_euvia   =  40
  Npt_max_euvib   =  40
  Npt_max_eit     =  40
  Npt_max_mk4     =  80
  Npt_max_kcor    =  80
  Npt_max_ucomp   =  40
  Npt_max_c2      = 150
  
; Default value in all arrays.
  default = -678.
  
; Define needed arrays:
; 1st Line Geometry and Bfield:
  Footpoint_Rad_A = fltarr(N_fl)
  Footpoint_Lon_A = fltarr(N_fl)
  Footpoint_Lat_A = fltarr(N_fl)
  Termpoint_Rad_A = fltarr(N_fl)
  Termpoint_Lon_A = fltarr(N_fl)
  Termpoint_Lat_A = fltarr(N_fl)

; 2nd Tomography Products: 
  if keyword_set(aia) then begin
     Npt_aia              = intarr(N_fl)             + default
     Ne_aia_A             = fltarr(N_fl,Npt_max_aia) + default
     Tm_aia_A             = fltarr(N_fl,Npt_max_aia) + default
     WT_aia_A             = fltarr(N_fl,Npt_max_aia) + default
     ldem_flag_aia_A      = fltarr(N_fl,Npt_max_aia) + default
     rad_aia_A            = fltarr(N_fl,Npt_max_aia) + default
     lat_aia_A            = fltarr(N_fl,Npt_max_aia) + default
     lon_aia_A            = fltarr(N_fl,Npt_max_aia) + default
     if keyword_set(trace_Bs) then begin
        s_aia_A            = fltarr(N_fl,Npt_max_aia) + default
        Br_aia_A           = fltarr(N_fl,Npt_max_aia) + default
        Bth_aia_A          = fltarr(N_fl,Npt_max_aia) + default
        Bph_aia_A          = fltarr(N_fl,Npt_max_aia) + default
        B_aia_A            = fltarr(N_fl,Npt_max_aia) + default
     endif
  endif
  
  if keyword_set(euvia) then begin
     Npt_euvia              = intarr(N_fl)               + default
     Ne_euvia_A             = fltarr(N_fl,Npt_max_euvia) + default
     Tm_euvia_A             = fltarr(N_fl,Npt_max_euvia) + default
     WT_euvia_A             = fltarr(N_fl,Npt_max_euvia) + default
     ldem_flag_euvia_A      = fltarr(N_fl,Npt_max_euvia) + default
      rad_euvia_A           = fltarr(N_fl,Npt_max_euvia) + default
      lat_euvia_A           = fltarr(N_fl,Npt_max_euvia) + default
      lon_euvia_A           = fltarr(N_fl,Npt_max_euvia) + default
      if keyword_set(trace_Bs) then begin
         s_euvia_A           = fltarr(N_fl,Npt_max_euvia) + default
         Br_euvia_A          = fltarr(N_fl,Npt_max_euvia) + default
         Bth_euvia_A         = fltarr(N_fl,Npt_max_euvia) + default
         Bph_euvia_A         = fltarr(N_fl,Npt_max_euvia) + default
         B_euvia_A           = fltarr(N_fl,Npt_max_euvia) + default
      endif
   endif
  
  if keyword_set(euvib) then begin
     Npt_euvib              = intarr(N_fl)               + default
     Ne_euvib_A             = fltarr(N_fl,Npt_max_euvib) + default
     Tm_euvib_A             = fltarr(N_fl,Npt_max_euvib) + default
     WT_euvib_A             = fltarr(N_fl,Npt_max_euvib) + default
     ldem_flag_euvib_A      = fltarr(N_fl,Npt_max_euvib) + default
      rad_euvib_A           = fltarr(N_fl,Npt_max_euvib) + default
      lat_euvib_A           = fltarr(N_fl,Npt_max_euvib) + default
      lon_euvib_A           = fltarr(N_fl,Npt_max_euvib) + default
      if keyword_set(trace_Bs) then begin
         s_euvib_A           = fltarr(N_fl,Npt_max_euvib) + default
         Br_euvib_A          = fltarr(N_fl,Npt_max_euvib) + default
         Bth_euvib_A         = fltarr(N_fl,Npt_max_euvib) + default
         Bph_euvib_A         = fltarr(N_fl,Npt_max_euvib) + default
         B_euvib_A           = fltarr(N_fl,Npt_max_euvib) + default
      endif
   endif

  if keyword_set(eit) then begin
     Npt_eit              = intarr(N_fl)             + default
     Ne_eit_A             = fltarr(N_fl,Npt_max_eit) + default
     Tm_eit_A             = fltarr(N_fl,Npt_max_eit) + default
     WT_eit_A             = fltarr(N_fl,Npt_max_eit) + default
     ldem_flag_eit_A      = fltarr(N_fl,Npt_max_eit) + default
      rad_eit_A           = fltarr(N_fl,Npt_max_eit) + default
      lat_eit_A           = fltarr(N_fl,Npt_max_eit) + default
      lon_eit_A           = fltarr(N_fl,Npt_max_eit) + default
      if keyword_set(trace_Bs) then begin
        s_eit_A           = fltarr(N_fl,Npt_max_eit) + default
        Br_eit_A          = fltarr(N_fl,Npt_max_eit) + default
        Bth_eit_A         = fltarr(N_fl,Npt_max_eit) + default
        Bph_eit_A         = fltarr(N_fl,Npt_max_eit) + default
        B_eit_A           = fltarr(N_fl,Npt_max_eit) + default
     endif
   endif
  
  if keyword_set(mk4) then begin
     Npt_mk4              = intarr(N_fl)             + default
     Ne_mk4_A             = fltarr(N_fl,Npt_max_mk4) + default
     rad_mk4_A            = fltarr(N_fl,Npt_max_mk4) + default
     lat_mk4_A            = fltarr(N_fl,Npt_max_mk4) + default
     lon_mk4_A            = fltarr(N_fl,Npt_max_mk4) + default
     if keyword_set(trace_Bs) then begin
        s_mk4_A           = fltarr(N_fl,Npt_max_mk4) + default
        Br_mk4_A          = fltarr(N_fl,Npt_max_mk4) + default
        Bth_mk4_A         = fltarr(N_fl,Npt_max_mk4) + default
        Bph_mk4_A         = fltarr(N_fl,Npt_max_mk4) + default
        B_mk4_A           = fltarr(N_fl,Npt_max_mk4) + default
     endif
  endif
  
  if keyword_set(kcor) then begin
     Npt_kcor              = intarr(N_fl)              + default
     Ne_kcor_A             = fltarr(N_fl,Npt_max_kcor) + default
     rad_kcor_A            = fltarr(N_fl,Npt_max_kcor) + default
     lat_kcor_A            = fltarr(N_fl,Npt_max_kcor) + default
     lon_kcor_A            = fltarr(N_fl,Npt_max_kcor) + default
     if keyword_set(trace_Bs) then begin
        s_kcor_A           = fltarr(N_fl,Npt_max_kcor) + default
        Br_kcor_A          = fltarr(N_fl,Npt_max_kcor) + default
        Bth_kcor_A         = fltarr(N_fl,Npt_max_kcor) + default
        Bph_kcor_A         = fltarr(N_fl,Npt_max_kcor) + default
        B_kcor_A           = fltarr(N_fl,Npt_max_kcor) + default
     endif
  endif
  
  if keyword_set(ucomp) then begin
     Npt_ucomp              = intarr(N_fl)               + default
     Ne_ucomp_A             = fltarr(N_fl,Npt_max_ucomp) + default
     rad_ucomp_A            = fltarr(N_fl,Npt_max_ucomp) + default
     lat_ucomp_A            = fltarr(N_fl,Npt_max_ucomp) + default
     lon_ucomp_A            = fltarr(N_fl,Npt_max_ucomp) + default
     if keyword_set(trace_Bs) then begin
        s_ucomp_A           = fltarr(N_fl,Npt_max_ucomp) + default
        Br_ucomp_A          = fltarr(N_fl,Npt_max_ucomp) + default
        Bth_ucomp_A         = fltarr(N_fl,Npt_max_ucomp) + default
        Bph_ucomp_A         = fltarr(N_fl,Npt_max_ucomp) + default
        B_ucomp_A           = fltarr(N_fl,Npt_max_ucomp) + default
     endif
  endif
  
  if keyword_set(lascoc2) then begin
     Npt_c2              = intarr(N_fl)            + default
     Ne_c2_A             = fltarr(N_fl,Npt_max_c2) + default
     rad_c2_A            = fltarr(N_fl,Npt_max_c2) + default
     lat_c2_A            = fltarr(N_fl,Npt_max_c2) + default
     lon_c2_A            = fltarr(N_fl,Npt_max_c2) + default
     if keyword_set(trace_Bs) then begin
        s_c2_A           = fltarr(N_fl,Npt_max_c2) + default
        Br_c2_A          = fltarr(N_fl,Npt_max_c2) + default
        Bth_c2_A         = fltarr(N_fl,Npt_max_c2) + default
        Bph_c2_A         = fltarr(N_fl,Npt_max_c2) + default
        B_c2_A           = fltarr(N_fl,Npt_max_c2) + default
     endif
  
  endif
           
  for i_fl = 0L,N_fl-1 do begin ; Start loop in fieldlines  
     initialized = 'no'
     readf,1,filename
     
   ; AIA
     if keyword_set(aia)   then begin
        file_aia   = filename+'_aia.out'
        if i_fl eq 0 then structure_filename = structure_filename + '_aia'
        if NOT keyword_set(trace_Bs) then $
           readcol,fl_dir+file_aia,x_l,y_l,z_l,rad_l,lat_l,lon_l,Ne_aia_l,Tm_aia_l,$
                   WT_aia_l, ldem_flag_aia_l, index_aia_l  ,FORMAT='D,D,D,D,D,D'
        if     keyword_set(trace_Bs) then $
           readcol,fl_dir+file_aia,x_l,y_l,z_l,s_l,Br_l,Bth_l,Bph_l,B_l,rad_l,lat_l,lon_l,Ne_aia_l,Tm_aia_l,$
                   WT_aia_l, ldem_flag_aia_l, index_aia_l  ,FORMAT='D,D,D,D,D,D,D,D,D,D,D'
        sample_fl, Ne_l=Ne_aia_l, index_l=index_aia_l, index_sampling_l=index_sampling_l
        if initialized eq 'no' then begin
           index_foot = where( rad_l eq min(rad_l))
           index_term = where( rad_l eq max(rad_l))
           footpoint_rad_A(i_fl) = rad_l(index_foot[0])
           footpoint_lat_A(i_fl) = lat_l(index_foot[0])
           footpoint_lon_A(i_fl) = lon_l(index_foot[0])
           termpoint_rad_A(i_fl) = rad_l(index_term[0])
           termpoint_lat_A(i_fl) = lat_l(index_term[0])
           termpoint_lon_A(i_fl) = lon_l(index_term[0])
           initialized = 'yes'
        endif
        indsamp = where(index_sampling_l eq 1)
        if indsamp[0] ne -1 then begin
           Nl_samp                           = n_elements(indsamp)
           Npt_aia        (i_fl)             = Nl_samp 
           Ne_aia_A       (i_fl,0:Nl_samp-1) = Ne_aia_l(indsamp)
           Tm_aia_A       (i_fl,0:Nl_samp-1) = Tm_aia_l(indsamp)
           WT_aia_A       (i_fl,0:Nl_samp-1) = WT_aia_l(indsamp)
           ldem_flag_aia_A(i_fl,0:Nl_samp-1) = ldem_flag_aia_l(indsamp)
           rad_aia_A      (i_fl,0:Nl_samp-1) = rad_l(indsamp)
           lat_aia_A      (i_fl,0:Nl_samp-1) = lat_l(indsamp)
           lon_aia_A      (i_fl,0:Nl_samp-1) = lon_l(indsamp)
           if keyword_set(trace_Bs) then begin
              s_aia_A     (i_fl,0:Nl_samp-1) =   S_l(indsamp)
              Br_aia_A    (i_fl,0:Nl_samp-1) =  Br_l(indsamp)
              Bth_aia_A   (i_fl,0:Nl_samp-1) = Bth_l(indsamp)
              Bph_aia_A   (i_fl,0:Nl_samp-1) = Bph_l(indsamp)
              B_aia_A     (i_fl,0:Nl_samp-1) =   B_l(indsamp)
           endif
        endif
     endif

; EUVI-A
     if keyword_set(euvia)  then begin
        file_euvia = filename+'_euvia.out'
        if i_fl eq 0 then structure_filename = structure_filename + '_euvia'
        if NOT keyword_set(trace_Bs) then $
           readcol,fl_dir+file_euvia,x_l,y_l,z_l,rad_l,lat_l,lon_l,Ne_euvia_l,Tm_euvia_l,$
                   WT_euvia_l, ldem_flag_euvia_l, index_euvia_l  ,FORMAT='D,D,D,D,D,D'
        if     keyword_set(trace_Bs) then $
           readcol,fl_dir+file_euvia,x_l,y_l,z_l,s_l,Br_l,Bth_l,Bph_l,B_l,rad_l,lat_l,lon_l,Ne_euvia_l,Tm_euvia_l,$
                   WT_euvia_l, ldem_flag_euvia_l, index_euvia_l  ,FORMAT='D,D,D,D,D,D,D,D,D,D,D'
        sample_fl, Ne_l=Ne_euvia_l, index_l=index_euvia_l, index_sampling_l=index_sampling_l
        if initialized eq 'no' then begin
           index_foot = where( rad_l eq min(rad_l))
           index_term = where( rad_l eq max(rad_l))
           footpoint_rad_A(i_fl) = rad_l(index_foot[0])
           footpoint_lat_A(i_fl) = lat_l(index_foot[0])
           footpoint_lon_A(i_fl) = lon_l(index_foot[0])
           termpoint_rad_A(i_fl) = rad_l(index_term[0])
           termpoint_lat_A(i_fl) = lat_l(index_term[0])
           termpoint_lon_A(i_fl) = lon_l(index_term[0])
           initialized = 'yes'
        endif
        indsamp = where(index_sampling_l eq 1)
        if indsamp[0] ne -1 then begin
           Nl_samp          = n_elements(indsamp)
           Npt_euvia         (i_fl)             = Nl_samp 
           Ne_euvia_A        (i_fl,0:Nl_samp-1) = Ne_euvia_l(indsamp)
           Tm_euvia_A        (i_fl,0:Nl_samp-1) = Tm_euvia_l(indsamp)
           WT_euvia_A        (i_fl,0:Nl_samp-1) = WT_euvia_l(indsamp)
           ldem_flag_euvia_A (i_fl,0:Nl_samp-1) = ldem_flag_euvia_l(indsamp)
           rad_euvia_A       (i_fl,0:Nl_samp-1) = rad_l(indsamp)
           lat_euvia_A       (i_fl,0:Nl_samp-1) = lat_l(indsamp)
           lon_euvia_A       (i_fl,0:Nl_samp-1) = lon_l(indsamp)
           if keyword_set(trace_Bs) then begin
              s_euvia_A      (i_fl,0:Nl_samp-1) =   S_l(indsamp)
              Br_euvia_A     (i_fl,0:Nl_samp-1) =  Br_l(indsamp)
              Bth_euvia_A    (i_fl,0:Nl_samp-1) = Bth_l(indsamp)
              Bph_euvia_A    (i_fl,0:Nl_samp-1) = Bph_l(indsamp)
              B_euvia_A      (i_fl,0:Nl_samp-1) =   B_l(indsamp)
           endif  
        endif
     endif
        
; EUVI-B
     if keyword_set(euvib)  then begin
        file_euvib = filename+'_euvib.out'
        if i_fl eq 0 then structure_filename = structure_filename + '_euvib'
        if NOT keyword_set(trace_Bs) then $
           readcol,fl_dir+file_euvib,x_l,y_l,z_l,rad_l,lat_l,lon_l,Ne_euvib_l,Tm_euvib_l,$
                   WT_euvib_l, ldem_flag_euvib_l, index_euvib_l  ,FORMAT='D,D,D,D,D,D'
        if     keyword_set(trace_Bs) then $
           readcol,fl_dir+file_euvia,x_l,y_l,z_l,s_l,Br_l,Bth_l,Bph_l,B_l,rad_l,lat_l,lon_l,Ne_euvib_l,Tm_euvib_l,$
                   WT_euvib_l, ldem_flag_euvib_l, index_euvib_l  ,FORMAT='D,D,D,D,D,D,D,D,D,D,D'
        sample_fl, Ne_l=Ne_euvib_l, index_l=index_euvib_l, index_sampling_l=index_sampling_l
        if initialized eq 'no' then begin
           index_foot = where( rad_l eq min(rad_l))
           index_term = where( rad_l eq max(rad_l))
           footpoint_rad_A(i_fl) = rad_l(index_foot[0])
           footpoint_lat_A(i_fl) = lat_l(index_foot[0])
           footpoint_lon_A(i_fl) = lon_l(index_foot[0])
           termpoint_rad_A(i_fl) = rad_l(index_term[0])
           termpoint_lat_A(i_fl) = lat_l(index_term[0])
           termpoint_lon_A(i_fl) = lon_l(index_term[0])
           initialized = 'yes'
        endif
        indsamp = where(index_sampling_l eq 1)
        if indsamp[0] ne -1 then begin
           Nl_samp          = n_elements(indsamp)
           Npt_euvib         (i_fl)             = Nl_samp 
           Ne_euvib_A        (i_fl,0:Nl_samp-1) = Ne_euvib_l(indsamp)
           Tm_euvib_A        (i_fl,0:Nl_samp-1) = Tm_euvib_l(indsamp)
           WT_euvib_A        (i_fl,0:Nl_samp-1) = WT_euvib_l(indsamp)
           ldem_flag_euvib_A (i_fl,0:Nl_samp-1) = ldem_flag_euvib_l(indsamp)
           rad_euvib_A       (i_fl,0:Nl_samp-1) = rad_l(indsamp)
           lat_euvib_A       (i_fl,0:Nl_samp-1) = lat_l(indsamp)
           lon_euvib_A       (i_fl,0:Nl_samp-1) = lon_l(indsamp)
           if keyword_set(trace_Bs) then begin
              s_euvib_A      (i_fl,0:Nl_samp-1) =   S_l(indsamp)
              Br_euvib_A     (i_fl,0:Nl_samp-1) =  Br_l(indsamp)
              Bth_euvib_A    (i_fl,0:Nl_samp-1) = Bth_l(indsamp)
              Bph_euvib_A    (i_fl,0:Nl_samp-1) = Bph_l(indsamp)
              B_euvib_A      (i_fl,0:Nl_samp-1) =   B_l(indsamp)
           endif  
        endif
     endif
     
; EIT
     if keyword_set(eit)  then begin
        file_eit = filename+'_eit.out'
        if i_fl eq 0 then structure_filename = structure_filename + '_eit'
        if NOT keyword_set(trace_Bs) then $
           readcol,fl_dir+file_eit,x_l,y_l,z_l,rad_l,lat_l,lon_l,Ne_eit_l,Tm_eit_l,$
                   WT_eit_l, ldem_flag_eit_l, index_eit_l  ,FORMAT='D,D,D,D,D,D'
        if     keyword_set(trace_Bs) then $
           readcol,fl_dir+file_eit,x_l,y_l,z_l,s_l,Br_l,Bth_l,Bph_l,B_l,rad_l,lat_l,lon_l,Ne_eit_l,Tm_eit_l,$
                   WT_eit_l, ldem_flag_eit_l, index_eit_l  ,FORMAT='D,D,D,D,D,D,D,D,D,D,D'
        sample_fl, Ne_l=Ne_eit_l, index_l=index_eit_l, index_sampling_l=index_sampling_l
        if initialized eq 'no' then begin
           index_foot = where( rad_l eq min(rad_l))
           index_term = where( rad_l eq max(rad_l))
           footpoint_rad_A(i_fl) = rad_l(index_foot[0])
           footpoint_lat_A(i_fl) = lat_l(index_foot[0])
           footpoint_lon_A(i_fl) = lon_l(index_foot[0])
           termpoint_rad_A(i_fl) = rad_l(index_term[0])
           termpoint_lat_A(i_fl) = lat_l(index_term[0])
           termpoint_lon_A(i_fl) = lon_l(index_term[0])
           initialized = 'yes'
        endif
        indsamp = where(index_sampling_l eq 1)
        if indsamp[0] ne -1 then begin
           Nl_samp          = n_elements(indsamp)
           Npt_eit         (i_fl)             = Nl_samp 
           Ne_eit_A        (i_fl,0:Nl_samp-1) = Ne_euvib_l(indsamp)
           Tm_eit_A        (i_fl,0:Nl_samp-1) = Tm_euvib_l(indsamp)
           WT_eit_A        (i_fl,0:Nl_samp-1) = WT_euvib_l(indsamp)
           ldem_flag_eit_A (i_fl,0:Nl_samp-1) = ldem_flag_euvib_l(indsamp)
           rad_eit_A       (i_fl,0:Nl_samp-1) = rad_l(indsamp)
           lat_eit_A       (i_fl,0:Nl_samp-1) = lat_l(indsamp)
           lon_eit_A       (i_fl,0:Nl_samp-1) = lon_l(indsamp)
           if keyword_set(trace_Bs) then begin
              s_eit_A      (i_fl,0:Nl_samp-1) =   S_l(indsamp)
              Br_eit_A     (i_fl,0:Nl_samp-1) =  Br_l(indsamp)
              Bth_eit_A    (i_fl,0:Nl_samp-1) = Bth_l(indsamp)
              Bph_eit_A    (i_fl,0:Nl_samp-1) = Bph_l(indsamp)
              B_eit_A      (i_fl,0:Nl_samp-1) =   B_l(indsamp)
           endif  
        endif
     endif
     
   ; Mk4
     if keyword_set(mk4)    then begin
        file_mk4   = filename+'_mk4.out'
        if i_fl eq 0 then structure_filename = structure_filename + '_mk4'
        if NOT keyword_set(trace_Bs) then $
           readcol,fl_dir+file_mk4,x_l,y_l,z_l,rad_l,lat_l,lon_l,Ne_mk4_l,$
                   index_mk4_l  ,FORMAT='D,D,D,D,D,D'
        if     keyword_set(trace_Bs) then $
           readcol,fl_dir+file_mk4,x_l,y_l,z_l,s_l,Br_l,Bth_l,Bph_l,B_l,rad_l,lat_l,lon_l,Ne_mk4_l,$
                   index_mk4_l  ,FORMAT='D,D,D,D,D,D,D,D,D,D,D'
        sample_fl, Ne_l=Ne_mk4_l, index_l=index_mk4_l, index_sampling_l=index_sampling_l
        if initialized eq 'no' then begin
           index_foot = where( rad_l eq min(rad_l))
           index_term = where( rad_l eq max(rad_l))
           footpoint_rad_A(i_fl) = rad_l(index_foot[0])
           footpoint_lat_A(i_fl) = lat_l(index_foot[0])
           footpoint_lon_A(i_fl) = lon_l(index_foot[0])
           termpoint_rad_A(i_fl) = rad_l(index_term[0])
           termpoint_lat_A(i_fl) = lat_l(index_term[0])
           termpoint_lon_A(i_fl) = lon_l(index_term[0])
           initialized = 'yes'
        endif
        indsamp = where(index_sampling_l eq 1)
        if indsamp[0] ne -1 then begin
           Nl_samp                    = n_elements(indsamp)
           Npt_mk4 (i_fl)             = Nl_samp 
           Ne_mk4_A(i_fl,0:Nl_samp-1) = Ne_mk4_l(indsamp)
           rad_mk4_A       (i_fl,0:Nl_samp-1) = rad_l(indsamp)
           lat_mk4_A       (i_fl,0:Nl_samp-1) = lat_l(indsamp)
           lon_mk4_A       (i_fl,0:Nl_samp-1) = lon_l(indsamp)
           if keyword_set(trace_Bs) then begin
              s_mk4_A      (i_fl,0:Nl_samp-1) =   S_l(indsamp)
              Br_mk4_A     (i_fl,0:Nl_samp-1) =  Br_l(indsamp)
              Bth_mk4_A    (i_fl,0:Nl_samp-1) = Bth_l(indsamp)
              Bph_mk4_A    (i_fl,0:Nl_samp-1) = Bph_l(indsamp)
              B_mk4_A      (i_fl,0:Nl_samp-1) =   B_l(indsamp)
           endif  
        endif
     endif

   ; KCOR
     if keyword_set(kcor)    then begin
        file_kcor   = filename+'_kcor.out'
        if i_fl eq 0 then structure_filename = structure_filename + '_kcor'
        if NOT keyword_set(trace_Bs) then $
           readcol,fl_dir+file_mk4,x_l,y_l,z_l,rad_l,lat_l,lon_l,Ne_kcor_l,$
                   index_kcor_l  ,FORMAT='D,D,D,D,D,D'
        if     keyword_set(trace_Bs) then $
           readcol,fl_dir+file_kcor,x_l,y_l,z_l,s_l,Br_l,Bth_l,Bph_l,B_l,rad_l,lat_l,lon_l,Ne_kcor_l,$
                   index_kcor_l  ,FORMAT='D,D,D,D,D,D,D,D,D,D,D'
        sample_fl, Ne_l=Ne_kcor_l, index_l=index_kcor_l, index_sampling_l=index_sampling_l
        if initialized eq 'no' then begin
           index_foot = where( rad_l eq min(rad_l))
           index_term = where( rad_l eq max(rad_l))
           footpoint_rad_A(i_fl) = rad_l(index_foot[0])
           footpoint_lat_A(i_fl) = lat_l(index_foot[0])
           footpoint_lon_A(i_fl) = lon_l(index_foot[0])
           termpoint_rad_A(i_fl) = rad_l(index_term[0])
           termpoint_lat_A(i_fl) = lat_l(index_term[0])
           termpoint_lon_A(i_fl) = lon_l(index_term[0])
           initialized = 'yes'
        endif
        indsamp = where(index_sampling_l eq 1)
        if indsamp[0] ne -1 then begin
           Nl_samp                    = n_elements(indsamp)
           Npt_kcor (i_fl)             = Nl_samp 
           Ne_kcor_A(i_fl,0:Nl_samp-1) = Ne_kcor_l(indsamp)
           rad_kcor_A       (i_fl,0:Nl_samp-1) = rad_l(indsamp)
           lat_kcor_A       (i_fl,0:Nl_samp-1) = lat_l(indsamp)
           lon_kcor_A       (i_fl,0:Nl_samp-1) = lon_l(indsamp)
           if keyword_set(trace_Bs) then begin
              s_kcor_A      (i_fl,0:Nl_samp-1) =   S_l(indsamp)
              Br_kcor_A     (i_fl,0:Nl_samp-1) =  Br_l(indsamp)
              Bth_kcor_A    (i_fl,0:Nl_samp-1) = Bth_l(indsamp)
              Bph_kcor_A    (i_fl,0:Nl_samp-1) = Bph_l(indsamp)
              B_kcor_A      (i_fl,0:Nl_samp-1) =   B_l(indsamp)
           endif  
        endif
     endif
     
   ; UCOMP
     if keyword_set(ucomp)    then begin
        file_ucomp   = filename+'_ucomp.out'
        if i_fl eq 0 then structure_filename = structure_filename + '_ucomp'
        if NOT keyword_set(trace_Bs) then $
           readcol,fl_dir+file_ucomp,x_l,y_l,z_l,rad_l,lat_l,lon_l,Ne_ucomp_l,$
                   index_ucomp_l  ,FORMAT='D,D,D,D,D,D'
        if     keyword_set(trace_Bs) then $
           readcol,fl_dir+file_ucomp,x_l,y_l,z_l,s_l,Br_l,Bth_l,Bph_l,B_l,rad_l,lat_l,lon_l,Ne_ucomp_l,$
                   index_ucomp_l  ,FORMAT='D,D,D,D,D,D,D,D,D,D,D'
        sample_fl, Ne_l=Ne_ucomp_l, index_l=index_ucomp_l, index_sampling_l=index_sampling_l
        if initialized eq 'no' then begin
           index_foot = where( rad_l eq min(rad_l))
           index_term = where( rad_l eq max(rad_l))
           footpoint_rad_A(i_fl) = rad_l(index_foot[0])
           footpoint_lat_A(i_fl) = lat_l(index_foot[0])
           footpoint_lon_A(i_fl) = lon_l(index_foot[0])
           termpoint_rad_A(i_fl) = rad_l(index_term[0])
           termpoint_lat_A(i_fl) = lat_l(index_term[0])
           termpoint_lon_A(i_fl) = lon_l(index_term[0])
           initialized = 'yes'
        endif
        indsamp = where(index_sampling_l eq 1)
        if indsamp[0] ne -1 then begin
           Nl_samp                    = n_elements(indsamp)
           Npt_ucomp (i_fl)             = Nl_samp 
           Ne_ucomp_A(i_fl,0:Nl_samp-1) = Ne_ucomp_l(indsamp)
           rad_ucomp_A       (i_fl,0:Nl_samp-1) = rad_l(indsamp)
           lat_ucomp_A       (i_fl,0:Nl_samp-1) = lat_l(indsamp)
           lon_ucomp_A       (i_fl,0:Nl_samp-1) = lon_l(indsamp)
           if keyword_set(trace_Bs) then begin
              s_ucomp_A      (i_fl,0:Nl_samp-1) =   S_l(indsamp)
              Br_ucomp_A     (i_fl,0:Nl_samp-1) =  Br_l(indsamp)
              Bth_ucomp_A    (i_fl,0:Nl_samp-1) = Bth_l(indsamp)
              Bph_ucomp_A    (i_fl,0:Nl_samp-1) = Bph_l(indsamp)
              B_ucomp_A      (i_fl,0:Nl_samp-1) =   B_l(indsamp)
           endif  
        endif
     endif
     
   ; LASCO-C2
     if keyword_set(lascoc2) then begin
        file_c2    = filename+'_lascoc2.out'
        if i_fl eq 0 then structure_filename = structure_filename + '_lascoc2'
        if NOT keyword_set(trace_Bs) then $
           readcol,fl_dir+file_c2,x_l,y_l,z_l,rad_l,lat_l,lon_l,Ne_c2_l,$
                   index_c2_l  ,FORMAT='D,D,D,D,D,D'
        if     keyword_set(trace_Bs) then $
           readcol,fl_dir+file_c2,x_l,y_l,z_l,s_l,Br_l,Bth_l,Bph_l,B_l,rad_l,lat_l,lon_l,Ne_c2_l,$
                   index_c2_l  ,FORMAT='D,D,D,D,D,D,D,D,D,D,D'
        sample_fl, Ne_l=Ne_c2_l, index_l=index_c2_l, index_sampling_l=index_sampling_l
        if initialized eq 'no' then begin
            index_foot = where( rad_l eq min(rad_l))
           index_term = where( rad_l eq max(rad_l))
           footpoint_rad_A(i_fl) = rad_l(index_foot[0])
           footpoint_lat_A(i_fl) = lat_l(index_foot[0])
           footpoint_lon_A(i_fl) = lon_l(index_foot[0])
           termpoint_rad_A(i_fl) = rad_l(index_term[0])
           termpoint_lat_A(i_fl) = lat_l(index_term[0])
           termpoint_lon_A(i_fl) = lon_l(index_term[0])
           initialized = 'yes'
        endif
        indsamp = where(index_sampling_l eq 1)
        if indsamp[0] ne -1 then begin
           Nl_samp                  = n_elements(indsamp)
           Npt_c2    (i_fl)         = Nl_samp        
           Ne_c2_A   (i_fl,0:Nl_samp-1) = Ne_c2_l
           rad_c2_A  (i_fl,0:Nl_samp-1) = rad_l(indsamp)
           lat_c2_A  (i_fl,0:Nl_samp-1) = lat_l(indsamp)
           lon_c2_A  (i_fl,0:Nl_samp-1) = lon_l(indsamp)
           if keyword_set(trace_Bs) then begin
              s_c2_A      (i_fl,0:Nl_samp-1) =   S_l(indsamp)
              Br_c2_A     (i_fl,0:Nl_samp-1) =  Br_l(indsamp)
              Bth_c2_A    (i_fl,0:Nl_samp-1) = Bth_l(indsamp)
              Bph_c2_A    (i_fl,0:Nl_samp-1) = Bph_l(indsamp)
              B_c2_A      (i_fl,0:Nl_samp-1) =   B_l(indsamp)
           endif
        endif
     endif                 
  endfor                        ; End loop in fieldlines    
  close,1

; POINTER-STRUCTURE  
; Create a pointer structure to store field line extraction information  
  trace_data = { N_fl:                ptr_new(N_fl)                ,$
                 footpoint_rad:       ptr_new(footpoint_rad_A)     ,$
                 footpoint_lat:       ptr_new(footpoint_lat_A)     ,$
                 footpoint_lon:       ptr_new(footpoint_lon_A)     ,$
                 termpoint_rad:       ptr_new(termpoint_rad_A)     ,$
                 termpoint_lat:       ptr_new(termpoint_lat_A)     ,$
                 termpoint_lon:       ptr_new(termpoint_lon_A)      }
  undefine,footpoint_rad_A
  undefine,footpoint_lat_A
  undefine,footpoint_lon_A
  undefine,termpoint_rad_A
  undefine,termpoint_lat_A
  undefine,termpoint_lon_A

if keyword_set(trace_Bs) then begin
; Read-in the leg-label of all fieldlines.
  leg_label_A = lonarr(N_fl)
  leglab=0L
  openr,1,fl_dir+'legs-label.dat'
  readf,1,N_fl
  for i_fl=0L,N_fl-1 do begin
     readf,1,leglab
     leg_label_A(i_fl)=leglab
  endfor
  close,1
; Read-in the leg-length of all fieldlines.
  leg_length_A = fltarr(N_fl)
  leglength=0.
  openr,1,fl_dir+'legs-length.dat'
  readf,1,N_fl
  for i_fl=0L,N_fl-1 do begin
     readf,1,leglength
     leg_length_A(i_fl)=leglength
  endfor
  close,1
; Read-in the leg footpoint Bfield of all fieldlines.
  leg_footbfield_A = fltarr(N_fl,3)
  xx='' & Br=0. & Bth=0. & Bph=0.
  openr,1,fl_dir+'legs-footpoint-Bfield.dat'
  readf,1,N_fl
  readf,1,xx
  for i_fl=0L,N_fl-1 do begin
     readf,1,Br,Bth,Bph
     leg_footbfield_A(i_fl,*)=[Br,Bth,Bph]
  endfor
  close,1
; Add arrays to structure
  trace_data = create_struct( trace_data ,$
               'leg_label'      , ptr_new(leg_label_A)     ,$
               'leg_length'     , ptr_new(leg_length_A)    ,$
               'leg_footbfield' , ptr_new(leg_footbfield_A) )
  undefine,leg_label_A
  undefine,leg_length_A
  undefine,leg_footbfield_A
endif

; Store into the structure traced tomographic resuls
  if keyword_set(aia) then begin
     trace_data = create_struct( trace_data ,$
              'Npt_max_aia'  ,ptr_new(       Npt_max_aia  ) ,$                        
                  'Npt_aia' , ptr_new(           Npt_aia  ) ,$            
                   'Ne_aia' , ptr_new(            Ne_aia_A) ,$
                   'Tm_aia' , ptr_new(            Tm_aia_A) ,$
                   'WT_aia' , ptr_new(            WT_aia_A) ,$
            'ldem_flag_aia' , ptr_new(     ldem_flag_aia_A) ,$
                  'rad_aia' , ptr_new(           rad_aia_A) ,$
                  'lat_aia' , ptr_new(           lat_aia_A) ,$
                  'lon_aia' , ptr_new(           lon_aia_A) )
     undefine,Npt_aia
     undefine,Ne_aia_A
     undefine,Tm_aia_A
     undefine,WT_aia_A
     undefine,ldem_flag_aia_A
     undefine,rad_aia_A
     undefine,lat_aia_A
     undefine,lon_aia_A
     if keyword_set(trace_Bs) then begin
        trace_data = create_struct( trace_data ,$
                                    's_aia', ptr_new(s_aia_A), $
                                    'B_aia', ptr_new(B_aia_A), $
                                    'Br_aia', ptr_new(Br_aia_A), $
                                    'Bth_aia', ptr_new(Bth_aia_A), $
                                    'Bph_aia', ptr_new(Bph_aia_A) )
        undefine,s_aia_A
        undefine,B_aia_A
        undefine,Br_aia_A
        undefine,Bth_aia_A
        undefine,Bph_aia_A
     endif
  endif

; EUVI-A 
  if keyword_set(euvia) then begin
     trace_data = create_struct( trace_data ,$
              'Npt_max_euvia'  ,ptr_new(       Npt_max_euvia  ) ,$                        
                  'Npt_euvia' , ptr_new(           Npt_euvia  ) ,$            
                   'Ne_euvia' , ptr_new(            Ne_euvia_A) ,$
                   'Tm_euvia' , ptr_new(            Tm_euvia_A) ,$
                   'WT_euvia' , ptr_new(            WT_euvia_A) ,$
            'ldem_flag_euvia' , ptr_new(     ldem_flag_euvia_A) ,$
                  'rad_euvia' , ptr_new(           rad_euvia_A) ,$
                  'lat_euvia' , ptr_new(           lat_euvia_A) ,$
                  'lon_euvia' , ptr_new(           lon_euvia_A) )
     undefine,Npt_euvia
     undefine,Ne_euvia_A
     undefine,Tm_euvia_A
     undefine,WT_euvia_A
     undefine,ldem_flag_euvia_A
     undefine,rad_euvia_A
     undefine,lat_euvia_A
     undefine,lon_euvia_A
     if keyword_set(trace_Bs) then begin
        trace_data = create_struct( trace_data ,$
                                    's_euvia', ptr_new(s_euvia_A), $
                                    'B_euvia', ptr_new(B_euvia_A), $
                                    'Br_euvia', ptr_new(Br_euvia_A), $
                                    'Bth_euvia', ptr_new(Bth_euvia_A), $
                                    'Bph_euvia', ptr_new(Bph_euvia_A) )
        undefine,s_euvia_A
        undefine,B_euvia_A
        undefine,Br_euvia_A
        undefine,Bth_euvia_A
        undefine,Bph_euvia_A
     endif
  endif

; EUVI-B
  if keyword_set(euvib) then begin
     trace_data = create_struct( trace_data ,$
              'Npt_max_euvib'  ,ptr_new(       Npt_max_euvib  ) ,$                        
                  'Npt_euvib' , ptr_new(           Npt_euvib  ) ,$            
                   'Ne_euvib' , ptr_new(            Ne_euvib_A) ,$
                   'Tm_euvib' , ptr_new(            Tm_euvib_A) ,$
                   'WT_euvib' , ptr_new(            WT_euvib_A) ,$
            'ldem_flag_euvib' , ptr_new(     ldem_flag_euvib_A) ,$
                  'rad_euvib' , ptr_new(           rad_euvib_A) ,$
                  'lat_euvib' , ptr_new(           lat_euvib_A) ,$
                  'lon_euvib' , ptr_new(           lon_euvib_A) )
     undefine,Npt_euvib
     undefine,Ne_euvib_A
     undefine,Tm_euvib_A
     undefine,WT_euvib_A
     undefine,ldem_flag_euvib_A
     undefine,rad_euvib_A
     undefine,lat_euvib_A
     undefine,lon_euvib_A
     if keyword_set(trace_Bs) then begin
        trace_data = create_struct( trace_data ,$
                                    's_euvib', ptr_new(s_euvib_A), $
                                    'B_euvib', ptr_new(B_euvib_A), $
                                    'Br_euvib', ptr_new(Br_euvib_A), $
                                    'Bth_euvib', ptr_new(Bth_euvib_A), $
                                    'Bph_euvib', ptr_new(Bph_euvib_A) )
        undefine,s_euvib_A
        undefine,B_euvib_A
        undefine,Br_euvib_A
        undefine,Bth_euvib_A
        undefine,Bph_euvib_A
     endif
  endif
 ; EIT 
 if keyword_set(eit) then begin
    trace_data = create_struct( trace_data ,$
              'Npt_max_eit'  ,ptr_new(       Npt_max_eit  ) ,$                        
                  'Npt_eit' , ptr_new(           Npt_eit  ) ,$            
                   'Ne_eit' , ptr_new(            Ne_eit_A) ,$
                   'Tm_eit' , ptr_new(            Tm_eit_A) ,$
                   'WT_eit' , ptr_new(            WT_eit_A) ,$
            'ldem_flag_eit' , ptr_new(     ldem_flag_eit_A) ,$
                  'rad_eit' , ptr_new(           rad_eit_A) ,$
                  'lat_eit' , ptr_new(           lat_eit_A) ,$
                  'lon_eit' , ptr_new(           lon_eit_A) )
     undefine,Npt_eit
     undefine,Ne_eit_A
     undefine,Tm_eit_A
     undefine,WT_eit_A
     undefine,ldem_flag_eit_A
     undefine,rad_eit_A
     undefine,lat_eit_A
     undefine,lon_eit_A
     if keyword_set(trace_Bs) then begin
        trace_data = create_struct( trace_data ,$
                                    's_eit', ptr_new(s_eit_A), $
                                    'B_eit', ptr_new(B_eit_A), $
                                    'Br_eit', ptr_new(Br_eit_A), $
                                    'Bth_eit', ptr_new(Bth_eit_A), $
                                    'Bph_eit', ptr_new(Bph_eit_A) )
        undefine,s_eit_A
        undefine,B_eit_A
        undefine,Br_eit_A
        undefine,Bth_eit_A
        undefine,Bph_eit_A
     endif
  endif

;MK4
 if keyword_set(mk4) then begin
    trace_data = create_struct( trace_data ,$
              'Npt_max_mk4'  ,ptr_new(       Npt_max_mk4  ) ,$                        
                  'Npt_mk4' , ptr_new(           Npt_mk4  ) ,$            
                   'Ne_mk4' , ptr_new(            Ne_mk4_A) ,$
                  'rad_mk4' , ptr_new(           rad_mk4_A) ,$
                  'lat_mk4' , ptr_new(           lat_mk4_A) ,$
                  'lon_mk4' , ptr_new(           lon_mk4_A) )
     undefine,Npt_mk4
     undefine,Ne_mk4_A
     undefine,rad_mk4_A
     undefine,lat_mk4_A
     undefine,lon_mk4_A
     if keyword_set(trace_Bs) then begin
        trace_data = create_struct( trace_data ,$
                                    's_mk4', ptr_new(s_mk4_A), $
                                    'B_mk4', ptr_new(B_mk4_A), $
                                    'Br_mk4', ptr_new(Br_mk4_A), $
                                    'Bth_mk4', ptr_new(Bth_mk4_A), $
                                    'Bph_mk4', ptr_new(Bph_mk4_A) )
        undefine,s_mk4_A
        undefine,B_mk4_A
        undefine,Br_mk4_A
        undefine,Bth_mk4_A
        undefine,Bph_mk4_A
     endif
  endif
 
;KCOR
  if keyword_set(kcor) then begin
     trace_data = create_struct( trace_data ,$
              'Npt_max_kcor'  ,ptr_new(       Npt_max_kcor  ) ,$                        
                  'Npt_kcor' , ptr_new(           Npt_kcor  ) ,$            
                   'Ne_kcor' , ptr_new(            Ne_kcor_A) ,$
                  'rad_kcor' , ptr_new(           rad_kcor_A) ,$
                  'lat_kcor' , ptr_new(           lat_kcor_A) ,$
                  'lon_kcor' , ptr_new(           lon_kcor_A) )
     undefine,Npt_kcor
     undefine,Ne_kcor_A
     undefine,rad_kcor_A
     undefine,lat_kcor_A
     undefine,lon_kcor_A
     if keyword_set(trace_Bs) then begin
        trace_data = create_struct( trace_data ,$
                                    's_kcor', ptr_new(s_kcor_A), $
                                    'B_kcor', ptr_new(B_kcor_A), $
                                    'Br_kcor', ptr_new(Br_kcor_A), $
                                    'Bth_kcor', ptr_new(Bth_kcor_A), $
                                    'Bph_kcor', ptr_new(Bph_kcor_A) )
        undefine,s_kcor_A
        undefine,B_kcor_A
        undefine,Br_kcor_A
        undefine,Bth_kcor_A
        undefine,Bph_kcor_A
     endif
  endif

; UCOMP
  if keyword_set(ucomp) then begin
     trace_data = create_struct( trace_data ,$
              'Npt_max_ucomp' , ptr_new(       Npt_max_ucomp  ) ,$                        
                  'Npt_ucomp' , ptr_new(           Npt_ucomp  ) ,$            
                   'Ne_ucomp' , ptr_new(            Ne_ucomp_A) ,$
                  'rad_ucomp' , ptr_new(           rad_ucomp_A) ,$
                  'lat_ucomp' , ptr_new(           lat_ucomp_A) ,$
                  'lon_ucomp' , ptr_new(           lon_ucomp_A) )
     undefine,Npt_ucomp
     undefine,Ne_ucomp_A
     undefine,rad_ucomp_A
     undefine,lat_ucomp_A
     undefine,lon_ucomp_A
     if keyword_set(trace_Bs) then begin
        trace_data = create_struct( trace_data ,$
                                    's_ucomp', ptr_new(s_ucomp_A), $
                                    'B_ucomp', ptr_new(B_ucomp_A), $
                                    'Br_ucomp', ptr_new(Br_ucomp_A), $
                                    'Bth_ucomp', ptr_new(Bth_ucomp_A), $
                                    'Bph_ucomp', ptr_new(Bph_ucomp_A) )
        undefine,s_ucomp_A
        undefine,B_ucomp_A
        undefine,Br_ucomp_A
        undefine,Bth_ucomp_A
        undefine,Bph_ucomp_A
     endif
  endif

; LASCO-C2
  if keyword_set(lascoc2) then begin
     trace_data = create_struct( trace_data ,$
                             'Npt_max_c2' , ptr_new(       Npt_max_c2  ) ,$                        
                                 'Npt_c2' , ptr_new(           Npt_c2  ) ,$            
                                  'Ne_c2' , ptr_new(            Ne_c2_A) ,$
                                  'rad_c2', ptr_new(           rad_c2_A) ,$
                                  'lat_c2', ptr_new(           lat_c2_A) ,$
                                  'lon_c2', ptr_new(           lon_c2_A) )
     undefine,Npt_c2
     undefine,Ne_c2_A
     undefine,rad_c2_A
     undefine,lat_c2_A
     undefine,lon_c2_A
     if keyword_set(trace_Bs) then begin
        trace_data = create_struct( trace_data ,$
                                    's_c2', ptr_new(s_c2_A), $
                                    'B_c2', ptr_new(B_c2_A), $
                                    'Br_c2', ptr_new(Br_c2_A), $
                                    'Bth_c2', ptr_new(Bth_c2_A), $
                                    'Bph_c2', ptr_new(Bph_c2_A) )
        undefine,s_c2_A
        undefine,B_c2_A
        undefine,Br_c2_A
        undefine,Bth_c2_A
        undefine,Bph_c2_A
     endif

  endif

; Perform fits:
  fit_trace_data, aia=aia, euvia=euvia, euvib=euvib, eit=eit,$
                  mk4=mk4, kcor=kcor, ucomp=ucomp, lascoc2=lascoc2,$
                  fl_dir=fl_dir
  
 ; Save structure in fl_dir:
  save, trace_data, filename = fl_dir + structure_filename + '_sampled.sav'

  print,'Output in:'
  print,fl_dir + structure_filename + '_sampled.sav'
  
  return
end
