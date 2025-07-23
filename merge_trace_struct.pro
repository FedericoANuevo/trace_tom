;
; PURPOSE: This code merges the tracing of tomographic products based
; on data from different instruments along given AWSoM fieldlines, and
; store all the information in a pointer structure.
;
; It also performs, if requested (/fits) analytical fits to the
; various tomographic products along individual field lines, including
; the fits as part of the structure.
;
; INPUTS:
; fir_fl and fl_list: STRINGS. directory where field lines are
; located, and filename of list of names of field lines filenames.
;
; FLAGS: one per allowed instrument, use those for which you want to
; merge results, provided in any order.
;
; HISTORY: V1.0 AMV & FAN, CLaSP, October 2023.
;          V1.1 AMV, IAFE, November 2023. Added Sampling, expanded to
;          all instruments and performed overall polishing.
;          V2.0 AMV, IAFE, January 2024. Added fit to tomographic resuls.
;          V2.1 AMV, CLASP, May 2024. pointers for each instrument are
;          defined only if its data is used.
;          V2.2 FAN, CLASP, May 2024. Added WT and LDEM_FLAG to the
;          traced result.
;          V2.4 AMV, CLaSP, May 2024, expanded to include closed FL,
;          added line-label to structure
;          V3.0 FAN, IAFE, JUL 2025, optimization of the used RAM memory.
;          Add the calculation of footpoints vectors (before in load_traced_data_structure)
pro merge_trace_struct, fl_dir=fl_dir, fl_list=fl_list, $
                        aia=aia, euvia=euvia, euvib=euvib, eit=eit, $
                        mk4=mk4, kcor=kcor, ucomp=ucomp, lascoc2=lascoc2, $
                        struture_filename=structure_filename,$
                        trace_Bs=trace_Bs

  if not keyword_set(fl_dir) or not keyword_set(fl_list) then STOP

; Set up filename for output structure:
  if not keyword_set(structure_filename) then structure_filename = fl_list
  structure_filename = structure_filename + '-tracing-structure-merge'

; Read the list with the field-lines  
  N_fl     = 0L
  filename = ''
  openr,1,fl_dir+fl_list
  readf,1,N_fl

; Maximum number of point along the fieldline
  Npt_max = 10100 ; Ask JUDIT for an optimal value.
; Default value in all arrays.
  default = -678.
  
; Define needed arrays:
; 1st Line Geometry and Bfield:
    x_A = dblarr(N_fl,Npt_max) + default
    y_A = dblarr(N_fl,Npt_max) + default
    z_A = dblarr(N_fl,Npt_max) + default
  rad_A = dblarr(N_fl,Npt_max) + default
  lat_A = dblarr(N_fl,Npt_max) + default
  lon_A = dblarr(N_fl,Npt_max) + default
 if keyword_set(trace_Bs) then begin
     s_A = dblarr(N_fl,Npt_max) + default
    Br_A = dblarr(N_fl,Npt_max) + default
   Bth_A = dblarr(N_fl,Npt_max) + default
   Bph_A = dblarr(N_fl,Npt_max) + default
     B_A = dblarr(N_fl,Npt_max) + default
 endif  
   Npt_v = intarr(N_fl)
; 2nd Tomography Products: 
  if keyword_set(aia) then begin
     Ne_aia_A             = fltarr(N_fl,Npt_max) + default
     Tm_aia_A             = fltarr(N_fl,Npt_max) + default
     WT_aia_A             = fltarr(N_fl,Npt_max) + default
     ldem_flag_aia_A      = fltarr(N_fl,Npt_max) + default
     index_aia_A          = intarr(N_fl,Npt_max) + default
     index_sampling_aia_A = intarr(N_fl,Npt_max) + default
  endif
  if keyword_set(euvia) then begin
     Ne_euvia_A             = fltarr(N_fl,Npt_max) + default
     Tm_euvia_A             = fltarr(N_fl,Npt_max) + default
     WT_euvia_A             = fltarr(N_fl,Npt_max) + default
     ldem_flag_euvia_A      = fltarr(N_fl,Npt_max) + default
     index_euvia_A          = intarr(N_fl,Npt_max) + default
     index_sampling_euvia_A = intarr(N_fl,Npt_max) + default
  endif
  if keyword_set(euvib) then begin
     Ne_euvib_A             = fltarr(N_fl,Npt_max) + default
     Tm_euvib_A             = fltarr(N_fl,Npt_max) + default
     WT_euvib_A             = fltarr(N_fl,Npt_max) + default
     ldem_flag_euvib_A      = fltarr(N_fl,Npt_max) + default
     index_euvib_A          = intarr(N_fl,Npt_max) + default
     index_sampling_euvib_A = intarr(N_fl,Npt_max) + default
  endif
  if keyword_set(eit) then begin
     Ne_eit_A             = fltarr(N_fl,Npt_max) + default
     Tm_eit_A             = fltarr(N_fl,Npt_max) + default
     WT_eit_A             = fltarr(N_fl,Npt_max) + default
     ldem_flag_eit_A      = fltarr(N_fl,Npt_max) + default
     index_eit_A          = intarr(N_fl,Npt_max) + default
     index_sampling_eit_A = intarr(N_fl,Npt_max) + default
  endif
  if keyword_set(mk4) then begin
     Ne_mk4_A    = fltarr(N_fl,Npt_max) + default
     index_mk4_A = intarr(N_fl,Npt_max) + default
     index_sampling_mk4_A = intarr(N_fl,Npt_max) + default
  endif
  if keyword_set(kcor) then begin
     Ne_kcor_A    = fltarr(N_fl,Npt_max) + default
     index_kcor_A = intarr(N_fl,Npt_max) + default
     index_sampling_kcor_A = intarr(N_fl,Npt_max) + default
  endif
  if keyword_set(ucomp) then begin
     Ne_ucomp_A    = fltarr(N_fl,Npt_max) + default
     index_ucomp_A = intarr(N_fl,Npt_max) + default
     index_sampling_ucomp_A = intarr(N_fl,Npt_max) + default
  endif  
  if keyword_set(lascoc2) then begin
     Ne_c2_A    = fltarr(N_fl,Npt_max) + default
     index_c2_A = intarr(N_fl,Npt_max) + default
     index_sampling_c2_A = intarr(N_fl,Npt_max) + default
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
           N_l         = n_elements(x_l)
           Npt_v(i_fl) = N_l
           initialized = 'yes'
           x_A   (i_fl,0:N_l-1) = x_l
           y_A   (i_fl,0:N_l-1) = y_l
           z_A   (i_fl,0:N_l-1) = z_l
           rad_A (i_fl,0:N_l-1) = rad_l
           lat_A (i_fl,0:N_l-1) = lat_l
           lon_A (i_fl,0:N_l-1) = lon_l
           if keyword_set(trace_Bs) then begin
             s_A (i_fl,0:N_l-1) =   s_l
            Br_A (i_fl,0:N_l-1) =  Br_l
           Bth_A (i_fl,0:N_l-1) = Bth_l
           Bph_A (i_fl,0:N_l-1) = Bph_l
             B_A (i_fl,0:N_l-1) =   B_l
           endif
        endif
        Ne_aia_A            (i_fl,0:N_l-1) = Ne_aia_l
        Tm_aia_A            (i_fl,0:N_l-1) = Tm_aia_l
        WT_aia_A            (i_fl,0:N_l-1) = WT_aia_l
        ldem_flag_aia_A     (i_fl,0:N_l-1) = ldem_flag_aia_l
        index_aia_A         (i_fl,0:N_l-1) = index_aia_l
        index_sampling_aia_A(i_fl,0:N_l-1) = index_sampling_l
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
           N_l         = n_elements(x_l)
           Npt_v(i_fl) = N_l
           initialized = 'yes'
           x_A   (i_fl,0:N_l-1) = x_l
           y_A   (i_fl,0:N_l-1) = y_l
           z_A   (i_fl,0:N_l-1) = z_l
           rad_A (i_fl,0:N_l-1) = rad_l
           lat_A (i_fl,0:N_l-1) = lat_l
           lon_A (i_fl,0:N_l-1) = lon_l
           if keyword_set(trace_Bs) then begin
             s_A (i_fl,0:N_l-1) =   s_l
            Br_A (i_fl,0:N_l-1) =  Br_l
           Bth_A (i_fl,0:N_l-1) = Bth_l
           Bph_A (i_fl,0:N_l-1) = Bph_l
             B_A (i_fl,0:N_l-1) =   B_l
           endif
        endif
        Ne_euvia_A            (i_fl,0:N_l-1) = Ne_euvia_l
        Tm_euvia_A            (i_fl,0:N_l-1) = Tm_euvia_l
        WT_euvia_A            (i_fl,0:N_l-1) = WT_euvia_l
        ldem_flag_euvia_A     (i_fl,0:N_l-1) = ldem_flag_euvia_l
        index_euvia_A         (i_fl,0:N_l-1) = index_euvia_l    
        index_sampling_euvia_A(i_fl,0:N_l-1) = index_sampling_l
     endif

   ; EUVI-B
     if keyword_set(euvib)  then begin
        file_euvib = filename+'_euvib.out'
        if i_fl eq 0 then structure_filename = structure_filename + '_euvib'
        if NOT keyword_set(trace_Bs) then $
        readcol,fl_dir+file_euvib,x_l,y_l,z_l,rad_l,lat_l,lon_l,Ne_euvib_l,Tm_euvib_l,$
                WT_euvib_l, ldem_flag_euvib_l, index_euvib_l  ,FORMAT='D,D,D,D,D,D'
        if     keyword_set(trace_Bs) then $
        readcol,fl_dir+file_euvib,x_l,y_l,z_l,s_l,Br_l,Bth_l,Bph_l,B_l,rad_l,lat_l,lon_l,Ne_euvib_l,Tm_euvib_l,$
                WT_euvib_l, ldem_flag_euvib_l, index_euvib_l  ,FORMAT='D,D,D,D,D,D,D,D,D,D,D'
        sample_fl, Ne_l=Ne_euvib_l, index_l=index_euvib_l, index_sampling_l=index_sampling_l
        if initialized eq 'no' then begin
           N_l         = n_elements(x_l)
           Npt_v(i_fl) = N_l
           initialized = 'yes'
           x_A   (i_fl,0:N_l-1) = x_l
           y_A   (i_fl,0:N_l-1) = y_l
           z_A   (i_fl,0:N_l-1) = z_l
           rad_A (i_fl,0:N_l-1) = rad_l
           lat_A (i_fl,0:N_l-1) = lat_l
           lon_A (i_fl,0:N_l-1) = lon_l
           if keyword_set(trace_Bs) then begin
             s_A (i_fl,0:N_l-1) =   s_l
            Br_A (i_fl,0:N_l-1) =  Br_l
           Bth_A (i_fl,0:N_l-1) = Bth_l
           Bph_A (i_fl,0:N_l-1) = Bph_l
             B_A (i_fl,0:N_l-1) =   B_l
           endif
        endif
        Ne_euvib_A            (i_fl,0:N_l-1) = Ne_euvib_l
        Tm_euvib_A            (i_fl,0:N_l-1) = Tm_euvib_l
        WT_euvib_A            (i_fl,0:N_l-1) = WT_euvib_l
        ldem_flag_euvib_A     (i_fl,0:N_l-1) = ldem_flag_euvib_l
        index_euvib_A         (i_fl,0:N_l-1) = index_euvib_l
        index_sampling_euvib_A(i_fl,0:N_l-1) = index_sampling_l
     endif
   ; EIT
     if keyword_set(eit)    then begin
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
           N_l         = n_elements(x_l)
           Npt_v(i_fl) = N_l  
           initialized = 'yes'
           x_A   (i_fl,0:N_l-1) = x_l
           y_A   (i_fl,0:N_l-1) = y_l
           z_A   (i_fl,0:N_l-1) = z_l
           rad_A (i_fl,0:N_l-1) = rad_l
           lat_A (i_fl,0:N_l-1) = lat_l
           lon_A (i_fl,0:N_l-1) = lon_l  
           if keyword_set(trace_Bs) then begin
             s_A (i_fl,0:N_l-1) =   s_l
            Br_A (i_fl,0:N_l-1) =  Br_l
           Bth_A (i_fl,0:N_l-1) = Bth_l
           Bph_A (i_fl,0:N_l-1) = Bph_l
             B_A (i_fl,0:N_l-1) =   B_l
           endif
        endif
        Ne_eit_A            (i_fl,0:N_l-1) = Ne_eit_l
        Tm_eit_A            (i_fl,0:N_l-1) = Tm_eit_l
        WT_eit_A            (i_fl,0:N_l-1) = WT_eit_l
        ldem_flag_eit_A     (i_fl,0:N_l-1) = ldem_flag_eit_l
        index_eit_A         (i_fl,0:N_l-1) = index_eit_l
        index_sampling_eit_A(i_fl,0:N_l-1) = index_sampling_l
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
           N_l         = n_elements(x_l)
           Npt_v(i_fl) = N_l  
           initialized = 'yes'
           x_A   (i_fl,0:N_l-1) = x_l
           y_A   (i_fl,0:N_l-1) = y_l
           z_A   (i_fl,0:N_l-1) = z_l
           rad_A (i_fl,0:N_l-1) = rad_l
           lat_A (i_fl,0:N_l-1) = lat_l
           lon_A (i_fl,0:N_l-1) = lon_l
           if keyword_set(trace_Bs) then begin
             s_A (i_fl,0:N_l-1) =   s_l
            Br_A (i_fl,0:N_l-1) =  Br_l
           Bth_A (i_fl,0:N_l-1) = Bth_l
           Bph_A (i_fl,0:N_l-1) = Bph_l
             B_A (i_fl,0:N_l-1) =   B_l
           endif
        endif
        Ne_mk4_A   (i_fl,0:N_l-1) = Ne_mk4_l
        index_mk4_A(i_fl,0:N_l-1) = index_mk4_l
        index_sampling_mk4_A(i_fl,0:N_l-1) = index_sampling_l
     endif

   ; KCOR
     if keyword_set(kcor)    then begin
        file_kcor   = filename+'_kcor.out'
        if i_fl eq 0 then structure_filename = structure_filename + '_kcor'
        if NOT keyword_set(trace_Bs) then $
           readcol,fl_dir+file_kcor,x_l,y_l,z_l,rad_l,lat_l,lon_l,Ne_kcor_l,$
                   index_kcor_l  ,FORMAT='D,D,D,D,D,D'
        if     keyword_set(trace_Bs) then $
           readcol,fl_dir+file_kcor,x_l,y_l,z_l,s_l,Br_l,Bth_l,Bph_l,B_l,rad_l,lat_l,lon_l,Ne_kcor_l,$
                   index_kcor_l  ,FORMAT='D,D,D,D,D,D,D,D,D,D,D'
        sample_fl, Ne_l=Ne_kcor_l, index_l=index_kcor_l, index_sampling_l=index_sampling_l
        if initialized eq 'no' then begin
           N_l         = n_elements(x_l)
           Npt_v(i_fl) = N_l  
           initialized = 'yes'
           x_A   (i_fl,0:N_l-1) = x_l
           y_A   (i_fl,0:N_l-1) = y_l
           z_A   (i_fl,0:N_l-1) = z_l
           rad_A (i_fl,0:N_l-1) = rad_l
           lat_A (i_fl,0:N_l-1) = lat_l
           lon_A (i_fl,0:N_l-1) = lon_l  
           if keyword_set(trace_Bs) then begin
             s_A (i_fl,0:N_l-1) =   s_l
            Br_A (i_fl,0:N_l-1) =  Br_l
           Bth_A (i_fl,0:N_l-1) = Bth_l
           Bph_A (i_fl,0:N_l-1) = Bph_l
             B_A (i_fl,0:N_l-1) =   B_l
           endif
        endif
        Ne_kcor_A   (i_fl,0:N_l-1) = Ne_kcor_l
        index_kcor_A(i_fl,0:N_l-1) = index_kcor_l        
        index_sampling_kcor_A(i_fl,0:N_l-1) = index_sampling_l
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
           N_l         = n_elements(x_l)
           Npt_v(i_fl) = N_l  
           initialized = 'yes'
           x_A   (i_fl,0:N_l-1) = x_l
           y_A   (i_fl,0:N_l-1) = y_l
           z_A   (i_fl,0:N_l-1) = z_l
           rad_A (i_fl,0:N_l-1) = rad_l
           lat_A (i_fl,0:N_l-1) = lat_l
           lon_A (i_fl,0:N_l-1) = lon_l  
           if keyword_set(trace_Bs) then begin
             s_A (i_fl,0:N_l-1) =   s_l
            Br_A (i_fl,0:N_l-1) =  Br_l
           Bth_A (i_fl,0:N_l-1) = Bth_l
           Bph_A (i_fl,0:N_l-1) = Bph_l
             B_A (i_fl,0:N_l-1) =   B_l
           endif
        endif
        Ne_ucomp_A   (i_fl,0:N_l-1) = Ne_ucomp_l
        index_ucomp_A(i_fl,0:N_l-1) = index_ucomp_l        
        index_sampling_ucomp_A(i_fl,0:N_l-1) = index_sampling_l
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
           N_l         = n_elements(x_l)
           Npt_v(i_fl) = N_l  
           initialized = 'yes'
           x_A   (i_fl,0:N_l-1) = x_l
           y_A   (i_fl,0:N_l-1) = y_l
           z_A   (i_fl,0:N_l-1) = z_l
           rad_A (i_fl,0:N_l-1) = rad_l
           lat_A (i_fl,0:N_l-1) = lat_l
           lon_A (i_fl,0:N_l-1) = lon_l  
           if keyword_set(trace_Bs) then begin
             s_A (i_fl,0:N_l-1) =   s_l
            Br_A (i_fl,0:N_l-1) =  Br_l
           Bth_A (i_fl,0:N_l-1) = Bth_l
           Bph_A (i_fl,0:N_l-1) = Bph_l
             B_A (i_fl,0:N_l-1) =   B_l
           endif
        endif
        Ne_c2_A   (i_fl,0:N_l-1) = Ne_c2_l
        index_c2_A(i_fl,0:N_l-1) = index_c2_l    
        index_sampling_c2_A(i_fl,0:N_l-1) = index_sampling_l
     endif
                              
  endfor                        ; End loop in fieldlines    
  close,1

; POINTER-STRUCTURE  
; Create a pointer structure to store field line extraction information  
  trace_data = { N_fl:    ptr_new(N_fl)    ,$
                 Npt_max: ptr_new(Npt_max) ,$
                 Npt_v:   ptr_new(Npt_v)   ,$
                 x:       ptr_new(x_A)     ,$
                 y:       ptr_new(y_A)     ,$
                 z:       ptr_new(z_A)     ,$
                 rad:     ptr_new(rad_A)   ,$
                 lat:     ptr_new(lat_A)   ,$
                 lon:     ptr_new(lon_A)    }
; undefine,N_fl
  undefine,Npt_max
  undefine,Npt_v
  undefine,x_A
  undefine,y_A
  undefine,z_A
  undefine,rad_A
  undefine,lat_A
  undefine,lon_A


; Determine the Rad, Lat and Lon of the footppoint of each field line.
; 1D Arrays: radial index corresponding to Rmin for each field line:
  irmin=intarr(N_fl)
  for i=0,N_fl-1 do irmin(i)=where(abs( (*trace_data.rad)(i,*) ) eq min(abs( (*trace_data.rad)(i,*) ) ) )
; 1D Arrays: Footpoint Lon and Lat for each field line:
  Footpoint_Rad_A = fltarr(N_fl)
  Footpoint_Lon_A = fltarr(N_fl)
  Footpoint_Lat_A = fltarr(N_fl)
  for ifl = 0,N_fl-1 do begin
     Footpoint_Rad_A(ifl) = (*trace_data.rad)(ifl,irmin[ifl])
     Footpoint_Lon_A(ifl) = (*trace_data.lon)(ifl,irmin[ifl])
     Footpoint_Lat_A(ifl) = (*trace_data.lat)(ifl,irmin[ifl])
  endfor

  trace_data = create_struct( trace_data       ,$
     'footpoint_rad', ptr_new(Footpoint_Rad_A),$
     'footpoint_lat', ptr_new(Footpoint_Lat_A),$
     'footpoint_lon', ptr_new(Footpoint_Lon_A) )

  undefine,Footpoint_Rad_A
  undefine,Footpoint_Lat_A
  undefine,Footpoint_Lon_A
  
  if keyword_set(trace_Bs) then begin
     trace_data = create_struct( trace_data,$
                 's'  ,   ptr_new(s_A)     ,$
                 'Br' ,   ptr_new(Br_A)    ,$
                 'Bth',   ptr_new(Bth_A)   ,$
                 'Bph',   ptr_new(Bph_A)   ,$
                 'B'  ,   ptr_new(B_A)     )
     undefine,s_A
     undefine,Br_A    
     undefine,Bth_A
     undefine,Bph_A
     undefine,B_A     
  endif

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
                   'Ne_aia' , ptr_new(            Ne_aia_A) ,$
                   'Tm_aia' , ptr_new(            Tm_aia_A) ,$
                   'WT_aia' , ptr_new(            WT_aia_A) ,$
            'ldem_flag_aia' , ptr_new(     ldem_flag_aia_A) ,$     
                'index_aia' , ptr_new(         index_aia_A) ,$
       'index_sampling_aia' , ptr_new(index_sampling_aia_A)  )
     undefine,Ne_aia_A
     undefine,Tm_aia_A
     undefine,WT_aia_A
     undefine,ldem_flag_aia_A
     undefine,index_aia_A
     undefine,index_sampling_aia_A     
  endif
  
  if keyword_set(euvia) then begin
     trace_data = create_struct( trace_data ,$
                 'Ne_euvia' , ptr_new(            Ne_euvia_A) ,$
                 'Tm_euvia' , ptr_new(            Tm_euvia_A) ,$
                 'WT_euvia' , ptr_new(            WT_euvia_A) ,$
          'ldem_flag_euvia' , ptr_new(     ldem_flag_euvia_A) ,$
              'index_euvia' , ptr_new(         index_euvia_A) ,$
     'index_sampling_euvia' , ptr_new(index_sampling_euvia_A)  )
     undefine,Ne_euvia_A
     undefine,Tm_euvia_A
     undefine,WT_euvia_A
     undefine,ldem_flag_euvia_A
     undefine,index_euvia_A
     undefine,index_sampling_euvia_A
  endif
  if keyword_set(euvib) then begin
     trace_data = create_struct( trace_data ,$
                 'Ne_euvib' , ptr_new(            Ne_euvib_A) ,$
                 'Tm_euvib' , ptr_new(            Tm_euvib_A) ,$
                 'WT_euvib' , ptr_new(            WT_euvib_A) ,$
          'ldem_flag_euvib' , ptr_new(     ldem_flag_euvib_A) ,$
              'index_euvib' , ptr_new(         index_euvib_A) ,$
     'index_sampling_euvib' , ptr_new(index_sampling_euvib_A)  )
     undefine,Ne_euvib_A
     undefine,Tm_euvib_A
     undefine,WT_euvib_A
     undefine,ldem_flag_euvib_A
     undefine,index_euvib_A
     undefine,index_sampling_euvib_A
  endif
  if keyword_set(eit) then begin
     trace_data = create_struct( trace_data ,$
                   'Ne_eit' , ptr_new(            Ne_eit_A) ,$
                   'Tm_eit' , ptr_new(            Tm_eit_A) ,$
                   'WT_eit' , ptr_new(            WT_eit_A) ,$
            'ldem_flag_eit' , ptr_new(     ldem_flag_eit_A) ,$
                'index_eit' , ptr_new(         index_eit_A) ,$
       'index_sampling_eit' , ptr_new(index_sampling_eit_A)  )
     undefine,Ne_eit_A
     undefine,Tm_eit_A
     undefine,WT_eit_A
     undefine,ldem_flag_eit_A
     undefine,index_eit_A
     undefine,index_sampling_eit_A
  endif
  if keyword_set(mk4) then begin
     trace_data = create_struct( trace_data ,$
                   'Ne_mk4' , ptr_new(            Ne_mk4_A) ,$
                'index_mk4' , ptr_new(         index_mk4_A) ,$
       'index_sampling_mk4' , ptr_new(index_sampling_mk4_A)  )
     undefine,Ne_mk4_A
     undefine,index_mk4_A
     undefine,index_sampling_mk4_A
  endif
  if keyword_set(kcor) then begin
     trace_data = create_struct( trace_data ,$
                 'Ne_kcor' , ptr_new(            Ne_kcor_A) ,$
              'index_kcor' , ptr_new(         index_kcor_A) ,$
     'index_sampling_kcor' , ptr_new(index_sampling_kcor_A)  )
     undefine,Ne_kcor_A
     undefine,index_kcor_A
     undefine,index_sampling_kcor_A
  endif
  if keyword_set(ucomp) then begin
     trace_data = create_struct( trace_data ,$
                 'Ne_ucomp' , ptr_new(            Ne_ucomp_A) ,$
              'index_ucomp' , ptr_new(         index_ucomp_A) ,$
     'index_sampling_ucomp' , ptr_new(index_sampling_ucomp_A)  )
     undefine,Ne_ucomp_A
     undefine,index_ucomp_A
     undefine,index_sampling_ucomp_A
  endif
  if keyword_set(lascoc2) then begin
     trace_data = create_struct( trace_data ,$
                   'Ne_c2' , ptr_new(            Ne_c2_A) ,$
                'index_c2' , ptr_new(         index_c2_A) ,$
       'index_sampling_c2' , ptr_new(index_sampling_c2_A)  )
     undefine,Ne_c2_A
     undefine,index_c2_A
     undefine,index_sampling_c2_A
  endif

; Perform fits:
  fit_trace_data, aia=aia, euvia=euvia, euvib=euvib, eit=eit,$
                  mk4=mk4, kcor=kcor, ucomp=ucomp, lascoc2=lascoc2,$
                  fl_dir=fl_dir, trace_data = trace_data 
  
 ; Save structure in fl_dir:
  save, trace_data, filename = fl_dir + structure_filename + '.sav'

  print,'Output in:'
  print,fl_dir + structure_filename + '.sav'
  
  return
end
