;;
; Calling sequence examples:
;
; expand_structure, cr_number = 2082, map_number = 1, /euvia, /lascoc2, /ldpower
; expand_structure, cr_number = 2082, map_number = 7, /euvia, /lascoc2
; expand_structure, cr_number = 2099, map_number = 1, /aia  , /lascoc2, /mk4, /ldpower
; expand_structure, cr_number = 2099, map_number = 7, /aia  , /lascoc2, /mk4
;;

pro expand_structure, map_number=map_number, cr_number=cr_number, $
                      aia = aia, euvia = euvia, euvib = euvib, eit = eit,$
                      mk4 = mk4, kcor = kcor, lascoc2 = lascoc2,$
                      npower = npower, dpower = dpower, linfit = linfit, ldpower = ldpower

  
; Set filename of the structure corresponding to selected CR_NUM and MAP_NUM.
;
  dir = './'
  if cr_number eq 2082 then begin
     if map_number eq 1 then begin
        dir                = '/data1/DATA/trace_tom_files/CR2082/field_lines_geometry_map1/'
        structure_filename = 'list.map1.new.txt-tracing-structure-merge_euvia_lascoc2_sampled.sav'
     endif  
     if map_number eq 7 then begin
        dir                = '/data1/DATA/trace_tom_files/CR2082/field_lines_geometry_map7/'
        structure_filename = 'list.map7.new.txt-tracing-structure-merge_euvia_lascoc2_sampled.sav'
     endif
  endif
  if cr_number eq 2099 then begin
     if map_number eq 1 then begin
        dir                = '/data1/DATA/trace_tom_files/CR2099/field_lines_geometry_map1/'
        structure_filename = 'list.map1.new.txt-tracing-structure-merge_aia_mk4_lascoc2_sampled.sav'
     endif  
     if map_number eq 7 then begin
        dir                = '/data1/DATA/trace_tom_files/CR2099/field_lines_geometry_map7/'
        structure_filename = 'list.map7.new.txt-tracing-structure-merge_aia_mk4_lascoc2_sampled.sav'
     endif
  endif
  if NOT keyword_set(structure_filename) then begin
     print, 'Invalid CR and/or MAP number, try again.'
     return
  endif

; Restore the structure with the traced-data  
  restore,dir+structure_filename

; Number of field-lines  
  N_fl = *trace_data.N_FL
;;

  if keyword_set(dpower)  then suffix = '_dpow.sav'
  if keyword_set(ldpower) then suffix = '_ldpow.sav'
  if keyword_set(linfit)  then suffix = '_linf.sav'
  if keyword_set(npower)  then suffix = '_npow.sav'
  
  if keyword_set(aia)     then Inst_list   = ['aia']
  if keyword_set(euvia)   then Inst_list   = ['euvia']
  if keyword_set(mk4)     then Inst_list   = [Inst_list,'mk4']
  if keyword_set(lascoc2) then Inst_list   = [Inst_list,'c2']
  
  if keyword_set(dpower) then  npar = 4
  if keyword_set(ldpower)then  npar = 4
  if keyword_set(linfit) then  npar = 2
  if keyword_set(npower) then  npar = n_elements(Inst_list)
  fit_par_A = fltarr(N_fl,npar)

  for i_fl =0, N_fl-1 do begin
; This module concatenate the densities of all instrument and fit it.
     if keyword_set(aia) then begin
        Nsamp   = (*trace_data.Npt_aia)(i_fl)
        rad_aia = reform((*trace_data.rad_aia)(i_fl,0:Nsamp-1))
        Ne_aia  = reform((*trace_data.Ne_aia) (i_fl,0:Nsamp-1))
        rad_concat = rad_aia
        Ne_concat  = Ne_aia
        inst_concat = 'aia' + strarr(n_elements(rad_aia))
     endif
     if keyword_set(euvia) then begin
        Nsamp     = (*trace_data.Npt_euvia)(i_fl)
        rad_euvia = reform((*trace_data.rad_euvia)(i_fl,0:Nsamp-1))
        Ne_euvia  = reform((*trace_data.Ne_euvia) (i_fl,0:Nsamp-1))
        rad_concat = rad_euvia
        Ne_concat  = Ne_euvia
        inst_concat = 'euvia' + strarr(n_elements(rad_euvia))
     endif   
     if keyword_set(mk4) then begin
        Nsamp   = (*trace_data.Npt_mk4)(i_fl)
        rad_mk4 = reform((*trace_data.rad_mk4)(i_fl,0:Nsamp-1))
        Ne_mk4  = reform((*trace_data.Ne_mk4) (i_fl,0:Nsamp-1))
        rad_concat = [rad_concat,rad_mk4]
        Ne_concat  = [Ne_concat,Ne_mk4]
        inst_concat = [inst_concat , 'mk4' + strarr(n_elements(rad_mk4))]
     endif
     if keyword_set(lascoc2) then begin
        Nsamp   = (*trace_data.Npt_c2)(i_fl)
        rad_c2  = reform((*trace_data.rad_c2)(i_fl,0:Nsamp-1))
        Ne_c2   = reform((*trace_data.Ne_c2) (i_fl,0:Nsamp-1))
        rad_concat = [rad_concat,rad_c2]
        Ne_concat  = [Ne_concat,Ne_c2]
        inst_concat = [inst_concat , 'c2' + strarr(n_elements(rad_c2))]
     endif
   
     if keyword_set(npower) then $
        power_concat_fit, Inst_list, rad_concat, Ne_concat, inst_concat, A, chisq, /weighted
     if keyword_set(dpower) then $
        double_power_concat_fit, rad_concat, Ne_concat, inst_concat, A, chisq, /weighted
     if keyword_set(ldpower) then $
        double_power_concat_fit, rad_concat, alog10(Ne_concat), inst_concat, A, chisq, /weighted
     if keyword_set(linfit) then begin
        index = where(rad_concat lt 5.5) 
        linear_fit, rad_concat(index)-1, alog10(Ne_concat(index)), A, r2, chisqr = chisq, /linfit_idl
     endif

     fit_par_A(i_fl,*) = A 
  endfor

  trace_data = create_struct( trace_data ,$
                              'fit_par', ptr_new(fit_par_A) )


  nname                     = strlen(structure_filename)
  structure_filename_expand = strmid(structure_filename,0,nname-4)+suffix

; Save expanded structure 
  save, trace_data, filename = dir + structure_filename_expand
  
  STOP
  return
end
 

