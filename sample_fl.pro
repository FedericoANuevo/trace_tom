;
; Purpose: Provided 1D-arrays Ne_l and index_l, this routine returns
; imdex_sampling_l, a 1D-array bof the same size as index_l, containg -1
; or +1. Values +1 indicate sampling indexes.
;
; HISTORY: V1.0 AMV, IAFE, November-2023.
; HISTORY: V1.0.1 FAN, IAFE, April-2024, goto,skip_sample added
; HISTORY: V1.0.1 FAN, ClaSP, May-2024, change GOTO for an IF
; HISTORY: V1.0.2 FAN, IAFE, Sept-2024, add goto ...
;
pro sample_fl, Ne_l=Ne_l, index_l=index_l, index_sampling_l=index_sampling_l
  
  index_sampling_l = lonarr(n_elements(index_l))
  il_range   = where(index_l ne -1)
  if il_range[0] ne -1 then begin
     il_min     = min(il_range)
     il_max     = max(il_range)
     il         = il_min
     while il le il_max do begin
        il1 = il
        while (index_l(il+1) eq index_l(il1)) do begin
           il=il+1
          ;next line: stop if the next pint is the last one.
           if il eq n_elements(index_l)-1 then goto,skip_while
        endwhile
        il2 = il
        ind_sampling_cell = fix(mean([il1,il2]))
        if Ne_l(ind_sampling_cell) gt 0. then index_sampling_l(ind_sampling_cell) = +1
        il = il2+1
     endwhile
  endif
  return
end
