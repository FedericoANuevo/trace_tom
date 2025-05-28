
goto,skip_tag_fullrange
; Highlight in red Box field lines for which tag_fullrange_A = +1 and there is a fit with scN < scN_crit
ifl_A = where(tag_box_A eq +1 AND tag_fullrange_A eq +1 AND fitflag_AIA_A eq +1 AND scN_fit_aia_A le scN_crit)
oplot,Footpoint_Lon_A(ifl_A),Footpoint_Lat_A(ifl_A),psym=4,th=2,color=red
; Compute and plot the average Ne(r) of the RED field lines
Ne_fit_aia_Avg = fltarr(n_elements(rad_fit_aia_A)) - 678.
index_selected_A = intarr(N_fl)
for ir=0,n_elements(rad_fit_aia_A)-1 do begin
   Nptsavg = 0
   for ifl=0,N_fl-1 do begin
      if tag_box_A[ifl] eq +1 AND tag_fullrange_A[ifl] eq +1 AND fitflag_AIA_A[ifl] eq +1 AND scN_fit_aia_A[ifl] le scN_crit then begin
         index_selected_A(ifl) = +1
         if Ne_fit_aia_A(ifl,ir) gt 0. then begin
            Nptsavg = Nptsavg + 1
            Ne_fit_aia_Avg(ir) = Ne_fit_aia_Avg(ir) + Ne_fit_aia_A(ifl,ir)
         endif
      endif      
   endfor ; ifl
   Ne_fit_aia_Avg(ir) = Ne_fit_aia_Avg(ir) / float(Nptsavg)
endfor ; ir
ifl_A = where(index_selected_A eq +1)
oplot,Footpoint_Lon_A(ifl_A),Footpoint_Lat_A(ifl_A),psym=4,th=2,color=red
skip_tag_fullrange:
