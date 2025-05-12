;-------------------------------------------------------------------------
;
; PURPOSE: This routine tests if the radsamp values cover in a
; reasonably uniform fashion the full range [radfit_min,radfit_max].
;
; HISTORY:  v1.0 AMV & FAN, CLaSP, May 2024.
;           v1.1 AMV, CLaSP, May 2025.
;                A much more flexible and simpler routine.
;
;-------------------------------------------------------------------------

pro test_coverage, radsamp=radsamp, radfit_min=radfit_min, radfit_max=radfit_max, covgflag=covgflag,$
                   aia=aia, euvia=euvia, euvib=euvib, eit=eit,$
                   mk4=mk4, kcor=kcor, lascoc2=lascoc2
  
  Nseg = 4 ; Too much? maybe Nseg = 3 is enough.
  DR   = (radfit_max-radfit_min)/float(Nseg)
 ;Set of Nseg+1 Ri values that bound Nseg equally long segments.
  Ri   = radfit_min + DR * findgen(Nseg+1)
 ;Test coverage.
  covgflag = 'no'
  if (where(radsamp ge Ri(0) and radsamp le Ri(1)))(0) ne -1 AND  $
     (where(radsamp gt Ri(1) and radsamp le Ri(2)))(0) ne -1 AND  $
     (where(radsamp gt Ri(2) and radsamp le Ri(3)))(0) ne -1 AND  $
     (where(radsamp gt Ri(3) and radsamp le Ri(4)))(0) ne -1 THEN covgflag = 'yes'
  
  return
end
