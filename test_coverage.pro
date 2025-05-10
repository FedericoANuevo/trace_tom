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
  if (where(radsamp gt Ri(0) and radsamp le Ri(1)))(0) ne -1 AND  $
     (where(radsamp gt Ri(1) and radsamp le Ri(2)))(0) ne -1 AND  $
     (where(radsamp gt Ri(2) and radsamp le Ri(3)))(0) ne -1 AND  $
     (where(radsamp gt Ri(3) and radsamp le Ri(4)))(0) ne -1 THEN covgflag = 'yes'

goto,SKIP_OLDCODE
  if keyword_set(aia) or keyword_set(euvia) or keyword_set(euvib)  or keyword_set(eit) then begin
     R0=1.00
     R1=1.10
     R2=1.15
     R3=1.20
     R4=1.25
     if (where(radsamp gt R0 and radsamp le R1))(0) ne -1 AND  $
        (where(radsamp gt R1 and radsamp le R2))(0) ne -1 AND  $
        (where(radsamp gt R2 and radsamp le R3))(0) ne -1 AND  $
        (where(radsamp gt R3 and radsamp le R4))(0) ne -1 THEN covgflag = 'yes'
  endif

  if keyword_set(mk4) then begin
     R0=1.1
     R1=1.2
     R2=1.3
     R3=1.4
     R4=1.5
     if (where(radsamp gt R0 and radsamp le R1))(0) ne -1 AND  $
        (where(radsamp gt R1 and radsamp le R2))(0) ne -1 AND  $
        (where(radsamp gt R2 and radsamp le R3))(0) ne -1 AND  $
        (where(radsamp gt R3 and radsamp le R4))(0) ne -1 THEN covgflag = 'yes'
  endif

  if keyword_set(lascoc2) then begin
     R0=2.5
     R1=3.0
     R2=4.0
     R3=5.0
     R4=6.0
     if (where(radsamp gt R0 and radsamp le R1))(0) ne -1 AND  $
        (where(radsamp gt R1 and radsamp le R2))(0) ne -1 AND  $
        (where(radsamp gt R2 and radsamp le R3))(0) ne -1 AND  $
        (where(radsamp gt R3 and radsamp le R4))(0) ne -1 THEN covgflag = 'yes'
  endif
SKIP_OLDCODE:
  
  return
end
