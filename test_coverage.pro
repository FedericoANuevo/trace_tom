;-------------------------------------------------------------------------
pro test_coverage, radsamp=radsamp, covgflag=covgflag, $
                   aia = aia, euvia = euvia, euvib = euvib, eit = eit,$
                   mk4 = mk4, kcor = kcor, lascoc2 = lascoc2

  covgflag = 'no'

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

  return
end
