;-------------------------------------------------------------------------
;
; PURPOSE: This routine tests if the radsamp values cover in a
; reasonably uniform fashion the full range [Irmin,Irmax] or [Irmin,rad_fl_max].
;
; HISTORY:  v1.0 AMV & FAN, CLaSP, May 2024.
;           v1.1 AMV, CLaSP, May 2025.
;                A much more flexible and simpler routine.
;
;-------------------------------------------------------------------------

pro test_coverage, radsamp=radsamp, radfit_min=radfit_min, radfit_max=radfit_max, rad_fl_max=rad_fl_max, covgflag=covgflag
  common tomgrid,nr,nt,np,rmin,rmax,Irmin,Irmax 

; Set rad_1:
  rad_1 = Irmin

; Set rad_2 upon type of fl (open, large-closed, small-closed)
; Open OR large-closed:
  if rad_fl_max gt Irmax then rad_2 = Irmax
; Small-closed
  if rad_fl_max le Irmax then rad_2 = radfit_max ; apex!

  Nseg = 3 ; Let us test this value.
  DR   = (rad_2-rad_1)/float(Nseg)
 ;Set of Nseg+1 Ri values that bound Nseg equally long segments.
  Ri   = rad_1 + DR * findgen(Nseg+1)
 ;Test coverage.
  covgflag = 'no'
  if (where(radsamp ge Ri(0) and radsamp le Ri(1)))(0) ne -1 AND  $
     (where(radsamp gt Ri(1) and radsamp le Ri(2)))(0) ne -1 AND  $
     (where(radsamp gt Ri(2) and radsamp le Ri(3)))(0) ne -1 AND  $
     n_elements(radsamp) ge 5 $
  THEN covgflag = 'yes'
  
  return
end
