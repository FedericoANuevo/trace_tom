;
; PURPOSE: Read-in tomgrid and set up fit-grid
;          Called by fit_trace_data.pro
;
; HISTORY: v1.0, AMV, CLaSP, 10-May-2025.
;

pro read_tomgrid_and_define_fitgrid,fl_dir=fl_dir,instr_string=instr_string
  common tomgrid,nr,nt,np,rmin,rmax,Irmin,Irmax
  common fitgrid,radmin_fit,radmax_fit,drad_fit,Npt_fit
  
   ;Set all the variables needed to read in the tomographic grid parameters.
    empty_string = '' & nr=0 & nt=0 & np=0 & rmin=0. & rmax=0. & Irmin=0. & Irmax=0.
  
   ;Read-in the tomographic computational ball grid parameters  
    openr,1,fl_dir+'tom.grid.'+instr_string+'.dat'
    readf,1,empty_string
    readf,1,nr,nt,np,rmin,rmax,Irmin,Irmax
    close,1

   ;Set radmin_fit and radmax_fit equal to the Instrumental FoV limits
    radmin_fit = Irmin
    radmax_fit = Irmax
    
   ;Set an adequate radial resolution for the fit
    if instr_string eq 'aia' or instr_string eq 'euvia' or instr_string eq 'euvib' or instr_string eq 'eit' or $
       instr_string eq 'mk4' or instr_string eq 'kcor'  or instr_string eq 'ucomp' then $
       drad_fit = 0.01 ;Rs
    if instr_string eq 'lascoc2' then $
       drad_fit = 0.10 ;Rs
    
   ;Compute de number of points in the fit
    Npt_fit = round((radmax_fit-radmin_fit)/drad_fit)

    return
end
