pro read_vlsrt,tom_dir,tom_file,nr,nt,np,rmin,rmax,Irmin,Irmax,rad,lat,lon,N_e  
; Path to use de x-tools
  !PATH = Expand_Path('+/data1/tomography/SolarTom_idl') + ':' + !PATH
; Create the Tom. Grid  
  grid_array1D,rmin,rmax,nr,nt,np,rad,lat,lon
; Read WL results
  xread,dir=tom_dir,file=tom_file,nr=nr,nt=nt,np=np,map=N_e
; Deprecate heights off the instrumental range
  irad = where( rad gt Irmin and rad lt Irmax)
  rad  = rad(irad)
  N_e  = N_e(irad,*,*)
  nr   = n_elements(rad)
  return
end
