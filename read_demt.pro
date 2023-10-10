pro read_demt,tom_dir,tom_file,nr,nt,np,rmin,rmax,Irmin,Irmax,rad,lat,lon,N_e,T_e  
; Path to use de x-tools
  !PATH = Expand_Path('+/data1/tomography/SolarTom_idl') + ':' + !PATH
; Create the Tom. Grid  
  grid_array1D,rmin,rmax,nr,nt,np,rad,lat,lon
; Define name of products to read   
  file_Ne     = 'Ne_'+tom_file
  file_Te     = 'Te_'+tom_file
; Read DEMT results
  xread,dir=tom_dir,file=file_Ne,nr=nr,nt=nt,np=np,map=N_e
  xread,dir=tom_dir,file=file_Te,nr=nr,nt=nt,np=np,map=T_e  
; Deprecate heights off the instrumental range
  irad = where( rad gt Irmin and rad lt Irmax)
  rad  = rad(irad)
  N_e  = N_e(irad,*,*)
  T_e  = T_e(irad,*,*)
  nr   = n_elements(rad)
  return
end
