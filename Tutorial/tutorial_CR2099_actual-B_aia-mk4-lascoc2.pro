pro mini_tutorial
  common data, N_fl, Npt_max, Npt_v, x_A, y_A, z_A, rad_A, lat_A, lon_A,$
     Ne_aia_A, Tm_aia_A, WT_aia_A, ldem_flag_aia_A, index_aia_A, index_sampling_aia_A,$
     Ne_euvia_A, Tm_euvia_A,  WT_euvia_A, ldem_flag_euvia_A, index_euvia_A, index_sampling_euvia_A,$
     Ne_euvib_A, Tm_euvib_A, WT_euvib_A, ldem_flag_euvib_A, index_euvib_A, index_sampling_euvib_A,$
     Ne_eit_A, Tm_eit_A, WT_eit_A, ldem_flag_eit_A, index_eit_A, index_sampling_eit_A,$
     Ne_mk4_A, index_mk4_A, index_sampling_mk4_A,$
     Ne_kcor_A, index_kcor_A, index_sampling_kcor_A,$
     Ne_c2_A, index_c2_A, index_sampling_c2_A,$
     rad_fit_aia_A, Ne_fit_aia_A, Tm_fit_aia_A, fitflag_aia_A,scN_fit_aia_A,scT_fit_aia_A,$
     rad_fit_euvia_A, Ne_fit_euvia_A, Tm_fit_euvia_A, fitflag_euvia_A,scN_fit_euvia_A,scT_fit_euvia_A,$
     rad_fit_euvib_A, Ne_fit_euvib_A, Tm_fit_euvib_A, fitflag_euvib_A,scN_fit_euvib_A,scT_fit_euvib_A,$
     rad_fit_eit_A, Ne_fit_eit_A, Tm_fit_eit_A, fitflag_eit_A,scN_fit_eit_A,scT_fit_eit_A,$
     rad_fit_c2_A, Ne_fit_c2_A, fitflag_c2_A,scN_fit_c2_A,$
     rad_fit_mk4_A, Ne_fit_mk4_A, fitflag_mk4_A,scN_fit_mk4_A,$
     N0_fit_aia_A,lN_fit_aia_A,T0_fit_aia_A,dTdr_fit_aia_A,$
     N0_fit_euvia_A,lN_fit_euvia_A,T0_fit_euvia_A,dTdr_fit_euvia_A,$
     N0_fit_euvib_A,lN_fit_euvib_A,T0_fit_euvib_A,dTdr_fit_euvib_A,$
     N0_fit_eit_A,lN_fit_eit_A,T0_fit_eit_A,dTdr_fit_eit_A,$
     N1_fit_aia_A,N2_fit_aia_A,p1_fit_aia_A,p2_fit_aia_A,$
     N1_fit_euvia_A,N2_fit_euvia_A,p1_fit_euvia_A,p2_fit_euvia_A,$
     N1_fit_euvib_A,N2_fit_euvib_A,p1_fit_euvib_A,p2_fit_euvib_A,$
     N1_fit_eit_A,N2_fit_eit_A,p1_fit_eit_A,p2_fit_eit_A,$
     N0_fit_mk4_A,lN_fit_mk4_A,$
     N1_fit_mk4_A,N2_fit_mk4_A,p1_fit_mk4_A,p2_fit_mk4_A,$
     N1_fit_c2_A,N2_fit_c2_A,p1_fit_c2_A,p2_fit_c2_A,$
     lN_fit_c2_A,$
     fit_F_Ne_aia,fit_F_Ne_mk4,fit_F_Ne_c2,$
     fit_F_Ne_euvia,fit_F_Ne_euvib,fit_F_eit_c2
  
; 1) Declare the DIR where the structure is located, and the filename.
  dir = './'
  structure_filename = 'CR2099_AWSoM-map1_tracing-structure-merge_aia_mk4_lascoc2.sav'
  structure_filename = 'CR2099_AWSoM-map7_tracing-structure-merge_aia_mk4_lascoc2.sav'

; 2) Load structure into memory and extract all available arrays from it.
  load_traced_data_structure, dir=dir, structure_filename=structure_filename, trace_data=trace_data, /aia, /mk4, /lascoc2

  goto,plots
 
print, 'Press SPACE BAR to continue.'
pause

; 3) See the full contents of the structure.
print
print, 'Bellow is the list of the full contents of the structure, named "trace_data", loaded into memory by the command:'
print
print, 'load_traced_data_structure, dir=dir, structure_filename=structure_filename, trace_data=trace_data, /aia, /mk4, /lascoc2'
print
print, 'For this particular structure only AIA, MK4 and C2 results exist.'
print
print,'-------------------------------'
print,'help, trace_data, /str'
print
print, 'Press SPACE BAR to continue.'
pause
print
help, trace_data, /str
print,'-------------------------------'
print

print, 'Press SPACE BAR to continue.'
pause

; 4) See the arrays extracted from the structure.
print
print, 'Besides loading into memory the structure above, the command "load_traced_data_structure...'
print, 'also creates in memory several variables and arrays with the results stored in the structure,'
print, 'so they are convenientely ready to be used.'
print
print, 'Let us see the output of the command "help" for all variables and arrays that MAY BE created from the structure:'
print
print, 'Press SPACE BAR to continue.'
pause
print
print,'-------------------------------'
help,N_fl, Npt_max, Npt_v, x_A, y_A, z_A, rad_A, lat_A, lon_A,$
     Ne_aia_A, Tm_aia_A, WT_aia_A, ldem_flag_aia_A, index_aia_A, index_sampling_aia_A,$
     Ne_euvia_A, Tm_euvia_A,  WT_euvia_A, ldem_flag_euvia_A, index_euvia_A, index_sampling_euvia_A,$
     Ne_euvib_A, Tm_euvib_A, WT_euvib_A, ldem_flag_euvib_A, index_euvib_A, index_sampling_euvib_A,$
     Ne_eit_A, Tm_eit_A, WT_eit_A, ldem_flag_eit_A, index_eit_A, index_sampling_eit_A,$
     Ne_mk4_A, index_mk4_A, index_sampling_mk4_A,$
     Ne_kcor_A, index_kcor_A, index_sampling_kcor_A,$
     Ne_c2_A, index_c2_A, index_sampling_c2_A,$
     rad_fit_aia_A, Ne_fit_aia_A, Tm_fit_aia_A, fitflag_aia_A,scN_fit_aia_A,scT_fit_aia_A,$
     rad_fit_euvia_A, Ne_fit_euvia_A, Tm_fit_euvia_A, fitflag_euvia_A,scN_fit_euvia_A,scT_fit_euvia_A,$
     rad_fit_euvib_A, Ne_fit_euvib_A, Tm_fit_euvib_A, fitflag_euvib_A,scN_fit_euvib_A,scT_fit_euvib_A,$
     rad_fit_eit_A, Ne_fit_eit_A, Tm_fit_eit_A, fitflag_eit_A,scN_fit_eit_A,scT_fit_eit_A,$
     rad_fit_c2_A, Ne_fit_c2_A, fitflag_c2_A,scN_fit_c2_A,$
     rad_fit_mk4_A, Ne_fit_mk4_A, fitflag_mk4_A,scN_fit_mk4_A,$
     N0_fit_aia_A,lN_fit_aia_A,T0_fit_aia_A,dTdr_fit_aia_A,$
     N0_fit_euvia_A,lN_fit_euvia_A,T0_fit_euvia_A,dTdr_fit_euvia_A,$
     N0_fit_euvib_A,lN_fit_euvib_A,T0_fit_euvib_A,dTdr_fit_euvib_A,$
     N0_fit_eit_A,lN_fit_eit_A,T0_fit_eit_A,dTdr_fit_eit_A,$
     N1_fit_aia_A,N2_fit_aia_A,p1_fit_aia_A,p2_fit_aia_A,$
     N1_fit_euvia_A,N2_fit_euvia_A,p1_fit_euvia_A,p2_fit_euvia_A,$
     N1_fit_euvib_A,N2_fit_euvib_A,p1_fit_euvib_A,p2_fit_euvib_A,$
     N1_fit_eit_A,N2_fit_eit_A,p1_fit_eit_A,p2_fit_eit_A,$
     N0_fit_mk4_A,lN_fit_mk4_A,$
     N1_fit_mk4_A,N2_fit_mk4_A,p1_fit_mk4_A,p2_fit_mk4_A,$
     N1_fit_c2_A,N2_fit_c2_A,p1_fit_c2_A,p2_fit_c2_A,$
     lN_fit_c2_A,$
     fit_F_Ne_aia,fit_F_Ne_mk4,fit_F_Ne_c2,$
     fit_F_Ne_euvia,fit_F_Ne_euvib,fit_F_eit_c2
print,'-------------------------------'
print

; 5) Some explanations follow.
print,'In this case only arrays wih filed line geometry, or pertaining AIA, Mk4 or C2, are defined.'
print
print,'In this case there are "N_fl=4" field lines that were traced between rmin=1 Rs and rmax=10 Rs. For each field line a maximum of "Npt_max=10100" points is allowed. This number will be optimized later on, depending on what we find to be a typical number of points for Judit field lines. For now we just set it as very large number. These are synthetic radial field lines with a uniform radial step of 0.0015 Rs, so that these four lines contain "Npt_v = [6000, 6000, 6000, 6000]" points (because (10-1)/0.0015=6000.'
print
print,'All variables named "name_A" are 2D arrays of dimensions "(N_fl,Npt_max)". For example: "rad_A(ifl,*)" contains the radial points of field line "ifl (=0,..,3 in this example)". Similarly, "Ne_aia_A(ifl,*)" contains the AIA-based DEMT result for Ne at points "rad_A(ifl,*)". Etc.'
print
print, 'Press SPACE BAR to continue.'
pause
print
print,'There is a KEY array, provided for each instrument, that deserves a separate explanation: "index_sampling_INSTRUMENT_A", with INSTRUMENT being in this case: "aia, mk4 or c2". For field line ifl, "index_sampling_INSTRUMENT_A(ifl,*)" contains values 0 and 1: the 1 values indicate a single point within a given tomographic cell grid that samples that cell: it is the mean point within the cell, and its value is 1 ONLY if there is a result for that INSTRUMENT in that cell. If its value is 0, that cell is off the reconstructed volume for that instrument, or the cell is within the volume reconstructed for that instrument but that specific cell has not been reconstructed (being either a ZDA for both VL-tomography and DEMT, or a location where the LDEM does not predict accurately all FBEs, in the case of DEMT only).' 
print
print, 'Press SPACE BAR to continue.'
pause

plots:

; Set up graph options.
Device, retain = 2, true_color = 24, decomposed = 0

window,0,xs=1000,ys=1000
!p.multi=[0,1,2]
loadct,0
!p.color=0
!p.background=255


; color table for fieldlines
ctbl = 12 ; 16-LEVEL

; 1D Arrays of radial index for Rmin and Rmax 
irmin=intarr(N_FL)
irmax=intarr(N_FL)
for i=0,N_FL-1 do irmax(i)=where(    rad_A(i,*)  eq max(    rad_A(i,*)) )
for i=0,N_FL-1 do irmin(i)=where(abs(rad_A(i,*)) eq min(abs(rad_A(i,*))))
; 1D Arrays of Footpoint and Terminal Lon and Lat
Footpoint_Lon = fltarr(N_FL)
Footpoint_Lat = fltarr(N_FL)
 Terminal_Lon = fltarr(N_FL)
 Terminal_Lat = fltarr(N_FL)
for ifl = 0,N_FL-1 do Footpoint_Lon(ifl) = lon_A(ifl,irmin[ifl])
for ifl = 0,N_FL-1 do Footpoint_Lat(ifl) = lat_A(ifl,irmin[ifl])
for ifl = 0,N_FL-1 do  Terminal_Lon(ifl) = lon_A(ifl,irmax[ifl])
for ifl = 0,N_FL-1 do  Terminal_Lat(ifl) = lat_A(ifl,irmax[ifl])

; ID groups of field lines by means of user-defined
; ranges of terminal Longitudes.
CritTermLon = [0.,100.,180.,270.,360.]
CritTermLon = [0.,100.,180.,270.,310.,360.]
Ngroups     = n_elements(CritTermLon)-1
; ID each field line.
line_groupID  = intarr(N_FL)
for ig=0,Ngroups-1 do begin
   ifl_A = where(Terminal_Lon gt CritTermLon(ig) AND Terminal_Lon lt CritTermLon(ig+1))
   line_groupID(ifl_A) = ig
endfor

; Use the 16-level color table to color each group:
colors = 16 + 190 * (indgen(Ngroups))/float(Ngroups-1)

; Lat/Lon plots of FootPoint and TerminalPoint
plot,lon_A,lat_A,xr=[0,360],yr=[-90,+90],xstyle=1,ystyle=1,/nodata,charsize=2,title=strmid(structure_filename,0,17)+'  r = 1.0 Rs',ytitle='Carrington Latitude [deg]'
loadct,ctbl
for ifl=0,N_FL-1 do oplot,[Footpoint_Lon(ifl)],[Footpoint_Lat(ifl)],psym=4,color=colors(line_groupID(ifl)),th=2
loadct,0
plot,lon_A,lat_A,xr=[0,360],yr=[-90,+90],xstyle=1,ystyle=1,/nodata,charsize=2,xtitle='Carrington Longitude [deg]',title='r = 23.68 Rs',ytitle='Carrington Latitude [deg]'
loadct,ctbl
for ifl=0,N_FL-1 do oplot,[Terminal_Lon(ifl)],[Terminal_Lat(ifl)],psym=4,color=colors(line_groupID(ifl)),th=2
loadct,0
!p.multi=0
print, 'Press SPACE BAR to see the plot.'
;pause
stop

irmin=intarr(N_FL)
irmax=intarr(N_FL)
for i=0,N_FL-1 do irmax(i)=where(    rad_A(i,*)  eq max(    rad_A(i,*)) )
for i=0,N_FL-1 do irmin(i)=where(abs(rad_A(i,*)) eq min(abs(rad_A(i,*))))
loadct,0
plot,lon_A,lat_A,xr=[0,360],yr=[-90,+90],xstyle=1,ystyle=1,/nodata,charsize=2,title='r = 1.0 Rs',ytitle='Carrington Latitude [deg]'
loadct,12
for ifl=0,N_FL-1 do oplot,[lon_A(ifl,irmin[ifl])],[lat_A(ifl,irmin[ifl])],psym=4,color=ifl*32
loadct,0
plot,lon_A,lat_A,xr=[0,360],yr=[-90,+90],xstyle=1,ystyle=1,/nodata,charsize=2,xtitle='Carrington Longitude [deg]',title='r = 23.68 Rs',ytitle='Carrington Latitude [deg]'
loadct,12
for ifl=0,N_FL-1 do oplot,[lon_A(ifl,irmax[ifl])],[lat_A(ifl,irmax[ifl])],psym=4,color=ifl*32
loadct,39
!p.multi=0
stop
print, 'Press SPACE BAR to see the plot.'
pause

print
print,'Let us see an example of result from the structure.'
print,'Say one wants to plot the AIA results along field line ifl=0, then one does this:'
print
print,'          ifl = 0'
print,'          tmp = reform(index_sampling_aia_A(ifl,*))'
print,' ind_samp_aia = where(tmp eq 1)'
print,' window, 0'
print,' plot,rad_A(ifl,ind_samp_aia),Ne_aia_A(ifl,ind_samp_aia)'
print, 'Press SPACE BAR to see the plot.'
pause
window,0
ifl=0
tmp = reform(index_sampling_aia_A(ifl,*))
ind_samp_aia = where(tmp eq 1)
plot,rad_A(ifl,ind_samp_aia),Ne_aia_A(ifl,ind_samp_aia),charsize=2,xtitle='r [Rsun]',title='AIA-DEMT Ne(r) [cm!U-3!N]',psym=4,th=4, /nodata, yr=[0,1.e8], ystyle=1, xr=[1,1.3], xstyle=1
loadct,12
Ne_fit_aia_avg = 0. * rad_fit_aia_A

for ifl=0,N_fl-1 do begin 
  print, 'Press SPACE BAR to plot next line.'
  pause
  tmp = reform(index_sampling_aia_A(ifl,*))
  ind_samp_aia = where(tmp eq 1)
  col = (ifl+1)*40
  oplot,rad_A(ifl,ind_samp_aia),Ne_aia_A(ifl,ind_samp_aia),psym=4,th=2,color=col
  if fitflag_AIA_A(ifl) eq +1. then begin
    oplot,rad_fit_aia_A,Ne_fit_aia_A(ifl,*),color=col
    Ne_fit_aia_avg = Ne_fit_aia_avg + reform(Ne_fit_aia_A(ifl,*))
    print,'Fit Score:',scN_fit_aia_A(ifl)
  endif
endfor
   print, 'Press SPACE BAR to plot average trend.'
   pause
  N_fits = n_elements( where(fitflag_AIA_A eq +1.) )
  Ne_fit_aia_avg = Ne_fit_aia_avg / float(N_fits)
  loadct,0
  oplot,rad_fit_aia_A,Ne_fit_aia_avg,th=4

   print, 'Press SPACE BAR to see the Te(r) plot for the same field lines.'
   pause

window,1
ifl=0
tmp = reform(index_sampling_aia_A(ifl,*))
ind_samp_aia = where(tmp eq 1)
MK = 1.e6 ; K
plot,rad_A(ifl,ind_samp_aia),Tm_aia_A(ifl,ind_samp_aia)/MK,charsize=2,xtitle='r [Rsun]',title='AIA-DEMT Te(r) [MK]',psym=4,th=4, /nodata, yr=[0.,2.], ystyle=1, xr=[1,1.3], xstyle=1
loadct,12
Tm_fit_aia_avg = 0. * rad_fit_aia_A
for ifl=0,N_fl-1 do begin
   print, 'Press SPACE BAR to plot next line.'
   pause
  tmp = reform(index_sampling_aia_A(ifl,*))
  ind_samp_aia = where(tmp eq 1)
  col = (ifl+1)*40
  oplot,rad_A(ifl,ind_samp_aia),Tm_aia_A(ifl,ind_samp_aia)/MK,psym=4,th=2,color=col
  if fitflag_AIA_A(ifl) eq +1. then begin
    oplot,rad_fit_aia_A,Tm_fit_aia_A(ifl,*)/MK,color=col
    Tm_fit_aia_avg = Tm_fit_aia_avg + reform(Tm_fit_aia_A(ifl,*))
    print,'Fit Score:',scT_fit_aia_A(ifl)
  endif
endfor
   print, 'Press SPACE BAR to plot average trend.'
   pause
  N_fits = n_elements( where(fitflag_AIA_A eq +1.) )
  Tm_fit_aia_avg = Tm_fit_aia_avg / float(N_fits)
  loadct,0
  oplot,rad_fit_aia_A,Tm_fit_aia_avg/MK,th=4

  
print
print, 'Press SPACE BAR to continue.'
pause
print
print,'Similarly, to plot the c2 results for the SAME field line, one does:'
print
print,'         ifl = 0'
print,'         tmp = reform(index_sampling_c2_A(ifl,*))'
print,' ind_samp_c2 = where(tmp eq 1)'
print,' window, 2'
print,' plot,rad_A(ifl,ind_samp_c2),Ne_c2_A(ifl,ind_samp_c2)'
print, 'Press SPACE BAR to see the plot.'
pause
c2:
ifl=0
tmp = reform(index_sampling_c2_A(ifl,*))
ind_samp_c2 = where(tmp eq 1)

window,2
 plot,rad_A(ifl,ind_samp_c2),Ne_c2_A(ifl,ind_samp_c2),charsize=2,xtitle='r [Rsun]',title='C2-SRT Ne(r) [cm!U-3!N]',psym=4,th=4,/nodata, xr=[2.5,6.0], xstyle=1, yr=[0,7.e4], ystyle=1
loadct,12
Ne_fit_c2_avg = 0. * rad_fit_c2_A
for ifl=0,N_fl-1 do begin
   print, 'Press SPACE BAR to plot next line.'
   pause
  tmp = reform(index_sampling_c2_A(ifl,*))
  ind_samp_c2 = where(tmp eq 1)
  col = (ifl+1)*40
  oplot,rad_A(ifl,ind_samp_c2),Ne_c2_A(ifl,ind_samp_c2),psym=4,th=2,color=col
  if fitflag_c2_A(ifl) eq +1. then begin
    oplot,rad_fit_c2_A,Ne_fit_c2_A(ifl,*),color=col
    Ne_fit_c2_avg = Ne_fit_c2_avg + reform(Ne_fit_c2_A(ifl,*))
    print,'Fit Score:',scN_fit_c2_A(ifl)
  endif
endfor

   print, 'Press SPACE BAR to plot average trend.'
   pause
  N_fits = n_elements( where(fitflag_c2_A eq +1.) )
  Ne_fit_c2_avg = Ne_fit_c2_avg / float(N_fits)
  loadct,0
  oplot,rad_fit_c2_A,Ne_fit_c2_avg,th=4

print
print, 'Press SPACE BAR to continue.'
pause
print
print,'Similarly, to plot the Mk4 results for the SAME field line, one does:'
print
print,'         ifl = 0'
print,'         tmp = reform(index_sampling_mk4_A(ifl,*))'
print,' ind_samp_mk4 = where(tmp eq 1)'
print,' window, 3'
print,' plot,rad_A(ifl,ind_samp_mk4),Ne_mk4_A(ifl,ind_samp_c2)'
print, 'Press SPACE BAR to see the plot.'
pause
lastgraph:
ifl=0
tmp = reform(index_sampling_mk4_A(ifl,*))
ind_samp_mk4 = where(tmp eq 1)
window,3
 plot,rad_A(ifl,ind_samp_mk4),Ne_mk4_A(ifl,ind_samp_mk4),charsize=2,xtitle='r [Rsun]',title='MK4-SRT Ne(r) [cm!U-3!N]',psym=4,th=4,/nodata, xr=[1.15,1.5], xstyle=1, yr=[0,8.e7], ystyle=1
loadct,12
Ne_fit_mk4_avg = 0. * rad_fit_mk4_A
for ifl=0,N_fl-1 do begin
   print, 'Press SPACE BAR to plot next line.'
   pause
  tmp = reform(index_sampling_mk4_A(ifl,*))
  ind_samp_mk4 = where(tmp eq 1)
  col = (ifl+1)*40
  oplot,rad_A(ifl,ind_samp_mk4),Ne_mk4_A(ifl,ind_samp_mk4),psym=4,th=2,color=col
  if fitflag_mk4_A(ifl) eq +1. then begin
    oplot,rad_fit_mk4_A,Ne_fit_mk4_A(ifl,*),color=col
    Ne_fit_mk4_avg = Ne_fit_mk4_avg + reform(Ne_fit_mk4_A(ifl,*))
    print,'Fit Chisq:',scN_fit_mk4_A(ifl)
  endif
endfor
   print, 'Press SPACE BAR to plot average trend.'
   pause
  N_fits = n_elements( where(fitflag_mk4_A eq +1.) )
  Ne_fit_mk4_avg = Ne_fit_mk4_avg / float(N_fits)
  loadct,0
  oplot,rad_fit_mk4_A,Ne_fit_mk4_avg,th=4

print
print,'Note that rad_A is NOT associated to an instrument.'
print
print,'Note that the aia, or lascoc2, sampled points are MUCH less than Npt_max=15000:'
print
help, ind_samp_aia, rad_A(ifl,ind_samp_aia),Ne_aia_A(ifl,ind_samp_aia)
help, ind_samp_c2, rad_A(ifl,ind_samp_c2),Ne_c2_A(ifl,ind_samp_c2)
print
print,'This concludes this mini-tutorial, to show the core part of our approach. We can pass you a plotting routine we are building that, based on this idea explained above, easily and flexibly combines any provided instruments on a single graph, computes average trends, etc. Let us first know what you think of the approach and data format.'
print

  stop
  return
end
