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
  structure_filename = 'CR2082_AWSoM-map1_tracing-structure-merge_euvia.sav'
  structure_filename = 'CR2082_AWSoM-map7_tracing-structure-merge_euvia.sav'
  
; 2) Load structure into memory and extract all available arrays from it.
  load_traced_data_structure, dir=dir, structure_filename=structure_filename, trace_data=trace_data, /euvia
  goto,plots
 
print, 'Press SPACE BAR to continue.'
pause

; 3) See the full contents of the structure.
print
print, 'Bellow is the list of the full contents of the structure, named "trace_data", loaded into memory by the command:'
print
print, 'load_traced_data_structure, dir=dir, structure_filename=structure_filename, trace_data=trace_data, /instruments....'
print
print, 'For this particular structure only EUVI-A results exist.'
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
print,'In this case only arrays wih filed line geometry, or pertaining EUVI-A, are defined,'
print
print,'In this case there are "N_fl=124" field lines that were traced in the height range'
print,'print, min(rad_A(where(rad_A gt 0.))), max(rad_A(where(rad_A gt 0.)))'
print, min(rad_A(where(rad_A gt 0.))), max(rad_A(where(rad_A gt 0.)))
print
print,'For each field line a maximum of "Npt_max=10100" points is allowed. This number will be optimized later on, depending on what we find to be a typical number of points for Judit field lines. For now we just set it as very large number.'
print
print,'All variables named "name_A" are 2D arrays of dimensions "(N_fl,Npt_max)". For example: "rad_A(ifl,*)" contains the radial points of field line "ifl (=0,..,3 in this example)". Similarly, "Ne_aia_A(ifl,*)" contains the AIA-based DEMT result for Ne at points "rad_A(ifl,*)". Etc.'
print
print, 'Press SPACE BAR to continue.'
pause
print
print,'There is a KEY array, provided for each instrument, that deserves a separate explanation: "index_sampling_INSTRUMENT_A", with INSTRUMENT being in this case: "euvia". For field line ifl, "index_sampling_INSTRUMENT_A(ifl,*)" contains values 0 and 1: the 1 values indicate a single point within a given tomographic cell grid that samples that cell: it is the mean point within the cell, and its value is 1 ONLY if there is a result for that INSTRUMENT in that cell. If its value is 0, that cell is off the reconstructed volume for that instrument, or the cell is within the volume reconstructed for that instrument but that specific cell has not been reconstructed (being either a ZDA for both VL-tomography and DEMT, or a location where the LDEM does not predict accurately all FBEs, in the case of DEMT only).' 
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

; Set a chisq threshold to select field lines:
chsq_threshold = 0.25

; Store average trends of Ne(r) and Tm(r) of each group in two 2D arrays.
Nrad_fit = n_elements(rad_fit_euvia_A)
Ne_avg = fltarr(Ngroups,Nrad_fit)
Tm_avg = fltarr(Ngroups,Nrad_fit)
for ig=0,Ngroups-1 do begin
   ;Get average trends for group ig.
   Nlines = 0
   Ne_fit_euvia_avg = fltarr(Nrad_fit)
   Tm_fit_euvia_avg = fltarr(Nrad_fit)
   for ifl = 0,N_FL-1 do begin
      if line_groupID(ifl) eq ig AND fitflag_EUVIA_A(ifl) eq +1. AND scN_fit_EUVIA_A(ifl) lt chsq_threshold AND scT_fit_EUVIA_A(ifl) lt chsq_threshold then begin
         Nlines = Nlines+1
         Ne_fit_euvia_avg = Ne_fit_euvia_avg + reform(Ne_fit_euvia_A(ifl,*))
         Tm_fit_euvia_avg = Tm_fit_euvia_avg + reform(Tm_fit_euvia_A(ifl,*))
      endif
   endfor
   Ne_avg(ig,*) = Ne_fit_euvia_avg / float(Nlines)
   Tm_avg(ig,*) = Tm_fit_euvia_avg / float(Nlines)
endfor

; Plot median fieldline and average trend of each group
window,1,xs=1000,ys=1000
!p.multi=[0,1,2]
ifl=0
tmp = reform(index_sampling_euvia_A(ifl,*))
ind_samp_euvia = where(tmp eq 1)
loadct,0
plot,rad_fit_euvia_A,Ne_avg(0,*),charsize=2,xtitle='r [Rsun]',title='EUVIA-DEMT Ne(r) [cm!U-3!N]',$
     psym=4,th=4, /nodata, yr=[0,3.e8], ystyle=1, xr=[1,1.3], xstyle=1
loadct,ctbl
for ig=0,Ngroups-1 do begin
   ; overplot Avg trend of group ig as a thick-line.
   oplot,rad_fit_euvia_A,Ne_avg(ig,*),color=colors(ig),th=4
   ; overplot Median line of group ig, data as diamonds and fit as a thin-line.
   ifl_A    = where(line_groupID eq ig AND fitflag_EUVIA_A eq +1. AND scN_fit_EUVIA_A lt chsq_threshold AND scT_fit_EUVIA_A lt chsq_threshold)
   ifl      = median(ifl_A)
   tmp      = reform(index_sampling_euvia_A(ifl,*))
   ind_samp = where(tmp eq 1)
   oplot,rad_A(ifl,ind_samp),Ne_euvia_A(ifl,ind_samp),color=colors(ig),psym=4 ; DEMT data points
   oplot,rad_fit_euvia_A,Ne_fit_euvia_A(ifl,*),color=colors(ig)               ; fit
endfor

loadct,0
MK = 1.e6 ; K
plot,rad_fit_euvia_A,Tm_avg(0,*)/MK,charsize=2,xtitle='r [Rsun]',title='EUVIA-DEMT Tm(r) [MK]',$
     psym=4,th=4, /nodata, yr=[0.,2.], ystyle=1, xr=[1,1.3], xstyle=1
loadct,ctbl
for ig=0,Ngroups-1 do begin
   ; overplot Avg trend of group ig as a thick-line.
   oplot,rad_fit_euvia_A,Tm_avg(ig,*)/MK,color=colors(ig),th=4
   ; overplot Median line of group ig, data as diamonds and fit as a thin-line.
   ifl_A    = where(line_groupID eq ig AND fitflag_EUVIA_A eq +1. AND scN_fit_EUVIA_A lt chsq_threshold AND scT_fit_EUVIA_A lt chsq_threshold)
   ifl      = median(ifl_A)
   tmp      = reform(index_sampling_euvia_A(ifl,*))
   ind_samp = where(tmp eq 1)
   oplot,rad_A(ifl,ind_samp),Tm_euvia_A(ifl,ind_samp)/MK,color=colors(ig),psym=4 ; DEMT data points
   oplot,rad_fit_euvia_A,Tm_fit_euvia_A(ifl,*)/MK,color=colors(ig)               ; fit
endfor
loadct,0


stop
  return
end


pro record_jpg,dir,filename

; set graph stuff
device, retain     = 2
device, true_color = 24
device, decomposed = 0

image24 = TVRD(True=1)
image2d = Color_Quan(image24, 1, r, g, b)
write_jpeg,dir+filename,image24,quality=100,true=1
return
end
