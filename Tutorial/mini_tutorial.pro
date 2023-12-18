pro mini_tutorial
  common data, N_fl, Npt_max, x_A, y_A, z_A, rad_A, lat_A, lon_A,$
     Ne_aia_A, Tm_aia_A, index_aia_A, index_sampling_aia_A,$
     Ne_euvia_A, Tm_euvia_A, index_euvia_A, index_sampling_euvia_A,$
     Ne_euvib_A, Tm_euvib_A, index_euvib_A, index_sampling_euvib_A,$
     Ne_eit_A, Tm_eit_A, index_eit_A, index_sampling_eit_A,$
     Ne_mk4_A, index_mk4_A, index_sampling_mk4_A,$
     Ne_kcor_A, index_kcor_A, index_sampling_kcor_A,$
     Ne_c2_A, index_c2_A, index_sampling_c2_A

; 1) Declare the DIR where the structure is located, and the filename.

dir = './'
structure_filename = 'list_synth.txt-tracing-structure-merge_aia_mk4_lascoc2_aia_mk4_lascoc2_aia_mk4_lascoc2_aia_mk4_lascoc2.sav'

; 2) Load structure into memory and extract all available arrays from it.

load_traced_data_structure, dir=dir, structure_filename=structure_filename, trace_data=trace_data,$
                            /aia, /mk4, /lascoc2

; 3) See the full contents of the structure.
print
print, 'Bellow is the list of the full contents of the structure, named "trace_data", loaded into memory by the command:'
print
print, 'load_traced_data_structure, dir=dir, structure_filename=structure_filename, trace_data=trace_data, /aia, /mk4, /lascoc2'
print
print, 'Everything that says *NullPointer* does not exist. For this particular structure only AIA, MK4 and C2 results exist.'
print
print,'-------------------------------'
print,'help, trace_data, /str'
print
help, trace_data, /str
print,'-------------------------------'
print

print, 'Press SPACE BAR to continue.'
pause

; 4) See the arrays extracted from the structure.
print
print, 'Let us see next the output of the command "help" for all variables contained in the header for which there is data:'
print
print,'-------------------------------'
help,N_fl, Npt_max, x_A, y_A, z_A, rad_A, lat_A, lon_A,$
     Ne_aia_A, Tm_aia_A, index_aia_A, index_sampling_aia_A,$
     Ne_mk4_A, index_mk4_A, index_sampling_mk4_A,$
     Ne_c2_A, index_c2_A, index_sampling_c2_A
print,'-------------------------------'
print


; 5) Some explanations follow.
print,'In this case there are "N_fl=4" field lines that were traced between rmin=1 Rs and rmax=10 Rs. These are synthetic radial field lines with a uniform radial step of 0.0015 Rs. For each field line a maximum of "Npt_max=15000" points is allowed. This number will be optimized later on, depending on what we find to be a typical number of points for Judit field lines. For now we just set it as very large number.'
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
ifl=0
tmp = reform(index_sampling_aia_A(ifl,*))
ind_samp_aia = where(tmp eq 1)
window,0
plot,rad_A(ifl,ind_samp_aia),Ne_aia_A(ifl,ind_samp_aia),charsize=2,xtitle='r [Rsun]',title='AIA-DEMT Ne(r) [cm!U-3!N]'

print
print, 'Press SPACE BAR to continue.'
pause
print
print,'Similarly, to plot the c2 results for the SAME field line, one does:'
print
print,'         ifl = 0'
print,'         tmp = reform(index_sampling_c2_A(ifl,*))'
print,' ind_samp_c2 = where(tmp eq 1)'
print,' window, 1'
print,' plot,rad_A(ifl,ind_samp_c2),Ne_c2_A(ifl,ind_samp_c2)'
print, 'Press SPACE BAR to see the plot.'
pause
ifl=0
tmp = reform(index_sampling_c2_A(ifl,*))
ind_samp_c2 = where(tmp eq 1)
window,1
plot,rad_A(ifl,ind_samp_c2),Ne_c2_A(ifl,ind_samp_c2),charsize=2,xtitle='r [Rsun]',title='C2-SRT Ne(r) [cm!U-3!N]'


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
