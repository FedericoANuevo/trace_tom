pro plot_lines_call, aia = aia, euvia = euvia, euvib = euvib, eit = eit, $
                     mk4 = mk4, kcor = kcor, lascoc2 = lascoc2, $
                     log = log

  structure_filename = 'list_synth.txt-tracing-structure-merge_aia_mk4_lascoc2_aia_mk4_lascoc2_aia_mk4_lascoc2_aia_mk4_lascoc2.sav'
  
  plot_lines, dir=dir, structure_filename=structure_filename, $
              aia = aia, euvia = euvia, euvib = euvib, eit = eit, $
              mk4 = mk4, kcor = kcor, lascoc2 = lascoc2, $
              log = log
  
  return
end
