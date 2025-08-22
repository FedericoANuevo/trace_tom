; Calling sequence examples:
  Ne3_analysis,/aia, /kcor, /ucomp, /closed, plot_filename_suffix='EQ-STR',/positparam, /load,latlimits=[-50.,50.],/histo,tit='Streamer'
  Ne3_analysis,/aia, /kcor, /open, plot_filename_suffix='NCH',/positparam,/load, r_max=1.195,latlimits=[ 50.,90.],/histo,tit='North-CH' 
  Ne3_analysis,/aia, /kcor, /open, plot_filename_suffix='SCH',/positparam,       r_max=1.195,latlimits=[-90.,-50],/histo,tit='South-CH' 
