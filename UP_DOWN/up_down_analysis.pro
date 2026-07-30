
pro up_down_analysis, LonLimits=LonLimits, LatLimits=LatLimits, $
                      plot_filename_suffix=plot_filename_suffix,$
                      aia=aia, euvi=euvi,$
                      open=open,closed=closed,$
                      positparam=positparam,$
                      only_loops=only_loops, loadstruct=loadstruct,$
                      histo = histo, tit = tit

  
  common datastructure, trace_data
  
  
; Create custom made symbol (psym=8) for scatter plots
  N=25
  A = FINDGEN(N) * (!PI*2/float(N-1))
  f=2.
  USERSYM, COS(A)/f, SIN(A)/f,/FILL 
;===============================================================================================
; Select the project to analyze:
  PROJECT_NAME = 'CR2223'
; Select structure to read:
  structure_filename='fdips_field_150X180X360_hmi.Synoptic_Mr_polfil.2223_prep.ubdat_fline-filenames_list.txt-tracing-structure-merge_aia_euvia_sampled.sav'
; Select dir where the structure is located (labeled after the selection of starting points)
  field_line_geometry_suffix_dir = '_aunifgrid_multirad-6h_2x2deg_HMI-PolFil/'
;===============================================================================================
  
  dir = '/data1/DATA/trace_tom_files/'+PROJECT_NAME+'/field_lines_geometry'+field_line_geometry_suffix_dir
  
; Load structure if so requested:
  if keyword_Set(loadstruct) then begin
     PRINT,'RESTORING THE POINTER-STRUCTURE'  
     restore, FILENAME = dir + structure_filename
     PRINT,'RESTORE COMPLETED'
  endif
  
  N_fl = *trace_data.N_fl


  
  midpoint_rad = fltarr(N_fl) -678.
  midpoint_lat = fltarr(N_fl) -678.
  midpoint_lon = fltarr(N_fl) -678.
; Calculate mid-point (~ 1.105 Rsun) of fieldlines  
  for ifl = 0L, N_fl-1 do begin
     Nsamp   = (*trace_data.Npt_aia)(ifl)
     if Nsamp ne -678. then begin
        radsamp = (*trace_data.rad_aia)(ifl,0:Nsamp-1)
        latsamp = (*trace_data.lat_aia)(ifl,0:Nsamp-1)
        lonsamp = (*trace_data.lon_aia)(ifl,0:Nsamp-1)
        index = where(radsamp gt 1.10 and radsamp lt 1.11)
        if index(0) ne -1 then begin
           midpoint_rad(ifl) = radsamp(index)
           midpoint_lat(ifl) = latsamp(index)
           midpoint_lon(ifl) = lonsamp(index)
        endif
     endif
  endfor
  

  
; Define plot filename suffix
  if not keyword_set(plot_filename_suffix) then plot_filename_suffix='footpoints-map'
; Default index line step and connect step, for plots only, this does


; Define box of footpoint Lat/Lot to analyze
  if not keyword_set(LonLimits) then LonLimits = [  0., 360.]
  if not keyword_set(LatLimits) then LatLimits = [-90.,+ 90.]

; Define a polarity tag
  tag_polarity_A = reform((*trace_data.leg_footBfield)(*,0)/abs((*trace_data.leg_footBfield)(*,0)))
  if ( where(finite(tag_polarity_A) eq 0) )(0) ne -1 then STOP ; assumes foot-Br is not ZERO
  
; Tag field lines with footpoints within the BOX., either open or closed.
  tag_box_A = fltarr(N_fl)
  index_box = where( (*trace_data.Footpoint_Lon ge LonLimits(0) AND *trace_data.Footpoint_Lon le LonLimits(1)) AND $
                       (*trace_data.Footpoint_Lat ge LatLimits(0) AND *trace_data.Footpoint_Lat le LatLimits(1)) )
  tag_box_A(index_box) = +1.

; Now, UNTAG CLOSED field lines whose other leg's footpoint is NOT wihin the BOX,
; as well as CLOSED field lines that do NOT comply with opposite polarity.
  if keyword_Set(only_loops) then begin
     ifl=0L
     while ifl le N_fl-2 do begin
        if (*trace_data.leg_label)(ifl) eq 0. then begin
           ifl=ifl+1
        endif else begin
           if (*trace_data.leg_label)(ifl) ne (*trace_data.leg_label)(ifl+1) then STOP ; !this should never happen.
           if (tag_box_A(ifl) eq +1 AND tag_box_A(ifl+1) eq  0) then tag_box_A(ifl  )=0
           if (tag_box_A(ifl) eq  0 AND tag_box_A(ifl+1) eq +1) then tag_box_A(ifl+1)=0
           if (tag_polarity_A(ifl) eq tag_polarity_A(ifl+1))    then tag_box_A(ifl:ifl+1)=0
           ifl=ifl+2
        endelse
     endwhile
  endif
    
; Create index arrays for closed and open field lines.
  ifl_closed_A = where(*trace_data.leg_label ne 0)
  ifl_open_A   = where(*trace_data.leg_label eq 0)
 
; Create index arrays for positive/negative footpoint Brad lines
  ifl_pos_A = where(tag_polarity_A eq +1.)
  ifl_neg_A = where(tag_polarity_A eq -1.)

; Set a maximum threshold for scN
  scN_crit = 0.10
  scT_crit = 0.10
  IF keyword_set(aia) then begin
; Independently for each instrument, tag field lines for which their
; DPL Ne-fit has all parameters positive     
     tag_posfit_aia_A = fltarr(N_Fl)
     index            = where(*trace_data.N1_fit_aia gt 0. AND *trace_data.N2_fit_aia gt 0. AND *trace_data.p1_fit_aia gt 0. AND *trace_data.p2_fit_aia gt 0.)
     tag_posfit_aia_A(index) = +1
; Independently for each instrument, index lines that:
; 1) are in the box, 2) have a fit, and 3) have a low chisq fit to Ne and Tm
     ifl_aia_A   = where(tag_box_A eq +1 AND *trace_data.fitflag_aia  eq +1 AND *trace_data.scN_fit_aia  le scN_crit AND *trace_data.scT_fit_aia  le scT_crit)
; Also independently for each instrument,
; filter OUT lines for which NOT all their DPL Ne-fit parameters are positive
     if keyword_set(positparam) then begin
        indpositparam_aia   = where(tag_posfit_aia_A   eq +1)  &  ifl_aia_A   = intersect(ifl_aia_A  ,indpositparam_aia  ) 
     endif
     tag_gradT_aia_A   = fltarr(N_Fl)
     index          = where(*trace_data.fitflag_aia  eq +1 and *trace_data.dTdr_fit_aia gt 0.)
     tag_gradT_aia_A(index)   = +1
     index          = where(*trace_data.fitflag_aia  eq +1 and *trace_data.dTdr_fit_aia lt 0.)
     tag_gradT_aia_A(index)   = -1
     if keyword_Set(only_loops) then begin
        ifl=0L
        while ifl le N_fl-2 do begin
           if (*trace_data.leg_label)(ifl) eq 0. then begin
              ifl=ifl+1
           endif else begin
              IF (*trace_data.fitflag_aia)(ifl) eq +1 and (*trace_data.fitflag_aia)(ifl+1) AND $
                 (*trace_data.dTdr_fit_aia)(ifl)*(*trace_data.dTdr_fit_aia)(ifl+1) lt 0.  then $
                    tag_gradT_aia_A(ifl:ifl+1)=0
              ifl=ifl+2
           endelse
        endwhile
     endif
  endif

  IF keyword_set(euvi) then begin
; Independently for each instrument, tag field lines for which their
; DPL Ne-fit has all parameters positive     
     tag_posfit_euvia_A = fltarr(N_Fl)
     index            = where(*trace_data.N1_fit_euvia gt 0. AND *trace_data.N2_fit_euvia gt 0. AND *trace_data.p1_fit_euvia gt 0. AND *trace_data.p2_fit_euvia gt 0.)
     tag_posfit_euvia_A(index) = +1
; Independently for each instrument, index lines that:
; 1) are in the box, 2) have a fit, and 3) have a low chisq fit to Ne and Tm
     ifl_euvi_A   = where(tag_box_A eq +1 AND *trace_data.fitflag_euvia  eq +1 AND *trace_data.scN_fit_euvia  le scN_crit AND *trace_data.scT_fit_euvia  le scT_crit)
; Also independently for each instrument,
; filter OUT lines for which NOT all their DPL Ne-fit parameters are positive
     if keyword_set(positparam) then begin
        indpositparam_euvia   = where(tag_posfit_euvia_A   eq +1)  &  ifl_euvi_A   = intersect(ifl_euvi_A  ,indpositparam_euvia  ) 
     endif
     tag_gradT_euvi_A   = fltarr(N_Fl)
     index          = where(*trace_data.fitflag_euvia  eq +1 and *trace_data.dTdr_fit_euvia gt 0.)
     tag_gradT_euvi_A(index)   = +1
     index          = where(*trace_data.fitflag_euvia  eq +1 and *trace_data.dTdr_fit_euvia lt 0.)
     tag_gradT_euvi_A(index)   = -1
     if keyword_Set(only_loops) then begin
        ifl=0L
        while ifl le N_fl-2 do begin
           if (*trace_data.leg_label)(ifl) eq 0. then begin
              ifl=ifl+1
           endif else begin
              IF (*trace_data.fitflag_euvia)(ifl) eq +1 and (*trace_data.fitflag_euvia)(ifl+1) AND $
                 (*trace_data.dTdr_fit_euvia)(ifl)*(*trace_data.dTdr_fit_euvia)(ifl+1) lt 0.  then $
                    tag_gradT_euvi_A(ifl:ifl+1)=0
              ifl=ifl+2
           endelse
        endwhile
     endif
  endif


  
;-----------------PLOTS SECTION--------------------------------------------------------------
; Define a few color codes.
  blue  = 100
  red   = 200
  green =  16

; Lat/Lon plots of FootPoints
  fig_dir = './'
; fig_dir = '~/Downloads/'
; fig_dir = dir 
  ps1,fig_dir+structure_filename+'_'+plot_filename_suffix+'.eps'
  np=1000
  !p.multi=[0,1,2]
  loadct,0
  !p.color=0
  !p.background=255
  csz=1
  


; First plot
  title = 'AIA: Physical location of leg at 1.105 Rsun'
  plot,*trace_data.Footpoint_lon,*trace_data.Footpoint_lat,xr=[0,360],yr=[-90,+90],xstyle=1,ystyle=1,/nodata,charsize=csz,$
       title=title,ytitle='Carrington Latitude [deg]',xtitle='Carrington Longitude [deg]',font=0
  if not keyword_set(open) and not keyword_set(closed) then oplot,midpoint_Lon              ,midpoint_Lat              ,psym=3
  if     keyword_set(open)                             then oplot,midpoint_Lon(ifl_open_A)  ,midpoint_Lat(ifl_open_A)  ,psym=3
  if                               keyword_set(closed) then oplot,midpoint_Lon(ifl_closed_A),midpoint_Lat(ifl_closed_A),psym=3

; Make an index array with the common elements of all ifl_INSTRUMENT_A
  ifl_A = indgen(N_fl,/LONG)    ; start index with ALL lines
  if keyword_set(aia)    then ifl_A = intersect(ifl_A,ifl_aia_A  )
; Intersect with CLOSED or OPEN, if so requested
  if keyword_set(closed) then ifl_A = intersect(ifl_A,ifl_closed_A)
  if keyword_set(open)   then ifl_A = intersect(ifl_A,ifl_open_A  )
; Create index arrays for positive/negative footpoint Brad lines
  ifl_pos_A = where(tag_gradT_aia_A   eq +1.)
  ifl_neg_A = where(tag_gradT_aia_A   eq -1.)

  loadct,12
; Color-highlight all footpoints indicated by ifl_A accordind to their
; polarity / or UP AND DOWN
  indxpos_A = intersect(ifl_A,ifl_pos_A)
  indxneg_A = intersect(ifl_A,ifl_neg_A)
  if n_elements(indxpos_A) gt 1 then $
     oplot,midpoint_lon(indxpos_A),midpoint_lat(indxpos_A),psym=3,th=2,color=red
  if n_elements(indxneg_A) gt 1 then $
     oplot,midpoint_lon(indxneg_A),midpoint_lat(indxneg_A),psym=3,th=2,color=blue 
  loadct,0
  print,'AIA'
  print,'fraction of Up   loops:',n_elements(indxpos_A)*1./( n_elements(indxpos_A)+n_elements(indxneg_A))
  print,'fraction of Down loops:',n_elements(indxneg_A)*1./( n_elements(indxpos_A)+n_elements(indxneg_A))

; Calculate a 1D vector with dTm/dr for histogram  
  index = where(abs(tag_gradT_aia_A) eq 1 and abs(*trace_data.dTdr_fit_aia)/1.e6 gt 0.5)
  index = intersect(index,ifl_A)
  gradT_aia = (*trace_data.dTdr_fit_aia)(index)/1.e6
  
  ; Second plot
  title = 'EUVI: Physical location of leg at 1.105 Rsun'
  plot,*trace_data.Footpoint_lon,*trace_data.Footpoint_lat,xr=[0,360],yr=[-90,+90],xstyle=1,ystyle=1,/nodata,charsize=csz,$
       title=title,ytitle='Carrington Latitude [deg]',xtitle='Carrington Longitude [deg]',font=0
  if not keyword_set(open) and not keyword_set(closed) then oplot,midpoint_Lon              ,midpoint_Lat              ,psym=3
  if     keyword_set(open)                             then oplot,midpoint_Lon(ifl_open_A)  ,midpoint_Lat(ifl_open_A)  ,psym=3
  if                               keyword_set(closed) then oplot,midpoint_Lon(ifl_closed_A),midpoint_Lat(ifl_closed_A),psym=3

; Make an index array with the common elements of all ifl_INSTRUMENT_A
  ifl_A = indgen(N_fl,/LONG)    ; start index with ALL lines
  if keyword_set(euvi)    then ifl_A = intersect(ifl_A,ifl_euvi_A  )
; Intersect with CLOSED or OPEN, if so requested
  if keyword_set(closed) then ifl_A = intersect(ifl_A,ifl_closed_A)
  if keyword_set(open)   then ifl_A = intersect(ifl_A,ifl_open_A  )
; Create index arrays for positive/negative footpoint Brad lines
  ifl_pos_A = where(tag_gradT_euvi_A   eq +1.)
  ifl_neg_A = where(tag_gradT_euvi_A   eq -1.)

  loadct,12
; Color-highlight all footpoints indicated by ifl_A accordind to their
; polarity / or UP AND DOWN
  indxpos_A = intersect(ifl_A,ifl_pos_A)
  indxneg_A = intersect(ifl_A,ifl_neg_A)
  if n_elements(indxpos_A) gt 1 then $
     oplot,midpoint_lon(indxpos_A),midpoint_lat(indxpos_A),psym=3,th=2,color=red
  if n_elements(indxneg_A) gt 1 then $
     oplot,midpoint_lon(indxneg_A),midpoint_lat(indxneg_A),psym=3,th=2,color=blue 
  loadct,0
  print,'EUVI'
  print,'fraction of Up   loops:',n_elements(indxpos_A)*1./( n_elements(indxpos_A)+n_elements(indxneg_A))
  print,'fraction of Down loops:',n_elements(indxneg_A)*1./( n_elements(indxpos_A)+n_elements(indxneg_A))
  !p.multi=0
; Calculate a 1D vector with dTm/dr for histogram    
  index = where(abs(tag_gradT_euvi_A) eq 1 and abs(*trace_data.dTdr_fit_euvia)/1.e6 gt 0.5)
  index = intersect(index,ifl_A)
  gradT_euvi = (*trace_data.dTdr_fit_euvia)(index)/1.e6

  

; Histograms of dTm/dr  
  xhisto2,gradT_aia ,comp_suffix='GradT_aia' ,sufijo=plot_filename_suffix, histo_x_tit='dTm/dr [MK/Rsun]',Nvals =50, dir_fig ='./',mini=-8.,maxi=8.,tit='AIA'
  xhisto2,gradT_euvi,comp_suffix='GradT_euvi',sufijo=plot_filename_suffix, histo_x_tit='dTm/dr [MK/Rsun]',Nvals =50, dir_fig ='./',mini=-8.,maxi=8.,tit='EUVI-A'
  
  PS2
  stop
  return
end

PRO ps1,archivo
set_plot,'ps'
device,filename=archivo,bits_per_pixel=8,/color,/encapsulated
return
end

PRO ps2
device,/close
set_plot,'x'
!p.multi=0
return
end


pro xhisto2,vector,comp_suffix=comp_suffix,sufijo=sufijo,tit=tit,histo_x_tit=histo_x_tit,Nvals=Nvals,cleanstat=cleanstat,dir_fig=dir_fig,mini=mini,maxi=maxi

  if not keyword_set(dir_fig) then dir_fig = '/data1/tomography/SolarTom_idl/Figures/'


  if not keyword_set(mini) then mini = min(vector)
  if not keyword_set(maxi) then maxi = max(vector)

  vector = vector > mini < maxi
  
;  mini      = min(vector)
;  maxi      = max(vector)
  delta     = (maxi-mini)/Nvals

  histo_vector = histogram(vector,binsize=delta,locations=xval)
  histo_vector = float(histo_vector) / float(n_elements(vector))

  index = where( vector gt mini + delta and vector lt maxi)
  v_cut = vector(index)

  ; i1= where( vector ge mini and vector le mini + delta)
  ; i2= where( vector ge maxi and vector le maxi + delta)
  ; stop
  
  if not keyword_set(cleanstat) then begin
     avg        =   mean(vector)
     med        = median(vector)
     stdv       =  stdev(vector)
     stdev_frac =  stdev(vector)/abs(avg)
  endif
  
  if keyword_set(cleanstat) then begin
     avg        =   mean(v_cut)
     med        = median(v_cut)
     stdv       =  stdev(v_cut)
     stdev_frac =  stdev(v_cut)/abs(avg) 
  endif

  print,'uncut - cut'
  print,'mean',mean(vector), mean(v_cut)
  print,'median',median(vector),median(v_cut)
  print,'st. dev.',stdev(vector),stdev(v_cut)

  
; redondeo de med y stdv a dos cifras decimales
  f    = 100.
  med  = round(float(med)*f)/f
  stdv = round(float(stdv)*f)/f

  if med ge 1. then  med_str =strmid(string(med),6,4)
  if med lt 1. then  med_str =strmid(string(med),5,4)
  stdv_str=strmid(string(stdv),5,4)
    
  cant       = long(n_elements(vector))
  name_fig   = 'comparison_'+comp_suffix+'_'+sufijo
  name_fig   =  (STRJOIN(STRSPLIT(name_fig, /EXTRACT,'.'), '_'))
  ps1,dir_fig+name_fig+'.eps'

  plot,xval,histo_vector ,font=0,xtitle=histo_x_tit,title=tit,linestyle=8,psym=10,thick=4,charsize=2.2,xrange=[mini-0.1,maxi+0.1],xstyle=1,yrange=[0.,0.16],ystyle=1
; xyouts,0.7*[1,1,1,1],0.98-[0.18,0.25,0.32,0.38],['m='+med_str,'','!9s!3='+stdv_str,''],/normal,charthick=1,Font=0,charsize=2.2
  !p.multi = 0
     
  return
end
