;
; PURPOSE: Performs analytical fits to tomographic products along each
; traced field line and stores results in the trace_data structure.
;
; HISTORY: V1.0, AMV, January 2024, IAFE.
;          V1.1, AMV & FAN, January 2024, IAFE. Added Theil-Sen and R² metric.
;          V1.2, AMV, January 2024, IAFE. Added Double-Power-Law for Mk4.
;                                         Added fits' parameters to output structure.
;                                         Added 2-PoweLaw fit for C2.
;          V1.3, AMV, January 2024, IAFE. Added DPL fit for EUV.
;                                         Added fitted function name
;                                         and automated parameter
;                                         storing in output structure.
;          V1.4, AMV, January 2024, IAFE. Added <lambda_N> for all
;                                         functional fits.
;          V1.4.1, FAN, April 2024, IAFE. Added goto,skip_fit_aia
;          and if indsamp[0] ne -1 then ...
;          v1.4.1, FAN, May 2024, ClaSP. Change GOTO by IF
;          v1.4.1, FAN, May 2024, ClaSP. Change scN_fit and scT_fit as
;          a function of chisqr
;          v1.4.2, FAN, May 2024, ClaSP. euvia keyword added.


pro fit_trace_data, aia=aia, euvia=euvia, euvib=euvib, eit=eit,$
                    mk4=mk4, kcor=kcor, ucomp=ucomp, lascoc2=lascoc2,$
                    fl_dir=fl_dir

    common to_fit_data, trace_data, $
     N_fl, Npt_max, Npt_v,$
     x_A, y_A, z_A, rad_A, lat_A, lon_A,$
     Ne_aia_A, Tm_aia_A, index_aia_A, index_sampling_aia_A,$
     Ne_euvia_A, Tm_euvia_A, index_euvia_A, index_sampling_euvia_A,$
     Ne_euvib_A, Tm_euvib_A, index_euvib_A, index_sampling_euvib_A,$
     Ne_eit_A, Tm_eit_A, index_eit_A, index_sampling_eit_A,$
     Ne_mk4_A, index_mk4_A, index_sampling_mk4_A,$
     Ne_kcor_A, index_kcor_A, index_sampling_kcor_A,$
     Ne_ucomp_A, index_ucomp_A, index_sampling_ucomp_A,$
     Ne_c2_A, index_c2_A, index_sampling_c2_A
    common radcrits, radcritA, radcritB
    common tomgrid,nr,nt,np,rmin,rmax,Irmin,Irmax
    common fitgrid,radmin_fit,radmax_fit,drad_fit,Npt_fit

   ;Set a default value for all the elements of all the [NAME]_A arrays
    default = -678.
    
    if keyword_set(aia) then begin
       read_tomgrid_and_define_fitgrid,fl_dir=fl_dir,instr_string='aia'
       fitflag_aia_A = fltarr(N_fl        ) + default
        N0_fit_aia_A = fltarr(N_fl        ) + default
        lN_fit_aia_A = fltarr(N_fl        ) + default
        N1_fit_aia_A = fltarr(N_fl        ) + default
        N2_fit_aia_A = fltarr(N_fl        ) + default
        p1_fit_aia_A = fltarr(N_fl        ) + default
        p2_fit_aia_A = fltarr(N_fl        ) + default
       scN_fit_aia_A = fltarr(N_fl        ) + default
        T0_fit_aia_A = fltarr(N_fl        ) + default
      dTdr_fit_aia_A = fltarr(N_fl        ) + default
       scT_fit_aia_A = fltarr(N_fl        ) + default
        Ne_fit_aia_A = fltarr(N_fl,Npt_fit) + default
        Tm_fit_aia_A = fltarr(N_fl,Npt_fit) + default
        rad_fit_aia_A = radmin_fit + drad_fit/2. + drad_fit * findgen(Npt_fit)
             dNe_dr_A = fltarr(Npt_fit) + default
       for ifl=0,N_fl-1 do begin
          tmp = reform(index_sampling_aia_A(ifl,*))
          ind_samp_aia = where(tmp eq 1)
          if ind_samp_aia[0] ne -1 then begin
             radsamp = reform(rad_A(ifl,ind_samp_aia)) ; Rsun
            ;Determine radsamp_max (which is the apex in the case of small closed loops,
            ;                       or even a smaller height if the apex is a ZDA.
             radsamp_max=max(radsamp)
             Nradsamp   =n_elements(radsamp)
            ;Sanity check of radsamp_max
             if radsamp(0) lt radsamp(Nradsamp-1) and radsamp_max ne radsamp(Nradsamp-1) then STOP
             if radsamp(0) gt radsamp(Nradsamp-1) and radsamp_max ne radsamp(0)          then STOP
            ;Determine the min and max rad over which we will actually evaluate the fit.
             radfit_min =                  min(rad_fit_aia_A)
             radfit_max = min([radsamp_max,max(rad_fit_aia_A)]) 
            ;Determine range of rad_fit_[instrument] over which we will actually evaluate the fit.
             range_fit = where(rad_fit_aia_A ge radfit_min AND rad_fit_aia_A le radfit_max)
            ;Test if there is proper coverage of actual sample data for a decent least squares fit.
             test_coverage, radsamp=radsamp, radfit_min=radfit_min, radfit_max=radfit_max, covgflag=covgflag, /aia
             if covgflag eq 'yes' then begin
                fitflag_aia_A(ifl) = +1.
                Nesamp = reform(Ne_aia_A(ifl,ind_samp_aia))
                Tmsamp = reform(Tm_aia_A(ifl,ind_samp_aia))
                goto,skip_aia_isohthermal_hydrostatic
                fit_F_Ne_aia  = 'IHS'
                linear_fit, 1./radsamp   , alog(Nesamp), AN, r2N, /linfit_idl
                scN_fit_aia_A(ifl)   = r2N             
                 N0_fit_aia_A(ifl)   = exp(AN[0]+AN[1])                                                                         ; cm-3
                 lN_fit_aia_A(ifl)   = 1./AN[1]                                                                                 ; Rsun
                 Ne_fit_aia_A(ifl,range_fit) = N0_fit_aia_A(ifl) * exp(-(1/lN_fit_aia_A(ifl))*(1.-1./rad_fit_aia_A(range_fit))) ; cm-3
                     dNe_dr_A    (range_fit) = reform(Ne_fit_aia_A(ifl,range_fit)) * float(-(1/lN_fit_aia_A(ifl))) / rad_fit_aia_A(range_fit)^2    ; cm-3 / Rsun
                 indsamp                     = where(dNe_dr_A lt 0. AND dNe_dr_A ne default)
                 v = abs(dNe_dr_A(indsamp)/reform(Ne_fit_aia_A(ifl,indsamp)))^(-1)
                 lN_fit_aia_A(ifl)   =  int_tabulated(rad_fit_aia_A(indsamp),v) / (max(rad_fit_aia_A(indsamp))-min(rad_fit_aia_A(indsamp))) ; cm-3 / Rsun
                 print,lN_fit_aia_A(ifl), float(mean(v)), float(median(v)), float(1./AN[1])
                skip_aia_isohthermal_hydrostatic:
               ;goto,skip_aia_double_power_law
                fit_F_Ne_aia  = 'DPL'
                double_power_fit, radsamp, Nesamp, A, chisq ;, /weighted
                scN_fit_aia_A(ifl)  = sqrt(chisq)/mean(Nesamp)
                 N1_fit_aia_A(ifl)   = A[0]                                                                          ; cm-3
                 p1_fit_aia_A(ifl)   = A[1]                                                                          ; dimensionless exponent of power law
                 N2_fit_aia_A(ifl)   = A[2]                                                                          ; cm-3
                 p2_fit_aia_A(ifl)   = A[3]                                                                          ; dimensionless exponent of power law
                 Ne_fit_aia_A(ifl,range_fit) = A[0] * rad_fit_aia_A(range_fit)^(-A[1]) + A[2] * rad_fit_aia_A(range_fit)^(-A[3]) ; cm-3
                     dNe_dr_A    (range_fit) = - A[1]*A[0] * rad_fit_aia_A(range_fit)^(-A[1]-1) - A[3]*A[2] * rad_fit_aia_A(range_fit)^(-A[3]-1) ; cm-3 / Rsun
                 indsamp                     = where(dNe_dr_A lt 0. AND dNe_dr_A ne default)
                if indsamp[0] ne -1 then begin
                   v = abs(dNe_dr_A(indsamp)/reform(Ne_fit_aia_A(ifl,indsamp)))^(-1)
                   lN_fit_aia_A(ifl)   =  int_tabulated(rad_fit_aia_A(indsamp),v) / (max(rad_fit_aia_A(indsamp))-min(rad_fit_aia_A(indsamp)))
                   ; print,lN_fit_aia_A(ifl), float(mean(v)), float(median(v))
                   ; stop
                endif
                skip_aia_double_power_law: 
                ;Linear fit to Te(r)
                linear_fit,    radsamp-1.,      Tmsamp , AT, r2T, /theil_sen, chisqr = chisqr
                 scT_fit_aia_A(ifl) = sqrt(chisqr)/mean(Tmsamp)                                                                    
                  T0_fit_aia_A(ifl) = AT[0]                                                                             ; K
                dTdr_fit_aia_A(ifl) = AT[1]                                                                             ; K/Rsun
                  Tm_fit_aia_A(ifl,range_fit) = T0_fit_aia_A(ifl) + dTdr_fit_aia_A(ifl) * (rad_fit_aia_A(range_fit)-1.) ; K
             endif ; covgflag = 'yes'
          endif
       endfor                   ; field lines loop.
       trace_data = create_struct( trace_data                           ,$
                                  'fitflag_aia',ptr_new( fitflag_aia_A) ,$
                                 'fit_F_Ne_aia',ptr_new(  fit_F_Ne_aia) ,$
                                  'scN_fit_aia',ptr_new( scN_fit_aia_A) ,$
                                  'scT_fit_aia',ptr_new( scT_fit_aia_A) ,$
                                  'rad_fit_aia',ptr_new( rad_fit_aia_A) ,$
                                   'Ne_fit_aia',ptr_new(  Ne_fit_aia_A) ,$
                                   'Tm_fit_aia',ptr_new(  Tm_fit_aia_A) ,$
                                   'T0_fit_aia',ptr_new(  T0_fit_aia_A) ,$
                                 'dTdr_fit_aia',ptr_new(dTdr_fit_aia_A) ,$
                                   'lN_fit_aia',ptr_new(  lN_fit_aia_A) )
       if fit_F_Ne_aia eq 'IHS' then $
          trace_data = create_struct( trace_data                        ,$
                                   'N0_fit_aia',ptr_new(  N0_fit_aia_A) )                                                                        
       if fit_F_Ne_aia eq 'DPL' then $
          trace_data = create_struct( trace_data                        ,$
                                    'N1_fit_aia',ptr_new( N1_fit_aia_A) ,$
                                    'N2_fit_aia',ptr_new( N2_fit_aia_A) ,$
                                    'p1_fit_aia',ptr_new( p1_fit_aia_A) ,$
                                    'p2_fit_aia',ptr_new( p2_fit_aia_A) )
       ;Add to structure the parameters of the tomography and fit grids.
       trace_data = create_struct( trace_data                                                      ,$
                                   'tomgrid_aia_hdr',ptr_new(['nr','nt','np','rmin','rmax','Irmin','Irmax'])   ,$
                                   'tomgrid_aia'    ,ptr_new([ nr , nt , np , rmin , rmax , Irmin , Irmax ])   ,$
                                   'fitgrid_aia_hdr',ptr_new(['radmin_fit','radmax_fit','drad_fit','Npt_fit']) ,$
                                   'fitgrid_aia'    ,ptr_new([ radmin_fit , radmax_fit , drad_fit , Npt_fit ]) )
    endif ; AIA

    if keyword_set(euvia) then begin
       read_tomgrid_and_define_fitgrid,fl_dir=fl_dir,instr_string='euvia'
       fitflag_euvia_A = fltarr(N_fl        ) + default
        N0_fit_euvia_A = fltarr(N_fl        ) + default
        lN_fit_euvia_A = fltarr(N_fl        ) + default
        N1_fit_euvia_A = fltarr(N_fl        ) + default
        N2_fit_euvia_A = fltarr(N_fl        ) + default
        p1_fit_euvia_A = fltarr(N_fl        ) + default
        p2_fit_euvia_A = fltarr(N_fl        ) + default
       scN_fit_euvia_A = fltarr(N_fl        ) + default
        T0_fit_euvia_A = fltarr(N_fl        ) + default
      dTdr_fit_euvia_A = fltarr(N_fl        ) + default
       scT_fit_euvia_A = fltarr(N_fl        ) + default
        Ne_fit_euvia_A = fltarr(N_fl,Npt_fit) + default
        Tm_fit_euvia_A = fltarr(N_fl,Npt_fit) + default
       rad_fit_euvia_A = radmin_fit + drad_fit/2. + drad_fit * findgen(Npt_fit)
       for ifl=0,N_fl-1 do begin      
          tmp = reform(index_sampling_euvia_A(ifl,*))
          ind_samp_euvia = where(tmp eq 1)
          if ind_samp_euvia[0] ne -1 then begin
             radsamp = reform(rad_A(ifl,ind_samp_euvia)) ; Rsun
             test_coverage, radsamp=radsamp, covgflag=covgflag, /euvia
             if covgflag eq 'yes' then begin
                fitflag_euvia_A(ifl) = +1.
                Nesamp = reform(Ne_euvia_A(ifl,ind_samp_euvia))
                Tmsamp = reform(Tm_euvia_A(ifl,ind_samp_euvia))
                goto,skip_euvia_isohthermal_hydrostatic
                fit_F_Ne_euvia  = 'IHS'
                linear_fit, 1./radsamp   , alog(Nesamp), AN, r2N, /linfit_idl
                scN_fit_euvia_A(ifl)  = r2N             
                N0_fit_euvia_A(ifl)   = exp(AN[0]+AN[1])                                                            ; cm-3
                lN_fit_euvia_A(ifl)   = 1./AN[1]                                                                    ; Rsun
                Ne_fit_euvia_A(ifl,*) = N0_fit_euvia_A(ifl) * exp(-(1/lN_fit_euvia_A(ifl))*(1.-1./rad_fit_euvia_A)) ; cm-3
                dNe_dr              = reform(Ne_fit_euvia_A(ifl,*)) * float(-(1/lN_fit_euvia_A(ifl))) / rad_fit_euvia_A^2 ; cm-3 / Rsun
                indsamp = where(rad_fit_euvia_A ge min(radsamp) and rad_fit_euvia_A le max(radsamp) AND dNe_dr lt 0.)
                v = abs(dNe_dr(indsamp)/reform(Ne_fit_euvia_A(ifl,indsamp)))^(-1)
                lN_fit_euvia_A(ifl)   =  int_tabulated(rad_fit_euvia_A(indsamp),v) / (max(rad_fit_euvia_A(indsamp))-min(rad_fit_euvia_A(indsamp))) ; cm-3 / Rsun
                print,lN_fit_euvia_A(ifl), float(mean(v)), float(median(v)), float(1./AN[1])
                skip_euvia_isohthermal_hydrostatic:
               ;goto,skip_euvia_double_power_law
                fit_F_Ne_euvia  = 'DPL'
                double_power_fit, radsamp, Nesamp, A, chisq ;, /weighted
                scN_fit_euvia_A(ifl)  = sqrt(chisq)/mean(Nesamp)
                N1_fit_euvia_A(ifl)   = A[0]                                                                          ; cm-3
                p1_fit_euvia_A(ifl)   = A[1]                                                                          ; dimensionless exponent of power law
                N2_fit_euvia_A(ifl)   = A[2]                                                                          ; cm-3
                p2_fit_euvia_A(ifl)   = A[3]                                                                          ; dimensionless exponent of power law
                Ne_fit_euvia_A(ifl,*) = A[0] * rad_fit_euvia_A^(-A[1]) + A[2] * rad_fit_euvia_A^(-A[3])                   ; cm-3
                dNe_dr              = - A[1]*A[0] * rad_fit_euvia_A^(-A[1]-1) - A[3]*A[2] * rad_fit_euvia_A^(-A[3]-1)   ; cm-3 / Rsun
                indsamp = where(rad_fit_euvia_A ge min(radsamp) and rad_fit_euvia_A le max(radsamp) AND dNe_dr lt 0.)
                if indsamp[0] ne -1 then begin
                   v = abs(dNe_dr(indsamp)/reform(Ne_fit_euvia_A(ifl,indsamp)))^(-1)
                   lN_fit_euvia_A(ifl)   =  int_tabulated(rad_fit_euvia_A(indsamp),v) / (max(rad_fit_euvia_A(indsamp))-min(rad_fit_euvia_A(indsamp)))
                   ; print,lN_fit_euvia_A(ifl), float(mean(v)), float(median(v))
                   ; stop
                endif
                skip_euvia_double_power_law: 
                ;Linear fit to Te(r)
                linear_fit,    radsamp-1.,      Tmsamp , AT, r2T, /theil_sen, chisqr = chisqr
                scT_fit_euvia_A(ifl)  = sqrt(chisqr)/mean(Tmsamp)                                                                    
                T0_fit_euvia_A(ifl)   = AT[0]                                                                         ; K
                dTdr_fit_euvia_A(ifl) = AT[1]                                                                         ; K/Rsun
                Tm_fit_euvia_A(ifl,*) = T0_fit_euvia_A(ifl) + dTdr_fit_euvia_A(ifl)       *      (rad_fit_euvia_A-1.) ; K
             endif                                                                                                    ; covgflag = 'yes'
          endif
       endfor                   ; field lines loop.
       trace_data = create_struct( trace_data                               ,$
                                  'fitflag_euvia',ptr_new( fitflag_euvia_A) ,$
                                 'fit_F_Ne_euvia',ptr_new(  fit_F_Ne_euvia) ,$
                                  'scN_fit_euvia',ptr_new( scN_fit_euvia_A) ,$
                                  'scT_fit_euvia',ptr_new( scT_fit_euvia_A) ,$
                                  'rad_fit_euvia',ptr_new( rad_fit_euvia_A) ,$
                                   'Ne_fit_euvia',ptr_new(  Ne_fit_euvia_A) ,$
                                   'Tm_fit_euvia',ptr_new(  Tm_fit_euvia_A) ,$
                                   'T0_fit_euvia',ptr_new(  T0_fit_euvia_A) ,$
                                 'dTdr_fit_euvia',ptr_new(dTdr_fit_euvia_A) ,$
                                   'lN_fit_euvia',ptr_new(  lN_fit_euvia_A) )
       if fit_F_Ne_euvia eq 'IHS' then $
          trace_data = create_struct( trace_data                        ,$
                                 'N0_fit_euvia',ptr_new(  N0_fit_euvia_A) )                                                                        
       if fit_F_Ne_euvia eq 'DPL' then $
          trace_data = create_struct( trace_data                        ,$
                                    'N1_fit_euvia',ptr_new( N1_fit_euvia_A) ,$
                                    'N2_fit_euvia',ptr_new( N2_fit_euvia_A) ,$
                                    'p1_fit_euvia',ptr_new( p1_fit_euvia_A) ,$
                                    'p2_fit_euvia',ptr_new( p2_fit_euvia_A) )

    endif ; EUVI-A

    if keyword_set(mk4) then begin
       read_tomgrid_and_define_fitgrid,fl_dir=fl_dir,instr_string='mk4'
       fitflag_mk4_A = fltarr(N_fl        ) + default
        N0_fit_mk4_A = fltarr(N_fl        ) + default
        lN_fit_mk4_A = fltarr(N_fl        ) + default
        N1_fit_mk4_A = fltarr(N_fl        ) + default
        N2_fit_mk4_A = fltarr(N_fl        ) + default
        p1_fit_mk4_A = fltarr(N_fl        ) + default
        p2_fit_mk4_A = fltarr(N_fl        ) + default
       scN_fit_mk4_A = fltarr(N_fl        ) + default
        Ne_fit_mk4_A = fltarr(N_fl,Npt_fit) + default
       rad_fit_mk4_A = radmin_fit + drad_fit/2. + drad_fit * findgen(Npt_fit)
       for ifl=0,N_fl-1 do begin      
          tmp = reform(index_sampling_mk4_A(ifl,*))
          ind_samp_mk4 = where(tmp eq 1)
          if ind_samp_mk4[0] ne -1 then begin
             radsamp = reform(rad_A(ifl,ind_samp_mk4)) ; Rsun
             test_coverage, radsamp=radsamp, covgflag=covgflag, /mk4
             if covgflag eq 'yes' then begin
                fitflag_mk4_A(ifl) = +1.
                Nesamp = reform(Ne_mk4_A(ifl,ind_samp_mk4))
                goto,skip_mk4_isohthermal_hydrostatic
                fit_F_Ne_mk4  = 'IHS'
                linear_fit, 1./radsamp   ,alog(Nesamp), AN, r2N, /linfit_idl
                scN_fit_mk4_A(ifl)   = r2N
                N0_fit_mk4_A(ifl)   = exp(AN[0]+AN[1])                                                                ; cm-3
                lN_fit_mk4_A(ifl)   = 1./AN[1]                                                                        ; Rsun
                Ne_fit_mk4_A(ifl,*) = N0_fit_mk4_A(ifl) * exp(-(1/lN_fit_mk4_A(ifl))*(1.-1./rad_fit_mk4_A))           ; cm-3
                dNe_dr              = reform(Ne_fit_mk4_A(ifl,*)) * float(-(1/lN_fit_mk4_A(ifl))) / rad_fit_mk4_A^2   ; cm-3 / Rsun
                indsamp = where(rad_fit_mk4_A ge min(radsamp) and rad_fit_mk4_A le max(radsamp) AND dNe_dr lt 0.)
                v = abs(dNe_dr(indsamp)/reform(Ne_fit_mk4_A(ifl,indsamp)))^(-1)
                lN_fit_mk4_A(ifl)   =  int_tabulated(rad_fit_mk4_A(indsamp),v) / (max(rad_fit_mk4_A(indsamp))-min(rad_fit_mk4_A(indsamp))) ; cm-3 / Rsun                  
                print,lN_fit_mk4_A(ifl), float(mean(v)), float(median(v)), float(1./AN[1])
                skip_mk4_isohthermal_hydrostatic:
                goto,skip_mk4_single_power_law
                fit_F_Ne_mk4  = 'SPL'
                linear_fit, alog(radsamp), alog(Nesamp), AN, r2N, /linfit_idl
                scN_fit_mk4_A(ifl)   = r2N
                N1_fit_mk4_A(ifl)   = exp(AN[0])                                               ; cm-3
                p1_fit_mk4_A(ifl)   =    -AN[1]                                                ; dimensionless exponent of power law
                Ne_fit_mk4_A(ifl,*) = N1_fit_mk4_A(ifl) * rad_fit_mk4_A^(-p1_fit_mk4_A(ifl))   ; cm-3
                indsamp = where(rad_fit_mk4_A ge min(radsamp) and rad_fit_mk4_A le max(radsamp)  AND dNe_dr lt 0.)
                v = abs(rad_fit_mk4_A(indsamp) / float(p1_fit_mk4_A(ifl))) 
                lN_fit_mk4_A(ifl)   = int_tabulated( rad_fit_mk4_A(indsamp), v) / (max(rad_fit_mk4_A(indsamp))-min(rad_fit_mk4_A(indsamp))) ; Rsun
                print,lN_fit_mk4_A(ifl), float(mean(v)), float(median(v))
                skip_mk4_single_power_law:
               ;goto,skip_mk4_double_power_law
                fit_F_Ne_mk4  = 'DPL'
                double_power_fit, radsamp, Nesamp, A, chisqr ;, /weighted
                scN_fit_mk4_A(ifl)   = sqrt(chisqr)/mean(Nesamp)
                N1_fit_mk4_A(ifl)   = A[0]                                                                          ; cm-3
                p1_fit_mk4_A(ifl)   = A[1]                                                                          ; dimensionless exponent of power law
                N2_fit_mk4_A(ifl)   = A[2]                                                                          ; cm-3
                p2_fit_mk4_A(ifl)   = A[3]                                                                          ; dimensionless exponent of power law
                Ne_fit_mk4_A(ifl,*) = A[0] * rad_fit_mk4_A^(-A[1]) + A[2] * rad_fit_mk4_A^(-A[3])                   ; cm-3
                dNe_dr              = - A[1]*A[0] * rad_fit_mk4_A^(-A[1]-1) - A[3]*A[2] * rad_fit_mk4_A^(-A[3]-1)   ; cm-3 / Rsun
                indsamp = where(rad_fit_mk4_A ge min(radsamp) and rad_fit_mk4_A le max(radsamp) AND dNe_dr lt 0.)
                if indsamp[0] ne -1 then begin
                   v = abs(dNe_dr(indsamp)/reform(Ne_fit_mk4_A(ifl,indsamp)))^(-1)
                   lN_fit_mk4_A(ifl)   = int_tabulated(rad_fit_mk4_A(indsamp), v) / (max(rad_fit_mk4_A(indsamp))-min(rad_fit_mk4_A(indsamp))) ; Rsun
                   ; print,lN_fit_mk4_A(ifl), float(mean(v)), float(median(v))
                   ; stop
                endif 
                skip_mk4_double_power_law:
             endif              ; covgflag = 'yes'
          endif
       endfor                   ; field lines loop.
       trace_data = create_struct( trace_data                          ,$
                                  'fitflag_mk4',ptr_new(fitflag_mk4_A) ,$
                                 'fit_F_Ne_mk4',ptr_new( fit_F_Ne_mk4) ,$
                                  'scN_fit_mk4',ptr_new(scN_fit_mk4_A) ,$
                                  'rad_fit_mk4',ptr_new(rad_fit_mk4_A) ,$
                                   'Ne_fit_mk4',ptr_new( Ne_fit_mk4_A) ,$
                                   'lN_fit_mk4',ptr_new(  lN_fit_mk4_A) )
       if fit_F_Ne_mk4 eq 'IHS' then $
          trace_data = create_struct( trace_data                        ,$
                                   'N0_fit_mk4',ptr_new(  N0_fit_mk4_A) )                                                                        
       if fit_F_Ne_mk4 eq 'SPL' then $
          trace_data = create_struct( trace_data                        ,$
                                    'N1_fit_mk4',ptr_new( N1_fit_mk4_A) ,$
                                    'p1_fit_mk4',ptr_new( p1_fit_mk4_A) )
       if fit_F_Ne_mk4 eq 'DPL' then $
          trace_data = create_struct( trace_data                        ,$
                                    'N1_fit_mk4',ptr_new( N1_fit_mk4_A) ,$
                                    'N2_fit_mk4',ptr_new( N2_fit_mk4_A) ,$
                                    'p1_fit_mk4',ptr_new( p1_fit_mk4_A) ,$
                                    'p2_fit_mk4',ptr_new( p2_fit_mk4_A) )
    endif ; Mk4
    
    if keyword_set(kcor) then begin
       read_tomgrid_and_define_fitgrid,fl_dir=fl_dir,instr_string='kcor'
       fitflag_kcor_A = fltarr(N_fl        ) + default
        N0_fit_kcor_A = fltarr(N_fl        ) + default
        lN_fit_kcor_A = fltarr(N_fl        ) + default
        N1_fit_kcor_A = fltarr(N_fl        ) + default
        N2_fit_kcor_A = fltarr(N_fl        ) + default
        p1_fit_kcor_A = fltarr(N_fl        ) + default
        p2_fit_kcor_A = fltarr(N_fl        ) + default
       scN_fit_kcor_A = fltarr(N_fl        ) + default
        Ne_fit_kcor_A = fltarr(N_fl,Npt_fit) + default
       rad_fit_kcor_A = radmin_fit + drad_fit/2. + drad_fit * findgen(Npt_fit)
       for ifl=0,N_fl-1 do begin      
          tmp = reform(index_sampling_kcor_A(ifl,*))
          ind_samp_kcor = where(tmp eq 1)
          if ind_samp_kcor[0] ne -1 then begin
             radsamp = reform(rad_A(ifl,ind_samp_kcor)) ; Rsun
             test_coverage, radsamp=radsamp, covgflag=covgflag, /kcor
             if covgflag eq 'yes' then begin
                fitflag_kcor_A(ifl) = +1.
                Nesamp = reform(Ne_kcor_A(ifl,ind_samp_kcor))
                fit_F_Ne_kcor  = 'DPL'
                double_power_fit, radsamp, Nesamp, A, chisqr ;, /weighted
                scN_fit_kcor_A(ifl)   = sqrt(chisqr)/mean(Nesamp)
                N1_fit_kcor_A(ifl)   = A[0]                                                                          ; cm-3
                p1_fit_kcor_A(ifl)   = A[1]                                                                          ; dimensionless exponent of power law
                N2_fit_kcor_A(ifl)   = A[2]                                                                          ; cm-3
                p2_fit_kcor_A(ifl)   = A[3]                                                                          ; dimensionless exponent of power law
                Ne_fit_kcor_A(ifl,*) = A[0] * rad_fit_kcor_A^(-A[1]) + A[2] * rad_fit_kcor_A^(-A[3])                   ; cm-3
                dNe_dr              = - A[1]*A[0] * rad_fit_kcor_A^(-A[1]-1) - A[3]*A[2] * rad_fit_kcor_A^(-A[3]-1)   ; cm-3 / Rsun
                indsamp = where(rad_fit_mk4_A ge min(radsamp) and rad_fit_mk4_A le max(radsamp) AND dNe_dr lt 0.)
                if indsamp[0] ne -1 then begin
                   v = abs(dNe_dr(indsamp)/reform(Ne_fit_kcor_A(ifl,indsamp)))^(-1)
                   lN_fit_kcor_A(ifl)   = int_tabulated(rad_fit_kcor_A(indsamp), v) / (max(rad_fit_kcor_A(indsamp))-min(rad_fit_kcor_A(indsamp))) ; Rsun
                   ; print,lN_fit_mk4_A(ifl), float(mean(v)), float(median(v))
                   ; stop
                endif 
             endif              ; covgflag = 'yes'
          endif
       endfor                   ; field lines loop.
       trace_data = create_struct( trace_data                          ,$
                                  'fitflag_kcor',ptr_new(fitflag_kcor_A) ,$
                                 'fit_F_Ne_kcor',ptr_new( fit_F_Ne_kcor) ,$
                                  'scN_fit_kcor',ptr_new(scN_fit_kcor_A) ,$
                                  'rad_fit_kcor',ptr_new(rad_fit_kcor_A) ,$
                                   'Ne_fit_kcor',ptr_new( Ne_fit_kcor_A) ,$
                                   'lN_fit_kcor',ptr_new(  lN_fit_kcor_A) )
       if fit_F_Ne_kcor eq 'DPL' then $
          trace_data = create_struct( trace_data                        ,$
                                    'N1_fit_kcor',ptr_new( N1_fit_kcor_A) ,$
                                    'N2_fit_kcor',ptr_new( N2_fit_kcor_A) ,$
                                    'p1_fit_kcor',ptr_new( p1_fit_kcor_A) ,$
                                    'p2_fit_kcor',ptr_new( p2_fit_kcor_A) )
    endif ; KCOR

    if keyword_set(ucomp) then begin
       read_tomgrid_and_define_fitgrid,fl_dir=fl_dir,instr_string='ucomp'
       fitflag_ucomp_A = fltarr(N_fl        ) + default
        N0_fit_ucomp_A = fltarr(N_fl        ) + default
        lN_fit_ucomp_A = fltarr(N_fl        ) + default
        N1_fit_ucomp_A = fltarr(N_fl        ) + default
        N2_fit_ucomp_A = fltarr(N_fl        ) + default
        p1_fit_ucomp_A = fltarr(N_fl        ) + default
        p2_fit_ucomp_A = fltarr(N_fl        ) + default
       scN_fit_ucomp_A = fltarr(N_fl        ) + default
        Ne_fit_ucomp_A = fltarr(N_fl,Npt_fit) + default
       rad_fit_ucomp_A = radmin_fit + drad_fit/2. + drad_fit * findgen(Npt_fit)
       for ifl=0,N_fl-1 do begin      
          tmp = reform(index_sampling_ucomp_A(ifl,*))
          ind_samp_ucomp = where(tmp eq 1)
          if ind_samp_ucomp[0] ne -1 then begin
             radsamp = reform(rad_A(ifl,ind_samp_ucomp)) ; Rsun
             test_coverage, radsamp=radsamp, covgflag=covgflag, /ucomp
             if covgflag eq 'yes' then begin
                fitflag_ucomp_A(ifl) = +1.
                Nesamp = reform(Ne_ucomp_A(ifl,ind_samp_ucomp))
                fit_F_Ne_ucomp  = 'DPL'
                double_power_fit, radsamp, Nesamp, A, chisqr ;, /weighted
                scN_fit_ucomp_A(ifl)   = sqrt(chisqr)/mean(Nesamp)
                N1_fit_ucomp_A(ifl)   = A[0]                                                                          ; cm-3
                p1_fit_ucomp_A(ifl)   = A[1]                                                                          ; dimensionless exponent of power law
                N2_fit_ucomp_A(ifl)   = A[2]                                                                          ; cm-3
                p2_fit_ucomp_A(ifl)   = A[3]                                                                          ; dimensionless exponent of power law
                Ne_fit_ucomp_A(ifl,*) = A[0] * rad_fit_kcor_A^(-A[1]) + A[2] * rad_fit_kcor_A^(-A[3])                   ; cm-3
                dNe_dr              = - A[1]*A[0] * rad_fit_ucomp_A^(-A[1]-1) - A[3]*A[2] * rad_fit_ucomp_A^(-A[3]-1)   ; cm-3 / Rsun
                indsamp = where(rad_fit_ucomp_A ge min(radsamp) and rad_fit_ucomp_A le max(radsamp) AND dNe_dr lt 0.)
                if indsamp[0] ne -1 then begin
                   v = abs(dNe_dr(indsamp)/reform(Ne_fit_ucomp_A(ifl,indsamp)))^(-1)
                   lN_fit_ucomp_A(ifl)   = int_tabulated(rad_fit_ucomp_A(indsamp), v) / (max(rad_fit_ucomp_A(indsamp))-min(rad_fit_ucomp_A(indsamp))) ; Rsun
                   ; print,lN_fit_mk4_A(ifl), float(mean(v)), float(median(v))
                   ; stop
                endif 
             endif              ; covgflag = 'yes'
          endif
       endfor                   ; field lines loop.
       trace_data = create_struct( trace_data                          ,$
                                  'fitflag_ucomp',ptr_new(fitflag_ucomp_A) ,$
                                 'fit_F_Ne_ucomp',ptr_new( fit_F_Ne_ucomp) ,$
                                  'scN_fit_ucomp',ptr_new(scN_fit_ucomp_A) ,$
                                  'rad_fit_ucomp',ptr_new(rad_fit_ucomp_A) ,$
                                   'Ne_fit_ucomp',ptr_new( Ne_fit_ucomp_A) ,$
                                   'lN_fit_ucomp',ptr_new( lN_fit_ucomp_A) )
       if fit_F_Ne_kcor eq 'DPL' then $
          trace_data = create_struct( trace_data                        ,$
                                    'N1_fit_ucomp',ptr_new( N1_fit_ucomp_A) ,$
                                    'N2_fit_ucomp',ptr_new( N2_fit_ucomp_A) ,$
                                    'p1_fit_ucomp',ptr_new( p1_fit_ucomp_A) ,$
                                    'p2_fit_ucomp',ptr_new( p2_fit_ucomp_A) )
    endif ; UCoMP

    if keyword_set(lascoc2) then begin
       read_tomgrid_and_define_fitgrid,fl_dir=fl_dir,instr_string='lascoc2'
       fitflag_c2_A = fltarr(N_fl        ) + default
        N1_fit_c2_A = fltarr(N_fl        ) + default
        N2_fit_c2_A = fltarr(N_fl        ) + default
        p1_fit_c2_A = fltarr(N_fl        ) + default
        p2_fit_c2_A = fltarr(N_fl        ) + default
        lN_fit_c2_A = fltarr(N_fl        ) + default
       scN_fit_c2_A = fltarr(N_fl        ) + default
        Ne_fit_c2_A = fltarr(N_fl,Npt_fit) + default
       rad_fit_c2_A = radmin_fit + drad_fit/2. + drad_fit * findgen(Npt_fit)
       for ifl=0,N_fl-1 do begin      
          tmp = reform(index_sampling_c2_A(ifl,*))
          ind_samp_c2 = where(tmp eq 1)
          if ind_samp_c2[0] ne -1 then begin
             radsamp = reform(rad_A(ifl,ind_samp_c2)) ; Rsun
             test_coverage, radsamp=radsamp, covgflag=covgflag, /lascoc2
             if covgflag eq 'yes' then begin
                fitflag_c2_A(ifl) = +1.
                Nesamp = reform(Ne_c2_A(ifl,ind_samp_c2))
                goto,skip_c2_single_power_law
                fit_F_Ne_c2 = 'SPL'
                linear_fit, alog(radsamp), alog(Nesamp), AN, r2N, /linfit_idl
                scN_fit_c2_A(ifl)  = r2N
                N1_fit_c2_A(ifl)   = exp(AN[0])                                           ; cm-3
                p1_fit_c2_A(ifl)   =    -AN[1]                                            ; dimensionless exponent of power law
                Ne_fit_c2_A(ifl,*) = N1_fit_c2_A(ifl) * rad_fit_c2_A^(-p1_fit_c2_A(ifl))  ; cm-3
                indsamp = where(rad_fit_c2_A ge min(radsamp) and rad_fit_c2_A le max(radsamp) AND dNe_dr lt 0.)
                v =  abs(rad_fit_c2_A(indsamp) / float(p1_fit_c2_A(ifl))) 
                lN_fit_c2_A(ifl)   = int_tabulated( rad_fit_c2_A(indsamp), v) / (max(rad_fit_c2_A(indsamp))-min(rad_fit_c2_A(indsamp))) ; Rsun
                print,lN_fit_c2_A(ifl), float(mean(v)), float(median(v))
                skip_c2_single_power_law:
               ;goto,skip_c2_double_power_law
                fit_F_Ne_c2  = 'DPL'
                double_power_fit, radsamp, Nesamp, A, chisqr ;, /weighted
                scN_fit_c2_A(ifl)   = sqrt(chisqr)/mean(Nesamp)
                N1_fit_c2_A(ifl)   = A[0]                                                                         ; cm-3
                p1_fit_c2_A(ifl)   = A[1]                                                                         ; dimensionless exponent of power law
                N2_fit_c2_A(ifl)   = A[2]                                                                         ; cm-3
                p2_fit_c2_A(ifl)   = A[3]                                                                         ; dimensionless exponent of power law
                Ne_fit_c2_A(ifl,*) = A[0] * rad_fit_c2_A^(-A[1]) + A[2] * rad_fit_c2_A^(-A[3])                    ; cm-3
                dNe_dr              = - A[1]*A[0] * rad_fit_c2_A^(-A[1]-1) - A[3]*A[2] * rad_fit_c2_A^(-A[3]-1)   ; cm-3 / Rsun
                indsamp = where(rad_fit_c2_A ge min(radsamp) and rad_fit_c2_A le max(radsamp) AND dNe_dr lt 0.)
                if indsamp[0] ne -1 then begin
                   v = abs(dNe_dr(indsamp)/reform(Ne_fit_c2_A(ifl,indsamp)))^(-1)
                   lN_fit_c2_A(ifl)   = int_tabulated(rad_fit_c2_A(indsamp), v) / (max(rad_fit_c2_A(indsamp))-min(rad_fit_c2_A(indsamp))) ; Rsun
                endif
                 ; print,lN_fit_c2_A(ifl), float(mean(v)), float(median(v))
                 ; stop
             skip_c2_double_power_law:
             endif              ; covgflag = 'yes'
          endif
       endfor                   ; field lines loop.
       trace_data = create_struct( trace_data                        ,$
                                  'fitflag_c2',ptr_new(fitflag_c2_A) ,$
                                 'fit_F_Ne_c2',ptr_new( fit_F_Ne_c2) ,$
                                  'scN_fit_c2',ptr_new(scN_fit_c2_A) ,$                                   
                                  'rad_fit_c2',ptr_new(rad_fit_c2_A) ,$
                                   'Ne_fit_c2',ptr_new( Ne_fit_c2_A) ,$
                                   'lN_fit_c2',ptr_new( lN_fit_c2_A) )
       if fit_F_Ne_c2 eq 'SPL' then $
          trace_data = create_struct( trace_data                     ,$
                                   'N1_fit_c2',ptr_new( N1_fit_c2_A) ,$
                                   'p1_fit_c2',ptr_new( p1_fit_c2_A) )
       if fit_F_Ne_c2 eq 'DPL' then $
          trace_data = create_struct( trace_data                     ,$
                                   'N1_fit_c2',ptr_new( N1_fit_c2_A) ,$
                                   'N2_fit_c2',ptr_new( N2_fit_c2_A) ,$
                                   'p1_fit_c2',ptr_new( p1_fit_c2_A) ,$
                                   'p2_fit_c2',ptr_new( p2_fit_c2_A) )
    endif ; LASCOC2

    return
 end
