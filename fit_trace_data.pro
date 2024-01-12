pro fit_trace_data, aia = aia, euvia = euvia, euvib = euvib, eit = eit,$
                    mk4 = mk4, kcor = kcor, lascoc2 = lascoc2

    common all_data, trace_data, $
     N_fl, Npt_max, Npt_v,$
     x_A, y_A, z_A, rad_A, lat_A, lon_A,$
     Ne_aia_A, Tm_aia_A, index_aia_A, index_sampling_aia_A,$
     Ne_euvia_A, Tm_euvia_A, index_euvia_A, index_sampling_euvia_A,$
     Ne_euvib_A, Tm_euvib_A, index_euvib_A, index_sampling_euvib_A,$
     Ne_eit_A, Tm_eit_A, index_eit_A, index_sampling_eit_A,$
     Ne_mk4_A, index_mk4_A, index_sampling_mk4_A,$
     Ne_kcor_A, index_kcor_A, index_sampling_kcor_A,$
     Ne_c2_A, index_c2_A, index_sampling_c2_A

    default = -678.
    
    if keyword_set(aia) then begin
       Npt_fit = 30 & rmin = 1.0 & rmax = 1.3
        dr_fit = (rmax-rmin)/float(Npt_fit)
       fitflag_aia_A = fltarr(N_fl        ) + default
        N0_fit_aia_A = fltarr(N_fl        ) + default
        lN_fit_aia_A = fltarr(N_fl        ) + default
       r2N_fit_aia_A = fltarr(N_fl        ) + default
        T0_fit_aia_A = fltarr(N_fl        ) + default
      dTdr_fit_aia_A = fltarr(N_fl        ) + default
       r2T_fit_aia_A = fltarr(N_fl        ) + default
        Ne_fit_aia_A = fltarr(N_fl,Npt_fit) + default
        Tm_fit_aia_A = fltarr(N_fl,Npt_fit) + default
         r_fit_aia_A = rmin + dr_fit/2. + dr_fit * findgen(Npt_fit)
       for ifl=0,N_fl-1 do begin      
          tmp = reform(index_sampling_aia_A(ifl,*))
          ind_samp_aia = where(tmp eq 1)
          rsamp = rad_A(ifl,ind_samp_aia) ; Rsun
          test_coverage, rsamp=rsamp, covgflag=covgflag, /aia
          if covgflag eq 'yes' then begin
             fitflag_aia_A(ifl) = +1.
             Nesamp = Ne_aia_A(ifl,ind_samp_aia)
             Tmsamp = Tm_aia_A(ifl,ind_samp_aia)
             linear_fit, 1./rsamp   , alog(Nesamp), AN, r2N, /linfit_idl
             linear_fit,    rsamp-1.,      Tmsamp , AT, r2T, /linfit_idl
             r2N_fit_aia_A(ifl)   = r2N
             r2T_fit_aia_A(ifl)   = r2T
              N0_fit_aia_A(ifl)   = exp(AN[0]+AN[1]) ; cm-3
              lN_fit_aia_A(ifl)   = 1./AN[1]         ; Rsun
              T0_fit_aia_A(ifl)   = AT[0]            ; K
            dTdr_fit_aia_A(ifl)   = AT[1]            ; K/Rsun
              Ne_fit_aia_A(ifl,*) = N0_fit_aia_A(ifl) * exp(-(1/lN_fit_aia_A(ifl))*(1.-1./r_fit_aia_A   )) ; cm-3
              Tm_fit_aia_A(ifl,*) = T0_fit_aia_A(ifl) + dTdr_fit_aia_A(ifl)       *      (r_fit_aia_A-1.)  ; K
          endif                 ; covgflag = 'yes'          
       endfor                   ; field lines loop.
       trace_data = create_struct( trace_data                          ,$
                                  'fitflag_aia',ptr_new(fitflag_aia_A) ,$
                                  'r2N_fit_aia',ptr_new(r2N_fit_aia_A) ,$
                                  'r2T_fit_aia',ptr_new(r2T_fit_aia_A) ,$
                                    'r_fit_aia',ptr_new(  r_fit_aia_A) ,$
                                   'Ne_fit_aia',ptr_new( Ne_fit_aia_A) ,$
                                   'Tm_fit_aia',ptr_new( Tm_fit_aia_A) )
    endif ; AIA


    if keyword_set(mk4) then begin
       Npt_fit = 40 & rmin = 1.1 & rmax = 1.5
        dr_fit = (rmax-rmin)/float(Npt_fit)
       fitflag_mk4_A = fltarr(N_fl        ) + default
        N0_fit_mk4_A = fltarr(N_fl        ) + default
        lN_fit_mk4_A = fltarr(N_fl        ) + default
       r2N_fit_mk4_A = fltarr(N_fl        ) + default
        Ne_fit_mk4_A = fltarr(N_fl,Npt_fit) + default
         r_fit_mk4_A = rmin + dr_fit/2. + dr_fit * findgen(Npt_fit)
       for ifl=0,N_fl-1 do begin      
          tmp = reform(index_sampling_mk4_A(ifl,*))
          ind_samp_mk4 = where(tmp eq 1)
          rsamp = rad_A(ifl,ind_samp_mk4) ; Rsun
          test_coverage, rsamp=rsamp, covgflag=covgflag, /mk4
          if covgflag eq 'yes' then begin
             fitflag_mk4_A(ifl) = +1.
             Nesamp = Ne_mk4_A(ifl,ind_samp_mk4)
             linear_fit, 1./rsamp   ,alog(Nesamp), AN, r2N, /linfit_idl
             r2N_fit_mk4_A(ifl)   = r2N
              N0_fit_mk4_A(ifl)   = exp(AN[0]+AN[1]) ; cm-3
              lN_fit_mk4_A(ifl)   = 1./AN[1]         ; Rsun
              Ne_fit_mk4_A(ifl,*) = N0_fit_mk4_A(ifl) * exp(-(1/lN_fit_mk4_A(ifl))*(1.-1./r_fit_mk4_A   )) ; cm-3
          endif                 ; covgflag = 'yes'          
       endfor                   ; field lines loop.
       trace_data = create_struct( trace_data                          ,$
                                  'fitflag_mk4',ptr_new(fitflag_mk4_A) ,$
                                  'r2N_fit_mk4',ptr_new(r2N_fit_mk4_A) ,$
                                    'r_fit_mk4',ptr_new(  r_fit_mk4_A) ,$
                                   'Ne_fit_mk4',ptr_new( Ne_fit_mk4_A) )
    endif ; MK4

    if keyword_set(lascoc2) then begin
       Npt_fit = 35 & rmin = 2.5 & rmax = 6.0
        dr_fit = (rmax-rmin)/float(Npt_fit)
         Rcrit = rmin
       fitflag_c2_A = fltarr(N_fl        ) + default
        N0_fit_c2_A = fltarr(N_fl        ) + default
         p_fit_c2_A = fltarr(N_fl        ) + default
       r2N_fit_c2_A = fltarr(N_fl        ) + default
        Ne_fit_c2_A = fltarr(N_fl,Npt_fit) + default
         r_fit_c2_A = rmin + dr_fit/2. + dr_fit * findgen(Npt_fit)
       for ifl=0,N_fl-1 do begin      
          tmp = reform(index_sampling_c2_A(ifl,*))
          ind_samp_c2 = where(tmp eq 1)
          rsamp = rad_A(ifl,ind_samp_c2) ; Rsun
          test_coverage, rsamp=rsamp, covgflag=covgflag, /lascoc2
          if covgflag eq 'yes' then begin
             fitflag_c2_A(ifl) = +1.
             Nesamp = Ne_c2_A(ifl,ind_samp_c2)
             linear_fit, alog(rsamp/Rcrit), alog(Nesamp), AN, r2N, /linfit_idl
             r2N_fit_c2_A(ifl)   = r2N
              N0_fit_c2_A(ifl)   = exp(AN[0]) ; cm-3
               p_fit_c2_A(ifl)   =    -AN[1]  ; dimensionless exponent of power law
              Ne_fit_c2_A(ifl,*) = N0_fit_c2_A(ifl) * (r_fit_c2_A / Rcrit)^(-p_fit_c2_A(ifl)) ; cm-3
          endif                 ; covgflag = 'yes'          
       endfor                   ; field lines loop.
       trace_data = create_struct( trace_data                        ,$
                                  'fitflag_c2',ptr_new(fitflag_c2_A) ,$
                                  'r2N_fit_c2',ptr_new(r2N_fit_c2_A) ,$                                   
                                    'r_fit_c2',ptr_new(  r_fit_c2_A) ,$
                                   'Ne_fit_c2',ptr_new( Ne_fit_c2_A) )
    endif ; LASCOC2

    return
 end

