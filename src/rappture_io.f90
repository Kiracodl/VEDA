!-----------------------------------------------------------------------
!             Copyright 2007 - 2012 by
!             John Melcher (jmelcher@purdue.edu), Daniel Kiracofe (drkiraco@purdue.edu)
!             Shuiqing Hu, Steven Johnson, Arvind Raman
!             MECHANICAL ENGINEERING, PURDUE UNIVERSITY
!             WEST LAFAYETTE, IN, 47907
!
!-----------------------------------------------------------------------
!
!VEDA is licensed under the Q Public License version 1.0 as reproduced in the LICENSE.txt file that
!accompanied this file, with the following addition:
!
!If VEDA (modified or unmodified), any portion thereof, or any derivative work (publicly distributed or not), contributes to any scientific or academic publication (including, but not limited to, journal articles, conference presentations or posters, or seminars) you must cite the original developers.
!
!Kiracofe, D.; Melcher, J. & Raman, A. "Gaining insight into the physics of dynamic atomic force microscopy in complex environments using the VEDA simulator Review of Scientific Instruments", (2012), 83, 013702
!
!You may also optionally cite the original VEDA article
!
! J. Melcher, S. Hu, A. Raman, "VEDA: A web-based virtual environment for dynamic atomic force microscopy" Review of Scientific Instruments 79, 061301 (2008).

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!file to put all rappture i/o in (with the intention that we could build a different executable 
!in the future that had all the same interals, but communicated with, say, matlab instead of rappture.

module rappture_io
  use params

  save

  real*8, allocatable :: outputbuffer(:,:,:,:)
  integer buffer_len4
  integer :: buffer_ndx  ! amp_index in main code can't be trusted.  that indexes the request time histories, but
!some might get skipped if there are jumps.  this indexes the time histories that actually occured.
  integer :: driver

  type time_history_settings
     logical :: want_voltage, want_TH_obs_defl, want_th_force_gap, want_surf_relax, want_contact_area, want_self_exc, want_FZ, want_acceleration, want_modal, want_PerfMetrics, want_SampleTopo
  end type time_history_settings

  type(time_history_settings), private :: th_set

contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!1
!first all of the reading routines, then all of the plot setup routines, then the writing routines

subroutine OpenInputFile(inFile)
  character(len=*), intent(in) :: inFile
  integer rp_lib

  driver = rp_lib(inFile)

  !fixme, really need a constructor. this is the best place I guess
  buffer_ndx = 0
end subroutine OpenInputFile

!unique inputs for the amplitude reduction tool
subroutine ReadAmpRedOpParameters( Asp_i, Asp_f, Asp_n,Asp_tolerance, iteration_method, ode_solver, z_debug_nm, A2_debug_nm, omega2_omega1_ratio, omega2_omega1_rounding, InputEcho)
  real*8, intent(out):: Asp_i, Asp_f, Asp_tolerance, z_debug_nm(100), A2_debug_nm(100), omega2_omega1_ratio
  integer, intent(out):: Asp_n, iteration_method, ode_solver, omega2_omega1_rounding
  character*1500, intent(inout) :: InputEcho
  character*200 strVal
  integer readGenericInteger
  integer status, i

  status = 0
  
  call read_value(driver, INPUT_PREFIX // "(op).number(Asp_i).current", status, Asp_i, InputEcho, "Asp_i "   )
  call read_value(driver, INPUT_PREFIX // "(op).number(Asp_f).current", status, Asp_f, InputEcho, "Asp_f "   )

  call read_value(driver, INPUT_PREFIX // "(sim).number(AmpRedIterTol).current", status, Asp_tolerance, InputEcho, "Asp_tol "   )

  call read_integer( driver,  INPUT_PREFIX // "(op).integer(Asp_n).current", Asp_n, InputEcho, "Asp_n ")

  call read_integer( driver,  INPUT_PREFIX // "(sim).group(backwards_compat).choice(ode_solver).current" , ode_solver , InputEcho, "ode_solver ")
  
  call read_integer(driver,INPUT_PREFIX//"(sim).group(backwards_compat).choice(iteration_method).current", iteration_method, InputEcho, "iteration_method ")

call read_value(driver, INPUT_PREFIX // "(op).number(omega2_omega1_ratio).current", status, omega2_omega1_ratio , InputEcho, "omega2_omega1_ratio "   )
  
  call read_integer(driver,INPUT_PREFIX//"(op).choice(omega1_omega2_rounding).current", omega2_omega1_rounding, InputEcho, "omega2_omega1_rounding ")

  
if ((iteration_method == 1) .or.  (iteration_method == 4)) then
     if (Asp_n > 100) call WriteFatalError("fixed z input only works for 99 Z points or less")
     
     strVal = rp_lib_get_wrap(driver, INPUT_PREFIX // "(sim).group(backwards_compat).string(debug_z_nm).current" )
     read( strVal, *, end=800, err=800) (Z_debug_nm(i),i=1,Asp_n)
  end if

  if (iteration_method ==1) then
     strVal = rp_lib_get_wrap(driver, INPUT_PREFIX // "(sim).group(backwards_compat).string(debug_A2_nm).current" )
     read( strVal, *, end=801, err=801) (A2_debug_nm(i),i=1,Asp_n)
  end if
  
  if (status > 0) call WriteFatalError( "Error reading Amplitude Reduction tool parameters. Please check values") 

return
  
  800  call WriteFatalError( "Could not parse debug Z values. Did you specify a comma separated list with the right number of values?")     
  801  call WriteFatalError( "Could not parse debug A2 values. Did you specify a comma separated list with the right number of values?")    
end subroutine ReadAmpRedOpParameters

subroutine ReadHydrodynamicFunction(fluid_rho, fluid_eta, cant_width, cant_len, InputEcho)
  real*8, intent(out) :: fluid_rho, fluid_eta, cant_width, cant_len
character*1500, intent(inout) :: InputEcho
  integer status

  status = 0
  
  call read_value(driver, INPUT_PREFIX // "(op).number(fluid_rho).current", status, fluid_rho, InputEcho, "fluid_rho "   )
  call read_value(driver, INPUT_PREFIX // "(op).number(fluid_eta).current", status, fluid_eta, InputEcho, "fluid_eta "   )
  call read_value(driver, INPUT_PREFIX // "(op).number(cant_width).current", status, cant_width, InputEcho, "cant_width "   )
  cant_width = cant_width * 1d-6;
  call read_value(driver, INPUT_PREFIX // "(op).number(cant_len).current", status, cant_len, InputEcho, "cant_len "   )
  cant_len = cant_len * 1d-6;

  if (status > 0) call WriteFatalError( "Error reading hydrodynamic function parameters. Please check values") 
end subroutine ReadHydrodynamicFunction

subroutine ReadControllerProperties( LineSpeed,  NoiseAmp, KP, KI, WantConstZ, Z_feedback_choice, modulation_type, transient_timeout_in, InputEcho)
  real*8, intent(out) :: LineSpeed,  NoiseAmp, KP, KI, transient_timeout_in
  integer, intent(out), optional :: Z_feedback_choice
  character*1500, intent(inout) :: InputEcho
  logical, intent(out) :: WantConstZ
  integer, intent(in) :: modulation_type

  real*8 :: SNratio
  character*100 strVal, tmpStr
  integer status, rp_units_convert_dbl, readGenericInteger

  status = 0

  call read_value(driver, INPUT_PREFIX // "(sim).number(Transient_timeout).current", status, transient_timeout_in, InputEcho, "transient_timeout " )
  if (status > 0) then
     !did not implement in all tools yet.  quick hack to default back to previous behavior
     status = 0
     transient_timeout_in = 0
  end if
  
  call read_value(driver, INPUT_PREFIX // "(op).number(LineSpeed).current", status, LineSpeed, InputEcho, "LineSpeed "   )
  
  z_feedback_choice = readGenericInteger(driver, INPUT_PREFIX // "(op).choice(Z_feedback_choice).current")
  write( tmpstr, *)   Z_feedback_choice
  InputEcho = trim(InputEcho) // "Z_feedback_choice " // trim( tmpstr) // char(10)

  call read_value(driver, INPUT_PREFIX // "(op).number(SNratio).current", status, SNratio, InputEcho, "SNratio "   )
  NoiseAmp = 10**(-SNratio/20)
  ! NoiseAmp is relative to the unconstrained deflection amplitude   

  
  strVal = rp_lib_get_wrap(driver,  INPUT_PREFIX // "(op).number(KP).current")
  InputEcho = trim(InputEcho) // "KP " // trim(strVal) // char(10)
  status = status +  rp_units_convert_dbl(strVal," ",KP)
  
  strVal = rp_lib_get_wrap(driver, INPUT_PREFIX // "(op).number(KI).current")
  InputEcho = trim(InputEcho) // "KI " // trim(strVal) // char(10)
  status = status +  rp_units_convert_dbl(strVal," ",KI)
  
  WantConstZ = daniel_get_boolean( INPUT_PREFIX // "(op).boolean(ConstZ).current")
  write( tmpstr, *)  WantConstZ
  InputEcho = trim(InputEcho) // "ConstZ " // trim(tmpstr) // char(10)
  
  if (status > 0) call WriteFatalError( "Error reading controller properties. Please check values") 
end subroutine ReadControllerProperties

subroutine ReadPeakForce(Zbase_Amp_dim, peakF_dim, want_FZ, InputEcho)
  real*8, intent(out) :: Zbase_Amp_dim, peakF_dim
  integer status, rp_units_convert_dbl
  logical, intent(out) :: want_FZ
  character*100 strVal
  character*1500, intent(inout) :: InputEcho
  status = 0

  strVal = rp_lib_get_wrap(driver,  INPUT_PREFIX // "(op).number(Zbase_amp).current")
  InputEcho = trim(InputEcho) // "Zbase_amp " // StrVal // char(10)
  status = status +  rp_units_convert_dbl(strVal," ",Zbase_amp_dim)
  Zbase_amp_dim = Zbase_amp_dim/1d9 
  
  strVal = rp_lib_get_wrap(driver,  INPUT_PREFIX // "(op).number(PeakF_dim).current")
  InputEcho = trim(InputEcho) // "PeakF_dim " // StrVal // char(10)
  status = status +  rp_units_convert_dbl(strVal," ",PeakF_dim)

  Want_FZ = daniel_get_boolean( INPUT_PREFIX // "(sim).group(timehist).boolean(Want_FZ).current")
 
  if (status > 0) call WriteFatalError( "Error reading Peak Force Values. Please check values") 
  

end subroutine ReadPeakForce

subroutine ReadConstrainedAmp(Amp_constr_dim, InputEcho)
  real*8, intent(out) :: Amp_constr_dim
  integer status, rp_units_convert_dbl
  character*1500, intent(inout) :: InputEcho
  character*100 strVal 
  status = 0
  strVal = rp_lib_get_wrap(driver,  INPUT_PREFIX // "(op).number(Amp_constr).current")
  InputEcho = trim(InputEcho) // "Amp_constr " // StrVal // char(10)
  status = status +  rp_units_convert_dbl(strVal," ",Amp_constr_dim)
 if (status > 0) call WriteFatalError( "Error reading constrained amplitude. Please check values") 

end subroutine ReadConstrainedAmp



subroutine ReadFMData(  fm_gain_k0, fm_gain_i0, fm_gain_d0, fm_gain_k1, fm_gain_i1,  fm_gain_d1, freq_shift_sp, fm_direct_control, z_feedback_choice, fm_want_noncontact,  WantCalcPLLGains, WantCalcAmpGains, want_pre_BPF, Want_NormFreqShift, InputEcho )
  real*8, intent(out) :: fm_gain_k0, fm_gain_i0, fm_gain_k1 , fm_gain_i1, freq_shift_sp,  fm_gain_d0,  fm_gain_d1
  integer, intent(in) :: z_feedback_choice
  logical, intent(out) ::  fm_direct_control, fm_want_noncontact,  WantCalcPLLGains,WantCalcAmpGains, want_pre_BPF,     Want_NormFreqShift
  character*1500, intent(inout) :: InputEcho

  integer :: rp_units_convert_dbl, readGenericInteger, status, tmp
  character*105:: strVal, tmpStr

  status = 0

  want_pre_BPF = daniel_get_boolean(INPUT_PREFIX // "(fm).boolean(want_pre_BPF).current")
  write(tmpStr, *) want_pre_BPF
  InputEcho = trim(InputEcho) // "want_pre_BPF " // tmpStr // char(10)

  fm_direct_control = daniel_get_boolean( INPUT_PREFIX // "(fm).boolean(fm_direct_control).current")
  write(tmpStr, *) fm_direct_control
  InputEcho = trim(InputEcho) // "fm_direct_control " // tmpStr // char(10)


  wantCalcPLLGains = daniel_get_boolean(INPUT_PREFIX // "(fm).boolean(WantCalcFMGains).current")
  write(tmpStr, *) wantCalcPLLGains
  InputEcho = trim(InputEcho) // "wantCalcPLLGains " // tmpStr // char(10)

  if (.not. wantCalcPLLGains) then
     strVal = rp_lib_get_wrap(driver, INPUT_PREFIX // "(fm).number(FmGainI0).current" )
     status = status +  rp_units_convert_dbl(strVal," ", fm_gain_i0  )
     InputEcho = trim(InputEcho) // "FmGainI0 = " // trim(strVal) // char(10)
     
     strVal = rp_lib_get_wrap(driver, INPUT_PREFIX // "(fm).number(FmGainK0).current" )
     status = status +  rp_units_convert_dbl(strVal," ", fm_gain_k0  )
     InputEcho = trim(InputEcho) // "FmGainK0 = " // trim(strVal) // char(10)

     strVal = rp_lib_get_wrap(driver, INPUT_PREFIX // "(fm).number(FmGainD0).current" )
     status = status +  rp_units_convert_dbl(strVal," ", fm_gain_d0  )
     InputEcho = trim(InputEcho) // "FmGainD0 = " // trim(strVal) // char(10)
  end if

  wantCalcAmpGains = daniel_get_boolean(INPUT_PREFIX // "(fm).boolean(WantCalcAmpGains).current")
  write(tmpStr, *) wantCalcAmpGains
  InputEcho = trim(InputEcho) // "wantCalcAmpGains " // tmpStr // char(10)

  if (.not. wantCalcAmpGains) then  
     strVal = rp_lib_get_wrap(driver, INPUT_PREFIX // "(fm).number(FmGainK1).current" )
     status = status +  rp_units_convert_dbl(strVal," ", fm_gain_k1  )
     InputEcho = trim(InputEcho) // "FmGainK1 = " // trim(strVal) // char(10)
     
     strVal = rp_lib_get_wrap(driver, INPUT_PREFIX // "(fm).number(FmGainI1).current" )
     status = status +  rp_units_convert_dbl(strVal," ", fm_gain_i1  )
     InputEcho = trim(InputEcho) // "FmGainI1 = " // trim(strVal) // char(10)

     strVal = rp_lib_get_wrap(driver, INPUT_PREFIX // "(fm).number(FmGainD1).current" )
     status = status +  rp_units_convert_dbl(strVal," ", fm_gain_d1  )
     InputEcho = trim(InputEcho) // "FmGainD1 = " // trim(strVal) // char(10)
  end if
     
  if ( Z_feedback_choice == PHASE_Z) then
     strVal = rp_lib_get_wrap(driver, INPUT_PREFIX // "(fm).number(phase_sp).current" )
     status = status +  rp_units_convert_dbl(strVal," ", freq_shift_sp  )
     InputEcho = trim(InputEcho) // "phase_sp = " // trim(strVal) // char(10)
  else
     strVal = rp_lib_get_wrap(driver, INPUT_PREFIX // "(fm).number(freq_shift_sp).current" )
     status = status +  rp_units_convert_dbl(strVal," ", freq_shift_sp  )
     InputEcho = trim(InputEcho) // "freq_shift_sp = " // trim(strVal) // char(10)
     freq_shift_sp = freq_shift_sp * 2 * pi
  end if

  tmp = readGenericInteger(driver, INPUT_PREFIX // "(fm).choice(FM_want_nc).current")
  write( tmpstr, *)  tmp
  InputEcho = trim(InputEcho) // "FM_want_nc " // trim( tmpstr) // char(10)
  fm_want_noncontact = (tmp == 1)

    Want_NormFreqShift = daniel_get_boolean( INPUT_PREFIX // "(fm).boolean(Want_NormFreqShift).current")

end subroutine ReadFMData

subroutine ReadTipData( Rtip_dim, Etip, Poisson_tip, tip_angle, tip_shape, InputEcho)
  real*8, intent(out) :: Rtip_dim, Etip, Poisson_tip, tip_angle
  integer, intent(out) :: tip_shape
  character*1500, intent(inout) :: InputEcho
  integer :: rp_units_convert_dbl, readGenericInteger
  integer :: status
  character*100:: strVal

  status = 0

  tip_shape =  readGenericInteger(driver,  INPUT_PREFIX // "(ts).choice(TipShape).current")
  write(strVal,*) tip_shape
  InputEcho = trim(InputEcho) // "tip shape = " // trim(strVal) // char(10)

  call read_value(driver, INPUT_PREFIX // "(ts).number(angle_tip).current" , status, tip_angle, InputEcho, "tip angle ")
  if (status == 1) then
     !backwards compat.
     tip_angle = 90
     tip_shape = PARABOLOID
     status = 0
  end if
  tip_angle = tip_angle * pi / 180d0

  call read_value(driver, INPUT_PREFIX // "(ts).number(Rtip).current" , status, Rtip_dim, InputEcho, "Rtip "   )  
  Rtip_dim = Rtip_dim/1d9
  !	Rtip_dim is the radius of the probe tip converted to (m)
  
  strVal = rp_lib_get_wrap(driver, INPUT_PREFIX // "(ts).number(Etip).current"  )
  InputEcho = trim(InputEcho) // "Etip = " // trim(strVal) // char(10)
  status = status +  rp_units_convert_dbl(strVal," ",Etip)
  Etip = Etip*1d9
  
  strVal = rp_lib_get_wrap(driver, INPUT_PREFIX // "(ts).number(Poisson_tip).current"  )
  InputEcho = trim(InputEcho) // "Poisson_tip = " // trim(strVal) // char(10)
  status = status +  rp_units_convert_dbl(strVal," ",Poisson_tip)
  
  if (status > 0) call WriteFatalError( "Could not read tip data.  Please check values.")

end subroutine ReadTipData

subroutine ReadSampleData( MatProp, path, InputEcho)
  
  use contactModels !this is ugly... really should move this routine inside contactModels. i.e. make it know how to read itself

  type(matl_prop), intent(out) :: MatProp
  character*3000, intent(inout) :: InputEcho

  character(len=*) :: path ! = "group(ts)" for substrate/normal and "group(feature).group(fp)" for feature
  character*200:: strVal
  character*500:: tmpstr
  integer :: status, rp_units_convert_dbl, readGenericInteger, i, readGenericInteger

  status = 0
  
  InputEcho = trim(InputEcho) // "mat properties for: " // path // char(10)
  
  MatProp%wantWLC = daniel_get_boolean(  path // ".boolean(WantWLC).current")
  write( tmpstr, *) MatProp%wantWLC
  InputEcho = trim(InputEcho) // "want WLC " // trim( tmpstr) // char(10)
  
  if (MatProp%wantWLC) then
     call read_value(driver, path // ".group(wlc).number(wlc_L0).current", status, MatProp%wlc_L0, InputEcho, "wlc_L0 "   )
     MatProp%wlc_L0 = MatProp%wlc_L0 * 1e-9
     call read_value(driver, path // ".group(wlc).number(wlc_Lp).current", status, MatProp%wlc_Lp, InputEcho, "wlc_Lp "   )
     MatProp%wlc_Lp = MatProp%wlc_Lp * 1e-9
     call read_value(driver, path // ".group(wlc).number(wlc_Fr).current", status, MatProp%wlc_Fr, InputEcho, "wlc_Fr ")
     MatProp%wlc_Fr = MatProp%wlc_Fr * 1e-9
  end if

  call read_value(driver, path // ".number(kts_R).current", status, MatProp%kts_R, InputEcho, "kts_R "   )
  call read_value(driver, path // ".number(kts_A).current", status, MatProp%kts_A, InputEcho, "kts_A "   )
		
  MatProp%CalcADMT =  readGenericInteger(driver, path // ".choice(CalcADMT).current")

  call read_value(driver, path // ".number(Fadhesion).current", status, MatProp%Fadhesion, InputEcho, "Fadhesion "   )
  MatProp%Fadhesion = MatProp%Fadhesion/1d9
           
  call read_value(driver, path // ".number(aDMT).current", status, MatProp%aDMT, InputEcho, "aDMT "   )
  MatProp%aDMT = MatProp%aDMT/1d9

!	A_hamaker is the hamaker consant user in the DMT model (J)
  call read_value(driver, path // ".number(A_hamaker).current", status, MatProp%A_hamaker_dim, InputEcho, "A_hamaker " )	

  MatProp%fts_model = readGenericInteger(driver, path // ".choice(fts_model).current") 	
  write( tmpstr, *) MatProp%fts_model
  InputEcho = trim(InputEcho) // "fts_model " // trim( tmpstr) // char(10)
  
  MatProp%Want_Tip_Squeeze =  daniel_get_boolean(path // ".group(nc).boolean(WantTipSqueeze).current")
  write( tmpstr, *) MatProp%Want_Tip_Squeeze
  InputEcho = trim(InputEcho) // "want tip squeeze " // trim( tmpstr) // char(10)
  
  if (MatProp%Want_Tip_Squeeze) then
     strVal = rp_lib_get_wrap(driver, path // ".group(nc).number(etaliquid).current"  )
     InputEcho = trim(InputEcho) // "etaliquid " // trim(strVal) // char(10)
     status = status +  rp_units_convert_dbl(strVal," ", MatProp%eta_liquid)
  else
     MatProp%eta_liquid=0
  end if
  
  MatProp%WantOscillatory =  daniel_get_boolean(  path // ".group(solvation).boolean(WantOscillatory).current")
  write( tmpstr, *) MatProp%WantOscillatory
  InputEcho = trim(InputEcho) // "want oscillatory " // trim( tmpstr) // char(10)


        MatProp%WantHydration =  daniel_get_boolean(  path // ".group(solvation).boolean(WantHydration).current")
        write( tmpstr, *) MatProp%WantHydration
        InputEcho = trim(InputEcho) // "want hydration " // trim( tmpstr) // char(10)


        if ( MatProp%WantOscillatory ) then
           strVal = rp_lib_get_wrap(driver, path // ".group(solvation).number(sigma_solvation).current"  )
           InputEcho = trim(InputEcho) // "sigma_solvation " // trim(strVal) // char(10)
           status = status +  rp_units_convert_dbl(strVal," ", MatProp%sigma_solvation_dim)
           MatProp%sigma_solvation_dim = MatProp%sigma_solvation_dim * 1d-9

           strVal = rp_lib_get_wrap(driver, path // ".group(solvation).number(rho_solvation).current"  )
           InputEcho = trim(InputEcho) // "rho_solvation " // trim(strVal) // char(10)
           status = status +  rp_units_convert_dbl(strVal," ", MatProp%rho)
        end if

        if ( MatProp%WantHydration ) then
           strVal = rp_lib_get_wrap(driver, path // ".group(solvation).number(lambda_solvation).current"  )
           InputEcho = trim(InputEcho) // "lambda_solvation " // trim(strVal) // char(10)
           status = status +  rp_units_convert_dbl(strVal," ", MatProp%lambda_solvation_dim)
           MatProp%lambda_solvation_dim = MatProp%lambda_solvation_dim * 1d-9

           strVal = rp_lib_get_wrap(driver, path // ".group(solvation).number(p_h).current"  )
           InputEcho = trim(InputEcho) // "p_h " // trim(strVal) // char(10)
           status = status +  rp_units_convert_dbl(strVal," ", MatProp%p_h)
        end if

        MatProp%Want_hyst_hydr =  daniel_get_boolean(  path // ".group(solvation).boolean(Want_hysteretic_hydration).current")
        write( tmpstr, *) MatProp%Want_hyst_hydr
        InputEcho = trim(InputEcho) // "Want_hsyteretic_hydration " // trim( tmpstr) // char(10)
       
        
        if ( MatProp%WantHydration .or.  MatProp%WantOscillatory .or. MatProp%Want_hyst_hydr ) then
           MatProp%WantContGrad = daniel_get_boolean(  path // ".group(solvation).boolean(cont_grad).current")
           write( tmpstr, *) MatProp%WantContGrad
           InputEcho = trim(InputEcho) // "cont grad " // trim( tmpstr) // char(10)
           
           MatProp%N_solv = 1  !fixme, this parameter should just be removed

           call read_value(driver, path // ".group(solvation).number(rtip_solv).current", status, MatProp%Rtip_solvation_dim, InputEcho,  "rtip_solv ")
           MatProp%Rtip_solvation_dim = MatProp%Rtip_solvation_dim * 1d-9

        end if

        if (MatProp%Want_hyst_hydr) then
           call read_value(driver, path // ".group(solvation).number(App_decay).current", status, MatProp%app_decay, InputEcho,  "app decay ")
           call read_value(driver, path // ".group(solvation).number(Ret_decay).current", status, MatProp%ret_decay, InputEcho,  "ret_decay ")
           call read_value(driver, path // ".group(solvation).number(app_scaling).current", status, MatProp%app_scaling, InputEcho,  "app_scaling ")
           call read_value(driver, path // ".group(solvation).number(cutoff_dist).current", status, MatProp%cutoff_dist, InputEcho,  "cutoff_dist ")
           MatProp%nside = readGenericInteger(driver, path // ".group(solvation).integer(nside).current")
           write( tmpstr, *)  MatProp%nside
           InputEcho = trim(InputEcho) // "nside " // trim( tmpstr) // char(10)
        end if


        if ( MatProp%fts_model == MAG_DIPOLE ) then
           call read_value(driver, path // ".group(md).number(mom_sample).current", status, MatProp%mom_sample, InputEcho,  "mom_sample ")
           call read_value(driver, path // ".group(md).number(mom_tip).current", status, MatProp%mom_tip, InputEcho,  "mom_tip ")
           call read_value(driver, path // ".group(md).number(mag_dia).current", status, MatProp%mag_dia, InputEcho,  "mag_dia ")
           MatProp%mag_dia = MatProp%mag_dia * 1e-9
           call read_value(driver, path // ".group(md).number(mag_delta).current", status, MatProp%mag_delta, InputEcho, "mag_delta ")
           MatProp%mag_delta = MatProp%mag_delta * 1e-9
        end if

        if (( MatProp%fts_model == ELECTROSTATIC_XU ) .or. ( MatProp%fts_model == ELECTROSTATIC_GIL) .or. ( MatProp%fts_model == ELECTROSTATIC_NONCONS)) then
           call read_value(driver, path // ".group(es).number(dc_bias_voltage).current", status, MatProp%dc_bias_voltage, InputEcho,  "dc_bias_voltage ")
           call read_value(driver, path // ".group(es).number(surface_pot).current", status, MatProp%surface_pot, InputEcho,  "surface_pot ")

           strVal = rp_lib_get_wrap(driver, path // ".group(es).number(ac_bias_voltage).current"  )
           status = status +  rp_units_convert_dbl(strVal," ",MatProp%ac_bias_voltage)
           InputEcho = trim(InputEcho) // "ac_bias_voltage " // trim(strVal) // char(10)

           strVal = rp_lib_get_wrap(driver, path // ".group(es).number(bias_freq).current"  )
           status = status +  rp_units_convert_dbl(strVal," ",MatProp%bias_freq_rads)
           InputEcho = trim(InputEcho) // "bias_freq " // trim(strVal) // char(10)
           MatProp%bias_freq_rads = MatProp%bias_freq_rads * 1e3 * 2 * pi
        end if

        if (( MatProp%fts_model == ELECTROSTATIC_GIL) .or. (MatProp%fts_model == ELECTROSTATIC_NONCONS)) then
           strVal = rp_lib_get_wrap(driver, path // ".group(es).number(theta_tip).current"  )
           InputEcho = trim(InputEcho) // "electrostatic theta tip " // trim(strVal) // char(10)
           status = status +  rp_units_convert_dbl(strVal," ", theta_tip)
           theta_tip = theta_tip * pi / 180d0

           strVal = rp_lib_get_wrap(driver, path // ".group(es).number(theta_lever).current"  )
           InputEcho = trim(InputEcho) // "electrostatic theta lever " // trim(strVal) // char(10)
           status = status +  rp_units_convert_dbl(strVal," ", theta_lever)
           theta_lever = theta_lever * pi / 180d0

           strVal = rp_lib_get_wrap(driver, path // ".group(es).number(es_h).current"  )
           InputEcho = trim(InputEcho) // "electrostatic height " // trim(strVal) // char(10)
           status = status +  rp_units_convert_dbl(strVal," ", es_h_dim)
           es_h_dim = es_h_dim * 1d-6

           strVal = rp_lib_get_wrap(driver, path // ".group(es).number(es_l).current"  )
           InputEcho = trim(InputEcho) // "electrostatic length " // trim(strVal) // char(10)
           status = status +  rp_units_convert_dbl(strVal," ", es_l_dim)
           es_l_dim = es_l_dim * 1d-6

           strVal = rp_lib_get_wrap(driver, path // ".group(es).number(es_w).current"  )
           InputEcho = trim(InputEcho) // "electrostatic width " // trim(strVal) // char(10)
           status = status +  rp_units_convert_dbl(strVal," ", es_w_dim)
           es_w_dim = es_w_dim * 1d-6
        end if

        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        !viscoelasticity starts here

        if   ((MatProp%fts_model == HERTZ) .or. ( MatProp%fts_model == DMT) ) then
           MatProp%VEchoice = readGenericInteger(driver, path // ".group(nc).choice(VEchoiceHertz).current")
        elseif (MatProp%fts_model == LINEAR) then
           MatProp%VEchoice = readGenericInteger(driver, path // ".group(nc).choice(VEchoiceLinear).current")
        elseif ( ( MatProp%fts_model == ATTARD_FOURIER_LSQ)  .or. (MatProp%fts_model == ATTARD_FOURIER_BAHRAM )) then
           MatProp%VEchoice = readGenericInteger(driver, path // ".group(nc).group(VEAttard).choice(VEChoiceAttard).current")
!           write(*,*) "foo: ", MatProp%VEchoice
        else
           MatProp%VEchoice = VEC_NONE
        end if
!backwards compatibility
        if (MatProp%VEchoice == 0)  MatProp%VEchoice = VEC_NONE
        write( tmpstr, *)  MatProp%VEchoice
        InputEcho = trim(InputEcho) // "VEchoice " // trim( tmpstr) // char(10)

        
        if (( MatProp%fts_model == ATTARD_FOURIER_BAHRAM) .or. ( MatProp%fts_model == ATTARD_FOURIER_LSQ)) then
             
           MatProp%N_attard_spatial = readGenericInteger( driver, path// ".group(nc).integer(N_attard_spatial).current")
           write( tmpstr, *)  MatProp%N_attard_spatial
           InputEcho = trim(InputEcho) // "N_attard " // trim( tmpstr) // char(10)

           ! call read_value(driver, path // ".group(nc).number(attard_stick_width).current", status, MatProp%attard_stick_width, InputEcho, "attard_stick_width "   )
           ! MatProp%attard_stick_width = MatProp%attard_stick_width * 1e-9

           ! call read_value(driver, path // ".group(nc).number(attard_stick_depth).current", status, MatProp%attard_stick_depth, InputEcho, "attard_stick_depth "   )
           ! MatProp%attard_stick_depth = MatProp%attard_stick_depth * 1e-12
           
           call read_value(driver, path // ".group(nc).number(attard_radial_extent).current", status, MatProp%attard_radial_extent, InputEcho, "attard_radial_extent "   )
           
           call read_value(driver, path // ".group(nc).number(attard_start).current", status, MatProp%attard_start_dim, InputEcho, "attard_start "   )
           MatProp%attard_start_dim = MatProp%attard_start_dim * 1e-9
           
           call read_value(driver, path // ".group(nc).number(attard_stop).current", status, MatProp%attard_stop_dim, InputEcho, "attard_stop "   )
           MatProp%attard_stop_dim = MatProp%attard_stop_dim * 1e-9
           
           MatProp%N_attard_fourier = readGenericInteger( driver, path// ".group(nc).integer(N_attard_fourier).current")
           write( tmpstr, *)  MatProp%N_attard_fourier
           InputEcho = trim(InputEcho) // "N_attard_fourier " // trim( tmpstr) // char(10)
           !call assert(MatProp%N_attard_fourier > 0, "MatProp%N_attard_fourier <= 0 ?")
        else
           MatProp%N_attard_fourier =0
        end if

        !this is more complicated than it needs to be, but this is the only way to make the gui look reasonable
        !and still have exact control over what options the user can pick.
        if  ((((MatProp%fts_model == DMT) .or. ( MatProp%fts_model == HERTZ)) .and. (MatProp%VEchoice /= VEC_NONE)) )  then
           MatProp%input_shear_modulus = (readGenericInteger( driver, path// ".group(nc).group(VEHertz).choice(modulus_type).current") == 2)

           call read_value(driver, path // ".group(nc).group(VEHertz).number(etasampleHertz).current", status, MatProp%etasample, InputEcho, "etasample "   )

           call read_value(driver, path // ".group(nc).group(VEHertz).number(threeelm_e2).current", status, MatProp%threeelm_e2, InputEcho, "threeelm_e2 "   )
           MatProp%threeelm_e2 = MatProp%threeelm_e2 * 1e9
           
           if ( MatProp%VEChoice == VEC_GENMAXWELL) then
              MatProp%N_gen_max = readGenericInteger(driver, path // ".group(nc).group(VEHertz).integer(gen_max_N).current")
              
              strVal = rp_lib_get_wrap(driver, path // ".group(nc).group(VEHertz).string(gen_max_E).current" )
              read( strVal, *, end=100, err=100) (MatProp%E_gen_max(i),i=1,MatProp%N_gen_max)
              inputEcho = trim(InputEcho) // "E_gen_max" // trim(strVal) // char(10)
              MatProp%E_gen_max = MatProp%E_gen_max * 1e9
           
              strVal = rp_lib_get_wrap(driver, path // ".group(nc).group(VEHertz).string(gen_max_tau).current" )
              read( strVal, *, end=101, err=101) (MatProp%tau_gen_max_dim(i),i=1,MatProp%N_gen_max)
              inputEcho = trim(InputEcho) // "tau_gen_max" // trim(strVal) // char(10)
           end if
        elseif  (( MatProp%fts_model == LIN_ATT)  .or. (MatProp%fts_model == LINEAR) ) then
           call read_value(driver, path // ".group(nc).number(etasampleLin).current", status, MatProp%etasample, InputEcho, "etasample ")
           call read_value(driver, path // ".group(nc).number(k2Lin).current", status, MatProp%k2, InputEcho, "k2 "   )
        elseif ((MatProp%fts_model == ATTARD_FOURIER_LSQ ) .or. ( MatProp%fts_model == ATTARD_FOURIER_BAHRAM)) then
           MatProp%input_shear_modulus = (readGenericInteger( driver, path// ".group(nc).group(VEAttard).choice(modulus_type).current") == 2)

!           write(*,*) MatProp%VEchoice
           
           if (MatProp%VEchoice ==  VEC_THREEELM_E1E2ETA) then           
              call read_value(driver, path // ".group(nc).group(VEAttard).number(etasampleHertz).current", status, MatProp%etasample, InputEcho, "etasample "   )
              call read_value(driver, path // ".group(nc).group(VEAttard).number(threeelm_e2).current", status, MatProp%threeelm_e2, InputEcho, "threeelm_e2 "   )
              MatProp%threeelm_e2 = MatProp%threeelm_e2 * 1e9
              call read_value(driver, path // ".group(nc).group(VEAttard).number(threeelm_e1).current", status, MatProp%threeelm_e1, InputEcho, "threeelm_e1 "   )
              MatProp%threeelm_e1 = MatProp%threeelm_e1 * 1e9
           elseif (MatProp%VEchoice ==  VEC_THREEELM_E0_EINF_TAU) then
              call read_value(driver, path // ".group(nc).group(VEAttard).number(threeelm_e0).current", status, MatProp%threeelm_e0, InputEcho, "threeelm_e0 "   )
              MatProp%threeelm_e0 = MatProp%threeelm_e0 * 1e9

!              write(*,*) "MatProp%threeelm_e0 ", MatProp%threeelm_e0
              
              call read_value(driver, path // ".group(nc).group(VEAttard).number(threeelm_einf).current", status, MatProp%threeelm_einf, InputEcho, "threeelm_einf "   )
              MatProp%threeelm_einf = MatProp%threeelm_einf * 1e9
              
              call read_value(driver, path // ".group(nc).group(VEAttard).number(threeelm_tau_stress).current", status, MatProp%threeelm_tau_stress_dim, InputEcho, "threeelm_tau_stress "   )

!              write(*,*) "MatProp%threeelm_tau_stress ", MatProp%threeelm_tau_stress_dim
           else
              call WriteFatalError('Unknown VE model for Attard, rappture_io')
           end if
           
        end if
        
        ! viscoelasticity ends here
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        
        if (( MatProp%fts_model == JKR) .or. ( MatProp%fts_model == DMT_DLVO)  .or. ( MatProp%fts_model == MORSE)  .or. ( MatProp%fts_model == LENNARD_JONES)) then
           !capillary with dlvo doesn't make sense (no cap in water).
           !cap w/ jkr disallowed for now b/c we don't implement it correctly
           MatProp%WantCapAd = .false.
        else
           MatProp%WantCapAd = daniel_get_boolean( path // ".group(nc).boolean(WantHYST).current")
           write(strVal, *) "WantCapAd ", MatProp%WantCapAd
           InputEcho = trim(InputEcho) // trim(strVal) // char(10)
        end if

        if (( MatProp%fts_model == JKR) .or. ( MatProp%fts_model == MORSE)  .or. ( MatProp%fts_model == LENNARD_JONES)) then
           MatProp%WantSurfHyst = .false.
        else
           MatProp%WantSurfHyst = daniel_get_boolean( path // ".group(nc).boolean(WantSurfHyst).current")           
           write(strVal, *) "WantSurfHyst ", MatProp%WantSurfHyst
           InputEcho = trim(InputEcho) // trim(strVal) // char(10)

           if (MatProp%WantSurfHyst) call read_value(driver, path // ".group(nc).number(SurfHystGamma).current", status, MatProp%SurfHyst_gamma, InputEcho, "Surf Hyst gamma " )
        end if

        MatProp%want_exp_dashpot =  daniel_get_boolean( path //".group(solvation).boolean(Want_exp_dashpot).current")

        if (MatProp%want_exp_dashpot) then
           strVal = rp_lib_get_wrap(driver, path // ".group(solvation).number(exp_dashpot_scale).current"  )
           status = status +  rp_units_convert_dbl(strVal," ",MatProp%exp_dashpot_scale)
           InputEcho = trim(InputEcho) // "exp_dashpot_scale " // trim(strVal) // char(10)

           strVal = rp_lib_get_wrap(driver, path // ".group(solvation).number(exp_dashpot_decay).current"  )
           status = status +  rp_units_convert_dbl(strVal," ",MatProp%exp_dashpot_decay)
           InputEcho = trim(InputEcho) // "exp_dashpot_decay " // trim(strVal) // char(10)
           MatProp%exp_dashpot_decay = MatProp%exp_dashpot_decay * 1e-9
        else
           MatProp%exp_dashpot_scale=0
           MatProp%exp_dashpot_decay=0
        end if
        
	
	if (MatProp%WantCapAd ) then
		strVal = rp_lib_get_wrap(driver, path // ".group(nc).number(D_0).current"  )
                InputEcho = trim(InputEcho) // "D_0 " // trim(strVal) // char(10)
        	status = status +  rp_units_convert_dbl(strVal," ", MatProp%D_0)
		MatProp%D_0 = MatProp%D_0/1d9

 		strVal = rp_lib_get_wrap(driver, path // ".group(nc).number(deltaE).current"  )
                InputEcho = trim(InputEcho) // "deltaE " // trim(strVal) // char(10)
        	status = status +  rp_units_convert_dbl(strVal," ", MatProp%deltaE)
		MatProp%deltaE = MatProp%deltaE/ electrons_per_coulomb                
	else
		MatProp%D_0 = 1d-9
		MatProp%deltaE = 0d0
	end if		


        call read_value(driver, path // ".number(Esample).current" , status, MatProp%Esample , InputEcho, "Esample "   )
 	MatProp%Esample = MatProp%Esample*1d9

!	Etip and Esample are the Young's modulus of the sample converted to (Pa)


        strVal = rp_lib_get_wrap(driver, path // ".number(Poisson_sample).current"  )   
        InputEcho = trim(InputEcho) // "Poisson_sample " // trim(strVal) // char(10)
        status = status +  rp_units_convert_dbl(strVal," ",MatProp%Poisson_sample)


!	DLVO Force parameters

        strVal = rp_lib_get_wrap( driver, path // ".number(KD).current")
        InputEcho = trim(InputEcho) // "KD " // trim(strVal) // char(10)
        status = status +  rp_units_convert_dbl(strVal," ",MatProp%KD_dim)
	MatProp%KD_dim = 1d6/MatProp%KD_dim

	strVal = rp_lib_get_wrap(driver, path // ".number(epsilon).current"  )
        InputEcho = trim(InputEcho) // "epsilon " // trim(strVal) // char(10)
        status = status +  rp_units_convert_dbl(strVal," ",MatProp%epsilon)

	strVal = rp_lib_get_wrap(driver, path // ".number(sigmat).current"  )
        InputEcho = trim(InputEcho) // "sigmat " // trim(strVal) // char(10)
        status = status +  rp_units_convert_dbl(strVal," ",MatProp%sigmat)

	strVal = rp_lib_get_wrap(driver, path // ".number(sigmas).current"  )
        InputEcho = trim(InputEcho) // "sigmas " // trim(strVal) // char(10)
        status = status +  rp_units_convert_dbl(strVal," ",MatProp%sigmas)
	

        !chadwick and the bottom edge correction models
        call read_value(driver, path // ".number(hs).current", status, MatProp%hs, InputEcho, "hs " )
        MatProp%hs = MatProp%hs*1d-9			
        

        if ( MatProp%fts_model == MORSE) then
           strVal = rp_lib_get_wrap(driver, path // ".number(MorseRc).current"  )
           InputEcho = trim(InputEcho) // "MorseRc " // trim(strVal) // char(10)
           status = status +  rp_units_convert_dbl(strVal," ",MatProp%MorseRc_dim)
           MatProp%MorseRc_dim= MatProp%MorseRc_dim*1e-9
           
           strVal = rp_lib_get_wrap(driver, path // ".number(MorseU0).current"  )
           InputEcho = trim(InputEcho) // "MorseU0 " // trim(strVal) // char(10)
           status = status +  rp_units_convert_dbl(strVal," ",MatProp%MorseU0)
           
           strVal = rp_lib_get_wrap(driver, path // ".number(MorseLambda).current"  )
           InputEcho = trim(InputEcho) // "MorseLambda " // trim(strVal) // char(10)
           status = status +  rp_units_convert_dbl(strVal," ",MatProp%MorseLambda_dim)
           MatProp%MorseLambda_dim= MatProp%MorseLambda_dim*1e-9
        elseif ((MatProp%fts_model == LENNARD_JONES) .or. (MatProp%fts_model == ATTARD_FOURIER_LSQ ) .or. ( MatProp%fts_model == ATTARD_FOURIER_BAHRAM)) then
           strVal = rp_lib_get_wrap(driver, path // ".number(LJ_r0).current"  )
           InputEcho = trim(InputEcho) // "LJ_r0 " // trim(strVal) // char(10)
           status = status +  rp_units_convert_dbl(strVal," ",MatProp%LJ_r0_dim)
           MatProp%LJ_r0_dim = MatProp%LJ_r0_dim * 1e-9
           
           strVal = rp_lib_get_wrap(driver, path // ".number(LJ_E0).current"  )
           InputEcho = trim(InputEcho) // "LJ_E0 " // trim(strVal) // char(10)
           status = status +  rp_units_convert_dbl(strVal," ",MatProp%LJ_E0)
        else
           MatProp%LJ_r0_dim=0
           MatProp%LJ_E0=0
        end if

        if ( MatProp%WantWLC .or. MatProp%WantOscillatory) then 
           call read_value(driver, path // ".number(temperature).current", status, MatProp%Temp, InputEcho, "Temperature " )
        end if

        if (MatProp%fts_model == CUSTOM_CONS) then
           strVal = rp_lib_get_wrap( driver, path // '.string(custom_filename).current')
           call readFtsTable_Custom_Cons(strVal, MatProp%N_custom, MatProp%gap_table_dim, MatProp%F_table_dim)
        else
           MatProp%N_custom = 0
        end if

        if (MatProp%fts_model == CUSTOM_MD) then
           strVal = rp_lib_get_wrap( driver, path // '.string(custom_filename).current')
           call readFtsTable_Custom_MD(strVal, MatProp%num_app_coeff, MatProp%num_ret_coeff, MatProp%num_x, MatProp%indentation_x, MatProp%ret_coeff, MatProp%app_coeff)
        else
           MatProp%N_custom = 0
        end if

     

        if (status > 0) call WriteFatalError( "Could not read sample and/or feature material properties.  Please check values. (path = " // path // ")" )

        return

100  call WriteFatalError( "Could not parse generalized maxwell modulii. Did you specify a comma separated list with the right number of values?")        
101  call WriteFatalError( "Could not parse generalized maxwell times. Did you specify a comma separated list with the right number of values?")        

      end subroutine ReadSampleData
      

subroutine read_value(driver, path, status, output_var, InputEcho, desc)
  character(len=*), intent(in) :: path, desc
  integer, intent(in) :: driver
  integer, intent(inout) :: status
  real*8, intent(out) :: output_var
  character*100:: strVal
  character*1500, intent(inout) :: InputEcho
  integer :: rp_units_convert_dbl

  strVal = rp_lib_get_wrap(driver, path)
  status = status +  rp_units_convert_dbl(strVal," ", output_var)
  InputEcho = trim(InputEcho) // desc // trim(strVal) // char(10)
end subroutine read_value

subroutine read_integer(driver, path, output_var, InputEcho, desc)
  character(len=*), intent(in) :: path, desc
  integer, intent(in) :: driver
  integer, intent(out) :: output_var
  character*100:: strVal
  character*1500, intent(inout) :: InputEcho
  integer, external :: readGenericInteger

  output_var =  readGenericInteger(driver, path)
  write(strVal,*) output_var
  InputEcho = trim(InputEcho) // desc // trim(strVal) // char(10)
end subroutine read_integer

subroutine ReadSimulationParameters(Zrange,Wanthist,WantHH,Want_Strob,Want_Impact,numpoinc,numHH,Z0,Zf,plotpnts,numincycle_dropdown, numincycle_direct,numHist,hist_cycles,xchoice, want_freqswp_sp_aprch, ASetPt, operating_mode, transient_allowance, want_hydrodynamic_function, Asp_f, freqshift_f, modulation_type, openmp_num_threads, InputEcho)

  implicit none
  logical, intent(out) :: Wanthist, WantHH, Want_Strob, Want_Impact,  want_freqswp_sp_aprch, want_hydrodynamic_function
  integer, intent(out) :: numHH, plotpnts, numincycle_dropdown, numincycle_direct, xchoice, numHist, hist_cycles, numpoinc,  Zrange, openmp_num_threads
  integer, intent(in) ::  operating_mode, modulation_type
  real*8, intent(out)  :: Z0, Zf, ASetPt, Asp_f, transient_allowance, freqshift_f
  character*1500, intent(inout) :: InputEcho

  character*100 strVal, tmpStr
  integer status, readGenericInteger, rp_units_convert_dbl, foo
  status = 0
  
  Zrange = readGenericInteger(driver, INPUT_PREFIX // "(op).choice(Zrange).current")
  write(tmpStr, *) Zrange
  InputEcho = trim(InputEcho) // "Zrange " // tmpStr // char(10)
  
  if ((operating_mode == FREQSWEEP) .or. (modulation_type == FREQUENCY)) then
     want_hydrodynamic_function =  daniel_get_boolean(INPUT_PREFIX // "(op).boolean(want_hydrodynamic_function).current")
  else
     want_hydrodynamic_function = .false.
  end if

  Wanthist =  daniel_get_boolean(INPUT_PREFIX // "(sim).group(timehist).boolean(Wanthist).current")

  WantHH =  daniel_get_boolean(INPUT_PREFIX // "(sim).boolean(WantHH).current")
  Want_Strob =  daniel_get_boolean(INPUT_PREFIX // "(sim).group(poinc).boolean(Want_Strob).current")
  Want_Impact =  daniel_get_boolean( INPUT_PREFIX // "(sim).group(poinc).boolean(Want_Impact).current")
		
  numHH = readGenericInteger(driver, INPUT_PREFIX // "(sim).integer(numHH).current")
  numpoinc = readGenericInteger(driver, INPUT_PREFIX // "(sim).group(poinc).integer(numpoinc).current")
		  
  if ((operating_mode /= SCAN) .or. (modulation_type == FORCE_VOL)) then
     strVal = rp_lib_get_wrap(driver, INPUT_PREFIX // "(op).number(Z0).current" )
     InputEcho = trim(InputEcho) // "Z0 " // StrVal // char(10)
     status = status +  rp_units_convert_dbl(strVal," ",Z0)
     Z0 = Z0/1d9
     
     if ((operating_mode == APPROACH ) .or. (operating_mode == APPROACH_STEP ) .or. (operating_mode == APPR_RET )) then 
        strVal = rp_lib_get_wrap(driver, INPUT_PREFIX // "(op).number(Zf).current" )
        InputEcho = trim(InputEcho) // "Zf " // StrVal // char(10)
        status = status +  rp_units_convert_dbl(strVal," ",Zf)
        Zf = Zf/1d9
        
        if (Zrange == ZRANGE_ASP) then
           strVal = rp_lib_get_wrap(driver, INPUT_PREFIX // "(op).number(Asp_f).current" )
           InputEcho = trim(InputEcho) // "Asp_f " // StrVal // char(10)
           status = status +  rp_units_convert_dbl(strVal," ",Asp_f)
        elseif (Zrange == ZRANGE_FREQSHIFT) then
           strVal = rp_lib_get_wrap(driver, INPUT_PREFIX // "(op).number(freqshift_f).current" )
           InputEcho = trim(InputEcho) // "freqshift_f " // StrVal // char(10)
           status = status +  rp_units_convert_dbl(strVal," ",freqshift_f)
        end if
     end if
  end if
  

  plotpnts = readGenericInteger( driver, INPUT_PREFIX // "(sim).integer(plotpnts).current")               
  if (plotpnts == 0) plotpnts = 1000 !backwards compat.
  

  numincycle_dropdown = readGenericInteger(driver, INPUT_PREFIX // "(sim).choice(numincycle).current")  
  numincycle_direct = readGenericInteger(driver, INPUT_PREFIX // "(sim).integer(numincycle).current")
 

   call read_integer( driver,  INPUT_PREFIX // "(sim).integer(openmp_num_threads).current", openmp_num_threads, InputEcho, "openmp_num_threads ")
  
  write( tmpStr, *) "plotpnts ", plotpnts, " numincycle_dropdown ", numincycle_dropdown, " numincycle_direct ", numincycle_direct
  InputEcho = trim(InputEcho) // trim(tmpStr) // char(10)

  if ((Wanthist) .or. (operating_mode == FIXED) ) then
     numHist = readGenericInteger(driver, INPUT_PREFIX // "(sim).group(timehist).integer(numHist).current")
          
     hist_cycles = readGenericInteger(driver, INPUT_PREFIX // "(sim).group(timehist).integer(Nhist).current") 
     write(tmpStr, *) "Hhist(cycles) ", hist_cycles
     InputEcho = trim(InputEcho) //  tmpStr // char(10)
  else
     numHist = 0
  end if
  
      if (operating_mode == FREQSWEEP) then
         want_freqswp_sp_aprch = daniel_get_boolean( INPUT_PREFIX // "(op).boolean(WantSetPt).current")
      else
         want_freqswp_sp_aprch = .false.   
      end if
      
      if ((operating_mode == SCAN) .or. (want_freqswp_sp_aprch)) then
         if ((modulation_type .ne. FORCE_VOL) .and.(modulation_type .ne. PEAK_FORCE))  then
            strVal = rp_lib_get_wrap(driver,  INPUT_PREFIX // "(op).number(Asp).current")
            InputEcho = trim(InputEcho) // "Asp " // StrVal // char(10)
            status = status +  rp_units_convert_dbl(strVal," ",ASetPt)
         end if
      end if
 
      xchoice = readGenericInteger(driver, INPUT_PREFIX // "(sim).choice(xchoice).current")
      if (status > 0) call WriteFatalError( "Error reading simulation parameters. Please check values1")  
      transient_allowance = -1d0
      call read_value(driver, INPUT_PREFIX // "(sim).number(transient).current", foo, transient_allowance, InputEcho, "transient_allowance " ) !no status on this one, backwards compatability
      
      if (status > 0) call WriteFatalError( "Error reading simulation parameters. Please check values2")
end subroutine ReadSimulationParameters
 
subroutine readForceVol(LineSpeed, F_ForceVol_dim, ForceVolSettleTime, z_feedback_choice, InputEcho)
  character*1500, intent(inout) :: InputEcho
  real*8, intent(out) :: F_ForceVol_dim, ForceVolSettleTime, LineSpeed
  character*100 strVal, tmpStr
  integer, intent(out) :: z_feedback_choice
  integer status
  integer :: rp_units_convert_dbl, readGenericInteger
  status = 0
   strVal = rp_lib_get_wrap(driver,  INPUT_PREFIX // "(op).number(F_ForceVol_dim).current")
   InputEcho = trim(InputEcho) // "F_ForceVol_dim " // StrVal // char(10)
   status = status +  rp_units_convert_dbl(strVal," ",F_ForceVol_dim)
   F_ForceVol_dim = F_ForceVol_dim / 1d9
   strVal = rp_lib_get_wrap(driver,  INPUT_PREFIX // "(op).number(ForceVolSettleTime).current")
   InputEcho = trim(InputEcho) // "ForceVolSettleTime " // StrVal // char(10)
   status = status +  rp_units_convert_dbl(strVal," ",ForceVolSettleTime)
   call read_value(driver, INPUT_PREFIX // "(op).number(LineSpeed).current", status, LineSpeed, InputEcho, "LineSpeed "   )
   z_feedback_choice = readGenericInteger(driver, INPUT_PREFIX // "(op).choice(Z_feedback_choice).current")
  
  write( tmpstr, *)   Z_feedback_choice
  InputEcho = trim(InputEcho) // "Z_feedback_choice " // trim( tmpstr) // char(10)

   if (status>0) call WriteFatalError("Could not read  force Colume parameters")

 end subroutine readForceVol
 

!get data about multiple higher harmonics
subroutine Read_HH( numHH, NHH, z_feedback_choice)
  integer, intent(in) ::  numHH
  integer, intent(in) :: z_feedback_choice
  real*8, intent(out) :: NHH(numHH) ! jtm - wanted to allow abitrary freq.
  character*100 strVal
  integer i
	
  strVal = rp_lib_get_wrap(driver, INPUT_PREFIX // "(sim).string(NHH1).current" )
  if (z_feedback_choice == MEAN_DEFL_Z) then
     if (NumHH>1) then
        read( strVal, *, end=994, err=994) (NHH(i),i=2,numHH)
        NHH(1) = 0
     elseif (NumHH == 1) then
        NHH = 0
     end if
 else  
    read( strVal, *, end=994, err=994) (NHH(i),i=1,numHH)
 end if
  return

994  call WriteFatalError( "Could not parse desired higher harmonics. Did you specify a comma separated list with the right number of values?")
end subroutine Read_HH

!get data about multiple time histories 
subroutine Read_TimeHist( numHist, Ahist, Wanthist_byA, Want_TH_obs_defl, want_TH_force_gap,do_fft,want_acceleration)
  integer, intent(in) ::  numHist
  real*8, intent(out) :: Ahist(numHist)
  logical, intent(out) :: Wanthist_byA, Want_TH_obs_defl, want_TH_force_gap, do_fft, want_acceleration
  character*100 strVal
  integer i

  Wanthist_byA = daniel_get_boolean( INPUT_PREFIX // "(sim).group(timehist).boolean(Wanthist_byA).current")

  do_fft=daniel_get_boolean (INPUT_PREFIX // "(sim).group(timehist).boolean(do_fft).current")

  Want_TH_obs_defl = daniel_get_boolean( INPUT_PREFIX // "(sim).group(timehist).boolean(Want_TH_obs_defl).current")
  Want_TH_force_gap = daniel_get_boolean( INPUT_PREFIX // "(sim).group(timehist).boolean(Want_TH_force_gap).current")
 
  want_acceleration = daniel_get_boolean( INPUT_PREFIX // "(sim).group(timehist).boolean(Want_TH_acceleration).current")
 
  strVal = rp_lib_get_wrap(driver, INPUT_PREFIX // "(sim).group(timehist).string(Ahist1).current" )
  read( strVal, *, end=993, err=993) (Ahist(i),i=1,numHist)
        
  return

993  call WriteFatalError( "Could not parse desired time history outputs. Did you specify a comma separated list with the right number of values?")

end subroutine Read_TimeHist


!New Subroutine added to read plot choice booleans
subroutine ReadPlotChoices(Want_A2,Want_P2,Want_AZ,Want_P1,Want_MF,Want_PF,Want_ED,Want_EP,Want_I,Want_CT, Want_Fourier, Want_ContTrans, Want_E_Anc, Want_ForceFourier, want_ev, want_virial, want_RMS)
  logical, intent(out) :: Want_AZ, Want_P1, Want_MF, Want_PF, Want_ED, Want_EP, Want_I, Want_CT, want_ev, Want_RMS
  logical, intent(out) :: Want_A2, Want_P2, Want_Fourier, Want_ContTrans, Want_E_Anc, Want_ForceFourier, want_virial


    Want_ev =           daniel_get_boolean( INPUT_PREFIX // "(sim).group(plots).boolean(Want_ev).current")
    Want_RMS =          daniel_get_boolean( INPUT_PREFIX // "(sim).group(plots).boolean(Want_RMS).current")
    Want_virial =       daniel_get_boolean( INPUT_PREFIX // "(sim).group(plots).boolean(Want_virial).current")
    Want_E_Anc =        daniel_get_boolean( INPUT_PREFIX // "(sim).group(plots).boolean(Want_E_Anc).current")
    Want_AZ =           daniel_get_boolean( INPUT_PREFIX // "(sim).group(plots).boolean(Want_AZ).current")
    Want_P1 =           daniel_get_boolean( INPUT_PREFIX // "(sim).group(plots).boolean(Want_P1).current")
    Want_MF =           daniel_get_boolean( INPUT_PREFIX // "(sim).group(plots).boolean(Want_MF).current")
    Want_PF =           daniel_get_boolean( INPUT_PREFIX // "(sim).group(plots).boolean(Want_PF).current")
    Want_ED =           daniel_get_boolean( INPUT_PREFIX // "(sim).group(plots).boolean(Want_ED).current")
    Want_EP =           daniel_get_boolean( INPUT_PREFIX // "(sim).group(plots).boolean(Want_EP).current")
    Want_I =            daniel_get_boolean( INPUT_PREFIX // "(sim).group(plots).boolean(Want_I).current")
    Want_CT =           daniel_get_boolean( INPUT_PREFIX // "(sim).group(plots).boolean(Want_CT).current")
    Want_A2 =           daniel_get_boolean( INPUT_PREFIX // "(sim).group(plots).boolean(Want_A2).current")
    Want_P2 =           daniel_get_boolean( INPUT_PREFIX // "(sim).group(plots).boolean(Want_P2).current")
    Want_Fourier =      daniel_get_boolean( INPUT_PREFIX // "(sim).group(plots).boolean(Want_Fourier).current")
    Want_ContTrans =    daniel_get_boolean( INPUT_PREFIX // "(sim).boolean(Want_ContTrans).current")
    Want_ForceFourier = daniel_get_boolean( INPUT_PREFIX // "(sim).group(plots).boolean(Want_ForceFourier).current")
end subroutine ReadPlotChoices

subroutine ReadPoincareParameters(Want_Strob,Want_Impact,numpoinc,Npoinc_X,Npoinc_Y)
  use Poincare, only : maxpoinc
  logical, intent(in) :: Want_Strob, Want_Impact
  integer, intent(in) :: numpoinc
  integer, intent(out) :: Npoinc_X(maxpoinc), Npoinc_Y(maxpoinc)

  character*100 strVal
  integer i
  
  if ((Want_Strob) .or. (Want_Impact)) then
	
     strVal = rp_lib_get_wrap(driver, INPUT_PREFIX // "(sim).group(poinc).string(Npoinc_X1).current" )
     read( strVal, *, end=995, err=995) (Npoinc_X(i),i=1,numpoinc)
	    
     strVal = rp_lib_get_wrap(driver, INPUT_PREFIX // "(sim).group(poinc).string(Npoinc_Y1).current" )
     read( strVal, *, end=995, err=995) (Npoinc_Y(i),i=1,numpoinc)
     
  end if
  
  return

995  call WriteFatalError( "Could not parse poincare data. Did you specify comma separated lists with one value per plot?")

end subroutine ReadPoincareParameters

subroutine ReadOperatingParameter( fexcite, exc_choice, sweepchoice, numModes, operating_mode, CalcInputK, AutoCalcOmega, modulation_type, AutoCalcAlpha, output_type, InputEcho)
  integer, intent(out) :: fexcite, exc_choice, numModes, operating_mode, sweepchoice
  integer, intent(out) :: modulation_type, output_type
  logical, intent(out) :: CalcInputK, AutoCalcOmega,  AutoCalcAlpha
  character*1500, intent(inout) :: InputEcho
  character*100 :: tmpStr
  integer :: readGenericInteger
  logical ::  AutoCalcChi

  !primary operating mode APPROACH = 1, FREQUENCY SWEEP = 2, SCAN = 3
  operating_mode = readGenericInteger(driver, INPUT_PREFIX // "(op).choice(operating_mode).current")

   !	  Choice of excitation 
  fexcite = readGenericInteger(driver, INPUT_PREFIX // "(op).choice(fexcite).current")


  modulation_type = readGenericInteger(driver, INPUT_PREFIX // "(op).choice(modulation_type).current")

  CalcInputK =  daniel_get_boolean( INPUT_PREFIX // "(op).boolean(CalcInputK).current")

  
  output_type = readGenericInteger(driver, INPUT_PREFIX // "(op).choice(output_type).current")
  !backwards compatibility
  if (output_type == 0) then
     AutoCalcChi = daniel_get_boolean( INPUT_PREFIX // "(op).boolean(AutoCalcChi).current")
     if (AutoCalcChi) then
        output_type = AUTO_CALC_CHI
     else
        output_type = MAN_CALC_CHI
     end if
  end if


  AutoCalcAlpha = daniel_get_boolean( INPUT_PREFIX // "(op).boolean(AutoCalcAlpha).current")
   
  AutoCalcOmega = daniel_get_boolean( INPUT_PREFIX // "(op).boolean(AutoCalcOmega).current")
  
  if ( operating_mode .ne. FREQSWEEP) then
     !(1) single mode, (2) bimodal
     exc_choice = readGenericInteger(driver, INPUT_PREFIX // "(op).choice(freqchoice).current")     
     write( tmpStr, *) "exc_choice = ", exc_choice     
  else
     sweepchoice = readGenericInteger(driver, INPUT_PREFIX // "(op).choice(freqchoice).current")
     exc_choice = SINGLE
     write( tmpStr, *) "freq_choice = ", sweepchoice
  end if
  
  numModes =  readGenericInteger(driver, INPUT_PREFIX // "(op).number(numModes).current")

  InputEcho = trim( InputEcho) // tmpStr // char(10)
  write(tmpStr, *) "operating mode ", operating_mode, " fexcite ", fexcite, " numModes ", numModes
  
  InputEcho = trim( InputEcho) // tmpStr // char(10)

end subroutine ReadOperatingParameter

subroutine ReadOpCondAndCantProps(modulation_type, operating_mode, exc_choice, fexcite,numModes, omegad_dim, omegad_start, omegad_stop, sweep_time, Ainitial_dim, AprchS_dim, mtip, want_AutoCalcTC, LockInOrder, LockInTC, gamma_drag, want_Abase_direct, Abase_input, F_input, want_autoCalcIC, Y_IC, Afluid, mstar_div_m, selfexc_gain, selfexc_phase, selfexc_saturation,&
Asample_dim, omegas_dim, wantSampleExc, InputEcho)
    
  integer, intent(in) ::  operating_mode, exc_choice, numModes, modulation_type
  integer, intent(inout) ::  fexcite
  real*8, intent(out) :: omegad_dim(maxExc), omegad_start, omegad_stop, sweep_time, Ainitial_dim(maxExc), AprchS_dim, mtip, LockInTC, gamma_drag, Abase_input(maxExc), F_input(maxModes), Y_IC(maxModes), mstar_div_m,  selfexc_gain, selfexc_phase, selfexc_saturation, Asample_dim, omegas_dim
  complex*8, intent(out) :: Afluid
  integer, intent(out) :: LockInOrder
  logical, intent(out) :: want_Abase_direct, want_autoCalcIC, want_AutoCalcTC, wantSampleExc
  character*1500, intent(inout) :: InputEcho
  
  logical :: want_Afluid
  character*100 strVal
  integer status, rp_units_convert_dbl, i, readGenericInteger
  real*8 :: tmp, Afluid_real, Afluid_imag

  status = 0
  if (modulation_type /= FORCE_VOL) then
     want_AutoCalcTC = daniel_get_boolean(INPUT_PREFIX // "(op).boolean(Want_AutoCalcTC).current")
  end if
  if (.not. want_AutoCalcTC) then
     LockInTC = readGenericInteger( driver, INPUT_PREFIX // "(op).choice(LockInTC).current")
     write(strVal, *) LockInTC
     InputEcho = trim(InputEcho) // "LockInTC (us) " // trim(strVal) // char(10)
     LockInTC = LockInTC / 1d6
  end if

  LockInOrder = readGenericInteger(driver,INPUT_PREFIX // "(op).choice(LockInOrder).current")  
  write(strVal, *) LockInOrder
  InputEcho = trim(InputEcho) // "LockInOrder " // trim(strVal) // char(10)

  strVal = rp_lib_get_wrap(driver, INPUT_PREFIX // "(op).number(mtip).current" )
  InputEcho = trim(InputEcho) // "mtip " // trim(strVal) // char(10)
  status = status +  rp_units_convert_dbl(strVal," ",mtip)



  if (exc_choice == SELFEXC) then
     strVal = rp_lib_get_wrap(driver, INPUT_PREFIX // "(op).group(selfexc).number(selfexc_gain).current" )
     InputEcho = trim(InputEcho) // "selfexc_gain " // trim(strVal) // char(10)
     status = status +  rp_units_convert_dbl(strVal," ",selfexc_gain)

     strVal = rp_lib_get_wrap(driver, INPUT_PREFIX // "(op).group(selfexc).number(selfexc_sat).current" )
     InputEcho = trim(InputEcho) // "selfexc_sat " // trim(strVal) // char(10)
     status = status +  rp_units_convert_dbl(strVal," ",selfexc_saturation)

     strVal = rp_lib_get_wrap(driver, INPUT_PREFIX // "(op).group(selfexc).number(selfexc_phase).current" )
     InputEcho = trim(InputEcho) // "selfexc_phase " // trim(strVal) // char(10)
     status = status +  rp_units_convert_dbl(strVal," ",selfexc_phase)
  end if

  !drive frequencies
  if ( operating_mode /= FREQSWEEP ) then 
  
     strVal = rp_lib_get_wrap(driver, INPUT_PREFIX // "(op).string(omegad_string).current" )
     InputEcho = trim(InputEcho) // "omegad " // trim(strVal) // char(10)

     if (exc_choice == BIMODAL) then
        read( strVal, *, end=996, err=996) (omegad_dim(i),i=1,2)
        omegad_dim = omegad_dim*2*pi*1d3
     else     		
        read( strVal, *, end=996, err=996) (omegad_dim(i),i=1,1)
        omegad_dim(1) = omegad_dim(1)*2*pi*1d3
     end if

  else
     strVal = rp_lib_get_wrap(driver, INPUT_PREFIX // "(op).number(omegad_start).current" )
     InputEcho = trim(InputEcho) // "omegad_start " // trim(strVal) // char(10)
     status = status +  rp_units_convert_dbl(strVal," ", tmp)
     omegad_start = tmp *2*pi*1d3 !rad/s
     
     strVal = rp_lib_get_wrap(driver, INPUT_PREFIX // "(op).number(omegad_stop).current" )
     InputEcho = trim(InputEcho) // "omegad_stop " // trim(strVal) // char(10)
     status = status +  rp_units_convert_dbl(strVal," ",tmp)
     omegad_stop = tmp *2*pi*1d3
     
     strVal = rp_lib_get_wrap(driver, INPUT_PREFIX // "(op).number(sweep_time).current" )
     InputEcho = trim(InputEcho) // "sweep_time " // trim(strVal) // char(10)
     status = status +  rp_units_convert_dbl(strVal," ", sweep_time)
     
  end if

  want_Afluid = daniel_get_boolean(INPUT_PREFIX // "(op).boolean(want_Afluid).current")
  if (want_Afluid) then

     strVal = rp_lib_get_wrap(driver, INPUT_PREFIX // "(op).number(Afluid_real).current" )
     InputEcho = trim(InputEcho) // "Afluid_real " // trim(strVal) // char(10)
     status = status +  rp_units_convert_dbl(strVal," ", Afluid_real)

     strVal = rp_lib_get_wrap(driver, INPUT_PREFIX // "(op).number(Afluid_imag).current" )
     InputEcho = trim(InputEcho) // "Afluid_imag " // trim(strVal) // char(10)
     status = status +  rp_units_convert_dbl(strVal," ", Afluid_imag)

     Afluid = CMPLX(Afluid_real, Afluid_imag)
     
     strVal = rp_lib_get_wrap(driver, INPUT_PREFIX // "(op).number(mstar_div_m).current" )
     InputEcho = trim(InputEcho) // "mstar_div_m " // trim(strVal) // char(10)
     status = status +  rp_units_convert_dbl(strVal," ", mstar_div_m)     
  else
     Afluid = 0d0
     mstar_div_m = 0d0
  end if
   
  !free amplitudes
!  want_Abase_direct = daniel_get_boolean(  INPUT_PREFIX // "(op).boolean(want_Abase_direct).current"))
  want_Abase_direct = daniel_get_boolean(INPUT_PREFIX // "(op).boolean(want_Abase_direct).current")
  if (want_Abase_direct) then
     strVal = rp_lib_get_wrap(driver, INPUT_PREFIX // "(op).string(Abase_direct).current" )
     if (fexcite == ACOUSTIC_IDEAL) then
        if (exc_choice == BIMODAL) then
           read( strVal, *, end=997, err=997) (Abase_input(i),i=1,2)
        else
           read( strVal, *, end=997, err=997) (Abase_input(i),i=1,1)
        end if
        inputEcho = trim(inputEcho) // "Abase (direct) = " // trim(strVal) // char(10)
     elseif (fexcite == ACOUSTIC_PIEZOR) then
        call assert(.false., 'fixme 1655')
     elseif (isMagnetic(fexcite)) then
        read( strVal, *, end=998, err=998) (F_input(i),i=1,numModes)
        inputEcho = trim(inputEcho) // "modal forces (direct) = " // trim(strVal) // char(10)
     elseif (fexcite == SAMPLE) then
        !for compatibility with force modulation mode.  ignore
        Abase_input = 0
     else
        call assert(.false., 'unhandled direct input')
     end if

     Abase_input = Abase_input / 1d9

     if (modulation_type == FREQUENCY) then
        Ainitial_dim(1) = readGenericDbl(status, INPUT_PREFIX // "(op).number(Ainitial).current", "Ainitial", InputEcho, 1d9) 
        !have to have a setpoint for FM controller and this is it.  probably could just calculate the free amp for the given inputs and use that, but this is easier
     else
        Ainitial_dim(1) = 1d-9  !have to have something for normalization
     end if
  else
     Ainitial_dim(1) = readGenericDbl(status, INPUT_PREFIX // "(op).number(Ainitial).current", "Ainitial", InputEcho, 1d9)
     if (exc_choice == BIMODAL) then
        Ainitial_dim(2) = readGenericDbl(status, INPUT_PREFIX // "(op).number(Ainitial2).current", "Ainitial2", InputEcho, 1d9)
     else
        Ainitial_dim(2)=0
     end if

     Abase_input(1) = 0
  end if

  if (fexcite == SAMPLE) then
     !for force mod tool
     fexcite = NO_EXC
     wantSampleExc = .true.
     omegas_dim =omegad_dim(1)
  else
     !for ssac tool
     wantSampleExc = daniel_get_boolean( INPUT_PREFIX // "(op).boolean(WantSampleExc).current")
     if (wantSampleExc) omegas_dim = readGenericDbl( status, INPUT_PREFIX // "(op).number(omegas).current", 'omegas', InputEcho, 2d0*pi*1d3)
  end if

  if (wantSampleExc) then
     Asample_dim = readGenericDbl( status, INPUT_PREFIX // "(op).number(Asample).current", "Asample", InputEcho, 1d9)
  else
     omegas_dim  = 0d0
     Asample_dim = 0d0
  end if

  if (.not. exc_choice == BIMODAL ) then
     Ainitial_dim(2) = 0
     Abase_input(2) = 0
  end if
  
  if (operating_mode == FIXED) then
     want_autoCalcIC = .not. daniel_get_boolean(  INPUT_PREFIX // "(op).boolean(want_direct_IC).current")
  else
     want_autoCalcIC = .true.
  end if
  if (.not. want_autoCalcIC) then
     strVal = rp_lib_get_wrap(driver, INPUT_PREFIX // "(op).string(IC).current" )
     read( strVal, *, end=999, err=999) (Y_IC(i),i=1,2*numModes)
  end if
!!!!!!!!!!!!

  if ((isOpModeApp(operating_mode)) .or. (modulation_type == FORCE_VOL) .or. (modulation_type == PEAK_FORCE)) then 
     strVal = rp_lib_get_wrap(driver, INPUT_PREFIX // "(op).number(gamma_drag).current" )
     InputEcho = trim(InputEcho) // "gamma_drag " // trim(strVal) // char(10)
     status = status +  rp_units_convert_dbl(strVal," ", gamma_drag)

     if (modulation_type /= PEAK_FORCE) then
        AprchS_dim = readGenericDbl(status, INPUT_PREFIX // "(op).number(AprchS).current", "Approach speed", InputEcho, 1d9)
     end if
  else
     AprchS_dim = 0
     gamma_drag = 0
  end if
  if (status > 0) call WriteFatalError( "Error reading operating conditions or cantilever properties. Please check values.")   
 
 return

996  call WriteFatalError( "Could not parse drive frequencies. Did you specify a comma separated list with two values?  Or perhaps you specified bimodal accidentally?")
997  call WriteFatalError( "Could not parse base amplitudes. Did you specify a comma separated list with two values?  Or perhaps you specified bimodal accidentally?")
998  call WriteFatalError( "Could not parse modal forces. Did you specify a comma separated list with one value per mode?")
999  call WriteFatalError( "Could not parse initial conditions. Did you specify a comma separated list with one displacement and one velocity per eigenmode?")

 end subroutine ReadOpCondAndCantProps


subroutine readFeatureProperties( FeatureType, HF, LF, LF2, SubsLen, WantTSCON, WantFProp, inputEcho)
  real*8, intent(out) :: HF, LF, LF2, SubsLen
  integer, intent(out) :: FeatureType
  logical, intent(out) :: WantTSCON, WantFprop
  character*1500, intent(inout) :: InputEcho

  character*100 strVal
  integer :: readGenericInteger
  integer status, rp_units_convert_dbl

  status = 0

  FeatureType = readGenericInteger(driver,INPUT_PREFIX // "(feature).choice(FeatureType).current")  
  !       step (1), trapizoid (2), sindusoid (3) 	
  	
  strVal = rp_lib_get_wrap(driver, INPUT_PREFIX // "(feature).number(HF).current")
  InputEcho = trim(InputEcho) // "HF " // trim(strVal) // char(10)
  status = status +  rp_units_convert_dbl(strVal," ",HF)
  HF = HF/1d9
  
  strVal = rp_lib_get_wrap(driver, INPUT_PREFIX // "(feature).number(LF).current")
  InputEcho = trim(InputEcho) // "LF " // trim(strVal) // char(10)
  status = status +  rp_units_convert_dbl(strVal," ",LF)
  LF = LF/1d9
  
  strVal = rp_lib_get_wrap(driver, INPUT_PREFIX // "(feature).number(LF2).current")
  InputEcho = trim(InputEcho) // "LF2 " // trim(strVal) // char(10)
  status = status +  rp_units_convert_dbl(strVal," ",LF2)
  LF2 = LF2/1d9
  
  strVal = rp_lib_get_wrap(driver, INPUT_PREFIX // "(sim).number(LS).current")
  InputEcho = trim(InputEcho) // "LS " // trim(strVal) // char(10)
  status = status +  rp_units_convert_dbl(strVal," ", SubsLen)
  SubsLen = SubsLen/1d9
 
  WantTSCON = daniel_get_boolean( INPUT_PREFIX // "(feature).boolean(WantTSCON).current")
  
  WantFprop = daniel_get_boolean( INPUT_PREFIX // "(feature).boolean(WantFprop).current")
  
  write(strVal, *) "Feature type ", FeatureType, "WantTSCon ", WantTSCon
  InputEcho = trim(InputEcho) // trim(strVal) // char(10)

  if (status > 0) call WriteFatalError( "Error reading feature properties. Please check values")

end subroutine readFeatureProperties

subroutine ReadSampleFreq( sample_freq, InputEcho)
  real*8, intent(out) :: sample_freq
  character*1500, intent(inout) :: InputEcho

  character*100 strVal
  integer status, rp_units_convert_dbl
  
  status = 0

  strVal =  rp_lib_get_wrap(driver,INPUT_PREFIX // "(op).number(sample_freq).current")
  InputEcho = trim(InputEcho) // "sample_freq Mhz " // trim(strVal) // char(10)
  status = status +  rp_units_convert_dbl(strVal," ", sample_freq)
  sample_freq = sample_freq * 1d6

  write( strVal, *)  "sample_freq ", sample_freq
  InputEcho = trim(InputEcho) // trim(strVal) // char(10)

  if (status > 0) call WriteFatalError( "Error reading sample frequency. Please check values") 
end subroutine ReadSampleFreq


subroutine readTriggeredFzMode(  jumpDeflLim, PiezoReverseTime_pct , Want_FzEig, InputEcho, fzshape, CurveRate_dim)
  character*1500, intent(inout) :: InputEcho
  real*8, intent(out) :: jumpDeflLim,  PiezoReverseTime_pct, CurveRate_dim
  logical, intent(out) :: Want_FzEig
  integer, intent(out) :: fzshape
  character*100 strVal
  integer :: status, rp_units_convert_dbl

  call read_integer(driver, INPUT_PREFIX // "(op).choice(fzshape).current" , fzshape , InputEcho, "fzshape")

  strVal = rp_lib_get_wrap(driver,  INPUT_PREFIX // "(op).number(PiezoReverseTime_pct).current")
  InputEcho = trim(InputEcho) // "PiezoReverseTime_pct " // trim(strVal) // char(10)
  status = rp_units_convert_dbl(strVal," ",PiezoReverseTime_pct)

  strVal = rp_lib_get_wrap(driver,  INPUT_PREFIX // "(op).number(jumpDeflLim).current")
  InputEcho = trim(InputEcho) // "jumpDeflLim " // trim(strVal) // char(10)
  status = rp_units_convert_dbl(strVal," ",jumpDeflLim)
  jumpDeflLim = jumpDeflLim / 1e9;

  Want_FzEig = daniel_get_boolean( INPUT_PREFIX // "(sim).group(plots).boolean(Want_FzEig).current")
  write( strVal, *) Want_FzEig
  InputEcho = trim(InputEcho) // "want_FzEig " // strVal // char(10)

  CurveRate_dim = readGenericDbl(status, INPUT_PREFIX // "(op).number(CurveRate).current", "Curve Rate", InputEcho, 1d0)
  
  if (status > 0) call WriteFatalError( "Error reading FZ curves properties. Please check values")   
  
end subroutine readTriggeredFzMode

subroutine readFtsTable_Custom_cons( filename, N_custom, gap_table_dim, F_table_dim)
  integer, intent(out) :: N_custom
  character(len=*), intent(in) :: filename
  real*8, allocatable, intent(inout) :: gap_table_dim(:), F_table_dim(:)
  integer :: i

  open( UNIT=10, FILE=filename)

  read(10, *) N_custom

!  write(*,*) N_custom

  allocate( gap_table_dim(N_custom))
  allocate( F_table_dim(N_custom))

  do i = 1, N_custom
     read(10, *) gap_table_dim(i), F_table_dim(i)
!     write(*,*)  gap_table_dim(i), F_table_dim(i)
  end do

  close(UNIT=10)
end subroutine readFtsTable_Custom_cons


subroutine readFtsTable_Custom_MD( filename, num_app_coeff,num_ret_coeff, num_x, indentation_x, ret_coeff, app_coeff )
  integer, intent(out) :: num_app_coeff, num_ret_coeff, num_x
  character(len=*), intent(in) :: filename
  real*8, allocatable, intent(inout) :: indentation_x(:), ret_coeff(:,:), app_coeff(:)
  integer :: i, j

  open( UNIT=10, FILE=filename)

  read(10, *) num_x, num_app_coeff, num_ret_coeff

  allocate( indentation_x(num_x))
  allocate( app_coeff(num_app_coeff))
  allocate( ret_coeff(num_ret_coeff, num_x))

  !assume same order as default matlab polyfit.  first value is x^n,  last value is x^0
  !x-coordinates are in nanometers (for now)
  if (num_app_coeff >0) read(10,*) ( app_coeff(i), i = 1,num_app_coeff)

  do i = 1, num_x
     !indentation_x, x_coordinates in meters, indentation negative
     read (10,*) indentation_x(i), (ret_coeff(j,i), j = 1,num_ret_coeff)
  end do

  close(UNIT=10)
end subroutine readFtsTable_Custom_MD


subroutine readCantModalProps(  numModes, fexcite, omegai,  Keq, Quality, CalcInputK, output_type, Chi, AutoCalcOmega, AutoCalcAlpha, alpha, osc_Q, osc_omega,phase_slope_dim, efficiency_slope_dim, want_nonideal_magnetic, want_custom_B, B_input, inputEcho)
  integer, intent(in) ::  numModes, fexcite, output_type
  logical, intent(in) :: CalcInputK, AutoCalcOmega, AutoCalcAlpha
  real*8, intent(out) :: omegai(numModes), Keq(numModes), Quality(numModes), Chi(numModes)
  real*8, intent(out) :: alpha(numModes), osc_Q, osc_omega, phase_slope_dim, efficiency_slope_dim, B_input(numModes)
  logical, intent(out) :: want_nonideal_magnetic, want_custom_B
  character*1500, intent(inout):: inputEcho
  character*100 strVal
  integer i
  integer status, rp_units_convert_dbl

  status = 0

  if (isMagnetic(fexcite)) then
     want_custom_B = daniel_get_boolean( INPUT_PREFIX // "(op).boolean(want_custom_B).current")
     if (want_custom_B) then
        strVal = rp_lib_get_wrap(driver, INPUT_PREFIX // "(op).string(B_input).current" )
        read( strVal, *, end=988, err=988) (B_input(i),i=1,numModes)
        inputEcho = trim(inputEcho) // "B = " // trim(strVal) // char(10)
     end if
  else
     want_custom_B = .false.
  end if


  if ( .not. AutoCalcAlpha) then
     strVal = rp_lib_get_wrap(driver, INPUT_PREFIX // "(op).string(alpha).current" )
     read( strVal, *, end=989, err=989) (alpha(i),i=1,numModes)
     inputEcho = trim(inputEcho) // "alpha = " // trim(strVal) // char(10)
  end if
  

  !natural frequencies !for now we'll keep all the labels as omega1, quality1, etc. so avoid having to change the gui simultaneously
  strVal = rp_lib_get_wrap(driver, INPUT_PREFIX // "(op).string(omega1).current" )
  if (AutoCalcOmega) then
     read( strVal, *, end=990, err=990) omegai(1)  
  else
     read( strVal, *, end=990, err=990) (omegai(i),i=1,numModes)  !parses a comma separated listed into the array values
  end if
  
  inputEcho = trim(inputEcho) // "omega = " // trim(strVal) // char(10)

  strVal = rp_lib_get_wrap(driver, INPUT_PREFIX // "(op).string(Keq1).current" )
  if (CalcInputK) then
     read( strVal, *, end=991, err=991) Keq(1)
  else
     read( strVal, *, end=991, err=991) (Keq(i),i=1,numModes)
  end if
  inputEcho = trim(inputEcho) // "Keq = " // trim(strVal) // char(10)
      
  if (output_type == MAN_CALC_CHI) then
     strVal = rp_lib_get_wrap(driver, INPUT_PREFIX // "(op).string(Chi2_n).current" )
     read( strVal, *, end=993, err=993) (Chi(i),i=1,numModes)
     inputEcho = trim(inputEcho) //  "Chi = " // trim(strVal) // char(10)
  else
     inputEcho = trim(inputEcho) //  "Chi = autocalc " // char(10)
  end if

  strVal = rp_lib_get_wrap(driver, INPUT_PREFIX // "(op).string(Quality1).current" )
  read( strVal, *, end=992, err=992) (Quality(i),i=1,numModes)
  inputEcho = trim(inputEcho) // "Q = " // trim(strVal) // char(10)
  
  if (fexcite == ACOUSTIC_PIEZOR) then
     strVal = rp_lib_get_wrap(driver, INPUT_PREFIX // "(op).number(osc_Q).current" )
     status = status +  rp_units_convert_dbl(strVal," ", osc_Q)
     inputEcho = trim(inputEcho) // "osc_Q = " // trim(strVal) // char(10)
     
     strVal = rp_lib_get_wrap(driver, INPUT_PREFIX // "(op).number(osc_omega).current" )
     status = status +  rp_units_convert_dbl(strVal," ", osc_omega)
     inputEcho = trim(inputEcho) // "osc_omega = " // trim(strVal) // char(10)
  end if

  if (isMagnetic(fexcite)) then
     want_nonideal_magnetic = daniel_get_boolean( INPUT_PREFIX // "(op).boolean(want_nonideal_magnetic).current")
  else
     want_nonideal_magnetic = .false.
  end if

  if ( isMagnetic(fexcite) .and. want_nonideal_magnetic) then
     strVal = rp_lib_get_wrap(driver, INPUT_PREFIX // "(op).number(phase_slope).current" )
     status = status +  rp_units_convert_dbl(strVal," ", phase_slope_dim)
     inputEcho = trim(inputEcho) // "phototh_phase_slope = " // trim(strVal) // char(10)
     
     strVal = rp_lib_get_wrap(driver, INPUT_PREFIX // "(op).number(efficiency_slope).current" )
     status = status +  rp_units_convert_dbl(strVal," ", efficiency_slope_dim)
     inputEcho = trim(inputEcho) // "phototh_efficiency_slope = " // trim(strVal) // char(10)
  end if
  

  
  return

988  call WriteFatalError( "Could not parse B. Did you specify a comma separated list with one value per eigenmode?")

989  call WriteFatalError( "Could not parse alpha. Did you specify a comma separated list with one value per eigenmode?")

990  call WriteFatalError( "Could not parse natural frequencies. Did you specify a comma separated list with one natural frequency per eigenmode?")

991  call WriteFatalError( "Could not parse stiffnesses. Did you specify a comma separated list with one stiffness per eigenmode?")

992   call WriteFatalError( "Could not parse Quality factors.  Did you specify a comma separated list with one Quality factor per eigenmode?")

993   call WriteFatalError( "Could not parse Chi values.  Did you specify a comma separated list with one Chi value for eigenmodes 2 - n?")
end subroutine readCantModalProps

!only used in forceViewer.f90
subroutine readForceViewerParameters( Z0, Zf, plotpnts, vel_model, tf, mov_avg_filt_len)
  integer, intent(out) :: vel_model, mov_avg_filt_len
  integer :: readGenericInteger, rp_units_convert_dbl, status
  character*100 :: strVal
  real*8 :: moving_avg_fraction 
  real*8, intent(out) :: tf, Z0, Zf
  integer*8, intent(out) :: plotpnts

  
  status = 0

  strVal = rp_lib_get_wrap(driver, INPUT_PREFIX // "(op).number(Z0).current" )
  !           InputEcho = trim(InputEcho) // "Z0 " // StrVal // char(10)
  status = status +  rp_units_convert_dbl(strVal," ",Z0)
  Z0 = Z0/1d9
  
  strVal = rp_lib_get_wrap(driver, INPUT_PREFIX // "(op).number(Zf).current" )
  status = status +  rp_units_convert_dbl(strVal," ",Zf)
  Zf = Zf/1d9

  plotpnts = readGenericInteger( driver, INPUT_PREFIX // "(op).integer(plotpnts).current")       

  vel_model = readGenericInteger( driver, INPUT_PREFIX // "(op).choice(vel_model).current")  

  if (vel_model > 1) then
     strVal = rp_lib_get_wrap(driver, INPUT_PREFIX // "(op).number(t).current" )
     status = status +  rp_units_convert_dbl(strVal," ",tf)
  end if


  if (vel_model == 4 ) then
     strVal = rp_lib_get_wrap(driver, INPUT_PREFIX // "(op).number(moving_avg_len).current" )
     status = status +  rp_units_convert_dbl(strVal," ",moving_avg_fraction)
     mov_avg_filt_len =  (moving_avg_fraction/100d0) * plotpnts
  end if
  
  if (status > 0) call WriteFatalError( "Could not read operating conditions. Please check values" )

end subroutine readForceViewerParameters

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!111
! setup plots


subroutine setup_hydro_plots()
  call  setupgenericplot( 'addM', 'added mass', 'hydrodynamic', 'freq', 'khz', 'm',    'g')
  call  setupgenericplot( 'addC', 'added viso', 'hydrodynamic', 'freq', 'khz', 'visc', '')
  call  setupgenericplot( 'dbgomega', 'omega', 'hydrodynamic',  'freq', 'khz', 'freq', '')
  call  setupgenericplot( 'dbgdamping', 'damping', 'hydrodynamic','freq','khz', 'nd', '')
end subroutine setup_hydro_plots

subroutine SetupDebugAMScan_Transient()
  call SetupGenericPlot( "Z0", "Transient Z0", "Transient", "t", "s", "Z", "nm")
  call SetupGenericPlot( "dbgAmp", "Transient Amp", "Transient", "t", "s", "Amp", "ND") 
end subroutine SetupDebugAMScan_Transient

subroutine SetupDebugForceModScan_Transient()
  call SetupGenericPlot( "Z0",      "Transient Z0", "Transient", "t", "s", "Z", "nm")
  call SetupGenericPlot( "dbgMeanDefl", "Transient Mean Defl", "Transient", "t", "s", "defl", "nm") 
  call SetupGenericPlot( "dbgDefl", "Transient Defl", "Transient", "t", "s", "defl", "nm") 
end subroutine SetupDebugForceModScan_Transient


subroutine SetupDebugContactScan_Transient()
  call SetupGenericPlot( "Z0",      "Transient Z0", "Transient", "t", "s", "Z", "nm")
  call SetupGenericPlot( "dbgDefl", "Transient Defl", "Transient", "t", "s", "defl", "nm") 
end subroutine SetupDebugContactScan_Transient

subroutine SetupDebugPeakForceScan_Transient()
  call SetupGenericPlot( "Z0",      "Transient Z0", "Transient Overview", "t", "s", "Z", "nm")
  call SetupGenericPlot( "dbgPeak", "Transient Force", "Transient Overview", "t", "s", "Force", "nN") 

  call SetupGenericPlot( "Z0_full",      "Transient Z0",    "Transient Full", "t", "s", "Z", "nm")
  call SetupGenericPlot( "dbgPeak_full", "Transient Force", "Transient Full", "t", "s", "Force", "nN") 
  call SetupGenericPlot( "d_full",       "Transient d",     "Transient Full", "t", "s", "Z", "nm")
  call SetupGenericPlot( "q1_full",      "Transient q1",    "Transient Full", "t", "s", "Z", "nm")
  call SetupGenericPlot( "dpr_full",     "Transient dpr",   "Transient Full", "t", "s", "Z", "um/s")
  
end subroutine SetupDebugPeakForceScan_Transient



subroutine SetupDebugFM_Transient()
!  call SetupGenericPlot( "thetaTrans",  "theta",  "Transient Lockin", "t", "(s)", "ND", "ND")
  call SetupGenericPlot( "lckR",  "Lock in r",  "Transient Lockin", "t", "(s)", "r", "ND")
  call SetupGenericPlot( "lckTh", "Lock in th", "Transient Lockin", "t", "(s)", "th", "deg")

  call SetupGenericPlot( "Z0",     "Z0",     "Transient", "t", "(s)", "Z", "nm")
  call SetupGenericPlot( "Zerr",   "Zerr",   "Transient", "t", "(s)", "Zerr", "")
  call SetupGenericPlot( "omegad", "Omegad", "Transient", "t", "(s)", "freq", "(kHz)")           
  call SetupGenericPlot( "DrvAmpTrans",   "Dr Amp",          "Transient", "t", "(s)", "amp", "ND")            
  call SetupGenericPlot( "PhaseErrTrans",    "Phase Err Trans",     "Transient", "t", "(s)", "deg", "ND")            
!  call SetupGenericPlot( "PhaseErrIntTrans", "Phase Err Int Trans", "Transient", "t", "(s)", "deg*t", "ND")            
  call SetupGenericPlot( "AmpErrTrans",    "Amp Err Trans",     "Transient", "t", "(s)", "ND", "ND")            
!  call SetupGenericPlot( "AmpErrIntTrans", "Amp Err Int Trans", "Transient", "t", "(s)", "ND*t", "ND")            

end subroutine SetupDebugFM_Transient

subroutine SetupScanPlots(fexcite,z_feedback_choice)
  character*50 :: strVal
  integer, intent(in) :: fexcite, z_feedback_choice

  strVal = "Measured topography"
  call SetupGenericPlot( "Topog", strVal, strVal, "X distance", "nm", "Height", "nm")

  call SetupGenericPlot( "Hsample",  "Sample height vs X distance", "Measured topography", "X distance", "nm", "Height", "nm")

  strVal = "Measurement error vs X distance"
  call SetupGenericPlot( "ErrorZ", strVal, strVal, "X distance", "nm",  "Measurement error", "nm")

  if ((fexcite .ne. NO_EXC) .and.(z_feedback_choice .ne. MEAN_DEFL_Z) .and. (z_feedback_choice .ne. MAX_FORCE)) then
     !fixme, AmpDist needs a better name.       
     strVal = "Amplitude error vs X distance"
     call SetupGenericPlot( "AmpDist", strVal, strVal, "X distance", "nm", "Amplitude error", "nm")
  elseif (z_feedback_choice .ne. MAX_FORCE) then
     strVal = "Deflection error vs X distance"
     call SetupGenericPlot( "DeflErr", strVal, strVal, "X distance", "nm", "Deflection error", "nm")
  end if
  
end subroutine SetupScanPlots


subroutine SetupPeakForcePlots()
  character*50 :: strVal1, strVal2

  strVal1 = "Peak Force vs X distance"
  call SetupGenericPlot( "PeakDistR", strVal1, strVal1, "X distance", "nm", "Peak Force", "nN")

  strVal2 = "Min Force vs X distance"
  call SetupGenericPlot( "PeakDistA", strVal2, strVal1, "X distance", "nm", "Peak Force", "nN")

end subroutine SetupPeakForcePlots

!might need to put scatter plot option back in to avoid confusing people
!formerly called "jump"
subroutine SetupFzCurvesNew( numModes, want_FzEig, wantAttard)
  integer, intent(in) ::  numModes
  logical, intent(in) :: want_FzEig, wantAttard
  integer :: i
  character*100 strVal
  character*2 tmpStr

  strVal = "Observed cantilever deflection vs Z distance"
  call SetupGenericPlot( "DeflZ_A", strVal, strVal, "Z distance", "nm", "Observed Deflection", "nm")
  call SetupGenericPlot( "DeflZ_R", strVal, strVal, "Z distance", "nm", "Observed Deflection", "nm")
  
  strVal = "Z distance vs time"
  call SetupGenericPlot( "ZT", strVal, strVal, "Time", "s", "Z distance", "nm")

!  strVal = "Z dot vs time"
!  call SetupGenericPlot( "ZdotT", strVal, strVal, "Time", "s", "Z velocity", "um/s")

  strVal =  "Tip-sample interaction force vs gap"
  call SetupGenericPlot( "FtsGap", strVal, strVal,"tip-sample gap", "nm", "Force", "nN")   

  strVal =  "Tip-sample interaction force vs Z distance"
  call SetupGenericPlot( "FtsZ_A", strVal, strVal,"Z distance", "nm", "Force", "nN")   
  call SetupGenericPlot( "FtsZ_R", strVal, strVal,"Z distance", "nm", "Force", "nN")   

  strVal = "Tip-sample gap vs Z distance"
  call SetupGenericPlot( "gapZ_A", strVal, strVal, "Z distance", "nm", "tip-sample gap", "nm")
  call SetupGenericPlot( "gapZ_R", strVal, strVal, "Z distance", "nm", "tip-sample gap", "nm")

  strVal = "Tip-sample gap vs time"
  call SetupGenericPlot( "gapT", strVal, strVal, "time", "s",  "Tip sample gap", "nm")

  strVal = "Derivative of Tip-sample gap vs time"
  call SetupGenericPlot( "dprT", strVal, strVal, "time", "s",  "Tip sample gap derivative", "um/s")
  
  strVal =  "Tip-sample interaction force vs time"
  call SetupGenericPlot( "FtsT", strVal, strVal, "time", "s", "Force", "nN")   

  strVal =  "Indentation vs Z distance"
  call SetupGenericPlot( "Indent_A", strVal, strVal,"Z distance", "nm", "Indentation", "nm", .true.)
  call SetupGenericPlot( "Indent_R", strVal, strVal,"Z distance", "nm", "Indentation", "nm", .true.)

  if ( want_FzEig) then
     do i = 1, numModes
        write(tmpStr, '(I1)') i
        strVal =  "Eigenmode responses vs time"
        call SetupGenericPlot( "Eig" // trim(tmpStr), strVal, strVal, "time", "s", "Modal Deflection", "nm")   
        call SetupGenericPlot( "EigVel" // trim(tmpStr), "Velocity " // strVal, "Velocity " // strVal, "time", "s", "Modal Velocity", "um / s")   
     end do
  end if

  if (wantAttard) then
     strVal =  "Attard Surface Position vs time"
     call SetupGenericPlot( "AttardSurfPos", strVal, "Tip-sample gap vs time", "time", "s", "Displacement", "nm")

     call start_movie(driver, "                ", "1", "Time (s)"  )     
  end if
     
end subroutine SetupFzCurvesNew


subroutine SetupPoincPlots(numpoinc,Npoinc_X,Npoinc_Y,Want_Strob,Want_Impact)
  use Poincare, only: maxpoinc
  integer, intent(in) :: numpoinc
  integer, intent(in) :: Npoinc_X(maxpoinc), Npoinc_Y(maxpoinc)
  logical, intent(in) :: Want_Strob, Want_Impact
  character*100 strVal
  character*1 Pnum
  character*3 Xvar
  character*8 Yvar, yunits, xunits
  integer :: ii

  if (( Want_Strob ).or.( Want_Impact )) then
     
     do ii = 1,numpoinc
        write (Pnum, '(I1)') ii
        if (Npoinc_X(ii) == -2 ) then
           Xvar = "Z"
           xunits = "nm"
        elseif (Npoinc_X(ii) == -1 ) then
           Xvar = "t"
           xunits = "s"
        elseif (Npoinc_X(ii) == 0 ) then
           Xvar = "d"
           xunits = "nm"
        elseif (Npoinc_X(ii) < 19 ) then
           ! 1-> 1, 3->2, 5->3... ii -> (ii-1)/2 ... 
           write( Xvar, '(A,I1)') 'q_', ((Npoinc_X(ii)+1)/2)
           xunits = "nm"
        else
           Xvar = "?"
           call assert(.false., 'cant handle more than 9 modes for poincare')
        end if
        
        if (Npoinc_Y(ii) == 0 ) then
!           Yvar = "d'"
           Yvar = "d/dt d"
           yunits = "um/s"
        else if (Npoinc_Y(ii) < 18 ) then
           if ( 2*(Npoinc_Y(ii)/2) == Npoinc_Y(ii)) then
              !even means velo
              write( Yvar, '(A,I1)') 'd/dt q_', (Npoinc_Y(ii)/2)          
              yunits = "um/s"
           else
              !odd means disp
              write( Yvar, '(A,I1)') 'q_', ((Npoinc_X(ii)+1)/2)
              yunits = "nm"
           end if
        else
           Yvar = "?"
           call assert(.false., 'cant handle more than 9 modes for poincare')
        end if
        
        if ( Want_Strob ) then
           strVal = "Stroboscopic Poincare (" // Pnum // ")"
           call SetupGenericPlot( "Poinc_Strob" // Pnum, strVal, strVal, Xvar, xunits, Yvar,yunits, .true.)
        end if
        
        if ( Want_Impact ) then
           strVal = "Impact Poincare (" // Pnum // ")"
           call SetupGenericPlot( "Poinc_Impact" // Pnum, strVal, strVal, Xvar, xunits, Yvar, yunits, .true.)
        end if
        
     end do
     
  end if

end subroutine SetupPoincPlots

subroutine SetupHHPlots(numHH,NHH,exc_choice, xlabel, xunits, Want_Fourier, useLockin, scatter)
  integer, intent(in) ::  numHH, exc_choice
  real*8, intent(in) :: NHH(numHH)
  logical, intent(in) :: Want_Fourier, useLockin, scatter
  character*100, intent(in) :: xlabel, xunits

  integer :: ii
  character*100 strVal, ylabel
  character*2 Snum
  character*5 HarmLabel

  
     do ii = 1,numHH       

        if ( ii >= 10) then  !got to be an easier way to do this
           write (Snum, '(I2)') ii
        else
           write (Snum, '(I1)') ii
        end if

        write (HarmLabel, '(F5.2)') NHH(ii)
        strVal = "Higher harmonic (" // HarmLabel // ") amplitude"
        
        if (exc_choice == BIMODAL) then
           ylabel = "A_" // HarmLabel // ",1"
        else
           ylabel = "A_" // HarmLabel
        end if

        if (NHH(ii) == 0) then
           strVal = "Mean Deflection (Zeroth Harmonic)"
           ylabel = "Mean deflection"
        end if

        if (Want_Fourier) call SetupGenericPlot( "AmpHH" // trim(Snum), strVal // " (Fourier integral)" , strVal, xlabel, xunits, trim(ylabel), "nm", scatter )
        if (useLockin) call SetupGenericPlot( "LockinHH" // Trim(Snum) // "R", strVal // " (Lockin)", strVal, xlabel, xunits, ylabel, "nm", scatter )

        if (numHH>0) then !dont want to output phase for zeroth harmonic
           strVal = "Higher harmonic (" // HarmLabel // ") phase"
        
           if (exc_choice == BIMODAL) then
              ylabel = "Phase_" // HarmLabel // ",1"
           else
              ylabel = "Phase_" // HarmLabel
           end if
          
           if (Want_Fourier) call SetupGenericPlot( "PhaseHH" // Trim(Snum), strVal // " (Fourier integral)" , strVal, xlabel, xunits, ylabel, "deg", scatter )
           if (useLockin) call SetupGenericPlot( "LockinHH" // Trim(Snum) // "Th", strVal // " (Lockin)" , strVal, xlabel, xunits, ylabel, "deg", scatter )
        end if
     end do

end subroutine SetupHHPlots

subroutine SetupMainOutputPlots(xchoice,exc_choice,Want_A2,Want_P2,Want_AZ,Want_P1,Want_MF,Want_PF,Want_ED,Want_EP,Want_Indent,Want_CT, Want_Fourier, Want_NumImpact, Want_E_Anc, Want_ForceFourier,  want_Virial, useLockIn, numModes, operating_mode, fexcite, want_ev, want_RMS, xlabel, xunits, scatter, AmpRedAvg)

  implicit none
  logical, intent(in) :: Want_AZ, Want_P1, Want_MF, Want_PF, Want_ED, Want_EP, Want_Indent, Want_CT,& 
                         Want_NumImpact , Want_A2, Want_P2, Want_Fourier, useLockIn, & 
                         Want_E_Anc, Want_ForceFourier, want_ev, want_virial, want_RMS, AmpRedAvg
  integer, intent(in) :: xchoice,  exc_choice, fexcite
  integer, intent(in) :: numModes, operating_mode
  character*100, intent(out) :: xlabel, xunits
  logical, intent(out) ::  scatter
  character*100 strVal, ylabel, yunits
  character*1 String1
  integer ii


  !do the first harmonic amplitude first.  it might use a different axis than all the rest

  if ((operating_mode == APPR_RET ) .or. (operating_mode == APPROACH )  .or. (operating_mode == APPROACH_STEP ) .or. (operating_mode == FIXED)) then
     if (xchoice == X_AMP_BIMODAL) then
        xlabel = "A_1,2/A_0,2"
        xunits = ""
     elseif (xchoice == X_MINGAP) then
        xlabel = "Min gap"
        xunits = "nm"        
     else
        xlabel = "Z distance"
        xunits = "nm"
     end if
  else if (operating_mode == FREQSWEEP) then
     xlabel = "Frequency"
     xunits = "kHz"
  else if ( operating_mode == SCAN) then
     xlabel = "X distance"
     xunits = "nm"        
  end if
  
  if (exc_choice == BIMODAL) then
     ylabel = "A_1,1"
  else
     ylabel = "A_1"
  end if

  if (AmpRedAvg) then
     strVal = "First harmonic amplitude"        
     call SetupGenericPlot( "Amp",       strVal // " (Target)",          strVal, xlabel, xunits, ylabel, "nm")
     call SetupGenericPlot( "Amp_iter",   strVal // " (Final Iteration)", strVal, xlabel, xunits, ylabel, "nm")
     
  elseif (Want_AZ) then
     strVal = "First harmonic amplitude"        

     if (want_Fourier) call SetupGenericPlot( "Amp", strVal // " (Fourier integral)", strVal, xlabel, xunits, ylabel, "nm")

     if (useLockIn) call SetupGenericPlot( "LockinR", "1st Harmonic Amplitude (Lock-in)", strVal, xlabel, xunits, ylabel, "nm")
     
  end if

  if (want_RMS) call SetupGenericPlot( "RMSAmp", "RMS Amplitude", "RMS Amplitude", xlabel, xunits, "A_RMS", "nm")

  !then the 2nd exc frequency if applicable.  it might use different units as well.
  if (exc_choice == BIMODAL) then
     if (xchoice == X_AMP) then
        xlabel = "A_1,1/A_0,1"
        xunits = ""
     elseif (xchoice == X_MINGAP) then
        xlabel = "Min gap"
        xunits = "nm"        
     else
        xlabel = "Z distance"
        xunits = "nm"
     end if

     strVal = "Second frequency amplitude"
     if (AmpRedAvg) then
        call SetupGenericPlot( "Amp2",      strVal // " (Final Iteration Input)",  strVal, xlabel, xunits, "A_1,2", "nm", scatter)        
        call SetupGenericPlot( "Amp2_iter", strVal // " (Final Iteration Output)", strVal, xlabel, xunits, "A_1,2", "nm", scatter)
     else     
        if (want_Fourier .and. Want_A2) then
           call SetupGenericPlot( "Amp2", strVal // " (Fourier Integral)", strVal, xlabel, xunits, "A_1,2", "nm", scatter)        
        end if
        if (Want_A2) then     
           call SetupGenericPlot( "LockinAmp2", strVal // " (Lockin)", strVal, xlabel, xunits, "A_1,2", "nm", scatter)
        end if
     end if

  end if
     
    !everything else will have same labels
    scatter = .false.
    if ( isOpModeApp(operating_mode)) then
       if (xchoice == X_AMP) then
!          scatter = .true. !per john's request
          if ( exc_choice == BIMODAL) then
             xlabel = "A_1,1/A_0,1"
             xunits = ""
          else
             xlabel = "A_1/A_0"
             xunits = ""
          end if
       elseif (xchoice == X_AMP_BIMODAL) then
!          scatter = .true. !per john's request
          xlabel = "A_1,2/A_0,2"
          xunits = ""
       elseif (xchoice == X_MINGAP) then
          xlabel = "Min gap"
          xunits = "nm"
       else
        xlabel = "Z distance"
        xunits = "nm"
       end if
    end if 
    
  if (Want_P1) then
     if (exc_choice == BIMODAL) then 
        strVal = "Phase_1,1"
     else
        strVal = "Phase_1"
     end if

     if (want_Fourier) then
        call SetupGenericPlot( "Phase", "First harmonic phase (Fourier integral)", "First harmonic phase", xlabel, xunits, strVal, "deg", scatter)
     end if

     if (useLockIn) call SetupGenericPlot( "LockinPhase", "First harmonic phase (Lock-in)", "First harmonic phase", xlabel, xunits, strVal, "deg", scatter)

  end if
    
    if (exc_choice == BIMODAL) then
	       
       if (want_Fourier .and. Want_P2) then
          strVal = "Second frequency phase"
          call SetupGenericPlot( "Phase2",  strVal // " (Fourier Integral)",strVal, xlabel, xunits, "Phase_1,2", "deg", scatter)
       end if
       if (Want_P2) then
          strVal = "Second frequency phase"
          call SetupGenericPlot( "LockinPhase2", strVal // " (Lockin)", strVal, xlabel, xunits, "Phase_1,2", "deg", scatter)
       end if

    end if

    if (Want_ForceFourier) then
       strVal = "1st harmonic of tip-sample force"
       call SetupGenericPlot( "FtsHarm", strVal, strVal, xlabel, xunits, "Force", "nN", scatter)
    end if

    if (Want_MF) then
       strVal = "Mean interaction forces"
       call SetupGenericPlot( "MeanForce", strVal, strVal, xlabel, xunits, "Mean force", "nN", scatter)
    end if
    
    if (Want_PF) then
       call SetupGenericPlot( "FPeakRep", "Peak repulsive force", "Peak interaction forces", xlabel, xunits, "Peak force", "nN", scatter)
       call SetupGenericPlot( "FPeakAtt", "Peak attractive force", "Peak interaction forces", xlabel, xunits, "Peak force", "nN", scatter)
    end if
     
    if (want_ev) then
       yunits = "eV / drive cycle"
    else
       yunits = "pW"
    end if
    
    if (Want_ED) then
       strVal = "Tip-sample energy dissipation"
                
       if (AmpRedAvg .and.  (exc_choice == BIMODAL)) then
                       !          (name, label,                 group, xlabel, xunits, ylabel,                    yunits,                       scatter
          call SetupGenericPlot( "Pts1", "Dissipation mode 1", strVal, xlabel, xunits, "Tip-sample dissipation", 'eV / first freq drive cycle', scatter)
          call SetupGenericPlot( "Pts2", "Dissipation mode 2", strVal, xlabel, xunits, "Tip-sample dissipation", 'eV / first freq drive cycle', scatter)
          call SetupGenericPlot( "Pts",  "Dissipation total",  strVal, xlabel, xunits, "Tip-sample dissipation", 'eV / first freq drive cycle', scatter)
       else
          call SetupGenericPlot( "Pts", strVal, strVal, xlabel, xunits, "Tip-sample dissipation", yunits, scatter)
       end if
    end if

    if (Want_Virial) then
       strVal = "Tip-sample virial"
       
       if (AmpRedAvg  .and.  (exc_choice == BIMODAL) ) then !amp red tool uses this
          call SetupGenericPlot( "virial1", "mode 1", strVal, xlabel, xunits, "Tip-sample virial (mode 1)", yunits, scatter)
          call SetupGenericPlot( "virial2", "mode 2", strVal, xlabel, xunits, "Tip-sample virial (mode 2)", yunits, scatter)
          call SetupGenericPlot( "virial",  "total", strVal, xlabel, xunits, "Tip-sample virial", yunits, scatter)
       else
          call SetupGenericPlot( "virial", strVal, strVal, xlabel, xunits, "Tip-sample virial", yunits, scatter)
       end if       
    end if


    if (Want_E_Anc) then
       strVal = "Ancyzkowski formula (for comparison only - may be inaccurate for liquids)"
       call SetupGenericPlot( "E_anc", strVal, strVal, xlabel, xunits, "Dissipation", yunits, scatter) 
    end if

    if (Want_EP) then
       if (isAcoustic(fexcite)) then
          call SetupGenericPlot( "Ebs", "Base Sample Work", "Base Sample Work", xlabel, xunits, "Base Sample Work", yunits, scatter)
       end if

        do ii=1,numModes

	  write (String1, '(I1)') ii
          strVal = "Energy propagation mode " // String1
          call SetupGenericPlot( "Eprop_" // String1, strVal, "Energy Propagation", xlabel, xunits, "Energy Propagation", yunits, scatter)
          strVal = "Fluid damping energy mode " // String1
          call SetupGenericPlot( "Edmp_" // String1, strVal, "Fluid damping energy", xlabel, xunits, "Fluid Damping Energy", yunits, scatter)
          strVal = "Drive input energy mode " // String1
          call SetupGenericPlot( "Edrv_" // String1, strVal, "Drive Input energy", xlabel, xunits, "Drive Input Energy", yunits, scatter)

	end do
      end if
            
      if (Want_Indent) then
        strVal = "Indentation"
        call SetupGenericPlot( "Indent", strVal, strVal, xlabel, xunits, "Indentation", "nm", scatter) 
      end if
      
      if (Want_CT) then
        strVal = "Contact time / drive cycle"
        call SetupGenericPlot( "tcontact", strVal, strVal, xlabel, xunits, strVal, "us", scatter) 
      end if

      if (Want_NumImpact) then
        strVal = "Number impacts per drive cycle"
        call SetupGenericPlot( "NumImpact" , strVal, strVal, xlabel, xunits, strVal, "ND", scatter) 
     end if

     
end subroutine SetupMainOutputPlots




subroutine SetupTimeHistoryPlots(value,xvalue, operating_mode, wantHist_byA, hist_cycles, numHist, numincycle, exc_choice, want_voltage, want_TH_obs_defl, want_th_force_gap, want_surf_relax, do_fft, want_contact_area, Want_FZ, want_acceleration, want_modal, want_PerfMetrics)  
  use data1, only : numModes
  implicit none
  
  logical, intent(in) :: wantHist_byA, want_voltage, want_TH_obs_defl, want_th_force_gap, want_surf_relax, want_contact_area, do_fft, Want_FZ, want_acceleration, want_modal, want_PerfMetrics
  integer, intent(in) ::  operating_mode, hist_cycles,  numHist, numincycle, exc_choice
  real*8, intent(in) :: value, xvalue
  character*16 SA
  character*100 strVal
  character*1 ST, SE
  integer i, status, numplots
  
  !remember this information so that we can output time histories later, even if called from WriteFatalError which does not
  !have that information handy
  th_set%want_voltage = want_voltage
  th_set%want_TH_obs_defl = want_TH_obs_defl
  th_set%want_th_force_gap = want_th_force_gap
  th_set%want_surf_relax = want_surf_relax
  th_set%want_contact_area = want_contact_area
  th_set%want_self_exc = (exc_choice == SELFEXC)
  th_set%want_FZ = want_FZ
  th_set%want_acceleration = want_acceleration
  th_set%want_modal = want_modal
  th_set%want_PerfMetrics = want_PerfMetrics
  th_set%want_sampletopo = (operating_mode == SCAN)
  
  buffer_len4 = hist_cycles * numincycle

  numplots = 2 + numModes

  buffer_ndx = buffer_ndx + 1
  
  if (((operating_mode == APPR_RET) .or.  (operating_mode == APPROACH)  .or.  (operating_mode == APPROACH_STEP)) .and. (.not. wantHist_byA)) then
     write (SA, '(A,F7.2,A)') "(Z=", value, ")"
  elseif (operating_mode == FREQSWEEP) then
     write (SA, '(A,F7.2,A)') "(F=", value, "kHz)"
  elseif (operating_mode == SCAN) then
     write (SA, '(A,F7.2,A)') "(X=", xvalue, ")" 
  elseif (operating_mode == FIXED) then
     SA = ""
  else
     write (SA, '(A,F7.3,A)') "(A/Ao=", value, ")"
  end if
  
  write (ST, '(I1)') buffer_ndx

  if (want_surf_relax) then
     if (do_fft) then
        call SetupGenericPlot(  "surf" // ST, "Surface Coord (center) " // SA, "Deflection and Force History " // SA, "Frequency", "harmonics", "Deflection", "nm")
     else
        call SetupGenericPlot(  "surf" // ST, "Surface Coord (center) " // SA, "Deflection and Force History " // SA, "t/T", "periods", "Deflection", "nm")
     end if
     numplots = numplots+1 

     call start_movie(driver, SA, ST, "Periods (t/T)")
  end if
  
  if (want_TH_obs_defl) then
     if (do_fft) then
        call SetupGenericPlot(  "defl" // ST, "Observed Deflection " // SA, "Deflection and Force History " // SA, "Frequency", "harmonics", "Observed Deflection", "nm")
     else
        call SetupGenericPlot(  "defl" // ST, "Observed Deflection " // SA, "Deflection and Force History " // SA, "t/T", "periods", "Observed Deflection", "nm")
     end if
     numplots = numplots+1 
  end if

  if (want_contact_area) then
     call SetupGenericPlot(  "cr" // ST, "Contact Radius " // SA, "Deflection and Force History " // SA, "t/T", "periods", "Contact Radius", "nm")
     numplots = numplots+1 
  end if

  
  if (want_th_force_gap) then
     call SetupGenericPlot(  "Fts" // ST, "Tip-sample force " // SA, "Deflection and Force History " // SA, "gap", "(nm)", "Tip-Sample Force", "nN")
  else
     numplots = numplots+1
     if (do_fft) then
        call SetupGenericPlot(  "gap" // ST, "Tip-sample gap (to undeformed surf) "   // SA, "Deflection and Force History " // SA, "Frequency", "harmonics", "Tip Sample Gap", "nm")  
        call SetupGenericPlot(  "Fts" // ST, "Tip-sample force " // SA, "Deflection and Force History " // SA, "Frequency", "harmonics", "Tip-Sample Force", "nN")
     else
        call SetupGenericPlot(  "gap" // ST, "Tip-sample gap (to undeformed surf) " // SA, "Deflection and Force History " // SA, "t/T", "periods", "Tip Sample Gap", "nm")  
        call SetupGenericPlot(  "Fts" // ST, "Tip-sample force " // SA, "Deflection and Force History " // SA, "t/T", "periods", "Tip-Sample Force", "nN")
     end if
  end if
   
  call SetupGenericPlot(  "Pspace" // ST, "Phase space " // SA, "Phase space " // SA, "d", "nm", "d/dt d", "um/s")

  if (want_modal) then
     do i =1,numModes
        write (SE, '(I1)') i                    
        call SetupGenericPlot(  "eig" // SE // ST, "Mode" // SE // " " // SA // ")" ,  "Modal response  " // SA, "t/T", "periods", "Response", "nm")
     end do
  end if
  
  if (exc_choice == SELFEXC) then
     call SetupGenericPlot(  "drvforce" // ST, "Modal Drive Force (mode 1) " // SA, "Deflection and Force History " // SA, "t/T", "periods", "Force", "nN")
     numplots = numplots+1
  end if

  if (want_voltage) then
     call SetupGenericPlot(  "voltage" // ST, "Voltage " // SA, "Deflection and Force History " // SA, "t/T", "periods", "Volts", "")
     numplots = numplots+1
  end if

  if (want_acceleration) then
     call SetupGenericPlot(  "TipAcceleration" // ST, "Tip Acceleration " // SA, "Deflection and Force History " // SA, "t/T", "periods", "mm/s^2", "")
     numplots = numplots+1
  end if


  if (want_FZ) then
     call SetupGenericPlot(  "FvsZ" // ST, "Force versus Z" // SA, "Force vs Z " // SA, "Z-distance", "(nm)", "Tip-Sample Force", "nN")
     numplots = numplots+1
  end if


  !this is really only for debugging.  average user will not care about this.  
  if (Want_PerfMetrics) then
     strVal = "Number calls to RES1 per output point"
     call SetupGenericPlot( "NumRES1" // ST, strVal // SA , strVal // SA, "t/T", "periods", "Count" , "")
     
     numplots = numplots+1
  end if


  if (th_set%want_sampletopo) then
     strVal = "Sample Topography"
     call SetupGenericPlot( "Topo" // ST, strVal // SA , "Deflection and Force History " // SA, "t/T", "periods", "Height" , "nm")
     
     numplots = numplots+1
  end if

     
  if (.not. allocated(outputbuffer) ) then
     allocate(outputbuffer( numplots, numHist, 2, buffer_len4), STAT=status)
     if (status > 0) then
        call WriteFatalError( "could not allocate output buffer");
     end if
  end if

     
end subroutine SetupTimeHistoryPlots

!just do these the slow way.  they are debug plots so we wont use them often
subroutine setup_debug_viscoelastic()
  call  setupgenericplot( 'visc_t1', 't1', 'Viscoelasticity Debug', 'time', 't/T', 't1', 'ndx')
end subroutine setup_debug_viscoelastic


subroutine SetupGenericPlot(name, label, group, xlabel, xunits, ylabel, yunits, scatter )
  character(len=*), intent(in) :: name, label, group, xlabel, xunits, ylabel, yunits
  logical, intent(in), optional :: scatter
         
  call rp_lib_put_str(driver, "output.curve(" //name// ").about.label", label, 0) 
  call rp_lib_put_str(driver, "output.curve(" //name// ").about.group", group, 0) 
  call rp_lib_put_str(driver, "output.curve(" //name// ").xaxis.label", xlabel, 0) 
  call rp_lib_put_str(driver, "output.curve(" //name// ").xaxis.units", xunits, 0) 
  call rp_lib_put_str(driver, "output.curve(" //name// ").yaxis.label", ylabel,0)
  call rp_lib_put_str(driver, "output.curve(" //name// ").yaxis.units", yunits, 0) 
  if (present(scatter)) then
     if (scatter) then
        call rp_lib_put_str(driver, "output.curve(" //name//").about.type", "scatter", 0)
     end if
  end if
end subroutine SetupGenericPlot


subroutine SetupFMPlots( xchoice, fm_direct_control, want_NormFreqShift)
  integer, intent(in) ::  xchoice
  logical, intent(in) :: fm_direct_control, want_NormFreqShift
  character*7 xlabel, ylabel

  if (xchoice == X_MINGAP) then
     xlabel = "Min gap"
  else
     xlabel = "Z dist"
  end  if

  if (want_NormFreqShift) then
     call SetupGenericPlot( "DriveFreq", "Freq Shift", "FM Feedback (delta f and drive amp vs X/Z)", xlabel, "nm", "Normalized Freq Shift", "ND")
  else
     call SetupGenericPlot( "DriveFreq", "Freq Shift", "FM Feedback (delta f and drive amp vs X/Z)", xlabel, "nm", "Freq Shift", "kHz")
  end if

  call SetupGenericPlot( "DriveAmp", "Drive Amp",   "FM Feedback (delta f and drive amp vs X/Z)", xlabel, "nm", "Drive Amp", "ND")

  call SetupGenericPlot( "AmpVsFreq", "Dr Amp vs Freq Shift", "FM Feedback (delta f vs drive amp)", "Freq Shift", "kHz", "Normalized Drive Amp", "ND")

  if (fm_direct_control) then
     ylabel = "deg"
  else
     ylabel = "ND"
  end if

  call SetupGenericPlot( "PhaseErr", "Phase Error",  "Phase Error", xlabel, "nm", "Phase Error", ylabel)
end subroutine SetupFMPlots


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! output routines

subroutine output_debug_viscoelastic(t, t1)
  real*8 :: t, t1
  character(len=*), parameter :: fmtStr = '(ES16.7E3,ES16.7E3,A)'
  character*150 :: tmpStr

  write(tmpStr,fmtStr) t, t1, char(10)
  call rp_lib_put_str(driver, "output.curve(visc_t1).component.xy", tmpStr, 1)
end subroutine output_debug_viscoelastic


subroutine Output_hydro_Plots( xaxis, added_m, added_c, omega, damping)
!subroutine OutputDebugPlots( xaxis, relaxed_d)
  real*8 :: xaxis, added_m, added_c, omega, damping
!  real*8, intent(in) :: xaxis, relaxed_d
  character(len=*), parameter :: fmtStr = '(ES16.7E3,ES16.7E3,A)'
  character*150 :: tmpStr

 write(tmpStr,fmtStr) xaxis,  added_m, char(10)
 call rp_lib_put_str(driver, "output.curve(addM).component.xy", tmpStr, 1)
 write(tmpStr,fmtStr) xaxis,  added_c, char(10)
 call rp_lib_put_str(driver, "output.curve(addC).component.xy", tmpStr, 1)
 write(tmpStr,fmtStr) xaxis,  omega, char(10)
 call rp_lib_put_str(driver, "output.curve(dbgomega).component.xy", tmpStr, 1)
 write(tmpStr,fmtStr) xaxis,  damping, char(10)
 call rp_lib_put_str(driver, "output.curve(dbgdamping).component.xy", tmpStr, 1)
end subroutine Output_hydro_Plots

subroutine OutputInputEcho( InputEcho1)
  character*3000, intent(in) :: InputEcho1

  call rp_lib_put_str(driver, "output.string(echo).about.label", "Echo of input parameters",0)
  call rp_lib_put_str(driver, "output.string(echo).current", InputEcho1,0)
end subroutine OutputInputEcho

!a separate tab to output various internal numbers that we might care about
subroutine OutputMiscParameters( numModes, F, phid, Keq, pad, numincycle, aDMT, Chi, omegai, B, fm_initial_phase, Z0, Abase, Y, beta, mu, Estar, damping, nstep, fm_amp_kp, fm_amp_ki, fm_amp_kd, fm_pll_kp, fm_pll_ki, fm_pll_kd, alpha, output_point_rate, delt, numofcycle, z_feedback_choice, minimum_relaxation, omegad_t_1)
!  use data1, only: NEQ
  use Nondimensionalization

  character(len = 210) :: s
  character(len=2500) :: Str
  real*8, intent(in) :: Keq(maxModes), aDMT, Chi(maxModes), omegai(maxModes), minimum_relaxation
  real*8, intent(in) :: F(maxExc, maxModes), phid(maxExc, maxModes), beta(maxModes), mu(maxModes), alpha(maxModes), delt
  real*8, intent(in) :: B(maxModes), fm_initial_phase, Z0, Abase(2), Y(:), Estar, damping(maxModes), fm_amp_kp, fm_amp_ki, fm_pll_kp, fm_pll_ki, fm_pll_kd, fm_amp_kd, omegad_t_1
  integer*8, intent(in) :: pad, nstep, output_point_rate, numofcycle
  integer, intent(in) ::  numModes,  numincycle, z_feedback_choice
  integer :: i



  call rp_lib_put_str(driver, "output.string(Misc).about.label","Misc. Internal values", 0) 
  
  Str = "Internal values for reference. Most users can ignore these."

  write(s, *) CHAR(10), 'svn revision $Revision:on date $Date: 2021-05-18 21:04:11 -0400 (Tue, 18 May 2021) $ in folder $HeadURL: https://nanohub.org/tools/veda/svn/trunk/src/rappture_io.f90 $'
  Str = trim(Str) // s

  write(s,  *) CHAR(10), "omega_scale ", NonDimenTime(1d0), CHAR(10)
  Str = trim(Str) // s


  write(s,  *) CHAR(10), "omegad(1) ", DimenFreq( omegad_t_1 ) , " (rad/s) " , DimenFreq( omegad_t_1/2d0 /pi ), " Hz, at last computed point ", CHAR(10)
  Str = trim(Str) // s


  write(s,  *) CHAR(10), "Z0: ", DimenLength(Z0 * 1d9), " (nm) at first computed point", CHAR(10) !char(10) is newline
  Str = trim(Str) // s

  write(s,  *) "delt: ", delt, " (nd), ", DimenTime(delt), " (dim) at last computed point", CHAR(10) !char(10) is newline
  Str = trim(Str) // s


  write(s,  *) "pad: ", pad, CHAR(10) !char(10) is newline
  Str = trim(Str) // s

  write(s,  *) "numofcycle: ", numofcycle, CHAR(10) !char(10) is newline
  Str = trim(Str) // s


  write(s,  *) "numincycle: ", numincycle, CHAR(10) !char(10) is newline
  Str = trim(Str) // s

  write(s,  *) "nstep: ", nstep, CHAR(10) !char(10) is newline
  Str = trim(Str) // s

  write(s,  *) "output_point_rate: ", output_point_rate, CHAR(10) !char(10) is newline
  Str = trim(Str) // s


  write(s, '(A,2E13.5,A)') "Abase (nm): ", DimenLength(Abase(1)*1d9),  DimenLength(Abase(2)*1d9), "At first computed point"
  Str = trim(Str) // s // CHAR(10)

  write(s, '(A,E12.5, A, E12.5,A, E12.5, A,A)') "F1: ", DimenForce(F(1,1)), ", ", DimenForce(F(1,2)), ", ", DimenForce(F(1,3)), " (N) At first computed point", CHAR(10) !fixme, should be arb # modes
  Str = trim(Str) // s

  write(s, '(A,E10.3, A, E10.3, A,A)') "F2: ", DimenForce(F(2,1)),  ", ",DimenForce(F(2,2)), " (N) At first computed point", CHAR(10) !fixme, should be arb # modes
  Str = trim(Str) // s


  write(s, '(100E10.3)')  (phid(1,i), i=1,numModes)
  Str =  trim(Str) // "phid1: " // s // " (rad) At first computed point" // CHAR(10) 

  write(s, '(100E10.3)')  (phid(2,i), i=1,numModes)
  Str =  trim(Str) // "phid2: " // s // " (rad) At first computed point" // CHAR(10) 

  write(s, '(A,E10.3,A,A)') "aDMT: ", DimenLength(aDMT* 1d9), " (nm) At first computed point", CHAR(10) !fixme, should be arb # modes
  Str = trim(Str) // s

  write(s, '(A,E12.5,A,A)') "fm_initial_phase: ", fm_initial_phase * 180d0 / pi, " (deg)", CHAR(10)
  Str = trim(Str) // s

  write( s, '(100E12.5)') (alpha(i) ,i=1,numModes)
  Str = trim(Str)  // "alpha: " // s // CHAR(10)

  write( s, '(100E12.5)') (mu(i) ,i=1,numModes)
  Str = trim(Str) // "mu: " // s

  write( s, '(100E12.5)') (beta(i) ,i=1,numModes)
  Str = trim(Str) // CHAR(10) // "beta: " // s 

  write( s, '(100E12.5)') (B(i) ,i=1,numModes)
  Str = trim(Str) // CHAR(10) // "B: " // s 

  write( s, '(100E14.7)') (omegai(i) / (2 * pi),i=1,numModes)
  Str = trim(Str) // CHAR(10) // CHAR(10) //"Omega_i (natural freq): " // s // "(Hertz)" // CHAR(10) 

  write( s, '(100E12.5)') (Chi(i),i=1,numModes)
  Str = trim(Str) // CHAR(10) // "Chi: " // s 

  write( s, '(100E12.5)') (Keq(i),i=1,numModes)
  Str = trim(Str) // CHAR(10) // "Keq: " // s

  write( s, '(E12.5)') Estar 
  Str = trim(Str) // CHAR(10) // "Estar (subs): " // s  

  write( s, '(100E10.3)') (Y(i), i=1,(2*numModes))
  Str = trim(Str) // CHAR(10) // "Y final: " // s 

  write( s, '(100E10.3)') (damping(i), i=1,numModes)
  Str = trim(Str) // CHAR(10) // "damping: " // s 

  write( s, '(A,E12.5,A,E12.5,A,E12.5,A)') "fm_amp_kp: ", fm_amp_kp, " fm_amp_ki: ", fm_amp_ki, " fm_amp_kd: ", fm_amp_kd
  Str = trim(Str) // CHAR(10) // s 
  write( s, '(A,E12.5,A,E12.5,A,E12.5,A)')"fm_pll_kp: ", fm_pll_kp,  " fm_pll_ki: ", fm_pll_ki,  " fm_pll_kd: ", fm_pll_kd
  Str = trim(Str) // CHAR(10) // s 

  write( s, *) z_feedback_choice
  Str = trim(Str) // CHAR(10) //  "z_feedback_choice: " // s  

  write( s, '(E12.5)') DimenTime(minimum_relaxation)
  Str = trim(Str) // CHAR(10) // "min relaxation time (subs): " // s  

  
  call rp_lib_put_str(driver, "output.string(Misc).current", Str, 0)  

end subroutine OutputMiscParameters

!drk: I tried replacing this point-by-point output with buffering the whole result into an array
!and outputting everything later.  It was only a 3% performance improvement.  Didn't think that was
!worth the extra complexity it would add to the code.  It does help a lot with time histories though,
!and we should probably do something to speed up poincare maps too.
 subroutine OutputMainResults(fexcite, avg_impacts_per_cycle, Want_NumImpact, FPeakAtt, Keq,&
       FPeakRep, MeanForce, Phase, Ets, Indent, tc, omegai_nd, Amp, modulation_type, &
       Want_AZ,Want_P1,Want_MF,Want_PF,Want_ED,Want_EP, Want_I,Want_CT, Want_Fourier, Want_E_Anc, want_virial, virial, Eprop, &
       Edamp, Edrive, Ebs, E_anc, numModes, LockInR,LockInTh, RMS, want_RMS, useLockIn, omegad_nd, want_ev, xaxis, xaxis_1, AmpRedAvg)
   use Nondimensionalization   

       integer, intent(in) :: numModes, modulation_type, fexcite
       real*8, intent(in) :: FPeakAtt, Keq(maxModes), FPeakRep,  &
       MeanForce, Phase, Ets, Indent,tc, omegai_nd(1), Amp, Eprop(numModes), Ebs, &
       E_anc, Edamp(numModes), Edrive(numModes), LockInR, LockInTh, avg_impacts_per_cycle,&
       omegad_nd, xaxis, xaxis_1, virial, RMS
       real*8 :: energy_units_conversion

       logical, intent(in) :: Want_AZ,Want_P1,Want_MF,Want_PF,Want_ED,Want_EP, Want_I,Want_CT, &
                              Want_E_Anc, Want_Fourier, useLockIn, Want_NumImpact, want_ev, want_virial,want_RMS, AmpRedAvg

! the program is very sensitive to roundoff error.  The 6th significant digit 
! is very often meaningless.  So there is no put in outputting 12 sig digits.  Let's be more realistic.
! also, the trailing 'E3' forces fortran to use 3 digits for the exponent (default is to use only 2, and then
! if the exponents is bigger, it drops the E and outputs a number like 1.0E-100 as 1.0-100, which confuses the
! heck out of C++
       character(len=*), parameter :: fmtStr = '(ES16.7E3,ES16.7E3,A)'

       real*8 :: tmp_ph
       character*1 :: String1
       character*50 tmpStr
       integer :: iprop
       

       !this one is different x-axis because if axis=Amp, then outputting Amp vs Amp isn't very useful.
       if (AmpRedAvg) then
             write(tmpStr,fmtStr) xaxis_1,  DimenLength(Amp*1d9), char(10)
             call rp_lib_put_str(driver, "output.curve(Amp).component.xy", tmpStr, 1)

             write(tmpStr,fmtStr) xaxis_1,  DimenLength(LockInR*1d9), char(10) !reusing variable name. not really LockIn
             call rp_lib_put_str(driver, "output.curve(Amp_iter).component.xy", tmpStr, 1)          
       elseif (Want_AZ) then
          if (want_Fourier) then
             write(tmpStr,fmtStr) xaxis_1,  DimenLength(Amp*1d9), char(10)
             call rp_lib_put_str(driver, "output.curve(Amp).component.xy", tmpStr, 1)
          end if
          
          if (useLockIn) then 
             write(tmpStr,fmtStr) xaxis_1,  DimenLength(LockInR*1d9), char(10)
             call rp_lib_put_str(driver, "output.curve(LockinR).component.xy", tmpStr, 1)
          end if

       end if

       
       if (want_RMS) then
          write(tmpStr,fmtStr) xaxis,  DimenLength(sqrt(RMS)*1d9), char(10)
          call rp_lib_put_str(driver, "output.curve(RMSAmp).component.xy", tmpStr, 1)
       end if


       if (Want_NumImpact) then
          write(tmpStr,fmtStr) xaxis,  avg_impacts_per_cycle, char(10)
          call rp_lib_put_str(driver, "output.curve(NumImpact).component.xy", tmpStr, 1)
       end if
       
       if (Want_P1) then     
          if (want_Fourier) then
             write(tmpStr,fmtStr) xaxis,  Phase*180/pi, char(10)
             call rp_lib_put_str(driver, "output.curve(Phase).component.xy",tmpStr, 1)

          end if
          
          if (useLockin) then
!make sure lockin +Fourier int wrap the same way         
             tmp_ph = 90d0 - LockInTh*180/pi
             if (tmp_ph > 180) tmp_ph = tmp_ph - 360
             write(tmpStr,fmtStr) xaxis, tmp_ph , char(10)
             call rp_lib_put_str(driver, "output.curve(LockinPhase).component.xy",tmpStr, 1)
          end if           
       end if
       
       if (Want_PF) then
          write(tmpStr,fmtStr) xaxis, DimenForce(FPeakAtt)*1d9, char(10)
          call rp_lib_put_str(driver, "output.curve(FPeakAtt).component.xy",tmpStr, 1)     
          
          write(tmpStr,fmtStr) xaxis,  DimenForce(FPeakRep)*1d9, char(10)
          call rp_lib_put_str(driver, "output.curve(FPeakRep).component.xy",tmpStr, 1) 
       end if
       
       if (Want_MF) then       
          write(tmpStr,fmtStr) xaxis,  DimenForce(MeanForce)*1d9, char(10)
          call rp_lib_put_str(driver, "output.curve(MeanForce).component.xy",tmpStr, 1)
       end if       

       if (want_ev) then
          energy_units_conversion = electrons_per_coulomb ! eV/drive cycle
       else
          energy_units_conversion = 1e12 * DimenFreq(omegad_nd) / (2d0 * pi) !pW
       end if
       
       if (Want_ED) then       
          write(tmpStr,fmtStr) xaxis,  Ets*energy_units_conversion, char(10)                  
          call rp_lib_put_str(driver, "output.curve(Pts).component.xy",tmpStr, 1)
       end if
       
       if (Want_EP) then
          if (isAcoustic(fexcite)) then
             write(tmpStr,fmtStr) xaxis,  Ebs*energy_units_conversion, char(10)
             call rp_lib_put_str(driver, "output.curve(Ebs).component.xy",tmpStr, 1)
          end if
          
         do iprop=1,numModes
	   write (String1, '(I1)') iprop
           write(tmpStr,fmtStr) xaxis,  Eprop(iprop)*energy_units_conversion, char(10)
	   call rp_lib_put_str(driver, "output.curve(Eprop_" // String1 // ").component.xy",tmpStr, 1)

           write(tmpStr,fmtStr) xaxis,  Edamp(iprop)*energy_units_conversion, char(10)
	   call rp_lib_put_str(driver, "output.curve(Edmp_" // String1 // ").component.xy",tmpStr, 1)

           write(tmpStr,fmtStr) xaxis,  Edrive(iprop)*energy_units_conversion, char(10)
	   call rp_lib_put_str(driver, "output.curve(Edrv_" // String1 // ").component.xy",tmpStr, 1)
	 end do
       end if
              
       if (Want_E_Anc) then 
          write(tmpStr,fmtStr) xaxis,  E_Anc*energy_units_conversion, char(10)
          call rp_lib_put_str(driver, "output.curve(E_anc).component.xy",tmpStr, 1)
       end if
          
       if (Want_virial) then 
          write(tmpStr,fmtStr) xaxis,  virial*energy_units_conversion, char(10)
          call rp_lib_put_str(driver, "output.curve(virial).component.xy",tmpStr, 1)
       end if

       if (Want_I) then 
        write(tmpStr,fmtStr) xaxis,  DimenLength(Indent*1d9), char(10)
        call rp_lib_put_str(driver, "output.curve(Indent).component.xy",tmpStr, 1)
       end if
       
       if (Want_CT) then       
        write(tmpStr,fmtStr) xaxis,  DimenTime(tc*1d6), char(10)
        call rp_lib_put_str(driver, "output.curve(tcontact).component.xy",tmpStr, 1)
       end if     
     end subroutine OutputMainResults
 
 
 
 ! New subroutine added to generate the impact Poincare plot(s)
 subroutine OutputPoincare(numpoinc, path, Y, Npoinc_X, Npoinc_Y, d, dpr,t)
   use data1, only: ComputeZdist_dim
   use timeAndcycle, only: timeSinceTrans
   use Poincare, only: maxpoinc
   use Nondimensionalization
   implicit none
   
   integer, intent(in) :: numpoinc,  Npoinc_X(maxpoinc), Npoinc_Y(maxpoinc)
   real*8, intent(in) ::  Y(:), d, dpr,t
   real*8 :: X_p, Y_p
   integer ii
   character*1 PlotNum
   character*50 tmpStr
   character(len=*), intent(in) :: path
   character(len=*), parameter :: fmtStr = '(ES15.7,ES15.7,A)'
   
   do ii=1,numpoinc
      if (Npoinc_X(ii) == -2) then
         X_p = computeZdist_dim(t)
      elseif ( Npoinc_X(ii) == -1 ) then
         X_p = DimenTime(timeSinceTrans(t))
      elseif ( Npoinc_X(ii) == 0 ) then
         X_p = DimenLength(d*1d9)
      else
         X_p = DimenLength(y( Npoinc_X(ii) ) *1d9)
      end if

      if ( Npoinc_Y(ii) == 0 ) then
         Y_p = DimenVelocity(dpr * 1d6)
      else
         if ( 2*(Npoinc_Y(ii)/2) == Npoinc_Y(ii)) then
            !even means velo         
            Y_p = DimenVelocity(y( Npoinc_Y(ii) ) * 1d6)
         else
            !odd means disp
            Y_p = DimenLength(y( Npoinc_Y(ii) ) * 1d9)
         end if
      end if      
      
      write (PlotNum, '(I1)') ii
      write(tmpStr,fmtStr) X_p, Y_p, char(10)
      call rp_lib_put_str(driver, "output.curve(Poinc_" // path // PlotNum // ").component.xy", tmpStr, 1)    
   end do
   
 end subroutine OutputPoincare
 

subroutine OutputHigherFrequencies( del_time, numHH, an, bn, xaxis_a, xaxis_p, label, Want_A2, Want_P2, Want_Fourier)
  use NonDimensionalization

  logical, intent(in) :: Want_A2, Want_P2, Want_Fourier
  integer, intent(in) :: label, numHH
  real*8, intent(in) :: del_time, an(numHH), bn(numHH), xaxis_a, xaxis_p

  real*8 :: AmpHH, PhaseHH
  character*100 tmpStr
  character*2 Snum
  integer ik

  ! calculate amp and phase from a and b (real and imag parts)
  do ik = 1, numHH
  
     if ( ik >= 10) then
        write (Snum, '(I2)') ik
     else
        write (Snum, '(I1)') ik
     end if
    
    if (label < 3) then
       AmpHH = DimenLength( 2d0 / (del_time)*sqrt(an(ik)**2+bn(ik)**2) )
    else
       AmpHH = DimenForce( 2d0 / (del_time)*sqrt(an(ik)**2+bn(ik)**2) )
    end if

    PhaseHH = atan2(bn(ik),an(ik))

    write(tmpStr, '(E15.7,E15.7,A)' ) xaxis_a, AmpHH*1d9, char(10)

    if (label == 3) then
       call rp_lib_put_str(driver,"output.curve(FtsHarm).component.xy", tmpStr, 1)
    elseif ((label == 2) .and. Want_A2  .and. Want_Fourier) then
        call rp_lib_put_str(driver,"output.curve(Amp2).component.xy", tmpStr, 1)
    elseif ((label == 1) .and. Want_Fourier) then
       call rp_lib_put_str(driver,"output.curve(AmpHH" // Trim(Snum) // ").component.xy", tmpStr, 1)	 
    end if

    write(tmpStr, '(E15.7,E15.7,A)' ) xaxis_p, PhaseHH*180/pi, char(10)

    if ((label == 2) .and. Want_P2 .and. Want_Fourier) then
         call rp_lib_put_str(driver, "output.curve(Phase2).component.xy", tmpStr, 1)
    elseif ((label == 1) .and. Want_Fourier) then
       call rp_lib_put_str(driver, "output.curve(PhaseHH" // Trim(Snum) // ").component.xy", tmpStr, 1)
    end if
    
  end do

end subroutine OutputHigherFrequencies

subroutine OutputHigherHarmonicLockins( del_time, numHH, HigherHarmLockin, xaxis)
  use NonDimensionalization
  use LockInDataType

  integer, intent(in) :: numHH
  real*8, intent(in) :: del_time, xaxis
  type(LockinData), intent(in) :: HigherHarmLockin(:)

  real*8 :: AmpHH, PhaseHH
  character*100 tmpStr
  character*2 Snum
  integer ik

  do ik = 1, numHH
  
     if ( ik >= 10) then
        write (Snum, '(I2)') ik
     else
        write (Snum, '(I1)') ik
     end if
     
    AmpHH = DimenLength( HigherHarmLockin(ik)%R *1d9 )
  
    PhaseHH =  90d0 -  HigherHarmLockin(ik)%Th *180d0 / pi

    write(tmpStr, '(E15.7,E15.7,A)' ) xaxis, AmpHH, char(10)

    call rp_lib_put_str(driver,"output.curve(LockinHH" // Trim(Snum) // "R).component.xy", tmpStr, 1) 
    
    write(tmpStr, '(E15.7,E15.7,A)' ) xaxis, PhaseHH, char(10)
    
    call rp_lib_put_str(driver, "output.curve(LockinHH" // Trim(Snum) // "Th).component.xy", tmpStr, 1)
    
  end do

end subroutine OutputHigherHarmonicLockins


subroutine OutputBimodalLockin( BimodalLockin, xaxis_a, xaxis_p, want_a2, want_p2)
  use NonDimensionalization
  use LockInDataType

  real*8, intent(in) :: xaxis_a, xaxis_p
  type(LockinData), intent(in) :: BimodalLockin
  logical, intent(in) :: want_a2, want_p2

  real*8 :: Amp, Phase
  character*100 tmpStr

    if (want_a2) then
       Amp = DimenLength(BimodalLockin%R *1d9)
       write(tmpStr, '(E15.7,E15.7,A)' ) xaxis_a, Amp, char(10)       
       call rp_lib_put_str(driver,"output.curve(LockinAmp2).component.xy", tmpStr, 1)	     
    end if

    if (want_p2) then
       Phase = 90d0 -  BimodalLockin%Th *180d0 / pi    
       write(tmpStr, '(E15.7,E15.7,A)' ) xaxis_p, Phase, char(10)
       call rp_lib_put_str(driver, "output.curve(LockinPhase2).component.xy", tmpStr, 1)
    end if
end subroutine OutputBimodalLockin


subroutine OutputBimodal_AmpRed( Amp2, Amp2_iter, Phase2, xaxis_a, xaxis_p, want_a2, want_p2, virial1, virial2, dissipation1, dissipation2)
  use NonDimensionalization


  real*8, intent(in) :: xaxis_a, xaxis_p, Amp2, Amp2_iter, Phase2,  virial1, virial2, dissipation1, dissipation2
  logical, intent(in) :: want_a2, want_p2

  real*8 :: Amp, Phase,  energy_units_conversion
  character*100 tmpStr

  character(len=*), parameter :: fmtStr = '(ES16.7E3,ES16.7E3,A)'

    if (want_a2) then
       Amp = DimenLength(Amp2 *1d9)
       write(tmpStr, fmtStr ) xaxis_a, Amp, char(10)       
       call rp_lib_put_str(driver,"output.curve(Amp2).component.xy", tmpStr, 1)

       Amp = DimenLength(Amp2_iter *1d9)
       write(tmpStr, fmtStr ) xaxis_a, Amp, char(10)       
       call rp_lib_put_str(driver,"output.curve(Amp2_iter).component.xy", tmpStr, 1)	     
    end if

    if (want_p2) then
       Phase = Phase2 *180d0 / pi    
       write(tmpStr, '(E15.7,E15.7,A)' ) xaxis_p, Phase, char(10)
       call rp_lib_put_str(driver, "output.curve(Phase2).component.xy", tmpStr, 1)
    end if


!    if (want_ev) then
       energy_units_conversion = electrons_per_coulomb ! eV/drive cycle
!    else
!       energy_units_conversion = 1e12 * DimenFreq(omegad_nd) / (2d0 * pi) !pW
!    end if
    
    
    write(tmpStr,fmtStr) xaxis_p,  virial1*energy_units_conversion, char(10)                  
    call rp_lib_put_str(driver, "output.curve(virial1).component.xy",tmpStr, 1)
    write(tmpStr,fmtStr) xaxis_p,  virial2*energy_units_conversion, char(10)                  
    call rp_lib_put_str(driver, "output.curve(virial2).component.xy",tmpStr, 1)

    write(tmpStr,fmtStr) xaxis_p,  dissipation1*energy_units_conversion, char(10)                  
    call rp_lib_put_str(driver, "output.curve(Pts1).component.xy",tmpStr, 1)
    write(tmpStr,fmtStr) xaxis_p,  dissipation2*energy_units_conversion, char(10)                  
    call rp_lib_put_str(driver, "output.curve(Pts2).component.xy",tmpStr, 1)
    
  end subroutine OutputBimodal_AmpRed



!formerly "jump"
subroutine OutputFzCurvesNew(u,Force,d,dpr,Zdist,Zdot, jumpState,t, y,want_FzEig, Indent, wantAttard, surf_u0, iout, cur_props )
  use Nondimensionalization
   use data1, only : numModes
   use TriggeredFzMode, only: TRIGGER_REVERSE
   use matl_prop_type_module  
   implicit none

   integer, intent(in) :: jumpState
   integer*8, intent(in) :: iout
   logical, intent(in) :: want_FzEig, wantAttard
   real*8, intent(in) ::  u,Force,d,dpr, Zdist, t, y(:), Indent, Zdot, surf_u0
   type(matl_prop), intent(in) :: cur_props
   
   real*8 defl
   integer :: i
   character*50 tmpStr
   character*2 dir, n
   character(len=*), parameter :: fmtStr = '(E16.7E3,E16.7E3,A)'   !see comment in outputMainResults
   
   if (jumpState == TRIGGER_REVERSE) then
      dir = "_R"
   else
      dir = "_A"
   end if

   defl = DimenLength(u*1d9)
   
   write(tmpStr,fmtStr) Zdist, defl, char(10)
   call rp_lib_put_str(driver, "output.curve(DeflZ" // dir // ").component.xy", tmpStr, 1) 

   write(tmpStr,fmtStr) Zdist,  DimenLength(d*1d9), char(10)
   call rp_lib_put_str(driver, "output.curve(gapZ" // dir // ").component.xy", tmpStr, 1) 

   write(tmpStr,fmtStr) Zdist,  DimenLength(Indent*1d9), char(10)
   call rp_lib_put_str(driver, "output.curve(Indent" // dir // ").component.xy",tmpStr, 1)

   
   write(tmpStr,fmtStr) T,  Zdist, char(10)
   call rp_lib_put_str(driver, "output.curve(ZT).component.xy", tmpStr, 1) 

   !doesn't work correctly for the new sine wave option on FZ curves, but this was mostly for debugging anyway,
   !so just got rid of it.   
!   write(tmpStr,fmtStr) T,  Zdot, char(10)
!   call rp_lib_put_str(driver, "output.curve(ZdotT).component.xy", tmpStr, 1) 

   
   write(tmpStr,fmtStr) T,  DimenLength(d*1d9), char(10)
   call rp_lib_put_str(driver, "output.curve(gapT).component.xy", tmpStr, 1) 

   write(tmpStr,fmtStr) T,  DimenVelocity(dpr*1d6), char(10)
   call rp_lib_put_str(driver, "output.curve(dprT).component.xy", tmpStr, 1) 
   
   write(tmpStr,fmtStr) DimenLength(d*1d9), DimenForce(Force)*1d9, char(10)
   call rp_lib_put_str(driver, "output.curve(FtsGap).component.xy", tmpStr, 1)

   write(tmpStr,fmtStr) Zdist, DimenForce(Force)*1d9, char(10)
   call rp_lib_put_str(driver, "output.curve(FtsZ" // dir // ").component.xy", tmpStr, 1)

   write(tmpStr,fmtStr) t, DimenForce(Force)*1d9, char(10)
   call rp_lib_put_str(driver, "output.curve(FtsT).component.xy", tmpStr, 1)
   
   if ( want_FzEig) then
      do i = 1, numModes
         write(n , '(I1)') i
         write(tmpStr,fmtStr) t, DimenLength(y(2*i-1) * 1d9), char(10) 
         call rp_lib_put_str(driver, "output.curve(Eig" // trim(n) // ").component.xy", tmpStr, 1)

         write(tmpStr,fmtStr) t, DimenVelocity(y(2*i) * 1d6), char(10)
         call rp_lib_put_str(driver, "output.curve(EigVel" // trim(n) // ").component.xy", tmpStr, 1)

      end do
   end if

   if (wantAttard) then
      write(tmpStr,fmtStr) t, DimenLength(surf_u0)*1d9, char(10)
      call rp_lib_put_str(driver, "output.curve(AttardSurfPos).component.xy", tmpStr, 1)

      if (d < cur_props%attard_start_nd) call output_movie_frame(Y, driver, d, iout , cur_props, 1, t )      
   end if
   
 end subroutine OutputFzCurvesNew


subroutine DebugFM_Transient( t, u, Z0, omegad, LockInR, LockInTh, DrvAmp, phase_err, phase_err_dt, amp_err, amp_err_dt, theta, Zerr)
  use NonDimensionalization

   real*8, intent(in) ::  t, u, omegad, phase_err, phase_err_dt, amp_err, amp_err_dt, theta
   real*8, intent(in) :: Z0, LockInR, LockInTh, DrvAmp, Zerr

   character*50 tmpStr
   character(len=*), parameter :: fmtStr = '(E16.7E3,E16.7E3,A)'   !see comment in ouputMainResults

!  write(tmpStr,fmtStr) t,  theta, char(10)
!  call rp_lib_put_str(driver, "output.curve(thetaTrans).component.xy", tmpStr, 1)    

  write(tmpStr,fmtStr) t,  DimenLength(LockInR*1d9), char(10)
  call rp_lib_put_str(driver, "output.curve(lckR).component.xy", tmpStr, 1)  

  write(tmpStr,fmtStr) t,  phase_err*180/Pi, char(10)
  call rp_lib_put_str(driver, "output.curve(PhaseErrTrans).component.xy", tmpStr, 1)  

!  write(tmpStr,fmtStr) t,  phase_err_dt*180/Pi, char(10)
!  call rp_lib_put_str(driver, "output.curve(PhaseErrIntTrans).component.xy", tmpStr, 1)  

  write(tmpStr,fmtStr) t,  amp_err, char(10)
  call rp_lib_put_str(driver, "output.curve(AmpErrTrans).component.xy", tmpStr, 1)  

!  write(tmpStr,fmtStr) t,  amp_err_dt, char(10)
!  call rp_lib_put_str(driver, "output.curve(AmpErrIntTrans).component.xy", tmpStr, 1)  

  write(tmpStr,fmtStr) t,  LockInTh*180/Pi, char(10)
  call rp_lib_put_str(driver, "output.curve(lckTh).component.xy", tmpStr, 1)  

  write(tmpStr,fmtStr) t,  DimenLength(Z0*1d9), char(10)
  call rp_lib_put_str(driver, "output.curve(Z0).component.xy", tmpStr, 1)  
 
  write(tmpStr,fmtStr) t, Zerr, char(10)
  call rp_lib_put_str(driver, "output.curve(Zerr).component.xy", tmpStr, 1)  
 
  write(tmpStr,fmtStr) t,  DimenFreq(omegad) / (1d3 * 2d0 * pi) , char(10)
  call rp_lib_put_str(driver, "output.curve(omegad).component.xy", tmpStr, 1)            
 
  write(tmpStr,fmtStr) t,  DrvAmp , char(10)
  call rp_lib_put_str(driver, "output.curve(DrvAmpTrans).component.xy", tmpStr, 1)            

  
end subroutine DebugFM_Transient

subroutine DebugAMScan_Transient( t,Z0, Amp)
  use NonDimensionalization

   real*8, intent(in) ::  t, Z0, Amp

   character*50 tmpStr
   character(len=*), parameter :: fmtStr = '(E16.7E3,E16.7E3,A)'   !see comment in ouputMainResults
    
   write(tmpStr,fmtStr) DimenTime(t),  Amp, char(10)
   call rp_lib_put_str(driver, "output.curve(dbgAmp).component.xy", tmpStr, 1)  
  
   write(tmpStr,fmtStr) DimenTime(t),  DimenLength(Z0*1d9), char(10)
   call rp_lib_put_str(driver, "output.curve(Z0).component.xy", tmpStr, 1)  

 end subroutine DebugAMScan_Transient


subroutine DebugContactScan_Transient( t,Z0, Defl)
  use NonDimensionalization
   real*8, intent(in) ::  t, Z0, Defl

   character*50 tmpStr
   character(len=*), parameter :: fmtStr = '(E16.7E3,E16.7E3,A)'   !see comment in ouputMainResults
    
   write(tmpStr,fmtStr) DimenTime(t), DimenLength(Defl*1e9), char(10)
   call rp_lib_put_str(driver, "output.curve(dbgDefl).component.xy", tmpStr, 1)  
  
   write(tmpStr,fmtStr) DimenTime(t),  DimenLength(Z0*1d9), char(10)
   call rp_lib_put_str(driver, "output.curve(Z0).component.xy", tmpStr, 1)  

 end subroutine DebugContactScan_Transient


subroutine DebugForceModScan_Transient( t,Z0, Defl, R)
  use NonDimensionalization
   real*8, intent(in) ::  t, Z0, Defl,R

   character*50 tmpStr
   character(len=*), parameter :: fmtStr = '(E16.7E3,E16.7E3,A)'   !see comment in ouputMainResults
    
   write(tmpStr,fmtStr) DimenTime(t), DimenLength(Defl*1e9), char(10)
   call rp_lib_put_str(driver, "output.curve(dbgDefl).component.xy", tmpStr, 1)  

   write(tmpStr,fmtStr) DimenTime(t), DimenLength(R*1e9), char(10)
   call rp_lib_put_str(driver, "output.curve(dbgMeanDefl).component.xy", tmpStr, 1)  
  
   write(tmpStr,fmtStr) DimenTime(t),  DimenLength(Z0*1d9), char(10)
   call rp_lib_put_str(driver, "output.curve(Z0).component.xy", tmpStr, 1)  

 end subroutine DebugForceModScan_Transient


subroutine DebugPeakForceScan_Transient( t,Z0, Force)
  use NonDimensionalization
   real*8, intent(in) ::  t, Z0, Force

   character*50 tmpStr
   character(len=*), parameter :: fmtStr = '(E16.7E3,E16.7E3,A)'   !see comment in ouputMainResults
    
   write(tmpStr,fmtStr) DimenTime(t), DimenLength(Force*1e9), char(10)
   call rp_lib_put_str(driver, "output.curve(dbgPeak).component.xy", tmpStr, 1)  
  
   write(tmpStr,fmtStr) DimenTime(t),  DimenLength(Z0*1d9), char(10)
   call rp_lib_put_str(driver, "output.curve(Z0).component.xy", tmpStr, 1)  

 end subroutine DebugPeakForceScan_Transient


 subroutine DebugPeakForceScan_FullTransient( t,d,dpr,q1, Z0, Force)
  use NonDimensionalization
   real*8, intent(in) ::  t,d,dpr, Z0, Force, q1

   character*50 tmpStr
   character(len=*), parameter :: fmtStr = '(E16.7E3,E16.7E3,A)'   !see comment in ouputMainResults
    
   write(tmpStr,fmtStr) DimenTime(t), DimenLength(Force*1e9), char(10)
   call rp_lib_put_str(driver, "output.curve(dbgPeak_full).component.xy", tmpStr, 1)  

   write(tmpStr,fmtStr) DimenTime(t),  DimenLength(Z0*1d9), char(10)
   call rp_lib_put_str(driver, "output.curve(Z0_full).component.xy", tmpStr, 1)  

   write(tmpStr,fmtStr) DimenTime(t),  DimenLength(d*1d9), char(10)
   call rp_lib_put_str(driver, "output.curve(d_full).component.xy", tmpStr, 1)  

   write(tmpStr,fmtStr) DimenTime(t),  DimenVelocity(dpr*1d6), char(10)
   call rp_lib_put_str(driver, "output.curve(dpr_full).component.xy", tmpStr, 1)  

   write(tmpStr,fmtStr) DimenTime(t),  DimenLength(q1*1d9), char(10)
   call rp_lib_put_str(driver, "output.curve(q1_full).component.xy", tmpStr, 1)  

   
 end subroutine DebugPeakForceScan_FullTransient



subroutine OutputScanData( xaxis, Ztrack, Amp, setpoint, ZF, ErrorZ, fexcite, u, z_feedback_choice)
  use NonDimensionalization
  real *8, intent(in) :: xaxis, Amp, setpoint, Ztrack, ZF, ErrorZ, u
  integer, intent(in) :: fexcite, z_feedback_choice
  character*50 tmpStr
  character(len=*), parameter :: fmtStr = '(E16.7E3,E16.7E3,A)'   !see comment in ouputMainResults
  
  if ((fexcite .ne. NO_EXC) .and.(z_feedback_choice .ne. MEAN_DEFL_Z).and. (z_feedback_choice .ne. MAX_FORCE)) then
     write(tmpStr, fmtStr ) xaxis, DimenLength((Amp-setpoint)*1d9), char(10)
     call rp_lib_put_str(driver, "output.curve(AmpDist).component.xy", tmpStr, 1) !fixme, needs a better name
  elseif (z_feedback_choice .ne. MAX_FORCE) then
     write(tmpStr, fmtStr ) xaxis, DimenLength((setpoint - u)*1d9), char(10)
     call rp_lib_put_str(driver, "output.curve(DeflErr).component.xy", tmpStr, 1) !fixme, needs a better name
  end if
  
  write( tmpStr, fmtStr ) xaxis, DimenLength(ZF*1d9), char(10)
  call rp_lib_put_str(driver, "output.curve(Hsample).component.xy",  tmpStr, 1) 
  write( tmpStr, fmtStr) xaxis, DimenLength(Ztrack*1d9), char(10)
  call rp_lib_put_str(driver, "output.curve(Topog).component.xy",  tmpStr, 1)

  write( tmpStr, fmtStr) xaxis,  DimenLength(ErrorZ*1d9), char(10)
  call rp_lib_put_str(driver, "output.curve(ErrorZ).component.xy", tmpStr, 1)
  
end subroutine OutputScanData

subroutine OutputPeakForceData( xaxis, PeakForceR,  PeakForceA )
  use NonDimensionalization
  real *8, intent(in) :: xaxis, PeakForceA, PeakForceR
  character*50 tmpStr
  character(len=*), parameter :: fmtStr = '(E16.7E3,E16.7E3,A)'   !see comment in ouputMainResults
  
  write(tmpStr, fmtStr ) xaxis, DimenForce(PeakForceR)*1e9, char(10)
  call rp_lib_put_str(driver, "output.curve(PeakDistR).component.xy", tmpStr, 1) !fixme, needs a better name

  write(tmpStr, fmtStr ) xaxis, DimenForce(PeakForceA)*1e9, char(10)
  call rp_lib_put_str(driver, "output.curve(PeakDistA).component.xy", tmpStr, 1) !fixme, needs a better name

end subroutine OutputPeakForceData


subroutine OutputFMPlots( xaxis, omegad, phase_err, Drive_Signal)
  real *8, intent(in) :: xaxis, omegad, Drive_signal, phase_err
  character*50 tmpStr

  character(len=*), parameter :: fmtStr = '(E16.7E3,E16.7E3,A)'   !see comment in ouputMainResults

!  write(tmpStr, fmtStr ) xaxis, freq_shift, char(10)
!  call rp_lib_put_str( "output.curve(FreqShift).component.xy", tmpStr, 1) 
!

  write(tmpStr, fmtStr ) omegad, Drive_signal, char(10)
  call rp_lib_put_str(driver, "output.curve(AmpVsFreq).component.xy", tmpStr, 1) 

  write(tmpStr, fmtStr ) xaxis, omegad, char(10)
  call rp_lib_put_str(driver, "output.curve(DriveFreq).component.xy", tmpStr, 1) 

  write(tmpStr, fmtStr ) xaxis, Drive_signal, char(10)
  call rp_lib_put_str(driver, "output.curve(DriveAmp).component.xy", tmpStr, 1) 

  write(tmpStr, fmtStr ) xaxis, phase_err * 180d0 / (2 * pi) , char(10)
  call rp_lib_put_str(driver, "output.curve(PhaseErr).component.xy", tmpStr, 1) 

end subroutine OutputFMPlots

  
subroutine SetupIterationStatsPlots(exc_choice)
  integer, intent(in) :: exc_choice
  
 !                       (name,   label,                group,                 xlabel,            xunits, ylabel, yunits, scatter )
  call  setupgenericplot( 'Z',     'Z',               'Iteration Statistics Z', 'Iteration Number', '',     'Z',     'nm', .false.)
  call  setupgenericplot( 'Asp_t', 'Setpoint target', 'Iteration Statistics Z', 'Iteration Number', '',     'A/A0',  'A/A0', .false.)
  call  setupgenericplot( 'Asp_e', 'Setpoint error',  'Iteration Statistics Z', 'Iteration Number', '',     'A/A0',  'A/A0', .false.)

  if (exc_choice == BIMODAL) then
     !                       (name,   label,                group,                 xlabel,            xunits, ylabel, yunits, scatter )
     call  setupgenericplot( 'Asp2',   'Asp2',            'Iteration Statistics A2', 'Iteration Number', '',     'A/A0',  'A/A0', .false.)
     call  setupgenericplot( 'Z_A2',   'Z',               'Iteration Statistics A2', 'Iteration Number', '',     'Z',  'nm', .false.)
     call  setupgenericplot( 'Asp2_e', 'Setpoint error',  'Iteration Statistics A2', 'Iteration Number', '',     'A/A0',  'A/A0', .false.)
  end if
     
end subroutine SetupIterationStatsPlots

subroutine OutputIterationStatsPlots_Z(iteration_number, Z, Asp_target, Asp_error)
  use Nondimensionalization
  integer, intent(inout) :: iteration_number
  real*8, intent(in) :: Z, Asp_target, Asp_error
  character(len=*), parameter :: fmtStr = '(I8,ES16.7E3,A)'
  character*50 tmpStr
  
  write(tmpStr,fmtStr) iteration_number,  DimenLength(Z*1d9), char(10)
  call rp_lib_put_str(driver, "output.curve(Z).component.xy",tmpStr, 1 )

  write(tmpStr,fmtStr) iteration_number, Asp_target, char(10)
  call rp_lib_put_str(driver, "output.curve(Asp_t).component.xy",tmpStr, 1)

  write(tmpStr,fmtStr) iteration_number, Asp_error, char(10)
  call rp_lib_put_str(driver, "output.curve(Asp_e).component.xy",tmpStr, 1)

  iteration_number = iteration_number+1
end subroutine OutputIterationStatsPlots_Z

subroutine OutputIterationStatsPlots_A2(iteration_number, Asp2, Z , Asp2_error)
  use Nondimensionalization
  integer, intent(inout) :: iteration_number
  real*8, intent(in) :: Asp2, Z, Asp2_error
  character(len=*), parameter :: fmtStr = '(I8,ES16.7E3,A)'
  character*50 tmpStr
  !fixme arguments need better names here.  Z is really A.  confusing... copied from the tapping mode and need a fresh routine
  write(tmpStr,fmtStr) iteration_number,  Asp2, char(10)
  call rp_lib_put_str(driver, "output.curve(Asp2).component.xy",tmpStr, 1 )

  write(tmpStr,fmtStr) iteration_number, DimenLength(Z*1e9) , char(10)
  call rp_lib_put_str(driver, "output.curve(Z_A2).component.xy",tmpStr, 1)

  
  write(tmpStr,fmtStr) iteration_number, Asp2_error, char(10)
  call rp_lib_put_str(driver, "output.curve(Asp2_e).component.xy",tmpStr, 1)

  iteration_number = iteration_number+1
end subroutine OutputIterationStatsPlots_A2


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!1111
! utility functions

!note: if the specified item is not found, and strVal was not blank beforhand,
!rappture simply leaves it intact.  it doesn't clear it out.  I spent several hours
!tracking a bug down for this reason.  so generate this wrapper function so that 
!we always call rp_lib_get with a fresh blank string
character*100  function rp_lib_get_wrap( driver, path)
  character (len=*), intent(in) :: path
  integer, intent(in) :: driver

  rp_lib_get_wrap = ""
  call rp_lib_get( driver, path, rp_lib_get_wrap)
end function rp_lib_get_wrap

real*8 function readGenericDbl( status, element_to_read, short_name, InputEcho, units_convert)
  integer, intent(inout) :: status
  character(len=*), intent(in) :: element_to_read, short_name
  character*1500, intent(inout) :: InputEcho
  real*8, intent(in) :: units_convert

  character*100 :: strVal
  integer :: rp_units_convert_dbl
  real*8 :: value
  
  strVal = rp_lib_get_wrap(driver, element_to_read )
  status = status +  rp_units_convert_dbl(strVal," ", value)
  InputEcho = trim(InputEcho) // short_name // " = " // trim(strVal) // char(10)
  readGenericDbl = value / units_convert
end function readGenericDbl


!rappture's boolen read returns true if the boolean was not found.  this is a big backward compability issue
!got bitten by this several times.  don't use rappture's function at all if you can help it!
logical function daniel_get_boolean( path )
  character*10 :: strval
  character(len=*), intent(in) :: path
  integer :: rp_lib_get_boolean, x  
  strVal = rp_lib_get_wrap(driver, path )
  if (len_trim( strVal) == 0) then
     daniel_get_boolean = .false.
     call assert( .false., 'rp_lib_get_boolean ' // path // ' not found')
  else
     !rappture's read boolean function actually returns an integer.  boo!
     x = rp_lib_get_boolean(driver, path)
     if (x == 1) then 
        daniel_get_boolean = .true.
     else if (x == 0) then
        daniel_get_boolean = .false.
     else 
        call assert( .false., 'bad value from rp_lib_get_boolean')
     end if
  end if
end function daniel_get_boolean

integer function readGenericInteger( element_to_read)
   character(len=*), intent(in) :: element_to_read
   readGenericInteger = readGenericInteger(driver, element_to_read)
end function
Anosike


integer function stringVal(Input_Prefix)
   character(len=*), intent(in) :: Input_Prefix
   stringVal = stringVal(Input_Prefix)
  strVal = rp_lib_get_wrap(driver,  INPUT_PREFIX // "")
  InputEcho = trim(InputEcho) // " " // StrVal // char(10)
  status = status +  rp_units_convert_dbl(strVal," "," ") 
end function

!Rappture is very inefficient when it comes to output.  Writing out 1000 points 
!one at a time is maybe 10 times slower than writing out all of the points in 
!one big string. doesn't make much difference to normal output, but time history
!data can get quite big. So, we buffer all of the results in arrays and then write
! them all out a the end.  This function generates the string from the arrays
  subroutine GenerateOutputString(n, xdata, ydata, buf)
    implicit none
    integer, parameter :: LineLength = 50
    character(len=*), parameter :: fmtStr = '(E16.7E3,E16.7E3,A)'   !see comment in ouputMainResults
    character(len=n * LineLength), intent(out) :: buf
    character(len=LineLength) tmp
    character(len=1), parameter :: newline = char(10)
    integer, intent(in) :: n
    real*8, intent(in) :: xdata(n), ydata(n)
    integer :: i
    
    do i = 1, n
       write(tmp,fmtStr)  xdata(i), ydata(i), newline
       buf( ((i-1)*LineLength+1):(i*LineLength) ) = tmp
    end do
  end subroutine GenerateOutputString


!how to add a new time history output:
!1) expand the outputbuffer in setupTimeHistoryPlots
!2) add a call to SetupGenericPlot in setupTimeHistoryPlots
!3) pass the new data into outputTimeHistory_buffered and store it in the array
!4) add a call in FlushTimeHistory to output it to rappture.
subroutine OutputTimeHistory_buffered(tau, u, dpr, Force,d, y, drvF, pnt, cur_props, contact_radius, Z0, yprime, attard_surf_coord, NumRES1, sample_topo, do_attard_movie_frame)
   use matl_prop_type_module
   use data1, only : numModes
   use Nondimensionalization

   implicit none
   type(matl_prop), intent(in) :: cur_props
   integer*8, intent(in) ::  pnt
   integer, intent(in) :: NumRES1
   real*8, intent(in) ::  tau, u, dpr, Force,d, drvF, contact_radius, Z0,  attard_surf_coord
   real*8, intent(in) :: y(:), yprime(:), sample_topo
   logical, intent(in) ::  do_attard_movie_frame

   integer :: i, k
   character(len=*), parameter :: fmtStr = '(E16.7E3,E16.7E3,A)'   !see comment in ouputMainResults

  k = 1

  if (th_set%want_contact_area) then
     outputbuffer( k,  buffer_ndx, 1, pnt) = tau  
     outputbuffer( k,  buffer_ndx, 2, pnt) = DimenLength(contact_radius*1d9)
     k= k+1
  end if

  if (th_set%want_surf_relax) then
     !for attard module only
     outputbuffer( k,  buffer_ndx, 1, pnt) = tau

     !old way. only good for the spatial discr
     !outputbuffer( k,  buffer_ndx, 2, pnt) = DimenLength(y(2*numModes+1)*1d9)

     outputbuffer( k,  buffer_ndx, 2, pnt) = DimenLength(attard_surf_coord *1d9)
     k= k+1
     
     if ( do_attard_movie_frame ) call output_movie_frame(Y, driver, d, pnt , cur_props, buffer_ndx, tau)
     
  end if

  if (th_set%want_TH_obs_defl) then
     outputbuffer( k,  buffer_ndx, 1, pnt) = tau  
     outputbuffer( k,  buffer_ndx, 2, pnt) = DimenLength(u*1d9) !defl
     k= k+1
  end if

  if (th_set%want_th_force_gap) then
     outputbuffer( k,  buffer_ndx, 1, pnt) = DimenLength(d*1d9)
     outputbuffer( k,  buffer_ndx, 2, pnt) = DimenForce(Force*1d9)
     k = k+1
  else
     outputbuffer( k,  buffer_ndx, 1, pnt) = tau  ! Fts
     outputbuffer( k,  buffer_ndx, 2, pnt) = DimenForce(Force*1d9)
     k = k+1
     outputbuffer( k,  buffer_ndx, 1, pnt) = tau  ! gap
     outputbuffer( k,  buffer_ndx, 2, pnt) = DimenLength(d*1d9)
     k=k+1
  end if

  outputbuffer( k,  buffer_ndx, 1, pnt) = DimenLength(d*1d9)  !Pspace
  outputbuffer( k,  buffer_ndx, 2, pnt) = DimenVelocity(dpr*1d6)
  k = k+1

  if (th_set%want_modal) then
     do i = 1, numModes
        outputbuffer( k,  buffer_ndx, 1, pnt) = tau  !eig
        outputbuffer( k,  buffer_ndx, 2, pnt) = DimenLength(y(2*i - 1)*1d9)
        k = k+1
     end do
  end if

  if (th_set%want_self_exc) then
      outputbuffer( k,  buffer_ndx, 1, pnt) = tau
      outputbuffer( k,  buffer_ndx, 2, pnt) = DimenForce(1d9* drvF)
      k=k+1
   end if

   if (th_set%want_voltage) then 
      outputbuffer( k,  buffer_ndx, 1, pnt) = tau
      outputbuffer( k,  buffer_ndx, 2, pnt) = Y(2*numModes+1)
      k =k+1
   end if

   if (th_set%want_acceleration) then 
      outputbuffer( k,  buffer_ndx, 1, pnt) = tau
      outputbuffer( k,  buffer_ndx, 2, pnt) = DimenAcceleration(sum(yprime(2:(2*numModes):2)))*1d3
      k =k+1
   end if

   !think this is only used for peak force
   if (th_set%want_FZ) then 
      outputbuffer( k,  buffer_ndx, 1, pnt) = Z0*1d9
      outputbuffer( k,  buffer_ndx, 2, pnt) = DimenForce(Force)*1d9
      k =k+1
   end if

   if (th_set%want_PerfMetrics) then 
      outputbuffer( k,  buffer_ndx, 1, pnt) = tau
      outputbuffer( k,  buffer_ndx, 2, pnt) = numRes1
      k =k+1
   end if

   if (th_set%want_SampleTopo) then 
      outputbuffer( k,  buffer_ndx, 1, pnt) = tau
      outputbuffer( k,  buffer_ndx, 2, pnt) = sample_topo
      k =k+1
   end if


end subroutine OutputTimeHistory_buffered

subroutine FlushTimeHistory(do_fft)
   use data1, only : numModes
   use timeAndCycle, only: numincycle,hist_cycles
   use params
   use fft

   integer :: j, i, k, ierr
   character*1  SE, ST
   character(len = buffer_len4 * 50)  tmp
   logical, intent(in) :: do_fft

   if (.not. allocated(outputbuffer) ) then
      call assert(.false., "tried to flush history but it was not allocated")
      return
   end if
   do j = 1,buffer_ndx

      k = 1
      write (ST, '(I1)') j

       if (th_set%want_contact_area) then
          call GenerateOutputString(buffer_len4, outputbuffer(k,j,1,:), outputbuffer(k,j,2,:), tmp)
          call rp_lib_put_str(driver, "output.curve(cr" // ST // ").component.xy", tmp,1)   
          k=k+1
       end if


       if (th_set%want_surf_relax) then
          if (do_fft) then
            call my_fft(outputbuffer(k,j,2,:), ierr)
            outputbuffer(k,j,1,:) = outputbuffer(k,j,1,:) * numincycle / hist_cycles
         end if

         
         call GenerateOutputString(buffer_len4, outputbuffer(k,j,1,:), outputbuffer(k,j,2,:), tmp)
         call rp_lib_put_str(driver, "output.curve(surf" // ST // ").component.xy", tmp,1)
         k = k+1
      end if


      if (th_set%want_TH_obs_defl) then

         if (do_fft) then
            call my_fft(outputbuffer(k,j,2,:), ierr)
            outputbuffer(k,j,1,:) = outputbuffer(k,j,1,:) * numincycle / hist_cycles
         end if

         call GenerateOutputString(buffer_len4, outputbuffer(k,j,1,:), outputbuffer(k,j,2,:), tmp)
         call rp_lib_put_str(driver, "output.curve(defl" // ST // ").component.xy", tmp,1)
         k = k+1
      end if

      if (th_set%want_th_force_gap) then
         call GenerateOutputString(buffer_len4, outputbuffer(k,j,1,:), outputbuffer(k,j,2,:), tmp)
         call rp_lib_put_str(driver, "output.curve(Fts" // ST // ").component.xy", tmp,1)
         k = k+1
      else
         if (do_fft) then
            call my_fft(outputbuffer(k,j,2,:), ierr)         
            outputbuffer(k,j,1,:) = outputbuffer(k,j,1,:) * numincycle / hist_cycles
         end if 

         call GenerateOutputString(buffer_len4, outputbuffer(k,j,1,:), outputbuffer(k,j,2,:), tmp)
         call rp_lib_put_str(driver, "output.curve(Fts" // ST // ").component.xy", tmp,1)
         k = k+1

         if (do_fft) then
            call my_fft(outputbuffer(k,j,2,:), ierr)         
            outputbuffer(k,j,1,:) = outputbuffer(k,j,1,:) * numincycle / hist_cycles
         end if 

         call GenerateOutputString(buffer_len4, outputbuffer(k,j,1,:), outputbuffer(k,j,2,:), tmp)
         call rp_lib_put_str(driver, "output.curve(gap" // ST // ").component.xy", tmp,1)
         k =k+1
      end if
      call GenerateOutputString(buffer_len4, outputbuffer(k,j,1,:), outputbuffer(k,j,2,:), tmp)
      call rp_lib_put_str(driver, "output.curve(Pspace" // ST // ").component.xy", tmp,1)   
      k=k+1

      if (th_set%want_modal) then
         do i =1, numModes
            write (SE, '(I1)') i
            call GenerateOutputString(buffer_len4, outputbuffer(k,j,1,:), outputbuffer(k,j,2,:), tmp)
            call rp_lib_put_str(driver, "output.curve(eig" // SE // ST // ").component.xy", tmp,1)   
            k=k+1
         end do
      end if
      
       if (th_set%want_self_exc) then
          call GenerateOutputString(buffer_len4, outputbuffer(k,j,1,:), outputbuffer(k,j,2,:), tmp)
          call rp_lib_put_str(driver, "output.curve(drvforce" // ST // ").component.xy", tmp,1)   
          k=k+1
       end if

       if (th_set%want_voltage) then
          call GenerateOutputString(buffer_len4, outputbuffer(k,j,1,:), outputbuffer(k,j,2,:), tmp)
          call rp_lib_put_str(driver, "output.curve(voltage" // ST // ").component.xy", tmp,1)   
          k=k+1
       end if

       if (th_set%want_acceleration) then
          call GenerateOutputString(buffer_len4, outputbuffer(k,j,1,:), outputbuffer(k,j,2,:), tmp)
          call rp_lib_put_str(driver, "output.curve(TipAcceleration" // ST // ").component.xy", tmp,1)   
          k=k+1
       end if


       if (th_set%want_FZ) then
          call GenerateOutputString(buffer_len4, outputbuffer(k,j,1,:), outputbuffer(k,j,2,:), tmp)
          call rp_lib_put_str(driver, "output.curve(FvsZ" // ST // ").component.xy", tmp,1)   
          k=k+1
       end if

       if (th_set%want_PerfMetrics) then
          call GenerateOutputString(buffer_len4, outputbuffer(k,j,1,:), outputbuffer(k,j,2,:), tmp)
          call rp_lib_put_str(driver, "output.curve(NumRES1" // ST // ").component.xy", tmp,1)   
          k=k+1
       end if


       if (th_set%want_sampletopo) then
          call GenerateOutputString(buffer_len4, outputbuffer(k,j,1,:), outputbuffer(k,j,2,:), tmp)
          call rp_lib_put_str(driver, "output.curve(Topo" // ST // ").component.xy", tmp,1)   
          k=k+1
       end if

   end do


   if (allocated(outputbuffer)) deallocate(outputbuffer)

 end subroutine FlushTimeHistory

    
    !can make a movie within rappture as a "sequence"
    !probably should be moved to rappture_io.f90
    !fixme, add the tip shape as a second curve
    subroutine output_movie_frame(Y, driver, d, framenumber, cur_props, buffer_ndx, t)
      use data1, only: numModes, NEQ
      use matl_prop_type_module
      use NonDimensionalization
      use Viscoelasticity_attard , only: compute_u_cur
     type(matl_prop), intent(in) :: cur_props
     real*8, intent(in) :: Y(:), d,t
     integer, intent(in) :: driver, buffer_ndx
     integer*8, intent(in) :: framenumber
     real*8 :: u(cur_props%N_attard_spatial)
     integer :: i, n1, N
     character*1 ST
          
     character*20 frameStr 
     character(len=*), parameter :: fmtStr = '(ES16.7E3,ES16.7E3,A)'
     character(len=*), parameter :: fmtStr1 = '(I8,ES16.7E3,A)'
     character(len = cur_props%N_attard_spatial * 50)  tmp

     N = cur_props%N_attard_spatial
     
     write (ST, '(I1)') buffer_ndx
     
     n1 = 2*numModes+1
     
     write(frameStr, '(E16.7)') t
     call rp_lib_put_str(driver, "output.sequence(AttardMovie" // ST // ").element(" // frameStr // ").index", frameStr, 0)

     call rp_lib_put_str(driver, "output.sequence(AttardMovie" // ST // ").element(" // frameStr // ").curve(Surf).xaxis.label", "Radius", 0)
     call rp_lib_put_str(driver, "output.sequence(AttardMovie" // ST // ").element(" // frameStr // ").curve(Surf).xaxis.units", "nm", 0)

     call rp_lib_put_str(driver, "output.sequence(AttardMovie" // ST // ").element(" // frameStr // ").curve(Surf).yaxis.label", "Tip / Surface coordinate", 0)
     call rp_lib_put_str(driver, "output.sequence(AttardMovie" // ST // ").element(" // frameStr // ").curve(Surf).yaxis.units", "nm", 0)

     call rp_lib_put_str(driver, "output.sequence(AttardMovie" // ST // ").element(" // frameStr // ").curve(Surf).about.label", "Sample Surface", 0)

     call rp_lib_put_str(driver, "output.sequence(AttardMovie" // ST // ").element(" // frameStr // ").curve(Tip).about.label", "Tip", 0)

     
     if (cur_props%N_attard_fourier == 0) then
        u = Y(n1:NEQ) !input is directly the deflection

        call GenerateOutputString(cur_props%N_attard_spatial, 1d9 * DimenLengthArray(cur_props%rr,N ), 1d9 * DimenLengthArray(u,N), tmp)                   
        call rp_lib_put_str(driver, "output.sequence(AttardMovie" // ST // ").element(" // frameStr // ").curve(Surf).component.xy", tmp, 1)


        call GenerateOutputString(cur_props%N_attard_spatial, 1d9 * DimenLengthArray(cur_props%rr,N ), 1d9 * DimenLengthArray(d + cur_props%tip_shape,N), tmp)           
        call rp_lib_put_str(driver, "output.sequence(AttardMovie" // ST // ").element(" // frameStr // ").curve(Tip).component.xy", tmp, 1)

     else
        if (.true.) then           
           u = compute_u_cur(Y(n1:NEQ), cur_props) !in this case the input is the fourier coeff
           
           call GenerateOutputString(cur_props%N_attard_spatial, 1d9 * DimenLengthArray(cur_props%rr,N ), 1d9 * DimenLengthArray(u,N), tmp)           
           call rp_lib_put_str(driver, "output.sequence(AttardMovie" // ST // ").element(" // frameStr // ").curve(Surf).component.xy", tmp, 1)
           
           
           call GenerateOutputString(cur_props%N_attard_spatial, 1d9 * DimenLengthArray(cur_props%rr,N), 1d9 * DimenLengthArray(d+cur_props%tip_shape,N), tmp)
           call rp_lib_put_str(driver, "output.sequence(AttardMovie" // ST // ").element(" // frameStr // ").curve(Tip).component.xy", tmp, 1) 
        else
           !direct output of the fourier coeff, for debugging only
           do i = 1, cur_props%N_attard_fourier
              write(tmp,fmtStr1) i , 1d9 * DimenLength(Y(n1+i)), char(10)
              
              call rp_lib_put_str(driver, "output.sequence(AttardMovie" // ST // ").element(" // frameStr // ").curve(Surf).component.xy", tmp, 1)
           end do
        end if
     end if
     
     
   end subroutine output_movie_frame

   
   subroutine start_movie(driver, SA, ST, time_units)
     integer, intent(in) :: driver
     character*16, intent(in) :: SA
     character(len=*), intent(in) :: time_units
     character*1, intent(in) :: ST

     call rp_lib_put_str(driver, "output.sequence(AttardMovie" // ST // ").about.label", "Surface Deformation Movie " // SA , 0)

     call rp_lib_put_str(driver, "output.sequence(AttardMovie" // ST // ").index.label",time_units , 0)
     
   end subroutine start_movie

 

end module rappture_io

!this is a hack. can't be inside the module because of dependencies issues
subroutine WriteFatalError(message)
  use timeHistory, only : Amp_index,  wantHist
  use rappture_io, only: driver,   FlushTimeHistory
  use params

  character(len=*), intent(in) :: message
 ! if (debugging) then
     !this is best for developers, so we can see debugging transients
     call rp_lib_put_str(driver, "output.string(ErrorMessage).about.label","ErrorMessage", 0) 
     call rp_lib_put_str(driver, "output.string(ErrorMessage).current", message, 0)

     if ((wantHist) .and. (Amp_index > 1)) call FlushTimeHistory(.false.) !fixme, should be able to know do_fft instead of hard code false
 
     call rp_result(driver)
     stop 0

!4/24/2012.  this used to work, but it appears to be broken in the newest version of rappture.
!i've filed a bug report and going back to the old way for now.  
 ! else
 !    !this behavior will be better for the average user.
 !    write(0,*) message !0 is stderr in fortran
 !    stop 1  !this tells rappture that we failed.
 ! end if

end subroutine WriteFatalError

