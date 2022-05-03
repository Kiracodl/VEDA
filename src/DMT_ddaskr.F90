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
! DMT_DDASKR PROGRAM: main file for VEDA 

!coding conventions:
!1) try to keep separate variables for the dimensional and non-dimensional quantities.
!e.g. use Ainitial_dim for the dimension quanity and then use Ainitial_nd for the non-dimensional quantity.
!helps to keep clear.  Note when I did this, the original variable Ainitial was deleted to make it easier to make
!sure that I had changed all references correctly.  variables without a suffix are most likely non-dimensional
!2) Prior to normalization, all parameters are in SI units MKS (meters-kilograms-seconds). most units
!are not SI in the GUI, the conversion is done in rappture_io.f90
!3) all oscillatory amplitudes are peak values, not peak-peak or RMS.
!4) prefix "want" indicates a logical (boolean) value

!wish list:
! should we supply a Jacobian routine to DDASKR?  Right now it is approximating
!derivatives with numerical differences.  We could trade two calls to RES1 for one call to JAC.
!may be faster.
    
      
program DMT_DDASKR
  use params
  use rappture_io
  use data0
  use data1
  use data2
  use contactModels
  use Viscoelasticity_ting
  use Scanning
  use Controller
  use Poincare
  use AppControl
  use TriggeredFzMode
  use FreqMod
  use LockIn
  use Approaching
  use ForcingCalculations
  use TimeAndCycle
  use BandPassFilter
  use HydrodynamicFunction
  use AntiAliasingFilter
  use Nondimensionalization
  use SelfExcitation
  use TimeHistory
  use tip_data
  use forcevolume
  use peakForceMod
  use checkFzCurves
#ifdef __INTEL_COMPILER
  use IFPORT !needed for rand() for ifort but not g95 or gfortran
#endif

  IMPLICIT none 

  !functions and subs
  EXTERNAL RES1, RT1, rk4, ab4
!  integer :: rk4, ab4

  character*100 inFile, xlabel, xunits

  integer Nprog2, Zrange, IDID,  INFO(20), LIW, LRW,  & 
       i,  xchoice,  numHH, z_feedback_choice, num_impacts, output_type, numincycle_dropdown, numincycle_direct, numincycle_ve, openmp_num_threads
  integer, parameter :: solver = SOLV_DDASKR

  integer :: JROOT(1)
  integer, allocatable :: IWORK(:)

  integer*8 :: Ihist,nstep,nstepFZ,iout_prevOutput, debug_skip, numofcycle
  integer*8, parameter :: progress_rate = 1000

  logical timeHC, checkh,  WantHH, want_autoCalcIC, scatter, want_TH_obs_defl, want_th_force_gap, do_fft, want_AutoCalcTC, want_FZ, want_acceleration
  logical Want_A2, Want_P2, Want_AZ, Want_P1, Want_MF, Want_PF, Want_ED, Want_Fourier, want_custom_B, flag, skip
  logical Want_E_Anc, Want_ForceFourier, want_ev, want_virial, want_RMS
  logical Want_I, Want_CT, Want_FzEig, Want_ContTrans, Want_NumImpact
  logical CalcInputK, AutoCalcOmega, AutoCalcAlpha, Want_EP, WantConstZ
  logical Want_Abase_direct, want_hydrodynamic_function, Want_NormFreqShift
  real*8 :: Abase_input(maxExc), F_input(maxModes), basevel, B_input(maxModes)
  real*8  Ainitial_nd(maxModes), accum3, accum4, transient_allowance

  REAL*8  omegad_start, omegad_stop, T, TOUT,   u, tfoo, &
       d, xaxis, xaxis_1, xaxis_2, FPeakAtt, MeanForce, dpr, dpr_old, min_d, & 
       Force, Force_old, Indent, FpeakRep, T_prevOutput, theta_prevOutput, Phase, omegad_dim(maxModes),   & 
       tc, theta_star, t_star, tau, omegas_dim, Asample_dim,  ASetPt, Noise, avg_impacts_per_cycle, &
       Phase_initial, new_omegad, fluid_rho, fluid_eta, cant_width, cant_len,  Y_IC(2*maxModes), F_ForceVol_dim, &
       Zbase_Amp_dim,  peakF_setpoint_dim, omegai_khz(maxModes), transient_timeout_in

  real*8, allocatable :: ATOL(:), RTOL(:), RWORK(:), Y(:), YPRIME(:)

  !real and imag parts of running fourier integrals
  Real*8 a1, b1, a2(1), b2(1),  F1r(1), F1i(1)
  real*8, allocatable :: NHH(:), an(:), bn(:)

  !note on a2,b2: we made 'an' and 'bn' to be arrays, a2 and b2 were still scalars, and they were being passed to the 
  !same routine (OutputHigherFrequencies).  to be technically correct, we can't mix arrays and scalars, so a2 and b2 have
  !to be arrays, just they are arrays of 1 element.

  character*3000 inputEcho, tmpStr

  real*8 rpar !dummy arg. required by ddaskr but we don't use
  integer ipar, jdum, psdum ! also not used
    
  ! Energy content paramters:	     
  REAL*8 Ets, Eprop(maxModes), Edamp(maxModes), Edrive(maxModes), Ebs, E_anc, Virial

  
  call getarg(1,inFile)

  call OpenInputFile(inFile)
  
  Checkh = .false.
  timeHC = .false.
  Ihist = 0
  Amp_index = 1
  InputEcho  = ""
  !getting a lot of weird floating point errors if these are left uninitialized
  omegai_khz = 0
  peakF_setpoint_dim = 0
  F_ForceVol_dim=0
  omegad_dim=0
!	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!	Parameters
!	Operating conditions and cantilever properties
!	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  call ReadOperatingParameter( fexcite, exc_choice, sweepchoice, numModes, &
       operating_mode, CalcInputK, AutoCalcOmega, modulation_type, AutoCalcAlpha, output_type, InputEcho)

  call ReadOpCondAndCantProps(modulation_type, operating_mode, exc_choice,fexcite,numModes, omegad_dim, &
       omegad_start, omegad_stop, sweep_time,Ainitial_dim,AprchS_dim, mtip, want_AutoCalcTC, & 
       LockInOrder, LockInTC, gamma_drag, want_Abase_direct, Abase_input,F_input, &
       want_autoCalcIC, Y_IC, Afluid, mstar_div_m, self_exc_gain, self_exc_phase, self_exc_saturation, Asample_dim, omegas_dim, wantSampleExc, InputEcho)

  if (numModes > maxModes) numModes = maxModes

  call readCantModalProps(  numModes, fexcite, omegai_khz,  Keq, Quality, CalcInputK, output_type, Chi, &
       AutoCalcOmega, AutoCalcAlpha, alpha, osc_Q, osc_omega_dim,phase_slope_dim, efficiency_slope_dim, want_nonideal_magnetic, &
       want_custom_B, B_input, inputEcho)
  
  if ((operating_mode == SCAN) .or. (modulation_type == FREQUENCY)) then
     if (modulation_type == PEAK_FORCE) then
        !samples exactly once per drive cycle, just used in calculating "pad"
        sample_freq_hz = omegad_dim(1) 
     else
        call ReadSampleFreq( sample_freq_hz, InputEcho)      
     end if
  end if
  
        if (operating_mode == SCAN) then
           if (modulation_type == FORCE_VOL) then 
              call ReadForceVol(LineSpeed_dim, F_ForceVol_dim, ForceVolSettleTime, z_feedback_choice, InputEcho)          
           else
              call ReadControllerProperties( LineSpeed_dim, NoiseAmp, KP, KI, WantConstZ, z_feedback_choice, modulation_type, transient_timeout_in, InputEcho)
              if (modulation_type == PEAK_FORCE) then
                 call ReadPeakForce(Zbase_Amp_dim, peakF_setpoint_dim, Want_FZ, InputEcho )
              else
                 Zbase_Amp_dim = 1d0 !for avoiding errors
              end if
           end if
           call readFeatureProperties( FeatType, HF, LF, LF2,SubsLen_dim, WantTSCON, WantFProp, InputEcho)
        else if (isFzCurves(operating_mode, fexcite)) then
           call readTriggeredFzMode(  triggerDeflLim, piezoReverseTime_pct, want_FzEig, InputEcho, fzshape, CurveRate_dim)
        end if

        
        call ReadTipData( Rtip_dim, Etip, Poisson_tip, tip_angle, tip_shape, InputEcho)
        call ReadSampleData(  substrate_props, "input.phase(ts)",  InputEcho)
        if (( operating_mode == SCAN) .and. WantFProp) then
           call ReadSampleData(  feature_props, "input.phase(feature).group(fp)",  InputEcho)
        end if        

        
        call ReadSimulationParameters(Zrange,Wanthist,WantHH,Want_Strob,Want_Impact,numpoinc,numHH,Z0,Zf,plotpnts,numincycle_dropdown, numincycle_direct ,numHist,hist_cycles,xchoice, want_freqswp_sp_aprch, ASetPt, operating_mode, transient_allowance, want_hydrodynamic_function, Asp_f, freqshift_f, modulation_type, openmp_num_threads, InputEcho)

#if defined(OPENMP)
        call OMP_set_num_threads(openmp_num_threads) !this is just for attard right now
#endif
        
        allocate( an( numHH))  !to avoid compiler errors, it's okay to allocate these with size 0
        allocate( bn( numHH))
        allocate( NHH( numHH))
        allocate( HigherHarmLockin(numHH))
	if (Wanthist) call Read_TimeHist( numHist, Ahist, Wanthist_byA, want_TH_obs_defl, want_th_force_gap, do_fft, want_acceleration)
	if (WantHH)   call Read_HH( numHH, NHH, z_feedback_choice)
           
	call ReadPoincareParameters(Want_Strob,Want_Impact,numpoinc,Npoinc_X,Npoinc_Y)
	call ReadPlotChoices(Want_A2,Want_P2,Want_AZ,Want_P1,Want_MF,Want_PF,Want_ED,Want_EP,Want_I,Want_CT, Want_Fourier, Want_ContTrans, Want_E_Anc, Want_ForceFourier, want_ev, want_virial, want_RMS )

        if (modulation_type == FREQUENCY) then
           call ReadFMData( fm_pll_kp, fm_pll_ki, fm_pll_kd, fm_amp_kp, fm_amp_ki,fm_amp_kd,  freq_shift_sp, fm_direct_control, z_feedback_choice, FM_want_noncontact, WantCalcPLLGains, WantCalcAmpGains, want_pre_BPF, Want_NormFreqShift, InputEcho )
        else
           want_pre_BPF = .false.
        end if
 
         !we already have a FatalError for this below, so this is redundant.  Plus, we are applying this check in all 
         !situations while it really only applies for force modulation -- Daniel
!        call assert(omegad_dim(1) .GT. omegai(1)*1.0d-1, 'Drive frequency is not less than 10 percent of natural frequency')

        if ( want_hydrodynamic_function ) then
           call ReadHydrodynamicFunction(fluid_rho, fluid_eta, cant_width, cant_len, InputEcho)
        end if

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!no more input to read after this point

        
        if (( Zrange == ZRANGE_ASP) .and. ( operating_mode == APPR_RET)) then
           call WriteFatalError('Approach/retract currently does not support approach to a specified setpoint.  please choose approach over a specified Z range instead')
           !better to do nothing than to do the wrong thing
        end if

        if (.not. wantHH) numHH = 0       
        if (operating_mode == SCAN) then
           if (z_feedback_choice == MEAN_DEFL_Z) then
              if (.not. wantHH) then
                 numHH = 1
                 wantHH = .true.
              else
                 numHH = numHH+1 !including zero
              end if
           end if
        end if

        if ((fexcite == NO_EXC) .or. (fexcite == ACOUSTIC_PEAK))  then
           Want_NumImpact = .false.
        else
           Want_NumImpact = .true.  !fixme make user input
        end if

        useLockin = (LockInOrder > 0)

	omegai = omegai_khz * 2.0*pi*1.0d3 ! internal value is rad/s.

        osc_omega_dim = osc_omega_dim * 2d0 * pi * 1e3
  
        call assert(omegai(1) > 0, 'omegai(1) > 0')
        
        if (operating_mode == SCAN) then
           if ( z_feedback_choice == MEAN_DEFL_Z) then
              if ( omegad_dim(1) .GT. omegai(1)*1.0d-1 ) then
                 call WriteFatalError('The drive frequency should be at most 10 percent of the cantilever natural frequency')
              end if
           end if
        end if
 
       call calcCantileverData(AutoCalcAlpha, numModes, mtip, alpha, beta, mu, fexcite, want_Custom_B, B, B_input, AutoCalcOmega, omegai, CalcInputK, Keq, invModalMass)
 
        call CalculateChi( Chi, output_type )	
	
      Vscan_dim = SubsLen_dim * LineSpeed_dim ! nm/s

        !this allows being more intelligent about number of points for Fz curves
	if ((fexcite == NO_EXC)) then
   	    omegad_dim(1) = omegai(1)
	end if 

      
!-------------------------------------------------------------


      if (operating_mode == FREQSWEEP) then
         call set_omega_scale( 0.5*(omegad_start+omegad_stop) )
      elseif (fexcite == NO_EXC) then
         if ( wantSampleExc) then
            call set_omega_scale( omegas_dim )
         else
            call set_omega_scale( omegai(1) )
         end if
      else
         call set_omega_scale( omegad_dim(1) )
      end if
      
      call set_a_scale( Ainitial_dim(1) )
	
        do i = 1, numModes
           omegai_nd(i) = NonDimenFreq(omegai(i))
        end do

        Ainitial_nd(1) = NonDimenLength( Ainitial_dim(1) ) 
	Ainitial_nd(2) = NonDimenLength( Ainitial_dim(2) ) 
                  
        Asample = NonDimenLength(Asample_dim)

        do i = 1, maxExc
           Abase_input(i) = NonDimenLength(Abase_input(i))
        end do

        if (Z0 < Zf)  AprchS_dim = -AprchS_dim  !reverse direction for retraction curves
        call SetApproachSpeedSetpoint( NonDimenVelocity(AprchS_dim) )

        F_ForceVol = NonDimenForce(F_ForceVol_dim)
        Rtip_nd = NonDimenLength(Rtip_dim)

        
        call NormalizeMatlProps( substrate_props,  Rtip_dim)
        if (( operating_mode == SCAN) .and. WantFProp) then
           !this has to be first, because thats where fts_ve_model gets set
           call NormalizeMatlProps( feature_props, Rtip_dim)
           
           if ((substrate_props%fts_ve_model == VEM_ATTARD_BASED) .or. (feature_props%fts_ve_model == VEM_ATTARD_BASED)) then
              !these restrictions are so that the number of equations in DDASKR does not change in the middle
              !could just change the GUI such that use doesn't have the option to pick differently? for now keep this
              if (substrate_props%fts_ve_model /= feature_props%fts_ve_model) then
                 call WriteFatalError('If attard model is used, it must be used for both feature and substrate.  Cannot switch.  Change both to attard model, or neither.')
              elseif (substrate_props%N_attard_fourier /= feature_props%N_attard_fourier) then
                 call WriteFatalError('If attard model is used, the number of fourier coefficents must be the same for both substrate and feature.  Change them to be the same')
              elseif ((substrate_props%N_attard_fourier==0) .and. (substrate_props%N_attard_spatial /= feature_props%N_attard_spatial)) then
                 call WriteFatalError('If spatial attard model is used, the number of spatial points must be the same for both substrate and feature.  Change them to be the same')
              end if
                 
           end if
           
        end if


        if (numincycle_dropdown == 0) then
           !expert mode, direct entry
           numincycle = numincycle_direct
        else
           !auto mode.  make some intelligent choices.  "standard speed/accuracy" = 750
           if (isFzCurves(operating_mode, fexcite)) then
              !10 points per natural period for single eigenmode 
              !For multiple eigenmode simulations, at least 3 points per highest frequency.
              !For viscoelasticity, at least a few points per relaxation time.
       
              numincycle = max( 10d0, 3d0 * omegai(numModes)/omegai(1) ) * ( numincycle_dropdown / 750d0)
              numincycle_ve = ceiling(2d0 * pi/ ViscoelasticModelsMaxDeltaT(substrate_props)) !feel like we shouldnt need 2 pi here but we do

              if (numincycle_ve > (2 * numincycle_dropdown)) then
                 write(tmpStr, *) "Viscoelastic relaxation time is very short. This will require a very small time step, resulting in a potentially very long simulation time. Check that you have specified material properties correctly. &
Remember that relaxation times more than 1 - 2 decades away from the contact time will not be physically meaningful. &
If you are sure you want to run this simulation as-is, then you can override this error message by picking 'expert mode' for 'Accuracy vs speed tradeoff' on the simulation tab, and setting 'Deflection points per cycle' to at least ", (numincycle_ve+1) 
                 call WriteFatalError(tmpStr)
              else
                 numincycle = max( numincycle, numincycle_ve)
              end if
                 

           elseif (modulation_type == PEAK_FORCE) then
              !at least 200 points per base motion cycle, at least 10 points per 1st mode cycle, at least 3 points per highest mode cycle
              numincycle = max( 200,  ceiling(10 * omegai(1)/omegad(t,1)), ceiling(3 * omegai(numModes)/omegad(t,1)) ) * ( numincycle_dropdown / 750d0)

              numincycle_ve = ceiling(2d0 * pi/ ViscoelasticModelsMaxDeltaT(substrate_props)) !feel like we shouldnt need 2 pi here but we do

              if (numincycle_ve > (1.5 * numincycle_dropdown)) then
                 write(tmpStr, *) "Viscoelastic relaxation time is very short. This will require a very small time step, resulting in a potentially very long simulation time. Check that you have specified material properties correctly. &
Remember that relaxation times more than 1 - 2 decades away from the contact time will not be physically meaningful. &
If you are sure you want to run this simulation as-is, then you can override this error message by picking 'expert mode' for 'Accuracy vs speed tradeoff' on the simulation tab, and setting 'Deflection points per cycle' to at least ", (numincycle_ve+1) 
                 call WriteFatalError(tmpStr)
              else
                 numincycle = max( numincycle, numincycle_ve)
              end if                            
           else
              numincycle = numincycle_dropdown              
           end if           
        end if

        
        !drk, 10/19/10.  was missing this 2pi earlier.  that made 
        !it all wrong.  It doesn't seem like normalization is the right place
        !to put it, because omega_scale is good for other time 
        !variable (at least, contact time comes out correctly)
        !probably this 2pi should really be inside the compute lockin
        !subroutine, but we don't want to multiply it out everytime so 
        !here is fine for now
        if (want_AutoCalcTC) then
           if (operating_mode == FREQSWEEP) then
              LockinTC = NonDimenTime( 10d0 / ( 0.5*(omegad_start+omegad_stop)))
           else
              LockInTC =  NonDimenTime(10d0/omegad_dim(1))
           end if
        else
           LockInTC =  NonDimenTime(LockInTC) / (2d0 * pi) ! lock-in time constant.
        end if
        triggerDeflLim = NonDimenLength(triggerDeflLim)

        if (operating_mode == FREQSWEEP) then
           omegad_start = NonDimenFreq(omegad_start)
           omegad_stop  = NonDimenFreq(omegad_stop)
        end if
        
	omegas = NonDimenFreq(omegas_dim)

        osc_omega_nd = NonDimenFreq(osc_omega_dim)
       
        Zbase_amp_nd = NondimenLength(Zbase_amp_dim)

        sweep_time = NonDimenTime(sweep_time)

        HF = NonDimenLength(HF)
        LF = NonDimenLength(LF)
        LF2 = NonDimenLength(LF2)
        SubsLen = NonDimenLength(SubsLen_dim)

        Vscan = NonDimenVelocity(Vscan_dim)

       !non-ideal magnetic
      !phase_slope_dim = degrees / kHz
      !phase_slope_nd  = radians / (  (rad/s) * omega_scale)
      !phase_slope_nd = phase_slope_dim * (pi / 180) *  omega_scale / (2 * pi)
      !effiency slope_dim = % / kHz = (1/100) / 1000 = fraction / (100 kHz)     
      phase_slope_nd = NonDimenTime( (phase_slope_dim / 1e3) * (1d0 / 180d0) / 2d0 )

      efficiency_slope_nd = NonDimenTime( (efficiency_slope_dim / 100e3) / (2 * pi) )

      CurveRate = NonDimenFreq( CurveRate_dim)
      
! done normalizing things
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!force coefficients

      subs_fc = computeForceCoeff(substrate_props)
      if (( operating_mode == SCAN) .and. WantFProp) then 
         ft_fc = computeForceCoeff(feature_props)
      end if


        cur_props => substrate_props	     	
        cur_fC => subs_fC
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


        ! ----------------------------------------------------------------------
        !  setup all ddaskr parameters
        ! ----------------------------------------------------------------------

        !NEQ_cant  = numModes*2
        NEQ_cant  = numModes*DOF
                
        if (cur_props%fts_ve_model == VEM_ATTARD_BASED) then           
           if (cur_props%N_Attard_fourier == 0) then
              NEQ = NEQ_cant + cur_props%N_Attard_spatial
           else
              NEQ = NEQ_cant + cur_props%N_Attard_fourier
           end if           
        else
           NEQ = NEQ_cant
        end if

        EQ_attard = NEQ_cant+1
        
        allocate(YPRIME(NEQ))
        allocate(Y(NEQ))

        if ((Want_Impact) .or. (Want_NumImpact) .or. haveHardHertzViscoelas(substrate_props) .or. (wantFProp .and. haveHardHertzViscoelas(feature_props))) then
           NRT = 1
        else
           NRT = 0!since we get the exact same answer without this, let's save ourselves ~10% computation time by not calculating it
        end if
               
        if (.not. want_autoCalcIC) then
           Y(1:NEQ_cant) = Y_IC(1:NEQ_cant)
        elseif ((fexcite /= NO_EXC) .and. (modulation_type /= PEAK_FORCE)) then        
              y = 0
              y(2)=omegad(0d0, 1) !fixme, this may be off by a phase angle? also doesn't consider higher modes  
        else
           y = 0 !all array elements
        end if

        YPRIME = 0
        isTransientOver = .false.
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!11


        do i = 1, numModes
           ! damping is the normalized damping coefficient=2*zeta*omega_i
           ! for large Q, Q ~= 1/ (2 * zeta).  for small Q, we assume that zeta is measured
           ! from a curve fit and then simply *defining* Q = 1/(2 * zeta).  Really need to be more rigorous about this 

           damping(i) = omegai_nd(i)/Quality(i)

           ! drag is the hydrodynamic drag applied to the Z pizeo motion.  in a previous version, the modal damping was 
           ! applied but that is not technically correct b/c the Z piezo may have much lower frequency content
           drag(i) = omegai_nd(i)**2 * beta(i) * omegai(i) * gamma_drag / Keq(i) 
           if (modulation_type == PEAK_FORCE) then
              drag_pf(i) = drag(i)
           else
              drag_pf(i) = 0
           end if
        end do


        call CalcZRange(operating_mode,modulation_type, ZRange, want_freqswp_sp_aprch, Z0,  Zf, Y, Zbase_Amp_nd)
	
        if (operating_mode == FREQSWEEP) then
           call set_omegad_startstop( omegad_start, omegad_stop )
        else
           call set_omegad( NonDimenFreqArray( omegad_dim, maxexc) ) 
        end if

!       Forcing amplitude calculations
        call CalcExcitAmp(output_type, omegai, fexcite, Quality, Chi, B, Keq, omegai_nd, & 
                     Ainitial_nd, F, Abase, phid, Abase_init, F1_init,numModes, beta, mu, mtip, &
                     want_Abase_direct, Abase_input, F_input, Afluid, mstar_div_m, fm_initial_phase, ZBase_Amp_nd, modulation_type)

        if (exc_choice == SELFEXC ) call InitSelfExc( omegai_nd(1), numincycle, fexcite )


        !controller(s) setup
        if (operating_mode == FIXED) then
           Z_Controller_On = .false.      
        else if  (isOpModeApp(operating_mode)) then           
           Z_Controller_On = .false.
           call triggeredFzModeInit()
           if (modulation_type == FREQUENCY) call SetInitialTransientPhase( TRANS_NUM)
           call SortTimeHistories(Z0 < Zf)
        else if (operating_mode == FREQSWEEP) then           
           sweep_rate = (omegad_stop - omegad_start) / sweep_time
           
           call assert(sweep_time > 0, 'sweep_time > 0')
           call assert(sweep_rate /= 0, 'sweep_rate /= 0')

           call SortTimeHistories(sweep_rate>0)
           Fhist = Ahist * NonDimenFreq(1d3 * 2d0 * Pi) !these are the normalized freq at which we should start outputting time history

           !guess on good controller properties.  may need to be more intelligent about this
           sample_freq_hz = 1e6
           KI = 0.002
           KP = 0.0001

           if (want_freqswp_sp_aprch) then            
              !long range forces require some settling time to avoid making the controller go nuts
              !also have a lock-in filter transient
              call SetInitialTransientPhase( TRANS1_CONT_OFF)
              Z_Controller_On = .false.

              call set_omegad( NonDimenFreqArray( omegai, maxExc) ) !approach at the 1st natural frequency
              call Init_Z_Controller(Z0, ASetPt, AMPLITUDE_Z, HF)               
           else
              call SetInitialTransientPhase( TRANS4_CONT_OFF) !skip right to the last phase
              Z_Controller_On = .false.
           end if

        else if (operating_mode == SCAN) then
           call SortTimeHistories(.true.)
           xHist = Ahist * NonDimenLength(1d0 / 1d9)
         
	   call Init_Scanning
    
           if (((fexcite == NO_EXC) .and. (modulation_type /= FORCE_VOL)) .or.(fexcite == ACOUSTIC_PEAK) ) then
              if (modulation_type == PEAK_FORCE) then

                 call Init_Z_Controller( Z0, peakF_setpoint_dim, MAX_FORCE, HF)

                 !12/22 try this
                 call Compute_Z_Position_PeakForce( cos(theta(t,1)), Z_Controller_On, Force,t )
                 
                 Z_Controller_On = .false.
                 call SetInitialTransientPhase( TRANS1 )
              
              else
                 call Init_Z_Controller(Z0, ASetPt, z_feedback_choice, HF)
                 Z_Controller_On = .true.
                 call SetInitialTransientPhase( TRANS1 )
                 
              end if
              
           elseif (modulation_type == AMPLITUDE)  then 
              call Init_Z_Controller(Z0, ASetPt, z_feedback_choice, HF)
              Z_Controller_On = .false.
              call SetInitialTransientPhase( TRANS0 )

              
           elseif (modulation_type == FREQUENCY) then
              if (z_feedback_choice == FREQ_SHIFT_Z) then
                 call Init_Z_Controller(Z0, NonDimenFreq(freq_shift_sp), z_feedback_choice, HF) !this has to be after calcexcitamp
              elseif (z_feedback_choice == PHASE_Z) then
                 call Init_Z_Controller(Z0, freq_shift_sp*pi/180d0, z_feedback_choice, HF)
              else
                 call Init_Z_Controller(Z0, freq_shift_sp, z_feedback_choice, HF)
              end if

              Z_Controller_On = .false.
              call SetInitialTransientPhase( TRANS1)
              if ( .not. FM_want_noncontact) then
                 KP = -KP
                 KI = -KI
              end if
           elseif (modulation_type == FORCE_VOL) then
              call init_ForceVolData(Z0)
              call ForceVolReset()
              Z_Controller_On = .false.
           else
              call assert(.false., 'unhandled case. scan')
           end if

           call assert(Vscan > 0, 'vscan > 0')
        else
           call assert( .false., 'unknown operating mode')
        end if
        if (modulation_type == FREQUENCY) then
           call FreqModInit()
           wanthist_byA = .false. ! no such thing as amplitude ratio for freq modulation
           if (exc_choice == BIMODAL) call WriteFatalError("bimodal FM is not supported yet.")
        end if

        if (want_hydrodynamic_function) then
           call InitHydrodynamicFunction(damping, omegai_nd, keq, fluid_rho, cant_width, cant_len, fluid_eta) !hardwired values for water
           
           !this overwrites previous values.  should be more elegant about this
           do i = 1,numModes
              damping(i) =  hydro_damping( DimenFreq(omegad(t,1)), i)
              omegai_nd(i) =  hydro_omega( DimenFreq(omegad(t,1)), i)
           end do           
        end if

        call ApproachStop()
           

        call CalcTimeCycleInfo(omegai, numModes, output_point_rate, nstep,Zf, pad, &
                   Z0, Quality,Vscan, SubsLen, sample_freq_hz, nstepFZ, triggerDeflLim,&
                   LockInTC, LockInOrder, transient_allowance, Want_Fourier, numofcycle, transient_timeout_in, fexcite)

        call InitBPF( dble(numincycle) / ( 2d0 * pi) )
        
        if ((operating_mode == SCAN) .and. (z_feedback_choice == MEAN_DEFL_Z)) then
           if ((LockInOrder < 3) .and. (LockInTC < NonDimenTime(10d0/omegad_dim(1)))) call WriteFatalError("Lock in time constant should be high enough so as to completely filter out the relatively high Mean deflection signal. Either increase the lockin time constant or increase the filter order (to have a higher roll-off)")
        else
           if ((LockInOrder > 0) .and. (LockInTC < (1/omegai_nd(1)))) call WriteFatalError(  "Lock In time constant must be such that filter bandwidth is less than the cantilever first natural frequency.  Check inputs and re-run.")
           
           if (modulation_type == FREQUENCY .and.(wantCalcPLLGains .or. wantCalcAmpGains)) then
              call CalculateFMGains(LockinOrder, LockinTc, Quality(1), omegai_nd(1), fexcite, osc_Q, osc_omega_nd, .true.)
           end if
        end if

        if ( operating_mode == APPROACH_STEP) call SetApproachSpeedSetpoint( 0d0) !fixme, a bit of a hack        
        !       ----------------------------------------------------------------------
!  setup Output plots
! ----------------------------------------------------------------------

        if (want_hydrodynamic_function) call Setup_hydro_Plots

        if (operating_mode == SCAN) then
           call SetupScanPlots(fexcite, z_feedback_choice)
           if (modulation_type == PEAK_FORCE) then
              call SetupPeakForcePlots()
           end if
           if ( Want_ContTrans) then
              if (fexcite == NO_EXC ) then
                 if (wantSampleExc) then
                    call SetupDebugForceModScan_Transient()
                 else
                    call SetupDebugContactScan_Transient()
                 end if
              elseif (modulation_type == AMPLITUDE) then
                 call SetupDebugAMScan_Transient()      
              elseif (modulation_type == PEAK_FORCE) then
                 call SetupDebugPeakForceScan_Transient()
              end if
           end if
        end if

        if ((operating_mode == FREQSWEEP) .and. want_freqswp_sp_aprch) call SetupDebugAMScan_Transient()

        if (isFzCurves(operating_mode, fexcite)) then
           call SetupFZcurvesNew(numModes, want_FzEig, substrate_props%fts_ve_model == VEM_ATTARD_BASED)
        end if
        
        call SetupMainOutputPlots(xchoice,exc_choice,Want_A2,Want_P2,Want_AZ,Want_P1,Want_MF,Want_PF,Want_ED,Want_EP,Want_I,Want_CT,Want_Fourier, Want_NumImpact, Want_E_Anc, Want_ForceFourier, want_virial, useLockIn, numModes, operating_mode, fexcite, want_ev, want_RMS, xlabel, xunits, scatter, .false.)

        if (WantHH) call SetupHHPlots(numHH,NHH,exc_choice, xlabel, xunits, Want_Fourier, useLockin, scatter)
        call SetupPoincPlots(numpoinc,Npoinc_X,Npoinc_Y,Want_Strob,Want_Impact)
        
        if ( modulation_type == FREQUENCY) then
           call SetupFMPlots( xchoice, fm_direct_control, want_NormFreqShift)
           if (Want_ContTrans) call SetupDebugFM_Transient()
        end if
        
! ----------------------------------------------------------------------
!  Initial conditions 
! ----------------------------------------------------------------------     
        IOUT = 0
        T = 0.0d0      
        TOUT = 0.0d0

        call setup_DDASKR_parameters(NEQ, DELT, LRW, LIW, INFO, RWORK, ATOL, RTOL, IWORK)
   
        call InitializeHystereticModels(Z0)
        call CheckViscoelasticModels(substrate_props, DELT)
        if (( operating_mode == SCAN) .and. WantFProp)  call CheckViscoelasticModels(feature_props, DELT)

        call AntiAlias_init()

        FPeakAtt = 0
        FPeakRep = 0
  
        call ResetCycleVariables(num_impacts, a1, b1, a2, b2, an, bn, F1r, F1i, Indent, MeanForce, Ets, Virial, Ebs, Eprop, Edamp,Edrive, Force_old,dpr_old, tc,min_d, numHH,numModes )

        Amp = 10 !need this, or keepGoing() gets confused before first output point.  time histories rewrittent to ignore Amp before the first output point.
	        

      call OutputMiscParameters( numModes, F, phid, Keq(1),  pad, numincycle, cur_props%aDMT, Chi, omegai, B, fm_initial_phase, Z0, Abase, Y, beta, mu, cur_props%Estar, damping, nstep, fm_amp_kp, fm_amp_ki, fm_amp_kd, fm_pll_kp, fm_pll_ki,fm_pll_kd, alpha, output_point_rate, delt, numofcycle, z_feedback_choice, minimum_tau(cur_props), omegad(t,1) )
              
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!1!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!start actual computation loop

      do while (  (.not. isTransientOver) .or.  (keepGoing( pointsSinceTrans(IOUT), nstep-1, t, XPos, SubsLen, GetApproachSpeedCurrent(t), Zrange, Amp, (DimenFreq(omegad(t,1))-omegad_dim(1)) / (2d0 * pi))))
         IOUT = IOUT + 1
         TOUT = TOUT + delt 
         
         IDID = 5 
         ! IDID = 5 indicates "The integration was successfully completed by finding one or more roots of R(T,Y,Y') at T."
         ! if IDID = 5 and t==tout then we got really lucky and got a root exactly where we'd requested an output.
         do while ( ((IDID == 1) .or. (IDID == 5)) .and. (T < TOUT))
            
            if (solver == SOLV_DDASKR) then
               CALL DDASKR (RES1,NEQ,T,Y,YPRIME,TOUT, INFO,RTOL,ATOL,IDID, & 
                     RWORK,LRW,IWORK,LIW,RPAR,IPAR, JDUM ,PSDUM, RT1,NRT,JROOT)
            elseif (solver == SOLV_LOCKIN_HACK) then
               ! if you pick this and then do a freq sweep, you see the freq resp of the lock-in. right now for debugging only
               Y(1) = sin(t)
               IDID = 5
               t = tout
            end if

            !IDID = ab4(t, Y, YPRIME, delt) !unused solvers.  see note below
            !IDID = rk4(t, Y,YPRIME, delt)
            if (IDID == -7) then
               if (sum(abs(Y)) < 1d-150) then
                  !this catches a nasty bug in approach/retract Fz curves for samples with no long range forces.
                  !after leaving contact with the sample, the cantilver motion will eventually decay to exactly zero
                  !analytically.  numerically it just gets *really* small, like 1e-168.  DDASKR still tries to find
                  !a solution though, and won't just let 1e-168 underflow to zero.  It spits out an error saying
                  !that the non-linear solver won't converge.  To fix it, just need to force Y to exactly zero
                  !and then restart.  
                  Y=0
                  YPRIME=0
                  !restart integration
                  info(1) = 0
                  cycle
               end if
            end if

            
            IF (IDID .LT. 0) GO TO 275 ! IDID < 0 indicates an error condition

            
            d   =  computeTipSampleGap(     t, y, NEQ_cant)  
            dpr =  computeTipSampleGapDeriv(t, y, NEQ_cant) !the output is not sensitive to roundoff errors here 
       
            

            call updateHystereticModels_outsideRES1(d, dpr,t, IOUT,  Y(EQ_attard : NEQ) )
            call Update_viscoelasticity_outsideRes1_Hertz( t, d - cur_props%aDMT, dpr, cur_props, Rtip_nd) 
            call Update_viscoelasticity_outsideRes1_linear( t, d - cur_props%aDMT, cur_props) 
            
            if ((IDID == 5) .and. (JROOT(1) == -1)) then
               !this is an intermediate point on the boundary between contact and non-contact
               !jroot(1) == -1 means a downward going point (instead of upwards).
                           
               if (Want_Impact .and. isTransientOver) then
                  call OutputPoincare(numpoinc,"Impact", Y, Npoinc_x, Npoinc_Y, d, dpr,t)   
               end if

               num_impacts = num_impacts + 1

            end if

         end do 

         if (output_type == INTERFEROMETER) then
            u = sum( y(1:NEQ_cant:DOF )) +  Abase(1)*cos(theta(t,1)) + Abase(2)*cos(theta(t,2))
         else
            !u is observed deflection from the photodiode:		
            u= sum( Chi(1:numModes) * y(1:(NEQ_cant):DOF )) ! chi(1) * y(1) + chi(2) * y(3) + chi(3) * y(5) + ...
         end if


         if (isnan(u) .or. ( abs(u) > 1e6 )) then
            call WriteFatalError('Sanity check failed.  Cantilever deflection is more than one million times the initial amplitude.  This usually indicates a non-converging solution.  Please submit a help ticket to the developers to resolve this issue.')
         end if


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!11
!lockin
         
         if ( isTransientOver )  then
            Noise = rand(0)*2*NoiseAmp-NoiseAmp 
         else
            Noise = 0 
         end if


!Previously solved Lock In filter equations inside RES1.  that was fine, 
!but it was very slow and overkill. a butterworth filter is not stiff or 
!non-smooth and we don't really need the same level of accuracy.  solving
!the equations here using forward euler is about 7 times faster than solving
!using DDASKR.  might consider using a slightly higher order method,
!or just directly discretizing these to digital filters, but 
!for now this is probably good.
         if (useLockIn) call ComputeLockin(MainLockin, delt, u, Noise, t, theta(t,1))

         call DoAdditionalLockins( t, delt, u+Noise , numHH, NHH, exc_choice, WantHH, useLockin, theta(t,1), theta(t,2), want_RMS) !move this earlier so that higher harmonics lockin transients settle
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! controller(s) (if applicable). 
        
         if (modulation_type == FREQUENCY) then            

            new_omegad = PLL_controller(MainLockin, IOUT, pad, t, delt)
           
            if ( PLL_on .and. (mod(IOUT, pad)==0)) then
               call setup_next_fm_point(t, new_omegad )
               if (want_hydrodynamic_function) then
                  do i = 1,numModes
                     damping(i) =  hydro_damping( DimenFreq(omegad(t,1)), i)
                     omegai_nd(i) =  hydro_omega( DimenFreq(omegad(t,1)), i)
                  end do
               end if
            end if

            !dont remember why, but we want to keep this outside the if statement below.
            !there's another check inside Amp_controller to see if the controller is on.
            call Amp_controller(MainLockin, IOUT, pad, delt)

            if (mod(IOUT, pad)==0 .and. ( PLL_On .or. Amp_Controller_On)) then  !fixme, combine this with the freq sweep rescaling? 
               if (fexcite == ACOUSTIC_IDEAL) then
                  Abase(1) = Abase_init * Drive_signal          
                  call CalcModalForcesIdealAcousticMode( numModes, Abase(1), omegad(t,1),omegai_nd, Quality,Keq, F(1,:), phid(1,:), beta,mu, mtip, Afluid, mstar_div_m)
               elseif (fexcite == ACOUSTIC_PIEZOR) then
                  call CalcModalForcesAcouticPiezoResMode( numModes, Abase(1), omegad(t,1),omegai_nd, Quality,Keq, F(1,:), phid(1,:), beta,mu, mtip, Abase_init * Drive_signal, Afluid, mstar_div_m)
               elseif (isMagnetic(fexcite)) then
                  call ScaleModalForcesMag(B, Keq, numModes, F(1,:), Abase(1), Phid, Drive_signal, omegad(t,1))
               else
                  call assert(.false., 'unhandled case. fm.')
               end if
            end if

         end if

         if (exc_choice == SELFEXC ) F = CalcModalForcesSelfExc(u, F1_init)
         if (modulation_type == PEAK_FORCE) then
            call runningForceCalcs(IOUT, IOUT_prevOutput, Force,  FPeakAtt,FPeakRep, MeanForce, d, dpr, t, y)
            call Compute_Z_Position_PeakForce( cos(theta(t,1)), Z_Controller_On, Force,t)
         end if

         if (Quality(1) > 1000) then
            debug_skip = 100 + Quality(1) / 10
         else
            debug_skip = 100
         end if
         

         if (Z_Controller_On) then                       
            call compute_Z_error( MainLockin%R, omegad(t,1), Drive_Signal,u,HigherHarmLockin(1)%R, MainLockin%Th, PeakForce )
            call Z_Controller(ApparentFeatHght(), ActualFeatHght(0d0), delt,  IOUT, pad, Z0, modulation_type)         
            if  (modulation_type == PEAK_FORCE) then
               PeakForce = 0
               PeakForceA = 0
               if ((Want_ContTrans) .and. (.not. isTransientOver)) then
                  call DebugPeakForceScan_Transient(t, Z0, DimenForce(PeakForce_output)*1e9)
               end if
            end if
         end if
         

           if ( Want_ContTrans .and. (.not. isTransientOver) .and. (modulation_type == PEAK_FORCE)) then
              call DebugPeakForceScan_FullTransient(t, d,dpr, Y(1), Z0, DimenForce(Force)*1e9)
           end if
            
         if ( Want_ContTrans .and. (.not. isTransientOver) .and. ( mod(IOUT,debug_skip) == 1) .and. (modulation_type == FREQUENCY)) then
            call DebugFM_Transient( DimenTime(t), u, Z0,  omegad(t,1), MainLockin%R, MainLockin%Th, Drive_signal, fm_phase_error, fm_phase_error_dt, fm_amp_error, int_amp_error, theta(t,1), Z_error ) 
         end if
         
         if (((operating_mode == SCAN) .and. Want_ContTrans) .or. ((operating_mode == FREQSWEEP) .and. want_freqswp_sp_aprch)) then
            if ((.not. isTransientOver) .and. ( mod(IOUT,debug_skip) == 1)) then
               if (fexcite == NO_EXC) then
                  if (wantSampleExc) then
                     call DebugForceModScan_Transient( t, Z0, u, HigherHarmLockin(1)%R ) 
                  else
                     call DebugContactScan_Transient( t, Z0, u ) 
                  end if
                elseif (modulation_type == AMPLITUDE) then
                   call DebugAMScan_Transient( t, Z0, MainLockin%R ) 
                end if
             end if
         end if
       
         if (operating_mode == SCAN) then      
            if (modulation_type == FORCE_VOL) then 
               call ComputeXPosForceVol(t)
            else
               call ComputeXPos(t) 
            end if
            if (wantFProp) then
               if ( onSubstrate(Xpos)) then
                  cur_props => substrate_props 
                  cur_fC => subs_fC             
               else 
                  cur_props => feature_props 
                  cur_fC => ft_fC                
               end if
            end if
         end if

         
!end controller
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!11
!!!!
!        if ((debugging) .and. (IOUT > 900000)) isTransientOver = .true.

         if ( .not. isTransientOver) then
            if (mod(IOUT, progress_rate) == 1)  then
               write(6,*) "=RAPPTURE-PROGRESS=>", transientPercent(t ) ," Computing", " transients..."
               flush(6)               
            end if

            call TransientStateMachine(t, Drive_signal,u, Quality(1), omegai_nd(1), fexcite)  !isTransientOver gets changed in here
            
            if (isTransientOver) then
!               if (debugging)  write(6,*) "transients now over!" // char(10)

               !this is a hack.  take it from the lock-in but later compare to fourier integral
               Phase_initial = MainLockin%Th

               call SetupFirstNonTransient( t, IOUT)
               
               theta_prevOutput = theta(t,1)
               T_prevOutput = t 
               IOUT_prevOutput = IOUT
               
               if (operating_mode == SCAN) then
                  call Center_Z_Controller(Z0)
                  if (WantConstZ) then
                     KI = 0
                     KP = 0
                  end if
               end if
               
            else
               cycle !go back to the top of the loop, skipping all outputs
            end if
         end if
               

! ----------------------------------------------------------------------
!  Calculated outputs
! ---------------------------------------------------------------------- 

         call DoMainFourierIntegrals( t, delt, u+Noise, numHH, NHH, exc_choice, WantHH, Want_A2, Want_p2, Want_Fourier, a1, b1, a2(1), b2(1), an, bn, Force, F1r, F1i, want_ForceFourier, theta(t,1), theta(t,2)) ! don't do this during transients (but lockins yes)

        if (d < min_d) min_d = d
        
	if (fexcite .ne. NO_EXC) then
        	if ((cur_props%aDMT-d).ge.Indent) then
           		Indent = cur_props%aDMT-d
        	end if
	else
	   	if (cur_props%aDMT .ge. d) then
           		Indent = cur_props%aDMT-d
        	end if
	end if
	
        tc = CalcContactTime(t, d, tc)

        if (modulation_type /= PEAK_FORCE) then
           !already did this earlier for peakforce mode in order to use it in the controller
           call runningForceCalcs(IOUT, IOUT_prevOutput, Force,  FPeakAtt,FPeakRep, MeanForce, d, dpr, t, y)
        end if

        
        !energy calculcations
        !Ets is trap rule.  others are just rectangle for now.
        Ets = Ets - ( Force*dpr + Force_old * dpr_old) / 2

        !this is a reasonable definition of the virial, but not a useful one.  this includes an offset due to Z position
        !which we don't want. 
        !Virial = Virial - Force *d
        !this includes only dynamic effects
        Virial = Virial - Force * sum( y(1:NEQ_cant:DOF ))


        basevel = - Abase(1) * omegad(t,1) * sin(theta(t,1)) - Abase(2) * omegad(t,2) * sin(theta(t,2))
        Ebs = Ebs + Force * basevel        

        
        do i=1,numModes
           Eprop(i) = Eprop(i) +  Force*(y(2*i))
           Edamp(i) = Edamp(i) +  UnscaledForce(damping(i) *  y(2*i)**2, i)           
           Edrive(i) = Edrive(i) +  F(1,i) * y(2*i) * cos(theta(t,1)+phid(1,i)) + F(2,i) * y(2*i) * cos(theta(t,2)+phid(2,i)) 
        end do
        
        Force_old = Force
        dpr_old = dpr       
        

        IF (Wanthist) THEN

           if ( .not. timeHC) then
              !if not currently outputting a time history, check to see if we should start the next time history output yet
              timeHC = startTimeHistory(want_Fourier, IOUT, output_point_rate, GetApproachSpeedCurrent(t), Amp, MainLockin%R, t, sweep_rate, XPos)              
           else
              !currently outputting a time history
              IF (Checkh .EQV. .false. ) THEN
                 !do this code once per time history output
                 Ihist = IOUT
                 theta_star = theta(t,1)
                 t_star = t
                 NumRes1 = 0
                 Checkh = .true.
                 call SetupTimeHistoryPlots(Ahist(Amp_index), XPos, operating_mode, wantHist_byA, hist_cycles,numHist, numincycle, exc_choice, substrate_props%fts_model == ELECTROSTATIC_NONCONS, want_TH_obs_defl, want_th_force_gap,  substrate_props%fts_ve_model == VEM_ATTARD_BASED , do_fft, haveHertzViscoelas(cur_props), want_FZ, want_acceleration, .true., .false.)
              END IF
              
              IF ((IOUT .LT. (Ihist+hist_cycles*numincycle)) .AND. (Amp_index.LE.numHist)) THEN                 
                 if ( modulation_type == FREQUENCY) then
                    tau = t - t_star
                 else
                    tau = (theta(t,1) - theta_star) / (2 * pi)
                 end if
                 
                 if (modulation_type == FORCE_VOL) then
                    if ((ForceVolState == FORCEVOL_APPR) .or. (ForceVolState == FORCEVOL_SETTLE) .or. (ForceVolState == FORCEVOL_REVERSE)) then
                       call OutputTimeHistory_buffered( tau, u, dpr, Force,d, y, F(1,1), IOUT - Ihist+1, cur_props, current_contact_radius, computeZdist_dim(t) , yprime, u0, NumRes1, ZTrack, .false. )
                    end if
                 else
                    call OutputTimeHistory_buffered( tau, u, dpr, Force,d, y, F(1,1), IOUT - Ihist+1, cur_props, current_contact_radius, computeZdist_dim(t), yprime, u0, NumRes1, ZTrack, .false. )     
                 end if

                 NumRes1 = 0
              ELSE
                 timeHC = .false.
                 Checkh = .false.
                 Amp_index = Amp_index + 1
              END IF
           end if ! timeHC
        END IF !wantHist

        
        if ((mod( INT(pointsSinceTrans(IOUT)),numincycle) .EQ. 1).and.( Want_Strob)) then
           call OutputPoincare(numpoinc, "Strob", Y, Npoinc_X, Npoinc_Y, d, dpr, t )
        end if
        
        if ((fexcite == NO_EXC).or. (fexcite == ACOUSTIC_PEAK)) then
           call AntiAlias_put( u,     AA_FZ_U)
           call AntiAlias_put( Force, AA_FZ_F)
           call AntiAlias_put( d,     AA_FZ_d)

           if (isFzCurves(operating_mode, fexcite)) then
              !fixme should anti-alias dpr but its mainly for debugging anywa
              !fixme, there's an off by one error in the antialias filter
              if ((mod(pointsSinceTrans(IOUT),nstepFZ).EQ.0) .AND. (pointsSinceTrans(IOUT).ne.0)) then
                 call OutputFzCurvesNew(AntiAlias_get(AA_FZ_U),AntiAlias_get(AA_FZ_F) ,AntiAlias_get(AA_FZ_d),dpr, computeZDist_dim(t), GetApproachSpeedCurrent(t), triggerState, DimenTime(timeSinceTrans(t)) , y, want_FzEig, Indent, substrate_props%fts_ve_model == VEM_ATTARD_BASED, u0, iout, cur_props )

                 
              end if
           end if
        end if
        
        if (modulation_type == FREQUENCY) then
           if (Want_NormFreqShift) then
              call AntiAlias_put( (DimenFreq(omegad(t,1)) -  omegad_dim(1)) /omegad_dim(1), AA_FREQ_OUT)
           else
              call AntiAlias_put(  (DimenFreq(omegad(t,1)) -  omegad_dim(1)) / (2d0 * pi) / 1e3, AA_FREQ_OUT)
           end if
           call AntiAlias_put( Drive_signal, AA_DRV_OUT)
        end if

        call AntiAlias_put( Ztrack, AA_Z_OUT)
        
        if (OutputThisStep( t,  output_point_rate, IOUT,  IOUT_prevOutput, theta(t,1), theta(t-delt,1), nstep, Z_Controller_On) ) then           
           avg_impacts_per_cycle = num_impacts / ( (theta(t,1) - theta_prevOutput) / (2d0 * pi) )

           MeanForce = Meanforce / (IOUT - IOUT_prevOutput)

           tc = tc /( (theta(t,1) - theta_prevOutput) / (2d0 * pi) )

           if (tc == 0d0) then
              !contact time of 0 for permanent contact really confuses people
              !old way, was incorrect for freq sweep and freq mod
!              if ( d < cur_props%aDMT ) tc = 2d0 * pi
              !new way, this is good for freq sweep and freq mod
               if ( d < cur_props%aDMT ) tc = (T - T_prevOutput) / ( (theta(t,1) - theta_prevOutput) / (2d0 * pi))
           end if
                           
           !old way was slightly wrong for freq sweep b/c delt is not constant (although it's close to 
           !constant over a given output point).  what we need to calculate is the total amount of time
           !that has passed since the last output point.  
           !Amp = 2/(numincycle*datacycle*delt)*sqrt(a1**2+b1**2) ! old way
           Amp = 2/(T - T_prevOutput)*sqrt(a1**2+b1**2) 
           
           Phase = atan2(b1,a1)	      
           
           !the 1d0 is ND stiffness Keq(1) / Keq(1). 2*Pi is needed to match some non-dimensionalization but I've lost track of what.
!           E_anc = (2 * pi) * 0.5 * ( 1d0 ) * omegai_nd(1) * ( Ainitial_nd(1) * Amp * sin(Pi/2 - MainLockin%Th + Phase_initial) - Amp ** 2) / Quality(1)
!           E_anc = (2 * pi) * 0.5 * ( 1d0 ) * omegai_nd(1) * ( Ainitial_nd(1) * Amp * sin(Phase) - Amp ** 2) / Quality(1)
           !5/7/2011.  originally only implemented the on resonance formula.  
!           E_anc = DimenEnergy((2d0 * pi) *0.5* ( 1d0 ) * omegai_nd(1) * ( Ainitial_nd(1) * Amp * sin(Phase) - (Amp ** 2) * omegad(t,1) ) / Quality(1))

           E_anc = 2d0 * pi * 0.5*  Keq(1)  * ( Ainitial_dim(1) * DimenLength(Amp) * sin(Phase) - (DimenLength(Amp) ** 2) * DimenFreq(omegad(t,1)) / omegai(1) ) / Quality(1)  ! J/drive cycle.  without the 2pi it's J/rad
           
           !Double check that these all got scaled right! Ets and Eprop used to be scaled by mode  1 
           ! and drive and damp used to be be i.  Due to the way that tip sample interaction force is calculated
           ! these are all scaled to J/drive cycle

           ! 9/7/2014: note that there is no factor of 2*pi on the virial.  Check the VEDA manual. Virial is defined slightly
           !differently than Ets.  Virial has a factor of (1/T) in front, but dissipation does not.  Basically, dissipation
           !is the *total* dissipation over the drive cycle, whereas virial is more like the *average* conservative work over the
           !cycle, not the sum of all of the conservative work over the cycle.  there was an erroneous factor of 2 * pi here  
           !previously.
           Virial = DimenEnergy(Virial)*delt / (theta(t,1) - theta_prevOutput) 
           Ets =    DimenEnergy(Ets)   *delt / ( (theta(t,1) - theta_prevOutput) / (2d0 * pi)) 
           Ebs =    DimenEnergy(Ebs)   *delt / ( (theta(t,1) - theta_prevOutput) / (2d0 * pi))

	   do i=1,numModes
             Edrive(i) = DimenEnergy(Edrive(i))* delt / ( (theta(t,1) - theta_prevOutput) / (2d0 * pi) )
             Edamp(i)  = DimenEnergy(Edamp(i)) * delt / ( (theta(t,1) - theta_prevOutput) / (2d0 * pi) )
             Eprop(i)  = DimenEnergy(Eprop(i)) * delt / ( (theta(t,1) - theta_prevOutput) / (2d0 * pi) )
	   end do

           !new method.  xaxis_1 = the choice for first harmonic first frequency
           !             xaxis_2 = the choice for first harmonic second frequency (bimodal)
           !             xaxis   = the choice for everything else.
           if  ((operating_mode == APPROACH) .or. (operating_mode == APPROACH_STEP) .or. (operating_mode == APPR_RET)) then
              if (xchoice == X_AMP) then
                 if (Want_Fourier) then 
                    xaxis = Amp
                 else
                    xaxis = MainLockin%R
                 end if
                 xaxis_1 = computeZdist_dim(t)
                 xaxis_2 = xaxis
              elseif (xchoice == X_AMP_BIMODAL) then
                 if (Want_Fourier) then 
                    xaxis = 2d0 / (T-T_prevOutput) * sqrt( a2(1) **2 + b2(1) ** 2) / Ainitial_nd(2)
                 else
                    xaxis = BimodalLockin%R / Ainitial_nd(2)
                 end if
                 xaxis_1 = xaxis
                 xaxis_2 = computeZdist_dim(t)
              elseif ( xchoice == X_ZDIST) then
                 xaxis = computeZdist_dim(t)
                 xaxis_1 = xaxis
                 xaxis_2 = xaxis
              elseif ( xchoice == X_MINGAP) then
                 xaxis = DimenLength(min_d  * 1d9)
                 xaxis_1 = xaxis
                 xaxis_2 = xaxis
              end if
           else if (operating_mode == APPROACH_TRIGGERED) then
              xaxis = computeZdist_dim(t)
           else if (operating_mode == FREQSWEEP) then
              xaxis = DimenFreq(omegad(t,1))  / (2d0 * pi) / 1e3 !back to kHz
              xaxis_1 = xaxis
              xaxis_2 = xaxis
           else ! scanning
              xaxis = DimenLength(Xpos *1d9)
              xaxis_1 = xaxis
              xaxis_2 = xaxis
           end if
           
           if (Want_ForceFourier) then
              call OutputHigherFrequencies(T-T_prevOutput,1,F1r,F1i,xaxis,xaxis,3,Want_A2,Want_P2,Want_Fourier)
           end if

           if (WantHH) then
              ! this outputs a user selected higher harmonic
              call OutputHigherFrequencies(T-T_prevOutput,numHH,an,bn,xaxis,xaxis,1,Want_A2,Want_P2,Want_Fourier)
              if (useLockin) call OutputHigherHarmonicLockins( T-T_prevOutput, numHH, HigherHarmLockin, xaxis)
           end if

           if (exc_choice == BIMODAL ) then
              call OutputBimodalLockin( BimodalLockin, xaxis_2, xaxis, want_a2, want_p2)
           end if           

           if ((operating_mode /= FREQSWEEP) .and. (exc_choice == BIMODAL )) then
              !this output the 2nd excitation frequency in bimodal schemes
              call OutputHigherFrequencies( T - T_prevOutput, 1, a2(1), b2(1), xaxis_2, xaxis, 2,Want_A2,Want_P2,Want_Fourier)
           end if
           

           if (modulation_type == FREQUENCY) then
              call OutputFMPlots(  xaxis, AntiAlias_get(AA_FREQ_OUT) ,fm_phase_error, AntiAlias_get(AA_DRV_OUT))
           end if
           
           call  OutputMainResults(fexcite, avg_impacts_per_cycle, Want_NumImpact, FPeakAtt,Keq,&
                FPeakRep, MeanForce, Phase, Ets, Indent, tc, omegai_nd, Amp,  & 
                modulation_type, Want_AZ,Want_P1,Want_MF,Want_PF,Want_ED,Want_EP,Want_I,Want_CT, &
                Want_Fourier, Want_E_Anc, want_virial, virial, Eprop, Edamp,Edrive, Ebs, E_Anc, numModes, MainLockin%R, &
                MainLockin%Th, RMS%R, want_RMS, useLockIn, omegad(t,1), want_ev, xaxis, xaxis_1, .false.)

           if (operating_mode == SCAN) then  
              if (modulation_type == PEAK_FORCE) then
                 call OutputScanData( xaxis, antiAlias_get(AA_Z_OUT), Amp, ASetPt, ActualFeatHght(0d0), ErrorZ, fexcite, u, z_feedback_choice)
                 call OutputPeakForceData( xaxis, PeakForce_output, PeakForceA_output )
              elseif (modulation_type /= FORCE_VOL) then
                 call OutputScanData( xaxis, antiAlias_get(AA_Z_OUT)  , Amp, ASetPt, ActualFeatHght(0d0), ErrorZ, fexcite, u, z_feedback_choice)                        
              else 
                 call OutputScanData( xaxis, ActualFeatHght(0d0) + ZForceVol_Output , Amp, ASetPt, ActualFeatHght(0d0), ZForceVol_Output, fexcite, u, z_feedback_choice) 
              end if
           end if

           

!           if (debugging) call OutputDebugPlots( xaxis, DimenLength(relaxed_d) * 1e9)
           if (want_hydrodynamic_function) call Output_hydro_Plots( xaxis, added_mass( DimenFreq(omegad(t,1)) ), added_viscosity( DimenFreq(omegad(t,1)) ), DimenFreq(hydro_omega( DimenFreq(omegad(t,1)),1)/(2e0*pi*1e3)), hydro_damping( DimenFreq(omegad(t,1)),1) )

           !done with outputs

           if (operating_mode ==  FREQSWEEP) then
              !update hydodynamic function for driving frequency changes.
              if (want_hydrodynamic_function) then
                 do i = 1,numModes
                    damping(i) =  hydro_damping( DimenFreq(omegad(t,1)), i)
                    omegai_nd(i) =  hydro_omega( DimenFreq(omegad(t,1)), i)
                 end do
              end if

              if ( sweepchoice == STAIRSTEP) then
                 call setup_next_stairstep_point(t)              

                 !for cont sweep this calc is done inside RES1 (which is annoying)
                 call RescaleModalForcesFreqSweep(fexcite, nummodes, Abase, t, omegai_nd, Quality, Keq, F, phid, beta, mu, mtip,Afluid, mstar_div_m, exc_choice, Abase_init, want_nonideal_magnetic, B)
              end if

           end if
                                 
!  Clear cyclic calculations
           call ResetCycleVariables(num_impacts, a1, b1, a2, b2, an, bn, F1r, F1i, Indent, MeanForce, Ets, Virial, Ebs, Eprop, Edamp,Edrive, Force_old, dpr_old, tc, min_d, numHH,numModes )
           IOUT_prevOutput = IOUT
           T_prevOutput = T
           theta_prevOutput = theta(t,1)
           
           ! Simulation progress bar              
           if ( Zrange == ZRANGE_ASP) then
              if ( operating_mode == APPR_RET) then
                 !i think this is correct, but untested.  computeTipSampleGap() is actually wrong for this case!
                 if ( GetApproachSpeedCurrent(t) < 0 ) then
                    Nprog2 = idint(100d0 *  (1d0 - Amp) / (1d0 - Asp_f) ) / 2
                 else
                    Nprog2 = 100 +  idint(100d0 *  (Amp- 1d0) / (1d0 - Asp_f) ) /2
                 end if
              else
                 Nprog2 = idint(100d0 *  (1d0 - Amp) / (1d0 - Asp_f) )
              end if
           elseif (modulation_type == FORCE_VOL) then
              Nprog2 = idint(1.0d2*Xpos/SubsLen)
           elseif (modulation_type == PEAK_FORCE) then
              Nprog2 = idint(1.0d2*Xpos/SubsLen)
           else
              Nprog2 = idint(1.0d2*timeSinceTrans(t)/SimulationDuration)
           end if
           
           write(6,*) "=RAPPTURE-PROGRESS=>",Nprog2," Performing simulation..."
           flush(6)
              
		
        end if !outputthisstep

! ----------------------------------------------------------------------
!  End of Outputs
        ! ----------------------------------------------------------------------
        if ( modulation_type == FORCE_VOL) then
           call updateForceVol(Force, t)         
        end if

	if ( operating_mode == FREQSWEEP ) call UpdateDelt(t)
        if (( operating_mode == APPR_RET) .or. ( operating_mode == APPROACH_TRIGGERED)) then
           call UpdateTriggeredFzMode(u, t)
        end if
       

 
 end do !iout
!
! PROBLEM COMPLETE.  PRINT FINAL STATISTICS.

275 CONTINUE ! this is where we end up for an error condition above

 call OutputInputEcho( inputEcho)

 if ((wantHist) .and. (Amp_index > 1))  call FlushTimeHistory(do_fft)

 IF ( IDID < 0) then
    if ((cur_props%fts_model == ELECTROSTATIC_XU) .or. (cur_props%fts_model == ELECTROSTATIC_GIL)) then
       d = computeTipSampleGap(t, y, NEQ_cant )                     

       if ( d <= 1e-4) then
          call WriteFatalError('Internal error.  Most likely cause is tip-crash in electrostatic model (i.e. the tip got too close and snapped in due to electrostatic forces).  Possible solutions are to keep the cantilever farther away from the sample (i.e. increase Z), increase the cantilever stiffness, or decrease the electrostatic forces. If that does not &
&help then contact developers for assistance')    
       end if
    end if
    
    !fall back generic error message
    call WriteFatalError( 'Internal error. Please contact developers for assistance (error code : DDASKR failed some tests)')
 end IF
      
 call rp_result(driver)

 STOP
END program DMT_DDASKR
      
! End of main function


! ----------------------------------------------------------------------
! Subroutines and external functions
! ----------------------------------------------------------------------


! ODEs subroutines



! from ddaskr docs: this is a subroutine which you provide to define the residual function G(t,y,y')
! of the differential/algebraic system.
!this function is  sensitive to roundoff error in dpr.  relative changes on the order of 1e-12 can
!make a relative change on the order of 1e-6 in the final output. 
SUBROUTINE RES1(T,Y,YPRIME,CJ,DELTA,IRES,RPAR,IPAR) 
  use data0
  use data1  
  use timeAndCycle
  use params
  use contactModels
  use Approaching, only: GetApproachSpeedCurrent
  use ForcingCalculations
  use Nondimensionalization, only: ScaledForce
  use Viscoelasticity_attard, only: compute_attard_derivative

 IMPLICIT none

  real*8, intent(in) :: Y(NEQ), YPRIME(NEQ), CJ, RPAR, IPAR, t
  real*8, intent(out) :: DELTA(NEQ), IRES

  logical :: fail
  integer  ii
  real*8   :: d, dpr, u, f_ts
  real*8, save :: sin1, cos1, cosp1(maxModes), cosp2(maxModes),last_t

  NumRes1 = NumRes1 + 1 !performance counter
  
  IRES = 0 !avoid compiler warning. means nothing to DDASKR	

  !this needs to be at the begining b/c phid can change here and that gets propagated through
  if ((operating_mode == FREQSWEEP) .and. (sweepchoice == CONTSWEEP)) then
     !for sweep, this has to be recalculated at every time point. note, for stairstep recalculated once per output point, but elsewhere     
     call RescaleModalForcesFreqSweep(fexcite, nummodes, Abase, t, omegai_nd, Quality, Keq, F, phid, beta, mu, mtip,Afluid, mstar_div_m, exc_choice, Abase_init,  want_nonideal_magnetic, B)
  end if

  !res1 frequently gets called twice in a row with the same t value.  since sin&cos are
  !very expensive, don't recompute if we don't need to.
  if ( t /= last_t) then
     if ((isAcoustic(fexcite)) .or. (modulation_type == PEAK_FORCE)) then
        sin1 = sin( theta(t,1))  !these calculations are for the base motion
        cos1 = cos( theta(t,1))
     end if
     do ii=1,numModes
        cosp1(ii) = cos(theta(t,1)+phid(1,ii)) !these calculations are for the modal forces
        if (exc_choice == BIMODAL) cosp2(ii) = cos(theta(t,2)+phid(2,ii))
     end do
  end if
  last_t = t
     	
  d   = computeTipSampleGap(t, y, NEQ_cant , cos1)   
  dpr = computeTipSampleGapDeriv(t,y, NEQ_cant , sin1)
  
  !  dpr = computeTipSampleGapDeriv_yprimeversion(t,y,yprime, NEQ_cant , sin1)
  
  call updateHystereticModels_insideRES1(d, dpr)  

   if (cur_props%fts_ve_model == VEM_ATTARD_BASED) then
      DELTA( EQ_attard:NEQ ) = compute_attard_derivative(d, dpr,t, Y(EQ_attard : NEQ), cur_props) - YPRIME(EQ_attard : NEQ)
   end if

   f_ts = Fts(d,dpr, fail, t, y)
  
   
   if (DOF==2) then 
      do ii = 1, numModes        
         DELTA(2*ii-1)=y(2*ii)-YPRIME(2*ii-1)		

         DELTA(2*ii) = -YPRIME(2*ii) -damping(ii)*y(2*ii) -omegai_nd(ii)**2*y(2*ii-1) + & 
              ScaledForce(f_ts + F(1,ii)* cosp1(ii)+F(2,ii)*cosp2(ii) ,ii)   + & 
              GetApproachSpeedCurrent(t) * drag(ii) + drag_pf(ii)*sin1           
      end do
   elseif (DOF==1) then
      !a hack to do just cx'+kx=F instead of mx''+cx'+kx=F.
      !not entirely sure that this is right.  should it really be omegai^2, or just k?
      do ii = 1, numModes
         DELTA(ii) = -YPRIME(ii)*damping(ii) -omegai_nd(ii)**2*y(ii) + & 
              ScaledForce(f_ts + F(1,ii)* cosp1(ii)+F(2,ii)*cosp2(ii) ,ii)   + & 
              GetApproachSpeedCurrent(t) * drag(ii) + drag_pf(ii)*sin1           
      end do
   end if
  
   
  if (fail) IRES = -1 !couldn't evaluate FTS at this point. tell DDASKR this
  
!hank's stuff
!  if (cur_props%fts_model == ELECTROSTATIC_NONCONS) then
!     DELTA(2*numModes+1) =  -yprime(2*numModes+1) + Fts_nc_es(d,dpr, y(2*numModes+1)) !something  
!  end if

end SUBROUTINE RES1



!Boundary subroutine.  DDASKR will compute an intermediate point at every zero of this function.
!note that the solver does not actually USE the intermediate point for anything.  It just gives
!us back that point, which we use for the impact poincare plots.  If impact poincare is not
!being used, this function has absolutely ZERO effect.
SUBROUTINE RT1 (N, T, Y, YPRIME, NRT, RVAL, RPAR, IPAR)
  use data1, only: NEQ,NEQ_cant, computeTipSampleGap
  use contactModels, only: cur_props

  IMPLICIT none
  INTEGER, intent(in) :: N, NRT, ipar
  real*8, intent(in) :: Y(NEQ), YPRIME(NEQ), T, RPAR
  real*8, intent(out) ::  RVAL(NRT)
  real*8 :: d
    
  d = computeTipSampleGap(t, y, NEQ_cant )
  RVAL(1) = d - cur_props%aDMT 
END SUBROUTINE RT1

    !for debugging
    subroutine assert( x, msg)
      logical :: x
      character (len=*) :: msg
      
      if (x .eqv. .FALSE. ) then
         write(6,*) "assertion failed ", msg
      end if
    end subroutine assert

