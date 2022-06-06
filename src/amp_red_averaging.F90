!-----------------------------------------------------------------------
! Amplitude modulated approach curves implemented by amplitude reduction (Averaging) method
! Daniel Kiracofe (daniel.kiracofe@gmail.com)
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


program amp_red_averaging
  use params
  use rappture_io
  use data0
  use data1
  use contactModels
!  use Viscoelasticity_ting
  use ForcingCalculations
  use TimeAndCycle, only: plotpnts, exc_choice,  operating_mode, numincycle, modulation_type, hist_cycles, set_omegad, omegad, delt, UpdateDelt
  use Nondimensionalization
  use TimeHistory
  use tip_data
  use amp_red
  use bisection
  use Fractions
  
  IMPLICIT none 

  real*8 :: Ainitial_nd_act(maxModes),Ainitial_nd_obs(maxModes) !actual vs observed deflection
  real*8 :: Phase1, Phase2, omegad_dim(maxModes)
  real*8 :: xaxis, xaxis_1, xaxis_2, omegai_khz(maxModes)
  real*8 :: Z1, error1, Z2, error2, Z_new,  Asp, Asp_i, Asp_f, Asp_err
  real*8 :: A2_1, A2_2, Asp2,  omega2_omega1_ratio
  real*8 :: FPeakRep, FPeakAtt, Indent, MeanForce, tc, E_anc
  real*8 :: Keq1_nd, Keq2_nd, Amp_iter, Amp2_iter, obs_defl
  integer ::  xchoice, Asp_n, Asp_ndx,  max_f1_cycles, max_f2_cycles,g,openmp_num_threads
  integer :: output_type , Nprog2, numincycle_direct, numincycle_dropdown
  integer :: i,k, iteration_method, ode_solver,  omega2_omega1_rounding
  integer*8 :: ii, longzero

  integer, parameter :: NONE=1, ITBRENT=2, BAHRAM=3,  ITBRENT_A2ONLY=4 !constants for iteration_method

  real*8 :: dz, dz_0, dzold, sign, signold !for matching Bahram's iteration method
  
  !fixme ForceFourier is being computed but axis looks off? 
  logical :: Want_I, Want_CT, Want_PF, Want_ED, Want_ev, Want_MF, Want_E_Anc, want_virial, want_th_obs_defl, want_th_force_gap,Want_P1, Want_P2 , Want_ForceFourier, Want_AZ, Want_A2, Want_AutoCalcAlpha, Want_CalcInputK,Want_Fourier, do_fft, scatter, trapezoidal_int

  
  character*3000 inputEcho
  character*400 tmp
  character*100 inFile, xlabel, xunits

  type(stack_type) :: time_history_stack
  
  real*8 :: adummy(maxModes)
  real*8 :: dummy
  logical :: ldummy
  integer :: idummy

!this is case 5
!  real*8, parameter, dimension(9) :: debug_Z_nm = [63.125,54.53125,46.640625,38.896484375,31.243896484375,23.743896484375,16.556396484375,9.447021484375,2.435302734375]

  
!  real*8, parameter, dimension(9) :: debug_Z_nm = [    18.898438,    16.825195,    14.764160,    12.687866,    10.687866,     8.747471,     6.821977,     4.915109,     3.002420]

  real*8 :: debug_Z_nm(100),  debug_A2_nm(100)
  
!    real*8, parameter, dimension(9) :: debug_Asp  = [0.9, 0.8, 0.7, 0.6, 0.5, 0.4, 0.3,0.2,0.1]
  
  adummy = 0d0
  longzero = 0
  
  ! read rappture inputs
  call getarg(1,inFile)

  call OpenInputFile(inFile)
  
  InputEcho  = "" 

!	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!	Parameters
!	Operating conditions and cantilever properties
!	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  call ReadOperatingParameter( fexcite, exc_choice,idummy , numModes, &
       operating_mode, Want_CalcInputK, ldummy , modulation_type, Want_AutoCalcAlpha , output_type, InputEcho)

  if (exc_choice == SINGLE) then
     numModes = 1
  else
     numModes = 2
  end if
  
  call ReadOpCondAndCantProps(modulation_type, operating_mode, exc_choice,fexcite,numModes, adummy, &
       dummy, dummy, dummy ,Ainitial_dim, dummy , mtip, ldummy , & 
       idummy , dummy , gamma_drag, ldummy ,adummy ,adummy, &
       ldummy , adummy , Afluid, mstar_div_m, dummy , dummy , dummy , dummy , dummy , wantSampleExc, InputEcho)
  
  if (numModes > maxModes) numModes = maxModes
  
  call ReadAmpRedOpParameters( Asp_i, Asp_f, Asp_n, Asp_tolerance, iteration_method, ode_solver, debug_Z_nm,debug_A2_nm, omega2_omega1_ratio, omega2_omega1_rounding, InputEcho)  

 
  !assuming that the custom_b/non-ideal magnetic, or the acoustic stuff does not really matter here.
  !exc frequency is fixed, and user will specify the free amplitude,
  !and force just is whatever it is to get that force. force doesn't enter into the amplitude reduction formula, been replaced by amp, K, Q.
  !base motion also does not enter into the amp red formula.  
  !that stuff matters a lot more for freq sweeps or FM.
  !
  !also setting autoCalOmega=true, to force only reading one mode (second mode will be calculated based on the ratio)
  call readCantModalProps(  numModes, fexcite, omegai_khz,  Keq, Quality, Want_CalcInputK, output_type, Chi, &
       .true. , Want_AutoCalcAlpha, alpha, osc_Q, osc_omega_dim, dummy , dummy , ldummy , &
       ldummy , adummy , inputEcho)

  omegai(1) = omegai_khz(1) * 2d0 * pi * 1000d0  
     
  call ReadTipData( Rtip_dim, Etip, Poisson_tip, tip_angle, tip_shape, InputEcho)
  call ReadSampleData(  substrate_props, "input.phase(ts)",  InputEcho)

  call ReadSimulationParameters(idummy ,Wanthist, ldummy , ldummy, ldummy , idummy , idummy ,Z0,Zf, idummy ,numincycle_dropdown,numincycle_direct ,numHist, idummy ,xchoice, ldummy , dummy , idummy , dummy , ldummy, dummy , dummy , idummy , openmp_num_threads, InputEcho)

#if defined(OPENMP)
    call OMP_set_num_threads(openmp_num_threads) !this is just for attard right now
#endif    
  
  if (numincycle_dropdown==0) then
     numincycle = numincycle_direct
  else
     numincycle = numincycle_dropdown
  end if
  
  if (Wanthist) then
     call Read_TimeHist( numHist, Ahist, Wanthist_byA, want_TH_obs_defl, want_th_force_gap, do_fft, ldummy )
     !Ahist is an array of the desired setpoint ratios.  make it into a stack (what we probably should have done in the original program)
     call init_stack(time_history_stack)
     call SortTimeHistories(.true.)
     do i = 1,numHist
        call push( Ahist(i), time_history_stack)
     end do
     Amp_index=1
  end if
     
  call ReadPlotChoices(Want_A2,Want_P2,Want_AZ,Want_P1,Want_MF,Want_PF,Want_ED, ldummy ,Want_I,Want_CT, Want_Fourier, ldummy , Want_E_Anc, Want_ForceFourier, want_ev, want_virial, ldummy )


  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!done reading inputs
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  
        ! initialize stuff, non-dimensionsionalize, setup tip sample models
        call calcCantileverData(Want_AutoCalcAlpha, numModes, mtip, alpha, beta, mu, fexcite, .false. , B, adummy , .false. , omegai, Want_CalcInputK, Keq, invModalMass)

        omegad_dim = omegai !excitation is always on original resonance

        
        call CalculateChi( Chi, output_type )	
	
        !warning: the iteration tolerances in Brent's method are somewhat dependent on this scaling
        !i.e. Z changing by 1e-9 would be considered "small" by that method. so need the scaling to make it work
        !but... using this scaling seems to make the rcond of Bahrams fourier matrices to be very bad
        !whereas my version is fine either way, and the spatial discretization doesn't change.
        if (.true.) then 
           call set_omega_scale( omegad_dim(1) )
           call set_a_scale( Ainitial_dim(1) )
        else
           call set_omega_scale( 1d0)
           call set_a_scale( 1d0)
        end if

			
        Ainitial_nd_obs(1) = NonDimenLength( Ainitial_dim(1) ) 
	Ainitial_nd_obs(2) = NonDimenLength( Ainitial_dim(2) ) 
 
        Ainitial_nd_act(1) = NonDimenLength( Ainitial_dim(1) ) / Chi(1)
	Ainitial_nd_act(2) = NonDimenLength( Ainitial_dim(2) ) / Chi(2)

 
!        do i = 1, maxExc
!           Abase_input(i) = NonDimenLength(Abase_input(i))
!        end do

        Rtip_nd = NonDimenLength(Rtip_dim)

        !for bimodal, want omegad(2) / omegad(1) to be a rational number with a small denominator
        !this is such that the drive signal repeats after not too long of a period.
        !round omegad2 here.  has no effect on single frequency exc.
        !this needs to happen before viscoelastic materials are initialized
        if (exc_choice == BIMODAL) then

           if ( omega2_omega1_rounding < 0) then
              !original idea.  this general results in a lower denominator, but its harder to
              !explain to the user what is going on.  this is labeled as "expert mode"
           
!              write(*,*) " omega2_omega1_rounding, ", omega2_omega1_rounding, omega2_omega1_ratio
              call fraction_limit_denominator(omega2_omega1_ratio, -omega2_omega1_rounding, max_f2_cycles, max_f1_cycles)
!              write(*,*) "max_f2_cycles f1 ",  max_f2_cycles, max_f1_cycles
           else
              !basic mode. easier to explain to user
              g = gcd( int( omega2_omega1_rounding * omega2_omega1_ratio + 0.5),  omega2_omega1_rounding)
              max_f2_cycles =  int( omega2_omega1_rounding * omega2_omega1_ratio + 0.5) / g
              max_f1_cycles = omega2_omega1_rounding   / g
           end if
           
           omegad_dim(2) = omegad_dim(1) * dble(max_f2_cycles) / dble(max_f1_cycles)
           !also correct omegai to match
           omegai(2) = omegad_dim(2)

           hist_cycles = max_f1_cycles
        else
           hist_cycles = 1
           max_f1_cycles = 1
           max_f2_cycles = 0
        end if
        
        do i = 1, numModes
           omegai_nd(i) = NonDimenFreq(omegai(i))
        end do
        

        do i = 1, numModes
           ! damping is the normalized damping coefficient=2*zeta*omega_i
           ! for large Q, Q ~= 1/ (2 * zeta).  for small Q, we assume that zeta is measured
           ! from a curve fit and then simply *defining* Q = 1/(2 * zeta).

           !Arvind didn't like this.  always allow even for low Q.
           ! if (Quality(i) <= 10) then
           !    call WriteFatalError('Error: To use this tool, Q must be at least 10 (>30 recommended).  For lower Q, use the Basic or Advanced AM Dynamic Approach Curves tools')
           ! end if
           damping(i) = omegai_nd(i)/Quality(i)
        end do

!        write(*,*) "omegad_dim ", omegad_dim
        
        call set_omegad( NonDimenFreqArray( omegad_dim, maxexc) )

!        write(*,*) "omegad ", omegad(0d0, 1), omegad(0d0, 2)

        call UpdateDelt(0d0)

        call NormalizeMatlProps( substrate_props,  Rtip_dim)
        
        call CheckViscoelasticModels(substrate_props, delt)

        !only applicable to this tool, so not included above
        if ( (3d0 * substrate_props%threeelm_tau_stress) > (1d0/ omegad(0d0, 1))  ) then
write(tmp,*)"Relaxation time tau=",DimenTime(substrate_props%threeelm_tau_stress)," is too large relative to the period of the excitation ", DimenTime(1d0/ omegad(0d0, 1)),". A key assumption of this tool is that the surface has time to completely relax between each tap. Either decrease the relaxation time (decrease eta) or decrease the excitation freq"
           call WriteFatalError(tmp)
        end if

       
        subs_fc = computeForceCoeff(substrate_props)
        cur_props => substrate_props	     	
        cur_fC => subs_fC


        !calculates drive forces, which are irrelevant here.
!        call CalcExcitAmp(output_type, omegai, fexcite, Quality, Chi, B, Keq, omegai_nd, & 
!                     Ainitial_nd, F, Abase, phid, Abase_init, F1_init,numModes, beta, mu, mtip, &
!                     .false. , Abase_input, adummy , Afluid, mstar_div_m, dummy, 0d0 , modulation_type)

        !the 2*numModes is backwards compatability to original VEDA which solved Attard's ODEs
        !with DDASKR. somewhat inconvenient but not ready to strip that out just yet.
        if ( cur_props%fts_ve_model == VEM_ATTARD_BASED) then
           if (cur_props%N_Attard_fourier == 0 ) then
              NEQ = numModes*2 + cur_props%N_Attard_spatial
           else
              NEQ = numModes*2 + cur_props%N_Attard_fourier
           end if
        else              
           NEQ = 0
        end if

        !write(*,*) "numModes NEQ ", NEQ
        
        call init_amp_red(ode_solver, max_f1_cycles, max_f2_cycles)
        
        call set_ainitial( Ainitial_nd_act )

        if ((iteration_method == ITBRENT_A2ONLY) .and. (exc_choice == SINGLE)) then
           iteration_method = NONE
        end if
        
        if (iteration_method == ITBRENT) then
           !could be more intelligent here, but should work
           Z1 = 1.05  * Ainitial_nd_act(1)
           Z2 = Asp_i * Ainitial_nd_act(1)
           error2 = nan()
        elseif (iteration_method == ITBRENT_A2ONLY) then
           !fixme initial guess calculated not hardcoded
           error2 = nan()
        else
           dz_0 = NonDimenLength(2.5d-9)
           dz   = dz_0
           sign=1d0
           Z_new =  Asp_i * Ainitial_nd_act(1) + dz_0
        end if

        cumulative_iteration = 1

        !removing this check per 11/4/2020 discussion
        ! if (exc_choice == BIMODAL) then
        !    if ( ( Keq(2) * Ainitial_nd_act(2)**2) >  ( Keq(1) * Ainitial_nd_act(1)**2) ) then
        !       call WriteFatalError('Amplitude at second drive frequency is too large relative to the amplitude at the first drive frequency.  A key assumption of this method is that the second frequency response is a small perturbation to the first frequency response.  Please lower A2')
        !    end if
        ! end if
        
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! done normalizing / initializing things.
! setup output plots
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

        call SetupMainOutputPlots(xchoice,exc_choice,Want_A2,Want_P2,Want_AZ,Want_P1,Want_MF,Want_PF,Want_ED , .false. ,Want_I,Want_CT, Want_Fourier, .false. , Want_E_Anc, Want_ForceFourier, want_virial, .false. , numModes, operating_mode, fexcite, want_ev, .false. , xlabel, xunits, scatter, .true.)

        call SetupIterationStatsPlots(exc_choice)
      
        
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! done with plot set.  Begin main computation loop.
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        
      ! loop over amplitude setpoints
      !if (Asp_i < Asp_f) call WriteFatalError('Final amplitude setpoint must be less than initial amplitude setpoint (i.e. only approach curves, cannot simulate retract)')

      if ((Asp_i == Asp_f) .and. (Asp_n > 1) ) call WriteFatalError('Final amplitude setpoint must be different than initial setpoint (unless number of points = 1)')
      
      do Asp_ndx = 1, Asp_n
         if (Asp_n == 1) then
            Asp = Asp_i
            delta_Asp = 0.1 !for sequential search we use the spacing between points, but for only 1 point there is no spacing
         else
            delta_Asp = ABS((Asp_i - Asp_f) / (Asp_n-1))
            Asp = Asp_i - delta_Asp * (Asp_ndx-1) 
         end if

         !write(*,*) "starting Asp: ", Asp, ' Z1 ', Z1, ' Z2 ', Z2
         call set_Asp1_target(Asp)
         
         !  inner  loop over Z ranges

         if (iteration_method == NONE) then
            !for debugging only, put Z on rails to exactly match Bahram's code
            Z_new = NonDimenLength(debug_Z_nm(Asp_ndx) *1e-9)
            !write(*,*) "Fixed Z(nm): ", debug_Z_nm(Asp_ndx), " A2(nm) ",  debug_A2_nm( Asp_ndx)
            Asp2 = NonDimenLength( debug_A2_nm( Asp_ndx) * 1e-9) / Ainitial_nd_act(2)
            call set_Asp2_target( Asp2 )
            Asp_err = amplitude_reduction_1_error(Z_new)
         elseif (iteration_method == ITBRENT_A2ONLY) then            
            !fixme, variable for initial guesses instead of hard coded
            !fixme, new stats function
            Z_new = NonDimenLength(debug_Z_nm(Asp_ndx) *1e-9)
            call set_Z_target(Z_new)
            !write(*,*) "Fixed Z(nm): ", Z_new
            A2_1 = 1.0
            A2_2 = 0.9
            call sequential_search_brent(amplitude_reduction_2_error , cumulative_iteration, A2_1  , A2_2  , error1, error2, Asp_tolerance, 0.01d0 , 10d0 ,  0d0 ,  OutputIterationStatsPlots_A2 ,  0d0, Asp2, Asp2_err )
         elseif (iteration_method == ITBRENT) then
            if (exc_choice == BIMODAL) then
               call sequential_search_brent(amplitude_reduction_1_error_bimodal , cumulative_iteration, Z1, Z2, error1, error2, Asp_tolerance, delta_Asp *  Ainitial_nd_act(1), Ainitial_nd_act(1) * 100. ,  -Ainitial_nd_act(1) * 10000. ,  OutputIterationStatsPlots_Z,  Asp, Z_new, Asp_err )
               
!               fixme need to retrieve Asp2_err data?
               Asp2 = Amp2 / Ainitial_nd_act(2)
            else                             
               call sequential_search_brent(amplitude_reduction_1_error , cumulative_iteration, Z1, Z2, error1, error2, Asp_tolerance, delta_Asp *  Ainitial_nd_act(1), Ainitial_nd_act(1) * 100. ,  -Ainitial_nd_act(1) * 10000. ,  OutputIterationStatsPlots_Z,  Asp, Z_new, Asp_err )
            end if
        elseif (iteration_method == BAHRAM) then
           if (exc_choice == BIMODAL) call WriteFatalError('Bahrams method not implemented for bimodal')
              
              
           do while(.true.)
              Asp_err = amplitude_reduction_1_error( Z_new)
              call OutputIterationStatsPlots_Z(cumulative_iteration, Z_new, Asp, Asp_err)
              !write(*,*) "bahram iteration: Z ", Z_new, " error", Asp_err
              
              signold=sign
              dzold=dz

              if ( abs(Asp_err)> Asp_tolerance ) then              
                 if (Asp_err > 0) then
                    sign=1
                    if (sign*signold>0) then 
                       dz=dzold
                    elseif (sign*signold<0) then
                       dz=dzold/2
                    end if
                    Z_new =Z_new -dz;
                 else
                    sign=-1
                    if (sign*signold>0) then
                       dz=dzold
                    elseif (sign*signold<0) then
                       dz=dzold/2
                    end if
                    Z_new=Z_new+dz                 
                 end if
              elseif (abs(Asp_err)<Asp_tolerance) then
                 exit
              end if
           
           end do
        else
           call WriteFatalError('Unknown iteration method.')
        end if

        !output the both the target amplitude AND the amplitude from the final iteration.
        !phase calculation is on the target
        !as long as iteration tolerance is tight they are the same.
        !Amp      = Asp          * Ainitial_nd(1)  !now use Amp1 from amp_red module
        
        !fixme but Amp1 is actual deflection, and we want to report observed deflection
        Amp_iter  = (Asp +Asp_err)  * Ainitial_nd_act(1)
        !fixme how to get Asp2_err data
        if (exc_choice == BIMODAL) Amp2_iter = (Asp2+Asp2_err) * Ainitial_nd_act(2)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
           
        Keq1_nd = NonDimenStiffness( Keq(1))
        Keq2_nd = NonDimenStiffness( Keq(2))
        !this is what Bahram had in macromolecules paper
        !           Phase = atan2(  ( 1/Quality(1) + dissipation / pi / Keq_nd / (Amp**2) )    ,  (-2d0 * virial / Keq_nd / (Amp**2)))
        !but if we swap sign on virial, to make it match original VEDA, then this needs to swap as well
        Phase1 = atan2(  ( 1/Quality(1) + dissipation1 / pi / Keq1_nd / (Amp1**2) )   ,  (2d0 * virial1 / Keq1_nd / (Amp1**2)))
        if (exc_choice == BIMODAL) Phase2 = atan2(  ( 1/Quality(2) + dissipation2 / pi / Keq2_nd / (Amp2**2) )   ,  (2d0 * virial2 / Keq2_nd / (Amp2**2)))
        
           !  xaxis_1 = the choice for first harmonic first frequency
           !  xaxis_2 = the choice for first harmonic second frequency (bimodal)
           !  xaxis   = the choice for everything else.
           if (xchoice == X_AMP) then           
              xaxis = Asp          
              xaxis_1 = DimenLength(Z_new*1e9)
              xaxis_2 = xaxis
           elseif (xchoice == X_AMP_BIMODAL) then              
              !              xaxis = BimodalLockin%R / Ainitial_nd(2)              
              xaxis = Asp2
              xaxis_1 = xaxis
              xaxis_2 = DimenLength(Z_new*1e9)
           elseif ( xchoice == X_ZDIST) then 
              xaxis = DimenLength(Z_new * 1e9)
              xaxis_1 = xaxis
              xaxis_2 = xaxis
           end if

           MeanForce = sum(Force) / real(numincycle * max_f1_cycles )
           FPeakRep  = maxval(Force)
!           if (FPeakRep  < 0)  FPeakRep = 0  !actually, this is not how the original VEDA tool worked... "repulsive" can be negative, so it really is just min and max forces
           FPeakAtt  = minval(Force)
!           if (FPeakAtt  > 0)  FPeakAtt = 0
           
           Indent = cur_props%aDMT-minval(gap)
           if (Indent < 0) Indent=0
           
           tc = count( (cur_props%aDMT-gap) >= 0) * delt / max_f1_cycles
           
           call  OutputMainResults(fexcite, 0d0 , .false., FPeakAtt,Keq,&
                FPeakRep, MeanForce, Phase1 ,DimenEnergy( dissipation1+dissipation2* omegad(0d0,2) / omegad(0d0,1)), Indent, tc, omegai_nd, Amp1 * Chi(1),  & 
                modulation_type, Want_AZ,Want_P1,Want_MF,Want_PF,Want_ED, .false. ,Want_I,Want_CT, &
                Want_Fourier, Want_E_Anc, want_virial, DimenEnergy( virial1+virial2) , adummy , adummy, adummy , 0d0, E_Anc, numModes, Amp_iter*Chi(1) , &
                0d0 , 0d0 , .false. , .false. , NonDimenFreq( omegad_dim(1))  , want_ev, xaxis, xaxis_1, .true.)
           
           if (exc_choice == BIMODAL) call OutputBimodal_AmpRed( Amp2 * Chi(2) , Amp2_iter  * Chi(2), Phase2, xaxis_2, xaxis, want_a2, want_p2, DimenEnergy(virial1),  DimenEnergy(virial2),  DimenEnergy(dissipation1),  DimenEnergy(dissipation2* omegad(0d0,2) /  omegad(0d0,1)))
           
           
           ! output time history here if applicable
           IF (Wanthist) THEN

               !check if this setpoint is in the list
              if ( length(time_history_stack) > 0) then
                 if ( Asp <= head(time_history_stack)) then
                    call pop(time_history_stack)
                    Amp_index = Amp_index + 1 !so writefatalerror knows to flush it
                    
                    call SetupTimeHistoryPlots( Asp, 0d0,  operating_mode, wantHist_byA, max_f1_cycles ,numHist, numincycle , exc_choice, substrate_props%fts_model == ELECTROSTATIC_NONCONS, want_TH_obs_defl, want_th_force_gap,substrate_props%fts_ve_model == VEM_ATTARD_BASED , do_fft, haveHertzViscoelas(cur_props), .false. , .false., .false., .false. )

                    !routine was setup to collect data one point at a time. we have it all at once, but want to reuse the existing routine
                    do ii=1,(numincycle*max_f1_cycles)
                       obs_defl = Amp1 * cos( time(ii) *  omegad(0d0, 1) ) * Chi(1) + Amp2 * cos( time(ii) *  omegad(0d0, 2) ) * Chi(2)

                       call OutputTimeHistory_buffered( DimenTime( time(ii)) * ( omegad_dim(1) / 2d0 /pi) , obs_defl , dpr(ii), Force(ii), gap(ii), Y(ii,:) , 0d0, ii, cur_props, current_contact_radius, Z_new, adummy, attard_surf_coord(ii), 0, 0d0, attard_enabled_history(ii) )
                    end do
                    
                 end if                 
              end if
           end IF

           if (iteration_method == ITBRENT) then
              Z1 = Z_new - 0.9 * delta_Asp * Ainitial_nd_act(1) !go slightly less than the next step would imply. don't want to go past
              Z2 = Z_new - 1.1 * delta_Asp * Ainitial_nd_act(1) !go slightly further than the next step would imply.  we always want to go past it
              error2 = nan()
           elseif (iteration_method == ITBRENT_A2ONLY) then
              error2 = nan()
           elseif (iteration_method == BAHRAM) then
              Z_new = Z_new - delta_Asp *  Ainitial_nd_act(1)
              dz=min(dz_0,5*dz)
           end if
           
           Nprog2 = idint( 100.d0 * Asp_ndx / Asp_n )
           write(6,*) "=RAPPTURE-PROGRESS=>",Nprog2," Performing simulation..."			
           flush(6)

        end do

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!        
  ! output data

   call OutputInputEcho( inputEcho)
   
   if (wantHist)   call FlushTimeHistory(do_fft)

   call OutputMiscParameters( numModes, F, phid, Keq(1), longzero, numincycle, cur_props%aDMT, Chi, omegai, B, 0d0, Z0, Abase, adummy, beta, mu, cur_props%Estar, damping, longzero , 0d0, 0d0, 0d0, 0d0, 0d0,0d0, alpha, longzero, 0d0, longzero, 0, minimum_tau(cur_props), omegad(0d0,1) )
   
   call rp_result(driver)

   STOP
end program amp_red_averaging


!Brent's method for root finding. slightly modified from original code at    
    !https://people.sc.fsu.edu/~jburkardt/f_src/brent/brent.f90
    ! biggest different is that we add a termination condition on the function value, in order to make sure
    ! we are always at least as good as Bahram's iteration
    !
    !tried to keep it generic, but the iteration statistics reporting means it is specialized to the amp red avg tool
    !if needed for something else maybe have a call back
    
!! ZERO seeks the root of a function F(X) in an interval [A,B].
!
!  Discussion:
!
!    The interval [A,B] must be a change of sign interval for F.
!    That is, F(A) and F(B) must be of opposite signs.  Then
!    assuming that F is continuous implies the existence of at least
!    one value C between A and B for which F(C) = 0.
!
!    The location of the zero is determined to within an accuracy
!    of 6 * MACHEPS * abs ( C ) + 2 * T.
!
!    Thanks to Thomas Secretin for pointing out a transcription error in the
!    setting of the value of P, 11 February 2013.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    11 February 2013
!
!  Author:
!
!    Original FORTRAN77 version by Richard Brent.
!    FORTRAN90 version by John Burkardt.
!
!  Reference:
!
!    Richard Brent,
!    Algorithms for Minimization Without Derivatives,
!    Dover, 2002,
!    ISBN: 0-486-41998-3,
!    LC: QA402.5.B74.
!
!  Parameters:
!
!    Input, real ( kind = 8 ) A, B, the endpoints of the change of 
!    sign interval.
!
!    input, real*8, fa_in, fb_in, the function value evaluated at the initial endpoints
!
!    Input, real ( kind = 8 ) MACHEP, an estimate for the relative machine
!    precision.
!
!    Input, real ( kind = 8 ) T, a positive error tolerance.
!
!    Input, external real ( kind = 8 ) F, the name of a user-supplied
!    function, of the form "FUNCTION F ( X )", which evaluates the
!    function whose zero is being sought.
!
!    Output, sb the estimated value of a zero of
!    the function F.  fb, value of the function at sb

recursive subroutine brent(f,  a, fa_in, b, fb_in, t, sb, fb, t_f, cumulative_iteration, Asp_for_stat, stat_func )
  use rappture_io
  
  implicit none

  interface
     subroutine stat_func(iteration_number, Z, Asp_target, Asp_error)
         integer, intent(inout) :: iteration_number
         real*8, intent(in) :: Z, Asp_target, Asp_error
     end subroutine stat_func
  end interface

  
  real*8, intent(in) :: a, b, fa_in, fb_in
  real*8, intent(out) :: sb, fb
  integer, intent(inout) :: cumulative_iteration
  real*8, intent(in) :: Asp_for_stat

  real*8 :: c,d,e,f, fa, fc, m, p, q, r, s, sa, t, t_f, tol
  real*8, parameter :: machep = 2.22e-16
  integer :: count

!
!  Make local copies of A and B.
!
  sa = a
  sb = b
  fa = fa_in
  fb = fb_in

  c = sa
  fc = fa
  e = sb - sa
  d = e

  count = 0
  
  
  do

     count = count + 1
     if (count > 50) then
        
        call WriteFatalError('Iteration failed to converge after 50 iterations. This may indicate a bistability or non-monotonic / discontinuous amplitude vs Z. The averaging tool may not be suitable for this simulation. Try the basic AM approach curves')
        return
     end if
     
    if ( abs ( fc ) < abs ( fb ) ) then

      sa = sb
      sb = c
      c = sa
      fa = fb
      fb = fc
      fc = fa

    end if

    tol = 2.0d0 * machep * abs ( sb ) + t
    m = 0.5d0 * ( c - sb )

    call stat_func(cumulative_iteration, sb, Asp_for_stat, fb)
    !write(*,*) "brent search Zl: ", sb, " Zr: ", c, " error@Zl: ", fb

    if ( abs ( m ) <= tol .or. fb == 0.0d0 )  then
       !write(*,*) "brent exiting on |m|<tol"
       exit
    elseif  (abs(fb) < t_f) then
       !write(*,*) "brent exiting on |fb|<t_f"
      exit
    end if


    
    if ( abs ( e ) < tol .or. abs ( fa ) <= abs ( fb ) ) then

      e = m
      d = e

    else

      s = fb / fa

      if ( sa == c ) then

        p = 2.0D+00 * m * s
        q = 1.0D+00 - s

      else

        q = fa / fc
        r = fb / fc
        p = s * ( 2.0D+00 * m * q * ( q - r ) - ( sb - sa ) * ( r - 1.0D+00 ) )
        q = ( q - 1.0D+00 ) * ( r - 1.0D+00 ) * ( s - 1.0D+00 )

      end if

      if ( 0.0D+00 < p ) then
        q = - q
      else
        p = - p
      end if

      s = e
      e = d

      if ( 2.0D+00 * p < 3.0D+00 * m * q - abs ( tol * q ) .and. &
        p < abs ( 0.5D+00 * s * q ) ) then
        d = p / q
      else
        e = m
        d = e
      end if

    end if

    sa = sb
    fa = fb

    if ( tol < abs ( d ) ) then
      sb = sb + d
    else if ( 0.0D+00 < m ) then
      sb = sb + tol
    else
      sb = sb - tol
    end if

    fb = f ( sb )

    if ( ( 0.0D+00 < fb .and. 0.0D+00 < fc ) .or. &
         ( fb <= 0.0D+00 .and. fc <= 0.0D+00 ) ) then
      c = sa
      fc = fa
      e = sb - sa
      d = e
    end if

  end do

!  brent = sb
end subroutine brent


recursive subroutine sequential_search_brent( func, cumulative_iteration, Z1, Z2, error1, error2, tol, search_step, max_Z1_limit, max_Z2_limit, stat_func, target_for_stat,  Z_new, err_new ) 
  integer, parameter :: max_iter = 100
  integer, intent(inout) :: cumulative_iteration
  real*8, intent(inout)  :: Z1, Z2, error1, error2, search_step
  real*8, intent(out) ::  Z_new, err_new
  real*8, intent(in) :: tol , search_step, target_for_stat, max_Z1_limit, max_Z2_limit
  real*8, external :: func
  integer :: k

  interface
     subroutine stat_func(iteration_number, Z, Asp_target, Asp_error)
         integer, intent(inout) :: iteration_number
         real*8, intent(in) :: Z, Asp_target, Asp_error
     end subroutine stat_func
  end interface

  
  !sequential search up to generat initial conditions for bisection. find two Z values that bound the target.
  !Z1 will have positive error, and Z2 will have negative error
  !always start with previous Z for Z1, except for initial point
  !guess a Z, evaluate amplitude reduction. if amp is lower than target, then increase (but we can use that for Z2)

  
  error1 = func( Z1)
  call stat_func(cumulative_iteration, Z1, target_for_stat , error1)  
  !write(*,*) "Z1 ", Z1, " error1 ", error1
  
  k = 1
  do while (( error1 < 0) .and. (abs(error1)> tol)) 
     k=k+1
     
     Z2 = Z1
     error2 = error1
     
     Z1 =  Z1 + (2**k) * search_step
     if  (Z1 > max_Z1_limit ) then
        call WriteFatalError('Seqential search failed. Z1 ridiculously high  ')
     end if
     
     if (k > max_iter) then
        call WriteFatalError('Sequential search failed to find a suitable Z1 value. It is possible that your requested final setpoint is too low (i.e. amplitude simply cannot go that low, even when the tip is permanently in contact with the sample).  Check values.')
     end if
     
     
     error1 = func(Z1)
     call stat_func(cumulative_iteration, Z1, target_for_stat , error1)       
     !write(*,*) "Z1 ", Z1, " error1 ", error1
  end do
  
  if (isnan(error2)) then
     if (abs(error1)> tol) then
        !dont have an estimate for error2 yet, and we will need one
        error2 = func( Z2)
        call stat_func(cumulative_iteration, Z2, target_for_stat , error2)       
     else
        !dont have an estimate for error2 yet, but will not need one, because error1 is already under tolerance
        !just set error2 to be high, in order to short circuit following checks
        error2 = 10. * tol
     end if        
  end if
  

  k = 1
  !write(*,*) error1, error2, tol
  do while (( error2 > 0) .and. ( min(abs(error1), abs(error2)) > tol))
     if (abs(error2) < abs(error1)) then
        !better estimate for Z1
        Z1 = Z2
        error1 = error2
     end if

     !putting the check here allows one trial right at the limit before failure. important for amplitude search where 0 is the min allowable
     if  (Z2 < max_Z2_limit ) then
        call WriteFatalError('Seqential search failed. Z2 ridiculously low. If your sample is very soft (kPa), this might represent a bug in the code ')
     end if

     !each sucessive step go twice as far
     Z2 = max( max_Z2_limit,  Z2 - (2**k) * search_step)
     
     k=k+1
     
     if (k >  max_iter) then
        call WriteFatalError('Sequential search failed to find a suitable Z2 value. It is possible that your requested final setpoint is too low (i.e. amplitude simply cannot go that low, even when the tip is permanently in contact with the sample).  Check values.')
     end if
     
     error2 =  func( Z2)
     call stat_func(cumulative_iteration, Z2, target_for_stat , error2)       
     !write(*,*) "Z2 ", Z2, " error2 ", error2
  end do
  
  !brent's method. usually faster than bisection (almost a factor of 2)
  !brent's method was written in terms of x tolerance, we hacked in a f(x) tolerance.  but it still needs an x tolerance to work with
  !make that one MUCH smaller so that the f(x) tolerance always governs
  if  ( min(abs(error1), abs(error2)) > tol) then
     call brent(func, Z1, error1, Z2, error2, tol/100, Z_new, err_new, tol , cumulative_iteration, target_for_stat,  stat_func )
  elseif (abs(error1) <= tol) then
     err_new = error1
     Z_new = Z1
  elseif (abs(error2) <= tol) then
     err_new = error2
     Z_new = Z2
  else
     call WriteFatalError('cant happen error')
  end if
  
  !write(*,*) "brent search finished Z_new ", Z_new, " error ", err_new
  !write(*,*) "----------------------------------"

end subroutine sequential_search_brent


    !for debugging
    subroutine assert( x, msg)
      logical :: x
      character (len=*) :: msg
      
      if (x .eqv. .FALSE. ) then
         write(6,*) "assertion failed ", msg
      end if
    end subroutine assert
