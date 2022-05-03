!-----------------------------------------------------------------------
!             Copyright 2007 - 2009 by
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
!If VEDA (modified or unmodified), any portion thereof, or any derivative work (publicly distributed or not), contributes to any scientific or academic publication (including, but not limited to, journal articles, conference presentations or posters, or seminars) you must cite the original developers:
!
!Kiracofe, D.; Melcher, J. & Raman, A. "Gaining insight into the physics of dynamic atomic force microscopy in complex environments using the VEDA simulator Review of Scientific Instruments", (2012), 83, 013702
!
!You may also optionally cite the original VEDA article
!
! J. Melcher, S. Hu, A. Raman, "VEDA: A web-based virtual environment for dynamic atomic force microscopy" Review of Scientific Instruments 79, 061301 (2008).
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!basic modules file
!the modules file was getting too big so need to split up.
!due to the way that fortran calculates modules, it is easy if we put in this file only
!modules that have no dependencies (no "use" statements).  we put all the other modules
!in the other file

module params

  !set to 1 for a hack to do just cx'+kx=F instead of mx''+cx'+kx=F.
  !set to 2 normally
  integer, parameter :: DOF=2
  
  character(len=*), parameter :: INPUT_PREFIX = 'input.phase'

  !named constants for tip shape
  integer, parameter :: PARABOLOID = 1, CONE = 2

  !named constants for solver
  integer, parameter :: SOLV_DDASKR = 1, SOLV_RK4 = 2, SOLV_LOCKIN_HACK = 3, SOLV_AB4 = 4, SOLV_FWDEULER = 5

  !named constants for output_type (used to be AutoCalcChi)
  integer, parameter :: AUTO_CALC_CHI = 1, MAN_CALC_CHI = 2, INTERFEROMETER = 3

  !named constant for Zrange
  integer, parameter :: ZRANGE_AUTO = 1, ZRANGE_MANUAL = 2, ZRANGE_ASP = 3, ZRANGE_FREQSHIFT = 4

  !phases of operating for freq sweep w/ approach to setpoint
  !1) TRANS1_CONT_OFF, a few drive cycles at resonance (needed if there are long range forces so that
  !                    the transients from that don't mess up the controller
  !2) TRANS2_CONT_ON, approach to setpoint with controller on at 1st nat freq
  !3) TRANS3_STABILIZE, stay there are make sure we are still within tolerance
  !3) TRANS4_CONT_OFF, a few drive cycles at sweep starting frequency, controller off
  !4) OPERATE    
  !the scanning tool essentially does states 2 -> 4
  !the approach curves essentially does states 3 -> 4  
  integer, parameter :: TRANS1_CONT_OFF = 1, TRANS2_CONT_ON = 2, TRANS3_STABILIZE =3, TRANS4_CONT_OFF = 4, OPERATE = 5

  !phases of operation for AM scanning 
  integer, parameter :: TRANS0 = 0, TRANS1 = 1, TRANS2 = 2

  !phases of operation for FM scanning
  integer, parameter :: TRANS_NUM = 0, TRANS_FTS = 1, TRANS2_PLL=2, TRANS3_AMP=3, TRANS4_Z_ON = 4, TRANS5_Z_STABL = 5, TRANS6_OP = 6
 
  logical, parameter :: debugging = .false. !set this for debugging outputs
  
  !named constants for frequency choice
  integer, parameter ::  SINGLE = 1, BIMODAL = 2, SELFEXC = 3
  
  !named constant for frequency sweep type
  integer, parameter :: CONTSWEEP = 3, STAIRSTEP = 4
  
  !named constants for operating_mode
  ! APPROACH           is a continuous one direction motion (untriggered)
  ! APPR_RET           is a continuous motion first down and then up (untriggered)
  ! APPROACH_STEP      is a set of discrete Z distances (untriggered)
  ! APPROACH_TRIGGERED is what we formerly called "jump" mode, which was just a triggered Fz curve anyway
  ! APPROACH_SINE      is approaching with a sinusoidal motion (untriggered only obviously)
  ! for freqsweep we just have one operating mode and then a sweep type variable.  for approach we
  ! have different operating modes.  in theory we should be consistent but oh well.
  integer, parameter :: APPROACH = 1, FREQSWEEP = 2, SCAN = 3, APPROACH_TRIGGERED = 4, FIXED=5, APPR_RET = 6, APPROACH_STEP = 7, APPROACH_SINE=8
  
  !named constants for modulation_type.  amplitude modulation also includes 
  !f-z curves and frequency sweep (which really don't have a modulation type)
  integer, parameter :: AMPLITUDE = 1, FREQUENCY = 2, FORCE_VOL = 3, PEAK_FORCE = 4

  ! named constants for fexcite
  integer, parameter :: ACOUSTIC_IDEAL = 1, MAGNETIC_LORENTZ = 2, NO_EXC = 3, SAMPLE=4, ACOUSTIC_PIEZOR = 5, MAGNETIC_TORQUE = 7, ACOUSTIC_PEAK = 8
  
  !named constants for feature type (scanning only)
  integer, parameter :: STEP = 1, TRAPEZOID = 2, SINUSOID = 3, CYLINDER = 4

  !named constants for xchoice
  integer, parameter ::  X_AMP = 1, X_ZDIST = 2, X_MINGAP = 3, X_AMP_BIMODAL = 4

  !named constants for z_feedback_choice
  !1 = tapping mode, 2 = standard frequency mode, 3 = drive modulation, 4 = contact mode
  integer, parameter :: AMPLITUDE_Z = 1, FREQ_SHIFT_Z = 2, DRIVE_AMP_Z = 3, DEFLECTION_Z = 4, MEAN_DEFL_Z = 5, PHASE_Z = 6, MAX_FORCE = 7

  !named constants for VEchoice
  integer, parameter :: VEC_NONE = 1, VEC_KELVINVOIGT = 2, VEC_MAXWELL = 3,VEC_KELVINVOIGT_LEGACY = 4, VEC_THREEELM_E1E2ETA = 5, VEC_GENMAXWELL = 6,  VEC_THREEELM_E0_EINF_TAU = 7

  !named constants for fts_ve_model.  
  integer, parameter :: VEM_HERTZ_BASED = 1, VEM_ATTARD_BASED = 2,  VEM_LINEAR_BASED = 3,  VEM_NONE = 4

  !max number of modes. theoretically could use dynamic allocation for all this, but 
  !9 modes is probably good for now.  if you change this, need to adjust output routines.
  !many are hard coded to single character for mode number. 
  integer, parameter :: maxModes=9 
  integer, parameter :: maxHist=9


  real*8, parameter :: pi=4.0d0*atan(1.0d0)
  real*8, parameter :: KB = 1.38062d-23 !boltzmann's constant. J/K
  
  real*8, parameter :: sqrt2 = sqrt(2d0)
  real*8, parameter :: butterworth4_1 = 0.7653668648, butterworth4_2 = 1.847759065

  !max number of driving frequencies (right now just bimodal)
  integer, parameter :: maxExc = 2
  
  !  Permittivity of free space		
  real *8, parameter :: epsilon0 = 8.8542d-12

  ! permeability of free space
  real* 8, parameter :: mu0 = 1.25663706d-6
  
  ! 1 eV = 1.60217653e-19 Joules, and 1 J = 6.24150974e18 eV, which is also the number of electrons per coulomb
  ! so multiply energy in joules by this to get energy in electron volts.
  real *8, parameter :: electrons_per_coulomb = 6.24150974d18 

  contains


  !simple macros for interpreting above choices.  put pure functions only here, no subroutines!
    pure logical function isAcoustic(f)
      integer, intent(in) :: f
      isAcoustic = ( ( f == ACOUSTIC_IDEAL) .or. (f == ACOUSTIC_PIEZOR))
    end function isAcoustic

    pure logical function isMagnetic(f)
      integer, intent(in) :: f
      isMagnetic = ( ( f == MAGNETIC_LORENTZ) .or. (f == MAGNETIC_TORQUE))
    end function isMagnetic
    
    pure logical function isOpModeApp(f)
      integer, intent(in) :: f
      isOpModeApp = ( (f == APPROACH) .or. (f == APPROACH_STEP) .or.  (f == APPR_RET) .or. ( f == APPROACH_TRIGGERED ) .or. ( f == APPROACH_SINE ) )
    end function isOpModeApp

    real*8 function nan() 
      character(4):: cinf(3)
      cinf = (/'NaN ','-Inf','+Inf'/)
      read(cinf,*) nan
    end function nan    
end module params

module matl_prop_type_module
  !named constants for fts_model
  integer, parameter :: LINEAR = 1, HERTZ = 2, DMT = 3, DMT_DLVO = 4, CHADWICK = 5,  LIN_ATT = 6, JKR = 7
  integer, parameter :: MORSE = 8, LENNARD_JONES = 9, ELECTROSTATIC_XU = 10, CUSTOM_CONS = 11
  integer, parameter :: ELECTROSTATIC_GIL = 13, ELECTROSTATIC_NONCONS = 14, MAG_DIPOLE = 16
  integer, parameter :: CUSTOM_MD = 17, ATTARD_FOURIER_BAHRAM = 18, HERTZ_BEC = 19,  ATTARD_FOURIER_LSQ = 20

  !named constants for calcADMT
  integer, parameter :: ENTER_FAD_ONLY = 1, ENTER_FAD_A0 = 2, ENTER_FAD_H = 3, ENTER_H_A0=4

  integer, parameter :: max_prony = 99

  type matl_prop
     integer :: fts_model, fts_ve_model
     
     !hertz, DMT, and misc.
     real*8 :: A_hamaker_dim, A_hamaker_nd, Estar, deltaE, D_0, Fadhesion, aDMT, Poisson_sample, Esample
     integer :: calcADMT

     !surf hysteresis
     logical :: wantSurfHyst
     real*8  :: SurfHyst_gamma
       
     !viscoelasticity
     !the star values are reduced (tip is always assumed rigid for viscoelasticity)
     integer :: VEchoice
     real*8 :: etasample, maxwell_t1, threeelm_e1e2,threeelm_tau_stress, threeelm_tau_stress_dim,threeelm_E2, threeelm_E2star, threeelm_surf_tau, k2
     real*8 :: threeelm_e1, threeelm_e0, threeelm_einf  !used only for attard
     real*8 :: E_gen_max(max_prony), E_gen_max_star(max_prony), tau_gen_max(max_prony), tau_gen_max_dim(max_prony)
     integer :: N_gen_max, N_attard_spatial, N_attard_fourier
     real*8, allocatable :: gen_maxwell_psi_lookup(:),  gen_maxwell_t_lookup(:)
     logical :: input_shear_modulus
     !attard
     real*8 :: attard_radial_extent, tau_creep
     real*8 :: attard_start_nd, attard_stop_nd, attard_start_dim, attard_stop_dim !only used for the amplitude reduction averaging tool
     real*8 ::  drr
     real*8, allocatable :: k(:, :), kt(:,:), tip_shape(:), rr(:),  bare_k(:, :)
     real*8, allocatable :: J(:,:),  b(:,:), C(:,:), CT(:,:), gamma(:,:)
     integer, allocatable :: ipiv(:)
     logical ::  use_bahram_version
     real*8 :: attard_E0, attard_Einfinity !these are reduced nondimensionalized young's modulus.  Input values are non-reduced, and might be shear or youngs, and might be e1,e2 or E0,Einf

     
     !linear
     real*8 :: kts_R, kts_A,Lo
     
     !chadwick
     real*8 :: hs
     
     !DLVO
     real*8 :: sigmas, sigmat, Kd_dim, Kd_nd, epsilon
     
     !electrostatics
     real*8 :: dc_bias_voltage, ac_bias_voltage, bias_freq_rads, surface_pot
     real*8 :: resis, epsilon_0, epsilon_d, g_d_dim, area

     !magnetic dipole
     real*8 :: mom_sample, mom_tip, mag_dia, mag_delta, mag_z0

     !lennard jones + Morse
     real*8 :: MorseLambda_dim, MorseRc_dim,LJ_r0_dim, MorseLambda_nd, MorseRc_nd,LJ_r0_nd, MorseU0,  LJ_E0
     !exp decaying viscosity
     real*8 ::  exp_dashpot_scale, exp_dashpot_decay
     
     !jkr
     real*8 :: acrit_JKR, dJKRcrit, a0_JKR, a_JKR, W_JKR
     
     !hydration + oscillatory forces
     logical wantOscillatory, wantHydration, wantContGrad, want_exp_dashpot
     integer :: N_solv
     
     real*8 :: sigma_solvation_nd, lambda_solvation_nd
     real*8 :: rho  ! water number density. 1/m^3
     real*8 :: p_h  !jeffrey does not give the value of his constant p_h.  I'm guessing 2d7 is good here.  think it should be overall repulsive.  
     real*8 :: sigma_solvation_dim, lambda_solvation_dim, Rtip_solvation_nd, Rtip_solvation_dim

     ! hysteretic solvation (uses preisach, not published)
     logical Want_hyst_hydr
     real*8 :: app_decay, ret_decay, app_scaling, cutoff_dist
     integer :: nside     
     logical wantPreisach
     
     !capillary
     logical wantCapAd
     
     !tip sample squeeze film
     logical :: want_tip_squeeze
     real*8 :: eta_liquid
     
    !worm like chain
     real*8 :: wlc_L0, wlc_Lp, wlc_Lr, wlc_Fr
     logical :: wantWLC
     
     real*8 :: Temp !deg kelvin

     !custom tip sample interaction model, direct input
     !the file for the gap table should be monotonically increasing.  indentations are first and then gaps. gap table should be dimensional (meters and Newtons)      
     integer :: N_custom
     real*8, allocatable :: gap_table_dim(:), F_table_dim(:)         

     !custom tip sample interaction fit from MD data.  
     integer num_ret_coeff, num_x, num_app_coeff
     real*8, allocatable :: indentation_x(:), ret_coeff(:,:), app_coeff(:)
  end type matl_prop

      type forceCoeff
         real*8  CvdW, CDMT,  Dts,     CCHAD,  CDLVO, D_exp
         real*8  Clin, C_att, C_CapAd,  Fad, C_Fjkr1, C_Fjkr2, C_deltJKR1, C_deltJKR2
         real*8  Coscillatory, Chydration, CMorse, CLJ, Celectrostatic_xu, Ces_gil_c1, Ces_gil_c2 
         real*8  Ces_gil_c5, tip_squeeze, kv_t0, vec_f !, Ces_nc_cap1, Ces_nc_dcdq1
         real*8  C_mag_dipole, hertz_bec
         real*8, allocatable :: gap_table_nd(:), F_table_nd(:)
      end type forceCoeff
end module matl_prop_type_module


!fortran is inconvenient about module dependencies.  need these to be separate
!so that Nondimensionalization can depend on this (because EVERYTHING depends
!on Nondimensionalization)
module data0
  use params
  real*8 :: omegai_nd(maxModes),Keq(maxModes), Omegai(maxModes), invModalMass(maxModes)
end module data0


  !this is a very simple antialiasing filter with cutoff freq 1/N
  !essentially an FIR filter with 2*N taps all equal to 1.
  !by all rights we should have an anti-alias filter on any output that is sampled and
  !on everything that has a controller on it.  for now, just have some of the ones where I've
  !noticed a problem.
 module AntiAliasingFilter
   integer, parameter :: naccum = 6
   real*8, private :: accum1(naccum),accum2(naccum)
   integer, private :: N1(naccum), N2(naccum)
   !enums
   integer, parameter :: AA_FREQ_OUT = 1, AA_DRV_OUT = 2, AA_Z_OUT =3, AA_FZ_U = 4, AA_FZ_F = 5, AA_FZ_d = 6
   logical :: firstCall

   contains

     subroutine AntiAlias_init()
       firstCall = .true.
       accum1=0
       accum2=0
       N1=0
       N2=0
     end subroutine AntiAlias_init

   subroutine AntiAlias_put(sample, num)
     real*8, intent(in) :: sample
     integer, intent(in) :: num
     accum1(num) = accum1(num) + sample
     N1(num) = N1(num) + 1
   end subroutine AntiAlias_put

   real*8 function AntiAlias_get(num)
     integer, intent(in) :: num

     if (firstCall) then
        firstCall = .false.
        accum2 = accum1 !prevent being 50% low on first call
        N2 = N1
     end if
     
     AntiAlias_get = (accum1(num) + accum2(num)) / DBLE(N1(num) + N2(num) )
     accum2(num) = accum1(num)
     accum1(num) = 0d0
     N2(num) = N1(num)
     N1(num) = 0
   end function AntiAlias_get
 end module AntiAliasingFilter

    module LockInDataType
      type :: LockInData
         real*8 :: LockInX, LockInY, LockInX_old, LockInY_old
         real*8 :: LockinXpr, LockInYpr, LockinXpr_old, LockInYpr_old

         real*8 :: LockIn2X,   LockIn2Y,   LockIn2X_old, LockIn2Y_old  !these needed for 4th order filter
         real*8 :: Lockin2Xpr, LockIn2Ypr, Lockin2Xpr_old, LockIn2Ypr_old
         
         real*8 :: R, Th
      end type LockInData
    end module LockInDataType


    ! this whole business is confusing because we have multiple different use cases.  there's approaching as the transient
    ! portion of a scan, there's approaching for fz curves, and there's approaching for dynamic curves.
  module Approaching
    real*8  AprchS_dim, CurveRate_dim, CurveRate
    real*8, private :: AprchS_cur , AprchS_setpoint
    real*8, private :: parabola_a, parabola_b, parabola_t0
    integer, private :: type
    integer, private, parameter :: LINEAR=1, PARABOLA=2
    contains

      
      subroutine SetApproachSpeedSetpoint(AprchS_in)
        real*8, intent(in) :: AprchS_in
        AprchS_setpoint = AprchS_in
        type = LINEAR
      end subroutine SetApproachSpeedSetpoint

      !this is for the rounded part of the rounded triangle
      subroutine SetApproachParabola(a,b,t0)
        real*8, intent(in) :: a,b,t0
        type = PARABOLA
        parabola_a=a
        parabola_b=b
        parabola_t0=t0
      end subroutine SetApproachParabola
      
      subroutine ApproachStop()
        AprchS_cur = 0
        type = LINEAR
      end subroutine ApproachStop

      subroutine ApproachStart()
        AprchS_cur  = AprchS_setpoint
        type = LINEAR
      end subroutine ApproachStart

      subroutine ApproachReverse()
        AprchS_cur = -AprchS_setpoint
        type = LINEAR
      end subroutine ApproachReverse

      real*8 function GetApproachSpeedCurrent(t)
        real*8, intent(in) :: t
        
        if (type == LINEAR) then
           GetApproachSpeedCurrent = AprchS_cur
        elseif (type == PARABOLA) then
           GetApproachSpeedCurrent = 2d0 * parabola_a * (t-parabola_t0) + parabola_b
        else
           call WriteFatalError('unknown approach type')
        end if
        
      end function GetApproachSpeedCurrent

      real*8 function GetApproachSpeedSetpoint()
        GetApproachSpeedSetpoint = AprchS_Setpoint
      end function GetApproachSpeedSetpoint

    end module Approaching
    
    
  module Poincare
    integer, parameter :: maxPoinc=9
    logical :: Want_Strob, Want_Impact
    integer :: X_impact(maxPoinc), Y_impact(maxPoinc),  Npoinc_X(maxPoinc), Npoinc_Y(maxPoinc), numpoinc
  end module Poincare

!this module keeps track of the velocity history which is needed for the convolution integrals
!in the viscoelasticity models.  this is not the cleanest stack implementation, but we want quick
!access to the whole array too, so its sort of a list type as well.
module Stack

  type :: stack_type
     real*8, allocatable :: data(:)
     integer current_loc
     integer buffer_len
  end type stack_type

contains
  subroutine init_stack(s)
    type(stack_type), intent(inout) :: s
    if (.not. allocated( s%data)) then
       s%buffer_len = 1000
!       write(*,*) "attempting to initial allocate"
       allocate(s%data(s%buffer_len))
!       write(*,*) allocated(s%data)
    end if
  end subroutine init_stack
  
    subroutine push(x,s)
      real*8, intent(in) :: x
      real*8, allocatable :: tmp(:)
    
      type(stack_type), intent(inout) :: s

!      write(*,*) s%buffer_len, allocated(s%data)
  
      
      if ((s%current_loc+1) >= s%buffer_len) then       
         if (s%buffer_len > 10000000 ) then
            !10 million points is somewhat arbitrary, but beyond this things will start to get slow.  could
            !potentially increase the limit if really needed
             call WriteFatalError( "History for viscoelasticity buffer exceeded (i.e. total time in contact exceeds the programs current capability for viscoelasticity.)  Please bug the developers to increase the capability")
          else
!            write(*,*) "attemp temp buffer allocate"
            allocate( tmp(s%buffer_len))
            tmp = s%data
!            write(*,*) "attempting to re-allocate", (s%buffer_len*5)
            deallocate( s%data)
            allocate( s%data(s%buffer_len * 5))
            s%data(1:s%buffer_len) = tmp
            s%buffer_len = s%buffer_len * 5
            deallocate(tmp)
         end if
      end if
      s%current_loc = s%current_loc + 1
      s%data( s%current_loc ) = x
    end subroutine push
  
    subroutine pop(s)
      type(stack_type), intent(inout) :: s
      s%current_loc = s%current_loc - 1
      call assert( s%current_loc >= 0, 'popped an empty stack')
    end subroutine pop
    
    real*8 function head(s)
      type(stack_type), intent(in) :: s
      !call assert(s%current_loc > 0, 'head empty stack')
      head = s%data(s%current_loc)
    end function head

    subroutine clear_stack(s)
      type(stack_type), intent(inout) :: s
      s%current_loc = 0
    end subroutine clear_stack

    pure integer function length(s)
      type(stack_type), intent(in) :: s
      length = s%current_loc
    end function length

end module Stack

module Interpolation
  contains

    !like extrapolate, but if outside the given range then clip to a saturated value (for computing exp, for which exp(-100) might as well be zero)
    !low and high refer to the numerical order (not necessarily table order)
    real*8 function interpolate_w_clip( x,y, xi, N, clip_low, clip_high)
      integer, intent(in) :: N
      real*8, intent(in) :: x(N), y(N), xi, clip_low, clip_high
      logical decreasing

      decreasing = ( x(2) < x(1))

      if (decreasing) then
         if ( xi > x(1) ) then
            interpolate_w_clip = clip_high
         elseif ( xi < x(N) ) then
            interpolate_w_clip = clip_low
         else
            interpolate_w_clip = interpolate( x, y, xi, N)
         end if
      else
         if ( xi < x(1) ) then
            interpolate_w_clip = clip_low
         elseif ( xi > x(N) ) then
            interpolate_w_clip = clip_high
         else
            interpolate_w_clip = interpolate( x, y, xi, N)
         end if
      end if


    end function interpolate_w_clip

    !assumption: x is monotonic
    real*8 function interpolate( x, y, xi, N)
      integer, intent(in) :: N
      real*8, intent(in) :: x(N), y(N), xi
      integer :: a, b, guess
      logical decreasing

      decreasing = ( x(2) < x(1))
            
      if (decreasing) then
         call assert( xi <= x(1), 'xi<=x(1) d')
         call assert( xi >= x(N), 'xi>=x(N) d')
      else
         call assert( xi <= x(N), 'xi<=x(N) i')
         call assert( xi >= x(1), 'xi>=x(1) i')
      end if

      !first, locate the two relevant indices by binary search
      a = 0
      b = N
      
      do while ( (b-a) > 1)
         guess = (a+b)/2
         if ((x(guess) > xi) .eqv. decreasing) then
            a = guess
         else
            b = guess
         end if            
      end do

      call assert( a > 0, "interp failed a > 0")
      call assert( b > 0, "interp failed b > 0")
      call assert( a <= N, "interp failed a <= N")
      call assert( b <= N, "interp failed b <= N")

      !then linear interpolation
      interpolate = y(a) + ( y(b) - y(a) ) * ( xi - x(a) ) / ( x(b) - x(a) )
    end function interpolate

    !could make this even faster by pre-computing dx and y(b)-y(a).  but 
    !this is good for now.
    !assumption x monotonically increasing and equally spaced
    real*8 function linear_table_lookup( x, y, xi, N)
      integer, intent(in) :: N
      real*8, intent(in) :: x(N), y(N), xi
      real*8 :: dx
      integer :: a,b
            
      dx = (x(2) - x(1))
      a = 1 + floor( xi /  dx ) 

      b=a+1

!      call assert( xi >= x(a), 'xi >= x(a)')
!      call assert( xi <= x(b), 'xi <= x(b)')
      
      linear_table_lookup = y(a) + ( y(b) - y(a) ) * ( xi - x(a) ) / dx
    end function linear_table_lookup

end module Interpolation

module Derivatives
contains

  !O(h^2) centered diff first order diffs at boundaries
  function centered_diff( x, n, dt)
    integer, intent(in) :: n
    real*8, intent(in) :: x(n), dt
    real*8 :: centered_diff(n)

    centered_diff(1) = (x(2) - x(1)) /dt
    centered_diff(2:n-1) = ( x(3:n) - x(1:n-2)) / (2d0 * dt)
    centered_diff(n) = ( x(n) - x(n-1)) / dt
  end function centered_diff


  !O(h^4) centered diff first order diffs at boundaries
  function centered_diff_4th( x, n, dt)
    integer, intent(in) :: n
    real*8, intent(in) :: x(n), dt
    real*8 :: centered_diff_4th(n)
    integer :: i

    centered_diff_4th(1) = (x(2) - x(1)) /dt
    centered_diff_4th(2) = (x(3) - x(2)) /dt

    do i = 3, n-2
       centered_diff_4th(i) = ( -x(i+2) + 8d0 * x(i+1) - 8d0 * x(i-1) + x(i-2) ) / (12d0 * dt)
    end do
    
    centered_diff_4th(n-1) = ( x(n-1) - x(n-2)) / dt
    centered_diff_4th(n)   = ( x(n) - x(n-1)) / dt
    
  end function centered_diff_4th
  
end module Derivatives

!implements a circular buffer, which we will use as a delay element for self-excitation
module CircularBuffer
  real*8, allocatable :: buffer(:)
  integer :: loc, len
  contains

    subroutine InitCircularBuffer( len_in)
      integer, intent(in) :: len_in    
      len = len_in
      allocate( buffer(0:(len-1)))
      loc = 0
    end subroutine InitCircularBuffer

    real*8 function ReadCircularBuffer()
      ReadCircularBuffer = buffer(loc)
    end function ReadCircularBuffer
    
    subroutine WriteCircularBuffer(x)
      real*8, intent(in) :: x
      buffer(loc) = x
      loc = mod((loc + 1), len)
    end subroutine  WriteCircularBuffer

  end module CircularBuffer

module Integrals
  contains

  !trapezoidal integration with uniform spacing
  real*8 function trapzu( x, n, dt)
    integer , intent(in) :: n
    real*8, intent(in) :: x(n), dt
    
    if (n <= 1) then
       trapzu = 0d0
    else
       trapzu = (dt / 2d0) * ( x(1) + x(n) + 2d0 * sum(x(2:(n-1))))
    end if
  end function trapzu



  !trapezoidal integration
  real*8 function trapz( x, n, dt, t, uniform_spacing)
    logical, intent(in) :: uniform_spacing
    integer , intent(in) :: n
    real*8, intent(in) :: x(n), dt, t(n)
    
    if (n <= 1) then
       trapz = 0d0
    else
       if (uniform_spacing) then
          !optimize for the common case
          trapz = (dt / 2d0) * ( x(1) + x(n) + 2d0 * sum(x(2:(n-1))))
       else
          !but allow for the non-uniform case too
          trapz = sum(  (x(1:n-1) + x(2:n)) * (t(2:n) - t(1:n-1))) / 2d0
       end if
    end if
  end function trapz

  !trapezoidal integration, possibly with uniform spacing but an additional last point that is not at the 
  !same spacing.  thought this would come up a lot with ddaskr and viscoelasticity.  actually not,
  !because I use the root finding, so the first point is never equally spaced either.  so really
  !the usefulness of this function is just that you don't have to append arrays
  real*8 function trapz_p1( x, n, dt, t, uniform_spacing, x2, dt2)
    integer , intent(in) :: n
    real*8, intent(in) :: x(n), dt, x2, dt2, t(n)
    logical, intent(in) :: uniform_spacing
    
    trapz_p1 = trapz(x,n,dt, t, uniform_spacing) +  ( x(n) + x2) * dt2 / 2d0
  end function trapz_p1

  !simpson's rule.  assumption: first point and last points are not equally spaced, but all others are.
!fixme: needs to be tested.  
  real*8 function simpson_p1( x, n, dt, t, x2, dt2)
    integer , intent(in) :: n
    real*8, intent(in) :: x(n), dt, x2, dt2, t(n)

    if (n <= 2) then
       !just use trapezoids.  won't improve anything here.
       simpson_p1 = trapz_p1( x, n, dt, t, .false. , x2, dt2)
    elseif (  (floor(n/2d0)*2) /= n) then
       !odd
       simpson_p1 = ( x(1) + x(2) ) * (t(2) - t(1)) / 2d0 + ( x(n) + x2 ) * dt2 / 2d0 + simpson(x(2:(n-1)), n-2, dt)
    else
       !even
       simpson_p1 = ( x(1) + x(2) ) * (t(2) - t(1)) / 2d0 + ( x(n) + x2 ) * dt2 / 2d0 + ( x(2) + x(3) ) * dt / 2d0 +  simpson(x(3:(n-1)), n-3, dt)
    end if

  end function simpson_p1


  !simpson's rule.  assumption: n is odd, equally spaced grid points.
  real*8 function simpson(x, n, dt)
    integer , intent(in) :: n
    real*8, intent(in) :: x(n), dt

    if (n <= 1) then
       simpson = 0
    elseif (n==3) then
       simpson = (dt / 3d0) * ( x(1) + 4d0 * x(2) + x(3) )
    else
       simpson = (dt / 3d0) * ( x(1) + 4d0 * sum(x(2:(n-1):2)) + 2d0 * sum(x(3:(n-2):2)) + x(n))
    end if
  end function simpson
    
end module Integrals

module LinearAlgebra
 contains
  subroutine solve_gauss_elim(A, b, ok)
    real*8, intent(inout), dimension(:,:) :: a 
    real*8, intent(inout), dimension(:) :: b  !The rhs vector on entry. Overwritten by solution.
    logical, intent(out) :: ok  !True after a successful entry and false otherwise.

    real*8 :: piv_tol !Smallest acceptable pivot.
    integer :: i !Row index
    integer :: j !Column index
    integer :: n !Matrix order
    real*8, dimension(size(b)) :: row
    real*8 :: element !Workspace variable.

    piv_tol = maxval(abs(a))*1.0e-10

    n = size(b)
    ok = ((size(a, 1) == n) .and. (size(a,2) == n))
    if (.not.ok) then
       return
    end if
    do j = 1, n
       ! Update elements in column j.
       do i = 1, j-1
          a(i+1:n, j) = a(i+1:n, j) - a(i,j)*a(i+1:n,i)
       end do
       ! Find pivot and check its size (using maxval just to obtain a scalar).
       i = maxval(maxloc(abs(a(j:n, j)))) + j - 1
       !maxval and maxloc
       if (abs(a(i,j)) < piv_tol) then
          ok=.false.
          return
       end if
       ! If necessary apply row interchange
       if (i/=j) then
          row = a(j,:); a(j,:) = a(i,:); a(i,:) = row
          element = b(j); b(j) = b(i); b(i) = element
       end if
       !Compute elements j+1 : n of j-th column.
       a(j+1:n,j) = a(j+1:n, j) / a(j,j)
    end do
    !Forward substitution
    do i = 1, n-1
       b(i+1:n) = b(i+1:n) - b(i)*a(i+1:n,i)
    end do
    !Back substitution
    do j = n, 1, -1
       b(j) = b(j)/a(j,j)
       b(1:j-1) = b(1:j-1) - b(j)*a(1:j-1,j)
    end do
  end subroutine solve_gauss_elim
end module LinearAlgebra

module Elasticity_Formulas
  contains
        real*8 pure function reducedModulus( E1, nu1, E2, nu2)
          real*8, intent(in) :: E1, nu1, E2, nu2
          reducedModulus  = 1/((1-nu1**2)/E1+ (1-nu2**2)/E2)   
        end function reducedModulus

        real*8 pure function reducedModulus1( E1, nu1)
          real*8, intent(in) :: E1, nu1
          reducedModulus1  = E1 / (1-nu1**2)
        end function reducedModulus1

        real*8 pure function ShearToReducedModulus( G, nu)          
          real*8, intent(in) :: G, nu
          ShearToReducedModulus = 2d0 * G / ( 1 - nu)
        end function ShearToReducedModulus

        real*8 pure function YoungsToShear( E, nu)
          real*8, intent(in) :: E, nu
          YoungsToShear = E / ( 2d0 * (1+nu))
        end function YoungsToShear
                
end module Elasticity_Formulas

module GridGeneration
contains
  !not the same as as matlab' logspace.  this is more intuitive to me
  real*8 function logspace( x1, x2, N) result(ls)
    real*8, intent(in) :: x1, x2
    integer, intent(in) :: N 
    integer i
    dimension  ls(N)

    ls = linspace( log10(x1), log10(x2), n)

    do i = 1,N
       ls(i) = 10**ls(i)
    end do

  end function logspace

  real*8 function linspace(x1, x2, N)
    real*8, intent(in) :: x1, x2
    integer,intent(in) :: N
    integer ::  i
    dimension linspace(N)

    do i = 1,N
       linspace(i) = x1 + (x2 - x1) * dble(i-1)/dble(N-1)
    end do
    
  end function linspace

  function quadspace(x0, xf, N, bias)
    real*8, intent(in) :: x0, xf, bias
    integer, intent(in) :: N
    real*8, dimension(N) :: quadspace

    quadspace = linspace(x0, xf, N) * (linspace( sqrt(1d0/bias), 1d0 ,N) ** 2)
  end function quadspace

end module GridGeneration

! I never finished this.  
!   module Photothermal

!   contains

!     !photothermal calculations reference Ramos JAP 06
!     pure real*8 function Phototh_DiffLen(K, beta, omega)
!       real*8, intent(in) :: K, beta, omega
!       Phototh_DiffLen = sqrt(( 2 * K) / (beta + sqrt(beta**2+omega**2)))
!     end function Phototh_DiffLen

!     pure real*8 function Phototh_WaveLen(K, beta, omega)
!       real*8, intent(in) :: K, beta, omega
!       Phototh_DiffLen = 2 * pi * sqrt(( 2 * K) / (-beta + sqrt(beta**2+omega**2)))
!     end function Phototh_DiffLen

!     !calc delta temp at a given location on the cantilever
!     pure complex*8 function Phototh_Del(K, beta, omega,x, x0)
!       real*8, intent(in) :: K, beta, omega,x, x0
!       real*8, del, lam
!       complex*8 

!       del = Phototh_DiffLen(K, beta, omega)
!       lam = Phototh_WaveLen(K, beta, omega)
!       !this is approximation.  see Ramos for full equation if we want more detail
!       real*8, parameter :: r = 1, theta = 0
!       if (x < x0) then
!           Phototh_DelT = Phototh_DelT0 * ( exp( CMPLX(1 / del, 2*pi/lam) * (x - x0)) + r * exp( CMPLX(1/del, 2*pi/lam) * (x + x0 - 2 * L) + CMPLX(0, theta)))
!        else
!           Phototh_DelT = Phototh_DelT0 * ( -exp( CMPLX(1 / del, 2*pi/lam) * (x - x0)) + r * exp( CMPLX(1/del, 2*pi/lam) * (x + x0 - 2 * L) + CMPLX(0, theta)))
!        end if
!      end function Phototh_DelT

!    end module Photothermal


module ellipticIntegrals
  !from numerical recipes in fortran 77
  private rf
  contains

    FUNCTION rf(x,y,z)
      real*8 rf,x,y,z,ERRTOL,TINY,BIG,THIRD,C1,C2,C3,C4
      PARAMETER (ERRTOL=.0025,TINY=1.5e-38,BIG=3.E37,THIRD=1./3.,C1=1./24.,C2=.1,C3=3./44.,C4=1./14.)
      !Computes Carlson’s elliptic integral of the first kind, RF (x, y, z). x, y, and z must be
      !nonnegative, and at most one can be zero. TINY must be at least 5 times the machine
      !underflow limit, BIG at most one fifth the machine overflow limit.
      REAL*8 alamb,ave,delx,dely,delz,e2,e3,sqrtx,sqrty,sqrtz,xt,yt,zt
      if(min(x,y,z).lt.0..or.min(x+y,x+z,y+z).lt.TINY.or.max(x,y,z).gt.BIG) call WriteFatalError( "invalid arguments in rf")
      xt=x
      yt=y
      zt=z
1     continue
      sqrtx=sqrt(xt)
      sqrty=sqrt(yt)
      sqrtz=sqrt(zt)
      alamb=sqrtx*(sqrty+sqrtz)+sqrty*sqrtz
      xt=.25*(xt+alamb)
      yt=.25*(yt+alamb)
      zt=.25*(zt+alamb)
      ave=THIRD*(xt+yt+zt)
      delx=(ave-xt)/ave
      dely=(ave-yt)/ave
      delz=(ave-zt)/ave
      if(max(abs(delx),abs(dely),abs(delz)).gt.ERRTOL)goto 1
      e2=delx*dely-delz**2
      e3=delx*dely*delz
      rf=(1.+(C1*e2-C2-C3*e3)*e2+C4*e3)/sqrt(ave)
      return
    END FUNCTION rf

!     FUNCTION rd(x,y,z)
!       REAL*8 rd,x,y,z,ERRTOL,TINY,BIG,C1,C2,C3,C4,C5,C6
!       PARAMETER (ERRTOL=.0015,TINY=1.e-25,BIG=4.5E21,C1=3./14.,C2=1./6.,C3=9./22.,C4=3./26.,C5=.25*C3,C6=1.5*C4)
!       !Computes Carlson’s elliptic integral of the second kind, RD(x, y, z). x and y must be
!       !nonnegative, and at most one can be zero. z must be positive. TINY must be at least twice
!       !the negative 2/3 power of the machine overflow limit. BIG must be at most 0.1 ×ERRTOL
!       !times the negative 2/3 power of the machine underflow limit.
!       REAL*8 alamb,ave,delx,dely,delz,ea,eb,ec,ed,ee,fac,sqrtx,sqrty, sqrtz,sum,xt,yt,zt
!       if(min(x,y).lt.0..or.min(x+y,z).lt.TINY.or.max(x,y,z).gt.BIG) call WriteFatalError("invalid arguments in rd")
!       xt=x
!       yt=y
!       zt=z
!       sum=0.
!       fac=1.
! 1     continue
!       sqrtx=sqrt(xt)
!       sqrty=sqrt(yt)
!       sqrtz=sqrt(zt)
!       alamb=sqrtx*(sqrty+sqrtz)+sqrty*sqrtz
!       sum=sum+fac/(sqrtz*(zt+alamb))
!       fac=.25*fac
!       xt=.25*(xt+alamb)
!       yt=.25*(yt+alamb)
!       zt=.25*(zt+alamb)
!       ave=.2*(xt+yt+3.*zt)
!       delx=(ave-xt)/ave
!       dely=(ave-yt)/ave
!       delz=(ave-zt)/ave
!       if(max(abs(delx),abs(dely),abs(delz)).gt.ERRTOL)goto 1
!       ea=delx*dely
!       eb=delz*delz
!       ec=ea-eb
!       ed=ea-6.*eb
!       ee=ed+ec+ec
!       rd=3.*sum+fac*(1.+ed*(-C1+C5*ed-C6*delz*ee)+delz*(C2*ee+delz*(-C3*ec+delz*C4*ea)))/(ave*sqrt(ave))
!       return
!     END FUNCTION rd

    !complete elliptic integral of the first kind (had accidentally used the integral of 2nd kind previously)
    real*8 function elliptic_k(k)
      real*8, intent(in):: k
      elliptic_k = rf(0d0, 1d0 - k**2, 1d0)
    end function elliptic_k
  end module ellipticIntegrals

module ZiegerNichols
  contains
    pure subroutine ZN_PI( k_ultimate, tau_ultimate, kp, ki, kd)
      real*8, intent(in) :: k_ultimate, tau_ultimate
      real*8, intent(out) :: kp, ki, kd
      Kp = k_ultimate / 2.2
      Ki = 1.2 * Kp / tau_ultimate 
      Kd = 0d0
    end subroutine ZN_PI

    pure subroutine ZN_PID( k_ultimate, tau_ultimate, kp, ki, kd)
      real*8, intent(in) :: k_ultimate, tau_ultimate
      real*8, intent(out) :: kp, ki, kd
      Kp = k_ultimate * 0.6
      Ki = 2d0 * Kp / tau_ultimate 
      Kd = kp * tau_ultimate / 8d0
    end subroutine ZN_PID
  end module ZiegerNichols

  module bisection
  contains
    !generic bisection search
    !tolerance is on the function value.  could make it more general purpose and put in Z as well
    subroutine bisection_search(f, Z1_in, error1_in, Z2_in, error2_in, f_tolerance, Z_out)
      real*8, intent(in) :: Z1_in, error1_in, Z2_in, error2_in, f_tolerance
      real*8, intent(out) :: Z_out
      real*8 :: f  !function that returns a real*8,
      
      real*8 :: error_new
      real*8 :: Z1, error1, Z2, error2
      integer :: count
      
      Z1 = Z1_in
      Z2 = Z2_in
      error1 = error1_in
      error2 = error2_in
      
      error_new = 2 * f_tolerance !also do at least one iteration

      count = 1
      
      do while ( abs(error_new) > f_tolerance) 
         Z_out = (Z1+Z2)/2
         error_new = f( Z_out)

         !write(*,*) 'bisection Z ', Z_out, 'error ', error_new
         
         if ( error_new * error1 < 0) then
            Z2 = Z_out
            error2 = error_new
         else
            Z1 = Z_out
            error1 = error_new
         end if
         count = count + 1
         if (count > 50) then

            call WriteFatalError('Error: Bisection failed to converge. This may indicate a bistability or non-monotonic / discontinuous amplitude vs Z. The averaging tool may not be suitable for this simulation. Try the basic AM approach curves')
!            write(*,*) "Error: Bisection failed to converge. continuing anyway"
            return
         end if
      end do
    end subroutine bisection_search

    
  end module bisection

  module ODE_solvers
    
  contains

  function fwdEuler(f, Y, t, dt)
        real*8, intent(in) :: t, dt
        real*8, intent(in) :: Y(:)

        interface
           function f(x,t)
             real*8, intent(in) :: x(:), t
            real*8, dimension(size(x)) :: f
          end function f
       end interface
        
        real*8, dimension(size(Y)) :: fwdEuler
     
        fwdEuler = Y + f(Y, t)*dt
        
      end function fwdEuler
      
    
    ! !generic runge-kutta 4th order solver. returns next Y
    ! f is a function that takes Y,t as input, returns Y'
      function rk4(f, Y, t, dt)
        real*8, intent(in) :: t, dt
        real*8, intent(in) :: Y(:)

        interface
           function f(x,t)
             real*8, intent(in) :: x(:), t
            real*8, dimension(size(x)) :: f
          end function f
       end interface
        
        real*8, dimension(size(Y)) :: rk4
        real*8, dimension(size(Y)) :: YP1, YP2, YP3, YP4
      
        YP1 = f(Y, t)
        YP2 = f(Y + 0.5d0 * YP1 * dt, t + 0.5d0 * dt)
        YP3 = f(Y + 0.5d0 * YP2 * dt, t + 0.5d0 * dt)
        YP4 = f(Y +         YP3 * dt, t +         dt)
        
        rk4 = Y + (YP1 + 2d0* YP2 + 2d0* YP3 + YP4)*dt/6d0
        
      end function rk4


! 7/7/2012.  I added these subroutines when I added the Attard model.  My thought at the time was that DDASKR would
! be too slow for that model. The reason being that DDASKR is an implicit method, there it needs to evaluate the jacobian
! (which we currently do by finite difference) and then does non-linear iterations.  This results in many function
! evaluations per time step.  I therefore thought that an explicit method like runge-kutta or adams-bashforth would
! be better because there are fewer function evaluations per time step (don't need the jacobian, etc).  However,
! to maintain stability with the explicit solvers, I needed to use very small time steps.  I had never actually tried
! solving the Attard model with DDASKR.  Turns out, it's really not all that slow.  Yes Attard model is slower than Hertz
! model, but for N_attard=50, it's not unreasonable, and I can use 10 times larger step with DDASKR and still maintain
! stability.  So, these routines are unnecessary, but I am leaving them in just in case you want to see them.

! !adams bashforth 4th order solver. This is same order as the runge-kutta solver but with four times less function evaluations
! !because the function evaluation is so expensive, this means the entire simulation is about four times faster.
! !it would be fairly easy to extent this to a fifth order method
! integer function ab4(t,Y,YPRIME, dt)
!   use data1, only : NEQ
!   real*8, intent(in) :: dt
!   real*8, intent(inout) :: Y(NEQ), YPRIME(NEQ), t
!   real*8 :: CJ, RPAR, IPAR, IRES
!   real*8, allocatable, save :: DELTA1(:), DELTA2(:), DELTA3(:), DELTA4(:)
!   integer, save :: number_calls = 0

!   number_calls = number_calls+1

!   if (number_calls == 1) then
!      !initialize. first call. use rk4 for first steps until we've built up four of them
!      allocate( DELTA1(NEQ))
!      allocate( DELTA2(NEQ))
!      allocate( DELTA3(NEQ))
!      allocate( DELTA4(NEQ))
     
!      call rk4(t,Y,YPRIME, dt)
!      DELTA4 =  YPRIME
!   elseif (number_calls == 2) then
!      call rk4(t,Y,YPRIME, dt)
!      DELTA3 =  YPRIME
!   elseif (number_calls == 3) then
!      call rk4(t,Y,YPRIME, dt)
!      DELTA2 =  YPRIME
!   else
!      YPRIME = 0
!      !adapt existing residual routine from ddaskr.  if we set YPRIME=0, then DELTA gives us the new Y'
!      call RES1(t, Y, YPRIME,CJ,DELTA1,IRES,RPAR,IPAR) 
     
!      Y = Y + dt * ( 55d0/24d0 * DELTA1 - 59d0/24d0 * DELTA2 + 37d0/24d0 * DELTA3 - 3d0/8d0 * DELTA4 )
     
!      DELTA4 = DELTA3 ! this is a lot of memory copying.  if that gets to be a bottleneck we can be smarter about it.
!      DELTA3 = DELTA2
!      DELTA2 = DELTA1
!   end if

!   t = t + dt
!   ab4 = IRES
! end function ab4

      
    end module ODE_solvers

    module Fractions
    contains

      integer function gcd(a_in, b_in)
        integer, intent(in) :: a_in,b_in
        integer :: a,b, old_a
        !Calculate the Greatest Common Divisor of a and b.
        !
        !Unless b==0, the result will have the same sign as b (so that when
        !b is divided by it, the result comes out positive).        
        a = a_in
        b = b_in
        do while( b > 0)
           old_a = a

           a = b
           b = mod(old_a,b)
        end do
        gcd = a
      end function gcd
      
      !this routine was based on the python function limit_denominator
      !https://svn.python.org/projects/python/trunk/Lib/fractions.py
      !the python version takes a fraction as in input, we take a real number,
      !and then round that to a fraction out of 1 million, and then remove gcd
      !
      
      subroutine fraction_limit_denominator( a_in, max_denom, p, q)
        real*8, intent(in) :: a_in
        integer, intent(in) :: max_denom

        integer, intent(out) :: p,q
        real*8 :: bound1, bound2
        integer :: p0,q0, p1, q1, n, d, q2, a, k, old_d, old_n, old_p0, g

        d = 1000000
        n = int( a_in * d)
        g = gcd(n, d);
        n = n / g
        d = d / g

        if (d <= max_denom) then
           !already satisfied!
           p = n
           q = d
           return
        end if
        
!        write(*,*) 'n d ', n,d
        call assert( max_denom > 1, ' max_denom > 1' )
        call assert( n > 0, 'n>0')
        
        p0=0
        q0=1
        p1=1
        q1=0

        do while(.true.)
!           write(*,*) 'n d ', n, d
           call assert( d > 0, 'd>0')
           
           a = n / d
           q2 = q0 + a * q1
           if (q2 > max_denom) then
              exit
           end if

           old_p0 = p0

           p0 = p1
           q0 = q1
           p1 = old_p0 + a * p1
           q1 = q2
           
           old_d = d
           old_n = n
           n = old_d
           d = old_n - a * old_d
        end do

        call assert( q1 > 0, 'q1>0')
        k = (max_denom - q0) / q1
        
        bound1 = (dble(p0+k*p1)) / ( dble(q0+k*q1))
        bound2 = dble(p1) / dble(q1)

!        write(*,*) 'bound1 ', bound1, ' bound2 ', bound2
        
        if (abs(bound2 - a_in) < abs(bound1 - a_in)) then
           p = p1
           q = q1
        else
           p=p0+k*p1
           q=q0+k*q1
        end if
        
      end subroutine fraction_limit_denominator
    end module Fractions

    !wanted to put this in the triggeredFzMode module, but had to pull it out for circular dependency reasons
    module checkFzCurves
      contains
        logical function isFzCurves(operating_mode, fexcite)
          use params
          integer, intent(in) :: operating_mode, fexcite
          isFzCurves = (((operating_mode == APPROACH) .OR. (operating_mode == APPR_RET) .or. ( operating_mode == APPROACH_TRIGGERED) .OR. (operating_mode == APPROACH_SINE ) ) .AND. (fexcite == NO_EXC))
        end function isFzCurves
      end module checkFzCurves
