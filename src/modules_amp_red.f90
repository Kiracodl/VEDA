
module amp_red

  real*8, private :: Asp1_target, Asp2_target,  Ainitial_nd(2), Z_target, Z,  Asp2_1, Asp2_2
  integer, private :: ode_solver, max_f1_cycles, max_f2_cycles

  integer :: cumulative_iteration
  real*8 :: Asp_tolerance, delta_Asp, Asp2_err
  real*8 :: virial1, dissipation1, virial2, dissipation2, Amp1, Amp2 !this will store the most recent calculation from each iteration

  !note gap has the same meaning as in rest of VEDA, gap is w.r.t. the undeformed surface
  !this set of variables stores the entire time history for output.  force may be zero for many or most of the points
  real*8, allocatable :: time(:), gap(:), dpr(:), force(:), attard_surf_coord(:), q1(:), q2(:), q1d(:), q2d(:)
  
  !this is for Attard.  putting here and storing all time points so we can use for movies
  real*8, allocatable :: Y(:,:) 
  logical, allocatable :: attard_enabled_history(:)
  
contains

  subroutine set_Asp1_target(Asp1_in)
    real*8, intent(in) :: Asp1_in
    Asp1_target = Asp1_in
  end subroutine set_Asp1_target

  !used for bimodal only
  subroutine set_Z_target(Z_in)    
    real*8, intent(in) :: Z_in
    !write(*,*) "setting Z target: ", Z_in
    Z_target = Z_in
  end subroutine set_Z_target

  
  subroutine set_Asp2_target(Asp2_in)
    real*8, intent(in) ::  Asp2_in    
    Asp2_target = Asp2_in
  end subroutine set_Asp2_target

  !this is actual deflection, not observed
  subroutine set_ainitial( Ainitial_nd_in)
    real*8, intent(in) :: Ainitial_nd_in(2)
    Ainitial_nd = Ainitial_nd_in
  end subroutine set_ainitial

  subroutine init_amp_red(ode_solver_in, max_f1_cycles_in, max_f2_cycles_in)    
    use TimeAndCycle, only: numincycle
    use data1, only:  NEQ
    integer, intent(in) :: ode_solver_in, max_f1_cycles_in, max_f2_cycles_in

    max_f1_cycles = max_f1_cycles_in
    max_f2_cycles = max_f2_cycles_in
    
    if (.not.(allocated(time))) allocate(time(numincycle*max_f1_cycles))
    if (.not.(allocated(gap))) allocate(gap(numincycle*max_f1_cycles))
    if (.not.(allocated(dpr))) allocate(dpr(numincycle*max_f1_cycles))

    if (.not.(allocated(q1))) allocate(q1(numincycle*max_f1_cycles))
    if (.not.(allocated(q2))) allocate(q2(numincycle*max_f1_cycles))
    if (.not.(allocated(q1d))) allocate(q1d(numincycle*max_f1_cycles))
    if (.not.(allocated(q2d))) allocate(q2d(numincycle*max_f1_cycles))
    
    if (.not.(allocated(force))) allocate(force(numincycle*max_f1_cycles))
    if (.not.(allocated(attard_surf_coord))) allocate(attard_surf_coord(numincycle*max_f1_cycles))

    if (.not.(allocated(Y))) allocate(  Y(numincycle*max_f1_cycles, NEQ) )

    if (.not.(allocated(attard_enabled_history))) allocate(  attard_enabled_history(numincycle*max_f1_cycles))
        
    !initial guess for first iteration. after that will reuse the final value from previous iteration
    Asp2_1 = 1.0
    Asp2_2 = 0.95
    
    
    ode_solver = ode_solver_in
  end subroutine init_amp_red

  !for bimodal outer loop.  input Z, returns Asp1 error
    real*8 function amplitude_reduction_1_error_bimodal(Z_new)
      use params
      use NonDimensionalization
      use rappture_io, only: OutputIterationStatsPlots_A2
      real*8, intent(in) :: Z_new
      real*8 :: error1, error2, Asp2 
      error1 = nan()
      error2 = nan()
      
      
      call set_Z_target(Z_new)
      call sequential_search_brent(amplitude_reduction_2_error , cumulative_iteration, Asp2_1, Asp2_2, error1, error2, Asp_tolerance, delta_Asp/2d0 , 10d0 ,  0d0 ,  OutputIterationStatsPlots_A2,  Z_new , Asp2, Asp2_err )
      call set_Asp2_target( Asp2)
      amplitude_reduction_1_error_bimodal =  amplitude_reduction_1_error( Z_new)
    end function amplitude_reduction_1_error_bimodal

  
  !fixme: would this be better conditioned if it was an Amp error not an Asp error? that's how Bahram did it
    
  !for bimodal innter loop.  input Asp2, returns Asp2 error
  !sign of error is such that Z iteration and A2 iteration move same way w.r.t. seq search
real*8 function amplitude_reduction_2_error( Asp2_in)
  use Nondimensionalization
  use params
  real*8, intent(in) :: Asp2_in
  real*8 :: new_Asp1,  new_Asp2
  
  call amplitude_reduction_iteration( Z_target, Asp1_target, Asp2_in, new_Asp1, new_Asp2 )
!  write(*,*) Asp2_in, new_Asp2
  amplitude_reduction_2_error = Asp2_in- new_Asp2
  !write(*,*)   "amplitude_reduction_2_error, Asp2_in: ",Asp2_in, "Asp1_target: ", Asp1_target, " Z_target: ", Z_target, " Asp2_out: ", new_Asp2, " error: ", amplitude_reduction_2_error
  !write(*,*)  "virial2 ", DimenEnergy(virial2* electrons_per_coulomb) , " dissipation2: ", DimenEnergy( dissipation2* electrons_per_coulomb) , " eV/mode2drivecycle"
end function amplitude_reduction_2_error

  
!for tapping mode, or called from within bimodal outer loop.  input Z, returns Asp1 error
real*8 function amplitude_reduction_1_error( Z_in)
  real*8, intent(in) :: Z_in
  real*8 :: new_Asp1,  new_Asp2
  
  call amplitude_reduction_iteration( Z_in, Asp1_target, Asp2_target, new_Asp1, new_Asp2 )
  Asp2_err                   =  new_Asp2 - Asp2_target
!write(*,*)   " Asp2_err, new_Asp2, Asp2_target:",   Asp2_err ,  new_Asp2 , Asp2_target
  amplitude_reduction_1_error = new_Asp1 - Asp1_target
end function amplitude_reduction_1_error


!inputs: Z & setpoint ratio 1 & 2 (for tapping mode leave 2 at 0)
!returns new setpoint ratio 1 & 2
subroutine amplitude_reduction_iteration( Z_in, Asp1_in, Asp2_in, Asp1_out, Asp2_out)
  use Integrals
  use params
  use contactModels !this is where we get cur_props!
  use data0, only: Keq
  use data1, only: Quality, numModes, NEQ
  use TimeAndCycle, only: omegad, numincycle, delt
  use ODE_solvers
  
  real*8, intent(in)  :: Z_in, Asp1_in, Asp2_in
  real*8, intent(out) :: Asp1_out, Asp2_out
  integer :: I_out !, I_comp
  logical :: fail
!  real*8 :: Ycheck1(NEQ),  Ycheck2(NEQ), Force_check, err_est  ! was for debugging mostly. removed.  
  real*8 :: Keq1_nd, Keq2_nd
  integer :: n1, x1, x2

  Amp1 = Asp1_in * Ainitial_nd(1)
  Amp2 = Asp2_in * Ainitial_nd(2)
  Z = Z_in ! so ODE yprime function can be only a function of time... oh how I wish FORTRAN had closures
  
  Force = 0
  attard_surf_coord=0
  Y = 0
  call InitializeHystereticModels(Z_in)
  
    
  !    based on Z, amp: generate time history of tip-sample gap
  do I_out = 1, numincycle*max_f1_cycles
     !generate exactly one cycle for tapping mode, multiple for bimodal
     time( I_out)  = (real( I_out-1)/numincycle) / (omegad(0d0, 1)/(2d0*pi))
     
     gap( I_out) = computeTipSampleGap(     time( I_out), q1(I_out) , q2(I_out) )
     dpr( I_out) = computeTipSampleGapDeriv(time( I_out), q1d(I_out), q2d(I_out) )
  end do

  
  do  I_out = 1,numincycle*max_f1_cycles

     attard_enabled_history(I_out) = attard_enabled
     
     !    generate time history of time-sample force (potentially output as time history)
     !  below is copied closely from force_viewer
          
     call updateHystereticModels_insideRES1(gap(I_out), dpr(I_out)) 

     !y here is just for doing attard.  it's the state of the individual discretization points. 
     
     Force(I_out) = Fts(gap(I_out), dpr(I_out), fail, time(I_out), Y(I_out, :)  )
     
     if (fail) exit     

     attard_surf_coord(I_out) = u0 !u0 is computed as a side effect of Fts.  so must be here, before the force check!

     !fixme reinstate force check
!      if  ((cur_props%fts_ve_model == VEM_ATTARD_BASED) .and. (ode_solver == SOLV_FWDEULER)) then
!            !getting a lot of false positives on the first few points after we turn the interaction on, which have  very low force levels
!            Force_check = Fts(gap(I_out), dpr(I_out), fail, time(I_out), Ycheck2)

!            err_est =  abs(100d0 *  ( Force_check - Force(I_out) ) / Force_check)

!            if ((Force(I_out) .ne. 0d0) .and. (Force_check .ne. 0d0) .and. (err_est > 100)) then
!               if (Force_check > ( Ainitial_nd(1) * Keq(1) / Quality(1) / 1000d0)) then
!                  call WriteFatalError('Convergence check failed.  Most likely cause is too large of time step.  Increase number of time steps. You may also try increasing number of spatial or fourier points, or switching to Runge-Kutta 4th order solver.')
! !              write(*,*) 'Convergence check failed.  Most likely cause is too large of time step.  Increase number of time steps. You may also try increasing number of spatial or fourier points, or switching to Runge-Kutta 4th order solver.'
! !              write(*,*) I_out, DimenLength(gap(I_out)), DimenVelocity(dpr(I_out)), Force(I_out), Force_check
!               end if
!            end if
!       end if
    

     
     !integrate the Attard ODEs here... for purposes of this file, it would be better if attard knows how to integrate itself
     !but that breaks the main VEDA code, and I'd still like to keep that open as a possibility.
     !the n1=2*numModes+1 is a bad hack to maintain compability with old version of VEDA when ODEs were integrated inside DDASKR
     !
     !note: this does need to go AFTER the force computation, not before.
     !otherwise the values of u and d are not computed at a consistent time step (and the LJ potential is very unforgiving
     !about being off by even a few angstroms)
     if ((cur_props%fts_ve_model == VEM_ATTARD_BASED) .and. (attard_enabled)) then
        n1 = 2 * numModes+1
        if (ode_solver == SOLV_FWDEULER) then
!           !this is mostly for debugging. compare 2 half steps to one full step. we might take it out for production
!           Ycheck1(n1:NEQ) = fwdEuler(ODE_yprime, Y(n1:NEQ),       time(I_out), delt/2d0)
!           Ycheck2(n1:NEQ) = fwdEuler(ODE_yprime, Ycheck1(n1:NEQ), time(I_out)+delt/2d0, delt/2d0)
           
           Y(I_OUT+1, n1:NEQ) = fwdEuler(ODE_yprime, Y(I_OUT, n1:NEQ), time(I_out), delt)           
        elseif (ode_solver == SOLV_RK4) then
           !RK4 permits approx a factor of 50 larger time step with same accuracy, so overall about 10x faster                      
           Y(I_OUT+1,n1:NEQ) = rk4(ODE_yprime, Y(I_OUT, n1:NEQ), time(I_out), delt)        
        end if
     else
        n1=1
     end if
     
     call updateHystereticModels_outsideRES1( gap(I_out), dpr(I_out),time(I_out), int( I_out,8 ),  Y(I_OUT+1,n1:NEQ) )
     call Update_viscoelasticity_outsideRes1_Hertz(  time(I_out), gap(I_out) - cur_props%aDMT, dpr(I_out), cur_props, Rtip_nd) 
     call Update_viscoelasticity_outsideRes1_linear( time(I_out), gap(I_out) - cur_props%aDMT, cur_props)

  end do    

  !    integrate force history to get virial and dissipation
  ! this adds up the whole vector, even though we know that large part of it is zeros
  ! could make it faster by keeping track of the non-zero code but that adds complexity and
  ! this is actually fairly fast compared to everything else in the loop
  dissipation1 =             -trapzu( Force * q1d , I_out-1 , delt) / max_f1_cycles
  if (Amp2>0) dissipation2 = -trapzu( Force * q2d , I_out-1 , delt)  / max_f2_cycles
  
!  if (dissipation1 < 0) then
!      call WriteFatalError("computed dissipation is negative. This likely indicates a numerical issue.  Try increasing the number of deflection points per cycle (accuracy vs speed tradeoff on simulation parameters tab). For Attard's model, try increasing number of spatial points, number of fourier basis functions, or extent of computational domain)")
!  end if

  
  !adding the negative sign here such that the results match the sign convention used
  !in the original VEDA implementation. virial gets squared in the amp red equation, so
  !it gives same answer either way.  have to also make corresponding change to phase
  virial1               = - (omegad(0d0, 1) / (2*pi)) * trapzu( Force * q1 , I_out-1 , delt) / max_f1_cycles
  if (Amp2 >0)  virial2 = - (omegad(0d0, 2) / (2*pi)) * trapzu( Force * q2 , I_out-1 , delt) / max_f2_cycles

  
  !based on virial and dissipation, get amplitude reduction  
  !this is setpoint ratio. Rajabifar, macromolecules, eq 4,
  !Keq is dimensionsal, whereas Amp is non-dimensional. force above is non-dimensional, therefore so is virial and dissipation
  Keq1_nd = NonDimenStiffness( Keq(1))
  Asp1_out = 1 /Quality(1)/sqrt(((1/Quality(1))+(dissipation1/pi/ Keq1_nd/(Amp1**2) ))**2+ (2*virial1/Keq1_nd/(Amp1**2))**2)

  if (Amp2>0) then
     Keq2_nd = NonDimenStiffness( Keq(2))
     Asp2_out = 1 /Quality(2)/sqrt(((1/Quality(2))+(dissipation2/pi/ Keq2_nd/(Amp2**2) ))**2+ (2*virial2/Keq2_nd/(Amp2**2))**2)
  else
     Asp2_out = 0d0
  end if
end subroutine amplitude_reduction_iteration



function ODE_yprime(Y, t)
  use Viscoelasticity_attard, only: compute_attard_derivative
  use contactModels, only: cur_props
  real*8, intent(in) :: t, Y(:)
  real*8, dimension(size(Y)) :: ODE_yprime
  real*8 :: gap, dpr

  gap = computeTipSampleGap(t)
  dpr = computeTipSampleGapDeriv(t)

  ODE_yprime = compute_attard_derivative(gap, dpr,t, Y, cur_props)
  
end function ODE_yprime

real*8 function computeTipSampleGap(t, q1_out, q2_out)
  use TimeAndCycle, only: omegad
  real*8, intent(in) :: t
  real*8, intent(out), optional :: q1_out, q2_out
  real*8 :: q1, q2
  q1 = Amp1 * cos( t *  omegad(0d0, 1) )
  if (Amp2>0) then
     q2 = Amp2 * cos( t *  omegad(0d0, 2) )
  else
     q2 = 0
  end if
  computeTipSampleGap   = Z + q1   +  q2
  if (present(q1_out)) then
     q1_out = q1
     q2_out = q2
  end if
end function computeTipSampleGap

real*8 function computeTipSampleGapDeriv(t, q1d_out, q2d_out)
  use TimeAndCycle, only: omegad
  real*8, intent(in) :: t
  real*8, intent(out), optional :: q1d_out, q2d_out
  real*8 :: q1d, q2d
  q1d = - omegad(0d0, 1) * Amp1 * sin( t *  omegad(0d0, 1) )
  if (Amp2>0) then
     q2d = - omegad(0d0, 2) * Amp2 * sin( t *  omegad(0d0, 2) )
  else
     q2d = 0
  end if
  computeTipSampleGapDeriv   = q1d+q2d
  if (present(q1d_out)) then
     q1d_out = q1d
     q2d_out = q2d
  end if
end function computeTipSampleGapDeriv


end module amp_red
