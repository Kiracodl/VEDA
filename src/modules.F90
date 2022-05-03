!-----------------------------------------------------------------------
!             Copyright 2007 - 2020 by
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
!modules file

!this has bit us so many times.  Just put all of it in ONE PLACE!
!
!	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!	Normalization
!	Spatial: L(ND) = L(m)/Ainitial(m)
!       Temporal: t(ND) = t(s) *  omega_scale(rad/s)      
!       Frequency: ND   = Freq (rad/s) / omega_scale (rad/s)
!	Forces: F (ND) = see below
!	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
module Nondimensionalization
  real*8, private :: omega_scale, A_scale

  contains
    !constructors 

    subroutine set_omega_scale(o)
      real*8, intent(in) :: o
      omega_scale = o
    end subroutine set_omega_scale

    subroutine set_a_scale(a)
      real*8, intent(in) :: a
      A_scale = a
    end subroutine set_a_scale

    !non-dimensionalizers

    pure real*8 function NonDimenTime( x )
      real*8, intent(in) :: x
      NonDimenTime = x * omega_scale
    end function NonDimenTime

    pure real*8 function NonDimenFreq( x )
      real*8, intent(in) :: x
      NonDimenFreq = x / omega_scale
    end function NonDimenFreq

    function NonDimenFreqArray( x, n )
      integer, intent(in) :: n
      real*8, intent(in) :: x(n)
      real*8 :: NonDimenFreqArray(n)
      integer :: i
      do i = 1,n
         NonDimenFreqArray(i) = NonDimenFreq(x(i))
      end do
    end function NonDimenFreqArray
      
    pure real*8 function NonDimenEnergy(x)
      real*8, intent(in) :: x
      NonDimenEnergy = NonDimenLength( NonDimenForce( x) )
    end function NonDimenEnergy

    pure real*8 function NonDimenInvLength( x )
      real*8, intent(in) :: x
      NonDimenInvLength = x * A_scale
    end function NonDimenInvLength

    pure real*8 function NonDimenLength( x )
      real*8, intent(in) :: x
      NonDimenLength = x / A_scale
    end function NonDimenLength

    pure real*8 function NonDimenLengthSq( x )
      real*8, intent(in) :: x
      NonDimenLengthSq = NonDimenLength( NonDimenLength( x) )
    end function NonDimenLengthSq

    pure real*8 function NonDimenInvLengthSq( x )
      real*8, intent(in) :: x
      NonDimenInvLengthSq = NonDimenInvLength( NonDimenInvLength( x) )
    end function NonDimenInvLengthSq

    pure real*8 function NonDimenVelocity(x)
      real*8, intent(in) :: x
      NonDimenVelocity =  NonDimenLength(x) / omega_scale
    end function NonDimenVelocity

    pure real*8 function NonDimenPressure(x)
      real*8, intent(in) :: x      

      NonDimenPressure = NonDimenForce( NonDimenInvLengthSq(x))
    end function NonDimenPressure


    !this used to both a non-dimensionalization and a scaling combined, which was bad.
    !non-dimension:  N = kg * m / s^2
    !now omegai_dim^2 / keq(i) = 1/mass_i,  the ith modal mass.  this is the scaling part of it, and it is dependent on the eigenmode
    !then we have units of m/s^2, which are non-dimensionalized by 1/ omega_scale^2 / A_scale, which is mode independent
    !the net result is omegai_nd^2 / a_scale / keq.
    !in the new version we split scaling and non-dimensionalization up into two functions
    pure real*8 function NonDimenForce(x)
      real*8, intent(in) :: x
      
      NonDimenForce = x / omega_scale**2 / A_scale 
    end function NonDimenForce

    !this scales out the ith modal mass.  this is used only from within in RES1. probably should just get rid of this
    !function and put invModalMass into RES1
    pure real*8 function ScaledForce(x,i)
       use data0, only: invModalMass
       real*8, intent(in) :: x
       integer, intent(in) :: i

       ScaledForce = x * invModalMass(i)
     end function ScaledForce

     pure real*8 function NonDimenStiffness(x)
       real*8, intent(in) :: x

       NonDimenStiffness = NonDimenInvLength(NonDimenForce(x))
     end function NonDimenStiffness

     !Dimensionalizers

     pure real*8 function DimenInvLength( x )
       real*8, intent(in) :: x
       DimenInvLength = x / A_scale
     end function DimenInvLength

    pure real*8 function DimenInvLengthSq( x )
      real*8, intent(in) :: x
      DimenInvLengthSq = DimenInvLength( DimenInvLength( x) )
    end function DimenInvLengthSq


     pure real*8 function DimenStiffness(x)
       real*8, intent(in) :: x

       DimenStiffness = DimenInvLength(DimenForce(x))
     end function DimenStiffness


     pure real*8 function DimenEnergy(x)
       real*8, intent(in) :: x
       DimenEnergy = DimenForce(DimenLength(x))
     end function DimenEnergy

     pure real*8 function DimenTime(x)
       real*8, intent(in) :: x
       DimenTime = x / omega_scale
     end function DimenTime

     pure real*8 function DimenLength(x)
       real*8, intent(in) :: x
       DimenLength = x * A_scale
     end function DimenLength

    function DimenLengthArray( x, n )
      integer, intent(in) :: n
      real*8, intent(in) :: x(n)
      real*8 :: DimenLengthArray(n)
      integer :: i
      do i = 1,n
         DimenLengthArray(i) = DimenLength(x(i))
      end do
    end function DimenLengthArray

     
     pure real*8 function DimenVelocity(x)
       real*8, intent(in) :: x
       DimenVelocity = DimenLength(x) * omega_scale   !oops! this was wrong until jun 23 2011
     end function DimenVelocity

     pure real*8 function DimenAcceleration(x)
       real*8, intent(in) :: x
       DimenAcceleration = DimenLength(x) * omega_scale * omega_scale
     end function DimenAcceleration

     pure real*8 function DimenFreq( x )
       real*8, intent(in) :: x
       DimenFreq = x * omega_scale
     end function DimenFreq

     !Forces: F (ND) = F (Newton) * omegai_nd(i)**2 /Ainitial/Keq(i)
     pure real*8 function DimenForce(x)
!       use data0
       real*8, intent(in) :: x      

       DimenForce = x * A_scale * omega_scale ** 2
     end function DimenForce

    pure real*8 function DimenPressure(x)
      real*8, intent(in) :: x      

      DimenPressure = DimenForce( DimenInvLengthSq(x))
    end function DimenPressure

     pure real*8 function UnscaledForce(x,i)
       use data0
       real*8, intent(in) :: x
       integer, intent(in) :: i

       UnscaledForce = x / omegai(i)**2 * Keq(i) 
     end function UnscaledForce

 end module NonDimensionalization


   module peakForceMod
     implicit none
     save
     real*8 :: Zbase_amp_nd, PeakForce, PeakForceA, PeakForce_output, PeakForceA_output
     real*8, private :: cos1_m1, cos1_m2
   contains

     subroutine Compute_Z_Position_PeakForce(  cos1, Z_Controller_On, Force,t )
       logical,intent(out) :: Z_Controller_On
       real*8, intent(in) :: cos1, Force,t

       PeakForce  = max( Force, PeakForce)
       PeakForceA = min( Force, PeakForceA)      

       !we may need to reset the PeakForce control variable before we output the value. so just a separate variable for output
       PeakForce_output = PeakForce
       PeakForceA_output = PeakForceA
       
       ! want the Z controller to be one for exactly one point at the top of the cycle.
       ! keep track of 3 successive values and look for peak
       if ( (cos1_m2 <= cos1_m1) .and. (cos1_m1 >= cos1)) then
          Z_Controller_On = .true.          
       else
          Z_Controller_On = .false.
       end if

       cos1_m2 = cos1_m1
       cos1_m1 = cos1
     end subroutine Compute_Z_Position_PeakForce

     
   end module peakForceMod

 
 module Controller
   use params
   implicit none
   save 

   logical Z_Controller_On
   real*8 :: KI, KP, ErrorZ, NoiseAmp, sample_freq_hz,  Ztrack, Z_error

   integer*8 :: pad
!removed the private scope for debugging
   real*8 :: Z00, setpoint, OffsetSetpoint, errordt !, Z01
   real*8, private :: approach_tol  !maybe make this a user inpuvt

   integer, private :: Z_feedback_choice

   integer, private :: Z_control_scheme
   integer, parameter :: Z_PI=1, Z_SCHED_PI=2
   integer, private :: sched, last_sched

   contains

     !this gets called when transient is over
     subroutine Center_Z_Controller(Z0)
       real*8, intent(in) :: Z0
       errordt = 0d0 !anti-windup.  
       Z00 = Z0 !fixme, better name
     end subroutine Center_Z_Controller


     subroutine Init_Z_Controller(Z0, setpoint_in, z_feedback_choice_in, HF)
       real*8, intent(in) :: setpoint_in, HF
       real*8, intent(inout) :: Z0
       integer, intent(in) :: z_feedback_choice_in
       errordt = 0d0  
       Z00 = Z0
       setpoint = setpoint_in 
       Z_feedback_choice = z_feedback_choice_in

       if ((z_feedback_choice == FREQ_SHIFT_Z) .or. (Z_feedback_choice == DRIVE_AMP_Z)  .or. (z_feedback_choice == DEFLECTION_Z)) then
          approach_tol = 0.01 !maybe make this a user input
       elseif (z_feedback_choice == MEAN_DEFL_Z) then
          approach_tol = 0.1 
       elseif (z_feedback_choice == MAX_FORCE) then
          approach_tol = 0.02
       elseif  (z_feedback_choice == AMPLITUDE_Z) then
          if (( HF > 0.5) .or. (HF == 0)) then
             approach_tol = 0.015 !maybe make this a user input
          else
             !feature is small.  a small offset will be more noticeable so take the time to get closer
             approach_tol = 0.005
          end if
       elseif (z_feedback_choice == PHASE_Z ) then
          approach_tol = 0.02
       else
          call assert(.false., 'unhandled case3')
       end if

     end subroutine Init_Z_Controller
          

     !we are assuming that the piezotube can respond essentially instantaneously to a controller command.
     !this is okay IF the bandwidth of the piezotube is much greater thabn the bandwidth of the lock-in.
     !may need to revisit this in the future.
     subroutine Z_Controller(AppFeatH, ActualFeatH, delt, IOUT, pad, Z0, modulation_type)
       real*8, intent(in) :: AppFeatH, ActualFeatH, delt
       real*8, intent(inout) :: Z0
       integer*8, intent(in) :: IOUT, pad
       integer, intent(in) :: modulation_type

 !      real*8, parameter :: up_scale = 0.5d0, dn_scale = 3d0
       real*8, parameter :: up_scale = 1d0, dn_scale = 1d0
 !      real*8, parameter :: up_scale = 2d0, dn_scale = 8d0
       if (modulation_type /= PEAK_FORCE) then
          if (mod(IOUT,pad) == 0) then         
          !delt is non-dimensional.  So really, what we ought to be doing is either non-dimensionalizing the controller KI gain,
          !OR, computing this integral with dimensional time.  However, that means we'd have to change a bunch of things
          !and the gains that everybody has gotten used to using might not work anymore.  So for now, this mean that if you change
          !omega_scale (e.g. by changing drive frequency) then you need to change your integral gain as well.  Daniel, feb 2012.
             errordt=errordt+Z_error*delt*pad
             
             Z_control_scheme = Z_PI

             if (Z_control_scheme ==  Z_PI) then
                Ztrack =  -Kp*Z_error-KI*errordt !this is output as "topog".  not used anywhere else
             elseif (Z_control_scheme ==  Z_SCHED_PI) then
             !I built this complicated gain scheduling algorithm, but it hasn't been published
             !yet.  this replicates sorta what some of the japanese guys have been doing
             !(I think Ando has something like this)
                if (Z_error < -0.9) then
                   sched = 1
                   !going up a step               
                elseif (Z_error > 0.9) then
                   sched = 2
                   !going down a step            
                else
                   sched = 3
                end if
                

                if (sched == 1) then
                   if (last_sched == 3) then
 !                  errordt = errordt / up_scale
                   elseif (last_sched == 2) then
 !                  errordt = errordt / up_scale * dn_scale
                      errordt = errordt * dn_scale
                   end if
                   
                   Ztrack =  -(Kp*up_scale)*Z_error-(KI)*errordt
                
                elseif (sched == 2) then
                   if (last_sched == 1) then
 !                  errordt = errordt / dn_scale * up_scale
                      errordt = errordt / dn_scale
                   elseif (last_sched == 3) then
                      errordt = errordt / dn_scale
                   end if

                   Ztrack =  -Kp*Z_error-(KI*dn_scale)*errordt

                elseif (sched == 3) then
                   if (last_sched == 1) then
 !                  errordt = errordt * up_scale
                   else if (last_sched == 2) then
                      errordt = errordt * dn_scale
                   end if
                   Ztrack =  -Kp*Z_error-KI*errordt
                end if
             
                last_sched = sched
             end if
             
             ErrorZ=   Ztrack - ActualFeatH 
             
             Z0 = Z00 + Ztrack - AppFeatH
          end if
          
       else
          !Peak Force
             errordt = errordt+Z_error*delt*pad         
             Z_control_scheme = Z_PI
             Ztrack =  -Kp*Z_error - KI*errordt
             ErrorZ =   Ztrack - ActualFeatH              
             Z0 = Z00 + Ztrack - AppFeatH
       end if


     end subroutine Z_Controller
     

      !this is used for freq mod.  don't know the exact frequency setpoint at the begining.
      !have to figure it out as we go. this eliminates omegaR
      subroutine Offset_setpoint( offset)
        real*8, intent(in) :: offset
        if (z_feedback_choice == FREQ_SHIFT_Z) then !fixme. this is a hack
           OffsetSetpoint =  offset + setpoint
        end if
      end subroutine Offset_setpoint

      !all of the information about the Z error is localized in this function so that we can
      !easily create a new imaging mode (such as feedback on the 1/4 harmonic) without changing
      !too much. 
      subroutine compute_Z_error( AmpC, omegad, Drive_signal, u, LockinHH1_R, phaseC, PeakForce)
        use Nondimensionalization, only : DimenForce
        real*8, intent(in) :: AmpC, omegad, Drive_signal, u, LockinHH1_R, phaseC, PeakForce

        if (Z_feedback_choice == AMPLITUDE_Z) then
           Z_error = AmpC - setpoint
        elseif (Z_feedback_choice == FREQ_SHIFT_Z) then
           Z_error = (omegad - OffsetSetpoint) 
        elseif (Z_feedback_choice == DRIVE_AMP_Z) then
           Z_error =  setpoint - Drive_signal
        elseif (Z_feedback_choice == PHASE_Z) then
           Z_error = setpoint - phaseC
        elseif (Z_feedback_choice == DEFLECTION_Z) then
           Z_error =  setpoint - u
           !fixme? low pass filtering on u for contact mode scan?
        elseif (Z_feedback_choice == MEAN_DEFL_Z) then
           Z_error = setpoint - LockinHH1_R          
        elseif (Z_feedback_choice == MAX_FORCE) then
           Z_Error = setpoint - DimenForce(PeakForce)*1d9
        else
           call assert(.false., 'unhandled case')
           Z_error = 0
        end if

      end subroutine compute_Z_error

      pure real*8 function ApproachPercent()
        ApproachPercent = 100d0 * ( 1  - abs(Z_error / setpoint))
      end function ApproachPercent     

      pure logical function isSampleEngageWithinTol()
        isSampleEngageWithinTol = (( abs(Z_error / setpoint)) < approach_tol)
        
      end function isSampleEngageWithinTol

    end module Controller


    module MovingAverageFilter
    contains
      !a bit of a hack. needed for forceViewer to create the rounded triangle profile.  doesn't touch
      !the begining or end of the data
      !is non-causual.  filters forward and then backward to avoid any phase delay
      function  doMovingAverageFilter( x, N, len) result( out)
        integer :: i
        integer, intent(in) :: len,N
        real*8, intent(in) :: x(N)
        real*8, dimension(N) :: out
        
        out(1:len) = x(1:len)
        do i = (1+len), (N-len)
           out(i) = sum(  x( (i-len):(i+len) ) ) / (2d0 * len+1)
        end do
        out( (N-len):N) = x( (N-len):N)                   
      end function doMovingAverageFilter
    end module MovingAverageFilter
    
   !used by FM, for matching specific other simulators and specific exp hardware
   !right now, just one instance, but could be generalized.
   !this is a 1 pole pair butterworth IIR filter, discretized by the bilinear
   !transform.   e.g. what you get in matlab from butter(1,w).
   !the center frequency is 1.  that makes the math easier, but if we decide to 
   !change omega_scale, this gets all messed (sorry!)
   module BandPassFilter
     integer, parameter, private :: N = 2
     real*8, private ::  a(0:N), b(0:N), x(0:N), y(0:N), Fs
     real*8, parameter :: Q = 2d0
     contains
     subroutine InitBPF( Fs_in )
       real*8, intent(in) :: Fs_in
       real*8 :: norm
       integer :: i

       Fs = Fs_in

       !1 pole pair filter
       b(0) = -2d0 * Fs
       b(1) = 0
       b(2) = 2d0 * Fs
       a(0) = Q + 4d0 * Q * Fs**2 - 2d0 * Fs
       a(1) = 2d0 * Q - 8d0 * Q * Fs**2
       a(2) =  4d0 * Q * Fs**2 + 2d0 * Fs + Q    
       ! 2 pole pair filter
       ! b(0) = 4d0 * Fs**2
       ! b(1) = 0d0
       ! b(2) = -8d0 * Fs**2
       ! b(3) = 0d0
       ! b(4) = 4d0 * Fs**2
       ! a(0) = 16d0*Q**2*Fs**4+4*Fs**2-8*sqrt(2d0)*Q*Fs**3+Q**2+8*Q**2*Fs**2-2*sqrt(2d0)*Q*Fs
       ! a(1) = 16d0*sqrt(2d0)*Q*Fs**3-4*sqrt(2d0)*Q*Fs+4*Q**2-64*Q**2*Fs**4
       ! a(2) = -8d0*Fs**2+96*Q**2*Fs**4-16*Q**2*Fs**2+6*Q**2
       ! a(3) = -64d0*Q**2*Fs**4+4*sqrt(2d0)*Q*Fs-16*sqrt(2d0)*Q*Fs**3+4*Q**2
       ! a(4) = 8d0*sqrt(2d0)*Q*Fs**3+16*Q**2*Fs**4+2*sqrt(2d0)*Q*Fs+Q**2+4*Fs**2+8*Q**2*Fs**2

       norm = a(N)
       do i = 0,N
          a(i) = a(i) / norm
          b(i) = b(i) / norm
       end do
       x = 0
       y = 0      

     end subroutine InitBPF

     complex*8 function bp_filter_gain( omegad)
       use params
       real*8, intent(in) :: omegad
       complex*8 :: num, den, z, j
       integer :: i
       num = 0
       den = 0
       j = CMPLX(0,1)

       z = exp( j * (1+omegad) / (Fs)  ) !do we need any 2 pis in here?
       do i = 0,N
          num = num + b(i) * (z**i)
          den = den + a(i) * (z**i)
       end do
       bp_filter_gain = num / den  
     end function bp_filter_gain
     
     real*8 function BPF(x_in)
       real*8, intent(in) :: x_in
       integer :: i      
       do i = 0,(N-1)
          x(i) = x(i+1)
          y(i) = y(i+1)
       end do
       x(N) = x_in
       y(N) = dot_product(b, x) -  dot_product(a(0:(N-1)), y(0:(N-1)))
       BPF = y(N)
     end function BPF

   end module BandPassFilter


  module LockIn
     use params
     use LockInDataType
     use BandPassFilter

     real*8 :: LockInTC  !for now all lockins have same tc and order.  this is only b/c
     !I don't feel like writing the i/o to read in a separate tc.  can be arbitary if you want.

     integer :: LockInOrder ! 0 = don't use lock in. 1st, 2nd, or 4th order filter is supported

     logical :: want_pre_BPF

     type(LockInData) :: MainLockin
     type(LockInData) :: BimodalLockin
     type(LockInData), allocatable :: HigherHarmLockin(:)
     type(LockInData) :: RMS  !no its not a lockin, but it was easier to hack it in here than write a separate module

     contains

       complex*8 function LockinGain( omegad )
         real*8, intent(in) :: omegad
         
         select case(LockInOrder)
         case (1)
            LockinGain = 1d0 / CMPLX(1,  omegad * LockinTC)
         case (2)
            LockinGain = 1d0 / (CMPLX( sqrt2 / 2d0, sqrt2/2d0 +  omegad * LockinTC) *CMPLX( sqrt2/2d0, -  sqrt2/2d0 + omegad * LockinTC)) !this appears to be correct
         case (4)
            LockinGain = 1d0 / (CMPLX( butterworth4_1/2d0, butterworth4_2/2d0 +  omegad * LockinTC) *CMPLX( butterworth4_1/2d0, -butterworth4_2/2d0 + omegad * LockinTC) * CMPLX( butterworth4_2/2d0, butterworth4_1/2d0 +  omegad * LockinTC) * CMPLX(butterworth4_2/2d0, -butterworth4_1/2d0 + omegad * LockinTC))
         case default
         end select
         
         if ( want_pre_BPF) then
            LockinGain = LockinGain * bp_filter_gain( omegad)
         end if
         
       end function LockinGain

       subroutine ComputeLockin(data, delt, u_in, Noise, t, theta_in, want_rms_in)

         type(LockInData), intent(inout) :: data
         real*8, intent(in) :: delt, u_in, Noise, t, theta_in
         real*8 :: u, ux, uy
         logical, optional :: want_rms_in
         logical :: want_rms

         if (.not. present(want_rms_in)) then
            want_rms = .false.
         else
            want_rms = want_rms_in
         end if

         if (want_pre_BPF) then
            u = BPF(u_in)
         else
            u = u_in
         end if

         if (want_rms) then
            !this is a bit of a misnomer to put the rms computation into the lockin modulues.
            !but we want to re-use the filter code.  so just hack it in here
            ux = (u+Noise)**2
            uy = 0
         else
            !previously these two lines were after the filter, but I don't know why.  this seems to make more sense.
            !also no reason for them to be part of the type data.
            ux = 2d0 * (u+Noise) * cos(theta_in)
            uy = 2d0 * (u+Noise) * sin(theta_in)
         end if

         ! butterworth filters.  use buttap in matlab to get the poles, multiply complex conjugates to get quadratics
         ! then convert to second order system chained together
         if ( LockInOrder == 1) then
            data%LockInX = data%LockInX_Old + delt * ( ux - data%LockInX_Old)/LockInTC
            data%LockInY = data%LockInY_Old + delt * ( uy - data%LockInY_Old)/LockInTC

            data%LockInX_Old = data%LockInX
            data%LockInY_Old = data%LockInY

         elseif (LockInOrder == 2) then
            data%LockInX   = data%LockInX_Old   + delt * data%LockInXpr_old
            data%LockInXpr = data%LockInXpr_Old + delt * ( ux - data%LockInX - sqrt2 * LockInTc * data%LockInXpr_Old) / LockInTC ** 2
            data%LockInY   = data%LockInY_Old   + delt * data%LockInYpr_old
            data%LockInYpr = data%LockInYpr_Old + delt * ( uy - data%LockInY - sqrt2 * LockInTc * data%LockInYpr_Old) / LockInTC ** 2
            data%LockInXpr_Old = data%LockInXpr
            data%LockInYpr_Old = data%LockInYpr
            data%LockInX_Old = data%LockInX
            data%LockInY_Old = data%LockInY
         elseif (LockInOrder == 4) then
            !4th order is like two 2nd orders chained together
            data%LockInX   = data%LockInX_Old   + delt * data%LockInXpr_old
            data%LockInXpr = data%LockInXpr_Old + delt * ( data%LockIn2X_Old - data%LockInX - butterworth4_1  * LockInTc * data%LockInXpr_Old) / LockInTC ** 2
            data%LockInY   = data%LockInY_Old   + delt * data%LockInYpr_old
            data%LockInYpr = data%LockInYpr_Old + delt * ( data%LockIn2Y_Old - data%LockInY - butterworth4_1 * LockInTc * data%LockInYpr_Old) / LockInTC ** 2
            data%LockInXpr_Old = data%LockInXpr
            data%LockInYpr_Old = data%LockInYpr
            data%LockInX_Old = data%LockInX
            data%LockInY_Old = data%LockInY

            data%LockIn2X   = data%LockIn2X_Old   + delt * data%LockIn2Xpr_old
            data%LockIn2Xpr = data%LockIn2Xpr_Old + delt * ( ux - data%LockIn2X - butterworth4_2 * LockInTc * data%LockIn2Xpr_Old) / LockInTC ** 2
            data%LockIn2Y   = data%LockIn2Y_Old   + delt * data%LockIn2Ypr_old
            data%LockIn2Ypr = data%LockIn2Ypr_Old + delt * ( uy - data%LockIn2Y - butterworth4_2 * LockInTc * data%LockIn2Ypr_Old) / LockInTC ** 2
            data%LockIn2Xpr_Old = data%LockIn2Xpr
            data%LockIn2Ypr_Old = data%LockIn2Ypr
            data%LockIn2X_Old = data%LockIn2X
            data%LockIn2Y_Old = data%LockIn2Y           
         else
            call assert(.false., "unknown filter order");
         end if

         data%Th = atan2( data%LockInX, data%LockInY) !this is cos / sin                  

         data%R  = sqrt( data%LockInX**2 + data%LockInY**2)
       end subroutine ComputeLockin

       subroutine DoAdditionalLockins( t, delt, u, numHH, NHH, exc_choice, WantHH, useLockin, theta1, theta2, want_RMS)
         use params

         integer, intent(in) :: numHH, exc_choice
         real*8, intent(in) :: t,delt, u, NHH(numHH), theta1, theta2
         logical, intent(in) :: WantHH, useLockin, want_RMS
         integer ij
         real*8 u_adj

         if (exc_choice == BIMODAL ) then
            if (useLockin) call ComputeLockin( BimodalLockin, delt, u, 0d0, t, theta2)
         end if

         if (WantHH) then
            do ij = 1, numHH
               if (NHH(ij) == 0) then
                  u_adj = u / 2
               else
                  u_adj = u
               end if

               if (useLockin) call ComputeLockin( HigherHarmLockin( ij), delt, u_adj, 0d0, t,  NHH(ij) * theta1)
            end do
         end if

         if (want_RMS) then
            call ComputeLockin( RMS, delt, u, 0d0, t, 0d0, .true.)
         end if
       end subroutine



       !instead of doing a whole fourier transform, we just do the integrals over the 
       !two or three harmonics we are intersted in.
       !this sub does a running sum of real and imaginary parts.

       !no, this isn't really a lockin, but we want it in a module somewhere so that fortran
       !will behave itself
       subroutine DoMainFourierIntegrals( t, delt, u, numHH, NHH, exc_choice, WantHH, Want_A2, Want_p2, WantFourier, a1, b1, a2, b2, an, bn,Force, F1r, F1i, want_ForceFourier, theta1, theta2)
         use params

         integer, intent(in) :: numHH, exc_choice
         real*8, intent(inout) :: a1, a2, b1, b2, an(numHH), bn(numHH), F1r(1), F1i(1)
         real*8, intent(in) :: t,delt, u, NHH(numHH), theta1, theta2, Force
         logical, intent(in) :: WantHH, WantFourier, Want_a2, want_p2, Want_ForceFourier
         integer ij
         real*8 u_adj

         if (wantFourier) then
            a1=a1+u*cos( theta1 )*delt 
            b1=b1+u*sin( theta1 )*delt           
         end if

         if (exc_choice == BIMODAL ) then
            if (wantFourier .and. (Want_A2 .or. Want_P2)) then
               a2=a2+u*cos( theta2 ) * delt
               b2=b2+u*sin( theta2 ) * delt		
            end if
         end if

         if (WantHH) then
            do ij = 1, numHH

               if (NHH(ij) == 0) then
                  u_adj = u / 2
               else
                  u_adj = u
               end if

               if (WantFourier) then
                  an(ij) = an(ij) + u_adj * cos( NHH(ij) * theta1 ) * delt
                  bn(ij) = bn(ij) + u_adj * sin( NHH(ij) * theta1 ) * delt
               end if

            end do
         end if

         if (want_ForceFourier) then
            F1r(1)=F1r(1)+Force*cos( theta1 )*delt 
            F1i(1)=F1i(1)+Force*sin( theta1 )*delt   
         end if
       end subroutine DoMainFourierIntegrals


   end module LockIn


   module FreqMod
     use params
     use LockIn
     implicit none
     !freq mod parameters.
     !time constants... keep in mind the non-dimensionalization by omega_scale
     ! we may want to scale based on the normalization
     real*8 :: fm_pll_kp, fm_pll_ki, fm_pll_kd, fm_amp_kp, fm_amp_ki, fm_amp_kd,  fm_initial_phase, fm_phase_error
     real*8 :: Drive_signal, freq_shift_sp, fm_amp_error, fm_phase_error_dt,  int_amp_error
     real*8, private :: prev_amp_error, prev_phase_error
     logical :: PLL_On, Amp_Controller_On, fm_direct_control, FM_want_noncontact
     logical :: WantCalcPLLGains,WantCalcAmpGains

     save

     contains

       real*8 function PLL_controller(MainLockin, IOUT, pad, t, delt)

         integer*8, intent(in) :: IOUT, pad
         real*8, intent(in) :: t, delt
         type(LockInData), intent(in) :: MainLockin
         real*8 :: new_omegad, phase_error_d

         if (fm_direct_control) then
            fm_phase_error = MainLockin%Th - fm_initial_phase           
         else
            fm_phase_error = MainLockin%R * sin( MainLockin%Th - fm_initial_phase)
         end if       

         if (mod(IOUT, pad)==0) then 
            !all three controllers run at same frequency for now.  many want to change in future

            !phase locked loop controller
            if (PLL_On) then                            

               phase_error_d = (fm_phase_error - prev_phase_error) / ( delt * dble(pad) )
               prev_phase_error = fm_phase_error

               fm_phase_error_dt = fm_phase_error_dt + fm_phase_error * delt * dble(pad) 

               PLL_Controller = 1d0 + fm_pll_kp * fm_phase_error + fm_pll_ki * fm_phase_error_dt + fm_pll_kd * phase_error_d

               if ( PLL_Controller < 0) then
                  call WriteFatalError('Drive frequency has become negative.  Check feedback gains')
               end if

            end if
         end if
       end function PLL_controller

       subroutine Amp_controller(MainLockin, IOUT, pad, delt)
         integer*8, intent(in) :: IOUT, pad
         type(LockInData), intent(in) :: MainLockin
         real*8, intent(in) :: delt
         real*8 :: amp_error_d

         !remember, amplitude is non-dimensionalized by initial amplitude
         if (fm_direct_control) then
            fm_amp_error = 1d0 - MainLockin%R 
         else
            fm_amp_error = 1d0 - MainLockin%R * cos( MainLockin%Th - fm_initial_phase)
         end if

         !all three controllers run at same frequency for now.  many want to change in future
         if (mod(IOUT, pad)==0) then

            if (Amp_Controller_On) then           
               !!previous version.  had forgotten to normalize by pad.  answer would change significantly
               !!if we change sampling rate.
               !int_amp_error = int_amp_error + fm_amp_error * delt

               !new version.  normalized by pad.  trapezoidal integration
               int_amp_error = int_amp_error + ((fm_amp_error + prev_amp_error) /2d0)  * delt * dble(pad)

               amp_error_d = (fm_amp_error - prev_amp_error) / ( delt * dble(pad) )

               Drive_signal = 1d0 + fm_amp_kp * fm_amp_error + fm_amp_ki * int_amp_error + fm_amp_kd * amp_error_d

               !with the new autocalculation of feedback gains, this is less necessary, and also is unnecessarily restrictive on the transient when the amp controller is first turned on.
 !              if (Drive_signal < 0d0) call WriteFatalError('Drive amplitude has become negative.  Check feedback gains')
               if (Drive_signal < 0d0) Drive_signal = 0d0
            end if

            prev_amp_error = fm_amp_error
         end if

       end subroutine Amp_controller

       !the transfer function for the slow time scale dynamics of the amplitude
     complex*8 function AmpXferFunc(omegad, Q, omegai_nd)
       real*8, intent(in) :: omegad, Q, omegai_nd
       AmpXferFunc = LockinGain(omegad) * ((omegai_nd / (2 * Q)) / CMPLX( omegai_nd / (2 * Q), omegad))       
     end function AmpXferFunc

     !see kilpatrick et al, RSI, 2009. we could maybe do better than this if we spent the time
     !to work it out, but this is not a bad start.
     subroutine CalculateFMGains(LockinOrder, LockinTc, Q, omegai_nd, fexcite, osc_Q, osc_omega_nd, want_pi)
       use ZiegerNichols
       real*8, intent(in) :: LockinTc, Q, omegai_nd, osc_Q, osc_omega_nd
       integer, intent(in) :: LockinOrder, fexcite
       real*8 :: a, b, guess, omegad_180, k_ultimate_amp, k_ultimate_pll, tau_ultimate
       logical, intent(in) :: want_pi

       complex*8 :: foo


      !start with amplitude controller
      !step 1, find the frequency at which the angle of the open loop
      !transfer function is -180 degree.  call it omega_180

      !simple bisection root finding
      a = 0d0;
      call assert(angle(AmpXferFunc(a,Q, omegai_nd)) > -Pi, 'a')
      call assert(CABS(AmpXferFunc(a,Q, omegai_nd) - 1e0) < 0.001, 'amp trans func != 1 at DC')


      select case (LockinOrder)
      case (1)
         !overall xfer function hits -180 at omega=infinity.  will go unstable due 
         !to numerical issues, or discrete time sampling, or non-linearities first.  
         !can't predict with this method.... although, we could if we use the 
         !the pre-bandpass filter
         call WriteFaTalError('Cannot autocalculate gains with 1st order lockin')
      case (2)
         !will wrap once at hit -270 (=+90) at infinity.  we can easily
         !find the -180 point
         b = max(10d0, 10d0/ LockinTC)
         !            call assert(angle( AmpXferFunc(b,Q, omegai_nd)) < -Pi, 'b')
      case(4)            
         !overall xfer will wrap again and hit -540 (=-90) at infinity. multi-valued
         !but lockin itself will only be single valued, so don't start out as far.
         b = 1d0 / LockinTC
      end select


      do while( (b -a) > 0.001)
         guess = (a+b)/2;         
         !essentially, we're looking for the wrap around point
         if ( angle(AmpXferFunc(guess,Q, omegai_nd))  > 0 ) then
            b = guess
         else
            a = guess
         end if

      end do
      omegad_180 = (a+b)/2

      !step 2, ultimate proportional gain is that which makes transfer function
      !have unity magnitude at omega_180.
      k_ultimate_amp = 1 / CABS( AmpXferFunc(omegad_180,Q, omegai_nd))
      !step 3, oscillation period at ultimate gain. previously had forgotten the 2 pi
      tau_ultimate = (2d0 * pi) / omegad_180 


      if (wantCalcAmpGains) then
         !step 4, calculate Ki and Kp from ultimate gain and period.
         !Zieger-Nichols.

         if (want_pi) then
            call ZN_PI( k_ultimate_amp, tau_ultimate, fm_amp_kp, fm_amp_ki, fm_amp_kd)
         else
            call ZN_PID( k_ultimate_amp, tau_ultimate, fm_amp_kp, fm_amp_ki, fm_amp_kd)
         end if

      end if

write (*,*) (k_ultimate_amp / 2.2), (1.2 * k_ultimate_amp / 2.2 / tau_ultimate)

      if (wantCalcPLLGains) then
       !PLL.  transfer function is actually the same, just scaled differently.
       !in particular, need the slope of the phase at natural frequency, which is 2 * Q / omega_i
       !at least for magnetic anyway.  may need to generalized for acoustic
       k_ultimate_pll = k_ultimate_amp / ( 2 * Q / omegai_nd )

       if (fexcite == ACOUSTIC_PIEZOR) then
          !actually, this is only if driving freq is close to the piezo resonance.
          !off resonance the slope is less.  so this is worst case I guess.
          k_ultimate_pll = k_ultimate_pll / ( 2 * osc_Q / osc_omega_nd)
       end if

       if (want_pi) then
          call ZN_PI( k_ultimate_pll, tau_ultimate, fm_pll_kp, fm_pll_ki, fm_pll_kd)
       else
          call ZN_PID( k_ultimate_pll, tau_ultimate, fm_pll_kp, fm_pll_ki, fm_pll_kd)
       end if

    end if
     end subroutine CalculateFMGains

     pure real*8 function angle( g)
       complex*8, intent(in) :: g
       angle = atan2( IMAG(g), REAL(g))
     end function angle

     subroutine FreqModInit()
       PLL_On = .false.
       Amp_Controller_On = .false.
       int_amp_error = 0
       fm_phase_error_dt = 0
       Drive_Signal = 1
     end subroutine FreqModInit

   end module FreqMod

   module timeAndCycle
     use params
     implicit none
     save

     integer :: plotpnts, exc_choice, sweepchoice, operating_mode, numincycle
     integer :: modulation_type, hist_cycles
     real*8  :: delt, maxModalRingdownTime, sweep_time, sweep_rate,  SimulationDuration
     logical :: isTransientOver, want_freqswp_sp_aprch, transient_suppress_Fts

     logical :: useLockIn ! false = old style fourier integrals, true = filtering lock-in

     integer :: current_stairstep_point  ! used for frequency sweeps only

     integer,   private :: transient_phase
     integer*8, private :: ntrans
     real*8,    private :: last_theta, last_time, prev_omegad
     real*8,    private :: cur_omegad(maxExc), ttrans
     real*8,    private :: omegad_start, omegad_stop ! used for frequency sweeps only

     real*8, private :: transient_timeout 

     contains

       logical function haveFeedbackController()
         haveFeedbackController = ((modulation_type == FREQUENCY) .or. (operating_mode == SCAN) .or. ((operating_mode == FREQSWEEP) .and. (want_freqswp_sp_aprch)))
       end function haveFeedbackController

       subroutine SetInitialTransientPhase( phase_in)
         integer :: phase_in
         transient_phase = phase_in
         transient_suppress_Fts = .false. !used only by FM
       end subroutine SetInitialTransientPhase


     !this function is called at every time step when transients are on.  It handles all of the operation of the transient
     !calculations, including sample engagement (if applicable), turning on the PLL and amplitude controllers (if applicable)
     !and letting numerical transients (settle).  the output variable 
     !this routine changes the variable isTransientOver to signal if transients should continue on the next time step or not
       subroutine TransientStateMachine(t, Drive_signal,u, Quality, omegai_nd, fexcite)
         use Controller, only: Z_Controller_On, isSampleEngageWithinTol, Offset_setpoint
         use LockIn, only: LockinOrder, LockInTC
         use FreqMod, only: Amp_Controller_On, PLL_On, fm_amp_kp, fm_amp_ki, CalculateFmGains, wantCalcAmpGains
         use Approaching, only: ApproachStart
         use NonDimensionalization

         real*8, intent(in) :: t, Drive_signal,u, Quality, omegai_nd       
         integer, intent(in) :: fexcite       
         character*1000 errStr
         
         if (((modulation_type == AMPLITUDE) .or. (modulation_type == PEAK_FORCE)) .and. (operating_mode == SCAN)) then
            select case (transient_phase)
               !transient is over when we get to within the tolerance AND stay there for some period of time
            case (TRANS0)  !let lock-in filter settle
               if ( t > 20d0 * LockInTC) then
                  transient_phase = TRANS1
                  Z_Controller_On = .true.
               end if
               isTransientOver = .false.
            case (TRANS1)            
               if (isSampleEngageWithinTol())  then 
                  transient_phase = TRANS2
                  if (modulation_type /= PEAK_FORCE) then
                     ttrans = t + max(50d0 * LockInTC,2d0 * 2d0*Quality / omegai_nd)
                  else
                     !maybe make this a user input?
                     !                     ttrans = t +  (100d0*2d0*Quality / omegai_nd)
                     ttrans = t +  (50d0*2d0*Quality / omegai_nd)
                  end if
               end if
               if ( t > transient_timeout) then
                  if (fexcite == NO_EXC) then
                     !contact mode
                     call WriteFatalError(  "Controller could not stabilize within alloted time." &            
                          // char(10) // "If the transient perctage is oscillating rapidly, you should decrease the gains." &
                          // char(10) // "If it appears stuck at a fixed value, you may need to increase the gains." & 
                          // char(10) // "Check examples for good values." &
                          // char(10) // "You can also choose to view a time history of the transients by checking the box in the simulation parameters tab to help debug.")
                  else
                     !am
                     errStr = "Controller could not stabilize within alloted time." &
                          // char(10) // "Suggestions: Check approach curve to verify that chosen setpoint is not in a bistable region." &
                          // char(10) // "If the transient perctage is oscillating rapidly, you may need to decrease the gains." &
                          // char(10) // "If it appears stuck at a fixed value, you may need to increase the gains." & 
                          // char(10) // "Also, check that lock-in time constant is not too small (bandwidth less than 1/10 cantilever resonance may be too small)." // char(10)

                     if (LockinOrder == 1) then
                        !the zieger nichols fails for 1st order lockin, so we cant suggest values.  probably not a great idea anyway
                        errStr = trim(errStr) // "you might have better luck with a 2nd or 4th order lockin rather than a 1st order lockin"
                     else
                        ! for hard surfaces, 1 nm in Z is 1 nm in amp, so actually ziegner nichols we use for the amp cont in FM
                        ! is a good guess (just scaled differently).  for soft surfaces that guesses too low
                        !this is a bit of a hack calling the FM routine.  should really write a more general purpose one.
                        wantCalcAmpGains = .true.
                        call CalculateFMGains(LockinOrder, LockinTc, Quality, omegai_nd, fexcite, 0d0, 0d0, .false.)
                        fm_amp_kp = fm_amp_kp / (2 * Quality / omegai_nd)
                        fm_amp_ki = fm_amp_ki / (2 * Quality / omegai_nd)
                        write(errStr,*) trim(errStr), "If you are using a hard sample (e.g. E > 10 GPa), good starting values for your parameters might be ",char(10), " kp = ", fm_amp_kp, " and ki = ", fm_amp_ki, char(10), "if you are using a softer sample, you might need to increase these a bit."
                     end if

                     errStr = trim( errStr) // char(10) // "Users of the advanced scanning tool can choose to view a time history of the transients by checking the box in the simulation parameters tab to help debug." &
                          // char(10) // "If you using the force modulation tool, make the excitation amplitude on the hard surface smaller."

                     call WriteFatalError(errStr)
                  end if
               end if
               !end if

               isTransientOver = .false.
            case (TRANS2)
               isTransientOver = .false.

               if (.not. isSampleEngageWithinTol())  then 
                  transient_phase = TRANS1 !try again!   
               else if (t > ttrans) then
                  isTransientOver = .true.                  
               end if               
            end select

         elseif (modulation_type == FREQUENCY) then
            !similar to am scanning, but we have to turn on multiple controllers and let
            !each stabilize
            isTransientOver = .false.
            select case( transient_phase)
            case (TRANS_NUM)
               transient_suppress_Fts = .true.  !need to calculate phase setpoint w/o tip-sample forces
               if ( t > ttrans) then
                  !               if (debugging) write(6,*) "turn on Fts at time ", DimenTime(t)
                  transient_phase = TRANS_FTS
                  ttrans = t + max( 20d0 * LockInTC, maxModalRingdownTime)
                  transient_suppress_Fts = .false. !now turn Fts back on
               end if
            case (TRANS_FTS)
               if ( t > ttrans) then
                  !               if (debugging) write(6,*) "turn on PLL at time ", DimenTime(t)
                  transient_phase = TRANS2_PLL
                  ttrans = t + max( 20d0 * LockInTC, maxModalRingdownTime)
                  PLL_On = .true.
               end if
            case (TRANS2_PLL)
               if ( t > ttrans) then
                  !               if (debugging) write(6,*) "turn on amp controller at time", DimenTime(t)

                  transient_phase = TRANS3_AMP
                  ttrans = t + max( 20d0 * LockInTC, maxModalRingdownTime)
                  Amp_Controller_On = .true.
               end if
            case(TRANS3_AMP)
               if ( t > ttrans) then

                  call Offset_setpoint( omegad(t,1) ) !this used to be omegaR

                  if (operating_mode == SCAN) then                  
                     !                  if (debugging) write(6,*) "turn on Z controller at time", t
                     transient_phase = TRANS4_Z_On
                     Z_Controller_On = .true.
                  else
                     call ApproachStart()
                     isTransientOver = .true.                  
                  end if
               end if
            case (TRANS4_Z_ON)
               if ( isSampleEngageWithinTol()) then 
                  !               if (debugging) write(6,*) "in tol at time ", t
                  transient_phase = TRANS5_Z_STABL
                  ttrans = t + 26d0 * LockInTC
               end if

               if ( t > transient_timeout) then
                  call WriteFatalError(  "The Z controller could not stabilize within alloted time." &
                       // char(10) // "If the transient perctage is oscillating rapidly, you may need to decrease the gains." &
                       // char(10) // "If it appears stuck at a fixed value, you may need to increase the gains." & 
                       // char(10) // "Also, check that lock-in time constant is not too small (bandwidth less than 1/10 cantilever resonance may be too small).")
               end if

            case (TRANS5_Z_STABL) 
               if (t > ttrans) then
                  if ( isSampleEngageWithinTol()) then 
                     transient_phase = TRANS6_OP
                     isTransientOver = .true.
                  else
                     transient_phase = TRANS4_Z_ON
                  end if
               end if
            case default
               call assert( .false., 'invalid state')
            end select
         else if ( isOpModeApp(operating_mode) .or. (operating_mode == FIXED)) then
            !for these cases, we just sit for a fixed amount of time       
            if (t > ttrans) then
               isTransientOver = .true.
               call ApproachStart()
            else
               isTransientOver = .false.
            end if

         else if (modulation_type == FORCE_VOL) then          
            if (t > ttrans) then
               isTransientOver = .true.
            else
               isTransientOver = .false.
            end if

         else if (operating_mode == FREQSWEEP) then
            !approach at resonance, then switch to start frequency and stabilize there
            isTransientOver = .false.
            select case (transient_phase)
            case (TRANS1_CONT_OFF) 
               if ( t > (ttrans )) then
                  if (debugging) write(*,*) "turn z controller on at time" , t
                  transient_phase = TRANS2_CONT_ON
                  Z_Controller_On = .true. 
                  if (debugging) write(*,*) "z controller on at time ", t
               end if
            case (TRANS2_CONT_ON)
               if (isSampleEngageWithinTol()) then
                  if (debugging) write(*,*) "in tolerance at time" , t
                  ttrans = t + 4d0 * LockInTC
                  transient_phase = TRANS3_STABILIZE
               else
                  if ( t > transient_timeout) then
                     call WriteFatalError(  "Controller could not stabilize within alloted time." &
                          // char(10) // "Suggestions: Check approach curve to verify that chosen setpoint is not in a bistable region")
                  end if
               end if
            case (TRANS3_STABILIZE)
               if (t > ttrans) then
                  if (isSampleEngageWithinTol( )) then  
                     transient_phase = TRANS4_CONT_OFF
                     last_theta = theta(t,1)
                     ttrans = t + ttrans 
                     call set_omegad_startstop(omegad_start, omegad_stop)               
                     last_time = t
                     call updateDelt(t)
                     Z_Controller_On = .false.
                     !fixme, update timeall?
                  else
                     !not yet stable. try again
                     transient_phase = TRANS2_CONT_ON
                  end if
               end if
            case (TRANS4_CONT_OFF)
               if (t > ttrans) then
                  last_theta = theta(t,1)
                  transient_phase = OPERATE
                  isTransientOver = .true.
               end if
            case (OPERATE)
               call assert( .false., 'should not have gotten called in this state!')
            case default
               call assert( .false., 'invalid state')
            end select
         end if
       end subroutine TransientStateMachine

       real*8 function timeSinceTrans(t)
         real*8, intent(in) :: t
         timeSinceTrans = t - ttrans 
       end function timeSinceTrans

     integer*8 function pointsSinceTrans(IOUT)
       integer*8, intent(in) :: IOUT
       pointsSinceTrans = IOUT - ntrans
     end function pointsSinceTrans

     subroutine UpdateDelt(t)
       real*8, intent(in) :: t
       delt = 1.0d0/numincycle*(2.0*pi)/ omegad(t,1) 
       call assert(delt > 0, 'delt > 0')
     end subroutine UpdateDelt

     subroutine SetupFirstNonTransient(t, IOUT)
       real*8, intent(in) :: t
       integer*8, intent(in) :: IOUT
       ttrans = t
       ntrans = IOUT
     end subroutine SetupFirstNonTransient



     subroutine CalcTimeCycleInfo(omegai, numModes, output_point_rate, nstep,Zf, pad, &
                Z0, Quality,Vscan, SubsLen, sample_freq_hz,nstepFZ, TriggerDeflLim,&
                LockInTC, LockInOrder, transient_allowance, want_fourier, numofcycle, transient_timeout_in, fexcite)
       use NonDimensionalization
       use Approaching
       use checkFzCurves
       integer, intent(in) :: numModes, LockInOrder, fexcite
       real*8, intent(inout) :: transient_allowance
       real*8, intent(in) :: omegai(numModes), Zf,Z0,Vscan, SubsLen, &
                             sample_freq_hz, Quality(maxModes), TriggerDeflLim, LockInTC, transient_timeout_in
       integer*8, intent(out) :: output_point_rate, nstep, pad, nstepFZ, numofcycle
       logical, intent(in) :: want_fourier

       real*8 ttrans_dim, ttrans_nd, lockin_allowance
       integer*8 transcycle
       integer i

       !to make the math easy, we insist that the time step of integration be an integer 
       !multiple of the sampling frequency. that way the controller's sampling of the waveform is 
       !just a decimation.  this may cause us to round up numincycle

       if (haveFeedbackController()) then
          pad = ceiling( numincycle / (  NonDimenFreq(sample_freq_hz * 2d0 * pi)))
          numincycle = nint( 2 * pi * NonDimenFreq(sample_freq_hz) * pad)
          call assert( pad > 0, 'pad > 0');
          call assert( numincycle > 0, 'numincycle > 0 (pad)')
       end if

 !Daniel Kiracofe 10/08. finds the time constants for all modes and picks the longest. also correct 
 !to use nat freq not drive freq.  gives slightly different answers than before for bimodal cases 
 !b/c we might not start at the same portion of the 2nd exc drive cycle.  

       if (transient_allowance == 0) then
          !this is a hidden flag.  should really be a user input
          lockin_allowance = 0d0
       else
          lockin_allowance = 2.5
       end if

       !Based on my results I think that 3*Q would be okay.  this is the default.  currently
       !only fixed point gui allows to override this
       if (transient_allowance < 0) transient_allowance = 3

       nstepFZ = 1 !prevent errors when in contact scan mode
       numofcycle = 0 !avoid compiler warning

       ttrans = 0
       do i = 1,numModes
          ttrans_dim = 2.0*pi/omegai(i)*transient_allowance*Quality(i) !dimensional time for transient to die down
          ttrans_nd  = NonDimenTime(ttrans_dim)
          if (ttrans_nd > ttrans) ttrans = ttrans_nd
       end do

       maxModalRingdownTime = ttrans

       !if the lock-in would take longer than the eigenmodes, use that as the limiting case.
       !be sure to give higher filter orders more time to settle.
       if ( (lockin_allowance*2*pi*LockInTC * sqrt(real(LockInOrder)) ) > ttrans) ttrans = (lockin_allowance*2*pi*LockInTC * sqrt(real(LockInOrder)))     

       if (haveFeedbackController() .and. (modulation_type /= PEAK_FORCE) ) then
          if ( ttrans*NonDimenFreq(sample_freq_hz) > 400000) then
             call WriteFatalError('Sampling frequency is excessively high compared to the bandwidth of the eigenmodes and the lockin filters.  This will lead to numerical instability.  Reduce sampling rate')
          elseif ( ttrans*NonDimenFreq(sample_freq_hz) < 5) then
             call WriteFatalError('Sampling frequency is very low compared to the bandwidth of the eigenmodes and the lockin filters.  This will lead to poor controller performance.  Increase sampling rate')

          end if
       end if

       !self exc may need much longer to stabilize?
       if ((exc_choice == SELFEXC) .and. (operating_mode /= FIXED)) ttrans = ttrans * 7

       !this rounds up to an integer number of drive cycles
       transcycle = ceiling( ttrans * omegad(0d0, 1) / (2.0 * pi) )
       call assert(transcycle > 0, 'transcycle > 0')        

       !recalculate the time to exactly match an integer number of cycles
       ttrans = transcycle * 2d0 * pi / omegad(0d0, 1)
       call assert(ttrans > 0, 'ttrans > 0')

       ntrans= transcycle*numincycle

       call assert(ntrans > 0, 'ntrans > 0')

       !	SimulationDuration refers to nondimensional time
       if (operating_mode == APPROACH_STEP) then
          numofcycle = plotpnts * Quality(1) * 5  !for now just 5*Q1 cycles per point.  might want to make a user input
          SimulationDuration = numofcycle * 2d0 * pi / omegad(0d0,1)  !recalc for integer number of cycles
          nstepFZ = ceiling(real(numincycle*numofcycle)/plotpnts)
       elseif ((operating_mode ==APPR_RET).or.(operating_mode == APPROACH).or.(operating_mode == APPROACH_SINE) ) then                              
          if ( Z0 == Zf) call WriteFatalError( "Initial and final Z separations are the same.  No simulation to perform!")

          if (isFzCurves(operating_mode, fexcite)) then
             !fz curves, and untriggered at that.  
             SimulationDuration = NonDimenTime( 1 / CurveRate_dim)
             
             if (operating_mode ==  APPROACH) then
                !approach only
                SimulationDuration =  SimulationDuration / 2d0
                call SetApproachSpeedSetpoint( (Z0-Zf) / SimulationDuration )
             elseif (operating_mode ==APPR_RET) then
                !approach retract
                call SetApproachSpeedSetpoint( 2*(Z0-Zf) / SimulationDuration )
             elseif  (operating_mode == APPROACH_SINE) then
                !nothing to do here actually
             end if
                          
          else
             !dynamic approach curves of some kind. these are always untriggered, and always used a speed
             !probably should be consistent with Fz curves, but haven't gotten there yet
             SimulationDuration = (Z0-Zf)/GetApproachSpeedSetpoint() !first estimate
             if (operating_mode == APPR_RET) SimulationDuration =  2.0 * SimulationDuration
          end if          
          
          numofcycle=ceiling(SimulationDuration*omegad(0d0, 1)/2.0/pi) !round up to an integer
          SimulationDuration = numofcycle * 2d0 * pi / omegad(0d0,1)  !recalc for integer number of cycles
          
          nstepFZ = ceiling(real(numincycle*numofcycle)/plotpnts)
          
       else if (operating_mode == APPROACH_TRIGGERED) then
          SimulationDuration = 2*(Z0+TriggerDeflLim)/GetApproachSpeedSetpoint() !this is guess and could be off by a lot for DLVO
          numofcycle=ceiling(SimulationDuration*omegad(0d0, 1)/2.0/pi) !round up to an integer
          SimulationDuration = numofcycle * 2d0 * pi / omegad(0d0,1)  !recalc for integer number of cycles
          nstepFZ = ceiling(real(numincycle*numofcycle)/plotpnts)

       else if (operating_mode == FREQSWEEP) then
          SimulationDuration = Sweep_time
          if (( operating_mode == FREQSWEEP) .and. (sweepchoice == CONTSWEEP)) then
             numofcycle = ceiling( omegad_start * Sweep_time + sweep_rate * (Sweep_time**2) / 2 /2d0/pi ) !unused for stairsteps
          end if

       else if ( operating_mode == SCAN) then
          SimulationDuration = SubsLen / Vscan
          numofcycle=ceiling( (SubsLen / Vscan) *omegad(0d0, 1)/2.0/pi) !round up to an integer number
          SimulationDuration = numofcycle * 2d0 * pi / omegad(0d0, 1)  !recalc for integer number of cycles 

       else if ( operating_mode == FIXED) then
          numofcycle = hist_cycles
          SimulationDuration = numofcycle * 2d0 * pi / omegad(0d0, 1)  !recalc for integer number of cycles  
       else
          call assert( .false., 'unknown mode')
       end if

       call assert( numofcycle > 0, 'numofcycle >0')

       call UpdateDelt(0d0)

       if (isOpModeApp(operating_mode)) then
          nstep = numofcycle*numincycle

          call assert(nstep > 0,'nstep > 0')
       elseif (operating_mode == FIXED) then
          nstep = numofcycle*numincycle + 2!there are some off-by-one or fence-post errors with fixed point.  this is a kludge
       end if

       !for modes with controllers, how to long to wait before giving up
       if (transient_timeout_in > 0) then
          transient_timeout = transient_timeout_in
       elseif (modulation_type == PEAK_FORCE) then
          transient_timeout = 9d3
       else
          transient_timeout = max( 9000d0, max(500d0 * LockInTC,20d0 * 2d0*Quality(1)/(NonDimenFreq(omegai(1)))))
       end if

       !ok.  now that we know what we are doing, how often to output data?
       if (Want_Fourier) then
          !we are not windowing the fourier integrals.  therefore, we can only output integer multiples of cycle.
          output_point_rate = ceiling( real(numofcycle)/plotpnts) * numincycle
       else
          !for modes like say contact mode, fz curves, (or even force modulation) we won't be using fourier integrals
          !and we won't want to lock ourselves into one output per drive cycle, which won't really even by defined
          !for static modes.  
          output_point_rate = ceiling( real(numofcycle * numincycle)/plotpnts) 
       end if

       if (modulation_type == PEAK_FORCE) then
          !for peak force scanning, outputs more frequently than once per cycle does not make sense
          output_point_rate = max( output_point_rate, numincycle)
       end if
       
       call assert(output_point_rate >= 0, 'ouput_point_rate >= 0')


     end subroutine CalcTimeCycleInfo


     ! t = time, i = index of drive freq (1 = first drive freq, 2 = second drive freq)
     real*8 function theta( t, i)
       real*8, intent(in) :: t
       integer, intent(in) :: i

       if (( modulation_type == AMPLITUDE) .or. (modulation_type == PEAK_FORCE)) then
          if (operating_mode .ne. FREQSWEEP) then
             !drive freq is constant, so integration yields linear (assume constant = 0)
             theta = cur_omegad(i) * t
          else 
             if ( (sweepchoice == CONTSWEEP) .and. ( isTransientOver) ) then
                !drive freq is linear function of time, so integration yields quadratic
                theta = last_theta + omegad_start * timeSinceTrans(t) + sweep_rate * (( timeSinceTrans(t) ) ** 2) / 2
             else 
                ! DDASKR use a backwards differentiation scheme.  that means that successive calls to theta may not
                ! have monotonically increasing time.  In particular, we might have successive calls that are on either
                ! side of the time where we switched from one frequency to another one.  So we need to keep track of
                ! enough information to calculate theta on either side of the most recent freq step
                if (t < last_time) then
                   theta = last_theta - ( last_time - t) * prev_omegad
                else
                   theta = last_theta + ( t - last_time) * cur_omegad(1)
                end if
             end if
          end if
       else
          !modulation == FREQUENCY.  this is a simple rectangular integration.  
          !fixme: be more accurate and use trapezoids or something
          if (t < last_time) then
             theta = last_theta - ( last_time - t) * prev_omegad
          else
             theta = last_theta + ( t - last_time) * cur_omegad(1)
          end if
       end if
     end function theta

     !this is the drive frequency.  it's just a constant for AM approach curves, 
     !but not constant for freq sweeps or FM
     real*8 function omegad(t, i)
       real*8, intent(in) :: t
       integer, intent(in) :: i

       if ((operating_mode == FREQSWEEP) .and. (sweepchoice == CONTSWEEP)) then
          !drive freq is linear function of time
          if ( isTransientOver ) then
             omegad = omegad_start + (omegad_stop - omegad_start) * timeSinceTrans(t) / sweep_time
          else
             omegad = cur_omegad(1) ! this is not necessarily the same as omegad_start anymore.
          end if
       else 
          omegad = cur_omegad(i)
       end if
     end function omegad

     subroutine set_omegad( omegad_in)
       real*8 omegad_in(maxExc)
       cur_omegad = omegad_in
     end subroutine

     subroutine set_omegad_startstop( o_start, o_stop)
       real*8 o_start, o_stop
       omegad_start = o_start
       cur_omegad(1) = omegad_start
       omegad_stop = o_stop
     end subroutine set_omegad_startstop

     subroutine setup_next_stairstep_point(t)
       real*8, intent(in) :: t
       real *8 :: sweep_fraction, tmp_theta
       integer point

       sweep_fraction = timeSinceTrans(t) / sweep_time !what percentage of the sweep is complete
       point =  floor(sweep_fraction * plotpnts)

       current_stairstep_point = point
       tmp_theta = theta(t,1) !

       prev_omegad = cur_omegad(1)
       cur_omegad(1) = omegad_start + sweep_fraction * (omegad_stop - omegad_start)

       last_theta = tmp_theta
       last_time = t     

     end subroutine setup_next_stairstep_point

     subroutine setup_next_fm_point(t, new_omegad)
       real*8, intent(in) :: t, new_omegad
       real*8 :: tmp_theta

       tmp_theta = theta(t,1)
       prev_omegad = cur_omegad(1)
       cur_omegad(1) = new_omegad

       last_theta = tmp_theta
       last_time = t
     end subroutine setup_next_fm_point

     integer function transientPercent(t)
       use Controller, only: ApproachPercent
       real*8, intent(in) :: t
       if (isOpModeApp(operating_mode) .or. (operating_mode == FIXED)) then
          transientPercent = idint(100d0*t/ttrans)
       else
          transientPercent =  ApproachPercent()  !not 100% correct for freq sweep b/c there is more ttrans after.  also may not be monotonic
       end if
     end function transientPercent

     logical function freqSweepKeepGoing(t)
       real*8, intent(in) :: t
       if ( omegad_stop > omegad_start) then
          freqSweepkeepGoing = ( omegad(t,1) < omegad_stop )
       else
          freqSweepkeepGoing = ( omegad(t,1) > omegad_stop )
       end if
     end function freqSweepKeepGoing
   end module timeAndCycle
   

   module ForceVolume
     implicit none
     save
     integer :: ForceVolState
     logical :: If_Output_ForceVolZ, t_scan_set, Force_Vol_FirstCall  
     real*8  :: F_ForceVol, ForceVolZStart, ForceVolZStop, ForceVolTStop, ForceVolTStart, ForceVolSettleTime, t_UptoScan, Zbase_relative, totalTimeApprRev, ZForceVol_Output      
     integer, parameter :: FORCEVOL_SCAN = 0, FORCEVOL_APPR = 1, FORCEVOL_SETTLE = 2, FORCEVOL_REVERSE = 3
                        
   contains
     
     subroutine init_ForceVolData(Z0)
       !need this subroutine to set the Zstart. Don't want to include it in the peakForceReset it needn't jave to be called multiple times.
       real*8, intent(in) :: Z0
       forceVolZstart = Z0
       t_UptoScan = 0
       t_scan_set = .false.
       Force_Vol_FirstCall = .true.
       totalTimeApprRev = 0
       ZForceVol_Output = 0
       !the call to the subroutine is done on only one instance from here. The rest of the reset calls are made in updatePeakForce subroutine
       call ForceVolReset()
     end subroutine init_ForceVolData

     subroutine ForceVolReset()
       ForceVolState = ForceVol_SCAN
       ForceVolTStop = 0
       ForceVolTStart = 0
     end subroutine ForceVolReset
     

     !to find the approached distance the time taken for the approach has to be explicitly 
     !differentiated from the time taken for the scanning     
     real*8 function timeForceVolApprRev(t)
       real*8, intent(in) :: t
       ! this function calculates the time taken for the approach alone and does not include the time taken for the scan.
       !Used for finding the approach distance in peakforce.
       timeForceVolApprRev = t - t_UptoScan     
     end function timeForceVolApprRev
     
     real*8 function ForceVolApproachedDistance(t, AprchS_cur)
       use Nondimensionalization
       real*8, intent(in) :: t, AprchS_cur
       if (ForceVolState == FORCEVOL_APPR) then
          ForceVolApproachedDistance = -AprchS_cur*timeForceVolApprRev(t)
       elseif (ForceVolState == FORCEVOL_SETTLE) then
          ForceVolApproachedDistance = ForceVolZStop
       elseif (ForceVolState == FORCEVOL_SCAN) then 
          ForceVolApproachedDistance = 0   
       else 
          ! FORCEVOL_REV
          ForceVolApproachedDistance = ForceVolZStop - AprchS_cur*( t - ForceVolTStart)
       end if         

     end function ForceVolApproachedDistance

     subroutine updateForceVol(Force, t)
      use Approaching
      use TimeAndCycle, only : isTransientOver
       real*8, intent(in) :: t, Force
       if ((isTransientOver) .and. (.not. t_scan_set)) then
          t_UptoScan = t
          t_scan_set = .true.
       end if
       if  (ForceVolState == FORCEVOL_APPR) then
          if (Force > F_ForceVol) then
             ForceVolZStop = ForceVolApproachedDistance(t, GetApproachSpeedCurrent(t) )
             ZForceVol_Output = ForceVolZStart + ForceVolZStop
             ForceVolTStop = t
             ForceVolTStart = t + ForceVolSettleTime       
             call ApproachStop()
             ForceVolState = FORCEVOL_SETTLE      
             If_Output_ForceVolZ = .true.
             if (Force_Vol_FirstCall) then
                Zbase_relative = ZForceVol_Output
                Force_Vol_FirstCall = .false.
             end if
             !ZPeak_output is actually the error measurement. ie. the difference between the actual topography and measured topography.
             !Zbase_relative is the first measured point taken as a reference.
             ZForceVol_Output = ZForceVol_Output - Zbase_relative
          end if
       elseif (( ForceVolState == FORCEVOL_SETTLE) .and.  (t > ForceVolTStart)) then
          ForceVolState = FORCEVOL_REVERSE
          call ApproachReverse()
          If_Output_ForceVolZ = .false.
       elseif ((ForceVolState == FORCEVOL_REVERSE) .and. ( ForceVolApproachedDistance(t, GetApproachSpeedCurrent(t) ) > 0 )) then
          call ForceVolReset()             
          ForceVolState = FORCEVOL_SCAN
          call ApproachStop()
          totalTimeApprRev = totalTimeApprRev + timeForceVolApprRev(t)
       elseif ( ForceVolState == FORCEVOL_SCAN) then
          t_uptoScan = t
       end if      
     end subroutine updateForceVol
     
   end module ForceVolume


   !perhaps poorly named.  untriggered FZ curves are handled here too, but with the "trigger" set to a z distance
   !that way we can accomodate the rounded triangle shape in both cases.
   module TriggeredFzMode
     use params
     implicit none
     integer :: triggerState, fzshape
     real*8  :: triggerDeflLim, triggerZStop, triggerTStop, triggerTStart, PiezoReverseTime_pct
     integer, parameter :: TRIGGER_FWD = 1, TRIGGER_SETTLE = 2, TRIGGER_REVERSE = 3, TRIGGER_ROUNDING=4
     integer, parameter, private :: TRIANGLE=1, ROUNDED_TRIANGLE=2
     real*8,private :: trigger_vel, PiezoReverseTime
     
   contains
               
     subroutine TriggeredFzModeInit()
       triggerState = TRIGGER_FWD
       triggerTStop = 0
       triggerTStart = 0
     end subroutine TriggeredFzModeInit

     real*8 function TriggeredFzModeApproachedDistance(t, AprchS_cur)
       use timeAndCycle, only: timeSinceTrans
       real*8, intent(in) :: t, AprchS_cur
       
       if (triggerState == TRIGGER_FWD) then
          TriggeredFzModeApproachedDistance = -AprchS_cur*timeSinceTrans(t)
       elseif (triggerState == TRIGGER_SETTLE ) then
          TriggeredFzModeApproachedDistance = triggerZStop
       elseif (triggerState == TRIGGER_REVERSE ) then
          TriggeredFzModeApproachedDistance = triggerZStop - AprchS_cur*( t - triggerTStart)
       elseif (triggerState == TRIGGER_ROUNDING) then
          TriggeredFzModeApproachedDistance = triggerZStop +  (t - triggerTstop)**2 * (  trigger_vel / PiezoReverseTime) - trigger_vel * (t - triggerTstop)
       end if
     end function TriggeredFzModeApproachedDistance

     subroutine updateTriggeredFzMode(u, t)
       use Approaching
       use TimeAndCycle, only: operating_mode,  SimulationDuration, timeSinceTrans
       use NonDimensionalization
       
       real*8, intent(in) :: u, t
          
       if  (triggerState == TRIGGER_FWD) then
          if ((operating_mode == APPR_RET) .and. (timeSinceTrans(t) > (SimulationDuration / 2d0)) .or. ((operating_mode == APPROACH_TRIGGERED) .and. (u > TriggerDeflLim))) then
             triggerZStop = triggeredFzModeApproachedDistance(t, GetApproachSpeedCurrent(t) )
             triggerTStop = t
             
             if (( fzshape == TRIANGLE) .or. (PiezoReverseTime_pct == 0)) then
                !triangle versions             
                triggerState = TRIGGER_REVERSE             
                triggerTStart = t
                call ApproachReverse()                          
             elseif (fzshape == ROUNDED_TRIANGLE) then
                !rounded version                          
                triggerState = TRIGGER_ROUNDING
                !we delay calculating until here, as for triggered curves we don't until we get here how long the
                !approach will take
                PiezoReverseTime = 2 * t * (PiezoReverseTime_pct / 100d0)
                triggerTStart = t + PiezoReverseTime
                trigger_vel = GetApproachSpeedCurrent(t)
                call SetApproachParabola( -trigger_vel/PiezoReverseTime, trigger_vel, t)
             end if
          end if
       elseif (triggerState == TRIGGER_ROUNDING) then
          if ( t > triggerTstart) then
             triggerState = TRIGGER_REVERSE
             call ApproachReverse()
          end if
       elseif (( triggerState == TRIGGER_SETTLE) .and.  (t > triggerTStart)) then
          !this state is currently unused
          triggerState = TRIGGER_REVERSE
          call ApproachReverse()          
       end if
       
     end subroutine updateTriggeredFzMode

   end module TriggeredFzMode


  module AppControl
    use params   
    use TimeAndCycle, only : isTransientOver
    implicit none
    real*8 Asp_f, freqshift_f

    contains


     !returns true if we should output right now, false otherwise
     !stairstep code decides point boundaries soley based on time. 
     !other amplitude mode codes put the boundary on an integer number of drive cycles     
     logical function OutputThisStep( t,  output_point_rate, IOUT, IOUT_prevOutput,theta, prev_theta, nstep, Z_Controller_On ) 
       use timeAndCycle, only: operating_mode, sweepchoice, timeSinceTrans, pointsSinceTrans, current_stairstep_point, modulation_type, numincycle, plotpnts, sweep_time
       use forcevolume, only : If_Output_ForceVolZ
       integer*8, intent(in) :: output_point_rate, IOUT, IOUT_prevOutput, nstep
       real*8, intent(in) :: t, prev_theta, theta
       logical, intent(in) :: Z_Controller_On
       real*8 sweep_fraction
       integer :: point

       if ( (operating_mode == FREQSWEEP) .and. (sweepchoice == STAIRSTEP)) then
          sweep_fraction = timeSinceTrans(t) / sweep_time !what percentage of the sweep is complete
          point =  floor(sweep_fraction * plotpnts)
          OutputThisStep = (point > current_stairstep_point)   
       elseif (operating_mode == FIXED) then
          OutputThisStep = (pointsSinceTrans(IOUT) ==  nstep-1)
 !originally wanted to output on drive cycle boundaries so that fourier integrals would still work.
 !but that makes anti-aliasing filtering harder b/c output rate is not constant.  so go back to other way.
 !      elseif (modulation_type == FREQUENCY) then
 !         !output on drive cycle boundaries
 !         OutputThisStep = ((cos(theta) > 0) .and. (cos(prev_theta) < 0) .and. ( IOUT - IOUT_prevOutput > numincycle*datacycle))
 !      elseif (modulation_type == PEAK_FORCE) then
 !            OutputThisStep = (If_Output_PeakForceZ) 
       elseif ((modulation_type == FORCE_VOL)) then
          OutputThisStep = ((mod( pointsSinceTrans(IOUT), 70*output_point_rate).EQ.0).and. ( pointsSinceTrans(IOUT).ne.0))
       elseif ((modulation_type == PEAK_FORCE)) then
          OutputThisStep = ( Z_Controller_On )  .and. (isTransientOver)
       else
          OutputThisStep = ((mod( pointsSinceTrans(IOUT), output_point_rate).EQ.0).and. ( pointsSinceTrans(IOUT).ne.0))
       end IF
     end function OutputThisStep


     logical function keepGoing( I, N,t, Xpos, SubsLen, AprchS_cur, Zrange, Amp, freqshift)
       use timeAndCycle, only: freqSweepKeepGoing, operating_mode
       use TriggeredFzMode, only: TriggeredFzModeApproachedDistance
       integer*8, intent(in) :: i, n
       integer, intent(in) :: Zrange
       real*8, intent(in) :: t, XPos,SubsLen, AprchS_cur, Amp, freqshift

       if (operating_mode == FREQSWEEP) then
          keepGoing = freqSweepKeepGoing(t)
       else if (operating_mode == SCAN) then
          keepGoing = XPos < SubsLen
       elseif  ((operating_mode == FIXED) .or. (operating_mode == APPR_RET) .or. (operating_mode == APPROACH_SINE) ) then
          keepGoing = (I < N)
       else if ((operating_mode == APPROACH ) .or.  (operating_mode == APPROACH_STEP )) then
          if ( Zrange == ZRANGE_ASP) then
             keepGoing = Amp > Asp_f
          elseif (Zrange == ZRANGE_FREQSHIFT) then
             if (freqshift_f > 0) then                           
                keepGoing = freqshift < freqshift_f
             else
                keepGoing = freqshift > freqshift_f
             end if
          else
             keepGoing = (I < N)
          end if
       else if (operating_mode == APPROACH_TRIGGERED) then
          keepGoing = (TriggeredFzModeApproachedDistance(t, AprchS_cur) <= 0)
       else
          call assert(.false., 'unhandled case1')
          keepGoing = .false.
       end if
     end function keepGoing

  end module AppControl


!this module is a collection of a bunch of misc stuff that doesn't fit anywhere else.  
!should probably find better homes for most of this stuff
   module data1
     use params
     implicit none
     save

     real*8 Z0, Zf, F(maxExc, maxModes), F1_init(maxModes), &
             Ainitial_dim(maxModes), mstar_div_m
     real*8 damping(maxModes),Quality(maxModes), Abase(maxExc), Abase_init
     real*8 phid(maxExc, maxModes), B(maxModes), beta(maxModes), alpha(maxModes), drag(maxModes), gamma_drag, drag_pf(maxModes)
     real*8 Asample, omegas ! subsurface parameters
     real*8 mu(maxModes), Chi(maxModes)
     real*8 Amp,mtip
     complex*8 :: Afluid
     logical wantSampleExc
     integer NEQ, NEQ_cant, EQ_attard , fexcite
     integer numModes ! number of modes in the current problem (versus the maximum possible in module params)
     integer ::  NRT  !number of boundaries. 1 = aDMT. boundary
     integer :: numRes1 !performance counter
     integer*8 :: IOUT,  output_point_rate

!     real*8, allocatable, private :: RWORK_save(:),  Y_save(:), YPRIME_save(:)
!     integer,allocatable, private :: IWORK_save(:)
!     real*8, private :: t_save

   contains
     
     subroutine setup_DDASKR_parameters(NEQ, DELT, LRW, LIW, INFO, RWORK, ATOL, RTOL, IWORK)
       integer, intent(in ) :: NEQ
       real*8, intent(in)   :: DELT
       integer, intent(out) :: INFO(20), LIW, LRW
       real*8, allocatable,  intent(out) :: ATOL(:), RTOL(:), RWORK(:)
       integer, allocatable, intent(out) :: IWORK(:)



       allocate(ATOL(NEQ))
       allocate(RTOL(NEQ))

       !john just used LRW=LIW=100*NEQ.  here is the proper computation based on the ddaskr documentation
       LRW =  60 + 9*NEQ + 3*NRT + NEQ**2
       LIW = 40 + NEQ

       allocate(RWORK(LRW))
       allocate(IWORK(LIW))       

        INFO = 0 !all array elem

        INFO(1) = 0      !     INDICATE START OF A NEW PROBLEM

        INFO(2) = 1         ! this indicates THAT RTOL AND ATOL ARE ARRAYS.  
                          ! EACH ENTRY OF RTOL AND ATOL MUST THEN BE DEFINED.

        INFO(5) = 0      !     TELL DDASKR THE JACOBIAN TYPE. 0 = no jacobian function provided.  use numerical derivatives

      !this is important! Do not change this unless you know what you are doing!
      !prevent ddaskr from taking too big of steps and thus missing 
      !when we switch things. outside of RES1.  If we don't do this, then anything that 
      !gets changed outside of RES1 might get completely ignored!  in some situations, DDASKR
      !may get fooled into thinking it can take quite long time steps even though we are
      !changing things outside of RES1 (turning the approach speed on for f-z curves after 
      !transient is over for example).  Note that contrary to what the VEDA manual or the original RSI article says, the
      !detection of a boundary in RT1 does not cause DDASKR to back up and re-integrate
      !any points. If it has integrated past the boundary then it just keeps going.  Therefore
      !this is necessary to limit how far DDASKR can step.
      !
      !note, this is definitely the safest option, but not necessarily the most efficient.
      !for example, in an approach curves example, with pure hertz contact, that starts in intermittent contact
      !with the sample, (and no feedback controller, no hysteretic forces) then this restriction is entirely unneccesary.

      INFO(7) = 1 
      RWORK(2) = DELT 
      
      ! SET INFO(11) = 1 IF DDASKR IS TO COMPUTE THE INITIAL YPRIME, AND
      ! GENERATE AN INITIAL GUESS FOR YPRIME.  OTHERWISE, SET INFO(11) = 0
      ! AND SUPPLY THE CORRECT INITIAL VALUE FOR YPRIME.
      INFO(11) = 0

      RTOL = 1d-9
      ATOL = 1d-9      

!see below
!      allocate(RWORK_save(LRW))
!      allocate(IWORK_save(LIW))
!      allocate(Y_save(NEQ))
!      allocate(YPRIME_save(NEQ))      
    end subroutine setup_DDASKR_parameters

!this was a hack.  didn't really help much.
   !  !this stores the complete state of ddaskr.  we can then later restore this state. essentially this allows us to "undo" 
   !  !one or more solution steps and re-solve them. was a hack for ting viscoelasticity to make it actually work the way john
   !  ! and shuiqing wanted it to work to begin with.
   !  subroutine save_DDASKR_state(t, Y, YPRIME, RWORK, IWORK)
   !   real*8,  intent(in) :: t, RWORK(:),  Y(:), YPRIME(:)
   !   integer, intent(in) :: IWORK(:)

   !   t_save = t
   !   RWORK_save = RWORK
   !   IWORK_save = IWORK
   !   Y_save = Y
   !   YPRIME_save = YPRIME
   ! end subroutine save_DDASKR_state

   !  subroutine restore_DDASKR_state(t, Y, YPRIME, RWORK, IWORK)
   !   real*8,  intent(out) :: t, RWORK(:),  Y(:), YPRIME(:)
   !   integer, intent(out) :: IWORK(:)

   !   t = t_save
   !   RWORK = RWORK_save
   !   IWORK = IWORK_save
   !   Y = Y_save
   !   YPRIME = YPRIME_save
   ! end subroutine restore_DDASKR_state


   !the Z distance is a quantity that is computed in a function every time that it is needed,
   !rather than being stored in a variable.  the advantage of a function is that it forces everything to be
   !in once place. you can't have multiple different routines manipulating the variable.  on the other hand
   !you need all of the necessary local variables available to make the computation all the time.
     real*8 function computeZdist_dim(t)
       use TriggeredFzMode, only: TriggeredFzModeApproachedDistance
       use Approaching, only: GetApproachSpeedCurrent, CurveRate
       use timeAndCycle, only : operating_mode, timeSinceTrans, SimulationDuration,   pointsSinceTrans, plotpnts, theta
       use NonDimensionalization
       real*8, intent(in) :: t


       !fixme, more or less duplicate of below routine.
       if (operating_mode == APPROACH) then
          computeZdist_dim = DimenLength((Z0 - GetApproachSpeedCurrent(t)*timeSinceTrans(t) )*1d9)
       elseif  (operating_mode == APPROACH_STEP) then
          computeZdist_dim = DimenLength(Z0 + ( pointsSinceTrans(IOUT) / output_point_rate ) * ( Zf - Z0 ) / plotpnts) * 1d9
       elseif (operating_mode == APPR_RET) then
          computeZdist_dim = DimenLength((Z0 + TriggeredFzModeApproachedDistance( t, GetApproachSpeedCurrent(t) ) ) *  1d9)
       else if (operating_mode == APPROACH_TRIGGERED) then
          computeZdist_dim = DimenLength((Z0 + TriggeredFzModeApproachedDistance( t, GetApproachSpeedCurrent(t) ) ) *  1d9)
       elseif (operating_mode == APPROACH_SINE ) then
          computeZdist_dim = DimenLength(  (Zf-Z0) * sin( pi * timeSinceTrans(t) * CurveRate) ) * 1d9
       elseif (fexcite == ACOUSTIC_PEAK) then
          computeZdist_dim = DimenLength( Z0 + Abase(1)*cos(theta(t,1) ))
       else
          call assert(.false., 'z dist called for unhandled case');
          computeZdist_dim =0 
       end if

     end function computeZdist_dim

     real*8 function computeTipSampleGap(t, y, neq, cos_in)
       use TriggeredFzMode, only: TriggeredFzmodeapproacheddistance
       use Approaching
       use forcevolume, only: ForceVolApproachedDistance
       use timeAndCycle, only: operating_mode, transient_suppress_Fts, timeSinceTrans, SimulationDuration, Modulation_type, theta, exc_choice, isTransientOver,  pointsSinceTrans, plotpnts
       
       integer, intent(in) :: neq
       real*8, intent(in) :: t,  y(neq)
       real*8, optional :: cos_in       

       real*8 :: ApproachedDistance, cos1

       if ( isTransientOver ) then
          !fixme, more or less duplicate of above routine.
          if (operating_mode == APPROACH) then
             ApproachedDistance = -GetApproachSpeedCurrent(t)*timeSinceTrans(t)
          elseif (operating_mode == APPROACH_STEP) then
             ApproachedDistance = ( pointsSinceTrans(IOUT) / output_point_rate ) * ( Zf - Z0 ) / plotpnts
          elseif (operating_mode == APPR_RET) then
             !fixme, is this correct for Zrange == ZRANGE_ASP????             
             ApproachedDistance = TriggeredFzModeApproachedDistance( t, GetApproachSpeedCurrent(t))
          elseif (operating_mode == APPROACH_TRIGGERED ) then            
             ApproachedDistance = TriggeredFzModeApproachedDistance( t, GetApproachSpeedCurrent(t))
          elseif (operating_mode == APPROACH_SINE ) then
             ApproachedDistance = (Zf-Z0) * sin( pi* timeSinceTrans(t) * CurveRate)
!             write(*,*) CurveRate, timeSinceTrans(t), ApproachedDistance
          elseif (modulation_type == FORCE_VOL) then
             ApproachedDistance = ForceVolApproachedDistance(t, GetApproachSpeedCurrent(t) )
          else
             ApproachedDistance = 0
          end if
       else
          if (( modulation_type == FREQUENCY) .and. (transient_suppress_Fts)) then
             !hack for FM.  want to calculate the phase setpoint in the absence of tip-sample forces
             ApproachedDistance = 1000
          else
             ApproachedDistance = 0
          end if
       end if

       
       ComputeTipSampleGap=Z0 + sum(Y(1:neq:DOF)) + ApproachedDistance ! sum y(1)+y(3)+etc
       !sines and cosines are expensive so only computed if needed.  Note the routine is so sensitive
       !to roundoff error, that skipping the addition of zero causes the final answer to change. 
       if (isAcoustic(fexcite) .or. fexcite == ACOUSTIC_PEAK) then
          if ( present(cos_in) ) then
             cos1 = cos_in
          else           
             cos1 = cos(theta(t,1))
          end if

          if (exc_choice == BIMODAL) then
             ComputeTipSampleGap = ComputeTipSampleGap + Abase(1)*cos1+Abase(2)*cos(theta(t,2))
          else             
             ComputeTipSampleGap = ComputeTipSampleGap + Abase(1)*cos1
          end if
       end if

       if (wantSampleExc) then
          ComputeTipSampleGap = ComputeTipSampleGap + Asample*cos(omegas*t)
       end if
     end function computeTipSampleGap

    !  !this was an idea I had to make something better or quicker, can't really remember what now
    !  !and I don't think it worked.
    !  real*8 function computeTipSampleGapDeriv_YPRIMEverion(t, y, yprime, neq, sin_in)
    !   use Approaching, only: GetApproachSpeedCurrent
    !   use timeAndCycle
    !    integer, intent(in) :: neq
    !    real*8, intent(in) :: y(neq), yprime(neq), t

    !    real *8 dpr, sin1

    !    real*8, optional :: sin_in


    ! !sines and cosines are very expensive. only compute if needed. 
    !    if (isAcoustic(fexcite) .or. (modulation_type == PEAK_FORCE)) then
    !       if ( present(sin_in) ) then 
    !          sin1 = sin_in
    !       else
    !          sin1 = sin(theta(t,1))
    !       end if

    !       if (exc_choice == BIMODAL) then
    !          dpr = -Abase(1)*omegad(t, 1)*sin1 -Abase(2)*omegad(t,2)*sin(theta(t,2)) + sum(Y(2:neq:2))
    !       else
    !          dpr = -Abase(1)*omegad(t, 1)*sin1  + sum(Y(2:neq:2))
    !       end if
    !    else
    !       dpr =  sum(YPRIME(1:neq:DOF))
    !    end if

    !    if (wantSampleExc) then
    !       dpr = dpr - Asample*omegas*sin(omegas*t)
    !    end if

    !    computeTipSampleGapDeriv_YPRIMEverion = dpr - GetApproachSpeedCurrent(t)        
       
    !  end function computeTipSampleGapDeriv_YPRIMEverion

     
     !derivative of above
     real*8 function computeTipSampleGapDeriv(t, y, neq, sin_in)
      use Approaching, only: GetApproachSpeedCurrent, CurveRate
      use timeAndCycle
       integer, intent(in) :: neq
       real*8, intent(in) :: y(neq), t

       real *8 dpr, sin1

       real*8, optional :: sin_in


    !sines and cosines are very expensive. only compute if needed. 
       if (isAcoustic(fexcite) .or. (modulation_type == PEAK_FORCE)) then
          if ( present(sin_in) ) then 
             sin1 = sin_in
          else
             sin1 = sin(theta(t,1))
          end if

          if (exc_choice == BIMODAL) then
             dpr = -Abase(1)*omegad(t, 1)*sin1 -Abase(2)*omegad(t,2)*sin(theta(t,2)) + sum(Y(2:neq:2))
          else
             dpr = -Abase(1)*omegad(t, 1)*sin1  + sum(Y(2:neq:2))
          end if
       else
          dpr =  sum(Y(2:neq:DOF))
       end if

       if (wantSampleExc) then
          dpr = dpr - Asample*omegas*sin(omegas*t)
       end if

       if (operating_mode == APPROACH_SINE) then
          !this is a hack, but getting this into GetApproachSpeedCurrent will be too hard
          dpr = dpr - (Zf-Z0) * cos( pi* timeSinceTrans(t) * CurveRate)  * pi* timeSinceTrans(t) * CurveRate
       end if
       
       computeTipSampleGapDeriv = dpr - GetApproachSpeedCurrent(t)        
       
     end function computeTipSampleGapDeriv

     subroutine CalcCantileverData(AutoCalcAlpha, numModes, mtip, alpha, beta, mu, fexcite, want_Custom_B, B, B_input, AutoCalcOmega, omegai, CalcInputK, Keq, invModalMass)
       
       logical, intent(in) :: AutoCalcAlpha, AutoCalcOmega, CalcInputK, want_Custom_B 
       real*8, intent(in) :: mtip
       real*8, dimension(numModes), intent(inout) :: alpha, beta, mu, B, keq, omegai, invModalMass,  B_input
       integer :: i
       integer, intent(in) :: fexcite, nummodes
       
       if ( AutoCalcAlpha ) then
          do i = 1, numModes
             alpha(i) = CalcAlpha(mtip, i)
          end do
       end if
       
       if (AutoCalcOmega) then 
          do i = 2,numModes 
             omegai(i) = omegai(1)* alpha(i)**2  / alpha(1)**2 
          end do
       end if


       do i = 1,numModes
          beta(i) = ModeShapeIntegral( alpha(i))
          mu(i) = ModeShapeSquaredIntegral( alpha(i))    
       end do
       
       if (isMagnetic(fexcite)) then
          if (want_Custom_B) then
             B = B_input
          elseif (fexcite == MAGNETIC_TORQUE) then
             B = 1d0
          elseif (fexcite == MAGNETIC_LORENTZ) then
             ! lorentz force excitation
             B = beta
          end if
       end if       
     
       if (CalcInputK) then
	  do i = 2,numModes
            Keq(i) = Keq(1)*Calcki(alpha(i),alpha(1))
         end do
      end if
      
      do i = 1,numModes
         invModalMass(i) = (omegai(i)**2) / Keq(i)
      end do
           
     end subroutine CalcCantileverData

     ! given a dispersion, integrate the square of the mode shape over the length of the cantilever.
     ! called this mu in the manual

     pure real*8 function ModeShapeSquaredIntegral(aa)
       !	calculate B

       real*8, intent(in) :: aa
       real*8 Sum2, xx, delta, Phi1, phi, uu
       integer klook, last

       Phi1 = sin(aa)-sinh(aa)
       Phi1=Phi1-((sin(aa)+sinh(aa))/(cos(aa)+cosh(aa))) *(cos(aa)-cosh(aa))

       Sum2 = 0
       last = 1000
       delta = 1/1d3
       do klook = 1,last     
          xx = (klook-0.5)/1d3
          uu = aa*xx
          phi = sin(uu)-sinh(uu)-((sin(aa)+sinh(aa))/(cos(aa) +cosh(aa)))*(cos(uu)-cosh(uu))
          phi = phi/Phi1
          Sum2=Sum2+((phi**2)*delta)	     
       end do

       ModeShapeSquaredIntegral = Sum2

       Return
     End function ModeShapeSquaredIntegral


     ! ******************************************    
     ! given a dispersion, find the integral of the mode shape over the length of the cantilever
     ! (used to be called BMag, but it's actually more general than that... I called this beta in the manual)
     pure real*8 function ModeShapeIntegral(aa)

       real*8, intent(in) :: aa
       real*8 Sum1, xx, delta, Phi1, phi, uu
       integer kloik, last


       Phi1=sin(aa)-sinh(aa)-((sin(aa)+sinh(aa))/(cos(aa) +cosh(aa)))*(cos(aa)-cosh(aa))

       Sum1 = 0
       last = 1000
       delta = 1/1d3
       do kloik = 1,last

          xx = (kloik-0.5)/1d3
          uu = aa*xx
          phi = sin(uu)-sinh(uu)-((sin(aa)+sinh(aa))/(cos(aa) +cosh(aa)))*(cos(uu)-cosh(uu))
          phi = phi/Phi1
          Sum1=Sum1+(phi*delta)

       end do

       ModeShapeIntegral = Sum1

     End function ModeShapeIntegral


 ! ******************************************    
 ! ************ Calculate ki Function start
 ! computes strain energy integral numerically. 
 ! could get an analytical expression for this, but it's 8 lines in
 ! maple.  not worth re-writing.
 ! ******************************************
 pure real*8 function Calcki(aa,aa1)
 !	calculate k2, k3, ...

   real*8, intent(in) :: aa, aa1
   real*8 Sum1, Sum2, xx, delta, PhiL
   real*8 Phi1L, phipp, phi1pp, uu, uu1, sigma, sigma1
   integer kloik, last

   sigma1 = ((sin(aa1)+sinh(aa1))/(cos(aa1)+cosh(aa1)))
   Phi1L = sin(aa1)-sinh(aa1)-sigma1*(cos(aa1)-cosh(aa1))

   sigma = ((sin(aa)+sinh(aa))/(cos(aa)+cosh(aa)))
   PhiL = sin(aa)-sinh(aa)-sigma*(cos(aa)-cosh(aa))

   Sum1 = 0
   Sum2 = 0
   last = 1000
   delta = 1/1d3
   do kloik = 1,last

       xx = (kloik-0.5)/1d3
       uu = aa*xx
       uu1 = aa1*xx
       phi1pp = (aa1**2)*(sigma1*(cos(uu1)+cosh(uu1))-sin(uu1)-sinh(uu1))
       phi1pp = phi1pp/Phi1L
       phipp = (aa**2)*(sigma*(cos(uu)+cosh(uu))-sin(uu)-sinh(uu))
       phipp = phipp/PhiL
       Sum1=Sum1+((phipp**2)*delta)
       Sum2=Sum2+((phi1pp**2)*delta)

   end do

   Calcki = Sum1/Sum2

   Return
 End function Calcki

 ! frequency dispersion relationship for the eigenvalues.
 ! newton-rhapson Calculation of alpha for the given tip mass. use zero-tip mass value as a starting guess
 pure real*8 function CalcAlpha(mtip, modeNum)
   use params
   integer, intent(in) :: modeNum  
   real*8, intent(in) :: mtip
   real*8 a1i, a1j, fa1, dfa1

   !1.87 = 1st root of cos(alpha) * cosh(alpha) = -1 for zero tip mass cantilever
   !4.69 = 2nd root, etc. !after first three, (n-0.5)*pi is very very close.  See melcher, hu, & raman, apl 2007.
   real*8, dimension(3), parameter :: guess = (/ 1.875, 4.694, 7.855/) 

   if ( modeNum <= 3) then
      a1i = guess(modeNum)
   else
      a1i = ( modeNum - 0.5) * pi
   end if

   fa1 = (cos(a1i)*cosh(a1i))+1+a1i*mtip*(cos(a1i)*sinh(a1i) -cosh(a1i)*sin(a1i))

   DO WHILE (abs(fa1).gt.(10d-5))
       dfa1=cos(a1i)*sinh(a1i)-sin(a1i)*(cosh(a1i) +2*a1i*mtip*sinh(a1i))
       a1j = a1i - (fa1/dfa1)
       a1i = a1j
       fa1=cos(a1i)*cosh(a1i)+1+a1i*mtip*(cos(a1i)*sinh(a1i) -cosh(a1i)*sin(a1i))
   END DO

   CalcAlpha = a1i


   Return
 End function CalcAlpha


 subroutine CalculateChi( Chi, output_type )
   integer :: i
   integer, intent(in) :: output_type
   real*8, intent(inout) :: Chi(maxmodes)

   if (output_type == AUTO_CALC_CHI) then
      Chi = 1 !assigns whole array so we don't get a divide by zero error in ampl red avg tool
      do i = 2,numModes
         Chi(i) = CalcChi2(alpha(1),alpha(i))
      end do
   elseif (output_type == INTERFEROMETER) then
      Chi = 1
   end if
 end subroutine CalculateChi

 pure real*8 function CalcChi2(a1a,a2a)

   real*8, intent(in) :: a1a, a2a
   real*8 y1y, y2y, sigma1, sigma2

   sigma1 = (sin(a1a) + sinh(a1a))/(cos(a1a) + cosh(a1a))
   sigma2 = (sin(a2a) + sinh(a2a))/(cos(a2a) + cosh(a2a))

   y1y = a1a*(cos(a1a)-cosh(a1a)-sigma1*(-sin(a1a)-sinh(a1a)))
   y1y = y1y/(sin(a1a)-sinh(a1a)-sigma1*(cos(a1a)-cosh(a1a)))

   y2y = a2a*(cos(a2a)-cosh(a2a)-sigma2*(-sin(a2a)-sinh(a2a)))
   y2y = y2y/(sin(a2a)-sinh(a2a)-sigma2*(cos(a2a)-cosh(a2a)))


   CalcChi2 = y2y/y1y

   Return
 End function CalcChi2


pure real*8 function sgn(x)
   real*8, intent(in) :: x
   if (x < 0) then
      sgn = -1d0
   else
      sgn = 1d0
   end if
 end function sgn

   end module data1





 !the preisach module is not currently used in any published code.
 !it allows simulation of hysteretic models.  we may return to this in the future

 !to do: change the matlab code to solve and generate the preisach model into a 2-d array.
 !think that is mostly done, just have to test it
 module PreisachModule
   use params

   integer :: N_Preisach  !number of relays
   real*8, private, allocatable :: preisach_alpha(:)  ! tip-sample gap where relays turn on (in d decreasing direction)
   real*8, private, allocatable :: preisach_beta(:)  ! tip-sample gap where relays turn off (in d decreasing direction)
   real*8, private, allocatable :: weights_dim(:,:)
   real*8, private :: force_coeff, force, force_prev, last_d
   logical, private, allocatable :: states(:,:)
   logical, private, allocatable :: states_new(:,:)
   integer, private :: last_alpha, last_beta

   contains

     !format: alpha, beta, weights should be in SI units (meters, newtons)
     subroutine readCustomPreisachTable( filename)
       character(len=*), intent(in) :: filename
       integer :: i, j

       open( UNIT=10, FILE=filename)

       read(10, *) N_Preisach           

       allocate( preisach_alpha(N_Preisach))
       allocate( preisach_beta(N_Preisach))

       allocate( weights_dim(N_Preisach,N_Preisach ))
       allocate( states(N_Preisach,N_Preisach ))
       allocate( states_new(N_Preisach,N_Preisach ))

       do i = 1, N_Preisach
          read(10, *) preisach_alpha(i), preisach_beta(i)
       end do

       do i = 1, N_Preisach
          do j = 1,i
             read(10,*) weights_dim(i,j)
          end do
       end do

       close(UNIT=10)
     end subroutine readCustomPreisachTable


     subroutine Normalize_Preisach()
       use Nondimensionalization
       integer :: i
       do i = 1, N_Preisach
          preisach_alpha(i) = NonDimenLength(preisach_alpha(i))
          preisach_beta(i)  = NonDimenLength(preisach_beta(i))
       end do

       force_coeff =  NonDimenForce( 1d0)

     end subroutine Normalize_Preisach

     subroutine Initialize_Preisach(Z0)
       real*8, intent(in) :: Z0
       integer :: i, j

       !this is a hack for testing.  we'll need to read in a bunch of states from a file or something      
       ! forall (i = 1:N_Preisach)
       !    alpha(i) = 1 - 0.02 * i
       !    beta(i) = 1.01 - 0.02 * i
       ! end forall
       ! weights = -0.01


       states = .false.

       !better just do the full thing here once.  the update functions 
 !     try to be efficient, but apparently miss something on the first call
       do i=1,N_Preisach
          do j = 1,i
             if ( Z0 < preisach_alpha(i)) states(i,j) = .true.
             if ( Z0 > preisach_beta(j))  states(i,j) = .false.
          end do
       end do

       force = sum( weights_dim, states)

     end subroutine Initialize_Preisach

     real*8 function Fts_Preisach()
       Fts_Preisach = force  * force_coeff
     end function Fts_Preisach

     !as with capillary adhesion and jkr, the updates inside res1 are all considered temporary
     !and relative to the previous state.  the updates outside res1 are permanent

     !we have to consider the relays in a 2-D plane.  Then we can have N^2 relays yet only need to
     !search through N alpha values to find which ones to turn on.
     subroutine UpdatePreisach_insideRES1(d)
       real*8, intent(in) :: d
       integer :: i, j, a1

!option 1.  faster, more messy
       !the assumption here is that the motion has been monotonic since the last call.  i.e.
       !we've either turned relays on, or turned relays off, but no cycling.  since we typically
       !compute 1000 points per drive cycle and drive near 1st or 2nd eigenmode, we'd need to get 
       !10 or 15 eigenmodes involved before this would be a bad assumption

      force = force_prev

       if (d < last_d) then
          !because we already have a good guess, linear search is much quicker than
          !binary search
          a1 = lin_find( preisach_alpha, d, last_alpha)

          do  i = max(last_alpha-1, 1), min(a1+1, N_Preisach)
             if (d < preisach_alpha(i)) then
                do j = 1,i
                   if (states(i,j) .eqv. .false.) then
                      force = force + weights_dim(i, j)
                   !running sum is okay here.  we do the full calc outside res1 so roundoff errors can't accumulate too much.
                   end if
                end do
             end if
          end do
       else
          a1 = lin_find( preisach_beta, d, last_beta)
         
          do  j = max(a1-1, 1), min(last_beta+1, N_Preisach)
             if (d > preisach_beta(j)) then
                do i = j, N_Preisach
                   if (states(i,j) .eqv. .true.) then
                      force = force - weights_dim(i, j)
                   end if
                end do
             end if
          end do
       end if


!option 2. this is painly slow, but definitely correct
     ! states_new = states
     ! do i=1,N_Preisach
     !    do j = 1,i
     !       if ((d < last_d) .and.  (d < preisach_alpha(i))) states_new(i,j) = .true.
     !       if ((d > last_d) .and.  (d > preisach_beta(j)))  states_new(i,j) = .false.
     !    end do
     ! end do
     ! force = sum( weights_dim, states_new )

    end subroutine UpdatePreisach_insideRES1

    
!standard binary search.  fixme, move to basic_modules.
    integer function find( table, target)
      real*8, intent(in) :: target, table(:)

      integer :: n, n_guess, a, b
    
      n = size(table)

      if ( target > table(1) ) then
         find = 1
         return
      elseif (target < table(n)) then
         find = n
         return
      end if

      a = 1
      b = n

      n_guess = (a+b) / 2


      do
         if ( table(n_guess) > target) then
            a = n_guess
         else  
            b = n_guess
         end if

         n_guess = (a+b) / 2
         if ( a+1 >= b) exit
      end do

      find = n_guess

    end function find

!the value we are looking for is often very near to the value we wanted on the last iteration
!there, instead of doing a binary search, its actually quicker to do a linear scan starting
!from the previous value
 integer function lin_find( table, target, guess)
      real*8, intent(in) :: target, table(:)
      integer, intent(in) :: guess

      integer :: n, n_guess
   
      n = size(table)

      call assert( guess >= 0, 'guess > 0')
      call assert( guess <= n, 'guess < n')

      if ( target >= table(1) ) then
         lin_find = 1
         return
      elseif (target <= table(n)) then
         lin_find = n
         return
      end if


      n_guess = guess

      do
         if (target < table(n_guess)) then
            if (target >= table(n_guess+1)) then
               lin_find = n_guess
               exit
            else 
               n_guess = n_guess + 1
            end if
         else
            n_guess = n_guess - 1
         end if
      end do

    end function lin_find

    subroutine UpdatePreisach_outsideRES1(d, IOUT)      
      real*8, intent(in) :: d
      integer*8, intent(in) :: IOUT
      integer*8, parameter :: rate = 1000
      integer :: i,j, a1, a2

! option 1. most straightforward, but slower
      ! do i=1,N_Preisach
      !    do j = 1,i
      !       if ((d < last_d) .and.  (d < preisach_alpha(i))) states(i,j) = .true.
      !       if ((d > last_d) .and.  (d > preisach_beta(j)))  states(i,j) = .false.
      !    end do
      ! end
      !do

!option 2. faster, but messy
       ! if (d < last_d) then
       !    do  i = 1, N_Preisach
       !       if (d < preisach_alpha(i)) then
       !          do j = 1,i
       !             states(i,j) = .true.                   
       !          end do
       !       else
       !          exit !alpha is monotonic, so if one fails, no need to test the rest
       !       end if
       !    end do
       ! else         
       !    do  j = N_Preisach, 1, -1
       !       if (d > preisach_beta(j)) then
       !          do i = j, N_Preisach
       !             states(i,j) = .false.                                      
       !          end do
       !       else
       !          exit
       !       end if
       !    end do
       ! end if

!option 3. fastest, but complicated
      force = force_prev
       if (d < last_d) then
          a1 = find( preisach_alpha, d )
          a2 = find( preisach_alpha, last_d)

          do  i = max(a2-1, 1), min(a1+1, N_Preisach)
             if (d < preisach_alpha(i)) then
                do j = 1,i
                   if (.not. states(i,j)) then
                      force = force + weights_dim(i,j)
                      states(i,j) = .true.                                      
                   end if
                end do
             end if
          end do
       else         
          a1 = find( preisach_beta, d )
          a2 = find( preisach_beta, last_d)
         
          do  j = max(a1-1, 1), min(a2+1, N_Preisach)

             if (d > preisach_beta(j)) then
                do i = j, N_Preisach
                   if (states(i,j)) then
                      states(i,j) = .false.           
                      force = force - weights_dim(i,j)                           
                   end if
                end do
             end if
          end do
       end if

!!! end options
!this sum accounts for essentially all of the time of the routine.
!yet, if we do without it, roundoff errors will eventually accumulate.
!so do it once every in awhile
       if (mod(IOUT, rate) == 0) then
          force = sum( weights_dim, states)  
       end if

      last_d = d

      force_prev = force

      last_alpha = find( preisach_alpha, last_d)
      last_beta = find( preisach_beta, last_d)

    end subroutine UpdatePreisach_outsideRES1

  end module PreisachModule

  module tip_data
      !tip
      real*8  Rtip_dim, Rtip_nd, Etip, Poisson_tip, tip_angle
      integer tip_shape
  end module tip_data

  module HertzModule
    use tip_data
    use params
    contains
      real*8 function Fts_Hertz(gap, CDMT)
        real*8, intent(in) :: gap, CDMT

   
        if (gap >= 0) then
           Fts_Hertz = 0
        else
           if (tip_shape == PARABOLOID) then
              Fts_Hertz = CDMT*(-gap)**(1.5d0)
           else
              Fts_Hertz = CDMT*(-gap)**(2d0)
           end if
        end if

      end function Fts_Hertz

        !find the maximum von mises stress in an elastic hertz contact (in order to see if we are yielding)
        !reference: Johnson, Contact mechanics.
        !fixme.  right now only works for spherical tip shape
        pure real*8 function Hertz_vonMises(d, F, nu)
          use params 

          real*8, intent(in) :: d, F, nu
          real*8 :: a, p0, z, sigma_r, sigma_th, sigma_z

          if (tip_shape /= PARABOLOID) then
             Hertz_vonMises =0d0
             return
          end if

          if (d > 0) then
             Hertz_vonMises = 0d0
          else
             !first find contact area a and max pressure p0 from indentation and force
             a = sqrt( -d * Rtip_nd)
             p0 = 3d0 * F / ( 2d0 * pi * a**2)

             !then find the z depth at which maximum stress occurs
             !the z depth is a complicated function of poisson's ratio. just work it out in maple and then curve fit
             z = a * ( .328013488739502d0*nu+.381494142661074d0+0.747980854650407d-2*nu**2 )
             

             !then find stresses at that depth
             sigma_r  = p0 * (-(1d0+nu)*(1d0-z*atan(a/z)/a)+0.50/(1d0+z**2/a**2))          
             sigma_th = sigma_r
             sigma_z  = -p0/(1d0+z**2/a**2)          
             

             !then find von mises
             Hertz_vonMises = sqrt( 0.5d0 * ( ( sigma_r-sigma_th)**2 + ( sigma_r-sigma_z)**2 +( sigma_z-sigma_th)**2))
             
          end if
        end function Hertz_vonMises


    end module HertzModule

!you are not expected to understand this module.  Take a look at the developers notes
!first and make sure you understand the reference matlab implementation.  then come back
!here and look at the FORTRAN code.
!
!this module implements the solution of Ting for Hertz contact.  attractive forces are hacked in
module Viscoelasticity_ting
  use Stack
  use params
  use matl_prop_type_module

  type(stack_type), private :: t_stack, d_stack, dpr_stack, deriv_force_stack
  real*8, private :: current_time, current_d, current_dpr
  integer,private :: direction
  integer, parameter :: INC = 1, DEC = 2, D_NONE = 0, LOST_CONTACT = 3

  integer, private :: t1_ndx,  last_t1_ndx, max_ndx_outsideRES1
  
  integer, private :: VEchoice
  real*8, private :: current_force

  !for determining if the sample is relaxed on the next tap
  real*8, private :: t_last_contact, d_last_contact
  real*8 :: relaxed_d

  real*8 :: current_contact_radius !bit of a hack, gets updated when force is computed.  so when we use it we have 
  !to assume that force has been computed recently enough that this is up to date.

! 10,000 points is not enough for a linearly spaced table when tau spans six decades.  because we force the table
! to go out to exp(-10) for the longest tau, the shortest has no resolution.
! e.g. if times span from 500s to 1e-4s, then table goes out to 5000s, 
! and delta_t = 0.5.  then the resolution on the shortest one is exp(-0 / 1e-4) = exp(0)
! straight to exp(-0.5 / 1e-4) = exp(-5000) = which might as well be zero.
! 
! only solution is to log space the table

  integer, private, parameter :: gen_maxwell_lookup_N = 50000
  integer*8 :: counter  !for debugging only

  real*8 :: storage_deltat

  real*8, allocatable :: psi_vec(:)
  integer :: psi_vec_len
  
  private pre_compute_genmaxwell_psi, pre_compute_exp, check_left_boundary, residual, psi, psi_sc,  deriv_force,  getViscoelasticForce_Hertz_KelvinVoigt

  contains

    logical function haveHertzViscoelas(cur_props)
      type(matl_prop), intent(in) :: cur_props
      haveHertzViscoelas = (cur_props%fts_ve_model == VEM_HERTZ_BASED).and. (cur_props%VEchoice /= VEC_NONE)
    end function haveHertzViscoelas


    logical function haveHardHertzViscoelas(cur_props)
      type(matl_prop), intent(in) :: cur_props
      haveHardHertzViscoelas = (cur_props%fts_ve_model == VEM_HERTZ_BASED).and.(( cur_props%VEchoice /= VEC_NONE) .and. (cur_props%VEchoice /= VEC_KELVINVOIGT))
    end function haveHardHertzViscoelas

    subroutine InitViscoelasticModels

    end subroutine InitViscoelasticModels
    
    subroutine NonDimensionalizeViscoelasticModels(props)
      use Elasticity_Formulas
      use NonDimensionalization
      use checkFZCurves
      use timeAndCycle, only : operating_mode, delt
      use data1, only : fexcite
      use data0, only: omegai_nd
      use Stack, only: init_stack

      type(matl_prop), intent(inout) :: props
      real*8 :: G1, G2
      integer :: i

      if ( .not. allocated( psi_vec) ) then
         psi_vec_len = 1000
         allocate( psi_vec( psi_vec_len))
      end if

      call init_stack(t_stack)
      call init_stack(d_stack)
      call init_stack(dpr_stack)
      call init_stack(deriv_force_stack)

      
      counter = 0

      direction = D_NONE
      
      if (props%fts_ve_model == VEM_HERTZ_BASED) then

         !tip is assumed to be rigid for viscoelasticity
         if (props%input_shear_modulus) then
            !a bit of a hack.  the same input box could be two different things
            props%Estar           = ShearToReducedModulus( props%Esample, props%Poisson_sample)
         else
            props%Estar           = reducedModulus1( props%Esample, props%Poisson_sample)
         end if
         
         if (props%VEchoice == VEC_MAXWELL) then
            if (props%input_shear_modulus) then
               G1 = props%Esample
            else
               G1 = YoungsToShear(props%Esample, props%Poisson_sample)
            end if
            
            props%maxwell_t1 =  NonDimenTime(props%etasample / G1)    
         elseif (props%VEchoice == VEC_THREEELM_E1E2ETA) then
            if (props%input_shear_modulus) then
               G1 = props%Esample 
               G2 = props%threeelm_E2
               props%threeelm_E2star = ShearToReducedModulus( props%threeelm_E2, props%Poisson_sample)
            else
               G1 = YoungsToShear(props%Esample , props%Poisson_sample)
               G2 = YoungsToShear(props%threeelm_E2, props%Poisson_sample)
               props%threeelm_E2star = reducedModulus1( props%threeelm_E2, props%Poisson_sample)
            end if
                        
            props%threeelm_surf_tau = NonDimenTime( props%etasample / G2)
            
            props%threeelm_tau_stress = NonDimenTime( props%etasample / ( G1 + G2))
            props%threeelm_e1e2 =  props%Esample / props%threeelm_E2
         elseif (props%VEchoice == VEC_GENMAXWELL) then
            
            do i = 1, (props%N_gen_max)
               if (props%input_shear_modulus) then
                  !kind of a hack.  "E" is really shear modulus here.
                  props%E_gen_max_star(i) = ShearToReducedModulus( props%E_gen_max(i), props%Poisson_sample)
               else                  
                  props%E_gen_max_star(i) = ReducedModulus1( props%E_gen_max(i), props%Poisson_sample)
               end if
               
               props%tau_gen_max(i) = NonDimenTime(props%tau_gen_max_dim(i))
            end do
            call  pre_compute_genmaxwell_psi( props )
         elseif (( props%VEchoice == VEC_KELVINVOIGT) .or.  ( props%VEchoice == VEC_KELVINVOIGT_LEGACY)) then
            !do nothing
         elseif (props%VEchoice /= VEC_NONE) then
            call assert(.false., 'unknwon ven model')
         end if
      elseif  (props%fts_ve_model == VEM_LINEAR_BASED) then !linear based
         if (props%VEchoice == VEC_THREEELM_E1E2ETA) then                   
            props%threeelm_tau_stress =  NonDimenTime(props%etasample / (props%kts_R+props%k2))
            call pre_compute_exp(props%threeelm_tau_stress, props )
         elseif (props%VEchoice == VEC_MAXWELL) then
            props%maxwell_t1 = NonDimenTime(props%etasample / props%kts_R)
            call pre_compute_exp(props%maxwell_t1, props )
         end if
      else
         call WriteFatalError('Shouldnt be here!')
      end if

      if ((isFzCurves(operating_mode, fexcite)) .or. (operating_mode == APPROACH_TRIGGERED)) then
         !storage_deltat = max(minimum_tau(props) / 10, delt)   ! original code. only problem was that delt had apparently not been computed yet.  not actually a problem for storage_deltat to be less than one time step.  it just results in storing every point

         !needs to be at least 10 points per relaxation time for accuracy. but also at least four points per oscillation cycle.
         storage_deltat = min( minimum_tau(props) / 10d0,( 1d0 / (omegai_nd(1)/(2d0*pi))) / 4d0 )
      else
         storage_deltat = 0d0
      endif


    end subroutine NonDimensionalizeViscoelasticModels




    real*8 function ViscoelasticModelsMaxDeltaT(props)
      type(matl_prop), intent(in) :: props
      
      if (props%fts_ve_model == VEM_HERTZ_BASED) then
         !tip is assumed to be rigid for viscoelasticity
         if (props%VEchoice == VEC_MAXWELL) then
            ViscoelasticModelsMaxDeltaT = props%maxwell_t1
         elseif (props%VEchoice == VEC_THREEELM_E1E2ETA) then
            !three elm will be a little more forgiving than maxwell,
            !b/c psi will relax to a finite value not zero.  still
            !if the relaxation time is very small, then why have it at all            
            ViscoelasticModelsMaxDeltaT =  10d0 * props%threeelm_tau_stress
         elseif (props%VEchoice == VEC_GENMAXWELL) then
            ViscoelasticModelsMaxDeltaT = maxval(props%tau_gen_max) * 10
         else
            call assert(.false., 'dont know how to compute max relaxation time')
            ViscoelasticModelsMaxDeltaT = 1d99
         end if
      elseif  (props%fts_ve_model == VEM_LINEAR_BASED)  then
         if (props%VEchoice == VEC_THREEELM_E1E2ETA) then
            ViscoelasticModelsMaxDeltaT = 10d0 * props%threeelm_tau_stress
         else
            call assert(.false., 'dont know how to compute max relaxation time')
            ViscoelasticModelsMaxDeltaT = 1d99
         end if
      elseif  (props%fts_ve_model == VEM_ATTARD_BASED)  then         
         if ((props%VEchoice == VEC_THREEELM_E1E2ETA) .or. (props%VEchoice == VEC_THREEELM_E0_EINF_TAU)) then
            ViscoelasticModelsMaxDeltaT =  props%threeelm_tau_stress            
         else
            call assert(.false., 'dont know how to compute max relaxation time')
            ViscoelasticModelsMaxDeltaT = 1d99
         end if
      elseif  (props%fts_ve_model == VEM_NONE)  then
         ViscoelasticModelsMaxDeltaT = 1d99
      else
         call assert(.false., 'dont know how to compute max relaxation time')
         ViscoelasticModelsMaxDeltaT = 1d99         
      end if
    end function ViscoelasticModelsMaxDeltaT

    
    !originally I had this with the non-dimensionalization above.  but at the time that I wanted to do that,
    !I had not calculated delt yet.  so here is just for checking the time constants
    subroutine CheckViscoelasticModels(props, delt)
      use NonDimensionalization
      real*8, intent(in) :: delt
      type(matl_prop), intent(in) :: props
      real*8 :: maxDeltaT
      character*1000 :: tmp

      maxDeltaT = ViscoelasticModelsMaxDeltaT(props)
      
      if (delt > maxDeltaT ) then
         write(tmp,*) "Solver time step ( = ", DimenTime(delt), ") is too large relative to the viscoelastic relaxation time.  Maximum allowable time step is (=", DimenTime(maxDeltaT), '). &
However, a smaller timestep is generally required for accurate results.   You can either decrease the time step (by increasing the number of deflection points per cycle &
on simulation tab) or adjust the material properties to increase the relaxation times (e.g. increase viscosity, or decrease youngs modulus).&
Keep in mind that typically, relaxation times far away from the contact time generally do not affect the tip-sample force.'
         call WriteFatalError(tmp)
      end if
    end subroutine CheckViscoelasticModels

    

  subroutine Update_viscoelasticity_outsideRes1_linear( t, d, cur_props)
    use NonDimensionalization
    real*8, intent(in) :: t, d
    type(matl_prop), intent(in) :: cur_props

    integer :: m

    if (.not. ((cur_props%fts_ve_model == VEM_LINEAR_BASED) .and. (cur_props%VEchoice == VEC_THREEELM_E1E2ETA .or. cur_props%VEchoice == VEC_MAXWELL))) return

    if (d > 0) then
       !not in contact.  clear anything out that might be there.
       call clear_stack(t_stack)
       call clear_stack(d_stack)

       current_force = 0d0
    else
       !in contact
       call push(t,t_stack)
       call push(d, d_stack)
    end if
    
  end subroutine Update_viscoelasticity_outsideRes1_linear

  real*8 function GetViscoelasticForce_linear(t, d, operating_mode, cur_props)
    use Integrals
    use NonDimensionalization
    use timeAndCycle, only : delt
    use Interpolation
    type(matl_prop), intent(in) :: cur_props
    real*8, intent(in) :: t, d
    integer, intent(in) :: operating_mode
    real*8 :: p2,  q, t0
    integer :: n, i

    if (( t == current_time) .and. (d == current_d)) then
       getViscoelasticForce_linear = current_force
    else
       current_time = t
       current_d = d
       
       !if we get here we're probably inside RES1.  do the full calculations but 
       !don't update the history stack

       if (d > 0) then
          getViscoelasticForce_linear = 0d0
       else
          n = length(t_stack)
          if (n == 0) then
             !fixme, something here!
             getViscoelasticForce_linear = 0 
          else
             t0 = t_stack%data(1) !the time we entered contact

             if (( cur_props%VEchoice == VEC_THREEELM_E1E2ETA) .or. (cur_props%VEchoice == VEC_MAXWELL )) then
                !original way as I got it from linear_3_elm.mw
                ! p = exp( (t_stack%data(1:n) - t0) / cur_props%threeelm_tau_stress)
                ! q = trapz( d_stack%data(1:n) * p(1:n), n,  delt, t_stack%data, .false.) 
                !  getViscoelasticForce_linear = - NonDimenStiffness(cur_props%kts_R) * ( d - NonDimenFreq(cur_props%kts_R / cur_props%etasample) * q * exp( -(t - t0) /  cur_props%threeelm_tau_stress))

                !after some rearranging and optimizing.  analytically the same, but numerically much better behaved and 30x faster.
                if (n > psi_vec_len) then
                   psi_vec_len = max( psi_vec_len * 5, 2 * n)
                   deallocate(psi_vec)
                   allocate( psi_vec( psi_vec_len))
                end if
                do i = 1,n
                   psi_vec(i) =  interpolate_w_clip(cur_props%gen_maxwell_t_lookup, cur_props%gen_maxwell_psi_lookup, -(t_stack%data(i)-t), gen_maxwell_lookup_N, 0d0, 1d0 )
                end do

                q = trapz( d_stack%data(1:n) * psi_vec(1:n), n,  delt, t_stack%data, .false.) !not necessarily equally spaced because of the root finding! 
                getViscoelasticForce_linear = - NonDimenStiffness(cur_props%kts_R) * ( d - NonDimenFreq(cur_props%kts_R / cur_props%etasample) * q)
!             elseif ( cur_props%VEchoice == VEC_MAXWELL ) then
!                !actually exact same code as 3 elm model, except k2 = 0, so we use maxwell_t1 instead of threelem_t1. but that has been
!                !factored out by the interpolation above.
!                !p = exp( (t_stack%data(1:n)-t) / cur_props%maxwell_t1)
!                q = trapz( d_stack%data(1:n) * p(1:n), n,  delt, t_stack%data, .false.) !not necessarily equally spaced because of the root finding! 
!                getViscoelasticForce_linear = - NonDimenStiffness(cur_props%kts_R) * ( d - NonDimenFreq(cur_props%kts_R / cur_props%etasample) * q)
             end if


             if (getViscoelasticForce_linear < 0) getViscoelasticForce_linear = 0
          end if
       end if
    end if

  end function GetViscoelasticForce_linear

!fixme, all of the different viscoelastic models are getting to be a pain to keep track of.
!need a better organization.
  subroutine Update_viscoelasticity_outsideRes1_Hertz( t, d, dpr, cur_props, Rtip_nd)
    use NonDimensionalization
    real*8, intent(in) :: t, d, dpr, Rtip_nd
    type(matl_prop), intent(in) :: cur_props
    real*8 :: relaxed_t, tau, G2

    if (.not. ((cur_props%fts_ve_model == VEM_HERTZ_BASED) .and. (( cur_props%VEchoice == VEC_GENMAXWELL) .or. ( cur_props%VEchoice == VEC_MAXWELL) .or. ( cur_props%VEchoice == VEC_THREEELM_E1E2ETA))))  return

    call assert( .not. isnan(d) )       

    if (d > 0) then
       !not in contact.  clear anything out that might be there.
       call clear_stack(t_stack)
       call clear_stack(d_stack)
       call clear_stack(dpr_stack)
       call clear_stack(deriv_force_stack) 
       current_force = 0d0
       direction = D_NONE
       t1_ndx = -1
       last_t1_ndx = -1
       max_ndx_outsideRES1 = -1
    elseif ((length(t_stack)>0) .and. (direction == DEC) .and. ( current_force == 0d0)) then
       !d < 0, but force is zero.  lost contact
       direction = LOST_CONTACT
       t_last_contact = t
       d_last_contact = d
    else
       !in contact

       !update the history stack
       if ( dpr < 0)  then
          !gap decreasing => indentation increasing => contact area is increasing.
          if (direction == D_NONE) then
             direction = INC
             
             !this is the first point in contact.  Check if the surface has fully relaxed from the previous tap
             if (cur_props%VEchoice == VEC_THREEELM_E1E2ETA) then
                relaxed_t = t - t_last_contact
                relaxed_d = d_last_contact * exp( - relaxed_t / cur_props%threeelm_surf_tau)
             end if

          elseif ((direction == LOST_CONTACT) .or. (direction == DEC)) then
             call WriteFatalError('Indentation reversal in viscoelastic contact is not implemented. Consult the manual section on viscoelasticity for details.');
          end if
       else
          !contact area is decreasing.  
          if (direction == INC) then
             !this is the first decreasing point.
             direction = DEC          
             max_ndx_outsideRES1 = length(t_stack)
          end if

       end if

       !used to always store every point.  But for f-z curves might need to  be store 1000 points/cycle for 1000 cycles,
       !even though most of that information was redundant.  So for f-z curves, only store at a reasonable rate (the
       !rate calculation is elsewhere).  for tapping mode, storage_delt=0.  this bit of a kludge here is needed b/c
       !fortran refuses to do short circuit evaluation.
       if (length(t_stack) > 100) then
          if ( (t-head(t_stack)) < storage_deltat) return
       end if
       call push(t, t_stack)
       call push(d, d_stack)
       call push(dpr, dpr_stack)    
       call push( deriv_force(d,dpr) , deriv_force_stack)     
    end if ! d>=0    
  end subroutine Update_viscoelasticity_outsideRes1_Hertz

  !this just includes the dependence on d, and dpr.  all of the constants are already included elsewhere
  pure real*8 function deriv_force( d, dpr)
    use tip_data
    real*8, intent(in) :: d, dpr

    if (tip_shape == PARABOLOID) then
       !f ~ d^(3/2), so f' ~ d' sqrt(d)
       deriv_force = dpr * sqrt(-d)
    else !cone
       !f ~ d^2, so     f' ~ d' d
       deriv_force = dpr * (-d)
    end if
  end function deriv_force

  !this is the elastic relation between indentation and contact area
  pure real*8 function contact_radius(d)
    use tip_data
    real*8, intent(in) :: d

    if (tip_shape == PARABOLOID) then
       contact_radius = sqrt( - d * Rtip_nd)
    else !cone
       contact_radius = - d * (2d0 / pi) / tan(tip_angle)
    end if
    
  end function contact_radius

  !this function is very complicated and I apologize for that... but it
  !must be this complicated in order to be computed quickly, and to 
  !be able to handle intermediate points from DDASKR.
  !
  !go look at the developer's notes and make sure you understand the reference
  !matlab implementation before trying to understand this code.
  real*8 function getViscoelasticForce_Hertz(t, d, dpr, cur_props, cur_fC, operating_mode)
    !d already has aDMT subtracted out
    use Integrals
    use Interpolation
    use timeAndCycle, only : delt
    use tip_data
    use HertzModule

    real*8, intent(in) :: t, d, dpr
    integer, intent(in) :: operating_mode
    type(matl_prop), intent(in) :: cur_props
    type(forceCoeff), intent(in) :: cur_fC
    integer :: i, search_ndx, n, max_ndx, slow_t1_ndx
    real*8 :: t1, residual_left, residual_right, tmp
    real*8 ::  dt, deriv_force_t1
!    real*8 :: foo !for debugging
    logical :: have_interpolation, no_history

!    counter = counter+1 !for debugging 

    have_interpolation = .false.
    
    if (direction == LOST_CONTACT) then
       getViscoelasticForce_Hertz = 0d0
    elseif (( t == current_time) .and. (d == current_d) .and. (dpr == current_dpr)) then
       !these calculations are very expensive, and sometimes we get called repeatedly with the same arguments.
       !memorize the result so we don't have to recalculate
       getViscoelasticForce_Hertz = current_force
    else
       current_time = t
       current_d = d
       current_dpr = dpr

       !if we get here we're probably inside RES1.  do the full calculations but 
       !don't update the history stack

       if (d > 0) then
          getViscoelasticForce_Hertz = 0d0
       else
          if ( cur_props%VEchoice == VEC_KELVINVOIGT) then
             current_force = getViscoelasticForce_Hertz_KelvinVoigt(t, d, dpr, cur_props, cur_fC)
          else
             n = length(t_stack)

             if ( dpr <= 0) then
                !indentation increasing, contact area increasing
                !can use elastic formula for contact area then calculate force

                no_history = .false.
                if (n == 0) then
                   no_history = .true.
                elseif (n == 1) then
                   if  (t == head(t_stack)) then  !fortran seems to lack short circuit evaluation, forcing this awkwardness
                      no_history = .true.
                   end if
                end if

                if (no_history) then
                   !no previous indentation history. guess at how long we've been in contact
                   !and do trapezoidal int from 0,0.
                   if (dpr == 0) then
                      dt = 0 !will give zero force.  what else can we do?
                   else
                      dt = d / dpr
                   end if
                   tmp = -psi_sc( 0d0, cur_props ) * deriv_force(d,dpr) * dt / 2
                else
                   tmp = trapz_p1( - psi( t- t_stack%data(1:n),n, cur_props) * deriv_force_stack%data(1:n), n, delt, t_stack%data, .false., -psi_sc( 0d0, cur_props ) * deriv_force(d,dpr), t - head(t_stack) )  !not necesarily equally spaced because of the root finding.

                   !thought this would be better, but actually makes the noise a little worse
!                   tmp = simpson_p1( - psi( t- t_stack%data(1:n),n, cur_props) * dpr_sqrt_d_stack%data(1:n), n, delt, t_stack%data, -psi_sc( 0d0, cur_props ) * sqrt(-d) * dpr, t - head(t_stack) )

                   !thought this would be better, but actually makes the noise a little worse
!                   tmp = simpson_p1( - psi( t- t_stack%data(1:n),n, cur_props) * dpr_sqrt_d_stack%data(1:n), n, delt, t_stack%data, -psi_sc( 0d0, cur_props ) * sqrt(-d) * dpr, t - head(t_stack) )
                end if

                call assert( tmp >= 0, 'tmp < 0')

                current_force =  cur_fC%vec_f * tmp

                current_contact_radius = contact_radius(d)

             elseif ( (dpr > 0) .and. ( n == 0)) then
                !indentation decreasing yet somehow no history of previous indentation.
                !implies coming in to and then out of contact within the same integration point
                !just call that no force I guess?             
                current_force = 0d0
                current_contact_radius = 0d0
             else
                if ( max_ndx_outsideRES1 == -1) then
                   max_ndx = n
                else
                   max_ndx = max_ndx_outsideRES1  !I don't remember what this was for?
                end if

                !calculation of psi gets very expensive.  precompute and re-use values is a must
                if (n > psi_vec_len) then
                   psi_vec_len = max(psi_vec_len * 5, n*2)
                   deallocate(psi_vec)
                   allocate( psi_vec( psi_vec_len))
                end if                
                psi_vec(1:n) = psi(t - t_stack%data(1:n), n, cur_props)

!                slow_t1_ndx = find_zero(psi_vec, t, dpr, cur_props, max_ndx)  !for debugging only. very slow!

                !indentation decreasing.  cannot use elastic formula for contact area.
                !we find contact area implictly by finding t1, the time at which the contact area 
                !previously took on the current value of contact area.
                !we do this by stepping through the stack and evaluating an integral at each point
                !the integral is zero precisely at t1, so try to find a zero crossing. The integrals are expensive to compute, so we
                !want to guess a value as close as possible.

                if (last_t1_ndx == -1) then
                   !this is the first decreasing point.  just go back one point.
                   search_ndx = max_ndx - 1
                else
                   !start new search where the last search stopped
                   !except if the last search stopped at 1, give ourselves a chance to go back to 2
                   !otherwise we might get stuck against the wall b/c of intermediate points within DDASKR
                   search_ndx = min( max(last_t1_ndx,2), max_ndx - 1)                   
                end if

                if (search_ndx <= 1) then
                   !if we start at search=1, then we can't check to the left.  already as far as we can go
                   t1_ndx = check_left_boundary(psi_vec, t, dpr, cur_props,  t1, have_interpolation)
                else
                   t1_ndx = -1

                   residual_left  = residual(psi_vec, search_ndx,  t, dpr, cur_props)
                   residual_right = residual(psi_vec, search_ndx+1,t, dpr, cur_props)
                end if

!fixme, might need some way to ensure we don't get stuck in infinite loop

                do while (t1_ndx == -1)
                   if (sign(1d0, residual_left) /= sign(1d0, residual_right)) then
                      ! if we stop here, that's nearest neighbor interpolation

                      t1_ndx = search_ndx

                      ! this interpolation should improve accuracy,
                      t1 = interpolate( (/residual_left, residual_right/) , (/t_stack%data(t1_ndx),t_stack%data(t1_ndx+1)/), 0d0, 2)
                      ! call assert( t1 > t_stack%data(t1_ndx-1), 't1 > t1_ndx-1')
                      ! call assert( t1 < t_stack%data(t1_ndx+1), 't1 > t1_ndx-1')
                      !call assert( t1_ndx >= 2, 't1_ndx<2')
                      !call assert( t1_ndx <= n-1, 't1_ndx>n-1')
                      deriv_force_t1 =  interpolate(  (/t_stack%data(t1_ndx),t_stack%data(t1_ndx+1)/), (/deriv_force_stack%data(t1_ndx),deriv_force_stack%data(t1_ndx+1)/), t1, 2)
                      if (t1 < t_stack%data(t1_ndx)) t1_ndx = t1_ndx - 1
                      have_interpolation = .true.

                   elseif ( abs(residual_left) < abs(residual_right)) then
                      !if left is smaller,  go left.  
                      if ( search_ndx <= 1) then
                         ! can't go any farther left
                         t1_ndx = check_left_boundary(psi_vec, t, dpr, cur_props, t1, have_interpolation)
                         if (have_interpolation) deriv_force_t1 =  interpolate(  (/t_stack%data(1),t_stack%data(2)/), (/deriv_force_stack%data(1),deriv_force_stack%data(2)/), t1, 2)
                      else
                         search_ndx = search_ndx - 1
                         residual_right = residual_left
                         residual_left = residual(psi_vec, search_ndx,t, dpr, cur_props)
                      end if
                   else
                      !else go right
                      if (search_ndx == max_ndx-1) then
                         !cannot go farther. !fixme. may want to interpolate here.
                         t1_ndx = max_ndx                                 
                      else
                         search_ndx = search_ndx + 1
                         residual_left = residual_right
                         residual_right = residual(psi_vec, search_ndx+1,t, dpr, cur_props)
                      end if
                   end if
                end do

             !   call assert( (slow_t1_ndx == t1_ndx) .or. (t1_ndx == 0) , 't1 search failed')

                last_t1_ndx = t1_ndx
                !force is found by same formula as for increasing contact area, but limits of integration are changed 
                if (t1_ndx > 0) then
                   if (have_interpolation) then                                                          
                      tmp = trapz_p1( -psi_vec(1:t1_ndx) * deriv_force_stack%data(1:t1_ndx), t1_ndx, delt, t_stack%data, .false., -psi_sc(t-t1, cur_props) * deriv_force_t1,  t1 -t_stack%data(t1_ndx))
                      current_contact_radius = contact_radius( interpolate( t_stack%data, d_stack%data, t1, n) )
                   else
                      tmp = trapz( -psi_vec(1:t1_ndx) * deriv_force_stack%data(1:t1_ndx), t1_ndx, delt, t_stack%data, .false.)
                      current_contact_radius = contact_radius( d_stack%data(t1_ndx) )
                   end if
                else
                   tmp = 0d0
                   current_contact_radius = 0d0
                end if
              
                current_force =  cur_fC%vec_f * tmp
              
             end if
          end if 

          getViscoelasticForce_Hertz = current_force
          call assert (getViscoelasticForce_Hertz >= 0d0, 'visc force negative')

       end if ! d>0
    end if ! t == current_time
  end function getViscoelasticForce_Hertz
  
  real*8 function getViscoelasticForce_Hertz_KelvinVoigt(t, d, dpr, cur_props, cur_fC)
    use tip_data
    use HertzModule
    real*8, intent(in) :: t, d, dpr
    type(matl_prop), intent(in) :: cur_props
    type(forceCoeff), intent(in) :: cur_fC
    real*8 d_to_find
    

    if (dpr <= 0) then
       !indentation increasing, contact area increasing
       !can use elastic formula for contact area then calculate force.
       !Dts is not the same as in ad-hoc though.  factors of 2/(1-nu).

       if (tip_shape == PARABOLOID) then
          getViscoelasticForce_Hertz_KelvinVoigt = Fts_Hertz(d, cur_fC%CDMT) - cur_fC%Dts*sqrt(-d)*dpr
       else
          call WriteFatalError('conical tip not implemented for KV viscoelasticity choice yet')
       end if
       current_contact_radius = contact_radius( d) 
    else
       !indentation decreasing, contact area decreasing.
       !cannot use elastic formula for contact area.
       !but, we don't need to do the whole integral business.
       !we can find solution in closed form

       d_to_find = cur_fC%kv_t0 * dpr + d
       if (d_to_find > 0) then
          !think this means that we will leave contact                   
          getViscoelasticForce_Hertz_KelvinVoigt = 0d0                   
       else
          !then compute force
          getViscoelasticForce_Hertz_KelvinVoigt = Fts_Hertz(d_to_find, cur_fC%CDMT)
       end if
       current_contact_radius = contact_radius( d_to_find) 
    end if
  end function getViscoelasticForce_Hertz_KelvinVoigt


  ! !use this for debuging only.  straightforward but very slow!
  ! integer function find_zero(psi_vec, t, dpr, cur_props, max_ndx)
  !   integer, intent(in) :: max_ndx
  !   type(matl_prop), intent(in) :: cur_props
  !   real*8,intent(in) ::  psi_vec(:), t, dpr
  !   integer :: i
  !   real*8 :: residual_table(max_len), abs_residual_table(max_len)
    
  !   do i=1,max_ndx
  !      residual_table(i)  =  residual(psi_vec, i, t, dpr, cur_props)
  !      abs_residual_table(i) = abs(residual_table(i))
  !   end do

  !   do i = 1,max_ndx-1
  !      if  (sign(1d0, residual_table(i)) /= sign(1d0, residual_table(i+1))) then
  !         find_zero = i
  !         return
  !      end if
  !   end do
          
  !   !if we get there, there is no zero crossing.  forced to just take the minimum absolute value
  !   find_zero = minloc( abs_residual_table(1:max_ndx),1 )
  ! end function find_zero


  !we backed up all the way to ndx=1.  can't go any farther.  check what to do next.
  integer function check_left_boundary(psi_vec, t, dpr, cur_props, t1, have_interpolation)
    use Interpolation
    real*8, intent(in) :: t, dpr
    type(matl_prop), intent(in) :: cur_props
    real*8 :: residual_left, residual_right, psi_vec(:)
    real*8, intent(out) :: t1
    logical, intent(out) :: have_interpolation 

    have_interpolation = .false.

    if ( length(t_stack) == 1) then
       !just one point on the stack.  can't properly tell what to do here without interpolating.
       !for now just assume that we're on that one point.
       check_left_boundary = 1
       return
    end if

    residual_left =  residual( psi_vec, 1, t, dpr, cur_props)
    residual_right = residual( psi_vec, 2, t, dpr, cur_props)
    
    if (sign(1d0,residual_left) /= sign(1d0,residual_right)) then
       !there is a sign change.  implies that the zero must be between these two points
       check_left_boundary = 1
       t1 = interpolate( (/residual_left, residual_right/), (/t_stack%data(1), t_stack%data(2)/), 0d0, 2)
!       call assert( t1 >= t_stack%data(1), 't1 > data(1)')
!       call assert( t1 <= t_stack%data(2), 't1 < data(2)')
       have_interpolation = .true.
    else
       !there is no sign change.  implies that there is no zero.  we must be leaving contact
       check_left_boundary = 0
    end if
  end function check_left_boundary

  real*8 function residual(psi_vec, t1_ndx, t, dpr, cur_props)
    use Integrals
    use timeAndCycle, only : delt
    type(matl_prop), intent(in) :: cur_props
    integer, intent(in) :: t1_ndx
    real*8, intent(in) :: dpr, t, psi_vec(:)
    integer :: n
    
    n = length(t_stack)
    call assert( t1_ndx <= n, 'residual: t1 > length')
    call assert( t1_ndx > 0, 't1 < 0')

    residual = trapz_p1( psi_vec(t1_ndx:n) * dpr_stack%data( t1_ndx:n) , n-t1_ndx+1, delt,  t_stack%data( t1_ndx:n), .false. , psi_sc(0d0, cur_props) * dpr, t - head(t_stack))

  end function residual

  !this is the stress relaxation function for the given material (e.g. the force at we would feel as a function
  !of time in response to a step change in indentation).  this is vector form
  function psi(t, n, cur_props )
    type(matl_prop), intent(in) :: cur_props
    integer, intent(in) :: n
    real*8, intent(in) :: t(n)
    real*8 :: psi(n)
    integer :: i
    
    do i = 1,n
       psi(i) = psi_sc(t(i), cur_props)
    end do    
  end function psi

  real*8 function psi_sc(t, cur_props )
    use Interpolation
    type(matl_prop), intent(in) :: cur_props
    real*8, intent(in) :: t
    
    if (cur_props%VEchoice == VEC_MAXWELL) then
       !this is maxwell
       !there is a leading scalar, but it's a constant so we factor it out of here and put it into the cur_fC 
       !later when we compute the forces.  for computing the residual we are looking for minimum anyway so
       !a constant factor is irrelevant
       psi_sc = exp(- t / cur_props%maxwell_t1 )    
    elseif (cur_props%VEchoice == VEC_THREEELM_E1E2ETA) then
       !again, the leading scalars is factored out
       !ps = (E1/(E1+E2)) * ( E2 + E1 * exp((tau-t)/T));
       !   = (E1/(E1+E2)) * E2 ( 1 + E1/E2 * exp((tau-t)/T))

       psi_sc = 1 + cur_props%threeelm_e1e2 * exp(  -t / cur_props%threeelm_tau_stress )    
    elseif (cur_props%VEchoice == VEC_GENMAXWELL) then
       !here we don't do any scaling
       !for large numbers of maxwell elements, this gets to be a huge resource hog.  50 - 60% of total
       !cpu time can be this line right here. so instead use the lookup table routine below. much faster
!       psi_sc =    cur_props%Estar +  dot_product(  cur_props%E_gen_max_star(1: cur_props%N_gen_max),  exp( -t /  cur_props%tau_gen_max(1: cur_props%N_gen_max)))

       if ( t >= cur_props%gen_maxwell_t_lookup(gen_maxwell_lookup_N)) then
          !t is so large that this might as well be fully relaxed
          psi_sc = cur_props%Estar
       else
          !we're log spacing the table now.  need a log lookup routine not a linear lookup routine.
!          psi_sc = linear_table_lookup( cur_props%gen_maxwell_t_lookup , cur_props%gen_maxwell_psi_lookup, t, gen_maxwell_lookup_N )
          psi_sc = interpolate( cur_props%gen_maxwell_t_lookup , cur_props%gen_maxwell_psi_lookup, t, gen_maxwell_lookup_N )
       end if
    else
       call WriteFatalError('unknown ve model')
    end if
  end function psi_sc

  real*8 function minimum_tau(cur_props)
    type(matl_prop), intent(inout) :: cur_props
    if (cur_props%VEchoice == VEC_GENMAXWELL) then
       minimum_tau = minval( cur_props%tau_gen_max(1: cur_props%N_gen_max))
    elseif (cur_props%VEchoice == VEC_THREEELM_E1E2ETA) then
       minimum_tau = cur_props%threeelm_tau_stress
    elseif (cur_props%VEchoice == VEC_MAXWELL) then
       minimum_tau = cur_props%maxwell_t1
    else
       minimum_tau = 0
    end if
  end function minimum_tau
  
  subroutine pre_compute_genmaxwell_psi( cur_props )
    use GridGeneration
    type(matl_prop), intent(inout) :: cur_props

    real*8 :: max_tau, min_tau
    integer :: i

    !need a separate look up table for each set of data
    max_tau = maxval( cur_props%tau_gen_max(1: cur_props%N_gen_max))
    min_tau = maxval( cur_props%tau_gen_max(1: cur_props%N_gen_max))

    allocate( cur_props%gen_maxwell_psi_lookup( gen_maxwell_lookup_N ) )
    allocate( cur_props%gen_maxwell_t_lookup( gen_maxwell_lookup_N ) )

    !generate the table so that for the largest possible tau, the table goes out to -10.  exp(-10)=4e-5 which is negligible.
    !cur_props%gen_maxwell_t_lookup(i)   = ( real(i-1) / real(gen_maxwell_lookup_N)) * 10 * max_tau
    !note may underflow
    
    cur_props%gen_maxwell_t_lookup(1) = 0
    cur_props%gen_maxwell_t_lookup(2:gen_maxwell_lookup_N) = logspace( 1e-6 * min_tau, 15 * max_tau, gen_maxwell_lookup_N-1) 
    
    do i = 1, gen_maxwell_lookup_N      
       cur_props%gen_maxwell_psi_lookup(i) =  cur_props%Estar +  dot_product(  cur_props%E_gen_max_star(1: cur_props%N_gen_max),  exp( - cur_props%gen_maxwell_t_lookup(i) /  cur_props%tau_gen_max(1: cur_props%N_gen_max)))
    end do

  end subroutine pre_compute_genmaxwell_psi

!fixme. work in progress. not finished yet.  
  !similar to above, but for piecewise linear contact based. just need to precompute a single exponential.
  !we will re-use the same variables, which is a bit of a hack
  subroutine pre_compute_exp(tau, cur_props )
    use GridGeneration
    real*8, intent(in) :: tau
    type(matl_prop), intent(inout) :: cur_props

    integer :: i

    allocate( cur_props%gen_maxwell_psi_lookup( gen_maxwell_lookup_N ) )
    allocate( cur_props%gen_maxwell_t_lookup( gen_maxwell_lookup_N ) )

    !generate the table so that the table goes out to -20 tau.  
    
    cur_props%gen_maxwell_t_lookup(1) = 0
    cur_props%gen_maxwell_t_lookup(2:gen_maxwell_lookup_N) = logspace( 1d-6 * tau, 20d0 * tau, gen_maxwell_lookup_N-1) 

    do i = 1, gen_maxwell_lookup_N
       cur_props%gen_maxwell_psi_lookup(i) = exp( -cur_props%gen_maxwell_t_lookup(i) /  tau)
    end do

  end subroutine pre_compute_exp

end module Viscoelasticity_ting




!this is the new module, includes the original spatial discretization (with the hysteretic stuff stripped out)
! also includes fourier discretization based on Bahram's works
! only 3 element model is implemented, not generalized maxwell (which I don't think every completely worked anyway)
!
module Viscoelasticity_attard
  use params
  use Nondimensionalization
  use matl_prop_type_module

      
  real*8 :: u0 ! surface deformation at center. with previous spatial discretization, could always grab Y(1), but now we might have to compute based on fourier coeff.  so just store for later use

  real*8 :: attard_start_time
  
  logical attard_enabled
  
  !profiler can't see private functions?
  private :: compute_new_udot, compute_new_adot_lsq, compute_new_adot_bahram

  
  contains

    !performance note: when attard is sudden turned on in the main code, it's a step change to the force (a very small step, but a step)
    !which causes ddaskr to take a lot of extra iterations (100 - 200 even). ramping the force in slowly does not help anything.
    !zero versus non-zero still makes DDASKR slow. Probably would have to go with a non-iterative solver to get any speedup.  
    subroutine CheckAttardStartStop(d,dpr,t, Y, cur_props)
      real*8,intent(in)    :: d, dpr, t
      real*8,intent(inout) :: Y(:) !Y is the surface dof from the ode solver, ie. either u_cur or a_cur, depending on whether we are doing spatial or fourier
      type(matl_prop), intent(in) :: cur_props
      
      if (attard_enabled) then
         !check if we should disable
         if (dpr > 0) then
            if ( (d-u0)  >  cur_props%attard_stop_nd) then
               attard_enabled = .false.               
               Y=0 !reset the surface. this is a BIG assumption.  fixme, do we need to check this?
            end if
         end if
      else
         !check if we should enable
         if (dpr < 0) then
            if (d <  cur_props%attard_start_nd) then
               attard_enabled = .true.
               attard_start_time = t
            end if
         end if
      end if

    end subroutine CheckAttardStartStop

!this didn't help    
!     real*8 function attardRampIn(t)
!       real*8, intent(in) :: t
! !      real*8, parameter :: decay_rate = 2 !arbitrary, but on non-dimensional time seems reasonable
!       if (t < attard_start_time) then
!          attardRampIn = 0  !shouldn't happen really
!       elseif (t == attard_start_time) then
!          attardRampIn  = 0
!       else
!          attardRampIn = exp( -1 / (t - attard_start_time))
!       end if      
!     end function attardRampIn


    
  subroutine initialize_attard( Rtip_nd, cur_props )
    use Derivatives
    use NonDimensionalization
    use GridGeneration
    use Elasticity_Formulas
    integer i,j
    real*8, intent(in) :: Rtip_nd
    type(matl_prop), intent(inout) :: cur_props
    real*8 :: G1,G2
    real*8, allocatable :: gamma1(:,:)       

    attard_enabled = .false.  

    attard_start_time = 0
    
    if (  cur_props%N_attard_fourier >  cur_props%N_attard_spatial) then
       call WriteFatalError('Error: the number of fourier basis functions for Attards model must be less or equal to than the number of spatial discretization points. please check inputs and try again')
    end if
    
    allocate( cur_props%k(cur_props%N_attard_spatial, cur_props%N_attard_spatial))
    
    allocate( cur_props%kt(cur_props%N_attard_spatial, cur_props%N_attard_spatial)) !transposed version
    
    allocate( cur_props%bare_k(cur_props%N_attard_spatial, cur_props%N_attard_spatial))
    
    allocate( cur_props%tip_shape(cur_props%N_attard_spatial))
    allocate( cur_props%rr(cur_props%N_attard_spatial))    
    allocate( cur_props%C(cur_props%N_attard_spatial, cur_props%N_attard_fourier) )

    allocate( cur_props%J(cur_props%N_attard_spatial,cur_props%N_attard_spatial))
    allocate( cur_props%b(cur_props%N_attard_spatial,1))
    allocate( cur_props%ipiv(cur_props%N_attard_spatial)) !work space for lapack routine

    
    cur_props%A_hamaker_nd  = NonDimenEnergy(  cur_props%A_hamaker_dim )

    
    !this duplicates a lot of code from NonDimensionalizeViscoelasticModel in Viscoelasticity_Ting, which I don't like
    !but it keeps each set of code manageable.  

    cur_props%attard_start_nd = NonDimenLength(cur_props%attard_start_dim)
    cur_props%attard_stop_nd = NonDimenLength(cur_props%attard_stop_dim)     
    
    if (cur_props%VEchoice == VEC_THREEELM_E1E2ETA) then
       if (cur_props%input_shear_modulus) then          
          cur_props%Estar           = ShearToReducedModulus(cur_props%threeelm_e1 ,cur_props%Poisson_sample)
          cur_props%threeelm_E2star = ShearToReducedModulus(cur_props%threeelm_E2,cur_props%Poisson_sample)

          !relaxation time is always defined w.r.t. shear modulus, although we may enter either shear or youngs
         cur_props%threeelm_tau_stress = NonDimenTime(cur_props%etasample / (cur_props%threeelm_e1 + cur_props%threeelm_E2))          
       else
          cur_props%Estar           = reducedModulus1(cur_props%threeelm_e1 ,cur_props%Poisson_sample)
          cur_props%threeelm_E2star = reducedModulus1(cur_props%threeelm_E2,cur_props%Poisson_sample)
          G1 = YoungsToShear(cur_props%Esample ,cur_props%Poisson_sample)
          G2 = YoungsToShear(cur_props%threeelm_E2,cur_props%Poisson_sample)
          cur_props%threeelm_tau_stress = NonDimenTime(cur_props%etasample / (G1+G2))
       end if
       
       cur_props%attard_E0        = NonDimenPressure( cur_props%Estar)
       cur_props%attard_Einfinity = NonDimenPressure(( cur_props%Estar * cur_props%threeelm_E2star)/( cur_props%Estar + cur_props%threeelm_E2star))
       
       call assert( cur_props%attard_Einfinity < cur_props%attard_E0, 'Einfinity < E0')
    elseif (cur_props%VEchoice == VEC_THREEELM_E0_EINF_TAU) then
       if (cur_props%input_shear_modulus) then
          cur_props%attard_E0        = NonDimenPressure(ShearToReducedModulus(cur_props%threeelm_E0,cur_props%Poisson_sample))
          cur_props%attard_Einfinity = NonDimenPressure(ShearToReducedModulus(cur_props%threeelm_Einf,cur_props%Poisson_sample))
       else
          cur_props%attard_E0        = NonDimenPressure(reducedModulus1(cur_props%threeelm_E0,  cur_props%Poisson_sample))
          cur_props%attard_Einfinity = NonDimenPressure(reducedModulus1(cur_props%threeelm_Einf,cur_props%Poisson_sample))
       end if

       cur_props%threeelm_tau_stress = NonDimenTime(cur_props%threeelm_tau_stress_dim)
       
    else
       call WriteFatalError("Unhandled constitutive equation for Attard model.  Only the three element model is implemented. please use that")
    end if


    cur_props%tau_creep = cur_props%threeelm_tau_stress * ( cur_props%attard_E0 / cur_props%attard_Einfinity) !time from creep and stress relax not the same!

    
    cur_props%rr = (/ linspace(0d0 , cur_props%attard_radial_extent*Rtip_nd, cur_props%N_attard_spatial ) /)
    !the above includes a point at r=zero.  But that leads to a singularity.  so offset all the points by one-half of the first interval
    cur_props%drr = cur_props%rr(2) - cur_props%rr(1)
    cur_props%rr = cur_props%rr + cur_props%drr / 2d0
    
    do i = 1,cur_props%N_attard_spatial         
       cur_props%tip_shape(i) = (cur_props%rr(i)**2 / (2d0 * Rtip_nd)); 
    end do
    
    
    cur_props%k = precompute_k_s(cur_props%rr, cur_props%rr, cur_props%drr, cur_props%N_attard_spatial, cur_props%bare_k)

    cur_props%kt = transpose(cur_props%k)

    !fourier transform matrix. 
    do i = 1,cur_props%N_attard_spatial
       do j = 1,cur_props%N_attard_fourier
          !bahram calls this cos_r
          cur_props%C(i,j) = cos( cur_props%rr(i) * (j-1) * pi / cur_props%rr(cur_props%N_attard_spatial) )
       end do
    end do

    
    cur_props%CT = transpose(cur_props%C)

    if (cur_props%fts_model == ATTARD_FOURIER_LSQ) then
       cur_props%use_bahram_version = .false.
    else
       cur_props%use_bahram_version = .true.

       allocate( cur_props%gamma( cur_props%N_attard_spatial, cur_props%N_attard_fourier)) !for bahrams's version
       allocate( gamma1( cur_props%N_attard_spatial, cur_props%N_attard_fourier)) !for bahrams's version

       !for bahram's versions
       !    gamma =  matmul(k, C)*drr*drr  !my k has already been mutliplied by radius, don't know why this isn't the same as the long way below

       gamma1 = matmul(cur_props%bare_k, cur_props%C)
       do j = 1, cur_props%N_attard_spatial
          cur_props%gamma(j,:) = cur_props%rr(j) * gamma1(j,:) * cur_props%drr * cur_props%drr
       end do

       deallocate( gamma1)
    end if
    
  end subroutine initialize_attard


  real*8 function Fts_attard_spatial(d,t, u_cur, cur_props)  
    use Integrals
    type(matl_prop), intent(in) :: cur_props
    real*8, intent(in) :: d, u_cur(cur_props%N_attard_spatial),t
    real*8 :: gap(cur_props%N_attard_spatial),  pressure(cur_props%N_attard_spatial)

    u0 = u_cur(1)

    if (.not. attard_enabled) then
       Fts_attard_spatial=0
       return
    end if
    
    gap = compute_gap(cur_props%tip_shape,d, u_cur, cur_props)
        
    pressure = compute_p(gap, cur_props) 

!    call assert( pressure(cur_props%N_attard_spatial) < 0.2 * pressure(1), 'radial extent too small?')


    Fts_attard_spatial =  2d0  * pi * trapzu(pressure * cur_props%rr, cur_props%N_attard_spatial, cur_props%drr )
       
    call assert( DimenForce( Fts_attard_spatial) < 1e3, 'attard Fts > 1e12 nN')    ! 1e12 nN = 1e3 N.  should never get this high
  end function Fts_attard_spatial


  real*8 function Fts_attard_fourier(d,t, a_cur, cur_props)  
    use Integrals
    type(matl_prop), intent(in) :: cur_props
    real*8, intent(in) :: d, a_cur(cur_props%N_attard_fourier),t
    real*8 :: u_cur(cur_props%N_attard_spatial)

    if (.not. attard_enabled) then
       Fts_attard_fourier=0
       return
    end if
    
    u_cur = compute_u_cur( a_cur, cur_props)    
    Fts_attard_fourier = Fts_attard_spatial(d,t, u_cur, cur_props)  
  end function Fts_attard_fourier

  function compute_u_cur( a_cur, cur_props)
    type(matl_prop), intent(in) :: cur_props
    real*8, intent(in) :: a_cur(cur_props%N_attard_fourier)
    real*8 :: compute_u_cur(cur_props%N_attard_spatial)
    integer :: i,j
    !expand spatial points in terms of their fourier components
    compute_u_cur = 0
    do i = 1,cur_props%N_attard_spatial
       do j = 1,cur_props%N_attard_fourier
          compute_u_cur(i) = compute_u_cur(i) + a_cur(j) * cur_props%C(i,j)
       end do
    end do
  end function compute_u_cur

  
  !this computes equation 16 in Rajabifar Macromolecules 2018
  !stores into module variables
  subroutine compute_J_b(d,dpr,t,u_cur, cur_props)
    use Integrals
    real*8, intent(in) :: d, dpr,t
    type(matl_prop), intent(inout) :: cur_props
      real*8 :: u_cur(cur_props%N_attard_spatial), uinfinity(cur_props%N_attard_spatial),  gap(cur_props%N_attard_spatial),  pressure(cur_props%N_attard_spatial), p_prime(cur_props%N_attard_spatial)

      integer :: j
      
      gap = compute_gap(cur_props%tip_shape,d, u_cur, cur_props)

      !$OMP PARALLEL SECTIONS
      !$OMP SECTION
      pressure = compute_p(gap, cur_props) 
      !$OMP SECTION
      p_prime = compute_p_prime(gap, cur_props)  !derivative of pressure
      !$OMP END PARALLEL SECTIONS
      
      ! do j = 1,cur_props%N_attard_spatial
      !       uinfinity(j) = -1d0/(cur_props%attard_Einfinity) * trapzu( pressure * cur_props%k(j,:), cur_props%N_attard_spatial, cur_props%drr) !vector multiply
      ! end do
      
      ! cur_props%b(:,1)  = (1/cur_props%tau_creep) * ( u_cur - uinfinity);
      
      ! do j = 1,cur_props%N_attard_spatial
      !    cur_props%J(j,:) = (1d0/cur_props%attard_E0) * cur_props%drr * p_prime * cur_props%k(j,:)

      !    !trapezoidal integration.  factor of 1/2 for first and last points, 1 for interior
      !    cur_props%J(j,1)                          = cur_props%J(j,1)/2d0
      !    cur_props%J(j,cur_props%N_attard_spatial) = cur_props%J(j,cur_props%N_attard_spatial)/2d0

         
      !    cur_props%b(j,1) = cur_props%b(j,1) + sum( cur_props%J(j,:) ) * dpr  !sum here, but we've already applied factors of 1/2 to J above, so it can effectively be trapezoidal
         
      !    cur_props%J(j,j) = cur_props%J(j,j) - 1d0
      ! end do

!this is a decent speedup.  maybe 40% overall faster with 3 threads
!$OMP PARALLEL DO
      do j = 1,cur_props%N_attard_spatial
         uinfinity(j) = -1d0/(cur_props%attard_Einfinity) * trapzu( pressure * cur_props%k(j,:), cur_props%N_attard_spatial, cur_props%drr) !vector multiply
                      
         cur_props%J(j,:) = (1d0/cur_props%attard_E0) * cur_props%drr * p_prime * cur_props%k(j,:)

         !trapezoidal integration.  factor of 1/2 for first and last points, 1 for interior
         cur_props%J(j,1)                          = cur_props%J(j,1)/2d0
         cur_props%J(j,cur_props%N_attard_spatial) = cur_props%J(j,cur_props%N_attard_spatial)/2d0
         
         cur_props%b(j,1) =  sum( cur_props%J(j,:) ) * dpr  !sum here, but we've already applied factors of 1/2 to J above, so it can effectively be trapezoidal
                  
         cur_props%J(j,j) = cur_props%J(j,j) - 1d0
      end do
!$OMP END PARALLEL DO
      
      cur_props%b(:,1)  = cur_props%b(:,1) +  (1/cur_props%tau_creep) * ( u_cur - uinfinity);
      
    end subroutine compute_J_b



!this computes equation 16 in Rajabifar Macromolecules 2018
!stores into module variables
    !the transposed way of doing it is faster (hits contiguous memory regions).
  subroutine compute_J_b_transposed(d,dpr,t,u_cur, cur_props)
    use Integrals
    real*8, intent(in) :: d, dpr,t
    type(matl_prop), intent(inout) :: cur_props
      real*8 :: u_cur(cur_props%N_attard_spatial), uinfinity(cur_props%N_attard_spatial),  gap(cur_props%N_attard_spatial),  pressure(cur_props%N_attard_spatial), p_prime(cur_props%N_attard_spatial)

      integer :: j
      
      gap = compute_gap(cur_props%tip_shape,d, u_cur, cur_props)

      !$OMP PARALLEL SECTIONS
      !$OMP SECTION
      pressure = compute_p(gap, cur_props) 
      !$OMP SECTION
      p_prime = compute_p_prime(gap, cur_props)  !derivative of pressure
      !$OMP END PARALLEL SECTIONS
      
!this is a decent speedup.  maybe 40% overall faster with 3 threads
!$OMP PARALLEL DO
      do j = 1,cur_props%N_attard_spatial
         uinfinity(j) = -1d0/(cur_props%attard_Einfinity) * trapzu( pressure * cur_props%kt(:,j), cur_props%N_attard_spatial, cur_props%drr) !vector multiply
                      
         cur_props%J(:,j) = (1d0/cur_props%attard_E0) * cur_props%drr * p_prime * cur_props%kt(:,j)

         !trapezoidal integration.  factor of 1/2 for first and last points, 1 for interior
         cur_props%J(1,j)                          = cur_props%J(1,j)/2d0
         cur_props%J(cur_props%N_attard_spatial,j) = cur_props%J(cur_props%N_attard_spatial,j)/2d0
         
         cur_props%b(j,1) =  sum( cur_props%J(:,j)) * dpr  !sum here, but we've already applied factors of 1/2 to J above, so it can effectively be trapezoidal
                  
         cur_props%J(j,j) = cur_props%J(j,j) - 1d0
      end do
!$OMP END PARALLEL DO
      
      cur_props%b(:,1)  = cur_props%b(:,1) +  (1/cur_props%tau_creep) * ( u_cur - uinfinity);
      
    end subroutine compute_J_b_transposed


    

    function compute_attard_derivative(d, dpr,t, u_cur, cur_props)
      real*8, intent(in) :: d, dpr, u_cur(:),t
      type(matl_prop), intent(inout) :: cur_props
      real*8, dimension( size(u_cur)) :: compute_attard_derivative
      
      if (.not. attard_enabled) then
         compute_attard_derivative=0
         return
      end if
    
      
      if (cur_props%N_attard_fourier == 0) then
         !transposed version is faster because it hits contiguous memory locations.        
!         compute_attard_derivative = compute_new_udot(d, dpr, u_cur, cur_props)
         compute_attard_derivative = compute_new_udot_transposed(d, dpr,t, u_cur, cur_props)
      else
         if (cur_props%use_bahram_version) then
            compute_attard_derivative = compute_new_adot_bahram(d, dpr,t, u_cur, cur_props)
         else
!            compute_attard_derivative = compute_new_adot_lsq(d, dpr, u_cur, cur_props)
            compute_attard_derivative = compute_new_adot_lsq_transposed(d, dpr,t, u_cur, cur_props)
         end if
      end if

    end function compute_attard_derivative





    function compute_new_udot_transposed(d, dpr,t, u_cur, cur_props)
      type(matl_prop), intent(inout) :: cur_props
#if defined(HAVE_LAPACK95MKL)
    use lapack95
#elif defined( HAVE_LAPACK)
    integer :: iwork(cur_props%N_attard_spatial)
    real*8 :: work(4*cur_props%N_attard_spatial)
    real*8, external :: dlange
#else
    use LinearAlgebra
#endif
    real*8 :: compute_new_udot_transposed(cur_props%N_attard_spatial), anorm, rcond
    real*8, intent(in) :: d, dpr, u_cur(cur_props%N_attard_spatial),t
    integer :: info
    logical ok
    character*16 :: istr
    
    call compute_J_b_transposed(d,dpr,t,u_cur, cur_props)

    
#if defined( HAVE_LAPACK95MKL)
    call WriteFatalError('not implemented')
#elif defined( HAVE_LAPACK)
     call dgetrf( cur_props%N_attard_spatial, cur_props%N_attard_spatial,  cur_props%J, cur_props%N_attard_spatial, cur_props%ipiv, info)
     if (info .ne. 0) then
        write(istr,*) info
        call WriteFatalError("error computing LU factorization. info = " // trim(istr) // " this may indicate a bug in the code, or a numerical issue.")
     end if
     
     anorm = dlange( '1',  cur_props%N_attard_spatial, cur_props%N_attard_spatial, cur_props%J, cur_props%N_attard_spatial, work)      
     call dgecon('1',  cur_props%N_attard_spatial,  cur_props%J,  cur_props%N_attard_spatial, anorm, rcond, work, iwork, info)
     if (rcond < 1e-10) write(*,*) "warning poorly conditioned matrix: rcond = ", rcond, "(spatial)"
     if (rcond < 1e-13) call WriteFatalError('Matrix is poorly conditioned. This may indicate a bug in the program. Or you may have input bad parameters. Trying increase # of points per cycle and/or # of spatial discretization points. This can also be caused by very high attractive forces on compliant samples. Increasing LJ equilibrium distance or decreasing Hamaker constant may help')
     
     call dgetrs( 'T', cur_props%N_attard_spatial, 1,  cur_props%J, cur_props%N_attard_spatial, cur_props%ipiv, cur_props%b,  cur_props%N_attard_spatial, info)
     if (info .ne. 0) call WriteFatalError("error computing solving matrix. this may indicate a bug in the code, or a numerical issue.")
#else
     call WriteFatalError('not implemented')     
#endif

     compute_new_udot_transposed =  cur_props%b(:,1) !this copy is not strictly necessary, but makes code easier to understand.

   end function compute_new_udot_transposed




    
    function compute_new_udot(d, dpr,t, u_cur, cur_props)
      type(matl_prop), intent(inout) :: cur_props
#if defined(HAVE_LAPACK95MKL)
    use lapack95
#elif defined( HAVE_LAPACK)
    integer :: iwork(cur_props%N_attard_spatial)
    real*8 :: work(4*cur_props%N_attard_spatial)
    real*8, external :: dlange
#else
    use LinearAlgebra
#endif
    real*8 :: compute_new_udot(cur_props%N_attard_spatial), anorm, rcond
    real*8, intent(in) :: d, dpr, u_cur(cur_props%N_attard_spatial),t
    integer :: info
    logical ok
    character*16 :: istr
    
    call compute_J_b(d,dpr,t,u_cur, cur_props)

    
#if defined( HAVE_LAPACK95MKL)
     !this is the LAPACK routine (intel math kernel library implementation)
     !three times faster than the generic routine below (at the expense of dragging in an 
     !extra library.  beware the math kernel libraries are not profiled.
     !these routines can also made multi-threaded on multi-core CPUs, leading to even more performance gains
     call dgetrf_f95( cur_props%J ,cur_props%ipiv)
     call dgetrs_f95( cur_props%J, cur_props%ipiv, cur_props%b)
#elif defined( HAVE_LAPACK)
     call dgetrf( cur_props%N_attard_spatial, cur_props%N_attard_spatial,  cur_props%J, cur_props%N_attard_spatial, cur_props%ipiv, info)
     if (info .ne. 0) then
        write(istr,*) info
        call WriteFatalError("error computing LU factorization. info = " // trim(istr) // " this may indicate a bug in the code, or a numerical issue.")
     end if
     
     anorm = dlange( '1',  cur_props%N_attard_spatial, cur_props%N_attard_spatial, cur_props%J, cur_props%N_attard_spatial, work)      
     call dgecon('1',  cur_props%N_attard_spatial,  cur_props%J,  cur_props%N_attard_spatial, anorm, rcond, work, iwork, info)
     if (rcond < 1e-10) write(*,*) "warning poorly conditioned matrix: rcond = ", rcond, "(spatial)"
     if (rcond < 1e-13) call WriteFatalError('Matrix is poorly conditioned. This may indicate a bug in the program. Or you may have input bad parameters. Trying increase # of points per cycle and/or # of spatial discretization points. This can also be caused by very high attractive forces on compliant samples. Increasing LJ equilibrium distance or decreasing Hamaker constant may help')
     
     call dgetrs( 'N', cur_props%N_attard_spatial, 1,  cur_props%J, cur_props%N_attard_spatial, cur_props%ipiv, cur_props%b,  cur_props%N_attard_spatial, info)
     if (info .ne. 0) call WriteFatalError("error computing solving matrix. this may indicate a bug in the code, or a numerical issue.")
#else
     !this is the generic routine we used in my ME614 class
     !can consume upwards of 90% of the program time.
     call solve_gauss_elim(cur_props%J , cur_props%b(:,1), ok)  !in place solve.  vector b is input, but also the output x.
#endif

     compute_new_udot =  cur_props%b(:,1) !this copy is not strictly necessary, but makes code easier to understand.

   end function compute_new_udot





   
 !    !this is bahrams's version   
 !    !u here is still the deflection
 !    !and a is fourier coeff
   function compute_new_adot_bahram(d, dpr,t, a_cur, cur_props)
     type(matl_prop), intent(in) :: cur_props
#if defined( HAVE_LAPACK95MKL)
     use lapack95
#elif defined( HAVE_LAPACK)
     integer :: iwork(cur_props%N_attard_fourier)
     real*8 :: work(4*cur_props%N_attard_fourier)
     real*8, external :: dlange
#else
     use LinearAlgebra
#endif
     real*8 :: compute_new_adot_bahram(cur_props%N_attard_fourier)

     real*8, intent(in) :: d, dpr, a_cur(cur_props%N_attard_fourier), t
     
     !fixme might be faster to pre-allocate all of these arrays once
     real*8 :: u_cur(cur_props%N_attard_spatial)
     real*8 :: uinfinity(cur_props%N_attard_spatial),  gap(cur_props%N_attard_spatial),  pressure(cur_props%N_attard_spatial), p_prime(cur_props%N_attard_spatial)
     real*8 :: A(cur_props%N_attard_fourier, cur_props%N_attard_fourier), AA(cur_props%N_attard_fourier), anorm, rcond
     real*8 :: Coef1(cur_props%N_attard_fourier), Coef2(cur_props%N_attard_fourier, cur_props%N_attard_fourier), temp(cur_props%N_attard_spatial, cur_props%N_attard_fourier), BB(cur_props%N_attard_fourier)
     
     logical ok
     character*20 :: istr

      integer :: j, info

      u_cur = compute_u_cur( a_cur, cur_props)

      gap = compute_gap(cur_props%tip_shape,d, u_cur, cur_props)
      
      pressure = compute_p(gap, cur_props) 
      
      p_prime = compute_p_prime(gap, cur_props)  !derivative of pressure
      
    
      Coef1   = cur_props%rr(cur_props%N_attard_spatial) *cur_props%attard_E0/2d0 !whole vector operation
      Coef1(1)= cur_props%rr(cur_props%N_attard_spatial) *cur_props%attard_E0
      Coef2 = 0 !whole array
      do j = 1,cur_props%N_attard_fourier
         Coef2(j,j)=Coef1(j)  !Coef2 = diag(Coef1)
      end do
      
      do j = 1,cur_props%N_attard_fourier
         temp(j,:)=  cur_props%C(j,:) * p_prime(j)
      end do
      A = matmul( transpose(cur_props%gamma), temp)
      
      BB= matmul(transpose(cur_props%gamma),pressure) !fixme, speedup possible here... gamma is only every used in transpose form, so transpose it once at initialization

      AA=(a_cur*Coef1/cur_props%tau_creep)+(cur_props%attard_E0/cur_props%tau_creep/cur_props%attard_Einfinity*BB)+A(:,1)*dpr
        
      A=A-Coef2

           
#if defined( HAVE_LAPACK95MKL)
     call dgetrf_f95( A ,cur_props%ipiv)
     call dgetrs_f95( A, cur_props%ipiv, AA)
#elif defined( HAVE_LAPACK)
     call dgetrf( cur_props%N_attard_fourier, cur_props%N_attard_fourier, A, cur_props%N_attard_fourier, cur_props%ipiv, info)

     !fixme, put some error checks on this
     anorm = dlange( '1',  cur_props%N_attard_fourier, cur_props%N_attard_fourier, A,  cur_props%N_attard_fourier, work)          
     call dgecon('1',  cur_props%N_attard_fourier,  A,  cur_props%N_attard_fourier, anorm, rcond, work, iwork, info)
     if (rcond < 1e-10) write(*,*) "warning poorly conditioned matrix: rcond = ", rcond, "(Bahram fourier)"
     if (rcond < 1e-13) call WriteFatalError('Matrix is poorly conditioned. This may indicate a bug in the program. Or you may have input bad parameters. Trying increase # of points per cycle and/or # of spatial discretization points. This can also be caused by very high attractive forces on compliant samples. Increasing LJ equilibrium distance or decreasing Hamaker constant may help')
     
     if (info .ne. 0) then
        write(istr,*) info
        call WriteFatalError("error computing LU factorization. info = " // trim(istr) // " this may indicate a bug in the code, or a numerical issue.")
     end if
     
     call dgetrs( 'N', cur_props%N_attard_fourier, 1, A , cur_props%N_attard_fourier, cur_props%ipiv, AA,  cur_props%N_attard_fourier, info)
     if (info .ne. 0) call WriteFatalError("error solving matrix.  this may indicate a bug in the code, or a numerical issue.")
#else
     call solve_gauss_elim(A , AA, ok)  !in place solve. 
#endif

     compute_new_adot_bahram =  AA !this copy is not strictly necessary, but makes code easier to understand.      
     
   end function compute_new_adot_bahram


    
  ! this is more of a least-squares fit.  matches the spatial results and Bahram's fourier results (when well conditioned).
  ! Slighty slower than Bahram's fourier results, but still faster than the spatial version.
  !   
  !u here is still the deflection
  !and a is fourier coeff
   function compute_new_adot_lsq(d, dpr,t, a_cur, cur_props)
     type(matl_prop), intent(inout) :: cur_props
#if defined( HAVE_LAPACK95MKL)
    use lapack95
#elif defined( HAVE_LAPACK)
    integer :: iwork(cur_props%N_attard_fourier)
    real*8 :: work(4*cur_props%N_attard_fourier)
    real*8, external :: dlange
#else
    use LinearAlgebra
#endif
    real*8 :: compute_new_adot_lsq(cur_props%N_attard_fourier)
    real*8, intent(in) :: d, dpr, a_cur(cur_props%N_attard_fourier),t
    
!fixme might be faster to pre-allocate all of these arrays once
    real*8 :: u_cur(cur_props%N_attard_spatial)
    real*8 :: JC(cur_props%N_attard_fourier, cur_props%N_attard_fourier), bC(cur_props%N_attard_fourier, 1), anorm, rcond
    
    integer :: info
    logical ok
    character*20 :: istr

    u_cur = compute_u_cur( a_cur, cur_props)
    
    call compute_J_b(d,dpr,t,u_cur, cur_props)
    
    !then here do matrix-matrix multiply of the spatial j times the fourier transform matrix
    if (cur_props%N_attard_spatial == cur_props%N_attard_fourier) then
       !square matrices
       JC = matmul( cur_props%J, cur_props%C)
       bC = cur_props%b
    else
       !this is what I originally had. can't rigorously justify it, although seems to match Bahram's code better than the below
       !both approaches have very poor matrix condition numbers
       JC =  matmul( cur_props%CT ,  matmul( cur_props%J, cur_props%C) )
       bC =  matmul( cur_props%CT, cur_props%b)

       !this is a correct least squares formulation. can rigorously justify it.  slightly different, but still doesnt match matlab.
!       JC =  matmul( matmul(CT, transpose(J)) ,  matmul( J, C) )      
!       bC =  matmul( CT , matmul( transpose(J), b))
    end if

    
#if defined( HAVE_LAPACK95MKL)
     !this is the LAPACK routine (intel math kernel library implementation)
     !three times faster than the generic routine below (at the expense of dragging in an 
     !extra library.  beware the math kernel libraries are not profiled.
     !these routines can also made multi-threaded on multi-core CPUs, leading to even more performance gains
     call dgetrf_f95( JC ,cur_props%ipiv)
     call dgetrs_f95( JC, cur_props%ipiv, bC)
#elif defined( HAVE_LAPACK)
     call dgetrf( cur_props%N_attard_fourier, cur_props%N_attard_fourier, JC, cur_props%N_attard_fourier, cur_props%ipiv, info)

     anorm = dlange( '1',  cur_props%N_attard_fourier, cur_props%N_attard_fourier, JC ,  cur_props%N_attard_fourier, work)          
     call dgecon('1',  cur_props%N_attard_fourier,  JC ,  cur_props%N_attard_fourier, anorm, rcond, work, iwork, info)
!     write(*,*) rcond
     if (rcond < 1e-10) write(*,*) "warning poorly conditioned matrix: rcond = ", rcond, "(daniel fourier)"
     if (rcond < 1e-13) call WriteFatalError('Matrix is poorly conditioned.  This may indicate a bug in the program.  Or you may have input a bad set of parameters.  Trying increase the number of points per cycle, the number of spatial discretization points, and/or the number of fourier basis functions')
     
     if (info .ne. 0) then
        write(istr,*) info
        call WriteFatalError("error computing LU factorization. info = " // trim(istr) // " this may indicate a bug in the code, or a numerical issue.")
     end if
     
     call dgetrs( 'N', cur_props%N_attard_fourier, 1,  JC, cur_props%N_attard_fourier, cur_props%ipiv, bC,  cur_props%N_attard_fourier, info)
     if (info .ne. 0) call WriteFatalError("error solving matrix.  this may indicate a bug in the code, or a numerical issue.")

#else
     !this is the generic routine we used in my ME614 class
     !can consume upwards of 90% of the program time.
     call solve_gauss_elim(JC , bC(:,1), ok)  !in place solve.  vector b is input, but also the output x.
#endif

     compute_new_adot_lsq =  bC(:,1) !this copy is not strictly necessary, but makes code easier to understand.

   end function compute_new_adot_lsq








 ! this is more of a least-squares fit.  matches the spatial results and Bahram's fourier results (when well conditioned).
  ! Slighty slower than Bahram's fourier results, but still faster than the spatial version.
  !   
  !u here is still the deflection
  !and a is fourier coeff
   function compute_new_adot_lsq_transposed(d, dpr,t, a_cur, cur_props)
     type(matl_prop), intent(inout) :: cur_props
#if defined( HAVE_LAPACK95MKL)
    use lapack95
#elif defined( HAVE_LAPACK)
    integer :: iwork(cur_props%N_attard_fourier)
    real*8 :: work(4*cur_props%N_attard_fourier)
    real*8, external :: dlange
#else
    use LinearAlgebra
#endif
    real*8 :: compute_new_adot_lsq_transposed(cur_props%N_attard_fourier)
    real*8, intent(in) :: d, dpr, a_cur(cur_props%N_attard_fourier),t

!fixme might be faster to pre-allocate all of these arrays once
    real*8 :: u_cur(cur_props%N_attard_spatial)
    real*8 :: JC(cur_props%N_attard_fourier, cur_props%N_attard_fourier), bC(cur_props%N_attard_fourier, 1), anorm, rcond
    
    integer :: info
    logical ok
    character*20 :: istr

    u_cur = compute_u_cur( a_cur, cur_props)
    
    call compute_J_b_transposed(d,dpr,t,u_cur, cur_props)
    
    !then here do matrix-matrix multiply of the spatial j times the fourier transform matrix
    if (cur_props%N_attard_spatial == cur_props%N_attard_fourier) then
       !square matrices
       JC = matmul( cur_props%J, cur_props%CT)
       bC = cur_props%b
    else
       !this is what I originally had. can't rigorously justify it, although seems to match Bahram's code better than the below
       !both approaches have very poor matrix condition numbers
       JC =  matmul( cur_props%CT ,  matmul( cur_props%J, cur_props%C) )
       bC =  matmul( cur_props%CT, cur_props%b)

       !this is a correct least squares formulation. can rigorously justify it.  slightly different, but still doesnt match matlab.
!       JC =  matmul( matmul(CT, transpose(J)) ,  matmul( J, C) )      
!       bC =  matmul( CT , matmul( transpose(J), b))
    end if

    
#if defined( HAVE_LAPACK95MKL)
    call WriteFatalError('not implemented')
#elif defined( HAVE_LAPACK)
     call dgetrf( cur_props%N_attard_fourier, cur_props%N_attard_fourier, JC, cur_props%N_attard_fourier, cur_props%ipiv, info)

     anorm = dlange( '1',  cur_props%N_attard_fourier, cur_props%N_attard_fourier, JC ,  cur_props%N_attard_fourier, work)          
     call dgecon('1',  cur_props%N_attard_fourier,  JC ,  cur_props%N_attard_fourier, anorm, rcond, work, iwork, info)
!     write(*,*) rcond
     if (rcond < 1e-10) write(*,*) "warning poorly conditioned matrix: rcond = ", rcond, "(daniel fourier)"
     if (rcond < 1e-13) call WriteFatalError('Matrix is poorly conditioned.  This may indicate a bug in the program.  Or you may have input a bad set of parameters.  Trying increase the number of points per cycle, the number of spatial discretization points, and/or the number of fourier basis functions')
     
     if (info .ne. 0) then
        write(istr,*) info
        call WriteFatalError("error computing LU factorization. info = " // trim(istr) // " this may indicate a bug in the code, or a numerical issue.")
     end if
     
     call dgetrs( 'T', cur_props%N_attard_fourier, 1,  JC, cur_props%N_attard_fourier, cur_props%ipiv, bC,  cur_props%N_attard_fourier, info)
     if (info .ne. 0) call WriteFatalError("error solving matrix.  this may indicate a bug in the code, or a numerical issue.")

#else
    call WriteFatalError('not implemented')
#endif

     compute_new_adot_lsq_transposed =  bC(:,1) !this copy is not strictly necessary, but makes code easier to understand.

   end function compute_new_adot_lsq_transposed




   
 function compute_p(gap, cur_props) !pressure
   use params
   type(matl_prop), intent(in) :: cur_props
   real*8, intent(in) :: gap(cur_props%N_attard_spatial)
   real*8 :: compute_p(cur_props%N_attard_spatial)

   compute_p = ( cur_props%A_hamaker_nd  /(6.d0 * pi * gap**3)) * ( (cur_props%LJ_r0_nd**6)      /(gap**6) - 1.d0)   
 end function compute_p

 function compute_p_prime(gap, cur_props)
   use params
   type(matl_prop), intent(in) :: cur_props
   real*8, intent(in) :: gap(cur_props%N_attard_spatial)
   real*8 :: compute_p_prime(cur_props%N_attard_spatial)

   compute_p_prime = ( 0.5d0 * cur_props%A_hamaker_nd *       (-3.d0 * cur_props%LJ_r0_nd**6 + gap**6)       / (pi * gap**10))

 end function compute_p_prime

 function compute_gap(tip_shape,d, u_cur, cur_props)
   type(matl_prop), intent(in) :: cur_props
   real*8, intent(in) :: tip_shape(cur_props%N_attard_spatial), d, u_cur(cur_props%N_attard_spatial)
   real*8 :: compute_gap(cur_props%N_attard_spatial)
   integer :: i
   
   do i = 1,cur_props%N_attard_spatial      
      compute_gap(i) = tip_shape(i) - u_cur(i) + d
   end do
 end function compute_gap
 

 ! this precomputes the kernel k, multiplied by s (because the kernel
 ! is only ever found multiplied by the radius).  for backwards compatibility
 !to bahrams' code, we also output just the kernal by itself
  function precompute_k_s(r, s, drr, N, bare_k)
    
  use ellipticIntegrals

  integer, intent(in) :: N
  real*8, intent(in) :: r(N), s(N), drr
  real*8, dimension(N,N) :: precompute_k_s
  real*8, dimension(N,N), intent(out) :: bare_k  !for comparison to bahram

  integer :: i,j
  real*8 ::  Inp, Inm

  do i = 1, N   
     do j = 1,N

        if (s(i) < r(j)) then
           !note the different definition of elliptic integral than matlab.  results in us doing sqrt of the argument
           !(pretty sure that's right).  
           precompute_k_s(i,j) = 4d0/(pi*r(j)) * elliptic_k(s(i) / r(j))
        elseif (s(i) > r(j)) then
           precompute_k_s(i,j) = 4d0/(pi*s(i)) * elliptic_k(r(j) / s(i))
        elseif (r(j) > 0) then            
           Inm = (drr/pi - drr**2/(4d0*pi*r(j))) * (1d0 + log((16d0*r(j)**2)/(r(j) * drr - drr**2/4d0)))
           Inp = (drr/pi) * log( (16d0*(r(j)+drr/2d0)**2) / (r(j) * drr + drr**2/4d0)) + 4d0 *r(j)/pi * log( (2d0* r(j)+drr)/(2d0 * r(j)+drr/2d0))    
           precompute_k_s(i,j) = (1d0 / ( r(j) * drr))  * (Inp + Inm)
        else
           precompute_k_s(i,j) = 0 !  r=s=0. this should not ever happen actually
        end if
       end do
    end do

    bare_k = precompute_k_s
        
    do i = 1,N
       precompute_k_s(i,:) = precompute_k_s(i,:) * s
    end do
    
  end function precompute_k_s

end module Viscoelasticity_attard


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


module CustomFts_MD
  use matl_prop_type_module
  use NonDimensionalization
  integer, parameter :: APP = 1, RET = 2, D_NONE = 0
  integer MD_state
  real*8 max_indent, md_reverse_ok, md_contact_start
  integer, parameter :: MD_DANIEL_1 = 1, MD_GABRIELA_1 = 2
  integer :: md_equation_type
  real*8, private, allocatable :: coeff(:)

  contains

    subroutine Initialize_Custom_MD
      MD_state = D_NONE
      !which fit function do we want to use?
      md_equation_type = MD_GABRIELA_1

      if (md_equation_type == MD_DANIEL_1) then
         md_reverse_ok = 0
         md_contact_start = 0
         if (.not.(allocated(coeff))) then
            allocate(coeff(5))
         end if
      else
         md_reverse_ok = 0.0
         if (.not.(allocated(coeff))) then
            allocate(coeff(8))
         end if
         !these could be read in from disk for more general situations
         coeff(4) = 0.376e-9 !dca
         coeff(5) = 0.0153e-18 !Ha
         coeff(6) = 2.5271e9 !Ea
         coeff(7) = 0.3e-9 !alpha
         coeff(8) = 6e-9 !Rtip  (12 was a diameter!!!)
         md_contact_start = NonDimenLength(4 * coeff(7) )
      end if

    end subroutine Initialize_Custom_MD

  !this function interpolates from molecular dynamic data
    real*8 function Fts_custom_md( d, mp)
      type(matl_prop), intent(in) :: mp
      real*8, intent(in) :: d

      if ( d > md_contact_start ) then
         Fts_custom_md = 0
      else
         if ( DimenLength(d) <  mp%indentation_x(1)) then
            call WriteFatalError( "indentation is below maximum MD data.  Cannot extrapolate.")
         end if

         if ( MD_state == APP) then
            Fts_custom_md = Fts_custom_md_app(d, mp)
         elseif (MD_state == RET) then
            Fts_custom_md = Fts_custom_md_ret(d, mp)
         else
            Fts_custom_md = 0d0
         end if
      end if
       
    end function Fts_custom_md

    real*8 function Fts_custom_md_app(d, mp)
      real*8, intent(in) :: d
      type(matl_prop), intent(in) :: mp
      real*8 b, dd
      integer :: i

      !right now matlab file coefficients are setup for inputs in nanometers, with
      !indentation negative.
      !and outputs are in nN

      select case(md_equation_type)
      case (MD_DANIEL_1)
         dd = DimenLength(d) * 1e9
         !polynomial.  horner's method is quicker (less multiplications)
         b = mp%app_coeff( 1 )
         do i = 2, mp%num_app_coeff         
            b = mp%app_coeff(i) + b * dd
         end do
         Fts_custom_md_app = NonDimenForce(b) * 1e-9
      case (MD_GABRIELA_1)
         dd = DimenLength(d)
         !attarctive part
         ! -Ha * Rtip / dshield(x-dca)^2
         b = -coeff(5) * coeff(8) / ( d_shield(dd-coeff(4), coeff(7))**2 )

         !repulsive part
         if ((dd-coeff(4)) < 0) b = b + (4./3.) * sqrt(coeff(8)) * coeff(6) * (coeff(4)-dd)**1.5         
         Fts_custom_md_app = NonDimenForce(b)
      end select


    end function Fts_custom_md_app

    real*8 function Fts_custom_md_ret(d, mp)
      type(matl_prop), intent(in) :: mp
      real*8, intent(in) :: d
      real*8 ddim
     integer :: i

      !to be completely generaly in handling the coefficients and potential functions, keep it all dimensional in here.  also currently coefficients are fit in nm


      select case(md_equation_type)
         case (MD_DANIEL_1)
            ddim = DimenLength( d ) * 1e9
            !this was foo opt 2 in 11-6-30-mdfit.  presented to strachan.  just a placeholder. will do something better 
            !fun =  x(1) * (d - x(2)).^2 .* (d <= x(5))  -  x(1) * (x(5) - x(2)).^2 .* (d <= x(5)) +  x(3) * exp( x(4) * x(5)) .* (d < x(5)) + x(3) * exp( x(4) * d) .* (d >= x(5)) ;      
            if (ddim < coeff(5)) then
               Fts_custom_md_ret = coeff(1) * (ddim- coeff(2))**2  -  coeff(1) * (coeff(5) - coeff(2))**2 + coeff(3) * exp( coeff(4) * coeff(5))
            elseif (ddim > coeff(5)) then
               Fts_custom_md_ret = coeff(3) * exp( coeff(4) * ddim)
            else ! ( ddim == coeff(5) ) 
               Fts_custom_md_ret = coeff(1) * (ddim- coeff(2))**2  -  coeff(1) * (coeff(5) - coeff(2))**2 + coeff(3) * exp( coeff(4) * ddim)
            end if

            Fts_custom_md_ret = NonDimenForce( Fts_custom_md_ret) * 1e-9

            case (MD_GABRIELA_1)                
               ddim = DimenLength( d )
               !attarctive part
               ! -Hr * Rtip / dshield(x-dcr)^2
               Fts_custom_md_ret = -coeff(3) * coeff(8) / d_shield(ddim-coeff(2), coeff(7))**2
               
               !repulsive part

               if  ((ddim-coeff(2)) < 0) then
                  Fts_custom_md_ret = Fts_custom_md_ret + (4./3.) * coeff(1) * (coeff(2)-ddim)**2
               end if
               Fts_custom_md_ret = NonDimenForce( Fts_custom_md_ret)
            end select


    end function Fts_custom_md_ret

    !aux function for gabriela
    pure real*8 function d_shield(r, alpha)
      real*8, intent(in) :: r, alpha
      if (r > 0) then
         d_shield = (r**3 + alpha**3) ** (1./3.)
      else
         d_shield = alpha
      end if
    end function d_shield

  subroutine UpdateCustomMD_OutsideRES1(d, dpr, mp)
    use Interpolation
    real*8, intent(in) :: d, dpr
    type(matl_prop), intent(in) :: mp
    integer :: i

    if (d > md_contact_start) then
       !not in contact.
       MD_state = D_NONE
       max_indent = 0d0
    else
       !d<=0 in contact
       if ( dpr < 0)  then
          !gap decreasing => indentation increasing.  approaching
          if (MD_state == D_NONE) then
             MD_state = APP        
             !this is the first point in contact.
          elseif (MD_state == RET) then
             call WriteFatalError('Indentation reversal with MD data is not implemented.');
          end if
       else
          if (MD_state == APP) then
             !this is the first decreasing point.
             MD_state = RET        
             max_indent = d

             select case (md_equation_type)
                case (MD_DANIEL_1)
                   do i = 1, mp%num_ret_coeff
                      !indentation_x is dimensional in meters.
                      !max_indent is non-dimensional
                      coeff(i) = interpolate( mp%indentation_x, mp%ret_coeff(i,:), DimenLength(max_indent), mp%num_x)
                   end do
                case (MD_GABRIELA_1)
                   if (max_indent > mp%indentation_x(mp%num_x)) then
                      !have not indented enough to have any reversal data.
                      !for now, just hack it back to approach
                      MD_state = APP
                   else
                      coeff(1) = interpolate( mp%indentation_x, mp%ret_coeff(1,:), DimenLength(max_indent), mp%num_x)
                      coeff(2) =  0.5789*DimenLength(max_indent) + 0.0609e-9 + coeff(4)  ! dcr [m]
                      coeff(3) = -0.0725e-9*DimenLength(max_indent) + 0.0184e-18        ! Hr [N*m]
                   end if
                end select
          end if
       end if
    end if ! d>=0

  end subroutine UpdateCustomMD_OutsideRES1

end module CustomFts_MD



module contactModels
      use params
      use PreisachModule
      use Viscoelasticity_ting
      use Viscoelasticity_attard
      use Nondimensionalization
      use matl_prop_type_module
      use CustomFts_MD
      use HertzModule !avoid circular dependencies
      implicit none
      save
      
      type(matl_prop), target :: substrate_props, feature_props
      type(matl_prop), pointer :: cur_props
      type(forceCoeff), target :: subs_fC, ft_fC
      type(forceCoeff), pointer :: cur_fC

      !electrostatics, gil
      real*8 :: theta_lever, theta_tip, es_c3, es_c4, es_delta, es_h, es_l, es_w,  es_h_dim, es_l_dim, es_w_dim

      real*8, private :: t1  !contact time calculations 
      logical, private :: check

      integer, private :: SurfHystMode, SurfHystMode_prev
      real*8, private, parameter :: SurfHyst_decay_dim = 0.1e-9
      real*8, private :: SurfHyst_decay_nd, SurfHyst_Y, SurfHyst_Y0, SurfHyst_Y_prev, SurfHyst_Y0_prev
      real*8 , private :: SurfHyst_d0, SurfHyst_d0_prev !the gap at which the most recent mode switch was made

      integer, private :: CapAdMode, CapAd_prev, WLC_mode, WLC_prev
      integer, private :: mJKR, mJKR_prev  ! These variables mark whether the JKR
      			       ! interaction force is active or if the surface forces have not yet
      			       ! begun.  The JKR force initiates when d = aDMT, then it breaks at some
      			       ! critical gap for which dJKRcrit > aDMT. 
  

      private Fts_CapAd, Fts_nc_vel_dep, Fts_cons, Fts_oscillatory, Fts_hydration, Fts_DMT

      contains

        !see dev notes for formulas. it's the cubic formula
        real*8 function WLC_ruppture_length(F_in, Lp, L0, temp)
          real*8, intent(in) :: F_in, Lp, L0, temp
          real*8 :: Fs, a, b,c,d, del
          complex*8 :: t, dL0, Q, dL1, dL2, QQ

          Fs = -F_in * Lp / (kB * temp)

          a = 1d0
          b = Fs - 9d0/4d0
          c = -2d0 * Fs + 1.5d0
          d = 1d0

           del = 18d0 * a * b * c * d - 4d0 * b**3 * d + b**2 * c**2 - 4d0 * a * c**3 - 27d0 * a**2 * d**2

          Q = 192.d0 * Fs**3 + 432.d0 * Fs**2 + 324.d0 * Fs + 405.d0
          t =  ( -144.d0 * Fs**2 - 108.d0 * Fs - 243.d0 - 64.d0 * Fs**3 + 12.d0 * sqrt(Q))
          dL0 = (1.d0/12.d0) * t**(1.d0/3.d0) - 12.d0 * (-Fs / 6.d0 - 1.d0/16.d0 - Fs**2 / 9.d0) * t**(-1.d0/3.d0) - Fs / 3.d0 + 0.d075

          dL1 = - (1.d0/24.d0) * t**(1.d0/3.d0) +  6.d0 * (-Fs / 6.d0 - 1.d0/16.d0 - Fs**2 / 9.d0) * t**(-1.d0/3.d0) - Fs / 3.d0 + 0.75d0 + CMPLX(0, sqrt(3.d0)/2.d0) * (  (1.d0/12.d0) *  t**(1.d0/3.d0) +  12.d0 * (-Fs / 6.d0 - 1.d0/16.d0 - Fs**2 / 9.d0) * t**(-1.d0/3.d0))

          dL2 = - (1.d0/24.d0) * t**(1.d0/3.d0) +  6.d0 * (-Fs / 6.d0 - 1.d0/16.d0 - Fs**2 / 9.d0) * t**(-1.d0/3.d0) - Fs / 3.d0 + 0.75d0 - CMPLX(0, sqrt(3.d0)/2.d0) * (  (1.d0/12.d0) *  t**(1.d0/3.d0) +  12.d0 * (-Fs / 6.d0 - 1.d0/16.d0 - Fs**2 / 9.d0) * t**(-1.d0/3.d0))

          !want the smallest of the three roots.  think this is usually dL1 but might change
          
          if (del < 0d0) then
             !one real root and two complex
             WLC_ruppture_length = dL0 * L0
          else
             !all three roots are real.  we want the smallest
             if (( abs(dL0) <= abs(dL1)) .and. ( abs(dL0) <= abs(dL2))) then
                WLC_ruppture_length = dL0 * L0
             elseif (( abs(dL1) <= abs(dL0)) .and. ( abs(dL1) <= abs(dL2))) then
                WLC_ruppture_length = dL1 * L0
             else
                WLC_ruppture_length = dL2 * L0
             end if
          end if
        end function WLC_ruppture_length
        
        subroutine InitializeHystereticModels(Z0)
          real*8, intent(in) :: Z0

          if (Z0 > cur_props%aDMT) then
             CapAdmode = 0
             mJKR = 0
             WLC_mode = 0
          else
             CapAdmode = 1
             mJKR = 1
             WLC_mode = 1
          end if

          SurfHyst_Y  = 0
          SurfHyst_Y0  = 0
          SurfHystMode = 0
          SurfHyst_d0 = Z0

          WLC_prev = WLC_mode
          mJKR_prev = mJKR
          CapAd_prev = CapAdmode
          SurfHystMode_prev = SurfHystMode

          if (cur_props%wantPreisach) then 
             call Initialize_Preisach(Z0)
          end if

          call Initialize_Custom_MD
        end subroutine InitializeHystereticModels

        !for outside res1, call with the _prev variable.  for inside res1, the current variables.
        !this allows us to do the same thing in both places w/o having to write it twice
        subroutine updateSurfHystModel(d, dpr, Y0, Y, d0, mode)
          real*8, intent(in) :: d,dpr
          real*8, intent(out) :: Y0, Y, d0
          integer, intent(inout) :: mode

             if ( (mode == 0) .and. (dpr > 0)) then
                !tip has been going down (force going up).  now we are starting to go back up

               if (d < cur_props%aDMT) then
                   !if we have fully decaying on to the down trajectory, then this is zero.
                   !if we have not fully reached the down traject, then there is non-zero.
                   Y0 =  Fts_nc_vel_ind(d)
                
                   !this is the difference between the current state and the desired up trajectory. 
                   !we need to decay this amount
                   Y = (-4 * pi * Rtip_nd * cur_props%SurfHyst_gamma) - Fts_nc_vel_ind(d)
                else
                   Y0 = Fts_nc_vel_ind(d) * (d**2 / cur_props%aDMT**2)
                   Y  = (-4 * pi * Rtip_nd * cur_props%SurfHyst_gamma ) - Fts_nc_vel_ind(d) * (d**2 / cur_props%aDMT**2)
                end if

                mode = 1
                d0 =  d
             elseif ( (mode == 1) .and. (dpr < 0)) then
                !tip has been doing up (force going down). now we start to go down.                
                if ( d < cur_props%aDMT) then
                   Y =  Fts_nc_vel_ind(d)
                else
                   Y =   Fts_nc_vel_ind(d) * (d**2 / cur_props%aDMT**2)
                end if
                Y0 = 0
                mode = 0
                d0 =  d
             end if
           end subroutine updateSurfHystModel

        !put anything that needs to be updated once per computed point here 
        !(not anything that should be computed for intermediate points though)
        subroutine updateHystereticModels_outsideRES1(d, dpr,t, IOUT, attardY)
          use params
          real*8, intent(in) :: d, dpr,t
          real*8, intent(out) :: attardY(:)
          integer*8, intent(in) ::IOUT

          
          if (cur_props%fts_model == JKR) then
             IF (d .LE. cur_props%aDMT) mJKR_prev = 1
             IF (d .GT. cur_props%dJKRcrit) mJKR_prev = 0
          elseif (cur_props%fts_model == CUSTOM_MD) then
             call UpdateCustomMD_OutsideRES1(d, dpr, cur_props)
          elseif (cur_props%fts_ve_model == VEM_ATTARD_BASED) then
             call CheckAttardStartStop(d, dpr,t, attardY, cur_props)
          end if

          if (cur_props%wantSurfHyst) then
             call updateSurfHystModel(d, dpr, SurfHyst_Y0_prev,  SurfHyst_Y_prev,  SurfHyst_d0_prev,  SurfHystMode_prev)
          end if

          if (cur_props%wantCapAd) then
             IF (d .LE. cur_props%aDMT) CapAd_prev = 1 
             IF (d .GT. cur_props%D_0) CapAd_prev = 0 
	  end if
             
          if (cur_props%WantPreisach) then 
             call UpdatePreisach_outsideRES1(d, IOUT)
          end if
          
          if (cur_props%wantWLC) then
             IF (d .LE. cur_props%aDMT) WLC_prev = 1 
             IF (d .GT. cur_props%wlc_Lr) WLC_prev = 0 
          end if

        end subroutine updateHystereticModels_outsideRES1

        subroutine updateHystereticModels_insideRES1(d, dpr)
          use params
          real*8, intent(in) :: d, dpr

!fixme, might want to put something in here for the MD stuff.  

          ! The calculation of CapAdMode must remain inside of RES1 so that each time Fts_CapAd is 
          ! called, the code will know whether or not the Hysteresis capillary force is acting. 
          ! The PrevCapAd integer is used as a reference to know whether the capillary force was 
          ! previously acting.  This PrevCapAd is used because RES1 can be called in 
          ! non-chronological order.  (kober's model using delta_E & D_0) 
          if (cur_props%wantCapAd) then 
             IF (d .LE. cur_props%aDMT) THEN 
                CapAdMode = 1 
             ELSE IF (d .GT. cur_props%D_0) THEN 
                CapAdMode = 0 
!this is what was in here before jun 29, 2011, but I honestly can't remember why it would have been like this.
!             ELSE IF ((d .GT. cur_props%aDMT) .AND. (CapAd_prev .EQ. 0)) THEN 
!                CapAdMode = 0 
!this makes a lot more sense on the face of it.  
             else
                CapAdMode = CapAd_prev
             END IF
          else
             CapAdMode = 0 !could be more elegant about this, but need to make sure that we handle subs->feat->subs correctly
          end if

          if (cur_props%wantSurfHyst) then
                SurfHystMode = SurfHystMode_prev
                SurfHyst_d0  = SurfHyst_d0_prev
                SurfHyst_Y   = SurfHyst_Y_prev 
                SurfHyst_Y0   = SurfHyst_Y0_prev 
                call updateSurfHystModel(d, dpr, SurfHyst_Y0,  SurfHyst_Y,  SurfHyst_d0,  SurfHystMode) !will overwrite above 4 lines if needed
          end if

          
          if (cur_props%fts_model == JKR) then
             IF (d .LE. cur_props%aDMT) THEN
                mJKR = 1
             ELSE IF (d .GT. cur_props%dJKRcrit) THEN
                mJKR = 0
!             ELSE IF ((d .GT. cur_props%aDMT) .AND. (mJKR_prev .EQ. 0)) THEN  !see above
!               mJKR = 0
             else
                mJKR = mJKR_prev
             END IF
          else
             mJKR = 0 !could be more elegant about this, but need to make sure that we handle subs->feat->subs correctly
          end if          

          if (cur_props%wantWLC) then
             IF (d .LE. cur_props%aDMT) then
                WLC_mode = 1 
             else IF (d .GT. cur_props%wlc_Lr)  then
                WLC_mode = 0 
             else
                WLC_mode = WLC_prev
             end IF
          end if


          if (cur_props%wantPreisach) then 
             call UpdatePreisach_insideRES1(d)
          end if
        end subroutine updateHystereticModels_insideRES1

        !technically, there is more stuff in here than just non-dimensionalization so its badly named
        subroutine NormalizeMatlProps( props, Rtip_dim)
          use data1, only: nummodes
          use data0
          use Elasticity_Formulas

          real*8, intent(in) :: Rtip_dim

          type(matl_prop), intent(inout) :: props
	  real*8 :: C1_d, C2_d, g_d
          integer :: i
          character*200 :: tmpstr

          if (props%fts_model .eq. MAG_DIPOLE) then
             props%mag_z0 = NonDimenLength(Rtip_dim/2d0 + props%mag_dia/2d0 + props%mag_delta)
          end if

          if ((props%fts_model .eq. DMT) .or. (props%fts_model .eq. DMT_DLVO)) then
             select case (props%calcADMT)
                case(ENTER_FAD_ONLY)
                   props%aDMT = 0.2e-9                   
                   props%A_Hamaker_dim = (6 * props%Fadhesion) * props%aDMT**2 / Rtip_dim 
                case( ENTER_FAD_A0)
                   props%A_Hamaker_dim = (6 * props%Fadhesion) * props%aDMT**2 / Rtip_dim 
                case (ENTER_FAD_H)
                   props%aDMT = sqrt( props%A_Hamaker_dim * Rtip_dim / (6 * props%Fadhesion)) !dimensional
                
                   if ( props%aDMT > 2e-9 ) then
                      call WriteFatalError( "Provided values for Hamaker constant, Tip radius, and van der Waals adhesion force lead to unrealistically high intermolecular distance." // char(10) // "The expected value is on the order of several Angstroms but calculated greater than 2 nm.  Check input values and try again.")
                   elseif (props%aDMT < 0.01e-9) then
                      call WriteFatalError( "Provided values for Hamaker constant, Tip radius, and van der Waals adhesion force lead to unrealistically low intermolecular distance. " // char(10) // " The expected value is on the order of several Angstroms but calculated less than 0.1 Angstrom.  Check input values and try again.")
                   end if

                case (ENTER_H_A0)
                   !nothing to do
                end select
                
             end if

          if  ((props%fts_model == DMT) .or.  (props%fts_model == DMT_DLVO)) then
             props%aDMT = NonDimenLength(props%aDMT)
          else
             props%aDMT = 0d0		
          end if

          if (props%VEchoice == VEC_NONE) then
             ! Estar is the effective Young's modulus (dimensional)         
             props%Estar           = reducedModulus( Etip, Poisson_tip, props%Esample, props%Poisson_sample)
          end if

          !lennard jones and also in Attard, so just always do it.  
          props%LJ_r0_nd = NonDimenLength(props%LJ_r0_dim)

          select case( props%fts_model) 
             case (LIN_ATT)
                props%Lo = props%aDMT + NonDimenLength(props%Fadhesion/props%kts_A)
             case (DMT_DLVO)
                props%KD_nd = NonDimenInvLength(props%KD_dim)
             case (MORSE)
                props%MorseRc_nd = NonDimenLength(props%MorseRc_dim)
                props%MorseLambda_nd = NonDimenLength(props%MorseLambda_dim)
             case (JKR)
                call assert(props%Estar > 0, 'Estar < 0')

                props%W_JKR = 2d0*props%Fadhesion/(3d0*pi*Rtip_dim)
                props%acrit_JKR = NonDimenLength(((    pi*(Rtip_dim**2)*props%W_JKR/(8d0*props%Estar))**(1d0/3d0)))
                props%a0_JKR    = NonDimenLength(((2d0*pi*(Rtip_dim**2)*props%W_JKR/     props%Estar)**(1d0/3d0)))
                C1_d = sqrt(NonDimenLength(2d0*pi*props%W_JKR*(props%acrit_JKR)/(props%Estar)))
                C2_d = (((props%acrit_JKR)**2)/Rtip_nd) - props%aDMT
                props%dJKRcrit = C1_d - C2_d
                props%a_JKR = props%a0_JKR                
             case (ELECTROSTATIC_GIL)
                es_h = NonDimenLength(es_h_dim)
                es_l = NonDimenLength(es_l_dim)
                es_w = NonDimenLength(es_w_dim)

                es_delta = Rtip_nd / (2 * tan( theta_tip / 2)**2)
                es_c4 = 2 * es_l * tan( theta_lever / 2)
                es_c3 = sin( theta_tip / 2) * ( es_h -2 * es_delta)

             case (ELECTROSTATIC_NONCONS)
                es_l = NonDimenLength(es_l_dim)
                es_w = NonDimenLength(es_w_dim)
                g_d = NonDimenLength(props%g_d_dim)
                props%area = es_l_dim*es_w_dim
              !  area = es_l*es_w
              !  a_dim
              !  a
              !  capac = (area*epsilon_0*epsilon_d)/(epsilon_0*g_d+epsilon_0*d) 
              !  dcdq = -(area*epsilon_0*epsilon_d^2)/(epsilon_0*g_d+epsilon_d*d)^2
           end select

           props%wantPreisach = .false.

           if (props%Want_hyst_hydr) then
              props%wantPreisach = .true.
              !fixme, replace this with a compiled matlab script on nanohub !  -nojvm -nosplash
              write( tmpstr, *) '../../src/preisach_generator/run_preisach.sh /opt/matlab ',  props%app_decay, props%ret_decay, props%app_scaling, props%nside, props%cutoff_dist
              
              !fixme, this gets long.  find a way to only re-run if the table has changed
              call system(tmpstr)
           end if
           

           if (props%wantPreisach) then
              call readCustomPreisachTable( 'preisach.txt')
              call Normalize_Preisach()
           end if

           props%SurfHyst_gamma = NonDimenInvLengthSq( NonDimenEnergy(  props%SurfHyst_gamma )) 
           
           props%D_0 = NonDimenLength(props%D_0)

           if (props%wantOscillatory) then
              props%Rtip_solvation_nd = NonDimenLength(props%Rtip_solvation_dim)           
              props%sigma_solvation_nd =   NonDimenLength(props%sigma_solvation_dim)
           else
              props%Rtip_solvation_nd =0d0
              props%sigma_solvation_nd=0d0
           end if

           if (props%wantHydration) then
              props%lambda_solvation_nd =  NonDimenLength(props%lambda_solvation_dim)
           else
              props%lambda_solvation_nd =  0d0
           end if

           if (props%want_exp_dashpot) then
              props%exp_dashpot_decay = NonDimenLength(props%exp_dashpot_decay)              
           end if
          
           SurfHyst_decay_nd = NonDimenLength( SurfHyst_decay_dim)

           !visoelastic non-conservative option is linked to conservative option

           if ((props%fts_model == ATTARD_FOURIER_BAHRAM) .or. (props%fts_model == ATTARD_FOURIER_LSQ)) then
              props%fts_ve_model = VEM_ATTARD_BASED
           elseif (props%VEchoice == VEC_NONE) then
              props%fts_ve_model = VEM_NONE
           else
              if ((props%fts_model == LINEAR) .or. (props%fts_model == LIN_ATT)) then
                 props%fts_ve_model = VEM_LINEAR_BASED
              elseif (( props%fts_model == HERTZ) .or. ( props%fts_model == DMT)) then
                 props%fts_ve_model = VEM_HERTZ_BASED
              else
                 props%fts_ve_model = VEM_NONE
              end if
           end if

           if ((props%fts_ve_model /= VEM_NONE) .and. (props%fts_ve_model /= VEM_ATTARD_BASED)) then
              call NonDimensionalizeViscoelasticModels(props)
           end if

           if (( props%fts_model == ELECTROSTATIC_XU ) .or. ( props%fts_model == ELECTROSTATIC_GIL) .or. ( props%fts_model == ELECTROSTATIC_NONCONS)) then
              props%bias_freq_rads = NonDimenFreq( props%bias_freq_rads )
           end if

          props%resis = 1e2
          props%epsilon_0 = 1*8.854e-12
          props%epsilon_d = 5*props%epsilon_0
          props%g_d_dim = 40e-9

          if (props%wantWLC) then
             props%wlc_Lr = nonDimenLength( WLC_ruppture_length( props%wlc_Fr, props%wlc_Lp, props%wlc_L0, props%Temp))

             props%wlc_L0 = nonDimenLength(props%wlc_L0)

             call assert (props%wlc_Lr <= props%wlc_L0, "worm like chain LR must be less than L0")
          end if
          !Lp gets non-dimen as part of the big force term

          if (props%fts_ve_model == VEM_ATTARD_BASED) then
             call initialize_attard(Rtip_nd, props )
          end if
        end subroutine NormalizeMatlProps


        real*8 function CalcContactTime(t, d, tc)  !keeps a running total of contact time for each output point
          real*8, intent(in) :: t,d, tc

           if ((check .eqv. .false. ) .AND. (d .lt. cur_props%aDMT)) then
              t1 = t ! just entered contact on this point
           end if

           if ((check .eqv. .true. ) .AND. (d.gt. cur_props%aDMT)) then
              CalcContactTime = tc + (t-t1)  !just left contact on this point
           else
              CalcContactTime = tc
           end if

           if (d.lt.cur_props%aDMT) then
              check = .true.
           else
              check = .false.			
           end if
         end function CalcContactTime
         
         recursive  function Fts_oscillatory(d) result(foo)
           real*8, intent(in) :: d !d has already had aDMT subtracted out
           real*8 :: foo, tmp, m, zk
           integer :: i

           !jeffrey et al. prB. 2004 assumed a spherical tip broken into N slices.  
           !this series does not appear to converge though!!
           !in the GUI we hard coded N=1 which makes this into just a decaying cosine
           !fixme: really should just take out the series
           if ((d < 0d0) .and. (.not. cur_props%WantContGrad)) then
              foo = Fts_oscillatory(0d0)
           else
              tmp = 0
              do i = 1,cur_props%N_solv
                 m = real(i)

                 zk = (m-1)*cur_props%Rtip_solvation_nd/cur_props%N_solv

                 tmp = tmp + (1d0-((cur_props%N_solv-m)/cur_props%N_solv)**2)*cos(2 *pi*(d+zk)/cur_props%sigma_solvation_nd)*exp(-(d+zk) / cur_props%sigma_solvation_nd)
              
              end do
              foo = cur_fC%Coscillatory*tmp
           end if

         end function Fts_oscillatory

         recursive  function Fts_hydration(d) result(foo)
           real*8, intent(in) :: d !d has already had aDMT subtracted out           
           real*8 :: foo, mm 
           integer :: i

           if ((d < 0d0) .and. (.not. cur_props%WantContGrad)) then
              foo = Fts_hydration(0d0)
           else
              foo = cur_fC%Chydration * exp(-d/ cur_props%lambda_solvation_nd)             
           end if
         end function Fts_hydration


         real*8 function tip_sample_voltage(t)
           real*8, intent(in) :: t

           if ( cur_props%fts_model == ELECTROSTATIC_NONCONS ) then
              call assert(.false., 'shouldnt be here')
           else
              tip_sample_voltage = cur_props%dc_bias_voltage - cur_props%surface_pot +  cur_props%ac_bias_voltage * cos(  cur_props%bias_freq_rads * t)
           end if
          
         end function tip_sample_voltage

         !  Tip-sample interaction function, conservative part
         real*8 function Fts_cons(gap, t,  fail)            
           real*8, intent(in) :: gap, t
           logical, intent(out) :: fail
           real*8 :: chi
           
           fail = .false.
           
           select case( cur_props%fts_model) 
           case (LINEAR)
              if (gap.lt.cur_props%aDMT) then
                 if ((cur_props%VEchoice == VEC_NONE) .or. (cur_props%VEchoice == VEC_KELVINVOIGT_LEGACY) .or. (cur_props%VEchoice == VEC_KELVINVOIGT)) then
                    Fts_cons = cur_fC%Clin*(-gap) 
                 else
                    Fts_cons = 0d0 !see comments under Hertz
                 end if
              else
                 Fts_cons = 0d0
              end if
	   case (LIN_ATT)
              if (gap.lt. 0d0) then
                 Fts_cons = cur_fC%Clin*(-gap) - cur_fC%Fad
              else if (gap.lt. cur_props%Lo) then
                 Fts_cons = cur_fC%C_att*(gap) - cur_fC%Fad
	      else
	         Fts_cons = 0d0
              end if
           case (CHADWICK)
              if (gap.lt. 0d0) then
                 Fts_cons=cur_fC%CCHAD*(-gap)**2
              else
                 Fts_cons = 0d0
              end if
           case (HERTZ_BEC)
              if (gap.lt. 0d0) then
                 if (tip_shape == PARABOLOID) then
                    chi = sqrt( Rtip_nd * (-gap) ) / NonDimenLength( cur_props%hs)                    
                    Fts_cons=cur_fC%hertz_bec * (-gap)**1.5d0 * (1d0+0.884d0*chi+0.781d0*chi**2+0.386d0*chi**3+0.0048d0*chi**4)
                 else !cone
                    chi = (-gap) *  tan(tip_angle) / NonDimenLength( cur_props%hs)
                    Fts_cons=cur_fC%hertz_bec *(-gap)**2 * (1d0 + 1.7795d0 * 2d0 * chi / pi**2   + 16d0 * 1.7795**2 * chi**2   )
                 end if
              else
                 Fts_cons = 0d0
              end if
           case (HERTZ)
              if ((cur_props%VEchoice /= VEC_NONE) .and. (cur_props%VEchoice /= VEC_KELVINVOIGT_LEGACY))  then
                 !for viscoelasticity, the conservative and nc parts all get lumped into 
                 !the same integral.  
                 Fts_cons = 0d0                 
              else
                 if (gap.ge. 0d0) then
                    Fts_cons = 0d0
                 else
                    Fts_cons = Fts_Hertz(gap, cur_fC%CDMT)
                 end if
              end if
           case (MORSE)
              if ( gap <= 0) then
                 fail = .true.  !zero can't happen for morse.
                 Fts_cons = 0d0
              else                 
                 Fts_cons = -cur_fC%CvdW/gap**2 -cur_fc%Cmorse * ( -exp(( cur_props%MorseRc_nd - gap) /cur_props%MorseLambda_nd) + exp(2 * ( cur_props%MorseRc_nd - gap) /cur_props%MorseLambda_nd))
              end if
           case (LENNARD_JONES)
              if ( gap <= 0) then
                 fail = .true.  !zero can't happen for lj
                 Fts_cons = 0d0
              else                 
                 Fts_cons = -cur_fC%CvdW/gap**2 -cur_fc%CLJ * ( ( cur_props%LJ_r0_nd / gap)**13 - (cur_props%LJ_r0_nd / gap)**7)
             end if
           case (ELECTROSTATIC_XU)
              if ( gap <= 0) then
                 fail = .true.  !avoid divide by zero                                  
              else
                 Fts_cons = tip_sample_voltage(t)**2 * cur_fC%Celectrostatic_xu / ( gap**2)
              end if
           case (ELECTROSTATIC_GIL)
              if ( gap <= 0) then
                 fail = .true. 
              else
                 Fts_cons =  -cur_fC%CvdW/gap**2 -cur_fC%Ces_gil_c5 * ( tip_sample_voltage(t)**2 * cur_fC%Ces_gil_c1 * log(( gap - es_delta + es_h)/(gap + es_delta)) -  tip_sample_voltage(t)**2 * cur_fC%Ces_gil_c1 * es_c3 * (gap - es_delta) / ((gap - es_delta + es_h)*(gap + es_delta)) +  tip_sample_voltage(t)**2 * cur_fC%Ces_gil_c2 / (( 1 + gap / es_h)*(1+ (gap + es_c4)/es_h)))
              end if
          case (CUSTOM_CONS)
             Fts_cons = Fts_custom_cons( gap)
          case (DMT_DLVO)
             Fts_cons = Fts_DMT(gap) + cur_fC%CDLVO*exp(-cur_props%KD_nd*gap)
          case (DMT)
             if ((cur_props%VEchoice == VEC_NONE) .or. (cur_props%VEchoice == VEC_KELVINVOIGT_LEGACY)) then
                Fts_cons = Fts_DMT(gap)
             else
                !for viscoelasticity, the in contact portion of this force is computed in the viscoelasticity mod
                if (gap >= cur_props%aDMT) then
                   Fts_cons = -cur_fC%CvdW/gap**2
                else
                   Fts_cons = -cur_fC%CvdW/(cur_props%aDMT**2)
                end if                
             end if
          case (ELECTROSTATIC_NONCONS)
             Fts_cons = 0d0
          case (MAG_DIPOLE)
             if ( (gap + cur_props%mag_z0) <= 0) then
                fail = .true.  !avoid divide by zero                                  
             else
                Fts_cons = cur_fC%C_mag_dipole / ( gap + cur_props%mag_z0)**4
             end if

          case default
              Fts_cons = 0d0
          end select

           if (cur_props%wantOscillatory) Fts_cons = Fts_cons + Fts_oscillatory(gap - cur_props%aDMT)
           if (cur_props%wantHydration) Fts_cons = Fts_cons + Fts_hydration(gap - cur_props%aDMT)
                      
                      
         End function Fts_cons


         real*8 function Fts_DMT(gap)
           real*8, intent(in) :: gap           
           if (gap >= cur_props%aDMT) then
              Fts_DMT = -cur_fC%CvdW/gap**2
           else
              Fts_DMT = -cur_fC%CvdW/(cur_props%aDMT**2) + Fts_Hertz(gap - cur_props%aDMT, cur_fC%CDMT)
           end if
         end function Fts_DMT
         

         real*8 function Fts_custom_cons( gap)          
           real*8, intent(in) :: gap

           integer i, N
 
           N = cur_props%N_custom
           
           if ( gap >= cur_fC%gap_table_nd(N) ) then
              !must extrapolate!
              Fts_custom_cons = cur_fC%F_table_nd(N) + (gap - cur_fC%gap_table_nd(N)) * (cur_fC%F_table_nd(N) - cur_fC%F_table_nd(N-1) ) / (cur_fC%gap_table_nd(N) - cur_fC%gap_table_nd(N-1))
           elseif ( gap <= cur_fC%gap_table_nd(1) ) then
              !must extrapolate
              Fts_custom_cons = cur_fC%F_table_nd(1) + (gap - cur_fC%gap_table_nd(1)) * (cur_fC%F_table_nd(2) - cur_fC%F_table_nd(1)) / (cur_fC%gap_table_nd(2) - cur_fC%gap_table_nd(1))
           else
              !interpolate
              !fixme, make this a binary search to be quicker
              do i = 1, N-1
                 if ( gap < cur_fC%gap_table_nd(i+1)) then
                    Fts_custom_cons = cur_fC%F_table_nd(i) + (gap - cur_fC%gap_table_nd(i)) * (cur_fC%F_table_nd(i+1) - cur_fC%F_table_nd(i) ) / (cur_fC%gap_table_nd(i+1) - cur_fC%gap_table_nd(i))
                    exit ! do loop
                 end if
              end do
           end if
         end function Fts_custom_cons
	 
	 real*8  function Fts_CapAd(gap)
           real*8 :: gap

           if (gap < cur_props%aDMT) then              
              Fts_CapAd = cur_fC%C_CapAd*(cur_props%aDMT-cur_props%D_0)
           else
              ! john + steven's original interpretation was that this equation below held
              ! whenever the force was on.  now belief that this equation only holds between
              ! aDMT and D_0, with it being a constant value below aDMT.  Kober is not really
              ! very clear on exactly what his model is. in the end may not make much
              ! difference...
              Fts_CapAd = cur_fC%C_CapAd*(gap-cur_props%D_0)
           end if
           
         End function Fts_CapAd
	 
	 
         real*8  function Fts_nc_vel_dep(d,dpr, t, y, fail)
           use data1, only: numModes, NEQ
           use timeAndCycle, only : operating_mode
           real*8, intent(in) :: d, dpr, t, y(:)
           logical, intent(out) :: fail

           fail = .false.

           if ((cur_props%VEchoice == VEC_NONE) .or. (cur_props%fts_ve_model == VEM_ATTARD_BASED)) then
              Fts_nc_vel_dep = 0d0
           else
              if (d < cur_props%aDMT) then

                 select case (cur_props%VEchoice)
                 case(VEC_KELVINVOIGT_LEGACY)
                    select case (cur_props%fts_ve_model)
                    case(VEM_HERTZ_BASED)
                       Fts_nc_vel_dep =  -cur_fC%Dts*sqrt(cur_props%aDMT-d)*dpr
                    case(VEM_LINEAR_BASED)
                       Fts_nc_vel_dep = -cur_fC%Dts*dpr 
                    end select
                 case(VEC_KELVINVOIGT)
                    select case (cur_props%fts_ve_model)
                    case(VEM_HERTZ_BASED)
                       Fts_nc_vel_dep = getViscoelasticForce_Hertz(t, d-cur_props%aDMT, dpr, cur_props, cur_fC, operating_mode)
                    case(VEM_LINEAR_BASED)
                       ! if ( cur_props%fts_model == LINEAR) then
                       !for the hertz model, there is a fairly complicated viscoelastic theory (see manual).
                       !for the linear model, we just impose the ad-hoc condition that the total
                       !force cannot go negative.  obviously, if the force goes negative then
                       !the tip should just leave contact b/c there is no attractive force to hold it on. 
                       !this is a different behavior than previous versions
                          if ( Fts_cons(d, t, fail)  -cur_fC%Dts*dpr > 0) then
                             Fts_nc_vel_dep = -cur_fC%Dts*dpr 
                          else
                             Fts_nc_vel_dep = - Fts_cons(d, t, fail)  !want total force zero.  cancel it out
                          end if
                       ! elseif  ( cur_props%fts_model == LIN_ATT) then
                       !    !in this situation, there could be a limited attractive force
                       !    if ( Fts_cons(d, t, mode, fail)  -cur_fC%Dts*dpr > -cur_fC%Fad) then
                       !       Fts_nc_vel_dep = -cur_fC%Dts*dpr 
                       !    else
                       !       Fts_nc_vel_dep = - Fts_cons(d, t, mode, fail) -cur_fC%Fad 
                       !    end if
                       !end if                     
                    case default
                       call assert(.false., 'unknown vem option')
                       Fts_nc_vel_dep = 0d0
                    end select
                 case (VEC_MAXWELL, VEC_THREEELM_E1E2ETA, VEC_GENMAXWELL )                          
                    if (cur_props%fts_ve_model == VEM_HERTZ_BASED) then
                       Fts_nc_vel_dep = getViscoelasticForce_Hertz(t, d-cur_props%aDMT, dpr, cur_props, cur_fC,  operating_mode)  
                    elseif  (cur_props%fts_ve_model == VEM_LINEAR_BASED) then
                       Fts_nc_vel_dep = GetViscoelasticForce_linear(t, d, operating_mode, cur_props)
                    end if
                 case default
                       call assert(.false., 'unknown vec model')
                       Fts_nc_vel_dep = 0d0
                 end select
              else
                 Fts_nc_vel_dep = 0d0
              end if
           end if

           if (cur_props%want_tip_squeeze) then
              if (d >= NonDimenLength(0.2d-9) ) then
                 Fts_nc_vel_dep = Fts_nc_vel_dep - cur_fC%tip_squeeze * dpr / d
              elseif (d > 0) then
                 !Fts_nc_vel_dep = Fts_nc_vel_dep - cur_fC%tip_squeeze * dpr / (6d0 * NonDimenLength(0.2d-9)- 5d0 * d)
                 Fts_nc_vel_dep = Fts_nc_vel_dep - cur_fC%tip_squeeze * dpr / (NonDimenLength(0.2d-9)) * (d / NonDimenLength(0.2d-9))
                 !discontinuous force does not play well with sovling diff eq.  smooth it out.
              end if
           end if

           if (cur_props%want_exp_dashpot) then
              if (d.ge.cur_props%aDMT) then
                 Fts_nc_vel_dep = Fts_nc_vel_dep - cur_fC%D_exp * exp( -(d - cur_props%aDMT) / cur_props%exp_dashpot_decay ) * dpr
              else
                 !prior to 5/7/2011, had this force drop instantly to zero aDMT.  leads to discontinuities though.  
                 !then, prior to 5/23/2012, had it decay smoothly to zero (exp). matched experimental data okay,
                 !but ad-hoc and hard to justify. per arvind's suggestion, just have the viscosity saturate.
                 !it's the same as the DMT assumption that the attractive forces continue to act in a zone
                 !just outside the contact area
                 Fts_nc_vel_dep = Fts_nc_vel_dep - cur_fC%D_exp * dpr
              end if
           end if


           if (cur_props%fts_ve_model == VEM_ATTARD_BASED) then
              if (cur_props%N_attard_fourier == 0) then
                 Fts_nc_vel_dep = Fts_nc_vel_dep + Fts_attard_spatial(d,t, Y(2*nummodes+1 : NEQ ), cur_props)
              else
                 Fts_nc_vel_dep = Fts_nc_vel_dep + Fts_attard_fourier(d,t, Y(2*nummodes+1 : NEQ ), cur_props)
              end if
           end if           
         End function Fts_nc_vel_dep

         real*8 function Fts_nc_vel_ind(d)

           real*8, intent(in) :: d

           real*8 :: JKR_err, f_1, df, delt


           if  ( cur_props%fts_model == JKR)  then
	      if (mJKR .gt. 0) then
		 JKR_err = 1d10
		 DO WHILE (JKR_err .GT. 1d-11)
		   f_1 = cur_props%aDMT-d-cur_fC%C_deltJKR1*cur_props%a_JKR**2+cur_fC%C_deltJKR2*sqrt(cur_props%a_JKR)
		   df = -2*cur_fC%C_deltJKR1*cur_props%a_JKR + cur_fC%C_deltJKR2/(2*sqrt(cur_props%a_JKR))
		   delt = f_1/df
		   cur_props%a_JKR = cur_props%a_JKR - f_1/df
		   JKR_err = abs(delt/cur_props%acrit_JKR)
		 END DO
		 Fts_nc_vel_ind = cur_fC%C_Fjkr1*cur_props%a_JKR**3 - cur_fC%C_Fjkr2*sqrt(cur_props%a_JKR**3)
	      else
                 Fts_nc_vel_ind = 0d0
              end if
           elseif (cur_props%fts_model == CUSTOM_MD) then
              Fts_nc_vel_ind = Fts_custom_MD(d, cur_props)
           else
              Fts_nc_vel_ind = 0d0
           end if

           if (cur_props%wantPreisach) then
              Fts_nc_vel_ind = Fts_nc_vel_ind + Fts_Preisach()
           end if

           if (cur_props%wantCapAd) then 
              Fts_nc_vel_ind = Fts_nc_vel_ind + CapAdMode*Fts_CapAd(d)
           end if

           if ((cur_props%wantSurfHyst) .and. (d < cur_props%aDMT)) then
              if (SurfHystMode == 0) then !going down. dpr < 0                 
                 Fts_nc_vel_ind = Fts_nc_vel_ind + SurfHyst_Y * exp( (d-SurfHyst_d0) / SurfHyst_decay_nd)
              elseif ((SurfHystMode == 1) .and. ((SurfHyst_Y /= 0) .or. (SurfHyst_Y0 /= 0))) then !going up. dpr > 0
                 Fts_nc_vel_ind = Fts_nc_vel_ind + SurfHyst_Y0 + SurfHyst_Y * (1 - exp(- (d-SurfHyst_d0) / SurfHyst_decay_nd))
              end if
           elseif  ((cur_props%wantSurfHyst) .and. (d > cur_props%aDMT)) then
              if     ((SurfHystMode == 0) .and. (SurfHyst_Y /= 0)) then !going down
!fixme, do not work for Hertz b/c aDMT = 0!
                 Fts_nc_vel_ind = Fts_nc_vel_ind +  (cur_props%aDMT**2 / d**2) * SurfHyst_Y * exp((d-SurfHyst_d0) / SurfHyst_decay_nd)
              elseif ((SurfHystMode == 1) .and. ((SurfHyst_Y /= 0) .or. (SurfHyst_Y0 /= 0))) then !going up
                 Fts_nc_vel_ind = Fts_nc_vel_ind + (cur_props%aDMT**2 / d**2) * ( SurfHyst_Y0 + SurfHyst_Y * (1d0-exp(- (d-SurfHyst_d0) / SurfHyst_decay_nd)))
              end if
           end if

           !previously was incorrectly in the "conservative" function
           if ((cur_props%wantWLC) .and. (WLC_mode == 1)) then                 
                 !fixme, pull out the constants to a fC
                 Fts_nc_vel_ind =  Fts_nc_vel_ind - NonDimenForce((kB * cur_props%Temp / cur_props%wlc_Lp) * ( 0.25d0 / (1d0 - 1d0 * d / cur_props%wlc_L0)**2 - 1/4d0 + d / cur_props%wlc_L0))
           end if

           
         end function Fts_nc_vel_ind

         ! Fts is now separated properly into cons and nc.  
         ! can we integrate over just Fts_nc for finding Ets?  might be more accurate...
         
         real*8 function Fts(d, dpr, fail,t, y)
           real*8, intent(in) :: d, dpr,t, y(:)
           logical, intent(out) :: fail
           logical :: fail1, fail2

           Fts = Fts_cons(d, t, fail1) + Fts_nc_vel_dep(d, dpr,t, y, fail2) + Fts_nc_vel_ind(d)
           fail = fail1 .or. fail2
         end function Fts

         type(forceCoeff) function computeForceCoeff( mp )
           use data0
           use Nondimensionalization
           use Elasticity_Formulas

            type(matl_prop) , intent(in) :: mp
                       
            real*8 :: G, prefactor
            integer :: j

           allocate( computeForceCoeff%gap_table_nd( mp%N_custom) )
           allocate( computeForceCoeff%F_table_nd( mp%N_custom) ) 
           
!see, Dissipation and oscillatory solvation forces in confined liquids studied by small-amplitude atomic force spectroscopy.  debeer, ende, mugule. 2010
           computeForceCoeff%tip_squeeze = 6d0 * pi * NonDimenForce( NonDimenTime(NonDimenInvLengthSq(mp%eta_liquid)) * Rtip_nd**2)

           computeForceCoeff%CvdW= NonDimenForce(NonDimenlength(mp%A_hamaker_dim)*Rtip_nd/6.0d0)
!    	   CvdW is the van der Waals force constant

!          CDMT is the Hertz contact constant, also used for DMT and DMT+DLVO models
           if (tip_shape == PARABOLOID) then
              computeForceCoeff%CDMT = NonDimenPressure(4d0/3d0* mp%Estar*sqrt(Rtip_nd))
           else
              computeForceCoeff%CDMT = NonDimenPressure(2d0 * mp%Estar * tan( tip_angle) / pi )
           end if


           computeForceCoeff%C_CapAd = NonDimenForce(NonDimenLength( 2d0 * mp%deltaE/(( mp%D_0**2))))

           if (mp%fts_model == JKR) then
              computeForceCoeff%C_Fjkr1 = NonDimenPressure(4d0*mp%Estar/(3d0*Rtip_nd))
              computeForceCoeff%C_Fjkr2 = NonDimenForce(NonDimenInvLength(2d0*sqrt(NonDimenInvLength(2d0*pi*mp%Estar*mp%W_JKR))))
              computeForceCoeff%C_deltJKR1 = 1/Rtip_nd
              computeForceCoeff%C_deltJKR2 = sqrt(NonDimenLength(2*pi*mp%W_JKR/(mp%Estar)))
           end if
        
           !10/31/10.  this are correct.  get the right energy dissipation versus 
           !an analytical formula
           if (mp%fts_ve_model == VEM_HERTZ_BASED) then            
              if (mp%VEchoice == VEC_KELVINVOIGT_LEGACY) then                 
                 ! etasample is the sample viscosity in Pa-s = (N/m^2)-s.  So Dts is a contact damping parameter. 
                 computeForceCoeff%Dts=  NonDimenPressure( NonDimenTime(sqrt(Rtip_nd) *mp%etasample))
              else
                 computeForceCoeff%Dts=  NonDimenPressure( NonDimenTime(sqrt(Rtip_nd) *mp%etasample)) * 2d0 / (1-mp%Poisson_sample)
              end if
           elseif (mp%fts_ve_model == VEM_LINEAR_BASED) then            
              !etasample is in N-s/m=kg/s, so that the below is normalization as well.  
              computeForceCoeff%Dts= NonDimenStiffness(NonDimenTime(mp%etasample))
           elseif ((mp%fts_ve_model == VEM_ATTARD_BASED) .or. (mp%fts_ve_model == VEM_NONE)) then
              !do nothing.
           else
              call assert(.false., "unknown nc model")
           end if

           if (tip_shape == PARABOLOID) then
              prefactor =  2d0 * sqrt(Rtip_nd)
           else
              prefactor =  4d0 * tan(tip_angle) / pi
           end if

           if (mp%VEchoice == VEC_MAXWELL) then
              computeForceCoeff%vec_f = prefactor * NonDimenPressure(mp%Estar)
           elseif (mp%VEchoice == VEC_THREEELM_E1E2ETA) then
              computeForceCoeff%vec_f = prefactor*NonDimenPressure((mp%Estar * mp%threeelm_E2star) / (mp%threeelm_E2star + mp%Estar))
           else
              !gen mawxwell, we include all of the E modulus into the psi directly, rather that pulling out here.              
              computeForceCoeff%vec_f =  prefactor* NonDimenForce(  NonDimenInvLengthSq(1d0) )
           end if
           
           
           G = YoungsToShear(mp%Esample, mp%Poisson_sample)
           computeForceCoeff%kv_t0 =  NonDimenTime( mp%etasample / (2d0 * G))

           computeForceCoeff%D_exp=  NonDimenStiffness(NonDimenTime(mp%exp_dashpot_scale))

           !custom
           do j = 1, mp%N_custom
              computeForceCoeff%gap_table_nd(j) = NonDimenLength(mp%gap_table_dim(j))
              computeForceCoeff%F_table_nd(j)   = NonDimenForce(mp%F_table_dim(j))
           end do


!	   Linear force coefficient
           computeForceCoeff%Clin = NonDimenStiffness(mp%kts_R)

!	   Linear attractive force coefficient
	   computeForceCoeff%C_att = NonDimenStiffness(mp%kts_A)
	   
!	   Scaled adhesion force
	   computeForceCoeff%Fad = NonDimenForce(mp%Fadhesion)

           computeForceCoeff%CDLVO = NonDimenForce( 4*pi*Rtip_dim*mp%sigmas*mp%sigmat /(mp%epsilon*epsilon0*mp%KD_dim))

	   computeForceCoeff%CCHAD = NonDimenPressure(2*pi*Rtip_dim*mp%Estar /(3d0 *mp%hs))

           if (tip_shape == PARABOLOID) then
              computeForceCoeff%hertz_bec = NonDimenPressure(16d0/9d0* mp%Esample*sqrt(Rtip_nd))
           else
              computeForceCoeff%hertz_bec = NonDimenPressure( 8d0 * mp%Esample * tan(tip_angle) / 3d0 / pi)
           end if

           if (mp%fts_model == MORSE) then
              computeForceCoeff%CMorse = NonDimenForce(NonDimenLength( -2d0 * mp%MorseU0 / (mp%MorseLambda_nd)))
           end if
           if (mp%fts_model == LENNARD_JONES) then
              computeForceCoeff%CLJ =    NonDimenForce(NonDimenLength(-12d0 * mp%LJ_E0   / (mp%LJ_r0_nd))      )
           end if

           !jeffrey is not really clear about the sign of the force, but from israelachvili think it should be attractive at d=0
           if (mp%wantOscillatory) then
              computeForceCoeff%Coscillatory = NonDimenForce(2 * pi * mp%Rtip_solvation_dim**2 * KB * mp%Temp * mp%rho)
           end if

           if (mp%wantHydration) then
              computeForceCoeff%Chydration = NonDimenForce(2 * pi * mp%Rtip_solvation_dim**2 * mp%p_h)
           end if

           computeForceCoeff%Celectrostatic_xu = NonDimenForce(-4d0 * pi * epsilon0 * ((mp%epsilon-1d0) / (mp%epsilon+1d0)) * Rtip_nd**2)

           if (mp%fts_model == ELECTROSTATIC_GIL) then
              computeForceCoeff%Ces_gil_c1 = 4 * pi * mp%epsilon * epsilon0  / (pi - theta_lever)**2
              computeForceCoeff%Ces_gil_c2 = tan( theta_lever)**2 * mp%epsilon * epsilon0 * es_l * es_w / (( theta_lever ** 2) * (es_h ** 2))
              computeForceCoeff%Ces_gil_c5 = NonDimenForce( 1d0)
           end if

           !11/9/2020 these are currently unused
!           computeForceCoeff%Ces_nc_cap1 = mp%area * mp%epsilon_0 * mp%epsilon_d
!           computeForceCoeff%Ces_nc_dcdq1 = -mp%area* mp%epsilon_0* mp%epsilon_d**2

           if (mp%fts_model == MAG_DIPOLE) then
              !need the right non-dim here.
              computeForceCoeff%C_mag_dipole = NonDimenForce( -1.5d0 * mu0 * NonDimenLengthSq( mp%mom_sample) * NonDimenLengthSq(mp%mom_tip) / pi)
           end if

         end function ComputeForceCoeff
         
       end module contactModels

!this is a number of routines that Hank started for his paper.  They were originally directly
!part of the contactModels module, but they were never finished and there is no
!current plan to continue working on them. Temporarily, I'm moving them here.
!long term if there is still no plan to finish these, they should be deleted
! module NonConsElectrostatic
!   contains
!          !this is hank's force
!          real*8 function Fts_electrostatic_noncons(gap, gappr, t,  V)
!            real*8, intent(in) :: gap, gappr, t, V
!            real*8 :: d_dim, dvdq, d_dot

!            d_dim = DimenLength(gap)
!            d_dot = DimenVelocity(gappr)

!            call assert(cur_props%resis*d_dot /= 0, 'd_dot R = 0')

!            dvdq = (cur_props%dc_bias_voltage/capac(d_dim))*dcdq(d_dim) - (V/capac(d_dim))*(1/(cur_props%resis*d_dot) + dcdq(d_dim))

!            Fts_electrostatic_noncons = NonDimenForce(-0.5D0*(dcdq(d_dim)*(cur_props%dc_bias_voltage - V)**2D0 - 2D0*capac(d_dim)*(cur_props%dc_bias_voltage-V)*dvdq))
!            !Fts_electrostatic_noncons = NonDimenForce( -(1/2)*(dcdq(d_dim)*(cur_props%dc_bias_voltage - V)**2 - 2*capac(d_dim)*(cur_props%dc_bias_voltage-V)*dvdq), M)
           
!          end function Fts_electrostatic_noncons

!          !this calculates dV/dt for RES1.  not actually a force!
!          real*8 function Fts_nc_es(d, dpr, V)
!            real*8, intent(in) :: d, dpr, V           
!            real*8 :: d_dim, d_dot
!            d_dim = DimenLength(d)
!            d_dot = DimenVelocity(dpr)
           
!            !hank: at Arvind's request I added a surface potential term to the other electrostatic models
!            !you may want to incorporate that in your models, but worry about getting it to work first
           
!            Fts_nc_es = dcdq(d_dim)*d_dot*cur_props%dc_bias_voltage/capac(d_dim) -(V/capac(d_dim))*(1/cur_props%resis+dcdq(d_dim)*d_dot)

!          end function Fts_nc_es

!          real*8 function capac(d_dim)
!            real*8, intent(in) :: d_dim
!            !this is a hack.  Ces_nc_cap1 is really not changing with each mode.
!            capac = cur_fC%Ces_nc_cap1/(cur_props%epsilon_0 * cur_props%g_d_dim + cur_props%epsilon_d*d_dim)
!          end function capac

!          real*8 function dcdq(d_dim)
!            real*8, intent(in) :: d_dim
!            dcdq = cur_fC%Ces_nc_dcdq1/(cur_props%epsilon_0 * cur_props%g_d_dim + cur_props%epsilon_d*d_dim)**2
!          end function dcdq

! end module NonConsElectrostatic


module Scanning
  use params
  implicit none
  save

  real*8 :: XPos, deltaX, Xpos_halt 
  real*8 :: HF, LF, LF2, SubsLen, Vscan
  real*8 :: LineSpeed_dim, Vscan_dim, SubsLen_dim
  real*8, private :: LS
  integer :: FeatType
  logical WantTSCON, WantFProp

  private TSCON

! !tip-sample convolution functions. 
   contains

     subroutine Init_Scanning()
       use TimeAndCycle, only : plotpnts
           XPos = 0d0
           deltaX = SubsLen/plotpnts
           XPos_halt = deltaX
           !definition of LS: substrate distance before and after feature           
           LS = (SubsLen - LF)/2

           if (LS < 0) call WriteFatalError(  "Substrate length must be longer than feature length.  Either increase subtrate length or decrease feature length.")

     end subroutine Init_Scanning

     subroutine computeXPos(t)
       use timeAndCycle, only : isTransientOver, timeSinceTrans
       real*8 , intent(in) :: t
       
       if ( isTransientOver ) then
          Xpos = Vscan * timeSinceTrans(t)
       else
          Xpos = 0
       end if
     end subroutine computeXPos

     subroutine computeXPosForceVol(t)
       use timeAndCycle, only : isTransientOver, timeSinceTrans
       use ForceVolume
       use Nondimensionalization
       use Approaching, only : ApproachStart
       real*8, intent(in) :: t             
       if (isTransientOver) then
          if (ForceVolState == FORCEVOL_SCAN) then
             XPos = VScan*timeSinceTrans( t_UptoScan - totalTimeApprRev )
             if ( Xpos > Xpos_halt ) then               
                ForceVolState = FORCEVOL_APPR
                t_UptoScan = t
                Xpos_halt = Xpos_halt + deltaX
                call ApproachStart()         
             end if
          end if
       else
          Xpos = 0
       end if        
     end subroutine computeXPosForceVol
     
    !thoughts for further improvement:
    !probably can use an adaptive step in theta, only need small steps around the edges.
!previous version got stuck in infinite loops.
    real*8 function TSCON()
      use params
      use contactModels, only : Rtip_nd
      real*8 :: theta, ht, h0
      integer, parameter :: N = 350
      integer :: i
      
      TSCON = 0d0
      
      Do i = -N, N
         theta = (pi / 2d0) * real(i) / real(N)   

         ht = ActualFeatHght(Theta)
         h0 = ActualFeatHght(0.0d0)
         
         if (ht > (h0+TSCON+Rtip_nd*(1-cos(Theta)))) then
            TSCON = ht - h0 - Rtip_nd*(1-cos(Theta))
         end if     
      end do
      
    END function TSCON
    
!this one used to be called Zs() in the old code    
    double precision function ApparentFeatHght()
      if (WantTSCON) then
         ApparentFeatHght = ActualFeatHght(0.0d0)+TSCON()
      else
         ApparentFeatHght = ActualFeatHght(0.0d0)
      end if
      
    END function ApparentFeatHght

    pure logical function OnSubstrate( X2 )
      real*8, intent(in) :: X2
      OnSubstrate = ((X2 <= LS).OR.(X2 >= (LS+LF)))
    end function OnSubstrate

!this one used to be called Zf in the old code
!  Theta = 0 is the actual sample height                
    double precision function ActualFeatHght(Angle)
      use contactModels, only : Rtip_nd
      real*8, intent(in) :: Angle
      
      real*8 :: X2, theta
      
      X2 = XPos + Rtip_nd*sin(Angle)
      
      if  (OnSubstrate(X2))  then         
         ActualFeatHght = 0.0d0
      else         
         if (FeatType .eq. STEP) then
            ActualFeatHght = HF
         else if (FeatType.eq. TRAPEZOID) then
            
            if ((X2-Ls).le.(LF-LF2)/2.0) then
               ActualFeatHght = (X2-Ls)*HF*2.0/(LF-LF2)
            else if ((X2-Ls).ge.(LF-(LF-LF2)/2.0)) then
               ActualFeatHght = HF-((X2-Ls)-(LF+LF2)/2.0)*HF*2.0/(LF-LF2)
            else
               ActualFeatHght = HF
            end if
            
         else if (FeatType.eq. SINUSOID) then
            ActualFeatHght = HF/2.0*(1.0+sin((X2-Ls)*2.0*pi/LF-pi/2))
         else if (FeatType .eq. CYLINDER) then
            theta = acos(  ((LF/2) - (X2-Ls)) / (LF/2)  ) 
            ActualFeatHght = (LF/2.0)*(1.0+sin(theta))
         else
            call assert( .false., 'unknown feature type')
            ActualFeatHght = 0.0d0
         end if
         
      end if
  
    End function ACTUALFEATHGHT
  end module Scanning



module timeHistory
  use params

  integer Amp_index, numHist 
  logical wanthist_byA, Wanthist
  real*8 Ahist(maxHist), Fhist(maxHist), xHist(maxHist)

  contains

    !this sorts Ahist
     subroutine SortTimeHistories(dir)
       logical, intent(in) :: dir

       real*8 tempHist
       integer i, i_A

       !	For retraction curves sort in reverse order
       tempHist = 0
       if ( dir ) then
          do i = 1, numHist
             do i_A = 1, numHist
                if ((Ahist(i_A) .lt. tempHist) .AND. (i_A > 1)) then
                   Ahist(i_A-1) = Ahist(i_A)
                   Ahist(i_A) = tempHist
                end if
                tempHist = Ahist(i_A)
             end do
          end do
       else
          do i = 1, numHist
             do i_A = 1, numHist
                if ((Ahist(i_A) .gt. tempHist) .AND. (i_A > 1)) then
                   Ahist(i_A-1)=Ahist(i_A)
                   Ahist(i_A) = tempHist
                end if
                tempHist = Ahist(i_A)
             end do
          end do
       end if

     end subroutine SortTimeHistories

      logical function startTimeHistory(want_Fourier, IOUT, output_point_rate, AprchS_set, Amp, MainLockinR, t, sweep_rate, XPos)
        use data1, only: computeZdist_dim
        use timeAndCycle, only: operating_mode, isTransientOver, pointsSinceTrans, omegad
        use Scanning, only: ComputeXPos
        logical, intent(in)    :: want_Fourier
        integer*8, intent(in)  :: IOUT, output_point_rate
        real*8, intent(in)     :: AprchS_set,  Amp, MainLockinR, t, sweep_rate, XPos

        real*8 Compamp

        
           if ( Amp_index > numHist) then
              !already did all of the time histories
              startTimeHistory = .false.
              return
           end if
           
           if (isOpModeApp(operating_mode)) then
              if (wantHist_byA) then !history by amplitude ratios
                 !(note, can't just check Amp == Ahist, b/c we might overshoot, and can't just do Amp<Ahist, 
                 !b/c we might be retracting.)
                 if (Want_Fourier) then
                    if ( INT(pointsSinceTrans(IOUT))< output_point_rate) then
                       !haven't computed first point yet.  There is no value in Amp
                       startTimeHistory = .false.
                       return
                    else
                       Compamp = sign(1d0,AprchS_set)*(Amp-Ahist(Amp_index))              
                    end if
                 else
                    Compamp = sign(1d0,AprchS_set)*(MainLockinR-Ahist(Amp_index))              
                 end if
                 
                 !steven: if the code could not find a time history corresponding to a certain Amplitude ratio, the code
                 !kept looking for that Amplitude ratio, even if other amplitude ratios were requested and actually occurred.  
                 !This change allows the code to search for the NEXT amplitude ratio if the current one was not found 
                 !(jumps can occur which would not allow a certain amplitude ratio to show up).  The Compamp .gt. -1.0d-2 is 
                 !used to find an amplitude ratio within a certain tolerance.
                 if (Compamp .lt. 0.0d0) then
                    if (Compamp .gt. -1.0d-2) then
                       startTimeHistory = .true.
                    else
                       startTimeHistory = .false.
                       Amp_index = Amp_index + 1
                    end if
                 else
                    startTimeHistory = .false.
                 end if
              else
                 startTimeHistory = sign(1d0,AprchS_set)*(computeZDist_dim(t) - Ahist(Amp_index) ) < 0d0
              end if
           else if (operating_mode == FREQSWEEP) then
              !raise a flag if we have passed the point where we want to start outputting time history
              startTimeHistory = sign(1d0,sweep_rate)* (Fhist(Amp_index) - omegad(t,1)) < 0d0
           else if (operating_mode == SCAN) then
              call ComputeXPos(t) 
              startTimeHistory = XPos > xHist(Amp_index)  
           else if (operating_mode == FIXED) then
              startTimeHistory = isTransientOver
           else
              call assert(.false., 'Unhandled time history case')
              startTimeHistory = .false.
           end if                 

      end function startTimeHistory

end module timeHistory


  module ForcingCalculations
    use params
    real*8 :: osc_Q, osc_omega_dim, osc_omega_nd
    real*8, private :: omegad_init, Fmag
    real*8 :: phase_slope_dim , efficiency_slope_dim, phase_slope_nd, efficiency_slope_nd
    logical want_nonideal_magnetic
  contains
    
    subroutine RescaleModalForcesFreqSweep(fexcite, nummodes, Abase, t, omegai_nd, Quality, Keq, F, phid, beta, mu, mtip,Afluid, mstar_div_m, exc_choice, Abase_init, want_nonideal_magnetic, B)
      use timeAndCycle, only : omegad
      integer, intent(in) :: fexcite, nummodes, exc_choice
      real*8, intent(in) :: t, omegai_nd(maxModes), Quality(maxModes), Keq(maxModes), beta(maxModes), mu(maxModes), mtip, mstar_div_m, Abase_init, B(maxModes)
      real*8, intent(inout) ::  Abase(maxExc)
      real*8, intent(out) :: F(maxExc,maxModes), phid(maxExc,maxModes)
      logical, intent(in) ::  want_nonideal_magnetic
      complex*8, intent(in) :: Afluid

      if ((fexcite==ACOUSTIC_IDEAL) .or. (fexcite == ACOUSTIC_PEAK)) then
         call CalcModalForcesIdealAcousticMode( numModes, Abase(1), omegad(t,1), omegai_nd, Quality,Keq, F(1,:), phid(1,:), beta,mu, mtip, Afluid, mstar_div_m)
         if (exc_choice == BIMODAL) call CalcModalForcesIdealAcousticMode( numModes, Abase(2), omegad(t,2), omegai_nd, Quality,Keq, F(2,:), phid(2,:), beta,mu, mtip, Afluid, mstar_div_m)
      elseif (fexcite == ACOUSTIC_PIEZOR) then
         call CalcModalForcesAcouticPiezoResMode( numModes, Abase(1), omegad(t,1), omegai_nd, Quality,Keq, F(1,:), phid(1,:), beta,mu, mtip, Abase_init, Afluid, mstar_div_m)
      elseif (isMagnetic(fexcite) .and. (want_nonideal_magnetic)) then
         call ScaleModalForcesMag(B, Keq, numModes, F(1,:), Abase(1), Phid, 1d0, omegad(t,1))
         !else don't need to rescale.
      end if
    end subroutine RescaleModalForcesFreqSweep
   
    subroutine CalcExcitAmp(output_type, omegai, fexcite, Quality, Chi, B, Keq,omegai_nd, &
         Ainitial_nd, F, Abase, phid,  Abase_init, F1_init, numModes, beta, mu, mtip, &
         want_Abase_direct, Abase_input, F_input, Afluid, mstar_div_m, fm_initial_phase, ZBase_Amp_nd, modulation_type)
      use params
      use timeAndCycle, only: operating_mode, omegad, exc_choice
      use NonDimensionalization

      implicit none
      integer, intent(in) :: output_type, fexcite, numModes, modulation_type
      real*8, intent(in) ::  omegai(maxModes), Quality(maxModes), Chi(maxModes), & 
         B(maxModes), Keq(maxModes), omegai_nd(maxModes), Ainitial_nd(maxModes), beta(maxModes),&
         Abase_input(maxExc), mu(maxModes), mtip, F_input(maxModes), mstar_div_m, Zbase_Amp_nd
      real*8, intent(out) :: F(maxExc, maxModes), Abase(maxExc), phid(maxExc, maxModes)
      real*8, intent(out) :: F1_init(maxModes), Abase_init, fm_initial_phase
      logical, intent(in) :: want_Abase_direct
      complex*8, intent(in) :: Afluid      

      integer :: i, j, numExc
      real *8 :: omegad_for_calc

      !defaults (sets all array elements)
      Abase = 0d0
      F = 0d0
      phid = 0d0

      if ( operating_mode == FREQSWEEP ) then
         omegad_init = omegai_nd(1)            
      else
         omegad_init = omegad(0d0, 1) !this is used for piezo resonance and non-ideal magnetic cases
      end if
            
      if (exc_choice == BIMODAL ) then
         numExc = 2
         if (isMagnetic(fexcite) .and. want_Abase_direct) call WriteFatalError("bimodal, magnetic, direct input is not supported. Either use acoustic drive, or input unconstrained amplitdes")         
      else
         numExc = 1
      end if

      do j = 1,numExc
         if ( operating_mode == FREQSWEEP ) then
            !initial drive freq is not a good choice for freq sweep
            omegad_for_calc = omegai_nd(1)            
         else
            omegad_for_calc = omegad(0d0, j)
         end if

         if (isMagnetic(fexcite)) then            
            ! Magnetic excitation            
            ! F(i,j) is the modal force on the jth mode due to ith excitation freq
            if ( .not. want_Abase_direct) then !fixme, this needs a better name
               call CalcFmag( omegad_for_calc, omegai, Quality, Chi, B, Keq, numModes, Ainitial_nd(j), fm_initial_phase )
               call ScaleModalForcesMag(B, Keq, numModes, F(j,:), Abase(j), Phid(j,:), 1d0, omegad_for_calc )
            else
               do i = 1,numModes
                  F(j,i) = NonDimenForce(F_input(i))
               end do
            end if
         elseif ((isAcoustic(fexcite)) .or. (fexcite == ACOUSTIC_PEAK)) then
            ! Acoustic or subsurface excitation		            
            Abase(j) = CalcBaseExcitation(output_type, omegad_for_calc, omegai,  Quality, Chi, numModes, Ainitial_nd(j), mtip, beta, mu, Afluid, mstar_div_m, fm_initial_phase)
            if (  want_Abase_direct) then
               !CalcBaseExcitation still need for fm_initial_phase, but throw away the Abase computation
               Abase = Abase_input
            end if
            
            if (modulation_type == PEAK_FORCE) then
               Abase(1) = Zbase_Amp_nd 
            end if
            
            if ((fexcite == ACOUSTIC_IDEAL) .or. (fexcite == ACOUSTIC_PEAK)) then
               call CalcModalForcesIdealAcousticMode( numModes, Abase(j), omegad_for_calc, omegai_nd, Quality, Keq, F(j,:), phid(j,:), beta, mu, mtip, Afluid, mstar_div_m)
            elseif (fexcite == ACOUSTIC_PIEZOR) then
               call CalcModalForcesAcouticPiezoResMode( numModes, Abase(j), omegad_for_calc, omegai_nd, Quality, Keq, F(j,:), phid(j,:), beta, mu, mtip, Abase(j), Afluid, mstar_div_m, fm_initial_phase)
            end if

         elseif (fexcite == NO_EXC) then
            !nothing to do
         else
            call WriteFatalError( 'unhandled excitation');
         end if

         !for self excitation overwrite any phid.  we'll implement any phase directly with feedback
         if (exc_choice == SELFEXC) phid = 0
                                 
      end do

       
      Abase_init = Abase(1) 
      F1_init = F(1,:)

    end subroutine CalcExcitAmp

!despite the name, this never really calculated photothermal per se.  it merely calculated a non-ideal excitation source
!where the excitation amplitude had some phase shift.  it's been absorbed into the non-ideal magnetic calculation.
!I started a proper calculation of photothermal in basic_model.f90, but never finished it.

    ! subroutine CalcInitialForcePhotothermal(omegad, omegai, Quality, Chi, B, Keq, numModes, Ainitial_nd,  F, Abase, Phid, fm_initial_phase)
    !   use Nondimensionalization
    !   use params
    !   implicit none
    !   integer, intent(in) :: numModes
    !   real*8, intent(in) ::  omegad, omegai(numModes), Quality(numModes),Chi(numModes), B(numModes), Keq(numModes), Ainitial_nd
    !   real*8, intent(out) :: F(numModes), Abase, Phid(numModes), fm_initial_phase
    !   real*8 Fpth, omegad_i(numModes)
    !   complex*8 Gtot
    !   integer i

    !   !fixme, need to calculate a proper B(i) for photothermal (assuming some spot position).
    !   !and then do this like we do magnetic mode.
    !   !for now we just use a hack and assume 100% of excitation goes to mode 1 and nothing to other modes

    !   omegad_i(1) = DimenFreq(omegad)/omegai(1)

    !   Gtot = Chi(1)*GAIN(omegad_i(1),Quality(1))/Keq(1)

    !   Fpth = DimenLength(Ainitial_nd)/CABS(Gtot)   !in a freq sweep these will be constant

    !   F(1) = NonDimenForce(Fpth)

    !   do i = 2,numModes
    !      Phid(i) = 0d0
    !      F(i) = 0d0
    !   end do

    !   Abase = 0d0

    !   fm_initial_phase = pi/2d0 + atan2(AIMAG(Gtot),REAL(Gtot))  

    ! end subroutine CalcInitialForcePhotothermal


    ! omegad   is normalized by omega_scale
    subroutine CalcFmag( omegad, omegai, Quality, Chi, B, Keq, numModes, Ainitial_nd, fm_initial_phase)     
      use params
      use Nondimensionalization
      integer, intent(in) :: numModes
      real*8, intent(in) :: omegad, omegai(numModes), Quality(numModes),Chi(numModes), B(numModes), Keq(numModes),  Ainitial_nd
      real*8, intent(out) :: fm_initial_phase
      real*8 omegad_i(numModes)
      complex*8 Gtot
      integer i

      !omegad is scaled by omega_scale in most of the program.  In this routine, we need it scaled by each
      !of the natural frequencies in turn.  So omegad_i(i) = omegad * omega_scale / omegai(i) is really
      !just the omegad_dimensional / omegai (dimensional). used to use variable omegad12 for that, but generalized to
      !arbitrary number of modes case

      Gtot = (0d0,0d0)
      do i = 1,numModes
         omegad_i(i) = DimenFreq(omegad)/omegai(i)
         Gtot = Gtot+Chi(i)*B(i)/Keq(i)*GAIN(omegad_i(i),Quality(i))
      end do

      Fmag = DimenLength(Ainitial_nd)/CABS(Gtot) !in a freq sweep these will be constant. make dim so can do a complete force non-dimen in one place below
      fm_initial_phase = pi/2d0 +  atan2(AIMAG(Gtot),REAL(Gtot))  

    end subroutine CalcFmag

    subroutine ScaleModalForcesMag(B, Keq, numModes, F, Abase, Phid, Drive_signal, omegad )
      use NonDimensionalization
      integer i
      integer, intent(in) :: numModes
      real*8, intent(in) :: B(numModes), Keq(numModes), Drive_signal, omegad
      real*8, intent(out) :: F(numModes), Abase, Phid(numModes)

      if (want_nonideal_magnetic) then
         do i = 1, numModes
            F(i) = NonDimenForce(Drive_signal*Fmag*B(i))  * ( 1 + efficiency_slope_nd *  (omegad_init - omegad ))
            Phid(i) = phase_slope_nd *  (omegad_init - omegad )
         end do
      else
         !ideal
         do i = 1,numModes
            Phid(i) = 0d0
            F(i) = NonDimenForce(Drive_signal*Fmag*B(i))
         end do
      end if

      Abase = 0d0
    end subroutine ScaleModalForcesMag


    !this is for acoustic mode. 
    !for low q and drive freq far away from natural frequency, it will look like this does not calculate the correct gain
    !but that is because you have to drive the base so hard off resonance that the relative deflection (what's measure) and the absolute deflection
    !may be very different and the tip & sample may start interacting at gaps way above what you think.
    real*8 function CalcBaseExcitation(output_type, omegad, omegai, Quality, Chi, numModes, Ainitial_nd, mtip, beta, mu, Afluid, mstar_div_m, fm_initial_phase)
      use Nondimensionalization
      use params
      implicit none
      integer, intent(in) :: output_type
      real*8, intent(in) :: omegad, omegai(numModes), Quality(numModes), Chi(numModes), Ainitial_nd, mtip, beta(numModes), mu(numModes), mstar_div_m
      complex*8, intent(in) :: Afluid
      real*8, intent(out) :: fm_initial_phase
      integer numModes
      real*8 omegad_i(numModes)
      complex*8 Gtot, tmp
      integer :: i

      Gtot = (0,0)
      do i = 1,numModes
         omegad_i(i) = DimenFreq(omegad) / omegai(i)
         Gtot=Gtot+Chi(i)*Gain_Aco_Num(omegad_i(i),Quality(i),beta(i),mu(i),mtip, Afluid, mstar_div_m)*GAIN(omegad_i(i),Quality(i))
      end do

      if (output_type == INTERFEROMETER) then
         if (CABS(Gtot) == -1d0) then
            call WriteFatalError('base motion is infinite.  check settings')
         end if

         tmp = Ainitial_nd / (1 + Gtot)
         CalcBaseExcitation = CABS(tmp)
         fm_initial_phase = atan2(AIMAG(tmp), REAL(tmp)) !fixme, double check this.
      else
         if (CABS(Gtot) == 0d0) then
            call WriteFatalError('base motion is infinite.  check fluid borne loading')
         end if

         CalcBaseExcitation  =  Ainitial_nd /CABS(Gtot) ! in a freq sweep, this should be constant
         fm_initial_phase = pi/2d0 +  atan2(AIMAG(Gtot),REAL(Gtot))  
      end if
    end function CalcBaseExcitation

    subroutine CalcModalForcesAcouticPiezoResMode(numModes, Abase, omegad, omegai_nd, Quality,Keq, F, phid, beta, mu, mtip, Abase_init, Afluid, mstar_div_m, fm_initial_phase)
      use params
      implicit none
      integer, intent(in) :: numModes
      real*8, intent(in) ::  omegad, omegai_nd(numModes), Quality(numModes),Keq(numModes), beta(numModes), mu(numModes), mtip, Abase_init, mstar_div_m
      complex*8, intent(in) :: Afluid
      real*8, intent(out) :: F(numModes), phid(numModes), Abase
      real*8, intent(inout), optional :: fm_initial_phase

      complex*8 :: g

      g = GAIN( omegad / osc_omega_nd, osc_Q) /  GAIN( omegad_init / osc_omega_nd, osc_Q) 

      !step 1, calculate new base amplitude based on drive frequency
      Abase = Abase_init * CABS(g)

      !step 2, use old routine for modal force and phase 
      call CalcModalForcesIdealAcousticMode(numModes, Abase, omegad, omegai_nd, Quality,Keq, F, phid, beta, mu, mtip, Afluid, mstar_div_m)

      !step 3, adjust phase
      phid = phid + atan2(AIMAG(g),REAL(g)) 
      if (present(fm_initial_phase)) then
         fm_initial_phase = fm_initial_phase + atan2(AIMAG(g),REAL(g)) 
      end if
    end subroutine CalcModalForcesAcouticPiezoResMode

    !for freq sweep code, this function needs to be called at every time step!
    subroutine CalcModalForcesIdealAcousticMode(numModes, Abase, omegad, omegai_nd, Quality, Keq, F, phid, beta, mu, mtip, Afluid, mstar_div_m)
      use NonDimensionalization
      implicit none
      integer, intent(in) :: numModes
      real*8, intent(in) :: Abase, omegad, omegai_nd(numModes), Quality(numModes), Keq(numModes), beta(numModes), mu(numModes), mtip, mstar_div_m
      complex*8, intent(in) :: Afluid
      real*8, intent(out) :: F(numModes), phid(numModes)
      complex*8 G
      real*8 omegad_i(numModes)
      integer i

      ! F is a (non-dimensional) modal force
      do i = 1,numModes
         omegad_i(i) = omegad/omegai_nd(i)
         G = Gain_Aco_Num(omegad_i(i),Quality(i),beta(i),mu(i),mtip, Afluid, mstar_div_m) * Keq(i)
         phid(i) = atan2(AIMAG(G),REAL(G))	
         F(i)= NonDimenForce(DimenLength(Abase)*CABS(G)) !in a freq sweep, this will not be constant.
      end do
     
    end subroutine CalcModalForcesIdealAcousticMode

    complex*8 function Gain(omega, Quality)
      real*8, intent(in):: omega, Quality

      Gain = 1/CMPLX(1-omega**2,omega/Quality)

    end function Gain

    ! numerator of the Acoustic mode transfer function
    ! also use to drive phases (phid)
    complex*8 function Gain_Aco_Num(omega,Quality,beta,mu,mtip, Afluid, mstar_div_m)

      real*8, intent(in):: omega, Quality, mu, beta, mtip
      real*8, intent(in) :: mstar_div_m !need a better name for this one.  ratio
      ! of added fluid mass to (cantilever mass + fluid mass).
      complex*8, intent(in) :: Afluid
      real*8 r1, r2, r3

      r1 = (beta+mtip)/(mu+mtip)
      r2 = beta/mu
      r3 = beta / (mu+mtip)

      !july 9, 2011.  use new convention for Afluid
      Gain_Aco_Num = omega**2*(r1+r3*mstar_div_m*Afluid) - CMPLX(0,1) * (omega/Quality*r2*(1+Afluid))

    end function Gain_Aco_Num

  end module ForcingCalculations

  module data2
    use params
  contains 
    subroutine runningForceCalcs(IOUT, IOUT_prevOutput, Force, FPeakAtt,FPeakRep, MeanForce, d, dpr, t, y)
      use contactModels, only : Fts

      integer*8, intent(in) :: IOUT, IOUT_prevOutput
      real*8, intent(inout) ::  FPeakAtt, FPeakRep, MeanForce
      real*8, intent(out) :: Force
      real*8, intent(in) :: d, dpr,t,y(:)
      integer :: i
      logical :: fail

      Force = Fts(d,dpr, fail, t, y)

      !	   Reset peak force parameters. this resets on the first step of a new outpoint point.  can't reset when we
      !          reset meanForce b/c we don't want to reset it to zero.
      if ((IOUT - IOUT_prevOutput) == 1) then
         FPeakAtt = Force
         FPeakRep = Force
      end if
      
      MeanForce=MeanForce+Force
      
      IF (Force .GT. FPeakRep) THEN
         FPeakRep=Force
      END IF
      
      IF (Force .LT. FPeakAtt) THEN
         FPeakAtt=Force
      END IF     
    end subroutine runningForceCalcs
    
    subroutine ResetCycleVariables(num_impacts, a1, b1, a2, b2, an, bn, F1r, F1i, Indent, MeanForce, Ets, Virial, Ebs, Eprop, Edamp,Edrive, Force_old, dpr_old, tc, min_d, numHH,numModes )
      integer, intent(in) :: numHH, numModes
      real*8, intent(out) :: a1, b1, a2(1), b2(1), an(numHH), bn(numHH), Indent, MeanForce, &
           Ets, Ebs, tc, Eprop(numModes), min_d, Edamp(numModes), Edrive(numModes), Force_old, &
           F1r(1), F1i(1), dpr_old, Virial
      integer, intent(out) :: num_impacts
      integer :: i

      tc = 0      
      Force_old = 0
      dpr_old = 0
      num_impacts = 0
      min_d = 1d99
      Ets = 0d0
      Ebs = 0d0
      Virial = 0d0
      MeanForce=0d0
      Indent = 0
      a1 = 0
      b1 = 0
      a2(1) = 0
      b2(1) = 0
      F1r(1) = 0
      F1i(1) = 0
      do i=1,numHH
         an(i) = 0
         bn(i) = 0
      end do
      
      do i=1,numModes
         Edrive(i) = 0
         Edamp(i) = 0
         Eprop(i) = 0
      end do
    end subroutine ResetCycleVariables

! Default Z range calculation
! Daniel's comments: the Z0 value is the value of the base of the Z piezo.In general, it is not the tip sample gap and it 
!does not include any static deflection due to long range forces.  Now for linear or hertz contact, there are no
!long range forces, so there is no static deflection. but with van der waals (if you're very close) or DLVO, there
!can be significant static deflection of the tip, which is not included in here. this is why I added the approach-to-setpoint
!option for freq sweeps.
    subroutine CalcZRange(operating_mode,modulation_type, ZRange, want_freqswp_sp_aprch, Z0,  Zf, y, Zbase_Amp_nd)
      use params
      use Nondimensionalization
      use ContactModels

      integer, intent(in) :: operating_mode, modulation_type, Zrange
      logical, intent(in) :: want_freqswp_sp_aprch
      real*8, intent(inout) :: Z0, Zf
      real*8, intent(in) :: y(:), Zbase_Amp_nd
      logical :: fail
      
      if ( (isOpModeApp(operating_mode) .and. (Zrange == ZRANGE_AUTO )) .or. &
           (( operating_mode == FREQSWEEP) .and. ( want_freqswp_sp_aprch))) then
!fixme, we need to do something intelligent for electrostatic models here.  
         if (cur_props%fts_model .ne. DMT_DLVO) then
            Z0 = 1d0+NonDimenLength(5d-9)
            Zf = 0d0
         else
            !Adjusts for the DLVO bias			
            !original method was starting too far neg for de-i water & 0.6N/m lever.  try this instead
            !may have a longer transient but hopefully it's more stable
            Z0 = 1d0-Fts(1d0,0d0, fail, 0d0,y)/4+NonDimenLength(5d-9)
            !Note: Fts is already normalized by Keq(1)		
            Zf = Z0-NonDimenLength(5d-9)-1d0
         end if
      else if (operating_mode == SCAN) then
         if (modulation_type == AMPLITUDE) then
            Z0 = 1d0 + cur_props%aDMT
         elseif (modulation_type == PEAK_FORCE) then
            Z0 =  cur_props%aDMT/2d0 + Zbase_Amp_nd           
         end if

         !at one point I had tried to guess a value of Z0 that would get closer (thus reducing
         !transient times, but we can't do that.  There may be a bistable region and if we start at the wrong place we may end
         !up on the wrong branch.  In the experiment you always start far away and approach to the setpoint so that's what we
         !do here.  note: these means we need longer transient times to let the controller find the right distance.
         Zf = NonDimenLength(Zf)
      else
         !just normalize the input values
         Zf = NonDimenLength(Zf)
         Z0 = NonDimenLength(Z0)
      end if
      
      if ( Z0 == cur_props%aDMT) then
         !for some reason, DDASKR does not like the initial conditions right here.  just fudge it a little
         Z0 = Z0 + 1e-6
      end if

    end subroutine CalcZRange
    

  end module data2

  module HydrodynamicFunction
    ! this module intended to scale damping and natural frequency according to
    ! hydrodynamic added mass and viscosity. i.e., the user still enters a 
    ! natural frequency and quality factor calibrated from experiment, but then
    ! this module predicts how that would change with a changing drive frequency (e.g. for FM)

    !fixme: I think there may be a bug in this function.  The 1st harmonic amplitude
    !does not match my matlab code on a freq sweep (near zero freq).  I think it might be because
    !of ScaledForce().  That is supposed to scale out the modal mass.  But with hydro
    !dynamic function, the modal mass is changing with drive freq.  Has that 
    !been correctly accounted for?
    use params
    use Nondimensionalization
   
    real*8, private :: rho, b, L, eta, damping0(maxModes), omegai0(maxModes), keq(maxModes)
    real*8, private :: added_viscosity0, added_mass0,  vacuum_mass
   
  contains
     
    !the fluid parameters are all dimensional
    !the viscosity and added mass is assumed to have been calculated near the wet natural frequency
    subroutine InitHydrodynamicFunction( damping_in, omegai_in, keq_in, rho_in, b_in, L_in, eta_in)
      real*8, intent(in) :: damping_in(maxModes), omegai_in(maxModes), keq_in(maxModes), rho_in, b_in, eta_in, L_in
      real*8 :: total_mass0, modal_mass0
      damping0 = damping_in
      omegai0 = omegai_in
      keq = keq_in
      rho= rho_in
      b = b_in
      L = L_in

      eta = eta_in
      added_viscosity0 = added_viscosity(DimenFreq(omegai0(1)))
      added_mass0 = added_mass(DimenFreq(omegai0(1)))

      modal_mass0 = keq(1) / (DimenFreq(omegai0(1)))**2
      total_mass0 = modal_mass0 *4 !fixme  1/alpha instead of 4
      vacuum_mass = total_mass0 - added_mass0
!      vacuum_mass = 2.5877d-12 * 4d0 !physical, not modal
      if ( vacuum_mass < 0) then
         call WriteFatalError('Error in hydrodynamic function calculation.  Entered values for cantilever modal parameters (k,Q, omega) are not consistent with physical parameters (b,L, rho, eta, etc).  Calcualted vacuum mass < 0.  check parameters')
      end if
    end subroutine InitHydrodynamicFunction

    real*8 function hydro_damping( omegad_dim, i)
      real*8, intent(in) :: omegad_dim
      integer, intent(in) :: i
      hydro_damping = damping0(i) * added_viscosity(omegad_dim) / added_viscosity0
!      hydro_damping = added_viscosity(omegad_dim) / ( vacuum_mass+ added_mass(omegad_dim)) / omega_scale
      call assert(hydro_damping > 0, 'hydrodamping < 0')
    end function hydro_damping

    real*8 function hydro_omega( omegad_dim, i)
      real*8, intent(in) :: omegad_dim
      integer, intent(in) :: i
      real*8 :: total_mass, modal_mass

      if (i == 1) then
         total_mass = vacuum_mass +  added_mass(omegad_dim)
         modal_mass = total_mass / 4 !fixme, alpha instead of 4
         hydro_omega = NonDimenFreq(sqrt( keq(1) / ( modal_mass )))
      else
         call WriteFatalError('hydro dynamic function not implemented for higher eigenmodes yet')
      end if
      call assert(hydro_omega > 0, 'hydrodomega < 0')
    end function hydro_omega

    real*8 function  added_mass(omegad_dim)
      real*8, intent(in) :: omegad_dim
      complex*8 :: g
      g = GammaRect(omegad_dim)
      added_mass = (pi/4) * rho * b**2 * L * real( g )      
      call assert( added_mass > 0, 'added mass < 0')
    end function added_mass

    real*8  function  added_viscosity(omegad_dim)
      real*8, intent(in) :: omegad_dim
      complex*8 :: g
      g =  GammaRect(omegad_dim) 
      added_viscosity =  (pi/4) * rho * b**2 * L * omegad_dim * AIMAG(g) 
      call assert( added_viscosity > 0, 'added visc < 0')
    end function added_viscosity

    complex*8 function GammaRect( omegad_dim)
      real*8, intent(in) :: omegad_dim
      GammaRect = GammaRectRe( rho * omegad_dim * b**2 / ( 4 * eta) )
      call assert( real(GammaRect) > 0, 're gammarect < 0')
      call assert( AIMAG(GammaRect) > 0, 'im gammarect < 0')
    end function GammaRect
    
   !ref sader, j. app. physic. 98    
   complex*8 function GammaRectRe(Re)
     real*8, intent(in) :: Re
     real*8 :: omega_r, omega_i, t
     complex*8 :: c

     call assert(Re > 0, 're < 0')

     t = log10(Re)

     omega_r = (0.913242 - 0.48274*t + 0.46842*t**2 - 0.12886*t**3+0.044055*t**4 - 0.0035117*t**5 + 0.00069085*t**6) * (1d0 / (1d0 - 0.56964*t + 0.48690 *t**2 - 0.13444 * t**3 + 0.045155 * t**4 - 0.0035862 * t**5 + 0.00069085*t**6))

     omega_i =   (-0.0241342 - 0.029256*t + 0.016294 *t**2 - 0.00010961*t**3 + 0.000064577*t**4 - 0.000044510*t**5) * (1/(1 - 0.59702*t + 0.55182*t**2 - 0.18357*t**3 +  0.079156*t**4 - 0.014369*t**5 + 0.0028361*t**6))

     c = cylinder_hydro( Re)
    GammaRectRe = CMPLX(omega_r, omega_i) * c

    call assert( real(GammaRectRe) > 0, 're gammarectre < 0')
    call assert( AIMAG(GammaRectRe) > 0, 'im gammarectre < 0')
  end function GammaRectRe

  complex*8 function cylinder_hydro(Re)
    real*8, intent(in) :: Re
    real*8 :: x, sr, si
!the definition of the hydrodynamic function is in terms of modified bessel functions.
!but those are really complicated and a simple library is not easily available.  so just
!do a polynomal fit in matlab for now 

!    complex*8 :: den, num,b0, s
!    complex*8, parameter :: i = CMPLX(0,1)
!    den =  sqrt(CMPLX(0, Re)) * BesselK(0, -i * sqrt(CMPLX(0,Re)))
!    num = 4 * i * besselk(1, -i * sqrt(CMPLX(0,Re)))
!    cylinder_hydro = 1d0  + num / den
    if (( Re < 1e-2) .or. (Re > 100)) then
       call WriteFatalError('hydrodynamic function fit not performed to this reynolds number')       
    end if
    
    x = log10(Re)
    sr = -2.386981e-05 * x**6  -2.735704e-04 * x**5  -9.839586e-04 * x**4 + 3.004040e-03 * x**3 + 7.688232e-02 * x**2  -3.979671e-01 * x + 5.985609e-01
    si = 8.640225e-05 * x**6  -1.037941e-05 * x**5  -1.893417e-03 * x**4  -8.419601e-04 * x**3 + 4.837697e-02 * x**2  -6.702764e-01 * x + 6.595996e-01

    cylinder_hydro = CMPLX(10 ** sr, 10 ** si)

    call assert( real(cylinder_hydro) > 0, 're cylinderhydro < 0')
    call assert( AIMAG(cylinder_hydro) > 0, 'im cylinderhydro < 0')
  end function cylinder_hydro


 end module HydrodynamicFunction

 module SelfExcitation
   use params
   use CircularBuffer
  real*8 :: self_exc_gain, self_exc_saturation, self_exc_phase
  integer :: delay
  contains
  subroutine InitSelfExc( omegai_nd, numincycle, fexcite)
    real*8, intent(in) :: omegai_nd
    integer, intent(in) :: numincycle, fexcite
    real*8 :: self_exc_phase_rad

    if (isAcoustic(fexcite)) then
       call WriteFatalError("self excitation does not work for acoustic drive right now")
    end if
    
    delay = floor( numincycle *  (self_exc_phase/360d0) / omegai_nd)
    call InitCircularBuffer(delay)
  end subroutine  InitSelfExc
    
  function CalcModalForcesSelfExc(u, F1_init)
    real*8 CalcModalForcesSelfExc(maxExc, maxModes)
    real*8, intent(in) :: u, F1_init(maxModes)

    if (delay == 0) then
       CalcModalForcesSelfExc(1,:) = F1_init *  min( max( ( u * self_exc_gain ), -self_exc_saturation), self_exc_saturation)       
    else
       CalcModalForcesSelfExc(1,:) = F1_init *  min( max( ( ReadCircularBuffer() * self_exc_gain ), -self_exc_saturation), self_exc_saturation)       
       call WriteCircularBuffer(u)
    end if

    CalcModalForcesSelfExc(2,:) = 0d0
  end function CalcModalForcesSelfExc
end module SelfExcitation
