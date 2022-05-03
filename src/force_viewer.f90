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
!If VEDA (modified or unmodified), any portion thereof, or any derivative work (publicly distributed or not), contributes to any scientific or academic publication (including, but not limited to, journal articles, conference presentations or posters, or seminars) you must cite the original developers:
!
!Kiracofe, D.; Melcher, J. & Raman, A. "Gaining insight into the physics of dynamic atomic force microscopy in complex environments using the VEDA simulator Review of Scientific Instruments", (2012), 83, 013702
!
!You may also optionally cite the original VEDA article
!
! J. Melcher, S. Hu, A. Raman, "VEDA: A web-based virtual environment for dynamic atomic force microscopy" Review of Scientific Instruments 79, 061301 (2008).
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!main file for forceViewer tool
!-----------------------------------------------------------------------
    
      
program DMT_DDASKR
  use rappture_io
  use data0
  use data1
  use data2
  use contactModels
  use AppControl
  use Approaching
  use NonDimensionalization
  use TimeAndCycle  
  use Viscoelasticity_ting
  use Viscoelasticity_attard
  use tip_data
  use HertzModule
  use Integrals
  use Derivatives
  use MovingAverageFilter
  
  IMPLICIT none

  character(len=*), parameter :: fmtStr = '(E16.7E3,E16.7E3,A)'

  character*100 inFile, strVal
  character*200 tmpStr

  integer vel_model, n1, mov_avg_filt_len

  integer*8 :: nstep, longzero

  real*8  omegad_dim(maxModes), delZ, tf, tf_dim, TotalTravel, t, gap2, dpr2
  real*8, allocatable ::  Force(:),  gap(:),  Gradient(:), Y(:),  YPRIME(:), dpr(:), YP2(:), YP3(:), YP4(:), q(:)

  real*8 virial, dissipation, center

  logical :: fail
  character*3000 inputEcho

  integer, parameter :: NONE = 1, TRIANGLE = 2, SINE = 3, ROUNDED_TRIANGLE=4
  integer, parameter :: FWD_EULER = 1, MOD_EULER = 2, RK4 = 3

  integer, parameter :: ATTARD_SOLVER = RK4
  
  call getarg(1,inFile)

  call OpenInputFile(inFile)  

  InputEcho  = ""

  call ReadForceViewerParameters( Z0, Zf, nstep, vel_model, tf_dim, mov_avg_filt_len)

  call ReadTipData( Rtip_dim, Etip, Poisson_tip, tip_angle, tip_shape,  InputEcho)  
  call ReadSampleData(  substrate_props, "input.phase(ts)",  InputEcho)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!no more input to read after this point

        if (( vel_model == TRIANGLE) .or. ( vel_model == ROUNDED_TRIANGLE)) then
           nstep = 2d0 * floor( nstep/2d0) !force even so we have v=0 exactly at a solved point.
        end if     

!	omegai(1) = 1d0 * (2.0*pi)*1.0d3 ! Have to have this for normalization later since
        omegai(1) = 10000d0 * 2. * pi
        !we copied stuff from main program.  But doesn't really matter.
        omegad_dim(1) = omegai(1)
        
        Ainitial_dim(1) = 10d-9 !intentionally make these not 1 in order to catch errors
          
        numModes = 1

	omegai_nd(1) = 1.0d0 !	Frequency ratios	

	Chi(1) = 1
            
        Keq(1) = 1.6 !intentionally make these not 1 in order to catch errors

        if (.true.) then
           call set_omega_scale( omegad_dim(1) )
           call set_a_scale( Ainitial_dim(1) )
        else
           call set_omega_scale(1d0)
           call set_a_scale( 1d0)
        end if

        tf = NonDimenTime(tf_dim)
        delt = tf / dble(nstep)

!	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!	Normalization
!	Spatial: L(ND) = L(m)/Ainitial(m)
!	Temporal: t(ND) = omegai(rad/s)*t(sec)
!	Prior to normalization, all parameters are in SI units MKS (meters-kilograms-seconds)
!	Forces: F (ND) = F (Newton) * omegai_nd(i)**2 /Ainitial/Keq(i)
!	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!	

        Rtip_nd = NonDimenLength(Rtip_dim)

        call NormalizeMatlProps( substrate_props, Rtip_dim)
       
! done normalizing things
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

        subs_fc = computeForceCoeff(substrate_props)
 
        cur_props => substrate_props	     	
        cur_fC => subs_fC

!write(*,*) cur_fC%C_mag_dipole

        Zf = NonDimenLength(Zf)
        Z0 = NonDimenLength(Z0)
! ----------------------------------------------------------------------
!  Output plots
! ----------------------------------------------------------------------
  strVal =  "Tip-sample interaction force vs gap"
  call SetupGenericPlot( "FtsGap", strVal, strVal,"tip-sample gap", "nm", "Force", "nN")   

  strVal =  "Tip-sample interaction stiffness vs gap (positive stiffness = positive freq. shift)"
  call SetupGenericPlot( "FtsGrad", strVal, strVal,"tip-sample gap", "nm", "Interaction Stiffess", "nN/nm")   


  if (( vel_model == SINE) .or. (vel_model == TRIANGLE) .or.  (vel_model == ROUNDED_TRIANGLE) ) then  
     strVal =  "Tip-sample relative velocity vs gap"
     call SetupGenericPlot( "ddot", strVal, strVal,"tip-sample gap", "nm", "tip-sample velocity", "um/s") 

     strVal =  "Tip-sample interaction force vs time"
     call SetupGenericPlot( "FtsTime", strVal, strVal,"time", "s", "Force", "nN")   

     strVal =  "Tip-sample gap versus time"
     call SetupGenericPlot( "d", strVal, strVal,"time", "s", "tip-sample gap", "nm")

     strVal =  "Tip-sample velocity versus time"
     call SetupGenericPlot( "dpr", strVal, strVal,"time", "s", "tip-sample gap", "nm") 
  end if


  if (( cur_props%fts_model == HERTZ) .and. (cur_props%VEchoice == VEC_NONE)) then
            call SetupGenericPlot( "Stress", "Maximum sample stress (von Mises)", "Maximum sample stress (von Mises)", "Tip-sample gap", "nm", "von Mises Stress", "GPa", .false.)
  end if

! ----------------------------------------------------------------------
!  Initial conditions 
! ----------------------------------------------------------------------     
  n1 = numModes*2+1  ! only used for attard but always calculated just in case
  if (cur_props%fts_ve_model == VEM_ATTARD_BASED) then           
     if ( cur_props%N_Attard_fourier==0) then
        NEQ = numModes*2 + cur_props%N_Attard_spatial
     else
        NEQ = numModes*2 + cur_props%N_Attard_fourier
     end if
  else
     NEQ = 2*numModes !just to prevent seg faults in places where this is read.  not actually used
  end if
  
        
        allocate(YPRIME(NEQ))
        allocate(YP2(NEQ)) !for attard only
        allocate(YP3(NEQ)) !for attard only
        allocate(YP4(NEQ)) !for attard only
        allocate(Y(NEQ))
  	
      call InitializeHystereticModels(Z0)
      call CheckViscoelasticModels(cur_props,delt )

      IOUT = 0
      
      if ( vel_model == NONE) then
         TotalTravel =  (Zf - Z0)
      else
         TotalTravel =  2 * (Zf - Z0)
      end if

      delZ = (1d0 / real(nstep)) * TotalTravel


      Y = 0
      YPRIME = 0

      call assert(nstep > 1, 'no points?')
      
      allocate( Force(nstep))
      allocate( Gradient(nstep))
      allocate( gap(nstep))
      allocate( q(nstep))
      allocate( dpr(nstep))


      if ( cur_props%fts_ve_model == VEM_ATTARD_BASED ) call start_movie(driver, "                ", "1", "Time (s)")
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!1!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!start actual computation loop

      !first compute all the gaps and times
      do IOUT = 1, nstep

         t = (IOUT - 1) * delt

         if ( vel_model == NONE) then
            gap(IOUT) = Z0 + real(IOUT) * delZ
         elseif (( vel_model == TRIANGLE) .or. ( vel_model == ROUNDED_TRIANGLE)) then
            if (IOUT < floor(nstep/2d0)) then
               gap(IOUT) = Z0 + real(IOUT) * delZ
            else
               gap(IOUT) = Z0 + real(nstep)/2d0 * delZ - real(IOUT-nstep/2) * delZ
            end if
         elseif ( vel_model == SINE) then
            gap(IOUT) = Z0 + ((Zf-Z0)/2) - ((Zf-Z0)/2) * cos( 2d0 * pi * real(IOUT-1) / real(nstep))
         else
            call assert(.false., 'unhandled case')
         end if

         if ( vel_model == TRIANGLE) then
            if (IOUT < floor(nstep/2d0) ) then
               dpr(IOUT) = 2 * (Zf - Z0) / tf
            elseif (IOUT == floor(nstep/2d0) ) then
               dpr(IOUT) = 0  !previously lacked this zero. adding it makes attard more well behaved
            else
               dpr(IOUT) = -2 * (Zf - Z0) / tf
            end if
         elseif ( vel_model == SINE) then
            dpr(IOUT) = 2d0 * pi * ((Zf-Z0)/2) * sin( 2d0 * pi * t / tf) / tf
         elseif ( vel_model == NONE) then
            dpr(IOUT) = 0d0
         end if
         !rounded triangle is handled outside the loop         
      end do

      if ( vel_model == ROUNDED_TRIANGLE) then
         !moving average filter here
         !this is not perfect, still a bit of a discontinuity in velocity.  need to do an IIR, or maybe
         !compute the velocity differently but it will work for now. 
         gap = doMovingAverageFilter( gap, size(gap) , mov_avg_filt_len)
         
         dpr = centered_diff_4th( gap , size(gap) , delt)
      end if
      
      
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!1!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      
      Force(1) = 0d0

      !second, compute forces.
      do IOUT = 2, nstep

         !fixme, need better control over frame rate
         if ((mod(IOUT,10)==0) .and. ( cur_props%fts_ve_model == VEM_ATTARD_BASED)) call output_movie_frame(Y, driver, gap(IOUT), IOUT, cur_props, 1, DimenTime(t) )

         
         t = (IOUT - 1) * delt

         call updateHystereticModels_insideRES1(gap(IOUT), dpr(IOUT)) 


         !integrate equations to next time step
         if ( cur_props%fts_ve_model == VEM_ATTARD_BASED )  then
            if (ATTARD_SOLVER == FWD_EULER) then
               !a bad hack to compute attard model using forward Euler for force viewer
               !would like to do this the way I did amp_red_avg, but would take a little
               !work to compute gap as function of time for the rounded rectangle profile. would
               !have to re-write
                YPRIME( n1:NEQ) = compute_attard_derivative(gap(IOUT-1), dpr(IOUT-1),t, Y(n1:NEQ), cur_props)
                Y( n1:NEQ) = Y( n1 :NEQ) + YPRIME( n1:NEQ) * delt
            elseif (ATTARD_SOLVER == MOD_EULER) then
                !modified euler
                YPRIME( n1:NEQ) = compute_attard_derivative(gap(IOUT-1), dpr(IOUT-1),t,      Y(n1:NEQ),                         cur_props)
                YP2( n1:NEQ)    = compute_attard_derivative(gap(IOUT),   dpr(IOUT),  t+delt, Y(n1:NEQ)+ YPRIME( n1:NEQ) * delt, cur_props)
                Y( n1:NEQ) = Y( n1 :NEQ) + (YPRIME(n1:NEQ)+ YP2(n1:NEQ))*delt/2d0
             else
                !RK 4
                !fixme, technically this just interpolates the sine wave, when we really should calculate it EXACTLY
               gap2 = (  gap(IOUT-1) + gap(IOUT) ) / 2d0
               dpr2 = (  dpr(IOUT-1) + dpr(IOUT) ) / 2d0

               YPRIME( n1:NEQ) = compute_attard_derivative(gap(IOUT-1), dpr(IOUT-1),t,            Y(n1:NEQ), cur_props)
               YP2( n1:NEQ)    = compute_attard_derivative(gap2,        dpr2,       t+0.5d0*delt, Y(n1:NEQ)+ 0.5d0 * YPRIME( n1:NEQ) * delt, cur_props)
               YP3( n1:NEQ)    = compute_attard_derivative(gap2,        dpr2,       t+0.5d0*delt, Y(n1:NEQ)+ 0.5d0 * YP2( n1:NEQ) * delt, cur_props)
               YP4( n1:NEQ)    = compute_attard_derivative(gap(IOUT),   dpr(IOUT),  t+delt,       Y(n1:NEQ)+ YP3( n1:NEQ) * delt, cur_props)
               
               Y( n1:NEQ) = Y( n1 :NEQ) + (YPRIME(n1:NEQ)+ 2d0* YP2(n1:NEQ)+2d0* YP3(n1:NEQ) +  YP4(n1:NEQ))*delt/6d0
            end if
         end if !attard

         !I've done this very confusingly.  in this file, I've implemented the ODE solver to start from
         !(gap(IOUT-1),   dpr(IOUT-1), etc, to compute the Y at IOUT.
         !thus force must be compute after the ODE solver.
         !in the amplitude reduction averaging program, I've done it starting from IOUT to compute IOUT+1
         !thus force must be computed BEOFRE the ODE solver.
         !probably should make them consistent
         Force(IOUT) = Fts(gap(IOUT), dpr(IOUT), fail, t, y)
         if (fail) exit

         
         write(tmpStr,fmtStr) DimenLength(gap(IOUT))*1d9, DimenForce(Force(IOUT))*1d9, char(10)
         call rp_lib_put_str(driver, "output.curve(FtsGap).component.xy", tmpStr, 1)
         
         if (( vel_model == SINE) .or. (vel_model == TRIANGLE)  .or. (vel_model == ROUNDED_TRIANGLE) ) then  
            write(tmpStr,fmtStr) DimenTime(t), DimenForce(Force(IOUT)*1d9), char(10)
            call rp_lib_put_str(driver, "output.curve(FtsTime).component.xy", tmpStr, 1)
            write(tmpStr,fmtStr) DimenLength(gap(IOUT))*1d9, DimenVelocity(dpr(IOUT))*1e6, char(10)
            call rp_lib_put_str(driver, "output.curve(ddot).component.xy", tmpStr, 1)
            write(tmpStr,fmtStr) DimenTime(t), DimenLength(gap(IOUT))*1d9, char(10)
            call rp_lib_put_str(driver, "output.curve(d).component.xy", tmpStr, 1)

            write(tmpStr,fmtStr) DimenTime(t), DimenLength(dpr(IOUT))*1d9, char(10)
            call rp_lib_put_str(driver, "output.curve(dpr).component.xy", tmpStr, 1)
         end if

         
         if (( cur_props%fts_model == HERTZ) .and. (cur_props%VEchoice == VEC_NONE)) then           
            write(tmpStr,fmtStr)  DimenLength(gap(IOUT))*1d9, DimenPressure(Hertz_vonMises(gap(IOUT), Force(IOUT), cur_props%Poisson_sample)/1e9) , char(10)
            call rp_lib_put_str(driver, "output.curve(Stress).component.xy", tmpStr, 1)
         end if

         call updateHystereticModels_outsideRES1(gap(IOUT), dpr(IOUT),t, IOUT, Y(n1:NEQ) ) 
         call Update_viscoelasticity_outsideRes1_Hertz( t,  gap(IOUT) - cur_props%aDMT, dpr(IOUT), cur_props, Rtip_nd) 
         call Update_viscoelasticity_outsideRes1_linear( t, gap(IOUT) - cur_props%aDMT, cur_props) 

         if (mod(IOUT, 100_8) == 1)  then
            write(6,*) "=RAPPTURE-PROGRESS=>", idint(100d0 * real(IOUT)/real(nstep)) ," Performing simulation..."
            flush(6)
         end if

      end do !iout

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!1!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!      
      
      !finally output results
      Gradient(1) = (Force(2) - Force(1) ) / (gap(2) - gap(1))
      do IOUT = 2,nstep-1
         if ((gap(IOUT+1)-gap(IOUT-1)) == 0) then
            Gradient(IOUT) = (Force(IOUT) - Force(IOUT-1) ) / ((gap(IOUT)-gap(IOUT-1)))
         else
            !2nd order central difference.
            Gradient(IOUT) = (Force(IOUT+1) - Force(IOUT-1) ) / ((gap(IOUT+1)-gap(IOUT-1)))
         end if
      end do
      Gradient(nstep) = ( Force(nstep) - Force(nstep-1) ) / (gap(nstep)-gap(nstep-1)) 
      
      do IOUT = 1,nstep
         !write kint = minus gradient.  sign convention is that a positive interaction force
         !yields a positive frequency shift.
         write(tmpStr,fmtStr) DimenLength(gap(IOUT))*1d9, -DimenStiffness(Gradient(IOUT)), char(10)
         call rp_lib_put_str(driver, "output.curve(FtsGrad).component.xy", tmpStr, 1)
      end do
          
     
      !see also comments in the DMT_ddaskr.f90 file. you would be tempted to use this definition
      !of the virial, but it contains and offset due to the Z position.  What we really
      !want to do is multiply by q, the eigenmode deflection.  We don't have that here
      !but assuming a perfect sine wave, we get calculate what it would be 
      !Virial = DimenEnergy( sum( Force * gap) * DimenTime(delt) ) * electrons_per_coulomb

      center = (Zf + Z0) / 2
      q = gap - center

      !note there is no factor of 2 pi here.  that was here earlier but it was wrong! (9/7/2014 drk)
      Virial = - DimenEnergy( trapzu(  Force * q, int(nstep) , delt)  ) * electrons_per_coulomb / tf



      if (VEL_MODEL /= NONE) then
         Dissipation = - DimenEnergy( trapzu(  Force * dpr, int(nstep) , delt)  ) * electrons_per_coulomb
         call rp_lib_put_str(driver, "output.string(VirialDiss).about.label","Virial and dissipation", 0) 

         write(tmpStr,  *) CHAR(10), "Virial: ", Virial, " (eV/cycle) ", CHAR(10), 'Dissipation: ', Dissipation, ' (eV/cycle)', char(10) !char(10) is newline
         call rp_lib_put_str(driver, "output.string(VirialDiss).current", tmpStr, 0) 
      else
         Dissipation = 0d0
         call rp_lib_put_str(driver, "output.string(VirialDiss).current", "Virial and dissipation are not calculated because velocity option has been set to 'None'.  To calculate virial and dissipation, select a 'non-conservative velocity option' on the operating conditions tab.", 0) 
      end if      



      call OutputInputEcho( inputEcho)

      longzero = 0 !keep gfortran happy
      call OutputMiscParameters( numModes, F, phid, Keq(1), longzero, numincycle, cur_props%aDMT, Chi, omegai, B, 0d0, Z0, Abase, Y, beta, mu, cur_props%Estar, damping, nstep, 0d0, 0d0, 0d0, 0d0, 0d0,0d0, alpha, longzero, delt, longzero, 0, minimum_tau(cur_props), 0d0)


 call rp_result(driver)
! if ( cur_props%fts_model == ATTARD) call stop_movie
 STOP
END program DMT_DDASKR
      
! End of main function



!added for debugging.  don't want to call in production code.
subroutine assert( x, msg)
  logical :: x
  character (len=*) :: msg
 
  if (x .eqv. .FALSE. ) then
    write(6,*) "assertion failed ", msg
 end if
end subroutine assert
