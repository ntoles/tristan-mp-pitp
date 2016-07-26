!
! User module
!
! This module contains functions that may be altered by the user to define a specific problem. 
! 
! These functions are called by tristanmainloop.F90 and should provide
! initialization of parameters specific to the problem being setup, 
! loading of initial plasma, injection of new plasma, and boundary conditions
! on particles and fields that are specific to user's problem.
! 
! Functions by name:
! read_input_user() -- parse the user-specific section of the input file for user parameters
! get_external_fields() -- called by mover to provide externally imposed EM fields if external_fields is on in input file
! init_EMfields_user() -- initializes EM fields 
! init_particle_distribution_user() -- loads initial particle distribution
! inject_particles_user() -- injects new particles on every time step as needed
! field_bc_user() -- applies user-specific boundary conditions to fields (beyond periodic or radiation BCs)
! particle_bc_user() -- enforces user-specific BCs on particles (e.g., reflecting walls)
! shift_domain_user() -- needed for some shock problems to remove empty space; ignore. 
!
#ifdef twoD 

module m_user

	use m_globaldata
	use m_system
	use m_aux
	use m_communications
	use m_fields
	use m_particles
	use m_inputparser
	use m_fparser
	use m_domain
	
#else

module m_user_3d

	use m_globaldata_3d
	use m_system_3d
	use m_aux_3d
	use m_communications_3d
	use m_fields_3d
	use m_particles_3d
	use m_inputparser_3d 
	use m_fparser_3d
	use m_domain_3d

#endif

	implicit none
		
	private

!-------------------------------------------------------------------------------
!	PARAMETERS
!-------------------------------------------------------------------------------

	
!-------------------------------------------------------------------------------
!	TYPE DEFINITIONS
!-------------------------------------------------------------------------------

!-------------------------------------------------------------------------------
!	VARIABLES
!-------------------------------------------------------------------------------

	real(sprec) :: temperature_ratio, sigma_ext, bz_ext0

!-------------------------------------------------------------------------------
!	INTERFACE DECLARATIONS
!-------------------------------------------------------------------------------
	
!-------------------------------------------------------------------------------
!	PUBLIC MODIFIERS
!-------------------------------------------------------------------------------

	public :: init_EMfields_user, init_particle_distribution_user, &
	inject_particles_user, read_input_user, field_bc_user, get_external_fields, &
	particle_bc_user, shift_domain_user

!-------------------------------------------------------------------------------
!	MODULE PROCEDURES AND FUNCTIONS
!-------------------------------------------------------------------------------

!
! External variables read from the input file earlier; they are available for use here
!
! sigma (electron sigma),
! Binit (B field based on sigma); can be reset here
! me, mi, qe, qi (mass and charge of the species, abs(qe)/me = 1
! ppc0 (fiducial particles per cell for ions and electrons together)
! btheta, bhi -- inclination angles of initial B field, in degrees
! delgam -- k T/ m_e c^2 -- fiducial normalized temperature of electrons 
!
! mx0,my0,mz0 -- real number of cells in each direction, including the ghost zones. Active domain is from 3 to mx0-3, inclusive, and same for y and z. This is global size of the grid.
! each core knows also about mx, my, mz, which is the size of the local grid on this core. This size includes 5 ghost zones. 
!
! iglob, jglob, kglob -- functions that return global i,j,k based on local i,j,k
! xglob, yglob, zglob -- functions that return global x,y,z based on local x,y,z

	contains

!-------------------------------------------------------------------------------
! 						subroutine read_input_user		
!									
! Reads any variables related to (or needed by) this module from section "problem" in input file
! 							
!-------------------------------------------------------------------------------

subroutine read_input_user()

	implicit none
	integer :: lextflds, luserpartbcs, lwall

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!CHANGE THE NAME ON THIS LINE IF YOU ARE CREATING A NEW USER FILE
!This helps to identify which user file is being used. 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	if(rank.eq.0) print *, "Using user file user_plasmaosc.F90"

!inputpar_getd_def -- read real parameter; input_geti_def -- read integer parameter

	call inputpar_getd_def("problem", "temperature_ratio", 1._sprec, Temperature_ratio)

	call inputpar_getd_def("problem","sigma_ext",0._sprec,sigma_ext)

	call inputpar_geti_def("problem","external_fields",0,lextflds)

	if(lextflds .eq. 1) then 
	   external_fields =.true.
	else
	   external_fields =.false.
	endif

	if(external_fields) bz_ext0 = sqrt((gamma0)*.5*ppc0*c**2*me*(1+me/mi)*sigma_ext)

	call inputpar_geti_def("problem","user_part_bcs",0,luserpartbcs)

	if(luserpartbcs .eq. 1) then 
	   user_part_bcs = .true.
	else
	   user_part_bcs = .false.
	endif

	call inputpar_geti_def("problem", "wall", 0, lwall)
	
	if (lwall==1) then
		wall=.true.
	else
		wall=.false.
	endif
	
	if(wall) user_part_bcs=.true.

 if(wall .and. rank .eq. 0) print *, "resetting wall to 0" !special for Weibel case
 if(wall) wall =.false.

!read_input_user is called last, so any parameter can be overwritten here
	
end subroutine read_input_user

!-------------------------------------------------------------
!     Compute external fields to be added to the mover. 
!     These fields do not evolve via Maxwell Eqs, but can depend on time
!-------------------------------------------------------------
	subroutine get_external_fields(x,y,z,ex_ext, ey_ext, ez_ext, bx_ext,by_ext,bz_ext, qm)
	
	real,intent(inout):: bx_ext,by_ext,bz_ext, ex_ext, ey_ext, ez_ext
	real, intent(in):: x,y,z
	real, optional:: qm
	ex_ext=0.
	ey_ext=0.
	ez_ext=0.
	bx_ext=0.
	by_ext=0.
	bz_ext=bz_ext0 !defined in the input reading part above  
        !external field can be directed in any way; only invoked for external_fields = 1
        !x,y,z come in as local coordinates, convert to global by calling xglob(x), yglob(y), zglob(z)	
	end subroutine get_external_fields

!-------------------------------------------------------------------------------
! 						subroutine init_EMfields_user		 
!												
! Sets the electromagnetic fields of any specific user purpose
!							
!-------------------------------------------------------------------------------

subroutine init_EMfields_user()
	
	! local variables
	
	integer :: i, j, k

! sqrt(sigma)=omega_c/omega_p = (qe B / gamma0 me c)/sqrt( qe 0.5 ppc0 (1+me/mi)/gamma0)  !4 pi=1, qe/me = 1 
! B = sqrt(gamma0*.5*ppc0*(1+me/mi)*c**2*(me)*sigma)
! (1+me/mi) term comes from omega_p^2 = omega_pe^2+omega_pi^2

       	Binit=sqrt(gamma0*ppc0*.5*c**2*(me*(1+me/mi))*sigma)

	!initialize B field to be set by Binit and the inclination angle 	
        !angles btheta and bphi are read earlier
	btheta=btheta/180.*pi
	bphi=bphi/180.*pi

	do  k=1,mz
		do  j=1,my
			do  i=1,mx
! can have fields depend on xglob(i), yglob(j), zglob(j) or iglob(i), jglob(j), kglob(k)
				
				bx(i,j,k)=0. !Binit*cos(btheta) 
				by(i,j,k)=0. !Binit*sin(btheta)*sin(bphi)
				bz(i,j,k)=0. !Binit*sin(btheta)*cos(bphi)

				ex(i,j,k)=0.
				ey(i,j,k)=0. !(-beta)*bz(i,j,k) 
				ez(i,j,k)=0. !-(-beta)*by(i,j,k)

			enddo
		enddo
	enddo

end subroutine init_EMfields_user


!-------------------------------------------------------------------------------
! 						subroutine init_particle_distribution_user()	
!											
! Sets the particle distrubtion for a user defined case
!
!-------------------------------------------------------------------------------

subroutine init_particle_distribution_user()

	implicit none

	! local variables
	
	integer :: i, n, direction
	real    gamma_drift, delgam_i, delgam_e, ppc, weight

	real(sprec) :: x1,x2,y1,y2,z1,z2

	pcosthmult=0 !if 0 the Maxwellian distribution corresponding to 
	             ! temperature is initialized in 2D (z temperature = 0), 1 for 3D distribution (z temp nonzero). 
                     ! 1 is default

	!set initial injection points 
	
	xinject=3. 
	xinject2=(mx0-2.) 

	! ------------------------------------------

!initialize plasma

	gamma_drift=0. !-gamma0 ! negative gamma_drift will send the plasma in the negative direction
	delgam_i=delgam  !delgam is read from input file in particles
	delgam_e=delgam*mi/me*Temperature_ratio

	x1=xinject  
	x2=xinject2 

	y1=3. !in global coordinates
	y2=my0-2.  
	z1=3.
	z2=mz0-2. !if 2D, it will reset automatically
	ppc=ppc0
	weight=1.

	direction=1  !drift along x, + or - determined by the sign of gamma_drift

!create a region of plasma with maxwellian distribution and drift  
	call inject_plasma_region(x1,x2,y1,y2,z1,z2,ppc,&
             gamma_drift,delgam_i,delgam_e,weight,direction)

!now go through all electrons and give them a nudge with a sinusoid.
!beta is velocity derived from gamma0 -- set in particles.F90. 
!gamma0 comes from input file. Use beta as amplitude of the kick
!
!number of ions and electrons on CPU: ions and lecs
!ions are located in p array from 1 to "ions". 
!electrons are in p array from maxhlf+1 to maxhlf+lecs
!maxhlf is the index of the midpoint of the array

 do n=maxhlf+1, maxhlf+lecs

!    p(n)%x = 3.+((mx-5.)/lecs)*(n-maxhlf)  !for quiet start to kill noise
!    p(n)%u = 0. !if you want to set temperature to 0. 
!    p(n)%v = 0. !
!    p(n)%w = 0. !

    p(n)%u = p(n)%u + beta*sin(2*pi/(mx0-5.) * xglob(p(n)%x) )

 enddo


 !now remove all ions; this makes ions immobile
       ions=0 

	call check_overflow()
	call reorder_particles()
	
end subroutine init_particle_distribution_user


!-------------------------------------------------------------------------------
! 				subroutine inject_particles_user()					 
!										
! Injects new particles in the simulation on every step. To skip injection, set ppc=0 below
!
!-------------------------------------------------------------------------------

subroutine inject_particles_user()

	implicit none
	real(sprec) :: x1,x2,y1,y2,z1,z2
	real delgam_i, delgam_e, injector_speed, ppc, betainj, gamma_drift, weight

	injectedions=0 !set counter to 0 here, so that all possible other injections on this step add up
	injectedlecs=0

!no injection is needed for weibel case
	
end subroutine inject_particles_user



!-------------------------------------------------------------------------------
! 				subroutine field_bc_user()	
!										
! Applies boundary conditions specific to user problem. 
! 
!-------------------------------------------------------------------------------

	subroutine field_bc_user()
	implicit none

!no special field BCs needed. Periodic conditions are taken care of in fieldboundaries

	end subroutine field_bc_user

!------------------------------------------------------------------------------

	subroutine particle_bc_user()
	implicit none
	real invgam, walloc, gammawall, betawall, gamma, x0, walloc0,xcolis,t0,t1
	integer n1,i0,i1, iter

	!loop over particles to check if they crossed special boundaries, like reflecting walls
	!outflow and periodic conditions are handled automatically in deposit_particles
	!
	!This routine is called after the mover and before deposit, thus allowing to avoid
        ! charge deposition behind walls. 
	
        if(1<0) then 	
	   do iter=1,2
	      if(iter.eq.1) then 
		 i0=1
		 i1=ions
	      else
		 i0=maxhlf+1
		 i1=maxhlf+lecs
	      endif
	      
	      do n1=i0,i1
		 
! no special conditions needed, this routine is not called if user_part_bcs=0 in input file
		 
	      enddo
	   enddo
 endif !skip loop
	end subroutine particle_bc_user

!-------------------------------------------------------------------------------
! 				subroutine shift_domain_user
!										
! shift fields and particles backward
! only used for some shock problems. Ignore but don't delete. 
!-------------------------------------------------------------------------------
 subroutine shift_domain_user
   implicit none

 end subroutine shift_domain_user
 
#ifdef twoD
end module m_user
#else
end module m_user_3d
#endif

