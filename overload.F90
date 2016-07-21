!
! Overload module
!
! This module contains all the branching that is done in tristan to either a default tristan function
! or a user defined function in the module m_user.
!
! A user may alter this module in the following way:
! 	- For every select case in this file add a line with a different number of your choice > 0
!	- Call your user supplied function whatever you like
!	- Supply this function either in m_user or in a different file (be sure to make tristan aware 
!	of this new file by including here the line "use m_..." and by including the file name in 
!	the Makefile in the correct place


#ifdef twoD 

module m_overload

	use m_globaldata
	use m_aux
	use m_particles
	use m_fields
	use m_user
	use m_domain
	use m_inputparser
	use m_communications
	
#else

module m_overload_3d

	use m_globaldata_3d
	use m_aux_3d
	use m_particles_3d
	use m_fields_3d
	use m_user_3d
	use m_domain_3d
	use m_inputparser_3d
	use m_communications_3d

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


!-------------------------------------------------------------------------------
!	INTERFACE DECLARATIONS
!-------------------------------------------------------------------------------
	
!-------------------------------------------------------------------------------
!	PUBLIC MODIFIERS
!-------------------------------------------------------------------------------

	! public functions

	public :: init_EMfields, init_particle_distribution, inject_particles, initialize_domain, &
	read_input_user, field_bc_user, shift_domain !, shift_ydomain 

!-------------------------------------------------------------------------------
!	MODULE PROCEDURES AND FUNCTIONS
!-------------------------------------------------------------------------------

	contains



!-------------------------------------------------------------------------------
! 						subroutine read_input_user					 
!																		
! Reads any variables related to (or needed by) this module
! 							
!-------------------------------------------------------------------------------

subroutine read_input()

	implicit none

	! local variables
	call read_input_user()

end subroutine read_input



!-------------------------------------------------------------------------------
! 						subroutine init_EMfields					 
!																		
! Sets the electromagnetic fields of any specific user purpose
!							
!-------------------------------------------------------------------------------

subroutine init_EMfields()

	implicit none

!	select case(caseinit)
!		case(0)
!			call init_EMfields_weibel()
!		case(1)
!			call init_EMfields_shock()
!		case default
!		   if(rank .eq. 0) print *, "No such caseinit defined in overload:",caseinit
!	end select

	call init_EMfields_user()

end subroutine init_EMfields



!-------------------------------------------------------------------------------
! 						subroutine init_particle_distribution					 
!																		
! Calls the correct function to set the initial particle distribution
! 							
!-------------------------------------------------------------------------------

subroutine init_particle_distribution()

	implicit none

!	select case(caseinit)
!		case(0)
!			call init_particle_distribution_weibel()
!		case(1)
!			call init_particle_distribution_shock()
!		case default
!		   if(rank .eq. 0) print *, "No such caseinit defined in overload:",caseinit
!	end select
	
	call init_particle_distribution_user()

end subroutine init_particle_distribution



!-------------------------------------------------------------------------------
! 						subroutine inject_particles					 
!																		
! Calls the correct function to inject the particles into the simulation box
! 							
!-------------------------------------------------------------------------------

subroutine inject_particles()
	
	implicit none


!	select case(caseinit)
!		case(0)
!			call inject_particles_weibel()
!		case(1)
!			call inject_particles_shock()
!		case default
!		   if(rank .eq. 0) print *, "No such caseinit defined in overload:",caseinit
!	end select
	
	call inject_particles_user()
	
end subroutine inject_particles


!-------------------------------------------------------------------------------
! 						subroutine shift_domain				 							
!-------------------------------------------------------------------------------

subroutine shift_domain()
	
	implicit none
	
	call shift_domain_user()
	
end subroutine shift_domain

!-------------------------------------------------------------------------------
! 						subroutine shift_domain				 							
!-------------------------------------------------------------------------------

!subroutine shift_ydomain()
	
!	implicit none

	
!	call shift_ydomain_user()
	
!end subroutine shift_ydomain



!-------------------------------------------------------------------------------
! 						subroutine field_bc_user 
!									
! Calls the correct function to apply problem-specific field BC
! Periodic BCs are handled in bc_e1 and bc_b1 in fieldboundaries.F90
!
!-------------------------------------------------------------------------------

subroutine field_bc()
	
	implicit none

!	select case(caseinit)
!		case(0)
!			call field_bc_weibel()
!		case(1)
!			call field_bc_shock()
!		case default
!		   if(rank .eq. 0) print *, "No such caseinit defined in overload:",caseinit
!	end select
	
	call field_bc_user()

end subroutine field_bc



!-------------------------------------------------------------------------------
! 						subroutine initialize_domain					 
!										
! Calls the correct function to set the m_domain variables correctly
! 							
!-------------------------------------------------------------------------------

subroutine initialize_domain()

	implicit none


!	select case(caseinit)
!		case(0)
!			call initialize_domain_default() !sets enlarge to 0 
!	end select
	
	
end subroutine initialize_domain

#ifdef twoD
end module m_overload
#else
end module m_overload_3d
#endif
