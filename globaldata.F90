!
! Global data module
!
! This module contains global variables, generaly needed throughout the code
! and that do not logically belong anywhere else. Please add variables (and functions
! or procedures) to this module only as needed.
!

#ifdef twoD 
module m_globaldata

	use m_system

#else
module m_globaldata_3d

	use m_system_3d

#endif

	implicit none
	
	public
	
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

!-------------------------------------------------------------------------------
!	MODULE PROCEDURES AND FUNCTIONS
!-------------------------------------------------------------------------------

!	contains


#ifdef twoD 
end module m_globaldata
#else
end module m_globaldata_3d
#endif
