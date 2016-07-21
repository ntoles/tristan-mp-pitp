!
! Aux module
!
! This module defines auxiliary functions / procedures usefull for at least two 
! different modules of the code
!
!

#ifdef twoD 

module m_aux

	use m_system
	use m_globaldata

#else

module m_aux_3d

	use m_system_3d
	use m_globaldata_3d

#endif


	implicit none
	
	private

!-------------------------------------------------------------------------------
!	PARAMETERS
!-------------------------------------------------------------------------------

#ifdef DEBUG
	logical, parameter :: debug=.true.
#else
	logical, parameter :: debug=.false.
#endif
	
!-------------------------------------------------------------------------------
!	TYPE DEFINITIONS
!-------------------------------------------------------------------------------

!-------------------------------------------------------------------------------
!	VARIABLES
!-------------------------------------------------------------------------------

	real(dprec) :: dseed
	integer :: lap, irestart, lapst !caseinit
	
!-------------------------------------------------------------------------------
!	INTERFACE DECLARATIONS
!-------------------------------------------------------------------------------
	
!-------------------------------------------------------------------------------
!	PUBLIC MODIFIERS
!-------------------------------------------------------------------------------

	! public functions

	public :: random, poisson

	! public variables

	public :: debug, lap, irestart, lapst, dseed ! caseinit,

!-------------------------------------------------------------------------------
!	MODULE PROCEDURES AND FUNCTIONS
!-------------------------------------------------------------------------------

	contains



!-------------------------------------------------------------------------------
! 						function random					 
!																		
! Random number generator between [0,1]
!							
!-------------------------------------------------------------------------------

real Function random(dseed)       
	
	implicit none
	
	real(dprec) :: DSEED
	integer ::I
	real(dprec) :: S2P31,S2P31M,SEED
	DATA              S2P31M/2147483647.D0/,S2P31/2147483648.D0/
	
	SEED = DSEED
	
	SEED = DMOD(16807.D0*SEED,S2P31M)
	random= SEED /S2P31
	DSEED = SEED
	
	return
	
end function



!-------------------------------------------------------------------------------
! 						subroutine poisson()					 
!																		
! Returns random numbers with a poisson distribution
!
!-------------------------------------------------------------------------------

real function poisson(numps)

	implicit none
	
	! dummy variables
	
	real(sprec) :: numps
	
	! local variables
	
	real(dprec) :: Lps, pps
	real kps, ups
	
	
	Lps=exp(-real(numps,8))
	kps=0
	pps=1
	do while ( pps .ge. Lps)
		kps=kps+1
		ups=random(dseed)
		pps=pps*ups
	enddo
	poisson=kps-1

end function poisson

#ifdef twoD
end module m_aux
#else
end module m_aux_3d
#endif
