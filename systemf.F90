!
! System module
!
! This module contains all system dependent functions, procedures and data  
! structures
!
!

#ifdef twoD 
module m_system
#else
module m_system_3d
#endif

	implicit none
	
!-------------------------------------------------------------------------------
!	PARAMETERS
!-------------------------------------------------------------------------------

    integer, parameter :: dprec = kind(1.0d0)
    integer, parameter :: sprec = kind(1.0e0)
	
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

	public :: dprec, sprec

!-------------------------------------------------------------------------------
!	MODULE PROCEDURES AND FUNCTIONS
!-------------------------------------------------------------------------------

	contains


!------------------------------------------------------------------------
! 						  subroutine mkdir							 
!																		
! Creates a directory tree as specified, making use of the c module system.c 
! 
!------------------------------------------------------------------------

subroutine mkdir(path, ierr)

	implicit none
   
   	! dummy variables
   
    character(len=*), intent(in) :: path
    integer, intent(out) :: ierr
    
    ! local variables
           
    character(len=len_trim(path)+1) :: lpath


    ! Make shure this is a proper C string
    lpath = trim(path)//char(0)
           
    call mkdir_f(lpath, ierr)
    
end subroutine mkdir



!------------------------------------------------------------------------
! 						  subroutine Delete							 
!																		
! Deletes a specified file making use of the c module system.c 
! 
!------------------------------------------------------------------------

subroutine Delete(fname, ierr)

	implicit none

   	! dummy variables
   
    character(len=*), intent(in) :: fname
    integer, intent(out) :: ierr
    
    ! local variables
           
    character(len=len_trim(fname)+1) :: lfname


    ! Make shure this is a proper C string
    lfname = trim(fname)//char(0)
           
    call remove_f(lfname, ierr)
	

end subroutine Delete



!------------------------------------------------------------------------
! 						  subroutine Move							 
!																		
! Moves a file making use of the c module system.c 
! 
!------------------------------------------------------------------------

subroutine Move(fnameprev, fnamefinal, ierr)

	implicit none

   	! dummy variables
   
    character(len=*), intent(in) :: fnameprev, fnamefinal
    integer, intent(out) :: ierr
    
    ! local variables
           
    character(len=len_trim(fnameprev)+1) :: lfnameprev
    character(len=len_trim(fnamefinal)+1) :: lfnamefinal


    ! Make sure this is a proper C string
    lfnameprev = trim(fnameprev)//char(0)
    lfnamefinal = trim(fnamefinal)//char(0)
           
    call move_f(lfnameprev, lfnamefinal, ierr)
	

end subroutine Move




#ifdef twoD 
end module m_system
#else
end module m_system_3d
#endif
