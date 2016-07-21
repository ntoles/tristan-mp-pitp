!
! Fortran input parser interface module. This module provides an easy way to call
! the par.c functions for parsing the input file.
!

#ifdef twoD 
module m_inputparser
#else
module m_inputparser_3d
#endif

	implicit none
	
	private
	
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
	
	interface inputpar_getd_def
		module procedure inputpar_getd_def_sprec
		module procedure inputpar_getd_def_dprec
	end interface
	
!-------------------------------------------------------------------------------
!	PUBLIC MODIFIERS
!-------------------------------------------------------------------------------

	public :: inputpar_open, inputpar_gets_def, inputpar_getd_def, inputpar_geti_def, &
			  inputpar_close

!-------------------------------------------------------------------------------
!	MODULE PROCEDURES AND FUNCTIONS
!-------------------------------------------------------------------------------

	contains


!------------------------------------------------------------------------
! 						  subroutine inputpar_open							 
!																		
! Opens the input file and reads all parameters into memory
! 
!------------------------------------------------------------------------

subroutine inputpar_open(filename)

	implicit none
   
   	! dummy variables
   
    character(len=*), intent(in) :: filename
    
    ! local variables
           
    character(len=len_trim(filename)+1) :: filename_c


    ! Make shure this is a proper C string
    filename_c = trim(filename)//char(0)
           
    call par_open_f(filename_c)
    
end subroutine inputpar_open



!------------------------------------------------------------------------
! 						  subroutine inputpar_gets_def							 
!																		
! Reads a string from the input file, or returns a default value 
! 
!------------------------------------------------------------------------

subroutine inputpar_gets_def(blockname, name, def, value)

	implicit none

   	! dummy variables
   
    character(len=*), intent(in) :: blockname, name, def
    character(len=*), intent(inout) :: value
    
    ! local variables
           
    character(len=len_trim(blockname)+1) :: blockname_c
    character(len=len_trim(name)+1) :: name_c
    character(len=len_trim(def)+1) :: def_c


    ! Make shure this is a proper C string
    
    blockname_c = trim(blockname)//char(0)
    name_c = trim(name)//char(0)
    def_c = trim(def)//char(0)
    
    call par_gets_def_f(blockname_c,name_c,def_c,value)

end subroutine inputpar_gets_def



!------------------------------------------------------------------------
! 						  subroutine inputpar_geti_def							 
!																		
! Reads an integer from the input file, or returns a default value
! 
!------------------------------------------------------------------------

subroutine inputpar_geti_def(blockname, name, def, value)

	implicit none

   	! dummy variables
   
    character(len=*), intent(in) :: blockname, name
    integer, intent(in) :: def
    integer, intent(out) :: value
    
    ! local variables
           
    character(len=len_trim(blockname)+1) :: blockname_c
    character(len=len_trim(name)+1) :: name_c

    
    ! Make shure this is a proper C string
    
    blockname_c = trim(blockname)//char(0)
    name_c = trim(name)//char(0)
           
    call par_geti_def_f(blockname_c, name_c, def, value)

end subroutine inputpar_geti_def



!------------------------------------------------------------------------
! 						  subroutine inputpar_getd_def							 
!																		
! Reads a double from the input file, or returns a default value
! 
!------------------------------------------------------------------------

subroutine inputpar_getd_def_sprec(blockname, name, def, value)

	implicit none

   	! dummy variables
   
    character(len=*), intent(in) :: blockname, name
    real(sprec), intent(in) :: def
    real(sprec), intent(out) :: value
    
    ! local variables
           
    character(len=len_trim(blockname)+1) :: blockname_c
    character(len=len_trim(name)+1) :: name_c
	real(dprec) :: readvalue
    
    ! Make shure this is a proper C string
    
    blockname_c = trim(blockname)//char(0)
    name_c = trim(name)//char(0)
           
    call par_getd_def_f(blockname_c, name_c, def, readvalue)
    
    value=readvalue

end subroutine inputpar_getd_def_sprec

!------------------------------------------------------------------------

subroutine inputpar_getd_def_dprec(blockname, name, def, value)

	implicit none

   	! dummy variables
   
    character(len=*), intent(in) :: blockname, name
    real(dprec), intent(in) :: def
    real(dprec), intent(out) :: value
    
    ! local variables
           
    character(len=len_trim(blockname)+1) :: blockname_c
    character(len=len_trim(name)+1) :: name_c

    
    ! Make shure this is a proper C string
    
    blockname_c = trim(blockname)//char(0)
    name_c = trim(name)//char(0)
           
    call par_getd_def_f(blockname_c, name_c, def, value)

end subroutine inputpar_getd_def_dprec



!------------------------------------------------------------------------
! 						  subroutine inputpar_close							 
!																		
! Closes the par module
! 
!------------------------------------------------------------------------

subroutine inputpar_close()

	implicit none
	
	call par_close_f()

end subroutine inputpar_close


#ifdef twoD 
end module m_inputparser
#else
end module m_inputparser_3d
#endif
