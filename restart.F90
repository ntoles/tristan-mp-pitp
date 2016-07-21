!
! Restart module
!
! This module is responsible for handling the restarts of the code 
! (starting the simulation from restart files from previous simulations)
!
!


#ifdef twoD 

module m_restart

	use m_globaldata
	use m_system
	use m_aux
	use m_communications
	use m_fields
	use m_particles
	use m_output
	use m_domain
	use m_inputparser
	
#else

module m_restart_3d

	use m_globaldata_3d
	use m_system_3d
	use m_aux_3d
	use m_communications_3d
	use m_fields_3d
	use m_particles_3d
	use m_output_3d
	use m_domain_3d
	use m_inputparser_3d

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

	public :: restart, initialize_restart, read_input_restart

!-------------------------------------------------------------------------------
!	MODULE PROCEDURES AND FUNCTIONS
!-------------------------------------------------------------------------------

	contains



!-------------------------------------------------------------------------------
! 						subroutine read_input_restart					 
!																		
! Reads any variables related to (or needed by) this module
! 							
!-------------------------------------------------------------------------------

subroutine read_input_restart()

	implicit none
	
	! local variables
	
	integer :: llaprestart, lnamedrestartint
		
	call inputpar_geti_def("restart", "irestart", 0, irestart)
	call inputpar_geti_def("restart", "intrestart", 2000, intrestart)
	call inputpar_geti_def("restart", "laprestart", 0, llaprestart)
	call inputpar_geti_def("restart", "lnamedrestartint", 1000000, lnamedrestartint)
	
	laprestart=llaprestart
	namedrestartint=lnamedrestartint
	
end subroutine read_input_restart



!-------------------------------------------------------------------------------
! 						subroutine initialize_restart					 
!
!	Initializes the module and gets restart parameters if restarting
!
!-------------------------------------------------------------------------------

subroutine initialize_restart()

	implicit none


	call SetRestartStrings()
	

end subroutine initialize_restart



!-------------------------------------------------------------------------------
! 						subroutine SetRestartStrings					 
!
!	Sets the strings needed in the module
!
!-------------------------------------------------------------------------------

subroutine SetRestartStrings()

	implicit none

	logical :: exst

	write(frootprt,"(a5,i3.3,a1)") "prtl.",rank,"."
	write(frootfld,"(a5,i3.3,a1)") "flds.",rank,"."
	
	!special numbering if we are on more than 10^3 processors
	if(size0 .le. 1000) then 
		write(rankchar,"(i3.3)")rank
	else
		write(rankchar,"(i4.4)")rank
	endif
	
	!use jobname from the file created by submit script
	inquire(file="jobname",exist=exst)
	
	if(exst) then 
		open(unit=10,file="jobname",form='formatted')
		read(10,*)jobname
		if(rank.eq.0)print *,rank, ": jobname=", jobname
	else
		jobname="default"
	endif
	
	
	frestartfld="restart/restflds."//trim(jobname)//"."// & 
	trim(rankchar)//".d"
	frestartprt="restart/restprtl."//trim(jobname)//"."// &
	trim(rankchar)//".d"
	if(laprestart .ne. 0) then
		write(lapchar,"(i6.6)")laprestart	
		frestartfldlap="restart/restflds."//trim(jobname) &
		//"."//"lap"//trim(lapchar)//"."//trim(rankchar)//".d"
		frestartprtlap="restart/restprtl."//trim(jobname) &
		//"."//"lap"//trim(lapchar)//"."//trim(rankchar)//".d"
	else
		frestartfldlap=frestartfld
		frestartprtlap=frestartprt
	endif	
	
	
	!special architectures that can use local disks

#ifdef LOCALSCRATCH
		frestartfldloc="/scratch/restflds."//trim(jobname)//"."// &
		trim(rankchar)//".d"
		frestartprtloc="/scratch/restprtl."//trim(jobname)//"."// &
		trim(rankchar)//".d"
		fenlargeloc="/scratch/fenlarge."//trim(jobname)//"."// &
		trim(rankchar)//".d"
#else
		frestartfldloc="restart/restflds."//trim(jobname)//"."// &
		trim(rankchar)//".d"
		frestartprtloc="restart/restprtl."//trim(jobname)//"."// &
		trim(rankchar)//".d"
		fenlargeloc="restart/fenlarge."//trim(jobname)//"."// &
		trim(rankchar)//".d"
#endif
	
	write(froottrq,"(a5,i3.3,a1)") "torq.",rank,"."
	

end subroutine SetRestartStrings



!-------------------------------------------------------------------------------
! 						subroutine restart					 
!																		
!
! 							
!-------------------------------------------------------------------------------

subroutine restart()

	implicit none
	
	! local variables

	integer :: dummy, ierr
	integer :: n, i,j,k, token, request,status(MPI_STATUS_SIZE),request1 &
	,status1(MPI_STATUS_SIZE)
	
		
	open(unit=7,file=frestartfldlap,form='unformatted')
	if(debug) print *, rank,": in restart, opening",frestartfldlap 
	rewind(7)
	read(7)mxrest,myrest,mzrest,(((ex(i,j,k),i=1,mxrest),j=1,myrest),k=1,mzrest) &
	,(((ey(i,j,k),i=1,mxrest),j=1,myrest),k=1,mzrest),(((ez(i,j,k),i=1 &
	,mxrest),j=1,myrest),k=1,mzrest),(((bx(i,j,k),i=1,mxrest),j=1,myrest),k=1 &
	,mzrest),(((by(i,j,k),i=1,mxrest),j=1,myrest),k=1,mzrest),(((bz(i,j,k),i &
	=1,mxrest),j=1,myrest),k=1,mzrest),dseed,lapst,xinject,xinject2,xinject3 &
	, leftwall, walloc !,split_E_ions,split_E_lecs
	close(7)
	print *,rank,": Read fields: lapst=",lapst,"xinj=",xinject &
	,"xinj2 =",xinject2 ,"walloc=",walloc,"dseed=", dseed
	
	lapst=lapst+1
	open(unit=8,file=frestartprtlap,form='unformatted')
	rewind(8)
	maxptl=maxptl0/size0
	maxhlf=maxptl/2
	
	print *, rank, ": in restart. Reading particles "
	! ions,lecs,maxptl,maxhlf,totalpartnum, ...
	read(8) ions,lecs,dummy,dummy, totalpartnum, &
	(p(n)%x,n=1,ions),(p(n)%x,n=maxhlf+1,maxhlf+lecs), &
	(p(n)%y,n=1,ions),(p(n)%y,n=maxhlf+1,maxhlf+lecs), &
	(p(n)%z,n=1,ions),(p(n)%z,n=maxhlf+1,maxhlf+lecs), &
	(p(n)%u,n=1,ions),(p(n)%u,n=maxhlf+1,maxhlf+lecs), &
	(p(n)%v,n=1,ions),(p(n)%v,n=maxhlf+1,maxhlf+lecs), &
	(p(n)%w,n=1,ions),(p(n)%w,n=maxhlf+1,maxhlf+lecs), &
	(p(n)%ch,n=1,ions),(p(n)%ch,n=maxhlf+1,maxhlf+lecs), &
	(p(n)%ind,n=1,ions),(p(n)%ind,n=maxhlf+1,maxhlf+lecs),   &
	(p(n)%proc,n=1,ions),(p(n)%proc,n=maxhlf+1,maxhlf+lecs), &
	(p(n)%splitlev,n=1,ions),(p(n)%splitlev,n=maxhlf+1,maxhlf+lecs)
	close(8)
	print *, rank,": read ions=",ions,"lecs=",lecs,"maxhlf" &
	,maxhlf
	print *,rank,": read particles, indexing"
	
	!            call index_particles()
	call reorder_particles()
	print *,rank,": done restart"
	
#ifdef MPI
!		goto 7321
!		if(modulo(rank,10).eq.0 .and. rank .ne. size0-1) then
!		token=rank
!		!              call MPI_Send(token,1,mpi_integer,rank+1,100
!		!     &                  ,MPI_COMM_WORLD,ierr)
!		
!		call mpi_isend(token,1, &
!		mpi_integer,rank+1,100,MPI_Comm_WORLD,request,ierr)
!		endif
!		7321      continue
		call mpi_barrier(MPI_COMM_WORLD,ierr)
	
#endif

end subroutine restart


#ifdef twoD
end module m_restart
#else
end module m_restart_3d
#endif
