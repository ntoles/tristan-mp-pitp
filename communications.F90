!
! Communications module
!
! Holds the communications (MPI) related data structures and procedures 
! 
!
!

#ifdef twoD 

module m_communications

	use m_system
	use m_aux
	use m_globaldata
	use m_inputparser
	
#else

module m_communications_3d

	use m_system_3d
	use m_aux_3d
	use m_globaldata_3d
	use m_inputparser_3d

#endif


	implicit none
	
	include "mpif.h"
	
	private

!-------------------------------------------------------------------------------
!	PARAMETERS
!-------------------------------------------------------------------------------
	
	integer, parameter :: tm_size=50

!-------------------------------------------------------------------------------
!	TYPE DEFINITIONS
!-------------------------------------------------------------------------------

!-------------------------------------------------------------------------------
!	VARIABLES
!-------------------------------------------------------------------------------

	integer :: rank, size0, sizex,sizey,sizez,ierr,statsize
	integer :: id1, id2, id3, id4, mpi_read, dst, mpi_dbl
	integer :: particletype
	real(dprec), dimension(tm_size) :: tm

!-------------------------------------------------------------------------------
!	INTERFACE DECLARATIONS
!-------------------------------------------------------------------------------

	interface init_communications
		module procedure init_communications
	end interface init_communications
	
	interface timer
		module procedure timer
	end interface timer


!-------------------------------------------------------------------------------
!	PUBLIC MODIFIERS
!-------------------------------------------------------------------------------

	! public functions (used in m_tris)
	
	public :: print_timers, init_communications, initialize_random_seed, timer, read_input_communications

	! public variables

	public :: rank, size0, sizex, sizey, sizez, statsize, mpi_read, mpi_dbl,MPI_COMM_WORLD, mpi_wtime, mpi_integer,&
			  mpi_sum, MPI_STATUS_SIZE, mpi_integer8, MPI_MAX, MPI_MIN, MPI_INFO_NULL, mpi_real,  &
			  particletype

	! Used by create_MPI_datatypes() in fieldboundaries.F90
	! Apparently this is needed to pass through variables from mpif.h to the rest of the code
	public :: MPI_ORDER_FORTRAN, tm

!-------------------------------------------------------------------------------
!	MODULE PROCEDURES AND FUNCTIONS
!-------------------------------------------------------------------------------

	contains



!-------------------------------------------------------------------------------
! 						subroutine read_input_communications					 
!																		
! Reads any variables related to (or needed by) this module
! 							
!-------------------------------------------------------------------------------

subroutine read_input_communications()

	implicit none
	
	call inputpar_geti_def("node_configuration", "sizex", 1, sizex)
	call inputpar_geti_def("node_configuration", "sizey", 1, sizey)

end subroutine read_input_communications



!-------------------------------------------------------------------------------
! 						subroutine init_communications					 
!																		
! Calls MPI init, and creates the MPI data structures for communication
! 							
!-------------------------------------------------------------------------------

subroutine init_communications()

	implicit none
	
	integer, dimension(0:1) :: oldtypes, blockcounts, offsets
	integer :: extent
	
	
	mpi_read=MPI_REAL 
	mpi_dbl=MPI_DOUBLE_PRECISION
	size0=1

	call MPI_Init(ierr)
	call MPI_Comm_rank(MPI_Comm_world, rank, ierr)
	call MPI_Comm_size(MPI_Comm_world, size0, ierr)

	statsize=MPI_STATUS_SIZE

	!set up description of the particle type fields
	!for double field

!	offsets(0)=0
!	oldtypes(0)=MPI_REAL8
!	blockcounts(0)=2

	!for real fields

!	call MPI_TYPE_EXTENT(MPI_REAL8,extent,ierr)

	offsets(0)=0
	oldtypes(0)= MPI_REAL
	blockcounts(0)=7

	!for integer ::fields

	call MPI_TYPE_EXTENT(MPI_REAL,extent,ierr)
	offsets(1)=7*extent+offsets(0)
	oldtypes(1)=MPI_INTEGER
	blockcounts(1)=2+1

	call MPI_TYPE_STRUCT(2,blockcounts,offsets,oldtypes,particletype,ierr)

	call MPI_TYPE_COMMIT(particletype, ierr)

	sizez=size0/sizey/sizex
	
#ifdef twoD
	sizez=1
#endif
	
	if(sizez*sizey*sizex .ne. size0) then
	   if(rank .eq. 0) print *, rank,":", "Error: sizex*sizey ne size0","size0=",size0,"sizez=",sizez,"sizey=",sizey
	   stop
	endif

	tm(:)=0.

#ifndef twoD
	if(sizex .ne. 1) then 
	   if(rank .eq. 0) print *, "In 3D, domain decomposition is not supported in x direction. Use sizex=1 in input file."
	   stop
	endif
#endif


end subroutine init_communications



!-------------------------------------------------------------------------------
! 						subroutine timer					 
!																		
! Calls the mpi timer to time the execution of the code
! 							
!-------------------------------------------------------------------------------

subroutine timer(i, tmstop)

	implicit none
	
	! dummy variables
	
	integer, intent(in) :: i
	logical, intent(in), optional :: tmstop
	

!	if (rank == 0) then
		if (present(tmstop)) then	! stop the timer
			tm(i)=mpi_wtime()-tm(i)
		else						! start the timer
			tm(i)=mpi_wtime()
		endif
!	endif
	
end subroutine timer



!-------------------------------------------------------------------------------
! 						subroutine initialize_random_seed					 
!																		
! Initializes the random seed for the random number generator 
!							
!-------------------------------------------------------------------------------

subroutine initialize_random_seed()

	implicit none

	dseed=123457.D0        
	dseed=dseed+rank

end subroutine initialize_random_seed



!-------------------------------------------------------------------------------
! 						subroutine print_timers					 
!																		
! Prints the timers
! 							
!-------------------------------------------------------------------------------

subroutine print_timers()
	implicit none
	real(dprec) tmax
	integer i,n
	real mean, minx, maxx, stddev
	character(20) timername
!	real(dprec) timtemp(size0*tm_size), arr(size0), time_arr(tm_size,size0) !, mean, stddev , maxval, minval
	real(sprec) tm4(tm_size),arr(size0), time_arr(tm_size,size0), summean !, mean, stddev , maxval, minval
	
	time_arr(:,:)=0.
	arr(:)=0.
	tm4(:)=0.

	tm4=real(tm,4) !send around single precision numbers

	call mpi_gather(tm4,tm_size,mpi_real,time_arr,tm_size,mpi_real,0,mpi_comm_world,ierr)

	if(rank == 0) then 
	   summean=0.
		write(*,*) "Timers for this step ----------------------------------"
		
		do n=1,8
		select case (n)
		   case (1)
		      arr=time_arr(2,:)+time_arr(4,:)+time_arr(9,:)	
		      timername="Field slv -->"
		   case (2)
		      arr=time_arr(3,:)
		      timername="Mover     -->" 
		   case (3)
		      arr=time_arr(5,:)
		      timername="Deposit   -->" 
		   case (4)
		      arr=time_arr(6,:)
		      timername="Part exch -->" 
		   case (5)
		      arr=time_arr(7,:)
		      timername="Currnt exch->" 
		   case (6)
		      arr=time_arr(8,:)
		      timername="Filter    -->" 
		   case (7)
		      arr=time_arr(10,:)
		      timername="Inject    -->" 
		   case (8)
		      arr=time_arr(1,:)
		      timername="=Total, sec=>"
		   case default
		   end select

		   mean=real(sum(arr)/(size0),4)
		   maxx=real(maxval(arr),4)
		   minx=real(minval(arr),4)
		   stddev=sqrt(sum((arr-mean)**2)/size0)
!		   if(n.ne.8)summean=summean+mean

		write(*,'(A13,ES9.2,1x,SP,A3,F7.1,S,A6,F5.1,A8,F5.1,A7,I4,I4, &
	A7,ES9.2)') timername,mean, "max",(maxx-mean)/max(1e-7,mean)*100, &
	"% min ",-1*(mean-minx)/max(1e-7,mean)*99.9, "% stddev", &
	min(stddev/mean*100,999.0), "% nodes",int(maxloc(arr)),int(minloc(arr)) !, " node0 ",arr(1)

	     enddo !n
	     
!	     print *, "sum of means", summean, "real total", time_arr(1,n)
!	     do n=1,size0,5
!		print *, "on node ", n, "timer sum=", sum(time_arr(2:10,n)), "real total", time_arr(1,n)
!	     enddo
	     
		write(*,*) "------------------------------------------------------"

	     endif !rank ? 0
 

	     if(rank.eq.0 .and. 1<0) then
		write(*,*) "Timers for this step ----------------------------------"
		write(*,*) "Field solve -->",real(tm(2)+tm(4)+tm(9),4),"   Mover -->",real(tm(3),4)
		write(*,*) "Deposit -->",real(tm(5),4),"   Particle exchange -->",real(tm(6),4)
		write(*,*) "Current Exchange -->",real(tm(7),4),"   Filter -->",real(tm(8),4)
		write(*,*) "Inject -->",real(tm(10),4),"   Mov window & Enl Domain -->",real(tm(11),4)
		write(*,*) "Reorder -->",real(tm(12),4),"   Redistribute x Domain -->",real(tm(31),4)
		write(*,*) "Redistribute y Domain -->",real(tm(32),4),"    Total time -->",real(tm(1),4)
		write(*,*) "------------------------------------------------------"
	     endif

end subroutine print_timers

#ifdef twoD
end module m_communications
#else
end module m_communications_3d
#endif
