!
! Initialize module
!
! This module reads the input file, and initializes all data structures in the 
! code
!
!


#ifdef twoD 

module m_initialize

	use m_system
	use m_aux
	use m_communications
	use m_particles
	use m_fields
	use m_fieldboundaries
	use m_domain
	use m_dynamic_domain
	use m_output
	use m_restart
	use m_overload
	use m_globaldata
	use m_inputparser
	
#else

module m_initialize_3d

	use m_system_3d
	use m_aux_3d
	use m_communications_3d
	use m_particles_3d
	use m_fields_3d
	use m_fieldboundaries_3d
	use m_domain_3d
	use m_dynamic_domain_3d
	use m_output_3d
	use m_restart_3d
	use m_overload_3d
	use m_globaldata_3d
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

	real(sprec) :: tcutoff
	integer :: nhist, nbins, nc
	character(len=32) :: input_file_name="input"
!-------------------------------------------------------------------------------
!	INTERFACE DECLARATIONS
!-------------------------------------------------------------------------------
	
!-------------------------------------------------------------------------------
!	PUBLIC MODIFIERS
!-------------------------------------------------------------------------------

	public :: initialize

!-------------------------------------------------------------------------------
!	MODULE PROCEDURES AND FUNCTIONS
!-------------------------------------------------------------------------------

	contains



!-------------------------------------------------------------------------------
! 						subroutine initialize					 
!																		
! Takes care of all calls necessary for the initialization of the code
! (reading input file, allocating memory, initializing data structures) 							
!-------------------------------------------------------------------------------
subroutine read_commandline_args()
  implicit none
  
  integer :: i
  character(len=32) :: arg, arg1

  do i = 1, command_argument_count()
     call get_command_argument(i,arg)
     
     select case (arg)
     case ('-i', '--input')
        call get_command_argument(i+1,arg1)
        input_file_name=trim(arg1)
     end select     
  end do
  
end subroutine read_commandline_args

subroutine initialize()

	implicit none

	! local variables
	
	integer :: ierr
	
        call read_commandline_args()

	call mkdir("output",ierr)
	call mkdir("restart",ierr)

	call set_default_values() ! sets default values for the main global variables

	call inputpar_open(input_file_name)	! parse the input file
	
	call read_input_communications()
		
	call init_communications() ! inits the MPI library and MPI data structures

        if(rank .eq. 0) print *, "Using input file ", input_file_name
	call read_input_file()		! have all modules read the required variables
	
	if(writetestlec)  call mkdir("output/tracking_elec",ierr)
	if(writetestion)  call mkdir("output/tracking_ion",ierr)
	
	if(rank.eq.0) time_beg=mpi_wtime()
	
	call initialize_restart()  ! Initializes Restart related variables
	
	call initialize_outputs()  ! Initializes variables for the outputs module
	
	call initialize_random_seed() 	!initialize random numbers on each rank

	call allocate_fields()  ! allocates memory for grid related variables (it DOES NOT SET the EM fields)

	call allocate_particles() ! allocates & initalizes particle-related variables
	
	call initialize_domain() ! Initalizes m_domain related variables
	
	if (irestart .ne. 1) call init_particle_distribution() ! sets the initial particle distribution
	
	call init_EMfields()		 ! Sets the initial EM fields
	
	if(debug)print *,rank,": before barrier"
	
	call mpi_barrier(MPI_COMM_WORLD,ierr)
	
	if(debug)print *,rank,": after barrier"
	
	if(irestart.eq.1) call restart()

	if(debug) print *,rank,":after init","step=",lap,"ions=",ions,"lecs=",lecs

	if(rank .eq. 0 ) then 
		call output_parameters_in_use()
	endif
	
	call inputpar_close()
	
end subroutine initialize



!-------------------------------------------------------------------------------
! 						subroutine set_default_values					 
!																		
! Sets default values for most of the global initialization variables in the code
! 							
!-------------------------------------------------------------------------------

subroutine set_default_values()

	implicit none
	
	time_cut=60.*15           ! in seconds, time left
	prt_first_lec=.true.
	prt_first_ion=.true.
	
	lapreorder=0
	
	mx0=100 
	my0=100+5 
	mz0=100+5 

#ifdef twoD        
	mz0=1
	betshock=.47
#else
	betshock=.33
#endif

	maxptl0=10e6 
	
	sizey=1
	
	periodicx=0
	periodicy=1
	periodicz=1
	
	highorder=.true. 
	
	corr=1.0 
	
	last=4000 
	lapst=1
	rewind(9)
	
	interval=100 
	torqint=20 !200 !100*100 !50 !100 !500 !100 !20
	
	graphics = 1
	ips= 1 
	nhist=last
	nbins=200
	nc=256
	ntimes= 4 
	
	irestart=0
	c=.45 
	sigma=0
	
	c_omp=10 
	
	gamma0=15
	delgam=1e-4 
	ppc0=4 

	btheta = 0
	bphi= 0.
	me=1.
	mi=1.
	
	writetestlec=.false.
	writetestion=.false.

	intrestlocal=-200 !negative means save no local restarts
	
!	leftclean=1*0   !need to move this to input file
!	rightclean=1*0 

end subroutine set_default_values



!-------------------------------------------------------------------------------
! 						subroutine read_input_time					 
!																		
! Reads any variables related to (or needed by) this module
! 							
!-------------------------------------------------------------------------------

subroutine read_input_time()

	implicit none

	call inputpar_geti_def("time", "last", 100, last)
	call inputpar_getd_def("time", "c", .45_sprec, c)
	call inputpar_getd_def("time", "timespan", 86400._dprec, timespan)
	

end subroutine read_input_time



!-------------------------------------------------------------------------------
! 						subroutine read_input_algorithm					 
!																		
! Reads any variables related to (or needed by) this module
! 							
!-------------------------------------------------------------------------------

subroutine read_input_algorithm()

	implicit none
	
	! local variables
	
	integer :: lhighorder
	
	call inputpar_geti_def("algorithm", "highorder", 1, lhighorder)

	if (lhighorder==1) then
		highorder=.true.
	else
		highorder=.false.
	endif

	call inputpar_getd_def("algorithm", "Corr", 1.025_sprec, corr)

	call inputpar_geti_def("algorithm", "ntimes", 4, ntimes)

end subroutine read_input_algorithm



!-------------------------------------------------------------------------------
! 						subroutine read_input_file					 
!																		
! Reads the input File
! 							
!-------------------------------------------------------------------------------

subroutine read_input_file()

	implicit none

	call read_input_time()

	call read_input_grid()
	
	call read_input_dynamic_grid()

	call read_input_algorithm()

	call read_input_restart()

	call read_input_output()

	call read_input_boundaries()

	call read_input_fields()

	call read_input_domain()

	call read_input_particles()

	call read_input_user()


end subroutine read_input_file


#ifdef twoD
end module m_initialize
#else
end module m_initialize_3d
#endif	
