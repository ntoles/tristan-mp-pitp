!
! Output module
!
! This module takes care of all the outputs from tristan
! (calls specific output procedures for each kind of output)
!
!


#ifdef twoD 

module m_output
  
	use m_system
	use m_aux
	use m_communications
	use m_fields
	use m_fieldboundaries
	use m_particles
	use m_domain
	use m_globaldata
	use m_inputparser
	use m_user
	
#else

module m_output_3d

	use m_system_3d
	use m_aux_3d
	use m_communications_3d
	use m_fields_3d
	use m_fieldboundaries_3d
	use m_particles_3d
	use m_domain_3d
	use m_globaldata_3d
	use m_inputparser_3d
	use m_user_3d
#endif

#ifdef HDF5
	use hdf5
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

	real(dprec) :: time_end, time_diff, time_beg, timespan

	integer(8) :: laprestart, namedrestartint
	integer :: intrestlocal, last, interval, torqint, istep1, graphics, ips, intrestart, &
			   pltstart, restart_status, stride, sele, seli, idx, idy, idz, istep
    ! save test particles
	integer :: dlaplec, dlapion, teststartlec, teststartion, testendlec, testendion
	! save fields frequently
			   
	! particle tracking		   
	logical :: writetestlec, writetestion
	
	character (len=5) :: indchar
	character (len=7) :: lapchar
	character (len=9) :: frootprt, frootfld, froottrq
	character (len=19) :: fnametrq
	character (len=100) :: fnameprt, fnamefld
	character (len=46) :: frestartfld, frestartprt
	character (len=34) :: frestartprtloc, frestartfldloc
	character (len=70) :: frestartprtlap
	character (len=10) :: jobname
    real(sprec), allocatable, dimension(:,:,:) :: temporary1,temporary2, &
    temporary3,temporary4,temporary5,temporary6,temporary7,temporary8, &
    temporary9,temporary10,temporary11,temporary12,temporary13,temporary14, &
    temporary15,temporary16

	real(sprec), allocatable, dimension(:,:,:) :: nav, navi
	type(prtnum), allocatable :: selprti(:), selprte(:) 
	
	real, allocatable, dimension(:,:,:) :: gam, pxm, pym, pzm

    logical :: xyzspec

!-------------------------------------------------------------------------------
!	INTERFACE DECLARATIONS
!-------------------------------------------------------------------------------
	
!-------------------------------------------------------------------------------
!	PUBLIC MODIFIERS
!-------------------------------------------------------------------------------

	! public functions
	
	public :: output_parameters_in_use, diagnostics, initialize_outputs, read_input_output

	! public variables 

	public :: frestartprtlap, frootprt, frootfld, froottrq, & 									   !used in restart module
			  last, interval, torqint, graphics, ips, c_omp, &     ! used in initialize module
			  intrestlocal, rankchar, jobname, frestartfld, frestartprt, laprestart, &  ! used in initialize module
			  lapchar, frestartprtloc, frestartfldloc, idx, idy, idz, timespan, time_beg, &        ! used in initialize module
			  intrestart, namedrestartint, pltstart, istep1, stride, istep						   ! used in initialize module

	public :: writetestlec, writetestion, dlaplec, dlapion, teststartlec, teststartion, &
			  testendlec, testendion
		

!-------------------------------------------------------------------------------
!	MODULE PROCEDURES AND FUNCTIONS
!-------------------------------------------------------------------------------

	contains



!-------------------------------------------------------------------------------
! 						subroutine read_input_output					 
!																		
! Reads any variables related to (or needed by) this module
! 							
!-------------------------------------------------------------------------------

subroutine read_input_output()

	implicit none
	
	! local variables 
	
	integer :: lwritetestion, lwritetestlec
	
	call inputpar_geti_def("output", "interval", 50, interval)
	call inputpar_geti_def("output", "torqint", 2000000, torqint)
	call inputpar_geti_def("output", "pltstart", 0, pltstart)

	call inputpar_geti_def("output", "istep", 2, istep)
	call inputpar_geti_def("output", "istep1", 2, istep1)
	call inputpar_geti_def("output", "stride", 100, stride)

	! save pre-selected test particles
	! dlap: save interval, initest: initial time step to save
	call inputpar_geti_def("output", "writetestlec", 0, lwritetestlec)
	call inputpar_geti_def("output", "dlaplec", 100, dlaplec)
	call inputpar_geti_def("output", "teststartlec", 0, teststartlec)
	call inputpar_geti_def("output", "testendlec", 0, testendlec)
	
	call inputpar_geti_def("output", "writetestion", 0, lwritetestion)
	call inputpar_geti_def("output", "dlapion", 1000, dlapion)
	call inputpar_geti_def("output", "teststartion", 0, teststartion)
	call inputpar_geti_def("output", "testendion", 0, testendion)
	
	if (lwritetestlec==1) then
		writetestlec=.true.
	else
		writetestlec=.false.
	endif
	
	if (lwritetestion==1) then
		writetestion=.true.
	else
		writetestion=.false.
	endif


end subroutine read_input_output



!-------------------------------------------------------------------------------
! 						subroutine initialize_outputs					 
!																		
! Initializes some variables used in the outputs module
!  							
!-------------------------------------------------------------------------------

subroutine initialize_outputs()

	idx=2!5 !averaging
	idy=2!5
	idz=2!5
	
#ifdef twoD
		idz=0
#endif
	
end subroutine initialize_outputs



!-------------------------------------------------------------------------------
! 						subroutine output_parameters_in_use					 
!																		
! Outputs the initial parameters read from the input file
!  							
!-------------------------------------------------------------------------------

subroutine output_parameters_in_use()

	implicit none

	print *, "Read the following parameters from the input file"
	print *, "mx0, my0, mz0", mx0, my0, mz0
	print *, "maxptl0", maxptl0
	print *, "sizex",sizex,"sizey", sizey !, "caseinit", caseinit
	print *, "c,sigma, c_omp, gamma0,delgam,ppc0,me,mi,freevar", c, &
	sigma,c_omp,  gamma0,  delgam ,ppc0,   me, mi ,dummyvar 
	print *, "irestart", irestart,"intrestart",intrestart,"laprest", &
	laprestart,"namedrestint", namedrestartint, intrestlocal
	print *,"periodicx,periodicy,periodicz",periodicx,periodicy,periodicz
	print *," wall,shiftint,shiftstart",wall,shiftint,shiftstart
	print *," Highorder", highorder
	print *, "Corr", corr
	print *, "Last", last
	print *, "interval, torqint", interval, torqint
	print *, "istep, istep1", istep, istep1
	print *, "ntimes", ntimes,"cleanfld",cleanfld,"cleanint",cleanint
	print *, "Graphics, ips", graphics, ips
	print *, "debug", debug
	print *, "btheta, bphi", btheta, bphi
	print *,"writetestlec", writetestlec
	print *,"writetestion", writetestion
	print *, "rank", rank, periodicx, radiationx, gamma0
	print *, "buffsize", buffsize

#ifdef HDF5
	print *,"hdf5 def"
#endif

	print *, "maxptl0",maxptl0
	print *, rank,":", "mz", mz, "mzall", mzall,"mz0",mz0
	print *, rank,":","mzlast",mzlast
	print *,rank,":","my", my, "myall", myall, "my0", my0
	print *, "mylast",mylast
	print *, "rank:", rank, frestartfldlap
	if(rank.eq.0) then
		print *, rank, ":", "qe=", qe
	endif

	print *,"Hello, world! I am ", rank, " of ", size0 
	print *,rank,":", frootprt, " ", frootfld," ",frestartfld
	
	if(rank .eq. 0) then
		print*, 'local mx=', mxl
	endif	
	

end subroutine output_parameters_in_use



!-------------------------------------------------------------------------------
! 						subroutine Diagnostics					 
!																		
! Manages all outputs to disk
!  							
!-------------------------------------------------------------------------------

subroutine Diagnostics()

	implicit none
	
        !user input
        ! if true then save spectra of single processors
        xyzspec=.false.
        !end user input

        ! save testparticles
		if(writetestlec) then			
			if(lap .ge. teststartlec .and. modulo((lap-pltstart),dlaplec) .eq.0 &
			   .and. lap .le. testendlec) then 					
					call write_test_e()
					prt_first_lec=.false.
			endif			
		endif
		
		if(writetestion) then			
			if(lap .ge. teststartion .and. modulo((lap-pltstart),dlapion) .eq.0 &
			   .and. lap .le. testendion) then 										
					call write_test_p()
					prt_first_ion=.false.
			endif			
		endif
					
		! add output, diagnostics, checkpointing, etc
		
		if(debug) print *,rank,": output", "step=", lap ,"step=",lap      

                call write_restart()
#ifdef HDF5
                call save_param()
                call output_tot()
!                call output_hug()
               call save_spectrum()
!              call save_spectrum_2d()
!                call save_momentumspec()
!                call save_velocityspec()
#endif

		if(debug) print *,rank,": b4 diagnost", "step=", lap

		call do_diagnostics()

end subroutine Diagnostics



!-------------------------------------------------------------------------------
! 						subroutine do_diagnostics					 
!																		
! Checks for overflow in the number of particles
!  							
!-------------------------------------------------------------------------------

subroutine do_diagnostics()

	implicit none
	
	! local variables
	
	integer ::request, add,status(MPI_STATUS_SIZE), n, ierr
	real dz
	real, allocatable :: templecsloc(:), tempionsloc(:), templecs(:,:) &
	, tempions(:,:)
	!temppar(:),tempparall(:)
	integer, allocatable :: pararr(:,:), pararrloc(:)
	integer ::ni,nl,totion,totlec, ni1, nl1, maxions, maxlecs,sizein &
	,sizeout
	integer ::maxlecsloc,maxionsloc
	integer*8 redparts(10), reduc(10)
	real norm
	
	
	redparts(1)=ions
	redparts(2)=lecs
!	redparts(3)=nionout
!	redparts(4)=nlecout
!	redparts(5)=injectedions
!	redparts(6)=injectedlecs
!	redparts(7)=receivedions
!	redparts(8)=receivedlecs
!	redparts(9)=nioneject
!	redparts(10)=nleceject
	
	call MPI_reduce(redparts,reduc,10,mpi_integer8,mpi_sum,0 &
	,MPI_Comm_World,ierr)
	
	if(rank.eq.0) print *, "lap", lap, reduc(1), reduc(2)
	
	if(ions.gt.maxhlf .or. lecs.gt.maxhlf) then 
		print *,rank,": OVERFLOW ions,lecs,maxhlf",ions,lecs,maxhlf
		stop
	endif
	
!	dz=.1

!	if(modulo(lap,torqint).eq. 0 .and. graphics .eq.1 ) then 
!		call meanq_fld_cur('tdens')
!	endif
	
end subroutine do_diagnostics

!-------------------------------------------------------------------------------
! 						subroutine save_spectrum				 
!																		
! Saves energy spectrum of particles, in x slices
!  							
!-------------------------------------------------------------------------------

subroutine save_spectrum()

	implicit none
	
	! local variables
	
	integer :: gambins, xbin, nbins, gbin,  ierr
	real gammin, gammax, gammax1, gammin1, dgam
	real dxslice, gam
	integer :: n, i, k, ind, procn
	integer :: iL, iR, jL, jR, mxmin, mxmax
	real vr, gvr, vx, vy, vz, uprime, vprime, wprime, gamprime
	character fname*20
	character(len=10) dsetname(20)

#ifdef HDF5
		integer(HID_T) :: file_id       ! File identifier 
		integer(HID_T) :: dset_id       ! Dataset identifier 
		integer(HID_T) :: dspace_id
		integer(HSIZE_T) data_dims2(2), data_dims1(1)
        integer(HSIZE_T) data_dims3(3)
#endif

	integer ::datarank, error
	
	real, allocatable:: spece(:,:), specp(:,:)
   	real, allocatable::	specein(:,:), specpin(:,:)	
	real, allocatable :: xgamma(:), xslice(:)
	
	integer :: status(MPI_STATUS_SIZE)
	
    !spectrum in the flow rest frame
	real, allocatable :: speceprime(:,:), specpprime(:,:)
	real, allocatable :: umean(:), vmean(:), wmean(:), numdens(:)	! mean flow speed, mean number density

	if(lap .ge.pltstart .and. modulo((lap-pltstart),interval) .eq.0)then
		
		gammin=1
		gammax=1
		gambins=100*2	! number of energy bins

		mxmin=3
		mxmax=mx0-2
		
		nbins=max((mxmax-mxmin)/100,1)
	
		allocate(spece(nbins,gambins),specp(nbins,gambins))
		allocate(xgamma(gambins),xslice(nbins))
		allocate(specpin(nbins,gambins),specein(nbins,gambins))

		!flow rest frame *********		
		allocate(speceprime(nbins,gambins),specpprime(nbins,gambins))
		allocate(umean(nbins),vmean(nbins),wmean(nbins))
		allocate(numdens(nbins))
		!********************		

		!bin size in each cpu. note that bin sizes are different in different cpus
		dxslice=1.*(mxmax-mxmin)/nbins
		do i=1, nbins
			xslice(i)=dxslice*(i-1)+.5*dxslice+mxmin
		enddo
		
		! local gamma min and max 
		do i=1,ions
			gam=sqrt(1+(p(i)%u**2+p(i)%v**2+p(i)%w**2)) 
			if(gam .gt. gammax) gammax=gam
			if(gam .lt. gammin) gammin=gam
			
		enddo
		do i=maxhlf+1,maxhlf+lecs
			gam=sqrt(1+(p(i)%u**2+p(i)%v**2+p(i)%w**2)) 
			if(gam .gt. gammax) gammax=gam
			if(gam .lt. gammin) gammin=gam
			
		enddo
		
		! global gamma min and max
		call mpi_allreduce(gammax,gammax1,1,mpi_read,mpi_max &
		,mpi_comm_world,ierr)
		gammax=gammax1 
		call mpi_allreduce(gammin,gammin1,1,mpi_read,mpi_min &
		,mpi_comm_world,ierr)
		gammin=gammin1
		
		gammin=max(gammin, 1.+1e-6)		
		dgam=(log10(gammax-1)-log10(gammin-1))/gambins		
			
		do i=1, gambins
			xgamma(i)=10.**((i-1.)*dgam+log10(gammin-1))
		enddo

		
		umean=0.
		vmean=0.
		wmean=0.
		numdens=0.
		specp=0.	
		specpprime=0.
			
       !(ion) flow mean velocity ***********
        do i=1,ions
		  xbin=int((p(i)%x+mxcum-mxmin)/dxslice+1)
		  if(xbin .ge. 1 .and. xbin .le. nbins) then
		    numdens(xbin)=numdens(xbin)+real(splitratio)**(1.-real(p(i &
					)%splitlev))*p(i)%ch
		    gam=sqrt(1+(p(i)%u**2+p(i)%v**2+p(i)%w**2)) 
		    umean(xbin)=umean(xbin)+(p(i)%u/gam)*real(splitratio)**(1.-real(p(i &
					)%splitlev))*p(i)%ch		
		    vmean(xbin)=vmean(xbin)+(p(i)%v/gam)*real(splitratio)**(1.-real(p(i &
					)%splitlev))*p(i)%ch		
			wmean(xbin)=wmean(xbin)+(p(i)%w/gam)*real(splitratio)**(1.-real(p(i &
					)%splitlev))*p(i)%ch				
		  endif
		enddo
		
		umean=umean/numdens
		vmean=vmean/numdens
		wmean=wmean/numdens

						
		!compute the spectrum
		if(debug) print *, "in spec, b4 ions"
		do i=1,ions
			xbin=int((p(i)%x+mxcum-mxmin)/dxslice+1)			
			if(xbin .ge. 1 .and. xbin .le. nbins) then
				gam=sqrt(1+(p(i)%u**2+p(i)%v**2+p(i)%w**2)) 								
				gbin=int((alog10(gam-1)-alog10(gammin-1))/dgam+1)
				if(gbin .ge. 1 .and. gbin .le. gambins ) then 				
					specp(xbin,gbin)=specp(xbin,gbin)+real(splitratio)**(1. &
					-real(p(i)%splitlev))*p(i)%ch                                                     
				endif
				!energy sepectrum in the flow rest frame
				!Lorentz transformation in the flow rest frame
				vx=umean(xbin)
				vy=vmean(xbin)
				vz=wmean(xbin)
				vr=sqrt(vx**2+vy**2+vz**2)
				gvr=1./sqrt(1.-vr**2)			
				uprime= -vx*gvr*gam+ &
					    (1+(gvr-1)*vx**2/vr**2)*p(i)%u+ &
						(gvr-1)*vx*vy/vr**2*p(i)%v+ &
						(gvr-1)*vx*vz/vr**2*p(i)%w
				vprime= -vy*gvr*gam+ &
						(gvr-1)*vx*vy/vr**2*p(i)%u+ &
						(1+(gvr-1)*vy**2/vr**2)*p(i)%v+ &
						(gvr-1)*vy*vz/vr**2*p(i)%w
				wprime= -vz*gvr*gam+ &
						(gvr-1)*vx*vz/vr**2*p(i)%u+ &
						(gvr-1)*vy*vz/vr**2*p(i)%v+ &
						(1+(gvr-1)*vz**2/vr**2)*p(i)%w
			    gamprime=sqrt(1+(uprime**2+vprime**2+wprime**2))
			    gbin=int((alog10(gamprime-1)-alog10(gammin-1))/dgam+1)	   
			    if(gbin .ge. 1 .and. gbin .le. gambins ) then 				
					specpprime(xbin,gbin)=specpprime(xbin,gbin)+real(splitratio)**(1. &
					-real(p(i)%splitlev))*p(i)%ch                                      
				endif
			endif !if(xbin .ge. 1 .and. xbin .le. nxbins) then
		enddo
				
		call mpi_allreduce(specp,specpin,nbins*gambins,mpi_read,mpi_sum &
		,mpi_comm_world,ierr)
        specp=specpin
		
        call mpi_allreduce(specpprime,specpin,nbins*gambins,mpi_read,mpi_sum &
		,mpi_comm_world,ierr)
        specpprime=specpin  

        do i=1, nbins
			specp(i,:)=specp(i,:)/xgamma(:)
			specpprime(i,:)=specpprime(i,:)/xgamma(:)
		enddo
        
        
		umean=0.
		vmean=0.
		wmean=0.
		numdens=0.
		spece=0.
		speceprime=0.
		
       !(electron) flow mean velocity ***********
        do i=maxhlf,maxhlf+lecs
		  xbin=int((p(i)%x+mxcum-mxmin)/dxslice+1)
		  if(xbin .ge. 1 .and. xbin .le. nbins) then
		    numdens(xbin)=numdens(xbin)+real(splitratio)**(1.-real(p(i &
					)%splitlev))*p(i)%ch
		    gam=sqrt(1+(p(i)%u**2+p(i)%v**2+p(i)%w**2)) 
		    umean(xbin)=umean(xbin)+(p(i)%u/gam)*real(splitratio)**(1.-real(p(i &
					)%splitlev))*p(i)%ch		
		    vmean(xbin)=vmean(xbin)+(p(i)%v/gam)*real(splitratio)**(1.-real(p(i &
					)%splitlev))*p(i)%ch		
			wmean(xbin)=wmean(xbin)+(p(i)%w/gam)*real(splitratio)**(1.-real(p(i &
					)%splitlev))*p(i)%ch				
		  endif
		enddo
		
		umean=umean/numdens
		vmean=vmean/numdens
		wmean=wmean/numdens        
        
        
		if(debug) print *, "in spec, b4 lecs"
		do i=maxhlf,maxhlf+lecs
			xbin=int((p(i)%x+mxcum-mxmin)/dxslice+1)
				
			if(xbin .ge. 1 .and. xbin .le. nbins) then
				gam=sqrt(1+(p(i)%u**2+p(i)%v**2+p(i)%w**2)) 								
				gbin=int((alog10(gam-1)-alog10(gammin-1))/dgam+1)
				if(gbin .ge. 1 .and. gbin .le. gambins ) then 				
					spece(xbin,gbin)=spece(xbin,gbin)+real(splitratio)**(1. &
					-real(p(i)%splitlev))*p(i)%ch                                             
				endif
				!energy sepectrum in the flow rest frame
				!Lorentz transformation in the flow rest frame
				vx=umean(xbin)
				vy=vmean(xbin)
				vz=wmean(xbin)
				vr=sqrt(vx**2+vy**2+vz**2)
				gvr=1./sqrt(1.-vr**2)
				uprime= -vx*gvr*gam+ &
					    (1+(gvr-1)*vx**2/vr**2)*p(i)%u+ &
						(gvr-1)*vx*vy/vr**2*p(i)%v+ &
						(gvr-1)*vx*vz/vr**2*p(i)%w
				vprime= -vy*gvr*gam+ &
						(gvr-1)*vx*vy/vr**2*p(i)%u+ &
						(1+(gvr-1)*vy**2/vr**2)*p(i)%v+ &
						(gvr-1)*vy*vz/vr**2*p(i)%w
				wprime= -vz*gvr*gam+ &
						(gvr-1)*vx*vz/vr**2*p(i)%u+ &
						(gvr-1)*vy*vz/vr**2*p(i)%v + &
						(1+(gvr-1)*vz**2/vr**2)*p(i)%w
			    gamprime=sqrt(1+(uprime**2+vprime**2+wprime**2))
			    gbin=int((alog10(gamprime-1)-alog10(gammin-1))/dgam+1)	   
			    if(gbin .ge. 1 .and. gbin .le. gambins ) then 				
					speceprime(xbin,gbin)=speceprime(xbin,gbin)+real(splitratio)**(1. &
					-real(p(i)%splitlev))*p(i)%ch                                           
				endif
			endif  !if(xbin .ge. 1 .and. xbin .le. nxbins) then
		enddo
		
		
        call mpi_allreduce(spece,specein,nbins*gambins,mpi_read,mpi_sum &
		,mpi_comm_world,ierr)
        spece=specein
        
        call mpi_allreduce(speceprime,specein,nbins*gambins,mpi_read,mpi_sum &
		,mpi_comm_world,ierr)
        speceprime=specein        

        do i=1, nbins
			spece(i,:)=spece(i,:)/xgamma(:)
			speceprime(i,:)=speceprime(i,:)/xgamma(:)
		enddo        
        
  	
		if(debug) print*, rank, "done spec calculation, start writing"
		
		if(rank .eq. 0) then     
			
			ind=(lap-pltstart)/interval
			
			write(indchar,"(i3.3)")ind
			
			fname="output/spect."//trim(indchar)
            print *,"",fname
			
			open(unit=11,file=fname,form='unformatted')
			close(11,status='delete')
			
			if(debug) print *,rank, ": in spec,starting save"
		
#ifdef HDF5
				!
				!  Initialize FORTRAN predefined datatypes
				!
                
				call h5open_f(error) 
				call h5fcreate_f(fname,H5F_ACC_TRUNC_F,file_id,error)
			
				dsetname(1)="spece"
				datarank=2		
				data_dims2(1)=nbins
				data_dims2(2)=gambins
				call h5screate_simple_f(datarank, data_dims2, dspace_id, error)
				call h5dcreate_f(file_id, dsetname(1), H5T_NATIVE_REAL &
				,dspace_id,dset_id,error)
				call h5dwrite_f(dset_id, H5T_NATIVE_REAL,spece, data_dims2, error)
				call h5dclose_f(dset_id,error)
				call h5sclose_f(dspace_id,error)
				
				dsetname(1)="specp"
				datarank=2				
				data_dims2(1)=nbins
				data_dims2(2)=gambins
				call h5screate_simple_f(datarank, data_dims2, dspace_id, error)
				call h5dcreate_f(file_id, dsetname(1), H5T_NATIVE_REAL &
				,dspace_id,dset_id,error)
				call h5dwrite_f(dset_id, H5T_NATIVE_REAL,specp, data_dims2, error)
				call h5dclose_f(dset_id,error)
				call h5sclose_f(dspace_id,error)
				
				!flow rest **********************
				dsetname(1)="specerest"
				datarank=2				
				data_dims2(1)=nbins
				data_dims2(2)=gambins		
				call h5screate_simple_f(datarank, data_dims2, dspace_id, error)				
				call h5dcreate_f(file_id, dsetname(1), H5T_NATIVE_REAL &
				,dspace_id,dset_id,error)				
				call h5dwrite_f(dset_id, H5T_NATIVE_REAL,speceprime, data_dims2, error)				
				call h5dclose_f(dset_id,error)
				call h5sclose_f(dspace_id,error)

				dsetname(1)="specprest"
				datarank=2				
				data_dims2(1)=nbins
				data_dims2(2)=gambins					
				call h5screate_simple_f(datarank, data_dims2, dspace_id, error)				
				call h5dcreate_f(file_id, dsetname(1), H5T_NATIVE_REAL &
				,dspace_id,dset_id,error)				
				call h5dwrite_f(dset_id, H5T_NATIVE_REAL,specpprime, data_dims2, error)				
				call h5dclose_f(dset_id,error)
				call h5sclose_f(dspace_id,error)	
				
				
				dsetname(2)="gamma"
				datarank=1			
				data_dims1(1)=gambins
				call h5screate_simple_f(datarank, data_dims1, dspace_id, error)
				call h5dcreate_f(file_id, dsetname(2), H5T_NATIVE_REAL &
				,dspace_id,dset_id,error)
				call h5dwrite_f(dset_id, H5T_NATIVE_REAL,xgamma,data_dims1,error)
				call h5dclose_f(dset_id,error)
				call h5sclose_f(dspace_id,error)
							
										
				dsetname(2)="xsl"
				datarank=1
				data_dims1(1)=nbins
				call h5screate_simple_f(datarank, data_dims1, dspace_id, error)
				call h5dcreate_f(file_id, dsetname(2), H5T_NATIVE_REAL &
				,dspace_id,dset_id,error)
				call h5dwrite_f(dset_id, H5T_NATIVE_REAL,xslice,data_dims1,error)
				call h5dclose_f(dset_id,error)
				call h5sclose_f(dspace_id,error) 
				
				dsetname(2)="gmin"
				datarank=1
				data_dims1(1)=1
				call h5screate_simple_f(datarank, data_dims1, dspace_id, error)
				call h5dcreate_f(file_id, dsetname(2), H5T_NATIVE_REAL &
				,dspace_id,dset_id,error)
				call h5dwrite_f(dset_id, H5T_NATIVE_REAL,gammin,data_dims1 &
				,error)
				
				call h5dclose_f(dset_id,error)
				call h5sclose_f(dspace_id,error) 	 
				
				dsetname(2)="gmax"
				datarank=1
				data_dims1(1)=1
				call h5screate_simple_f(datarank, data_dims1, dspace_id, error)
				call h5dcreate_f(file_id, dsetname(2), H5T_NATIVE_REAL &
				,dspace_id,dset_id,error)
				call h5dwrite_f(dset_id, H5T_NATIVE_REAL,gammax,data_dims1 &
				,error)
				call h5dclose_f(dset_id,error)
				call h5sclose_f(dspace_id,error)      
														
				
				call h5fclose_f(file_id, error)
				call h5close_f(error)


				
#endif
		endif !if rank eq 0
		
		call MPI_BARRIER(MPI_COMM_WORLD,ierr)
		
		if(allocated(xgamma)) deallocate(xgamma)
		if(allocated(xslice)) deallocate(xslice)
		if(allocated(specp)) deallocate(specp)
		if(allocated(spece)) deallocate(spece)
		
		if(allocated(specpin)) deallocate(specpin)
		if(allocated(specein)) deallocate(specein)
		
		if(allocated(speceprime)) deallocate(speceprime)		
		if(allocated(specpprime)) deallocate(specpprime)
        
		if(allocated(umean)) deallocate(umean)
		if(allocated(vmean)) deallocate(vmean)
		if(allocated(numdens)) deallocate(numdens)
	
		if(debug) print *, rank, "done spec"
		
					
	endif ! if lap is right for output
	
	
	call MPI_BARRIER(MPI_COMM_WORLD,ierr)
	
end subroutine save_spectrum

!-------------------------------------------------------------------------------
! 						subroutine save_spectrum_2d					 
!																		
! Saves energy spectrum of particles, in x slices
!  							
!-------------------------------------------------------------------------------

subroutine save_spectrum_2d()

	implicit none
	
	! local variables
	
	integer :: gambins, xbin, ybin, nxbins, nybins, gbin,  ierr
	real gammin, gammax, gammax1, gammin1, dgam
	real dxslice, dyslice, gam
	integer :: n, i, k, ind, procn
	integer :: iL, iR, jL, jR, mxmin, mxmax, mymin, mymax
	real vr, gvr, vx, vy, vz, uprime, vprime, wprime, gamprime
	character fname*20
	character(len=10) dsetname(20)

#ifdef HDF5
		integer(HID_T) :: file_id       ! File identifier 
		integer(HID_T) :: dset_id       ! Dataset identifier 
		integer(HID_T) :: dspace_id
		integer(HSIZE_T) data_dims(2), data_dims1(1)
        integer(HSIZE_T) data_dims3(3)
#endif

	integer ::datarank, error
	
	real, allocatable:: spece(:,:,:), specp(:,:,:)
   	real, allocatable::	specein(:,:,:), specpin(:,:,:)	
	real, allocatable :: xgamma(:), xslice(:), yslice(:)
	
	integer :: status(MPI_STATUS_SIZE)
	
    !spectrum in the flow rest frame
	real, allocatable :: speceprime(:,:,:), specpprime(:,:,:)!, &
!						 specelprime(:,:,:), specplprime(:,:,:)
	real, allocatable :: umean(:,:), vmean(:,:), wmean(:,:), numdens(:,:)	! mean flow speed, mean number density

	if(lap .ge.pltstart .and. modulo((lap-pltstart),interval) .eq.0)then
		
		gammin=1
		gammax=1
		gambins=100*2	! number of energy bins

		mxmin=int((mx0-5)*.5)-1200+3
		mxmax=int((mx0-5)*.5)+400+3
		mymin=3
		mymax=my0-2
		!mymin=20
		!mymax=int((my0-5)*.8)+3
		
		nxbins=max((mxmax-mxmin)/30,1)
		nybins=max((mymax-mymin)/30,1)
	
		allocate(spece(nxbins,nybins,gambins),specp(nxbins,nybins,gambins))
		allocate(xgamma(gambins),xslice(nxbins),yslice(nybins))
		allocate(specpin(nxbins,nybins,gambins),specein(nxbins,nybins,gambins))

		!flow rest frame *********		
		allocate(speceprime(nxbins,nybins,gambins),specpprime(nxbins,nybins,gambins))
!		allocate(specelprime(nxbins,nybins,gambins),specplprime(nxbins,nybins,gambins))
		allocate(umean(nxbins,nybins),vmean(nxbins,nybins),wmean(nxbins,nybins))
		allocate(numdens(nxbins,nybins))
		!********************		

		!bin size in each cpu. note that bin sizes are different in different cpus
		dxslice=1.*(mxmax-mxmin)/nxbins
		do i=1, nxbins
			xslice(i)=dxslice*(i-1)+.5*dxslice+mxmin
		enddo
		
		dyslice=1.*(mymax-mymin)/nybins
		do i=1, nybins
			yslice(i)=dyslice*(i-1)+.5*dyslice+mymin
		enddo
			 
		! local gamma min and max 
		do i=1,ions
			gam=sqrt(1+(p(i)%u**2+p(i)%v**2+p(i)%w**2)) 
			if(gam .gt. gammax) gammax=gam
			if(gam .lt. gammin) gammin=gam
			
		enddo
		do i=maxhlf+1,maxhlf+lecs
			gam=sqrt(1+(p(i)%u**2+p(i)%v**2+p(i)%w**2)) 
			if(gam .gt. gammax) gammax=gam
			if(gam .lt. gammin) gammin=gam
			
		enddo
		
		! global gamma min and max
		call mpi_allreduce(gammax,gammax1,1,mpi_read,mpi_max &
		,mpi_comm_world,ierr)
		gammax=gammax1 
		call mpi_allreduce(gammin,gammin1,1,mpi_read,mpi_min &
		,mpi_comm_world,ierr)
		gammin=gammin1
		
		gammin=max(gammin, 1.+1e-6)		
		dgam=(log10(gammax-1)-log10(gammin-1))/gambins		
			
		do i=1, gambins
			xgamma(i)=10.**((i-1.)*dgam+log10(gammin-1))
		enddo

		
		umean=0.
		vmean=0.
		wmean=0.
		numdens=0.
		specp=0.	
		specpprime=0.
!		specplprime=0.
			
       !(ion) flow mean velocity ***********
        do i=1,ions
        if(p(i)%ind .gt. 0) then
		  xbin=int((p(i)%x+mxcum-mxmin)/dxslice+1)
		  ybin=int((p(i)%y+mycum-mymin)/dyslice+1)
		  if(xbin .ge. 1 .and. xbin .le. nxbins .and. &
			 ybin .ge. 1 .and. ybin .le. nybins) then
		    numdens(xbin,ybin)=numdens(xbin,ybin)+p(i)%ch
		    gam=sqrt(1+(p(i)%u**2+p(i)%v**2+p(i)%w**2)) 
		    umean(xbin,ybin)=umean(xbin,ybin)+(p(i)%u/gam)*p(i)%ch
		    vmean(xbin,ybin)=vmean(xbin,ybin)+(p(i)%v/gam)*p(i)%ch
			wmean(xbin,ybin)=wmean(xbin,ybin)+(p(i)%w/gam)*p(i)%ch		
		  endif
		endif
		enddo
		
		umean=umean/numdens
		vmean=vmean/numdens
		wmean=wmean/numdens

						
		!compute the spectrum
		if(debug) print *, "in spec, b4 ions"
		do i=1,ions
		  if(p(i)%ind .gt. 0) then
			xbin=int((p(i)%x+mxcum-mxmin)/dxslice+1)
			ybin=int((p(i)%y+mycum-mymin)/dyslice+1)
			
			if(xbin .ge. 1 .and. xbin .le. nxbins .and. &
			   ybin .ge. 1 .and. ybin .le. nybins) then
				gam=sqrt(1+(p(i)%u**2+p(i)%v**2+p(i)%w**2)) 								
				gbin=int((alog10(gam-1)-alog10(gammin-1))/dgam+1)
!				if(gbin .ge. 1 .and. gbin .le. gambins ) then 				
!					specp(xbin,ybin,gbin)=specp(xbin,ybin,gbin)+p(i)%ch                                         
!				endif
				!energy sepectrum in the flow rest frame
				!Lorentz transformation in the flow rest frame
				vx=umean(xbin,ybin)
				vy=vmean(xbin,ybin)
				vz=wmean(xbin,ybin)
				vr=sqrt(vx**2+vy**2+vz**2)
				gvr=1./sqrt(1.-vr**2)			
				uprime= -vx*gvr*gam+ &
					    (1+(gvr-1)*vx**2/vr**2)*p(i)%u+ &
						(gvr-1)*vx*vy/vr**2*p(i)%v+ &
						(gvr-1)*vx*vz/vr**2*p(i)%w
				vprime= -vy*gvr*gam+ &
						(gvr-1)*vx*vy/vr**2*p(i)%u+ &
						(1+(gvr-1)*vy**2/vr**2)*p(i)%v+ &
						(gvr-1)*vy*vz/vr**2*p(i)%w
				wprime= -vz*gvr*gam+ &
						(gvr-1)*vx*vz/vr**2*p(i)%u+ &
						(gvr-1)*vy*vz/vr**2*p(i)%v+ &
						(1+(gvr-1)*vz**2/vr**2)*p(i)%w
			    gamprime=sqrt(1+(uprime**2+vprime**2+wprime**2))
			    gbin=int((alog10(gamprime-1)-alog10(gammin-1))/dgam+1)	   
			    if(gbin .ge. 1 .and. gbin .le. gambins ) then 				
					specpprime(xbin,ybin,gbin)=specpprime(xbin,ybin,gbin)+p(i)%ch                                       
				endif
!				if(gbin .ge. 1 .and. gbin .le. gambins .and. p(i)%ind .lt. 0) then 				
!					specplprime(xbin,ybin,gbin)=specplprime(xbin,ybin,gbin)+p(i)%ch                                       
!				endif
			endif
		  endif
		enddo
						
!        call mpi_allreduce(specp,specpin,nxbins*nybins*gambins,mpi_read,mpi_sum &
!		,mpi_comm_world,ierr)
!       specp=specpin
         
        call mpi_allreduce(specpprime,specpin,nxbins*nybins*gambins,mpi_read,mpi_sum &
		,mpi_comm_world,ierr)
        specpprime=specpin    
        
!        call mpi_allreduce(specplprime,specpin,nxbins*nybins*gambins,mpi_read,mpi_sum &
!		,mpi_comm_world,ierr)
!        specplprime=specpin  
        
		umean=0.
		vmean=0.
		wmean=0.
		numdens=0.
		spece=0.
		speceprime=0.
!		specelprime=0.
		
       !(electron) flow mean velocity ***********
        do i=maxhlf,maxhlf+lecs
        if(p(i)%ind .gt. 0) then
		  xbin=int((p(i)%x+mxcum-mxmin)/dxslice+1)
		  ybin=int((p(i)%y+mycum-mymin)/dyslice+1)
		  if(xbin .ge. 1 .and. xbin .le. nxbins .and. &
			 ybin .ge. 1 .and. ybin .le. nybins) then
		    numdens(xbin,ybin)=numdens(xbin,ybin)+p(i)%ch
		    gam=sqrt(1+(p(i)%u**2+p(i)%v**2+p(i)%w**2)) 
		    umean(xbin,ybin)=umean(xbin,ybin)+(p(i)%u/gam)*p(i)%ch
		    vmean(xbin,ybin)=vmean(xbin,ybin)+(p(i)%v/gam)*p(i)%ch
			wmean(xbin,ybin)=wmean(xbin,ybin)+(p(i)%w/gam)*p(i)%ch		
		  endif
		endif
		enddo
		
		umean=umean/numdens
		vmean=vmean/numdens
		wmean=wmean/numdens        
        
        
		if(debug) print *, "in spec, b4 lecs"
		do i=maxhlf,maxhlf+lecs
		 if(p(i)%ind .gt. 0) then
			xbin=int((p(i)%x+mxcum-mxmin)/dxslice+1)
			ybin=int((p(i)%y+mycum-mymin)/dyslice+1)
				
			if(xbin .ge. 1 .and. xbin .le. nxbins .and. &
			   ybin .ge. 1 .and. ybin .le. nybins) then
				gam=sqrt(1+(p(i)%u**2+p(i)%v**2+p(i)%w**2)) 								
				gbin=int((alog10(gam-1)-alog10(gammin-1))/dgam+1)
!				if(gbin .ge. 1 .and. gbin .le. gambins ) then 				
!					spece(xbin,ybin,gbin)=spece(xbin,ybin,gbin)+p(i)%ch                                         
!				endif
				!energy sepectrum in the flow rest frame
				!Lorentz transformation in the flow rest frame
				vx=umean(xbin,ybin)
				vy=vmean(xbin,ybin)
				vz=wmean(xbin,ybin)
				vr=sqrt(vx**2+vy**2+vz**2)
				gvr=1./sqrt(1.-vr**2)
				uprime= -vx*gvr*gam+ &
					    (1+(gvr-1)*vx**2/vr**2)*p(i)%u+ &
						(gvr-1)*vx*vy/vr**2*p(i)%v+ &
						(gvr-1)*vx*vz/vr**2*p(i)%w
				vprime= -vy*gvr*gam+ &
						(gvr-1)*vx*vy/vr**2*p(i)%u+ &
						(1+(gvr-1)*vy**2/vr**2)*p(i)%v+ &
						(gvr-1)*vy*vz/vr**2*p(i)%w
				wprime= -vz*gvr*gam+ &
						(gvr-1)*vx*vz/vr**2*p(i)%u+ &
						(gvr-1)*vy*vz/vr**2*p(i)%v + &
						(1+(gvr-1)*vz**2/vr**2)*p(i)%w
			    gamprime=sqrt(1+(uprime**2+vprime**2+wprime**2))
			    gbin=int((alog10(gamprime-1)-alog10(gammin-1))/dgam+1)	   
			    if(gbin .ge. 1 .and. gbin .le. gambins ) then 				
					speceprime(xbin,ybin,gbin)=speceprime(xbin,ybin,gbin)+p(i)%ch                                       
				endif
!				if(gbin .ge. 1 .and. gbin .le. gambins .and. p(i)%ind .lt. 0) then 				
!					specelprime(xbin,ybin,gbin)=specelprime(xbin,ybin,gbin)+p(i)%ch                                       
!				endif				
			endif
		  endif
		enddo
		
		
!        call mpi_allreduce(spece,specein,nxbins*nybins*gambins,mpi_read,mpi_sum &
!		,mpi_comm_world,ierr)
 !       spece=specein
        
        call mpi_allreduce(speceprime,specein,nxbins*nybins*gambins,mpi_read,mpi_sum &
		,mpi_comm_world,ierr)
        speceprime=specein        

!        call mpi_allreduce(specelprime,specein,nxbins*nybins*gambins,mpi_read,mpi_sum &
!		,mpi_comm_world,ierr)
!       specelprime=specein
        
  	
		if(debug) print*, rank, "done spec calculation, start writing"
		
		if(rank .eq. 0) then     
			
			ind=(lap-pltstart)/interval
			
			write(indchar,"(i3.3)")ind
			
			fname="output/spect."//trim(indchar)
            print *,"",fname
			
			open(unit=11,file=fname,form='unformatted')
			close(11,status='delete')
			
			if(debug) print *,rank, ": in spec,starting save"
		
#ifdef HDF5
				!
				!  Initialize FORTRAN predefined datatypes
				!
                
				call h5open_f(error) 
				call h5fcreate_f(fname,H5F_ACC_TRUNC_F,file_id,error)
			
!				dsetname(1)="spece"
!				datarank=3		
!				data_dims3(1)=nxbins
!				data_dims3(2)=nybins
!				data_dims3(3)=gambins
!				call h5screate_simple_f(datarank, data_dims3, dspace_id, error)
!				call h5dcreate_f(file_id, dsetname(1), H5T_NATIVE_REAL &
!				,dspace_id,dset_id,error)
!				call h5dwrite_f(dset_id, H5T_NATIVE_REAL,spece, data_dims3, error)
!				call h5dclose_f(dset_id,error)
!				call h5sclose_f(dspace_id,error)
				
!				dsetname(1)="specp"
!				datarank=3				
!				data_dims3(1)=nxbins
!				data_dims3(2)=nybins
!				data_dims3(3)=gambins
!				call h5screate_simple_f(datarank, data_dims3, dspace_id, error)
!				call h5dcreate_f(file_id, dsetname(1), H5T_NATIVE_REAL &
!				,dspace_id,dset_id,error)
!				call h5dwrite_f(dset_id, H5T_NATIVE_REAL,specp, data_dims3, error)
!				call h5dclose_f(dset_id,error)
!				call h5sclose_f(dspace_id,error)
				
				!flow rest **********************
				dsetname(1)="specerest"
				datarank=3				
				data_dims3(1)=nxbins
				data_dims3(2)=nybins
				data_dims3(3)=gambins		
				call h5screate_simple_f(datarank, data_dims3, dspace_id, error)				
				call h5dcreate_f(file_id, dsetname(1), H5T_NATIVE_REAL &
				,dspace_id,dset_id,error)				
				call h5dwrite_f(dset_id, H5T_NATIVE_REAL,speceprime, data_dims3, error)				
				call h5dclose_f(dset_id,error)
				call h5sclose_f(dspace_id,error)

				dsetname(1)="specprest"
				datarank=3				
				data_dims3(1)=nxbins
				data_dims3(2)=nybins
				data_dims3(3)=gambins					
				call h5screate_simple_f(datarank, data_dims3, dspace_id, error)				
				call h5dcreate_f(file_id, dsetname(1), H5T_NATIVE_REAL &
				,dspace_id,dset_id,error)				
				call h5dwrite_f(dset_id, H5T_NATIVE_REAL,specpprime, data_dims3, error)				
				call h5dclose_f(dset_id,error)
				call h5sclose_f(dspace_id,error)	
				
!				dsetname(1)="specelrest"
!				datarank=3				
!				data_dims3(1)=nxbins
!				data_dims3(2)=nybins
!				data_dims3(3)=gambins		
!				call h5screate_simple_f(datarank, data_dims3, dspace_id, error)				
!				call h5dcreate_f(file_id, dsetname(1), H5T_NATIVE_REAL &
!				,dspace_id,dset_id,error)				
!				call h5dwrite_f(dset_id, H5T_NATIVE_REAL,specelprime, data_dims3, error)				
!				call h5dclose_f(dset_id,error)
!				call h5sclose_f(dspace_id,error)

!				dsetname(1)="specplrest"
!				datarank=3				
!				data_dims3(1)=nxbins
!				data_dims3(2)=nybins
!				data_dims3(3)=gambins					
!				call h5screate_simple_f(datarank, data_dims3, dspace_id, error)				
!				call h5dcreate_f(file_id, dsetname(1), H5T_NATIVE_REAL &
!				,dspace_id,dset_id,error)				
!				call h5dwrite_f(dset_id, H5T_NATIVE_REAL,specplprime, data_dims3, error)				
!				call h5dclose_f(dset_id,error)
!				call h5sclose_f(dspace_id,error)
									
				dsetname(2)="gamma"
				datarank=1			
				data_dims1(1)=gambins
				call h5screate_simple_f(datarank, data_dims1, dspace_id, error)
				call h5dcreate_f(file_id, dsetname(2), H5T_NATIVE_REAL &
				,dspace_id,dset_id,error)
				call h5dwrite_f(dset_id, H5T_NATIVE_REAL,xgamma,data_dims1,error)
				call h5dclose_f(dset_id,error)
				call h5sclose_f(dspace_id,error)
							
										
				dsetname(2)="xsl"
				datarank=1
				data_dims1(1)=nxbins
				call h5screate_simple_f(datarank, data_dims1, dspace_id, error)
				call h5dcreate_f(file_id, dsetname(2), H5T_NATIVE_REAL &
				,dspace_id,dset_id,error)
				call h5dwrite_f(dset_id, H5T_NATIVE_REAL,xslice,data_dims1,error)
				call h5dclose_f(dset_id,error)
				call h5sclose_f(dspace_id,error) 
				
				dsetname(2)="ysl"
				datarank=1
				data_dims1(1)=nybins
				call h5screate_simple_f(datarank, data_dims1, dspace_id, error)
				call h5dcreate_f(file_id, dsetname(2), H5T_NATIVE_REAL &
				,dspace_id,dset_id,error)
				call h5dwrite_f(dset_id, H5T_NATIVE_REAL,yslice,data_dims1,error)
				call h5dclose_f(dset_id,error)
				call h5sclose_f(dspace_id,error) 
							 
				
				dsetname(2)="gmax"
				datarank=1
				data_dims1(1)=1
				call h5screate_simple_f(datarank, data_dims1, dspace_id, error)
				call h5dcreate_f(file_id, dsetname(2), H5T_NATIVE_REAL &
				,dspace_id,dset_id,error)
				call h5dwrite_f(dset_id, H5T_NATIVE_REAL,gammax,data_dims1 &
				,error)
				call h5dclose_f(dset_id,error)
				call h5sclose_f(dspace_id,error)      
														
				
				call h5fclose_f(file_id, error)
				call h5close_f(error)


				
#endif
		endif !if rank eq 0
		
		call MPI_BARRIER(MPI_COMM_WORLD,ierr)
		
		if(allocated(xgamma)) deallocate(xgamma)
		if(allocated(xslice)) deallocate(xslice)
		if(allocated(yslice)) deallocate(yslice)
		if(allocated(specp)) deallocate(specp)
		if(allocated(spece)) deallocate(spece)
		
		if(allocated(specpin)) deallocate(specpin)
		if(allocated(specein)) deallocate(specein)
		
		if(allocated(speceprime)) deallocate(speceprime)		
		if(allocated(specpprime)) deallocate(specpprime)
!		if(allocated(specelprime)) deallocate(specelprime)		
!		if(allocated(specplprime)) deallocate(specplprime)
        
		if(allocated(umean)) deallocate(umean)
		if(allocated(vmean)) deallocate(vmean)
		if(allocated(numdens)) deallocate(numdens)
	
		if(debug) print *, rank, "done spec"
		
					
	endif ! if lap is right for output
	
	
	call MPI_BARRIER(MPI_COMM_WORLD,ierr)
	
end subroutine save_spectrum_2d


!-------------------------------------------------------------------------------
! 						subroutine save_momentumspec					 
!																		
! Saves momentum spectrum of particles, in x slices
!  							
!-------------------------------------------------------------------------------

subroutine save_momentumspec()

  implicit none

  real ishockall, ishock, dxslice
  real maxdens, mindens, middens
  integer mombins, xbin, mbin, mbinneg, ierr
  real mommin, mommax, mommin1, mommax1, dmom, mom
  logical logscale
  integer n, i, ind, var, procn
  character fname*80
  character(len=20) dsetname(20)
  
#ifdef HDF5
  INTEGER(HID_T) :: file_id ! File identifier 
  INTEGER(HID_T) :: dset_id ! Dataset identifier 
  integer(HID_T) :: dspace_id
  INTEGER(HSIZE_T) data_dims(2), data_dims1(1)
  integer(HSIZE_T) data_dims3(3)
#endif
  integer datarank, error
  
  real, allocatable :: spece(:,:), specp(:,:), specin(:,:)
  real, allocatable :: speceb(:,:), specpb(:,:)
  real, allocatable :: logmom(:), xslice(:)
  real, allocatable :: dens1(:), densin(:),  dens(:), dens2(:)

  ! spectra sliced in x, y, z (each processor saves one) 
  real, allocatable:: tspece(:,:,:), tspecp(:,:,:)
  real, allocatable:: tspeceb(:,:,:), tspecpb(:,:,:)
  integer :: status(MPI_STATUS_SIZE)
  
  logscale=.true.
  mombins=100*2
  dxslice=100
  
  IF(lap.ge.pltstart.and.modulo((lap-pltstart),interval).eq.0) THEN
     
     do var=1,3             !loop over dimensions
        
        select case (var)
        case (1)
           if(rank .eq. 0.and.debug) PRINT *, "Doing x momenta"
        case (2)
           if(rank .eq. 0.and.debug) PRINT *, "Doing y momenta"
        case (3)
           if(rank .eq. 0.and.debug) PRINT *, "Doing z momenta"
        end select
        
        !     Find momenta with maximum and minimum magnitudes            
        mommax=0.
        mommin=1e8
        
        do i=1,ions
           select case (var)
           case (1) 
              mom=p(i)%u
           case (2)
              mom=p(i)%v
           case (3)
              mom=p(i)%w
           end select
           
           if(abs(mom).gt. mommax) mommax=abs(mom)
           if(abs(mom).lt. mommin) mommin=abs(mom)
        enddo
        
        do i=maxhlf+1,maxhlf+lecs
           select case (var)
           case (1) 
              mom=p(i)%u
           case (2)
              mom=p(i)%v
           case (3)
              mom=p(i)%w
           end select
           
           if(abs(mom).gt. mommax) mommax=abs(mom)
           if(abs(mom).lt. mommin) mommin=abs(mom)
        enddo
        
        
        !     Compare vals from each processor and determine global maxima and minima.
        call mpi_allreduce(mommax,mommax1,1,mpi_read,mpi_max &
             ,mpi_comm_world,ierr)
        mommax=mommax1
        call mpi_allreduce(mommin,mommin1,1,mpi_read,mpi_min &
             ,mpi_comm_world,ierr)
        mommin=mommin1
        
        !     checks
        IF(mommin.eq.0.and.mommax.eq.0) THEN
           if(rank .eq. 0.and.debug) PRINT *, "No interesting momenta"
        ENDIF
        mommin=max(mommin,10.d-5)
        IF(mommax.eq.0.)THEN
           mommax=10.d-4
        ENDIF
        
        !     Compute logarithmic binsize.      
        IF(logscale)dmom=(log10(mommax)-log10(mommin))/(mombins/2.)  
        IF(dmom.eq.0)THEN
           dmom=10.**(-5.)
        ENDIF
        
        !     Allocate appropriate storage for logarithmic histogram.
        allocate(spece(max(int(mx0/dxslice),1),mombins))
        allocate(specp(max(int(mx0/dxslice),1),mombins))
        allocate(speceb(max(int(mx0/dxslice),1),mombins))
        allocate(specpb(max(int(mx0/dxslice),1),mombins))
        allocate(specin(max(int(mx0/dxslice),1),mombins))
        allocate(logmom(mombins))
        allocate(xslice(max(int(mx0/dxslice),1)))
        allocate(dens(max(int(mx0/10),1)))
        allocate(dens2(max(int(mx0/10),1)))
        allocate(dens1(max(int(mx0/dxslice),1)))
        allocate(densin(max(int(mx0/dxslice),1)))

        if (xyzspec) then
           allocate(tspece(max(int(mx0/dxslice),1),mombins,size0))
           allocate(tspecp(max(int(mx0/dxslice),1),mombins,size0))
           allocate(tspeceb(max(int(mx0/dxslice),1),mombins,size0))
           allocate(tspecpb(max(int(mx0/dxslice),1),mombins,size0))
        endif

        !     Compute values of logarithmic bins
        do i=1,mombins/2
           !               logmom(i+mombins/2)=10.**((i-1)*dmom)*mommin
           !               logmom(mombins/2+1-i)=-(10.**((i-1)*dmom)*mommin)
           logmom(i+mombins/2)=10.**((i-0.5)*dmom)*mommin
           logmom(mombins/2+1-i)=-(10.**((i-0.5)*dmom)*mommin)
        enddo
        
        !     Set default spectra to 0. 
        spece(:,:)=0.
        specp(:,:)=0.
        speceb(:,:)=0.
        specpb(:,:)=0.
        specin(:,:)=0.
        dens1(:)=0.
        densin(:)=0.
        dens(:)=0.
        dens2(:)=0.

        if (xyzspec) then
           tspece(:,:,:)=0.
           tspecp(:,:,:)=0.
           tspeceb(:,:,:)=0.
           tspecpb(:,:,:)=0.
        endif

        !     find where is the shock
        maxdens=0
        mindens=1e8
        do n=1,ions
           i=int(p(n)%x)               
           dens(max(min(int(i/10.+1),mx0/10),1))= &
                dens(max(min(int(i/10.+1),mx0/10),1))+p(n)%ch
        enddo
        call mpi_allreduce(dens,dens2,max(int(mx0/10),1),mpi_read,mpi_sum &
             ,mpi_comm_world,ierr)            
        dens=dens2/size0
        
        maxdens=maxval(dens)
        do i=1,max(int(mx0/10),1)
           if(dens(i) .gt. maxdens/10.) then
              if(dens(i).lt.mindens)mindens=dens(i)
           endif
        enddo
        middens=minval(abs(dens-(maxdens+mindens)/2.))            
        do i=1,max(int(mx0/10),1)
           if(abs(dens(i)-(maxdens+mindens)/2.) .eq. middens) then
              ishock=i*10.-5
           endif
        enddo
        
        !     define spatial bins
        do i=1,max(int(mx0/dxslice),1)
           xslice(i)=(i-1)*dxslice+dxslice/2.
        enddo
        
        if(debug) print *, rank, ": spec calculation"
        
        !     Compute the spectrum for ions
        do i=1,ions
           
           !     Compute the x bin
           xbin=int(p(i)%x/dxslice+1)
           
           !     Compute momentum bin
           mbin=mombins/2+1
           if(xbin .ge. 1 .and. xbin .le. max(int(mx/dxslice),1)) then
              select case (var)
              case (1) 
                 mom=p(i)%u
              case (2)
                 mom=p(i)%v
              case (3)
                 mom=p(i)%w
              end select
              
              !     Boost momentum within limits (logscale)
              if(abs(mom).lt.mommin)then
                 if(mom.ge.0)then
                    mom=mommin 
                 else
                    mom=-mommin
                 endif
              endif
              
              !     Compute the logarithmically binned spectra.
              if(mom.lt.0)then
                 if(logscale) mbinneg=int((alog10(abs(mom))- &
                      alog10(mommin))/dmom+1)
                 mbin=mombins/2-mbinneg+1
              endif
              if(mom.gt.0)then
                 if(logscale) mbin=int((alog10(abs(mom))- &
                      alog10(mommin))/dmom+1)+mombins/2
              endif
              
              if(mbin .ge. 1 .and. mbin .le. mombins ) then 
                 dens1(xbin)=dens1(xbin)+ p(i)%ch
                 specp(xbin,mbin)=specp(xbin,mbin)+ p(i)%ch

                 if (p(i)%ind .lt. 0) then 
                    specpb(xbin,mbin)=specpb(xbin,mbin)+ p(i)%ch
                 endif

              endif
           endif !xbin
        enddo !over ions
        
        !     Compute the spectrum for electrons.
        do i=maxhlf+1,maxhlf+lecs

           !     Compute the x bin
           xbin=int(p(i)%x/dxslice+1)
           
           !     Compute momentum bin
           mbin=mombins/2+1
           if(xbin .ge. 1 .and. xbin .le. max(int(mx/dxslice),1)) then
              select case (var)
              case (1) 
                 mom=p(i)%u
              case (2)
                 mom=p(i)%v
              case (3)
                 mom=p(i)%w
              end select
              
              !     Boost momentum within limits
              if(abs(mom).lt.mommin)then
                 if(mom.ge.0)then
                    mom=mommin 
                 else
                    mom=-mommin
                 endif
              endif
              
              !     Compute the logarithmic bin.
              if(mom.lt.0) then
                 if(logscale) mbinneg=int((alog10(abs(mom))- &
                      alog10(mommin))/dmom+1)
                 mbin=mombins/2-mbinneg+1
              endif
              if(mom.gt.0)then
                 if(logscale) mbin=int((alog10(abs(mom))- &
                      alog10(mommin))/dmom+1)+mombins/2 
              endif
              
              !     Compute the logarithmic spectra.
              if(mbin .ge. 1 .and. mbin .le. mombins ) then 
                 dens1(xbin)=dens1(xbin)+ p(i)%ch
                 spece(xbin,mbin)=spece(xbin,mbin)+ p(i)%ch

                 if (p(i)%ind .lt. 0) then 
                    speceb(xbin,mbin)=speceb(xbin,mbin)+ p(i)%ch
                 endif
              endif
           endif  !xbin       
           
        enddo !over lecs
        
        call MPI_BARRIER(MPI_COMM_WORLD,ierr)
        if(debug) print *, rank, "reducing spectra"
        
        !     now reduce all to rank 0, which will save everything
        ! xyz spectra
        if (xyzspec) then
           !lecs	
           do procn=0,size0-1                   
              if (rank .eq. 0) then                      
                 if (procn .ne. 0) then !receive from procn                                
                    call MPI_Recv(specin,max(int(mx0/dxslice),1)*mombins &
                         ,mpi_read,procn,1,MPI_COMM_WORLD,status,error)
                    tspece(:,:,procn+1)=specin(:,:)
                    specin(:,:)=0.
                 else           !procn eq 0
                    tspece(:,:,1)=spece(:,:)               
                 endif
              else !rank ne 0                      
                 if (procn .eq. rank) then
                    call MPI_Send(spece,max(int(mx0/dxslice),1)*mombins &
                         ,mpi_read,0,1,MPI_COMM_WORLD,error)          
                 endif
              endif !rank ne 0                   
           enddo!procn     
           !ions
           do procn=0,size0-1                   
              if (rank .eq. 0) then                      
                 if (procn .ne. 0) then !receive from procn                                
                    call MPI_Recv(specin,max(int(mx0/dxslice),1)*mombins &
                         ,mpi_read,procn,1,MPI_COMM_WORLD,status,error)
                    tspecp(:,:,procn+1)=specin(:,:)
                    specin(:,:)=0.
                 else           !procn eq 0
                    tspecp(:,:,1)=specp(:,:)               
                 endif
              else !rank ne 0                      
                 if (procn .eq. rank) then
                    call MPI_Send(specp,max(int(mx0/dxslice),1)*mombins &
                         ,mpi_read,0,1,MPI_COMM_WORLD,error)          
                 endif
              endif !rank ne 0                   
           enddo!procn 
           !beam lecs
           do procn=0,size0-1                   
              if (rank .eq. 0) then                      
                 if (procn .ne. 0) then !receive from procn                                
                    call MPI_Recv(specin,max(int(mx0/dxslice),1)*mombins &
                         ,mpi_read,procn,1,MPI_COMM_WORLD,status,error)
                    tspeceb(:,:,procn+1)=specin(:,:)
                    specin(:,:)=0.
                 else           !procn eq 0
                    tspeceb(:,:,1)=speceb(:,:)               
                 endif
              else !rank ne 0                      
                 if (procn .eq. rank) then
                    call MPI_Send(speceb,max(int(mx0/dxslice),1)*mombins &
                         ,mpi_read,0,1,MPI_COMM_WORLD,error)          
                 endif
              endif !rank ne 0                   
           enddo!procn
           ! beam ions
           do procn=0,size0-1                   
              if (rank .eq. 0) then                      
                 if (procn .ne. 0) then !receive from procn                                
                    call MPI_Recv(specin,max(int(mx0/dxslice),1)*mombins &
                         ,mpi_read,procn,1,MPI_COMM_WORLD,status,error)
                    tspecpb(:,:,procn+1)=specin(:,:)
                    specin(:,:)=0.
                 else           !procn eq 0
                    tspecpb(:,:,1)=specpb(:,:)               
                 endif
              else !rank ne 0                      
                 if (procn .eq. rank) then
                    call MPI_Send(specpb,max(int(mx0/dxslice),1)*mombins &
                         ,mpi_read,0,1,MPI_COMM_WORLD,error)          
                 endif
              endif !rank ne 0                   
           enddo!procn
        endif
        ! end xyz spectra
        

        !     Sum the logarithmic spectra and make the final spectra global.
        call mpi_allreduce(spece,specin,max(int(mx0/dxslice),1) &
             *mombins,mpi_read,mpi_sum,mpi_comm_world,ierr)
        
        spece=specin
        specin(:,:)=0.
        
        call mpi_allreduce(specp,specin,max(int(mx0/dxslice),1) &
             *mombins,mpi_read,mpi_sum,mpi_comm_world,ierr)
        
        specp=specin
        specin(:,:)=0.
        
        call mpi_allreduce(speceb,specin,max(int(mx0/dxslice),1) &
             *mombins,mpi_read,mpi_sum,mpi_comm_world,ierr)
        
        speceb=specin
        specin(:,:)=0.
        
        call mpi_allreduce(specpb,specin,max(int(mx0/dxslice),1) &
             *mombins,mpi_read,mpi_sum,mpi_comm_world,ierr)
        
        specpb=specin
        
        call mpi_allreduce(dens1,densin,max(int(mx0/dxslice),1) &
             ,mpi_read,mpi_sum,mpi_comm_world,ierr)
        
        dens1=densin/size0
        
        !     Convert from dN/dlogp to dN/dp.
        do i=1,max(int(mx0/dxslice),1)
           spece(i,:)=spece(i,:)/abs(logmom)
           specp(i,:)=specp(i,:)/abs(logmom)
           speceb(i,:)=speceb(i,:)/abs(logmom)
           specpb(i,:)=specpb(i,:)/abs(logmom)
           if (xyzspec) then
              do procn=0,size0-1
                 tspece(i,:,procn+1)=tspece(i,:,procn+1)/abs(logmom)
                 tspecp(i,:,procn+1)=tspecp(i,:,procn+1)/abs(logmom)
                 tspeceb(i,:,procn+1)=tspeceb(i,:,procn+1)/abs(logmom)
                 tspecpb(i,:,procn+1)=tspecpb(i,:,procn+1)/abs(logmom)
              enddo
           endif
        enddo
        
        call MPI_BARRIER(MPI_COMM_WORLD,ierr)
        
        if(debug) print *, rank, "done spec calc,now writing"
        
        if(rank .eq. 0) then 
           ind=(lap-pltstart)/interval
           write(indchar,"(i3.3)")ind               
           fname="output/momentum."//trim(indchar)               
           if(var.eq.1) then 
              if(rank.eq.0)PRINT *, "", fname
              open(unit=11,file=fname,form='unformatted')
              close(11,status='delete')
           endif
#ifdef HDF5
           
           !     Initialize FORTRAN predefined datatypes.               
           !     Create storage file to save spectra to.               
           IF(var.eq.1)THEN 
              CALL h5open_f(error) 
              CALL h5fcreate_f(fname,H5F_ACC_TRUNC_F,file_id,error)
           ENDIF
           
           !xyz spectra
           if (xyzspec) then
              select case (var)
              case (1)
                 dsetname(1)="tpxelogsp"
              case (2)
                 dsetname(1)="tpyelogsp"
              case (3)
                 dsetname(1)="tpzelogsp"
              end select
              
              datarank=3				
              data_dims3(1)=max(int(mx0/dxslice),1)
              data_dims3(2)=mombins
              data_dims3(3)=size0
              
              call h5screate_simple_f(datarank, data_dims3, dspace_id, error)
              
              call h5dcreate_f(file_id, dsetname(1), H5T_NATIVE_REAL &
                   ,dspace_id,dset_id,error)
              
              call h5dwrite_f(dset_id, H5T_NATIVE_REAL,tspece, data_dims3, error)
              
              call h5dclose_f(dset_id,error)
              call h5sclose_f(dspace_id,error)
              
              select case (var)
              case (1)
                 dsetname(1)="tpxplogsp"
              case (2)
                 dsetname(1)="tpyplogsp"
              case (3)
                 dsetname(1)="tpzplogsp"
              end select
              
              datarank=3				
              data_dims3(1)=max(int(mx0/dxslice),1)
              data_dims3(2)=mombins
              data_dims3(3)=size0
              
              call h5screate_simple_f(datarank, data_dims3, dspace_id, error)
              
              call h5dcreate_f(file_id, dsetname(1), H5T_NATIVE_REAL &
                   ,dspace_id,dset_id,error)
              
              call h5dwrite_f(dset_id, H5T_NATIVE_REAL,tspecp, data_dims3, error)
              
              call h5dclose_f(dset_id,error)
              call h5sclose_f(dspace_id,error)
              
              select case (var)
              case (1)
                 dsetname(1)="tpxeblogsp"
              case (2)
                 dsetname(1)="tpyeblogsp"
              case (3)
                 dsetname(1)="tpzeblogsp"
              end select
              
              datarank=3				
              data_dims3(1)=max(int(mx0/dxslice),1)
              data_dims3(2)=mombins
              data_dims3(3)=size0
				
              call h5screate_simple_f(datarank, data_dims3, dspace_id, error)
              
              call h5dcreate_f(file_id, dsetname(1), H5T_NATIVE_REAL &
                   ,dspace_id,dset_id,error)
              
              call h5dwrite_f(dset_id, H5T_NATIVE_REAL,tspeceb, data_dims3, error)
              
              call h5dclose_f(dset_id,error)
              call h5sclose_f(dspace_id,error)
              
              select case (var)
              case (1)
                 dsetname(1)="tpxpblogsp"
              case (2)
                 dsetname(1)="tpypblogsp"
              case (3)
                 dsetname(1)="tpzpblogsp"
              end select
              
              datarank=3				
              data_dims3(1)=max(int(mx0/dxslice),1)
              data_dims3(2)=mombins
              data_dims3(3)=size0
              
              call h5screate_simple_f(datarank, data_dims3, dspace_id, error)
              
              call h5dcreate_f(file_id, dsetname(1), H5T_NATIVE_REAL &
                   ,dspace_id,dset_id,error)
              
              call h5dwrite_f(dset_id, H5T_NATIVE_REAL,tspecpb, data_dims3, error)
              
              call h5dclose_f(dset_id,error)
              call h5sclose_f(dspace_id,error)
           endif
           ! end xyz spectrum
           
           !     Write logarithmic spectra.
           select case (var)
           case (1)
              dsetname(1)="pxelogsp"
           case (2)
              dsetname(1)="pyelogsp"
           case (3)
              dsetname(1)="pzelogsp"
           end select
           
           datarank=2          
           data_dims(1)=max(int(mx0/dxslice),1)
           data_dims(2)=mombins
           CALL h5screate_simple_f(datarank,data_dims,dspace_id, &
                error)               
           call h5dcreate_f(file_id, dsetname(1), H5T_NATIVE_REAL &
                ,dspace_id,dset_id,error)
           CALL h5dwrite_f(dset_id, H5T_NATIVE_REAL,spece,data_dims &
                , error)
           call h5dclose_f(dset_id,error)
           call h5sclose_f(dspace_id,error)
           
           select case (var)
           case (1)
              dsetname(1)="pxplogsp"
           case (2)
              dsetname(1)="pyplogsp"
           case (3)
              dsetname(1)="pzplogsp"
           end select
           
           datarank=2               
           data_dims(1)=max(int(mx0/dxslice),1)
           data_dims(2)=mombins               
           CALL h5screate_simple_f(datarank, data_dims, dspace_id,&
                error)
           call h5dcreate_f(file_id, dsetname(1), H5T_NATIVE_REAL &
                ,dspace_id,dset_id,error)               
           CALL h5dwrite_f(dset_id, H5T_NATIVE_REAL,specp,data_dims &
                , error)               
           call h5dclose_f(dset_id,error)
           call h5sclose_f(dspace_id,error)      

           select case (var)
           case (1)
              dsetname(1)="pxpblogsp"
           case (2)
              dsetname(1)="pypblogsp"
           case (3)
              dsetname(1)="pzpblogsp"
           end select
           
           datarank=2               
           data_dims(1)=max(int(mx0/dxslice),1)
           data_dims(2)=mombins               
           CALL h5screate_simple_f(datarank, data_dims, dspace_id, &
                error)
           call h5dcreate_f(file_id, dsetname(1), H5T_NATIVE_REAL &
                ,dspace_id,dset_id,error)               
           CALL h5dwrite_f(dset_id, H5T_NATIVE_REAL,specpb,data_dims &
                , error)               
           call h5dclose_f(dset_id,error)
           call h5sclose_f(dspace_id,error) 
           
           select case (var)
           case (1)
              dsetname(1)="pxeblogsp"
           case (2)
              dsetname(1)="pyeblogsp"
           case (3)
              dsetname(1)="pzeblogsp"
           end select
           
           datarank=2               
           data_dims(1)=max(int(mx0/dxslice),1)
           data_dims(2)=mombins               
           CALL h5screate_simple_f(datarank, data_dims, dspace_id, &
                error)
           call h5dcreate_f(file_id, dsetname(1), H5T_NATIVE_REAL &
                ,dspace_id,dset_id,error)               
           CALL h5dwrite_f(dset_id, H5T_NATIVE_REAL,speceb,data_dims &
                , error)               
           call h5dclose_f(dset_id,error)
           call h5sclose_f(dspace_id,error)
           
           select case (var)
           case (1)
              dsetname(2)="pxbin"
           case (2)
              dsetname(2)="pybin"
           case (3)
              dsetname(2)="pzbin"
           end select
           
           datarank=1              
           data_dims1(1)=mombins
           CALL h5screate_simple_f(datarank, data_dims1,dspace_id, &
                error)               
           call h5dcreate_f(file_id, dsetname(2), H5T_NATIVE_REAL &
                ,dspace_id,dset_id,error)               
           CALL h5dwrite_f(dset_id, H5T_NATIVE_REAL,logmom, &
                data_dims1,error)               
           call h5dclose_f(dset_id,error)
           call h5sclose_f(dspace_id,error)      
           
           select case (var)
           case (1)
              dsetname(1)="pxlim"
           case (2)
              dsetname(1)="pylim"
           case (3)
              dsetname(1)="pzlim"
           end select
           
           datarank=1               
           data_dims1(1)=2               
           CALL h5screate_simple_f(datarank, data_dims1,dspace_id, & 
                error)               
           call h5dcreate_f(file_id, dsetname(1), H5T_NATIVE_REAL &
                ,dspace_id,dset_id,error)               
           CALL h5dwrite_f(dset_id,H5T_NATIVE_REAL,[mommin, mommax], &
                data_dims1, error)               
           call h5dclose_f(dset_id,error)
           call h5sclose_f(dspace_id,error)
           
           select case (var)
           case (1)
              dsetname(2)="dpx"
           case (2)
              dsetname(2)="dpy"
           case (3)
              dsetname(2)="dpz"
           end select
           
           datarank=1               
           data_dims1(1)=1
           CALL h5screate_simple_f(datarank, data_dims1, dspace_id, &
                error)               
           call h5dcreate_f(file_id, dsetname(2), H5T_NATIVE_REAL &
                ,dspace_id,dset_id,error)               
           CALL h5dwrite_f(dset_id, H5T_NATIVE_REAL,dmom,data_dims1 &
                ,error)               
           call h5dclose_f(dset_id,error)
           call h5sclose_f(dspace_id,error)      
           
           if(var .eq. 1) then
              dsetname(2)="xshock"
              datarank=1                  
              data_dims1(1)=1
              CALL h5screate_simple_f(datarank, data_dims1, &
                   dspace_id, error)                  
              call h5dcreate_f(file_id, dsetname(2), H5T_NATIVE_REAL &
                   ,dspace_id,dset_id,error)                  
              CALL h5dwrite_f(dset_id, H5T_NATIVE_REAL,ishock*1., &
                   data_dims1,error)                  
              call h5dclose_f(dset_id,error)
              call h5sclose_f(dspace_id,error)      
              
              dsetname(2)="xsl"               
              datarank=1               
              data_dims1(1)=max(int(mx0/dxslice),1)
              CALL h5screate_simple_f(datarank,data_dims1,dspace_id, &
                   error)               
              call h5dcreate_f(file_id, dsetname(2), H5T_NATIVE_REAL &
                   ,dspace_id,dset_id,error)               
              CALL h5dwrite_f(dset_id, H5T_NATIVE_REAL,xslice, &
                   data_dims1,error)
              
              call h5dclose_f(dset_id,error)
              call h5sclose_f(dspace_id,error) 
              
              dsetname(2)="dens"
              datarank=1
              data_dims1(1)=max(int(mx0/dxslice),1)
              CALL h5screate_simple_f(datarank, data_dims1, & 
                   dspace_id, error)                  
              call h5dcreate_f(file_id, dsetname(2), H5T_NATIVE_REAL &
                   ,dspace_id,dset_id,error)                  
              CALL h5dwrite_f(dset_id, H5T_NATIVE_REAL,dens1, &
                   data_dims1,error)                  
              call h5dclose_f(dset_id,error)
              call h5sclose_f(dspace_id,error)      
              
              dsetname(2)="dens0" 
              datarank=1                  
              data_dims1(1)=max(int(mx0/10),1)
              CALL h5screate_simple_f(datarank, data_dims1, &
                   dspace_id, error)                  
              call h5dcreate_f(file_id, dsetname(2), H5T_NATIVE_REAL &
                   ,dspace_id,dset_id,error)                  
              CALL h5dwrite_f(dset_id, H5T_NATIVE_REAL,dens, &
                   data_dims1,error)                  
              call h5dclose_f(dset_id,error)
              call h5sclose_f(dspace_id,error)            
           endif !var eq 1
           
           !     Close main file.!
           If (var.eq.3) THEN
              CALL h5fclose_f(file_id, error)
              CALL h5close_f(error)
           ENDIF
#endif
        endif               !if rank eq 0
        if(debug) print *, rank, "done writing, now cleaning"
        
        call MPI_BARRIER(MPI_COMM_WORLD,ierr)
        if(debug) print *, rank, "in spec deallocating"
        if(debug) print *, rank,"al spece?", allocated(spece)
        call MPI_BARRIER(MPI_COMM_WORLD,ierr)
        if(debug) print *, rank,"al specp?", allocated(specp)
        call MPI_BARRIER(MPI_COMM_WORLD,ierr)
        if(debug) print *, rank,"al dens?", allocated(dens)
        
        if(allocated(logmom)) deallocate(logmom)
        call MPI_BARRIER(MPI_COMM_WORLD,ierr)
        if(debug) print *, rank, "removed logmom"
        
        if(allocated(xslice)) deallocate(xslice)
        call MPI_BARRIER(MPI_COMM_WORLD,ierr)
        if(debug) print *, rank, "removed xslice"
        
        if(allocated(spece)) deallocate(spece)
        call MPI_BARRIER(MPI_COMM_WORLD,ierr)
        if(debug) print *, rank, "deallocated spece"
        
        if(allocated(specp)) deallocate(specp)
        call MPI_BARRIER(MPI_COMM_WORLD,ierr)
        if(debug) print *, rank, "removed specp"

        if(allocated(speceb)) deallocate(speceb)
        call MPI_BARRIER(MPI_COMM_WORLD,ierr)
        if(debug) print *, rank, "deallocated speceb"
        
        if(allocated(specpb)) deallocate(specpb)
        call MPI_BARRIER(MPI_COMM_WORLD,ierr)
        if(debug) print *, rank, "removed specpb"
        
        if(allocated(specin))deallocate(specin)
        call MPI_BARRIER(MPI_COMM_WORLD,ierr)
        if(debug) print *, rank, "removed specin"
        
        if(debug) print *, rank, "about to deallocate dens"
        if(allocated(dens)) deallocate(dens)
        call MPI_BARRIER(MPI_COMM_WORLD,ierr)
        if(debug) print *, rank, "deallocated dens"            
        
        if(allocated(dens1)) deallocate(dens1)
        call MPI_BARRIER(MPI_COMM_WORLD,ierr)
        if(debug) print *, rank, "removed dens1"
        
        if(allocated(dens2)) deallocate(dens2)
        call MPI_BARRIER(MPI_COMM_WORLD,ierr)
        if(debug) print *, rank, "removed dens2"
        
        if(allocated(densin)) deallocate(densin)
        call MPI_BARRIER(MPI_COMM_WORLD,ierr)
        if(debug) print *, rank, "removed densin"

        if (xyzspec) then
           if(allocated(tspece)) deallocate(tspece)
           call MPI_BARRIER(MPI_COMM_WORLD,ierr)
           if(debug) print *, rank, "removed tspece"
        
           if(allocated(tspecp)) deallocate(tspecp)
           call MPI_BARRIER(MPI_COMM_WORLD,ierr)
           if(debug) print *, rank, "removed tspecp"
           
           if(allocated(tspeceb)) deallocate(tspeceb)
           call MPI_BARRIER(MPI_COMM_WORLD,ierr)
           if(debug) print *, rank, "removed tspeceb"
           
           if(allocated(tspecpb)) deallocate(tspecpb)
           call MPI_BARRIER(MPI_COMM_WORLD,ierr)
           if(debug) print *, rank, "removed tspecpb"
        endif
        
        if(debug) print *, rank, "done momentum spec"
        
     end do                 ! var
  endif                     ! if lap is right for output
  
  call MPI_BARRIER(MPI_COMM_WORLD,ierr)
  
end subroutine save_momentumspec


!-------------------------------------------------------------------------------
! 						subroutine write_restart					 
!																		
! Saves restart files at user-specified times
!  							
!-------------------------------------------------------------------------------
subroutine write_restart()
  
  implicit none
  
  ! local variables  
  integer ::sendbuf,ierr,n
  logical exst1
  integer ::token,request,status(MPI_STATUS_SIZE),request1 &
       ,status1(MPI_STATUS_SIZE)
  
  !save restart file
  !naming convention: write to rest(flds/prtl).jobname.rank.d
  !if the interval is met, rename the file into rest(flds/prtl).jobname.lapNNNNN.rank.d
	
  sendbuf=0
  if(rank .eq. 0) then 
     inquire(file="write_restart",exist=exst1)
     if(exst1)sendbuf=1
  endif

  call mpi_bcast(sendbuf,1,mpi_integer,0,mpi_comm_world,ierr)
  exst1=.false.

  if(sendbuf .eq. 1) exst1=.true.
  
  if(modulo(lap,1) .eq.0 .or. exst1) then
     
#ifdef LOCALRESTART
     if(intrestlocal .gt. 0 .and. modulo(lap,intrestlocal) .eq.0) then 
        open(unit=7,file=frestartfldloc, form='unformatted')
        rewind 7
        if(rank.eq.0)print *, rank, ": writing local restart file"
        write(7)mx,my,mz,ex,ey,ez,bx,by,bz,dseed,lap,xinject &
             ,xinject2,xinject3,leftwall !,split_E_ions,split_E_lecs 
        close(7)

        open(unit=8,file=frestartprtloc, form='unformatted')
        rewind 8
        if(rank.eq.0)print *, rank,": wrote fields, writing local particles"

        write(8)ions,lecs,maxptl,maxhlf, totalpartnum, &
             (p(n)%x,n=1,ions),(p(n)%x,n=maxhlf+1,maxhlf+lecs), &
             (p(n)%y,n=1,ions),(p(n)%y,n=maxhlf+1,maxhlf+lecs), &
             (p(n)%z,n=1,ions),(p(n)%z,n=maxhlf+1,maxhlf+lecs), &
             (p(n)%u,n=1,ions),(p(n)%u,n=maxhlf+1,maxhlf+lecs), &
             (p(n)%v,n=1,ions),(p(n)%v,n=maxhlf+1,maxhlf+lecs), &
             (p(n)%w,n=1,ions),(p(n)%w,n=maxhlf+1,maxhlf+lecs), &
             (p(n)%ch,n=1,ions),(p(n)%ch,n=maxhlf+1,maxhlf+lecs), &
             (p(n)%ind,n=1,ions),(p(n)%ind,n=maxhlf+1,maxhlf+lecs), &
             (p(n)%proc,n=1,ions),(p(n)%proc,n=maxhlf+1,maxhlf+lecs), &
             (p(n)%splitlev,n=1,ions), &
             (p(n)%splitlev,n=maxhlf+1,maxhlf+lecs)
        
        close(8)
     endif ! if(modulo(lap,intrestlocal) .eq.0)
7645 continue
#endif
     !------------------------------
     
     call MPI_BARRIER(MPI_COMM_WORLD,ierr)
		
     !     check if I have time to save the restart file
     if(modulo(lap,intrestart) .eq. 0  .or.exst1) then
        if(rank.eq.0) then
           time_end=mpi_wtime()
           time_diff=time_end-time_beg !elapsed time
           if(debug)print *,"times for restart",real(timespan,4), &
                real(time_diff,4),real(time_cut,4) 
           if(timespan-time_diff.gt.time_cut)then 
              restart_status=1
           else
              restart_status=0
           endif

	   if(exst1) restart_status=1

        endif
        
        call MPI_BCAST(restart_status,1,MPI_INTEGER,0, &
             MPI_COMM_WORLD,ierr)               
        if(debug)print *,"rank, restart_status",rank,restart_status
     endif
     
     call MPI_BARRIER(MPI_COMM_WORLD,ierr)                
                     
     if((modulo(lap,intrestart) .eq. 0  .or.exst1) .and. restart_status .eq. 1) then 
        
        !install a semaphore
        
#ifdef MPI
        if(rank .gt. 0) then 
           token=rank-1
           !wait for the message from the previous rank
           call mpi_irecv(token,1, &
                mpi_integer,rank-1,100,MPI_Comm_WORLD,request,ierr)
           call mpi_wait(request,status,ierr)
           
           if (debug)print *,"rank ",rank," received the token from ",rank-1
           
           if(modulo(rank,15).ne.0 .and. rank .ne. size0-1) then
              token=rank
              call mpi_isend(token,1, &
                   mpi_integer,rank+1,100,MPI_Comm_WORLD,request1,ierr)
              call mpi_wait(request1,status1,ierr)
           endif
        endif
3004    continue
#endif
        
        open(unit=7,file=frestartfld, form='unformatted')
        close(7,status='delete')
        
        open(unit=7,file=frestartfld, form='unformatted')
        rewind 7
        if(rank.eq.0)print *,"writing restart file"
        write(7)mx,my,mz,ex,ey,ez,bx,by,bz,dseed,lap,xinject &
             ,xinject2,xinject3,leftwall,walloc !,split_E_ions,split_E_lecs 
        close(7)
        
        open(unit=8,file=frestartprt, form='unformatted')
        close(8,status='delete')
        
        open(unit=8,file=frestartprt, form='unformatted')
        rewind 8
        if(rank .eq. 0) print *,"wrote fields, writing particles"

         ! jaehong 9/2013; inludes totalpartnum.
        write(8)ions,lecs,maxptl,maxhlf, totalpartnum, &
             (p(n)%x,n=1,ions),(p(n)%x,n=maxhlf+1,maxhlf+lecs), &
             (p(n)%y,n=1,ions),(p(n)%y,n=maxhlf+1,maxhlf+lecs), &
             (p(n)%z,n=1,ions),(p(n)%z,n=maxhlf+1,maxhlf+lecs), &
             (p(n)%u,n=1,ions),(p(n)%u,n=maxhlf+1,maxhlf+lecs), &
             (p(n)%v,n=1,ions),(p(n)%v,n=maxhlf+1,maxhlf+lecs), &
             (p(n)%w,n=1,ions),(p(n)%w,n=maxhlf+1,maxhlf+lecs), &
             (p(n)%ch,n=1,ions),(p(n)%ch,n=maxhlf+1,maxhlf+lecs), &
             (p(n)%ind,n=1,ions),(p(n)%ind,n=maxhlf+1,maxhlf+lecs), &
             (p(n)%proc,n=1,ions),(p(n)%proc,n=maxhlf+1,maxhlf+lecs), &
             (p(n)%splitlev,n=1,ions), &
             (p(n)%splitlev,n=maxhlf+1,maxhlf+lecs)
        
        close(8)
        !           if(rank.eq.0) 
        if(rank.eq.0)print *,"done writing restart file"
               
        if(modulo(lap,namedrestartint) .eq. 0) then
           if(rank.eq.0)print *, rank, ": renaming restart"
           write(lapchar,"(i6.6)")lap	
           frestartfldlap="restart/restflds."//trim(jobname) &
                //"."//"lap"//trim(lapchar)//"."//trim(rankchar)//".d"
           frestartprtlap="restart/restprtl."//trim(jobname) &
                //"."//"lap"//trim(lapchar)//"."//trim(rankchar)//".d"
           
           
           open(unit=7,file=frestartfldlap, form='unformatted')
           close(7,status='delete')
           
           open(unit=7,file=frestartfldlap, form='unformatted')
           rewind 7
           if(rank.eq.0)print *,"writing laprestart file"
           write(7)mx,my,mz,ex,ey,ez,bx,by,bz,dseed,lap,xinject &
                ,xinject2,xinject3,leftwall,walloc !,split_E_ions,split_E_lecs 
           close(7)
           
           open(unit=8,file=frestartprtlap, form='unformatted')
           close(8,status='delete')
           
           open(unit=8,file=frestartprtlap, form='unformatted')
           rewind 8
           if(rank.eq.0)print *,"wrote fields, writing particles"
           
           write(8)ions,lecs,maxptl,maxhlf, totalpartnum, &
                (p(n)%x,n=1,ions),(p(n)%x,n=maxhlf+1,maxhlf+lecs), &
                (p(n)%y,n=1,ions),(p(n)%y,n=maxhlf+1,maxhlf+lecs), &
                (p(n)%z,n=1,ions),(p(n)%z,n=maxhlf+1,maxhlf+lecs), &
                (p(n)%u,n=1,ions),(p(n)%u,n=maxhlf+1,maxhlf+lecs), &
                (p(n)%v,n=1,ions),(p(n)%v,n=maxhlf+1,maxhlf+lecs), &
                (p(n)%w,n=1,ions),(p(n)%w,n=maxhlf+1,maxhlf+lecs), &
                (p(n)%ch,n=1,ions),(p(n)%ch,n=maxhlf+1,maxhlf+lecs), &
                (p(n)%ind,n=1,ions),(p(n)%ind,n=maxhlf+1,maxhlf+lecs), &
                (p(n)%proc,n=1,ions),(p(n)%proc,n=maxhlf+1,maxhlf+lecs), &
                (p(n)%splitlev,n=1,ions), &
                (p(n)%splitlev,n=maxhlf+1,maxhlf+lecs)
           
           close(8)
           !           if(rank.eq.0) 
           if(rank.eq.0)print *,"done writing restart file"
        endif
        
        !send the message to the next rank, if not the last
        
#ifdef MPI
        if(modulo(rank,15).eq.0 .and. rank .ne. size0-1) then
           token=rank
           call mpi_isend(token,1, &
                mpi_integer,rank+1,100,MPI_Comm_WORLD,request,ierr)
        endif
2121    continue
        
        call mpi_barrier(MPI_COMM_WORLD,ierr)
#endif

        if(rank .eq. 0 .and. exst1) then 
           open(unit=11,file="write_restart",form='unformatted')
           close(11,status='delete')
        endif
        
        !after a barrier can delete the temporary files in /tmp
        !            open(unit=7,file=frestartfldloc, form='unformatted')
        !            close(7,status='delete')
        !            open(unit=7,file=frestartprtloc, form='unformatted')
        !            close(7,status='delete')
        
     endif
  endif          
  
end subroutine write_restart

!-------------------------------------------------------------------------------
! 						subroutine save_param					 
!																		
! Saves parameters of the simulation run
!  							
!-------------------------------------------------------------------------------
subroutine save_param()

  implicit none

  character(len=20) fparam
  integer ::ind, datarank, error, i
  integer,dimension(30)::intdata
  real,dimension(30)::realdata
  character(len=10) dsetnamei(30),dsetnamer(30)
#ifdef HDF5
  integer(HSIZE_T) dimsf(3)
  integer(HID_T) :: file_id ! File identifier 
  integer(HID_T) :: dset_id(40) ! Dataset identifier 
  integer(HID_T) :: filespace(40) ! Dataspace identifier in file 

#endif
!  integer(HID_T) :: memspace ! Dataspace identifier in memory
!  integer(HID_T) :: plist_id ! Property list identifier
  
  if(lap .ge.pltstart .and. modulo((lap-pltstart),interval) .eq.0) then

     ind=(lap-pltstart)/interval
     write(indchar,"(i3.3)")ind
     !save parameter file
     !          goto 608
     fparam="output/param."//trim(indchar)
     
     if(rank .eq. 0) then
        print *,"",fparam
        open(unit=12,file=fparam,form='unformatted')
        close(12,status='delete')
        call h5open_f (error)
        call h5fcreate_f(fparam, H5F_ACC_TRUNC_F, file_id, error)
        datarank=1
        dimsf(1)=1
        
        !first pack them into an array
        intdata(1)=int(mx0)
        intdata(2)=my0
        intdata(3)=mz0
        intdata(4)=1 !caseinit
        intdata(5)=interval
        intdata(6)=torqint
        intdata(7)=pltstart
        intdata(8)=istep
        intdata(9)=istep1
        intdata(10)=stride
        intdata(11)=ntimes
        intdata(12)=0
!        if(cooling) intdata(12)=1
        intdata(12)=1
        intdata(13)=sizex
        intdata(14)=sizey
		intdata(15)=dlaplec
		intdata(16)=dlapion
		intdata(17)=teststartlec
		intdata(18)=teststartion
		intdata(19)=testendlec
		intdata(20)=testendion
			
        realdata(1)=c
        realdata(2)=sigma
        realdata(3)=c_omp
        realdata(4)=gamma0
        realdata(5)=delgam
        realdata(6)=ppc0
        realdata(7)=me
        realdata(8)=mi
        realdata(9)=dummyvar
        realdata(10)= 0 !acool
        realdata(11)=lap*(c/c_omp)
        realdata(12)=qi
        realdata(13)=btheta
        realdata(14)=bphi
        realdata(15)=xinject2
        realdata(16)=walloc
        
        dsetnamei(1)="mx0"
        dsetnamei(2)="my0"
        dsetnamei(3)="mz0"
        dsetnamei(4)="caseinit"
        dsetnamei(5)="interval"
        dsetnamei(6)="torqint"
        dsetnamei(7)="pltstart"
        dsetnamei(8)="istep"
        dsetnamei(9)="istep1"
        dsetnamei(10)="stride"
        dsetnamei(11)="ntimes"
        dsetnamei(12)="cooling"
        dsetnamei(13)="sizex"
        dsetnamei(14)="sizey"
        dsetnamei(15)="dlaplec"
        dsetnamei(16)="dlapion"
        dsetnamei(17)="teststartlec"
		dsetnamei(18)="teststartion"
		dsetnamei(19)="testendlec"
		dsetnamei(20)="testendion"

        
        dsetnamer(1)="c"
        dsetnamer(2)="sigma"
        dsetnamer(3)="c_omp"
        dsetnamer(4)="gamma0"
        dsetnamer(5)="delgam"
        dsetnamer(6)="ppc0"
        dsetnamer(7)="me"
        dsetnamer(8)="mi"
        dsetnamer(9)="dummy"
        dsetnamer(10)="acool"
        dsetnamer(11)="time"
        dsetnamer(12)="qi"
        dsetnamer(13)="btheta"
        dsetnamer(14)="bphi"
        dsetnamer(15)="xinject2"
        dsetnamer(16)="walloc"
        
         do i=1,20
           if(debug)print *, rank,": dataseti ",i
           call h5screate_simple_f(datarank, dimsf, filespace(1), &
                error)
           call h5dcreate_f(file_id, dsetnamei(i), H5T_NATIVE_INTEGER, &
                filespace(1),dset_id(1), error)
           call h5dwrite_integer_1(dset_id(1),H5T_NATIVE_INTEGER &
                ,intdata(i),dimsf,error)
           call h5dclose_f(dset_id(1), error)
           call h5sclose_f(filespace(1), error)
        enddo
        
        do i=1,16
           if(debug)print *,rank, ": datasetr ",i
           call h5screate_simple_f(datarank, dimsf, filespace(1), &
                error)
           call h5dcreate_f(file_id, dsetnamer(i), H5T_NATIVE_REAL, &
                filespace(1),dset_id(1), error)
           call h5dwrite_real_1(dset_id(1),H5T_NATIVE_REAL,realdata(i) &
                ,dimsf,error)
           call h5dclose_f(dset_id(1), error)
           call h5sclose_f(filespace(1), error)
        enddo   
        
        datarank=1
        dimsf=sizex
        dsetnamei(22)="mx"
        
 	    call h5screate_simple_f(datarank, dimsf, filespace(1), &
                error)
        call h5dcreate_f(file_id, dsetnamei(22), H5T_NATIVE_INTEGER, &
                filespace(1),dset_id(1), error)
        call h5dwrite_integer_1(dset_id(1),H5T_NATIVE_INTEGER &
                ,mxl(1:sizex),dimsf,error)
        call h5dclose_f(dset_id(1), error)
        call h5sclose_f(filespace(1), error)
              
        
        datarank=1
        dimsf=sizey
        dsetnamei(23)="my"
        
 	    call h5screate_simple_f(datarank, dimsf, filespace(1), &
                error)
        call h5dcreate_f(file_id, dsetnamei(23), H5T_NATIVE_INTEGER, &
                filespace(1),dset_id(1), error)
        call h5dwrite_integer_1(dset_id(1),H5T_NATIVE_INTEGER &
                ,myl(1:(sizey-1)*sizex+1:sizex),dimsf,error)
        call h5dclose_f(dset_id(1), error)
        call h5sclose_f(filespace(1), error)
           
           
        call h5fclose_f(file_id,error)
        call h5close_f(error)
        
     endif
608  continue
     
  endif !if (lap .ge. pltstart)
  
end subroutine save_param


!-------------------------------------------------------------------------------
! 						subroutine output_tot					 
!																		
! Saves fluid (flds.tot.) and particle (prtl.tot.) data 
!  							
!-------------------------------------------------------------------------------

subroutine output_tot()
  
  implicit none
  
  ! local variables
  integer ::ind,kglob,kloc,kstart,kfinish,nz,info,n,i,j,k,ierr
  integer ::istart,ifinish,jstart,jfinish,nx,ny,iv,mz1,my1,mx1
  real bx0,by0,bz0,ex0,ey0,ez0,bxd,byd,bzd,exd,eyd,ezd,curx0, &
       cury0,curz0
  integer ::error ! Error flags
  integer ::nvars, midvars
  integer ::all_ions(size0), all_lecs(size0)
  integer*8 all_ions_lng(size0), all_lecs_lng(size0)
  character(len=6) :: dsetname(40),varname  
  integer ::datarank
  integer :: ifindi,ifind
#ifdef HDF5
  integer(HID_T) :: file_id ! File identifier 
  integer(HID_T) :: dset_id(40) ! Dataset identifier 
  integer(HID_T) :: filespace(40) ! Dataspace identifier in file 
  integer(HID_T) :: memspace ! Dataspace identifier in memory
  integer(HID_T) :: plist_id ! Property list identifier
  integer(HSIZE_T) dimsf(3),dimsfions(1), dimsflecs(1)
  integer(HSIZE_T) dimsfi(7),dimsfiions(7), dimsfilecs(7)
  integer(HSIZE_T), DIMENSION(3) :: count  
  integer(HSSIZE_T), DIMENSION(3) :: offset
  integer(HSIZE_T), DIMENSION(1) :: countpart  
  integer(HSSIZE_T), DIMENSION(1) :: offsetpart
#endif 
  integer ::FILE_INFO_TEMPLATE !for MPI file info
  integer ::procn
  integer ::mx1list(size0),my1list(size0),mz1list(size0) 
  real*4, allocatable :: temporary0(:,:,:),temporary0_vec(:)  
  real*4,allocatable:: temporary_vec(:)
  integer*8, allocatable:: selecti_vec(:),selecte_vec(:) 
  integer :: status(MPI_STATUS_SIZE)
  integer ::tmp1,tmp1i,tmp1e
  logical :: saverank
  integer ::token,request,status2(MPI_STATUS_SIZE),request1 &
      ,status1(MPI_STATUS_SIZE)
  logical :: writecurr
 

  if(lap .ge.pltstart .and. modulo((lap-pltstart),interval) .eq.0)then

     !timestamp of the saved output
     ind=(lap-pltstart)/interval	    
     write(fnamefld,"(a16,i3.3)") "output/flds.tot.",ind
     write(fnameprt,"(a16,i3.3)") "output/prtl.tot.",ind     
     write(indchar,"(i3.3)")ind     
     if(rank.eq.0.and.debug)write(*,*) rank,": ",fnamefld,'  ',fnameprt
     
     !just in case destroy the file
     if(rank .eq. 0) then 
        open(unit=11,file=fnamefld,form='unformatted')
        close(11,status='delete')
        open(unit=12,file=fnameprt,form='unformatted')
        close(12,status='delete')
     endif
     
     call MPI_BARRIER(MPI_COMM_WORLD,ierr)     

     
     !specify FILE_INFO_TEMPLATE for parallel output
#ifndef serIO 
#ifdef MPI     
     call MPI_Info_create(FILE_INFO_TEMPLATE,error)
     call MPI_Info_set(FILE_INFO_TEMPLATE, "access_style",  &
          "write_once",error)
     call MPI_Info_set(FILE_INFO_TEMPLATE,  &
          "collective_buffering", "true",error)
     call MPI_Info_set(FILE_INFO_TEMPLATE, "cb_block_size", &
          "4194304",error)
     call MPI_Info_set(FILE_INFO_TEMPLATE, "cb_buffer_size", &
          "16777216",error)
     call MPI_Info_set(FILE_INFO_TEMPLATE, "cb_nodes", &
          "1",error)
#endif
#endif


     !FLDS.TOT. saving
     if (rank.eq.0)print *,"",fnamefld
     
     !user input: array of variables to be saved for flds.tot.
     !usage:
     !-the first entry should have the largest number of characters
     !(if needed, supplemented with empty spaces)
     !to be chosen from:
     !dens,densi: density and density of ions
     !ex,ey,ez,bx,by,bz,jx,jy,jz:fields and currents
     !v3x,v3xi,v3y,v3yi,v3z,v3zi:3-velocity
     !v4x,v4xi,v4y,v4yi,v4z,v4zi:4-velocity (momentum)
     !ken,keni:particle total energy
#ifdef twoD
  if (sigma .eq. 0. .and. pcosthmult .eq. 0) then 
     nvars=15
     dsetname(1:nvars)=(/&
          'dens  ','densi '&
!          ,'bdens ','bdensi'&
          ,'v3x   ','v3xi  ','v3y   ','v3yi  '&
!          ,'v4x   ','v4xi  ','v4y   ','v4yi  ','ken   ','keni  '&
          ,'ex    ','ey    ','ez    ','bx    ','by    ','bz    ','jx    ','jy    ','jz    '&
       /)
  else
  nvars=17!27!19
  dsetname(1:nvars)=(/&
       'dens  ','densi '&
!       ,'bdens ','bdensi'&
       ,'v3x   ','v3xi  ','v3y   ','v3yi  ','v3z   ','v3zi  '&
!       ,'v4x   ','v4xi  ','v4y   ','v4yi  ','v4z   ','v4zi  ','ken   ','keni  '&
       ,'ex    ','ey    ','ez    ','bx    ','by    ','bz    ','jx    ','jy    ','jz    '&
       /)
  endif
#else
  nvars=17
  dsetname(1:nvars)=(/&
       'dens  ','densi '&
!       ,'bdens ','bdensi'&
       ,'v3x   ','v3xi  ','v3y   ','v3yi  ','v3z   ','v3zi  '&
!       ,'v4x   ','v4xi  ','v4y   ','v4yi  ','v4z   ','v4zi  ','ken   ','keni  '&
       ,'ex    ','ey    ','ez    ','bx    ','by    ','bz    ','jx    ','jy    ','jz    '&
       /)
#endif

     
     !compute global size of fluid arrays to be saved
     dimsf(1)=mx0/istep
     dimsf(2)=my0/istep  
#ifndef twoD
     dimsf(3)=mz0/istep
#else
     dimsf(3)=1
#endif 
     dimsfi(1:3)=dimsf(1:3)
     dimsfi(4:7)=0
     datarank=3
     
#ifndef twoD     
     !computing my1 and mz1 for each processor
     kstart=3
     kfinish=mz-3
     if(rank/(sizex*sizey) .eq. 0) kstart=1
     if(rank/(sizex*sizey) .eq. sizez-1) kfinish=mz
110  if(modulo(kstart+mzcum,istep) .eq. 0) goto 120
     kstart=kstart+1
     goto 110
120  continue
     
130  if(modulo(kfinish+mzcum,istep) .eq. 0) goto 140
     kfinish=kfinish-1
     goto 130
140  continue
     
     mz1=(kfinish-kstart)/istep+1
#else
     kstart=1
     kfinish=1
     mz1=1
#endif
     
     jstart=3
     jfinish=my-3
     if(modulo(rank,sizex*sizey)/sizex .eq. 0) jstart=1
     if(modulo(rank,sizex*sizey)/sizex .eq. sizey-1) jfinish=my

210  if(modulo(jstart+mycum,istep) .eq. 0) goto 220	
     jstart=jstart+1
     goto 210
220  continue

230  if(modulo(jfinish+mycum,istep) .eq. 0)  goto 240
     jfinish=jfinish-1
     goto 230
240  continue      
     
     my1=(jfinish-jstart)/istep+1
     
     istart=3
     ifinish=mx-3
     if(modulo(rank,sizex) .eq. 0) istart=1
     if(modulo(rank,sizex) .eq. sizex-1) ifinish=mx
     
310  if(modulo(istart+int(mxcum),istep) .eq. 0) goto 320	
     istart=istart+1
     goto 310
320  continue

330  if(modulo(ifinish+int(mxcum),istep) .eq. 0)  goto 340
     ifinish=ifinish-1
     goto 330
340  continue
     
     mx1=(ifinish-istart)/istep+1
     
  
! determine whether the current processor should save any cell
     if (kfinish .lt. kstart) then
        kfinish=1
        kstart=1
        mz1=0
     endif
     if (jfinish .lt. jstart) then
        jfinish=1
        jstart=1
        my1=0
     endif
     
     if (ifinish .lt. istart) then
        ifinish=1
        istart=1
        mx1=0
     endif     
     
     !allocate temporary1 on every processor
     if(debug) print *, rank,"in output b4 allocate"
     allocate(temporary1(max(mx1,1),max(my1,1),max(mz1,1)))


     !getting information about the size of temporary1
#ifdef MPI
     if(debug) print *, rank,": in output, b4 allgather", mx1
     call mpi_allgather(mx1,1,mpi_integer,mx1list,1,mpi_integer &
          ,mpi_comm_world, error)
     if(debug) print *, rank,": in output, b4 allgather", my1
     call mpi_allgather(my1,1,mpi_integer,my1list,1,mpi_integer &
          ,mpi_comm_world, error)
     if(debug) print *, rank,": in output, b4 allgather", mz1
     call mpi_allgather(mz1,1,mpi_integer,mz1list,1,mpi_integer &
          ,mpi_comm_world, error) 
#ifdef serIO       
     if (rank .eq. 0) then
        allocate(temporary0(maxval(mx1list),maxval(my1list), &
             maxval(mz1list)))
        temporary0(:,:,:)=0.
     endif
#endif
#endif

!decide if I need to write currents to disk
     writecurr=.false.
     do iv=1,nvars
        varname=trim(dsetname(iv))
        if (varname .eq. 'jx' .or. varname .eq. 'jy' .or. varname .eq. 'jz' ) writecurr=.true. 
     enddo
     
     if (writecurr) then

        !install a semaphore        
#ifdef MPI
        if(rank .gt. 0) then 
           token=rank-1
           !wait for the message from the previous rank
           call mpi_irecv(token,1, &
                mpi_integer,rank-1,100,MPI_Comm_WORLD,request,ierr)
           call mpi_wait(request,status2,ierr)
           
           if (debug)print *,"rank ",rank," received the token from ",rank-1
           
           if(modulo(rank,15).ne.0 .and. rank .ne. size0-1) then
              token=rank
              call mpi_isend(token,1, &
                   mpi_integer,rank+1,100,MPI_Comm_WORLD,request1,ierr)
              call mpi_wait(request1,status1,ierr)
           endif
        endif
3004    continue
#endif

     !writing currents (all procs)
     open(unit=7,file=fenlargeloc, form='unformatted')
     rewind 7
     write(7) (((curx(i,j,k),i=1,mx),j=1,my),k=1,mz), &
          (((cury(i,j,k),i=1,mx),j=1,my),k=1,mz), &
          (((curz(i,j,k),i=1,mx),j=1,my),k=1,mz)
     close(7)

        !send the message to the next rank, if not the last        
#ifdef MPI
        if(modulo(rank,15).eq.0 .and. rank .ne. size0-1) then
           token=rank
           call mpi_isend(token,1, &
                mpi_integer,rank+1,100,MPI_Comm_WORLD,request,ierr)
        endif
2121    continue
        
        call mpi_barrier(MPI_COMM_WORLD,ierr)
#endif

     endif !writecurr

     !
     !  Initialize FORTRAN predefined datatypes
     !
#ifdef serIO
     if (rank .eq. 0) then
        call h5open_f(error)
     endif
#else
     call h5open_f(error) 
#endif
   
     ! 
     ! Setup file access property list if parallel I/O access.
     !   
#ifndef serIO
#ifdef MPI
     call h5pcreate_f(H5P_FILE_ACCESS_F, plist_id, error)
     !      call H5Pset_sieve_buf_size_f(plist_id, 262144,error) 
     !      call H5Pset_alignment_f(plist_id, 524288, 262144,error)
     call h5pset_fapl_mpio_f(plist_id, mpi_comm_world, &
          FILE_INFO_TEMPLATE,error)
#endif
#endif
     
     !
     ! Create the file
     !
#ifdef MPI
#ifdef serIO
     if (rank .eq. 0) then
        call h5fcreate_f(fnamefld, H5F_ACC_TRUNC_F, file_id, error, &
             h5p_default_f,h5p_default_f)
     endif
#else
     call h5fcreate_f(fnamefld, H5F_ACC_TRUNC_F, file_id, error, &
          access_prp = plist_id)
     call h5pclose_f(plist_id, error)
#endif
#else
     call h5fcreate_f(fnamefld, H5F_ACC_TRUNC_F, file_id, error)
#endif  

     !
     ! Create the dataspace for the dataset
     !
#ifdef serIO
     if (rank .eq. 0) then
        do i=1,nvars
           call h5screate_simple_f(datarank, dimsf, filespace(i), error)
        enddo
     endif
#else
     do i=1,nvars
        call h5screate_simple_f(datarank, dimsf, filespace(i), error)
     enddo
#endif
     
     !
     ! Loop over variables
     !
     do iv=1,nvars

        varname=trim(dsetname(iv))

        !
        ! Define the array
        !        
        !pre-processing of some field quantities
        !total density and ion density
        if(varname.eq.'dens' ) call meanq_fld_cur('tdens')
        if(varname.eq.'densi') call meanq_fld_cur('idens')
        if(varname.eq.'bdens' ) call meanq_fld_cur('btden')
        if(varname.eq.'bdensi') call meanq_fld_cur('biden')
        !velocities, momenta and energy
        if(varname.eq.'v3x'  ) call meanq_fld_cur('tbetx')
        if(varname.eq.'v3y'  ) call meanq_fld_cur('tbety')
        if(varname.eq.'v3z'  ) call meanq_fld_cur('tbetz') 
        if(varname.eq.'v3xi' ) call meanq_fld_cur('ibetx')
        if(varname.eq.'v3yi' ) call meanq_fld_cur('ibety')
        if(varname.eq.'v3zi' ) call meanq_fld_cur('ibetz')
        if(varname.eq.'v4x'  ) call meanq_fld_cur('tmomx')
        if(varname.eq.'v4y'  ) call meanq_fld_cur('tmomy')
        if(varname.eq.'v4z'  ) call meanq_fld_cur('tmomz')
        if(varname.eq.'v4xi' ) call meanq_fld_cur('imomx')
        if(varname.eq.'v4yi' ) call meanq_fld_cur('imomy')
        if(varname.eq.'v4zi' ) call meanq_fld_cur('imomz')
        if(varname.eq.'ken'  ) call meanq_fld_cur('tener')
        if(varname.eq.'keni' ) call meanq_fld_cur('iener')     
        !electric currents
        if(varname.eq.'jx'.or.varname.eq.'jy'.or.varname.eq.'jz') then 
!        if(varname.eq.'jx'.or.varname.eq.'jy') then 
           !     restore currents from file
           open(unit=7,file=fenlargeloc, form='unformatted')
           rewind 7
           if(debug.and.rank.eq.0) print *, rank,":reading back current file"
           read(7)(((curx(i,j,k),i=1,mx),j=1,my),k=1,mz), &
                (((cury(i,j,k),i=1,mx),j=1,my),k=1,mz), &
                (((curz(i,j,k),i=1,mx),j=1,my),k=1,mz)
           if(debug .and. rank.eq.0)print *, rank,"resetting arrays"           
        endif
        
        ! define temporary1 for every processor

        do k=kstart, kfinish, istep 
           nz=(k-kstart)/istep+1
           do j=jstart, jfinish,istep 
              ny=(j-jstart)/istep+1
              do i=istart, ifinish,istep 
              	 nx=(i-istart)/istep+1
  
              	 select case (varname)
                 case('dens','densi')
                    temporary1(nx,ny,nz)=real(curx(i,j,k),4)
                 case('ex')
                    call interpolfld(real(i),real(j), &
                         real(k),bx0,by0,bz0,ex0,ey0,ez0,bxd, &
                         byd,bzd,exd,eyd,ezd) 
                    temporary1(nx,ny,nz)=real(ex0+exd,4) 
                 case('ey')
                    call interpolfld(real(i),real(j), &
                         real(k),bx0,by0,bz0,ex0,ey0,ez0,bxd, &
                         byd,bzd,exd,eyd,ezd) 
                    temporary1(nx,ny,nz)=real(ey0+eyd,4) 
                 case('bz')
                    call interpolfld(real(i),real(j), &
                         real(k),bx0,by0,bz0,ex0,ey0,ez0,bxd, &
                         byd,bzd,exd,eyd,ezd) 
                    temporary1(nx,ny,nz)=real(bz0+bzd,4)
                 case('ez')
                    call interpolfld(real(i),real(j), &
                         real(k),bx0,by0,bz0,ex0,ey0,ez0,bxd, &
                         byd,bzd,exd,eyd,ezd) 
                    temporary1(nx,ny,nz)=real(ez0+ezd,4) 
                 case('bx')
                    call interpolfld(real(i),real(j), &
                         real(k),bx0,by0,bz0,ex0,ey0,ez0,bxd, &
                         byd,bzd,exd,eyd,ezd) 
                    temporary1(nx,ny,nz)=real(bx0+bxd,4)
                 case('by')
                    call interpolfld(real(i),real(j), &
                         real(k),bx0,by0,bz0,ex0,ey0,ez0,bxd, &
                         byd,bzd,exd,eyd,ezd) 
                    temporary1(nx,ny,nz)=real(by0+byd,4)
                 case('jx')
                    !                   call interpolcurrent(real(i),real(j), &
                    !                        real(k),curx0, cury0, curz0)
                    curx0=curx(i,j,k)
                    temporary1(nx,ny,nz)=real(curx0,4)
                 case('jy')
                    !                   call interpolcurrent(real(i),real(j), &
                    !                        real(k),curx0, cury0, curz0)
                    cury0=cury(i,j,k)
                    temporary1(nx,ny,nz)=real(cury0,4)
                 case('jz')
                    !                   call interpolcurrent(real(i),real(j), &
                    !                        real(k),curx0, cury0, curz0)
                    curz0=curz(i,j,k)
                    temporary1(nx,ny,nz)=real(curz0,4)
                 case('ken','v3x','v3y','v3z','v4x','v4y','v4z',&
                      'keni','v3xi','v3yi','v3zi','v4xi','v4yi','v4zi')
                    temporary1(nx,ny,nz)=real(curx(i,j,k),4)
                 end select
 
              enddo
           enddo
        enddo
        

#ifdef serIO
        !
        ! Serial writing
        !
        if (rank.eq.0.and.debug) print *,"serial writing of .tot. files"

        do procn=0,size0-1

           if (rank .eq. 0) then
              if (procn .eq. 0) then
                 call h5dcreate_f(file_id,varname,H5T_NATIVE_REAL, &
                      filespace(iv),dset_id(iv), error)
                 call h5dclose_f(dset_id(iv), error)
              endif
           endif
           
           ! Determine which processors should write to file
           saverank=.true.
           if (mx1list(procn+1)*my1list(procn+1)*mz1list(procn+1) .eq. 0) then
              saverank=.false.
           endif

           ! write only if saverank eq. .true.
           if (saverank) then
              
           if (rank .eq. 0) then
              
              if (procn .ne. 0) then !receive from procn                                
#ifdef MPI
                 call MPI_Recv(temporary0(1:mx1list(procn+1),1:my1list(procn+1), &
                      1:mz1list(procn+1)),mx1list(procn+1)*my1list(procn+1) &
                      *mz1list(procn+1),mpi_real,procn,1, &
                      MPI_COMM_WORLD,status,error)
#endif
              else           !procn eq 0
                 temporary0(1:mx1,1:my1,1:mz1)=temporary1               
              endif
              
              !     computing the offset
              count(1) = mx1list(procn+1)
              count(2) = my1list(procn+1)
              count(3) = mz1list(procn+1)
              
              offset(1) = 0
              offset(2) = 0
              offset(3) = 0
#ifndef twoD
              !     z offset
              if(procn/(sizex*sizey) .gt. 0) then 
                 offset(3)=sum(mz1list(1:procn/(sizex*sizey)*(sizex*sizey):(sizex*sizey)))
              endif
#else
              offset(3)=0
#endif
              !     y offset
              if(modulo(procn,sizex*sizey)/sizex .gt. 0) then 
                 offset(2)=sum(my1list(1:modulo(procn,sizex*sizey)/sizex*sizex:sizex))
              endif
              
              !     x offset
              if(modulo(procn,sizex) .gt. 0) then 
                 offset(1)=sum(mx1list((procn/sizex)*sizex+1:procn))
              endif   
              
              !hyperslab selection
#ifdef MPI
              ! this is necessay only if I close the filespace above      
              !               call h5dget_space_f(dset_id(iv), filespace(iv), error)               
              call h5sselect_hyperslab_f(filespace(iv),H5S_SELECT_SET_F &
                   ,offset,count,error)
#endif               
              
              !memory space
#ifdef MPI
              call h5screate_simple_f(datarank, count, memspace, error) 
#endif              
              
              !create or open the dataset
              call h5dopen_f(file_id,varname,dset_id(iv),error)
              

              !writing
              if(debug)print *,rank,": before h5dwrite","filespace", &
                   filespace(iv),"memspace",memspace
              
#ifdef MPI
              call h5dwrite_f(dset_id(iv), H5T_NATIVE_REAL,  &
                   temporary0(1:mx1list(procn+1),1:my1list(procn+1), &
                   1:mz1list(procn+1)),dimsfi,error,mem_space_id= &
                   memspace,file_space_id = filespace(iv))
            !               call h5dwrite_real_1(dset_id(iv), H5T_NATIVE_REAL, 
              !     &                   temporary0(:,:,:),dimsfi,error,file_space_id = 
              !     &                   filespace(iv), mem_space_id=memspace)
#else
              call h5dwrite_f(dset_id(iv),H5T_NATIVE_REAL, &
                   temporary(:,:,:) ,dimsfi,error)
#endif 
              
              !closing dataset and memory space
              if(debug)print *, rank, ": ", iv, "closing d"
              call h5dclose_f(dset_id(iv), error)
              if(debug)print *, rank, ": ", "closing m"
              call h5sclose_f(memspace, error)

           else !rank ne 0
              
              if (procn .eq. rank) then
#ifdef MPI
                 call MPI_Send(temporary1,mx1*my1 &
                      *mz1,mpi_real,0,1, &
                      MPI_COMM_WORLD,error) 
#endif             
              endif
              
           endif !rank ne 0
              
           endif !saverank .eq. .true.
           
        enddo!procn
        
#else 
! finished serial writing

        !
        ! Parallel writing
        !     
        if (rank.eq.0.and.debug) print *,"parallel writing of .tot. files"
           
        ! Determine which processors should write to file
        saverank=.true.
        if (my1*mz1 .eq. 0) then
           saverank=.false.
        endif

        !
        ! Create the dataset with default properties
        !
        call h5dcreate_f(file_id, varname, H5T_NATIVE_REAL, &
             filespace(iv),dset_id(iv), error)
        call h5sclose_f(filespace(iv), error)
        call h5dget_space_f(dset_id(iv), filespace(iv), error)

        !
        ! Create property list for collective dataset write
        !
#ifdef MPI
        call h5pcreate_f(H5P_DATASET_XFER_F, plist_id, error) 
        call h5pset_dxpl_mpio_f(plist_id, H5FD_MPIO_COLLECTIVE_F, error)
#endif
           
        !
        ! Each process defines dataset in memory and writes it to the hyperslab
        ! in the file. 
        !
        count(1) = mx1
        count(2) = my1
        count(3) = mz1
  
        offset(1) = 0
        offset(2) = 0
        offset(3) = 0
#ifndef twoD
        !z offset
       if(rank/(sizex*sizey) .gt. 0) then 
             offset(3)=sum(mz1list(1:rank/(sizex*sizey)*(sizex*sizey):(sizex*sizey)))
       endif
#else
        offset(3)=0
#endif
       !     y offset
       if(modulo(rank,sizex*sizey)/sizex .gt. 0) then 
             offset(2)=sum(my1list(1:modulo(rank,sizex*sizey)/sizex*sizex:sizex))
       endif
              
       !     x offset
       if(modulo(rank,sizex) .gt. 0) then 
            offset(1)=sum(mx1list((rank/sizex)*sizex+1:rank))
       endif 
        
        ! write only if saverank .eq. .true.
        if (saverank) then

#ifdef MPI
        call h5screate_simple_f(datarank, count, memspace, error) 
#endif
        
        ! 
        ! Select hyperslab in the file.
        !
#ifdef MPI              
        call h5sselect_hyperslab_f(filespace(iv),H5S_SELECT_SET_F,offset &
             ,count, error)
#endif
        
        if(debug)print *,rank,": before h5dwrite","filespace",filespace(iv &
             ),memspace
        
#ifdef MPI
        call h5dwrite_real_1(dset_id(iv), H5T_NATIVE_REAL, temporary1(:,: &
             ,:),dimsfi,error,file_space_id = filespace(iv), mem_space_id &
             =memspace,xfer_prp = h5p_default_f)
#else
        call h5dwrite_f(dset_id(iv),H5T_NATIVE_REAL,TEMPORARY1(:,: &
             ,:) ,dimsfi,error)
#endif      
        endif ! saverank .eq. .true.
     
#endif 
 !finished parallel writing    
        

        !
        !Closing dataset and dataspace    
        !
!#ifdef serIO    
!        if (rank .eq. 0) then
!           if(debug)print *, rank, ": ", iv, "closing d"
!           call h5dclose_f(dset_id(iv), error)
!#ifdef MPI
!           if(debug)print *, rank, ": ", iv, "closing s" 
!           call h5sclose_f(filespace(iv), error)
!#endif
!        endif
!#else
!        if(debug)print *, rank, ": ", iv, "closing d"
!        call h5dclose_f(dset_id(iv), error)
!#ifdef MPI
!        if(debug)print *, rank, ": ", iv, "closing s" 
!        call h5sclose_f(filespace(iv), error)     
!#endif
!#endif

!     enddo!nvars


        !
        !Closing dataset and dataspace    
        !
#ifdef serIO    
        if (rank .eq. 0) then
#ifdef MPI
           if(debug)print *, rank, ": ", iv, "closing s" 
           call h5sclose_f(filespace(iv), error)
#endif
        endif
#else
        if(debug)print *, rank, ": ", iv, "closing d"
        call h5dclose_f(dset_id(iv), error)
#ifdef MPI
        if(debug)print *, rank, ": ", iv, "closing s" 
        call h5sclose_f(filespace(iv), error)     
#endif
#endif

     enddo!nvars


     
     
     !
     ! Closing memory space, file and hdf5
     !
!#ifdef serIO
!     if (rank .eq. 0) then
!#ifdef MPI
!        if (saverank) then
!           if(debug)print *, rank, ": ", "closing m"
!           call h5sclose_f(memspace, error)
!        endif ! saverank .eq. .true.
!#endif
!        if(debug)print *, rank, ": ", "closing f"
!        call h5fclose_f(file_id, error)
!        if(debug)print *, rank, ": ", "closing h5"
!        call h5close_f(error)
!        deallocate(temporary0)
!     endif
!#else
!#ifdef MPI
!     if (saverank) then
!        if(debug)print *, rank, ": ", "closing m"
!        call h5sclose_f(memspace, error)
!     endif !saverank .eq. .true.
!     if(debug)print *, rank, ": ", "closing p"
!     call h5pclose_f(plist_id, error)
!#endif
!     if(debug)print *, rank, ": ", "closing f"
!     call h5fclose_f(file_id, error)
!     if(debug)print *, rank, ": ", "closing h5"
!     call h5close_f(error)    
!#endif
     
!     if(debug) print *, rank,": finished writing .tot. fields"
     
!     deallocate(temporary1)



     !
     ! Closing memory space, file and hdf5
     !
#ifdef serIO
     if (rank .eq. 0) then
        if(debug)print *, rank, ": ", "closing f"
        call h5fclose_f(file_id, error)
        if(debug)print *, rank, ": ", "closing h5"
        call h5close_f(error)
        deallocate(temporary0)
     endif
#else
#ifdef MPI
     if (saverank) then
        if(debug)print *, rank, ": ", "closing m"
        call h5sclose_f(memspace, error)
     endif !saverank .eq. .true.
     if(debug)print *, rank, ": ", "closing p"
     call h5pclose_f(plist_id, error)
#endif
     if(debug)print *, rank, ": ", "closing f"
     call h5fclose_f(file_id, error)
     if(debug)print *, rank, ": ", "closing h5"
     call h5close_f(error)    
#endif
     
     if(debug) print *, rank,": finished writing .tot. fields"
     
     deallocate(temporary1)



!---------------------------------------------------------------
!     goto 127
     !--------------------------------------
     ! PRTL.TOT. saving: now write particles
     !--------------------------------------
     if (rank.eq.0)print *,"",fnameprt

     !user input: number of variables to be saved for prtl.tot.
     nvars=21
     midvars=nvars/2

     !user input: array of variables to be saved for prtl.tot.
     !usage:
     !-the first entry should have the largest number of characters
     !(if needed, supplemented with empty spaces)
     !-same number of entries for electrons and ions
     !to be chosen from:
     !xi,yi,zi,xe,ye,ze:position
     !ui,vi,wi,ue,ve,we:4-velocity (momentum)
     !chi,che:weight
     !gammai,gammae:energy
     !indi,inde:index
     !proci,proce:processor
     dsetname(1:nvars)=(/&
          'gammai','xi    ','yi    ','zi    '&
          ,'ui    ','vi    ','wi    '&
          ,'chi   ','indi  ','proci '&
          ,'gammae','xe    ','ye    ','ze    '&
          ,'ue    ','ve    ','we    '&
          ,'che   ','inde  ','proce ','time  '&
          /)
     
     !find out number of ions and lecs
      datarank=1  
#ifdef MPI
     tmp1=max(ions,1)
     call mpi_allgather(tmp1,1,mpi_integer,all_ions,1 &
          ,mpi_integer,mpi_comm_world, error)
     tmp1=max(lecs,1)
     call mpi_allgather(tmp1,1,mpi_integer,all_lecs,1, &
          mpi_integer,mpi_comm_world, error)
     if(rank.eq.0.and.debug)print *,"all_ions",all_ions,tmp1 
     all_ions_lng=all_ions
     all_lecs_lng=all_lecs
#else
     all_ions_lng=ions
     all_lecs_lng=lecs
#endif    
     where(all_ions .lt. stride)all_ions=stride
     where(all_lecs .lt. stride)all_lecs=stride
     where(all_ions_lng .lt. stride)all_ions_lng=stride
     where(all_lecs_lng .lt. stride)all_lecs_lng=stride
     dimsfions(1)=sum(all_ions_lng/stride)
     dimsflecs(1)=sum(all_lecs_lng/stride)
     dimsfiions(1)=dimsfions(1)
     dimsfilecs(1)=dimsflecs(1)
     dimsfiions(2:7)=0
     dimsfilecs(2:7)=0
     if(rank .eq. 0 .and. debug) then 
        print *, "all_ions", all_ions_lng
        print *, "all_lecs", all_lecs_lng
        print *,"dimsions", sum(all_ions_lng/stride)
        print *,"dimslecs", sum(all_lecs_lng/stride)
     endif
        ! allocate arrays
     allocate(temporary_vec(max(max(ions/stride,lecs/stride),1)))       
#ifdef serIO
     if (rank .eq. 0) then
        allocate(temporary0_vec(max(maxval(all_ions_lng/stride), &
             maxval(all_lecs_lng/stride))))
        temporary0_vec(:)=0.
     endif
#endif  
     !
     !  Initialize FORTRAN predefined datatypes
     !
#ifdef serIO
     if (rank .eq. 0) then
        call h5open_f(error)
     endif
#else
     call h5open_f(error) 
#endif
     
     ! 
     ! Setup file access property list with parallel I/O access.
     !     
#ifdef MPI
#ifndef serIO
     call h5pcreate_f(H5P_FILE_ACCESS_F, plist_id, error)
     !      call H5Pset_sieve_buf_size_f(plist_id, 262144,error) 
     !      call H5Pset_alignment_f(plist_id, 524288, 262144,error)
     call h5pset_fapl_mpio_f(plist_id, mpi_comm_world, &
          FILE_INFO_TEMPLATE,error)
#endif
#endif

     !
     ! Create the file
     ! 
#ifdef MPI
#ifdef serIO
     if (rank .eq. 0) then
        call h5fcreate_f(fnameprt, H5F_ACC_TRUNC_F, file_id, error, &
             h5p_default_f,h5p_default_f)
     endif
#else
     call h5fcreate_f(fnameprt, H5F_ACC_TRUNC_F, file_id, error, &
          access_prp = plist_id)
     call h5pclose_f(plist_id, error)
#endif
#else
     call h5fcreate_f(fnameprt, H5F_ACC_TRUNC_F, file_id, error)
#endif      
          
     !
     ! Create the data space for the  dataset. 
     !     
#ifdef serIO
     if (rank .eq. 0) then
        do i=1,midvars
           call h5screate_simple_f(datarank, dimsfions(1), filespace(i), &
                error)
        enddo
        do i=midvars+1,nvars
           call h5screate_simple_f(datarank, dimsflecs(1), filespace(i), &
                error)
        enddo
     endif
#else
     do i=1,midvars
        call h5screate_simple_f(datarank, dimsfions(1), filespace(i), &
             error)
     enddo
     do i=midvars+1,nvars
        call h5screate_simple_f(datarank, dimsflecs(1), filespace(i), &
             error)
     enddo
#endif
     
     !
     ! Create property list for collective dataset write
     !
#ifndef serIO
#ifdef MPI
     call h5pcreate_f(H5P_DATASET_XFER_F, plist_id, error) 
     call h5pset_dxpl_mpio_f(plist_id, H5FD_MPIO_COLLECTIVE_F, error)
#endif
#endif
     
     
     !looping nvars
     do i=1,nvars

        varname=trim(dsetname(i))
        
        !define temporary_vec for every proc
        temporary_vec=0.
        !IONS
        if(varname.eq.'xi') temporary_vec(1:ions/stride)=p(1:ions:stride)%x &
        	 +mxcum
        if(varname.eq.'yi') temporary_vec(1:ions/stride)=p(1:ions:stride)%y &
             +mycum
        if(varname.eq.'zi') temporary_vec(1:ions/stride)=p(1:ions:stride)%z &
        	 +mzcum
        if(varname.eq.'ui') temporary_vec(1:ions/stride)=p(1:ions:stride)%u
        if(varname.eq.'vi') temporary_vec(1:ions/stride)=p(1:ions:stride)%v
        if(varname.eq.'wi') temporary_vec(1:ions/stride)=p(1:ions:stride)%w    
        if(varname.eq.'chi') temporary_vec(1:ions/stride)=p(1:ions:stride)%ch
        if(varname.eq.'gammai') temporary_vec(1:ions/stride)=sqrt(1. &
             +(p(1:ions:stride)%u**2+p(1:ions:stride)%v**2 &
             +p(1:ions:stride)%w**2))
        if(varname.eq.'indi') temporary_vec(1:ions/stride)=real(p(1:ions:stride) &
             %ind,4)
        if(varname.eq.'proci') temporary_vec(1:ions/stride)=real(p(1:ions:stride) &
             %proc,4)
        !ELECTRONS        
        if(varname.eq.'xe') temporary_vec(1:lecs/stride)=p(maxhlf+1:maxhlf &
             +lecs:stride)%x+mxcum
        if(varname.eq.'ye') temporary_vec(1:lecs/stride)=p(maxhlf+1:maxhlf &
             +lecs:stride)%y+mycum
        if(varname.eq.'ze') temporary_vec(1:lecs/stride)=p(maxhlf+1:maxhlf &
             +lecs:stride)%z+mzcum
        if(varname.eq.'ue') temporary_vec(1:lecs/stride)=p(maxhlf+1:maxhlf &
             +lecs:stride)%u
        if(varname.eq.'ve') temporary_vec(1:lecs/stride)=p(maxhlf+1:maxhlf &
             +lecs:stride)%v
        if(varname.eq.'we') temporary_vec(1:lecs/stride)=p(maxhlf+1:maxhlf &
             +lecs:stride)%w
        if(varname.eq.'che') temporary_vec(1:lecs/stride)=p(maxhlf+1:maxhlf+lecs:stride)%ch
        if(varname.eq.'gammae') temporary_vec(1:lecs/stride)=sqrt(1. +(p(maxhlf &
             +1:maxhlf+lecs:stride)%u**2+p(maxhlf+1:maxhlf +lecs:stride &
             )%v**2+p(maxhlf +1:maxhlf+lecs:stride)%w**2))
        if(varname.eq.'inde')temporary_vec(1:lecs/stride)=real(p(maxhlf+1:maxhlf &
             +lecs:stride)%ind,4)
        if(varname.eq.'proce')temporary_vec(1:lecs/stride)=real(p(maxhlf+1:maxhlf &
             +lecs:stride)%proc,4)
        
#ifdef serIO
        !
        ! Serial writing of prtl.tot.
        !        

        do procn=0,size0-1
           
           if (rank .eq. 0) then
              
              if (procn .ne. 0) then !receive from procn                                
#ifdef MPI
                 if(i .le. midvars) then
                    call MPI_Recv(temporary0_vec(1:all_ions_lng(procn+1)/stride) &
                         ,all_ions_lng(procn+1)/ &
                         stride,mpi_real,procn,1, &
                         MPI_COMM_WORLD,status,error)
                 else
                    call MPI_Recv(temporary0_vec(1:all_lecs_lng(procn+1)/stride) &
                         ,all_lecs_lng(procn+1)/ &
                         stride,mpi_real,procn,1, &
                         MPI_COMM_WORLD,status,error)
                 endif
#endif
              else           !procn eq 0
                 if (i .le. midvars) then
                    temporary0_vec(1:all_ions_lng(1)/stride)= &
                         temporary_vec(1:all_ions_lng(1)/stride)   
                 else
                    temporary0_vec(1:all_lecs_lng(1)/stride)= &
                         temporary_vec(1:all_lecs_lng(1)/stride)  
                 endif         
              endif
              
              !define count and offset
              if(i.le.midvars)countpart = max(all_ions_lng(procn+1)/ &
                   stride,1)
              if(i.gt.midvars)countpart = max(all_lecs_lng(procn+1)/ &
                   stride,1)
              
              offsetpart = 0
              if(procn .gt. 0) then 
                 if(i .le. midvars) then
                    offsetpart = sum(all_ions_lng(1:procn)/stride)
                 else
                    offsetpart = sum(all_lecs_lng(1:procn)/stride)
                 endif
              endif
              
              
              !hyperslab selection
#ifdef MPI
              ! this is necessary only if I close the filespace above      
              !               call h5dget_space_f(dset_id(i), filespace(i), error)               
              call h5sselect_hyperslab_f(filespace(i),H5S_SELECT_SET_F &
                   ,offsetpart,countpart,error)
#endif               
              
              if(debug)print *,rank,": selected hyperslab"
              
              
              !memory space
#ifdef MPI
              call h5screate_simple_f(1, countpart, memspace, error) 
#endif
                            
              !dataset creation/opening
              if (procn .eq. 0) then
                 call h5dcreate_f(file_id,varname,H5T_NATIVE_REAL, &
                      filespace(i),dset_id(i), error)
              else 
                 call h5dopen_f(file_id,varname,dset_id(i),error)
              endif
              
              
              !writing
#ifdef MPI
              if(i .le. midvars) then 
                 call h5dwrite_real_1(dset_id(i), H5T_NATIVE_REAL, &
                      temporary0_vec(1:all_ions_lng(procn+1)/stride),  &
                      dimsfiions,error,file_space_id = filespace(i),  &
                      mem_space_id=memspace)
              else 
                 call h5dwrite_real_1(dset_id(i), H5T_NATIVE_REAL, &
                      temporary0_vec(1:all_lecs_lng(procn+1)/stride),  &
                      dimsfilecs,error,file_space_id = filespace(i),  &
                      mem_space_id=memspace)
              endif
              
#else
              if(i .le. midvars) then 
                 call h5dwrite_f(dset_id(i), H5T_NATIVE_REAL, &
                      temporary_vec(1:max(ions/stride,1)),  &
                      dimsfiions,error)
              else 
                 call h5dwrite_f(dset_id(i), H5T_NATIVE_REAL, &
                      temporary_vec(1:max(lecs/stride,1)),  &
                      dimsfilecs,error)
              endif
#endif
              
              
              !closing dataset and memory space
              if(debug)print *, rank, ": ", i, "closing d"
              call h5dclose_f(dset_id(i), error)
              if(debug)print *, rank, ": ", i, "closing m"
              call h5sclose_f(memspace, error)
              
           else !rank ne 0
              if (procn .eq. rank) then
#ifdef MPI
                 if(i .le. midvars) then
                    call MPI_Send(temporary_vec,all_ions_lng(procn+1)/ &
                         stride,mpi_real,0,1, &
                         MPI_COMM_WORLD,error)
                 else
                    call MPI_Send(temporary_vec,all_lecs_lng(procn+1)/ &
                         stride,mpi_real,0,1, &
                         MPI_COMM_WORLD,error)
                 endif
#endif       
              endif
           endif               !rank ne 0
           
        enddo!procn
!else (finished serial writing)        
#else 
        !
        ! Parallel writing of prtl.tot.
        !

        !
        ! Create the dataset with default properties.
        !        
        call h5dcreate_f(file_id,varname, H5T_NATIVE_REAL, &
             filespace(i),dset_id(i), error)
        call h5sclose_f(filespace(i), error)

        !
        ! Each process defines dataset in memory and writes it to the hyperslab
        ! in the file. 
        !
        if(i.le.midvars)countpart = max(ions/stride,1)
        if(i.gt.midvars)countpart = max(lecs/stride,1)
        
        offsetpart = 0
        if(rank .gt. 0) then 
           if(i .le. midvars) then
              offsetpart = sum(all_ions_lng(1:rank)/stride)
           else
              offsetpart = sum(all_lecs_lng(1:rank)/stride)
           endif
        endif
        
#ifdef MPI
        call h5screate_simple_f(1, countpart, memspace, error)         
        call h5dget_space_f(dset_id(i), filespace(i), error)        
        call h5sselect_hyperslab_f(filespace(i),H5S_SELECT_SET_F &
             ,offsetpart,countpart, error)        
#endif
        
#ifdef MPI
        if(i .le. midvars) then 
           call h5dwrite_real_1(dset_id(i), H5T_NATIVE_REAL, &
                temporary_vec(1:max(ions/stride,1)), dimsfiions,error &
                ,file_space_id = filespace(i), mem_space_id=memspace &
                ,xfer_prp = h5p_default_f)
        else 
           call h5dwrite_real_1(dset_id(i), H5T_NATIVE_REAL, &
                temporary_vec(1:max(lecs/stride,1)), dimsfilecs,error &
                ,file_space_id = filespace(i), mem_space_id=memspace &
                ,xfer_prp = h5p_default_f)
        endif
        
#else
        if(i .le. midvars) then 
           call h5dwrite_f(dset_id(i), H5T_NATIVE_REAL, &
                temporary_vec(1:max(ions/stride,1)), dimsfiions,error)
        else 
           call h5dwrite_f(dset_id(i), H5T_NATIVE_REAL, &
                temporary_vec(1:max(lecs/stride,1)), dimsfilecs,error)
        endif
#endif
        
	if(debug)print *, rank, ": ", i, "closing m"
#ifdef MPI
        call h5sclose_f(memspace,error)
#endif
        if(debug)print *, rank, ": ", i, "closing d"
        call h5dclose_f(dset_id(i), error)
!endif finished parallel writing
#endif 
       
        
        !
        ! Close dataspaces.
        !
        if(debug)print *, rank, ": ", i, "closing s" 
#ifdef MPI
#ifdef serIO    
        if (rank .eq. 0) then
           call h5sclose_f(filespace(i), error)
        endif
#else
        call h5sclose_f(filespace(i), error)
#endif
#endif
        
     enddo!nvars

     
#ifdef serIO
     if (rank .eq. 0) then
        if(debug)print *, rank, ": ", "closing f"
        call h5fclose_f(file_id, error)
        if(debug)print *, rank, ": ", "closing h5"
        call h5close_f(error)
        deallocate(temporary0_vec)
     endif
#else
     if(debug)print *, rank, ": ", "closing p"
#ifdef MPI
     call h5pclose_f(plist_id, error)
#endif
     if(debug)print *, rank, ": ", "closing f"
     call h5fclose_f(file_id, error)
     if(debug)print *, rank, ": ", "closing h5"
     call h5close_f(error)
#endif
     
     if(debug) print *, rank,": finished writing particles"
     
     deallocate(temporary_vec)

127  continue
     
#ifdef MPI
#ifndef serIO
     call MPI_Info_free(FILE_INFO_TEMPLATE,error)
#endif
#endif
		
     
  endif !if lap ge pltstart
	
	
end subroutine output_tot



!-------------------------------------------------------------------------------
! 						subroutine output_hug					 
!														
! Save flds.hug. files for fluid data with different spatial resolution  					
!
!-------------------------------------------------------------------------------
subroutine output_hug()
 
  implicit none
  
  ! local variables
  integer ::ind,kglob,kloc,kstart,kfinish,nz,info,n,i,j,k,ierr
  integer ::jstart,jfinish,ny,iv,mz1,my1,mx1
  real bx0,by0,bz0,ex0,ey0,ez0,bxd,byd,bzd,exd,eyd,ezd,curx0, &
       cury0,curz0
  integer ::error ! Error flags
  integer ::nvars !, midvars
!  integer ::all_ions(size0), all_lecs(size0)
!  integer*8 all_ions_lng(size0), all_lecs_lng(size0)
  character(len=6) :: dsetname(40), varname  
  integer ::datarank
#ifdef HDF5
  integer(HID_T) :: file_id ! File identifier 
  integer(HID_T) :: dset_id(40) ! Dataset identifier 
  integer(HID_T) :: filespace(40) ! Dataspace identifier in file 
  integer(HID_T) :: memspace ! Dataspace identifier in memory
  integer(HID_T) :: plist_id ! Property list identifier
  integer(HSIZE_T) dimsf(3)
  integer(HSIZE_T) dimsfi(7)
  integer(HSIZE_T), DIMENSION(3) :: count  
  integer(HSSIZE_T), DIMENSION(3) :: offset
#endif 
  integer ::FILE_INFO_TEMPLATE !for MPI file info
  integer ::procn
  integer ::my1list(size0),mz1list(size0) 
  real*4, allocatable :: temporary0(:,:,:)!, temporary0_vec(:)   
  integer :: status(MPI_STATUS_SIZE)
  integer ::istep0
  logical :: saverank
  integer ::token,request,status2(MPI_STATUS_SIZE),request1 &
      ,status1(MPI_STATUS_SIZE)
  logical :: writecurr
  
  !user input: array of variables to be saved for flds.hug.
  !usage:
  !-the first entry should have the largest number of characters
  !(if needed, supplemented with empty spaces)
  !to be chosen from:
  !dens,densi: density and density of ions
  !ex,ey,ez,bx,by,bz,jx,jy,jz:fields and currents
  !v3x,v3xi,v3y,v3yi,v3z,v3zi:3-velocity
  !v4x,v4xi,v4y,v4yi,v4z,v4zi:4-velocity (momentum)
  !ken,keni:particle total energy
#ifdef twoD
  if (sigma .eq. 0. .and. pcosthmult .eq. 0) then 
     nvars=4
     dsetname(1:nvars)=(/&
          'dens  ','ex    ','ey    ','bz    '/)
  else
     nvars=7
     dsetname(1:nvars)=(/&
          'dens  ','ex    ','ey    ','ez    ','bx    ','by    ','bz    '/)
  endif
#else
  nvars=7
  dsetname(1:nvars)=(/&
       'dens  ','ex    ','ey    ','ez    ','bx    ','by    ','bz    '/)
#endif
  
  istep0=istep

  if(lap .ge.pltstart .and. modulo((lap-pltstart),torqint) .eq.0)then

     istep=istep1
     
     !timestamp of the saved output
     ind=(lap-pltstart)/torqint	    
     write(fnamefld,"(a16,i3.3)") "output/flds.hug.",ind     
     write(indchar,"(i3.3)")ind     
     if (rank.eq.0)print *,"",fnamefld
     
     !just in case destroy the file
     if(rank .eq. 0) then 
        open(unit=11,file=fnamefld,form='unformatted')
        close(11,status='delete')
     endif
     
     call MPI_BARRIER(MPI_COMM_WORLD,ierr)     

     
     !specify FILE_INFO_TEMPLATE for parallel output
#ifndef serIO 
#ifdef MPI     
     call MPI_Info_create(FILE_INFO_TEMPLATE,error)
     call MPI_Info_set(FILE_INFO_TEMPLATE, "access_style",  &
          "write_once",error)
     call MPI_Info_set(FILE_INFO_TEMPLATE,  &
          "collective_buffering", "true",error)
     call MPI_Info_set(FILE_INFO_TEMPLATE, "cb_block_size", &
          "4194304",error)
     call MPI_Info_set(FILE_INFO_TEMPLATE, "cb_buffer_size", &
          "16777216",error)
     call MPI_Info_set(FILE_INFO_TEMPLATE, "cb_nodes", &
          "1",error)
#endif
#endif


     !FLDS.HUG. saving

     !compute global size of fluid arrays to be saved
     dimsf(1)=mx0/istep
     dimsf(2)=my0/istep  
#ifndef twoD
     dimsf(3)=mz0/istep
#else
     dimsf(3)=1
#endif 
     dimsfi(1:3)=dimsf(1:3)
     dimsfi(4:7)=0
     datarank=3
     
     
     !computing my1 and mz1 for each processor
     kstart=3
     kfinish=mz-3
     if(rank/sizey .eq. 0) kstart=1
     if(rank/sizey .eq. sizez-1) kfinish=mz
     
1110  if(modulo(kstart+(rank/sizey)*(mzall-5)-1*0,istep) .eq. 0) goto 1120
     kstart=kstart+1
     goto 1110
1120  continue
     
1130  if(modulo(kfinish+(rank/sizey)*(mzall-5)-1*0,istep) .eq. 0) goto 1140
     kfinish=kfinish-1
     goto 1130
1140  continue
     
     mz1=(kfinish-kstart)/istep+1
#ifdef twoD
     kstart=1
     kfinish=1
     mz1=1
#endif
     
     jstart=3
     jfinish=my-3
     if(modulo(rank,sizey) .eq. 0) jstart=1
     if(modulo(rank,sizey) .eq. sizey-1) jfinish=my
     
1210  if(modulo(jstart+modulo(rank,sizey)*(myall-5)-1*0,istep) .eq. 0) goto 1220
     jstart=jstart+1
     goto 1210
1220  continue
     
1230  if(modulo(jfinish+modulo(rank,sizey)*(myall-5)-1*0,istep) .eq. 0)  goto 1240
     jfinish=jfinish-1
     goto 1230
1240  continue
     
     my1=(jfinish-jstart)/istep+1
  
! determine whether the current processor should save any cell
     if (kfinish .lt. kstart) then
        kfinish=1
        kstart=1
        mz1=0
     endif
     if (jfinish .lt. jstart) then
        jfinish=1
        jstart=1
        my1=0
     endif
   
     
     !allocate temporary1 on every processor
     if(debug) print *, rank,"in output b4 allocate"
     allocate(temporary1(mx0/istep,max(my1,1),max(mz1,1)))


     !getting information about the size of temporary1
#ifdef MPI
     if(debug) print *, rank,": in output, b4 allgather", my1
     call mpi_allgather(my1,1,mpi_integer,my1list,1,mpi_integer &
          ,mpi_comm_world, error)
     if(debug) print *, rank,": in output, b4 allgather", mz1
     call mpi_allgather(mz1,1,mpi_integer,mz1list,1,mpi_integer &
          ,mpi_comm_world, error) 
#ifdef serIO       
     if (rank .eq. 0) then
        allocate(temporary0(mx0/istep,maxval(my1list), &
             maxval(mz1list)))
        temporary0(:,:,:)=0.
     endif
#endif
#endif

!decide if I need to write currents to disk
     writecurr=.false.
     do iv=1,nvars
        varname=trim(dsetname(iv))
        if (varname .eq. 'jx' .or. varname .eq. 'jy') writecurr=.true. 
     enddo
     
     if (writecurr) then

        !install a semaphore        
#ifdef MPI
        if(rank .gt. 0) then 
           token=rank-1
           !wait for the message from the previous rank
           call mpi_irecv(token,1, &
                mpi_integer,rank-1,100,MPI_Comm_WORLD,request,ierr)
           call mpi_wait(request,status2,ierr)
           
           if (debug)print *,"rank ",rank," received the token from ",rank-1
           
           if(modulo(rank,15).ne.0 .and. rank .ne. size0-1) then
              token=rank
              call mpi_isend(token,1, &
                   mpi_integer,rank+1,100,MPI_Comm_WORLD,request1,ierr)
              call mpi_wait(request1,status1,ierr)
           endif
        endif
3004    continue
#endif

     !writing currents (all procs)
     open(unit=7,file=fenlargeloc, form='unformatted')
     rewind 7
     write(7) (((curx(i,j,k),i=1,mx),j=1,my),k=1,mz), &
          (((cury(i,j,k),i=1,mx),j=1,my),k=1,mz), &
          (((curz(i,j,k),i=1,mx),j=1,my),k=1,mz)
     close(7)

        !send the message to the next rank, if not the last        
#ifdef MPI
        if(modulo(rank,15).eq.0 .and. rank .ne. size0-1) then
           token=rank
           call mpi_isend(token,1, &
                mpi_integer,rank+1,100,MPI_Comm_WORLD,request,ierr)
        endif
2121    continue
        
        call mpi_barrier(MPI_COMM_WORLD,ierr)
#endif

    endif !writecurr
     
     !
     !  Initialize FORTRAN predefined datatypes
     !
#ifdef serIO
     if (rank .eq. 0) then
        call h5open_f(error)
     endif
#else
     call h5open_f(error) 
#endif
     
     ! 
     ! Setup file access property list if parallel I/O access.
     !     
#ifndef serIO
#ifdef MPI
     call h5pcreate_f(H5P_FILE_ACCESS_F, plist_id, error)
     !      call H5Pset_sieve_buf_size_f(plist_id, 262144,error) 
     !      call H5Pset_alignment_f(plist_id, 524288, 262144,error)
     call h5pset_fapl_mpio_f(plist_id, mpi_comm_world, &
          FILE_INFO_TEMPLATE,error)
#endif
#endif
     
     !
     ! Create the file
     !
#ifdef MPI
#ifdef serIO
     if (rank .eq. 0) then
        call h5fcreate_f(fnamefld, H5F_ACC_TRUNC_F, file_id, error, &
             h5p_default_f,h5p_default_f)
     endif
#else
     call h5fcreate_f(fnamefld, H5F_ACC_TRUNC_F, file_id, error, &
          access_prp = plist_id)
     call h5pclose_f(plist_id, error)
#endif
#else
     call h5fcreate_f(fnamefld, H5F_ACC_TRUNC_F, file_id, error)
#endif  

     !
     ! Create the dataspace for the dataset
     !
#ifdef serIO
     if (rank .eq. 0) then
        do i=1,nvars
           call h5screate_simple_f(datarank, dimsf, filespace(i), error)
        enddo
     endif
#else
     do i=1,nvars
        call h5screate_simple_f(datarank, dimsf, filespace(i), error)
     enddo
#endif
     
     !
     ! Loop over variables
     !
     do iv=1,nvars

        varname=trim(dsetname(iv))

        !
        ! Define the array
        ! 
        !pre-processing of some field quantities
        !total density and ion density
        if(varname.eq.'dens' ) call meanq_fld_cur('tdens')
        if(varname.eq.'densi') call meanq_fld_cur('idens')
        !velocities, momenta and energy
        if(varname.eq.'v3x'  ) call meanq_fld_cur('tbetx')
        if(varname.eq.'v3y'  ) call meanq_fld_cur('tbety')
        if(varname.eq.'v3z'  ) call meanq_fld_cur('tbetz') 
        if(varname.eq.'v3xi' ) call meanq_fld_cur('ibetx')
        if(varname.eq.'v3yi' ) call meanq_fld_cur('ibety')
        if(varname.eq.'v3zi' ) call meanq_fld_cur('ibetz')
        if(varname.eq.'v4x'  ) call meanq_fld_cur('tmomx')
        if(varname.eq.'v4y'  ) call meanq_fld_cur('tmomy')
        if(varname.eq.'v4z'  ) call meanq_fld_cur('tmomz')
        if(varname.eq.'v4xi' ) call meanq_fld_cur('imomx')
        if(varname.eq.'v4yi' ) call meanq_fld_cur('imomy')
        if(varname.eq.'v4zi' ) call meanq_fld_cur('imomz')
        if(varname.eq.'ken'  ) call meanq_fld_cur('tener')
        if(varname.eq.'keni' ) call meanq_fld_cur('iener')
        !electric currents
        if(varname.eq.'jx'.or.varname.eq.'jy'.or.varname.eq.'jz') then 
!        if(varname.eq.'jx'.or.varname.eq.'jy') then
           !     restore currents from file
           open(unit=7,file=fenlargeloc, form='unformatted')
           rewind 7
           if(debug.and.rank.eq.0) print *, rank,":reading back current file"
           read(7)(((curx(i,j,k),i=1,mx),j=1,my),k=1,mz), &
                (((cury(i,j,k),i=1,mx),j=1,my),k=1,mz), &
                (((curz(i,j,k),i=1,mx),j=1,my),k=1,mz)
           if(debug .and. rank.eq.0)print *, rank,"resetting arrays"           
        endif
        
        ! define temporary1 for every processor

        do k=kstart, kfinish, istep 
           nz=(k-kstart)/istep+1
           do j=jstart, jfinish,istep 
              ny=(j-jstart)/istep+1
              do i=1,mx/istep
       
                 select case (varname)
                 case('dens','densi')
                    temporary1(i,ny,nz)=real(curx(i*istep,j,k),4)
                 case('ex')
                    call interpolfld(real(i*istep),real(j), &
                         real(k),bx0,by0,bz0,ex0,ey0,ez0,bxd, &
                         byd,bzd,exd,eyd,ezd) 
                    temporary1(i,ny,nz)=real(ex0+exd,4)
                 case('ey')
                    call interpolfld(real(i*istep),real(j), &
                         real(k),bx0,by0,bz0,ex0,ey0,ez0,bxd, &
                         byd,bzd,exd,eyd,ezd) 
                    temporary1(i,ny,nz)=real(ey0+eyd,4)
                 case('bz')
                    call interpolfld(real(i*istep),real(j), &
                         real(k),bx0,by0,bz0,ex0,ey0,ez0,bxd, &
                         byd,bzd,exd,eyd,ezd) 
                    temporary1(i,ny,nz)=real(bz0+bzd,4)
                 case('ez')
                    call interpolfld(real(i*istep),real(j), &
                         real(k),bx0,by0,bz0,ex0,ey0,ez0,bxd, &
                         byd,bzd,exd,eyd,ezd) 
                    temporary1(i,ny,nz)=real(ez0+ezd,4)
                 case('bx')
                    call interpolfld(real(i*istep),real(j), &
                         real(k),bx0,by0,bz0,ex0,ey0,ez0,bxd, &
                         byd,bzd,exd,eyd,ezd) 
                    temporary1(i,ny,nz)=real(bx0+bxd,4)
                 case('by')
                    call interpolfld(real(i*istep),real(j), &
                         real(k),bx0,by0,bz0,ex0,ey0,ez0,bxd, &
                         byd,bzd,exd,eyd,ezd) 
                    temporary1(i,ny,nz)=real(by0+byd,4)
                 case('jx')
                    !                   call interpolcurrent(real(i*istep),real(j), &
                    !                        real(k),curx0, cury0, curz0)
                    curx0=curx(i*istep,j,k)
                    temporary1(i,ny,nz)=real(curx0,4)
                 case('jy')
                    !                   call interpolcurrent(real(i*istep),real(j), &
                    !                        real(k),curx0, cury0, curz0)
                    cury0=cury(i*istep,j,k)
                    temporary1(i,ny,nz)=real(cury0,4)
                 case('jz')
                    !                   call interpolcurrent(real(i*istep),real(j), &
                    !                        real(k),curx0, cury0, curz0)
                    curz0=curz(i*istep,j,k)
                    temporary1(i,ny,nz)=real(curz0,4)
                 case('ken','v3x','v3y','v3z','v4x','v4y','v4z',&
                      'keni','v3xi','v3yi','v3zi','v4xi','v4yi','v4zi')
                    temporary1(i,ny,nz)=real(curx(i*istep,j,k),4)
                 end select
                 
              enddo
           enddo
        enddo
        

#ifdef serIO
        !
        ! Serial writing
        !
        if (rank.eq.0.and.debug) print *,"serial writing for .hug."

        do procn=0,size0-1

           if (rank .eq. 0) then
              if (procn .eq. 0) then
                 call h5dcreate_f(file_id,varname,H5T_NATIVE_REAL, &
                      filespace(iv),dset_id(iv), error)
                 call h5dclose_f(dset_id(iv), error)
              endif
           endif
           
           ! Determine which processors should write to file
           saverank=.true.
           if (my1list(procn+1)*mz1list(procn+1) .eq. 0) then
              saverank=.false.
           endif

           ! write only if saverank eq. .true.
           if (saverank) then
              
           if (rank .eq. 0) then
              
              if (procn .ne. 0) then !receive from procn                                
#ifdef MPI
                 call MPI_Recv(temporary0(1:mx0/istep,1:my1list(procn+1), &
                      1:mz1list(procn+1)),mx0/istep*my1list(procn+1) &
                      *mz1list(procn+1),mpi_real,procn,1, &
                      MPI_COMM_WORLD,status,error)
#endif
              else           !procn eq 0
                 temporary0(1:mx0/istep,1:my1,1:mz1)=temporary1               
              endif
              
              !     computing the offset
              count(1) = dimsf(1)
              count(2) = my1list(procn+1)
              count(3) = mz1list(procn+1)
              
              offset(1) = 0
              offset(2) = 0
              offset(3) = 0
              !     z offset
              if(procn/sizey .gt. 0) then 
                 offset(3)=sum(mz1list(1:(procn/sizey)*sizey:sizey))
              endif
#ifdef twoD 
              offset(3)=0
#endif
              !     y offset
              if(modulo(procn,sizey) .gt. 0) then 
                 offset(2)=sum(my1list((procn/sizey)*sizey+1:procn))
              endif
              
              !hyperslab selection
#ifdef MPI
              ! this is necessay only if I close the filespace above      
              !               call h5dget_space_f(dset_id(iv), filespace(iv), error)               
              call h5sselect_hyperslab_f(filespace(iv),H5S_SELECT_SET_F &
                   ,offset,count,error)
#endif               
              
              !memory space
#ifdef MPI
              call h5screate_simple_f(datarank, count, memspace, error) 
#endif              
              
              !create or open the dataset
              call h5dopen_f(file_id,varname,dset_id(iv),error)
              

              !writing
              if(debug)print *,rank,": before h5dwrite","filespace", &
                   filespace(iv),"memspace",memspace
              
#ifdef MPI
              call h5dwrite_f(dset_id(iv), H5T_NATIVE_REAL,  &
                   temporary0(1:mx0/istep,1:my1list(procn+1), &
                   1:mz1list(procn+1)),dimsfi,error,mem_space_id= &
                   memspace,file_space_id = filespace(iv))
              !               call h5dwrite_real_1(dset_id(iv), H5T_NATIVE_REAL, 
              !     &                   temporary0(:,:,:),dimsfi,error,file_space_id = 
              !     &                   filespace(iv), mem_space_id=memspace)
#else
              call h5dwrite_f(dset_id(iv),H5T_NATIVE_REAL, &
                   temporary(:,:,:) ,dimsfi,error)
#endif 
              
              !closing dataset and memory space
              if(debug)print *, rank, ": ", iv, "closing d"
              call h5dclose_f(dset_id(iv), error)
              if(debug)print *, rank, ": ", "closing m"
              call h5sclose_f(memspace, error)

           else !rank ne 0
              
              if (procn .eq. rank) then
#ifdef MPI
                 call MPI_Send(temporary1,mx0/istep*my1 &
                      *mz1,mpi_real,0,1, &
                      MPI_COMM_WORLD,error) 
#endif             
              endif
              
           endif !rank ne 0
              
           endif !saverank .eq. .true.
           
        enddo!procn
        
#else 
! finished serial writing

        !
        ! Parallel writing
        !     
        if (rank.eq.0.and.debug) print *,"parallel writing of .hug. files"
           
        ! Determine which processors should write to file
        saverank=.true.
        if (my1*mz1 .eq. 0) then
           saverank=.false.
        endif

        !
        ! Create the dataset with default properties
        !
        call h5dcreate_f(file_id, varname, H5T_NATIVE_REAL, &
             filespace(iv),dset_id(iv), error)
        call h5sclose_f(filespace(iv), error)
        call h5dget_space_f(dset_id(iv), filespace(iv), error)

        !
        ! Create property list for collective dataset write
        !
#ifdef MPI
        call h5pcreate_f(H5P_DATASET_XFER_F, plist_id, error) 
        call h5pset_dxpl_mpio_f(plist_id, H5FD_MPIO_COLLECTIVE_F, error)
#endif
           
        !
        ! Each process defines dataset in memory and writes it to the hyperslab
        ! in the file. 
        !
        count(1) = dimsf(1)
        count(2) = my1
        count(3) = mz1
        
        offset(1) = 0
        offset(2) = 0
        offset(3) = 0
        !z offset
        if(rank/sizey .gt. 0) then 
           offset(3) = sum(mz1list(1:(rank/sizey)*sizey:sizey))
        endif
#ifdef twoD 
        offset(3)=0
#endif
        !now get the y offset
        if(modulo(rank,sizey) .gt. 0) then 
           offset(2)=sum(my1list((rank/sizey)*sizey+1:rank))
        endif
        
        ! write only if saverank .eq. .true.
        if (saverank) then

#ifdef MPI
        call h5screate_simple_f(datarank, count, memspace, error) 
#endif
        
        ! 
        ! Select hyperslab in the file.
        !
#ifdef MPI              
        call h5sselect_hyperslab_f(filespace(iv),H5S_SELECT_SET_F,offset &
             ,count, error)
#endif
        
        if(debug)print *,rank,": before h5dwrite","filespace",filespace(iv &
             ),memspace
        
#ifdef MPI
        call h5dwrite_real_1(dset_id(iv), H5T_NATIVE_REAL, temporary1(:,: &
             ,:),dimsfi,error,file_space_id = filespace(iv), mem_space_id &
             =memspace,xfer_prp = h5p_default_f)
#else
        call h5dwrite_f(dset_id(iv),H5T_NATIVE_REAL,TEMPORARY1(:,: &
             ,:) ,dimsfi,error)
#endif      
        endif ! saverank .eq. .true.
     
#endif 
 !finished parallel writing    
        

        !
        !Closing dataset and dataspace    
        !
#ifdef serIO    
        if (rank .eq. 0) then
#ifdef MPI
           if(debug)print *, rank, ": ", iv, "closing s" 
           call h5sclose_f(filespace(iv), error)
#endif
        endif
#else
        if(debug)print *, rank, ": ", iv, "closing d"
        call h5dclose_f(dset_id(iv), error)
#ifdef MPI
        if(debug)print *, rank, ": ", iv, "closing s" 
        call h5sclose_f(filespace(iv), error)     
#endif
#endif

     enddo!nvars
     
     
     !
     ! Closing memory space, file and hdf5
     !
#ifdef serIO
     if (rank .eq. 0) then
        if(debug)print *, rank, ": ", "closing f"
        call h5fclose_f(file_id, error)
        if(debug)print *, rank, ": ", "closing h5"
        call h5close_f(error)
        deallocate(temporary0)
     endif
#else
#ifdef MPI
     if (saverank) then
        if(debug)print *, rank, ": ", "closing m"
        call h5sclose_f(memspace, error)
     endif !saverank .eq. .true.
     if(debug)print *, rank, ": ", "closing p"
     call h5pclose_f(plist_id, error)
#endif
     if(debug)print *, rank, ": ", "closing f"
     call h5fclose_f(file_id, error)
     if(debug)print *, rank, ": ", "closing h5"
     call h5close_f(error)    
#endif
     
     if(debug) print *, rank,": finished writing .hug. fields"
     
     deallocate(temporary1)

     !-----------------------------------------
     !for huge files no need to write particles
     !-----------------------------------------

  endif !torqint

  istep=istep0

end subroutine output_hug


!-------------------------------------------------------------------------------
! 						subroutine write_test_e				 
! Jaehong, updated 9/2013																		
! Track preselected electrons (from restart in former simulation) 
! Serially collected HDF files are created in output/tracking/ directory.
!  							
!-------------------------------------------------------------------------------

subroutine write_test_e()

 	 implicit none
 	 
 	 integer i,n,n1,statusfile,procn
 	 real bx0,by0,bz0,ex0,ey0,ez0,dummy,gam
 	 character fnametst*20
 	 integer num1,proc1 	 
 	 ! number of detected particles in each rank
 	 integer nsave
 	 integer, allocatable :: nsavelist(:)
 	 ! number of detected particles (in all ranks)
 	 integer nsave0
 	 ! particle data (local)
 	 real, allocatable :: lecdata(:,:)
 	 ! particle data (global)
	 real, allocatable :: lecdata0(:,:)
	 integer, allocatable :: selectpart(:)
	 
	 integer istart, ifinal
	 character (len=50) :: fname
	 character(len=10) dsetname(20)
	 integer ::datarank, error, ierr
	 integer :: status(MPI_STATUS_SIZE)
	 
	 
#ifdef HDF5
	 integer(HID_T) :: file_id       ! File identifier 
	 integer(HID_T) :: dset_id       ! Dataset identifier 
	 integer(HID_T) :: dspace_id
	 integer(HSIZE_T) data_dims(2)
#endif
	 
	 allocate(nsavelist(size0))  
	
 	 if(prt_first_lec) then  !do it only at the first time 
    
 	 	!  count the number of selected particles in "select.testprt.dat"
 	 	write(fnametst,"(a17)") "select.test.e.dat"
 	 	open(unit=43, file=fnametst,form='formatted')
 	 	if(rank .eq. 0) print *,"loading particles from ",fnametst
 	 	sele=0
 	 	! read selected partcles; index, processor, gamma
301  	read(43,*,err=398,end=399) num1, proc1, gam
      	! number of electrons
        sele=sele+1
        goto 301
398  	print *,"error reading the test.sel. file"
399  	continue
     	close(43)
     	     	
     	allocate(selprte(sele))
     
     	open(unit=43, file=fnametst,form='formatted')     
     	sele=0
401  	read(43,*,err=498,end=499) num1,proc1,gam
			sele=sele+1
        	selprte(sele)%ind=num1
        	selprte(sele)%proc=proc1
        goto 401
498  	print *,"error reading the test.sel. file"
499  	continue

    	if(rank .eq. 0) print*, "selected electrons : ", sele
         	
  endif ! only for the first time
  
  
  		 allocate(selectpart(sele))
     	       	     
     	 nsave=0 ! number of detected electrons
     	 do n=maxhlf+1,maxhlf+lecs
     	 	if(p(n)%ind .ne. 0) then
     	 		do n1=1,sele
     	 			! detect condition
       				if(p(n)%ind .eq. selprte(n1)%ind & 
       					.and. p(n)%proc .eq. selprte(n1)%proc) then
        				nsave=nsave+1    
         				selectpart(nsave)=n
         			endif
         		enddo
         	endif
         enddo
     
      	allocate(lecdata(nsave,16))

      	call mpi_allgather(nsave,1,mpi_integer, &
                   nsavelist,1,mpi_integer,mpi_comm_world,ierr)
        call mpi_allreduce(nsave,nsave0,1,mpi_integer,mpi_sum &
			,mpi_comm_world,ierr)         
              
        if(nsave0 .ne. 0) then
        ! particle data (global)
        	allocate(lecdata0(nsave0,16))
        else
        	allocate(lecdata0(1,16))
        endif
	
        do i=1, nsave      
           	n=selectpart(i)
	       	call interpolfld(real(p(n)%x),real(p(n)%y),p(n)%z,bx0,by0, &
                  	 bz0,ex0,ey0,ez0,dummy,dummy,dummy,dummy,dummy,dummy)
     	
            lecdata(i,1)=lap
            lecdata(i,2)=p(n)%ind
            lecdata(i,3)=p(n)%proc
            lecdata(i,4)=p(n)%x+mxcum
            lecdata(i,5)=p(n)%y+mycum
            lecdata(i,6)=p(n)%z+mzcum
            lecdata(i,7)=p(n)%u
            lecdata(i,8)=p(n)%v
            lecdata(i,9)=p(n)%w
            lecdata(i,10)=sqrt(1.+(p(n)%u**2+p(n)%v**2+p(n)%w**2))	
            lecdata(i,11)=ex0
            lecdata(i,12)=ey0
            lecdata(i,13)=ez0
            lecdata(i,14)=bx0
            lecdata(i,15)=by0
            lecdata(i,16)=bz0
            
        enddo      

      ! serial writing	
      	  
      	if(nsave0 .ne. 0) then  
     
      	do procn=0, size0-1
           
           	if (rank .eq. 0) then
           	   
           		if (procn .ne. 0) then !receive from procn         
           		
                  	istart=sum(nsavelist(1:procn))+1
                  	ifinal=istart+nsavelist(procn+1)-1                   
                         
                    call MPI_Recv(lecdata0(istart:ifinal,1:16) &
                         ,nsavelist(procn+1)*16,mpi_real,procn,1, &
                         MPI_COMM_WORLD,status,ierr)    
                
                else  !procn eq 0          	                   	    	    
                          	 
                    lecdata0(1:nsave,1:16)=lecdata(1:nsave,1:16)              
              	endif
              
            else !rank ne 0
            	if (procn .eq. rank) then        
              	  	call MPI_Send(lecdata(1:nsave,1:16), &
                 	 	nsave*16,mpi_real,0,1,MPI_COMM_WORLD,ierr)                 	 	
                endif
              
            endif   !rank ne 0
	
       enddo !procn 
       
          
     else	!if nsvai0 .ne. 0
     	   nsave0=1
       	   if (rank .eq. 0) then
       	   	   lecdata0(:,:)=0
       	   endif
     endif
	   
       call mpi_barrier(MPI_COMM_WORLD,ierr)       
    	   
       
       if(rank .eq. 0) then 	      	   
       	   
       	   !******************************************************
       	   !****** 	save partcle tracking
       	   !******************************************************

       	   write(fname,"(a31,i7.7)") "./output/tracking_elec/testprt.",lap
       	   print *,"",fname        
        
       	   !   Initialize FORTRAN predefined datatypes
			
       	   call h5open_f(error) 
       	   call h5fcreate_f(fname,H5F_ACC_TRUNC_F,file_id,error)
	
   		   dsetname(1)="elec"
		   datarank=2
		
		   data_dims(1)=nsave0
		   data_dims(2)=16
		
		   call h5screate_simple_f(datarank, data_dims, dspace_id, error)
		   call h5dcreate_f(file_id, dsetname(1), H5T_NATIVE_REAL &
				,dspace_id,dset_id,error)			
		   call h5dwrite_f(dset_id, H5T_NATIVE_REAL,lecdata0, data_dims, error)
		   call h5dclose_f(dset_id,error)
		   call h5sclose_f(dspace_id,error)
	   
		   call h5fclose_f(file_id, error)
		   call h5close_f(error)
	   
	   endif ! rank 0
	
	   if(allocated(selectpart)) deallocate(selectpart)
	   if(allocated(lecdata)) deallocate(lecdata)
	   if(allocated(lecdata0)) deallocate(lecdata0)

	   
	   call mpi_barrier(MPI_COMM_WORLD,ierr)	 	  
	      
  
end subroutine write_test_e

!---------------------------------------------------------------------


subroutine write_test_p()

 	 implicit none
 	 
 	 integer i,n,n1,statusfile,procn
 	 real bx0,by0,bz0,ex0,ey0,ez0,dummy,gam
 	 character fnametst*20
 	 integer num1,proc1 	 
 	 ! number of detected particles in each rank
 	 integer nsavi
 	 integer, allocatable :: nsavilist(:)
 	 ! number of detected particles (in all ranks)
 	 integer nsavi0
 	 ! particle data (local)
 	 real, allocatable :: iondata(:,:)
 	 ! particle data (global)
	 real, allocatable :: iondata0(:,:)
	 integer, allocatable :: selectpart(:)
	 
	 integer istart, ifinal
	 character (len=50) :: fname
	 character(len=10) dsetname(20)
	 integer ::datarank, error, ierr
	 integer :: status(MPI_STATUS_SIZE)
	 
	 
#ifdef HDF5
	 integer(HID_T) :: file_id       ! File identifier 
	 integer(HID_T) :: dset_id       ! Dataset identifier 
	 integer(HID_T) :: dspace_id
	 integer(HSIZE_T) data_dims(2)
#endif
	 
	 allocate(nsavilist(size0))  

 
 	 if(prt_first_ion) then  !first time 
    
 	 	!    count the number of selected particles in "select.testprt.dat"
 	 	write(fnametst,"(a17)") "select.test.p.dat"
 	 	open(unit=43, file=fnametst,form='formatted')
 	 	if(rank .eq. 0) print *,"loading particles from ",fnametst
 	 	seli=0
 	 	! read selected partcles; index, processor, gamma
301  	read(43,*,err=398,end=399) num1, proc1, gam
   		! number of ions
        seli=seli+1
        goto 301
398  	print *,"error reading the test.sel. file"
399  	continue
     	close(43)
     	     	
     	allocate(selprti(seli))
       
     	open(unit=43, file=fnametst,form='formatted')     
     	seli=0
401  	read(43,*,err=498,end=499) num1,proc1,gam
    	   	seli=seli+1
        	selprti(seli)%ind=num1
        	selprti(seli)%proc=proc1 
        goto 401
498  	print *,"error reading the test.sel. file"
499  	continue

		if(rank .eq. 0) print*, "selected ions : ", seli
     
  endif ! only for the first time
  
  		allocate(selectpart(seli))
        
     	nsavi=0	     ! number of detected ions
     	do n=1,ions
     		if(p(n)%ind .ne. 0) then 
     			do n1=1,seli
     				! detect condition
     				if(p(n)%ind .eq. selprti(n1)%ind & 
     					.and. p(n)%proc .eq. selprti(n1)%proc) then
     					nsavi=nsavi+1    
     					selectpart(nsavi)=n           
     				endif
     			enddo
     		endif
     	enddo
  
		allocate(iondata(nsavi,16))	

	 	call mpi_allgather(nsavi,1,mpi_integer, &
               nsavilist,1,mpi_integer,mpi_comm_world,ierr)     
        call mpi_allreduce(nsavi,nsavi0,1,mpi_integer,mpi_sum &
			,mpi_comm_world,ierr)  
			
        if(nsavi0 .ne. 0) then
        	allocate(iondata0(nsavi0,16))         
        else
        	allocate(iondata0(1,16))
        endif
        
        
        do i=1, nsavi
	    	n=selectpart(i)
	    
	    	call interpolfld(real(p(n)%x),real(p(n)%y),p(n)%z,bx0,by0, &
                   bz0,ex0,ey0,ez0,dummy,dummy,dummy,dummy,dummy,dummy)
     	    
            iondata(i,1)=lap
            iondata(i,2)=p(n)%ind
            iondata(i,3)=p(n)%proc
            iondata(i,4)=p(n)%x+mxcum
            iondata(i,5)=p(n)%y+mycum
            iondata(i,6)=p(n)%z+mzcum
            iondata(i,7)=p(n)%u
            iondata(i,8)=p(n)%v
            iondata(i,9)=p(n)%w
            iondata(i,10)=sqrt(1.+(p(n)%u**2+p(n)%v**2+p(n)%w**2))
            iondata(i,11)=ex0
            iondata(i,12)=ey0
            iondata(i,13)=ez0
            iondata(i,14)=bx0
            iondata(i,15)=by0
            iondata(i,16)=bz0
            
        enddo  	           
      
        
        ! serial writing	        
        
      if(nsavi0 .ne. 0) then
      	
      	do procn=0, size0-1
           
           	if (rank .eq. 0) then
           	   
           		if (procn .ne. 0) then !receive from procn    
              
              	  	! ions
              	  	istart=sum(nsavilist(1:procn))+1
              	  	ifinal=istart+nsavilist(procn+1)-1              	  	
              	  	
              	  	call MPI_Recv(iondata0(istart:ifinal,1:16) &
                         ,nsavilist(procn+1)*16,mpi_real,procn,1, &
                         MPI_COMM_WORLD,status,ierr)                            
                 
                else  !procn eq 0              	      
               	    iondata0(1:nsavi,1:16)=iondata(1:nsavi,1:16)              	    
              	  
              	endif
              
            else !rank ne 0
            	if (procn .eq. rank) then
            		call MPI_Send(iondata(1:nsavi,1:16), &
                 	 	nsavi*16,mpi_real,0,1,MPI_COMM_WORLD,ierr)                  	            	 	
                endif
              
            endif   !rank ne 0
	
       enddo !procn          	   
       
     else	!if nsvai0 .ne. 0
     	   nsavi0=1
       	   if (rank .eq. 0) then
       	   	   iondata0(:,:)=0
       	   endif
     endif
      	   	   
          
     call mpi_barrier(MPI_COMM_WORLD,ierr)

       
       if(rank .eq. 0) then 	      	   
       	   
       	   !******************************************************
       	   !****** 	save partcle tracking
       	   !******************************************************

       	   write(fname,"(a30,i7.7)") "./output/tracking_ion/testprt.",lap
       	   print *,"",fname        
        
       	   !   Initialize FORTRAN predefined datatypes
			
       	   call h5open_f(error) 
       	   call h5fcreate_f(fname,H5F_ACC_TRUNC_F,file_id,error)
	
       	   dsetname(1)="ion"
       	   datarank=2
		
       	   data_dims(1)=nsavi0
       	   data_dims(2)=16
		
       	   call h5screate_simple_f(datarank, data_dims, dspace_id, error)
       	   call h5dcreate_f(file_id, dsetname(1), H5T_NATIVE_REAL &
				,dspace_id,dset_id,error)			
		   call h5dwrite_f(dset_id, H5T_NATIVE_REAL,iondata0, data_dims, error)
		   call h5dclose_f(dset_id,error)
		   call h5sclose_f(dspace_id,error)
   
	
		   call h5fclose_f(file_id, error)
		   call h5close_f(error)
	   
	   endif ! rank 0
	
	   if(allocated(selectpart)) deallocate(selectpart)
	   if(allocated(iondata)) deallocate(iondata)
	   if(allocated(iondata0)) deallocate(iondata0)
	   
	   call mpi_barrier(MPI_COMM_WORLD,ierr)	 	  
	      
  
end subroutine write_test_p


!-------------------------------------------------------------------------------
! 						subroutine interpolfld					 
!														
! To interpolate fields to output grid
!  							
!-------------------------------------------------------------------------------

subroutine interpolfld(x,y,z,bx0,by0,bz0,ex0,ey0,ez0,bx_ext,by_ext,bz_ext,ex_ext,ey_ext,ez_ext)
	
	implicit none
	
	! dummy variables
	
	real(sprec) :: x,y,z,bx0,by0,bz0,ex0,ey0,ez0,bx_ext,by_ext,bz_ext,ex_ext,ey_ext,ez_ext
	
	! local variables
	
	real dx, dy, dz, f, g
	real u0,v0,w0,u1,v1,w1
	integer ::l,i,j,k

	i=x
	dx=x-i
	j=y
	dy=y-j
	k=z
	dz=z-k
	if(i.eq.1)i=2
	if(i.eq.mx)i=mx-1
	if(j.eq.1)j=2
	if(j.eq.my)j=my-1
	if(k.eq.1)k=2
	if(k.eq.mz)k=mz-1
	
#ifdef twoD
		k=1
		dz=0
		iz=0
#endif
	
	l=i+iy*(j-1)+iz*(k-1)
	
	! Field interpolations are tri-linear (linear in x times linear in y
	! times linear in z). This amounts to the 3-D generalisation of "area
	! weighting". A modification of the simple linear interpolation formula
	!           f(i+dx) = f(i) + dx * (f(i+1)-f(i))
	! is needed since fields are recorded at half-integer ::locations in certain
	! dimensions: see comments and illustration with the Maxwell part of this
	! code. One then has first to interpolate from "midpoints" to "gridpoints"
	! by averaging neighbors. Then one proceeds with normal interpolation.
	! Combining these two steps leads to:
	!   f at location i+dx  = half of f(i)+f(i-1) + dx*(f(i+1)-f(i-1))
	! where now f(i) means f at location i+1/2. The halving is absorbed
	! in the final scaling.
	!   E-component interpolations:

	f=ex(l,1,1)+ex(l-ix,1,1)+dx*(ex(l+ix,1,1)-ex(l-ix,1,1))
	f=f+dy*(ex(l+iy,1,1)+ex(l-ix+iy,1,1)+dx*(ex(l+ix+iy,1,1)-ex(l-ix &
	+iy,1,1))-f)
	g=ex(l+iz,1,1)+ex(l-ix+iz,1,1)+dx*(ex(l+ix+iz,1,1)-ex(l-ix+iz,1 &
	,1))
	g=g+dy* &
	(ex(l+iy+iz,1,1)+ex(l-ix+iy+iz,1,1)+dx*(ex(l+ix+iy+iz,1,1) &
	-ex(l-ix+iy+iz,1,1))-g)
	ex0=(f+dz*(g-f))*(.5)
	
	!   ------------
	f=ey(l,1,1)+ey(l-iy,1,1)+dy*(ey(l+iy,1,1)-ey(l-iy,1,1))
	f=f+dz*(ey(l+iz,1,1)+ey(l-iy+iz,1,1)+dy*(ey(l+iy+iz,1,1)-ey(l-iy &
	+iz,1,1))-f)
	g=ey(l+ix,1,1)+ey(l-iy+ix,1,1)+dy*(ey(l+iy+ix,1,1)-ey(l-iy+ix,1 &
	,1))
	g=g+dz* &
	(ey(l+iz+ix,1,1)+ey(l-iy+iz+ix,1,1)+dy*(ey(l+iy+iz+ix,1,1) &
	-ey(l-iy+iz+ix,1,1))-g)
	ey0=(f+dx*(g-f))*(.5)
	!   ------------
	f=ez(l,1,1)+ez(l-iz,1,1)+dz*(ez(l+iz,1,1)-ez(l-iz,1,1))
	f=f+dx*(ez(l+ix,1,1)+ez(l-iz+ix,1,1)+dz*(ez(l+iz+ix,1,1)-ez(l-iz &
	+ix,1,1))-f)
	g=ez(l+iy,1,1)+ez(l-iz+iy,1,1)+dz*(ez(l+iz+iy,1,1)-ez(l-iz+iy,1 &
	,1))
	g=g+dx* &
	(ez(l+ix+iy,1,1)+ez(l-iz+ix+iy,1,1)+dz*(ez(l+iz+ix+iy,1,1) &
	-ez(l-iz+ix+iy,1,1))-g)
	ez0=(f+dy*(g-f))*(.5)
	!   ---------
	!   B-component interpolations:
	f=bx(l-iy,1,1)+bx(l-iy-iz,1,1)+dz*(bx(l-iy+iz,1,1)-bx(l-iy-iz,1 &
	,1))
	f=bx(l,1,1)+bx(l-iz,1,1)+dz*(bx(l+iz,1,1)-bx(l-iz,1,1))+f+dy &
	* (bx(l+iy,1,1)+bx(l+iy-iz,1,1)+dz*(bx(l+iy+iz,1,1)-bx(l+iy &
	-iz,1,1))-f)
	g=bx(l+ix-iy,1,1)+bx(l+ix-iy-iz,1,1)+dz*(bx(l+ix-iy+iz,1,1) &
	-bx(l+ix-iy-iz,1,1))
	g=bx(l+ix,1,1)+bx(l+ix-iz,1,1)+dz*(bx(l+ix+iz,1,1)-bx(l+ix-iz,1 &
	,1))+g+dy*(bx(l+ix+iy,1,1)+bx(l+ix+iy-iz,1,1)+dz*(bx(l+ix &
	+iy+iz,1,1)-bx(l+ix+iy-iz,1,1))-g)
	bx0=(f+dx*(g-f))*(.25)
	!   ------------
	f=by(l-iz,1,1)+by(l-iz-ix,1,1)+dx*(by(l-iz+ix,1,1)-by(l-iz-ix,1 &
	,1))
	f=by(l,1,1)+by(l-ix,1,1)+dx*(by(l+ix,1,1)-by(l-ix,1,1))+f+dz &
	* (by(l+iz,1,1)+by(l+iz-ix,1,1)+dx*(by(l+iz+ix,1,1)-by(l+iz &
	-ix,1,1))-f)
	g=by(l+iy-iz,1,1)+by(l+iy-iz-ix,1,1)+dx*(by(l+iy-iz+ix,1,1)-by(l &
	+iy-iz-ix,1,1))
	g=by(l+iy,1,1)+by(l+iy-ix,1,1)+dx*(by(l+iy+ix,1,1)-by(l+iy-ix,1 &
	,1))+g+dz*(by(l+iy+iz,1,1)+by(l+iy+iz-ix,1,1)+dx*(by(l+iy &
	+iz+ix,1,1)-by(l+iy+iz-ix,1,1))-g)
	by0=(f+dy*(g-f))*(.25)
	!   ------------
	f=bz(l-ix,1,1)+bz(l-ix-iy,1,1)+dy*(bz(l-ix+iy,1,1)-bz(l-ix-iy,1 &
	,1))
	f=bz(l,1,1)+bz(l-iy,1,1)+dy*(bz(l+iy,1,1)-bz(l-iy,1,1))+f+dx &
	* (bz(l+ix,1,1)+bz(l+ix-iy,1,1)+dy*(bz(l+ix+iy,1,1)-bz(l+ix &
	-iy,1,1))-f)
	g=bz(l+iz-ix,1,1)+bz(l+iz-ix-iy,1,1)+dy*(bz(l+iz-ix+iy,1,1)-bz(l &
	+iz-ix-iy,1,1))
	g=bz(l+iz,1,1)+bz(l+iz-iy,1,1)+dy*(bz(l+iz+iy,1,1)-bz(l+iz-iy,1 &
	,1))+g+dx*(bz(l+iz+ix,1,1)+bz(l+iz+ix-iy,1,1)+dy*(bz(l+iz &
	+ix+iy,1,1)-bz(l+iz+ix-iy,1,1))-g)
	bz0=(f+dz*(g-f))*(.25)
	
	!        call deutschsub((x),(y),(z),
	!     &   bxd,byd,bzd,exd,eyd,ezd)
	
	bx_ext=0
	by_ext=0
	bz_ext=0
	ex_ext=0
	ey_ext=0
	ez_ext=0

	if(external_fields) then
	   call get_external_fields(x,y,z,ex_ext,ey_ext,ez_ext,bx_ext,by_ext,bz_ext)
	endif


end subroutine interpolfld



!-------------------------------------------------------------------------------
! 						subroutine interpolcurrent					 
!														
! To interpolate currents to output grid 
! Not used 							
!
!-------------------------------------------------------------------------------

subroutine interpolcurrent(x,y,z,curx0, cury0, curz0)

	implicit none
	
	! dummy variables
	
	real(sprec) :: x, y, z, curx0, cury0, curz0
	
	! local variables
	
	real dx, dy, dz, f, g
	integer ::l,i,j,k
	
	
	i=x
	dx=x-i
	j=y
	dy=y-j
	k=z
	dz=z-k
	
	if(i.eq.1)i=2
	if(i.eq.mx)i=mx-1
	if(j.eq.1)j=2
	if(j.eq.my)j=my-1
	if(k.eq.1)k=2
	if(k.eq.mz)k=mz-1
	
#ifdef twoD
		k=1
		dz=0
		iz=0
#endif
	
	l=i+iy*(j-1)+iz*(k-1)
	
	! Field interpolations are tri-linear (linear in x times linear in y
	! times linear in z). This amounts to the 3-D generalisation of "area
	! weighting". A modification of the simple linear interpolation formula
	!           f(i+dx) = f(i) + dx * (f(i+1)-f(i))
	! is needed since fields are recorded at half-integer ::locations in certain
	! dimensions: see comments and illustration with the Maxwell part of this
	! code. One then has first to interpolate from "midpoints" to "gridpoints"
	! by averaging neighbors. Then one proceeds with normal interpolation.
	! Combining these two steps leads to:
	!   f at location i+dx  = half of f(i)+f(i-1) + dx*(f(i+1)-f(i-1))
	! where now f(i) means f at location i+1/2. The halving is absorbed
	! in the final scaling.
	!   E-component interpolations:
	
	f=curx(l,1,1)+curx(l-ix,1,1)+dx*(curx(l+ix,1,1)-curx(l-ix,1,1))
	f=f+dy*(curx(l+iy,1,1)+curx(l-ix+iy,1,1)+dx*(curx(l+ix+iy,1,1) &
	-curx(l-ix+iy,1,1))-f)
	g=curx(l+iz,1,1)+curx(l-ix+iz,1,1)+dx*(curx(l+ix+iz,1,1)-curx(l &
	-ix+iz,1,1))
	g=g+dy* (curx(l+iy+iz,1,1)+curx(l-ix+iy+iz,1,1)+dx*(curx(l+ix+iy &
	+iz,1,1) -curx(l-ix+iy+iz,1,1))-g)
	curx0=(f+dz*(g-f))*(.5)
	
	!   ------------
	f=cury(l,1,1)+cury(l-iy,1,1)+dy*(cury(l+iy,1,1)-cury(l-iy,1,1))
	f=f+dz*(cury(l+iz,1,1)+cury(l-iy+iz,1,1)+dy*(cury(l+iy+iz,1,1) &
	-cury(l-iy+iz,1,1))-f)
	g=cury(l+ix,1,1)+cury(l-iy+ix,1,1)+dy*(cury(l+iy+ix,1,1)-cury(l &
	-iy+ix,1,1))
	g=g+dz* (cury(l+iz+ix,1,1)+cury(l-iy+iz+ix,1,1)+dy*(cury(l+iy+iz &
	+ix,1,1) -cury(l-iy+iz+ix,1,1))-g)
	cury0=(f+dx*(g-f))*(.5)
	!   ------------
	f=curz(l,1,1)+curz(l-iz,1,1)+dz*(curz(l+iz,1,1)-curz(l-iz,1,1))
	f=f+dx*(curz(l+ix,1,1)+curz(l-iz+ix,1,1)+dz*(curz(l+iz+ix,1,1) &
	-curz(l-iz+ix,1,1))-f)
	g=curz(l+iy,1,1)+curz(l-iz+iy,1,1)+dz*(curz(l+iz+iy,1,1)-curz(l &
	-iz+iy,1,1))
	g=g+dx* (curz(l+ix+iy,1,1)+curz(l-iz+ix+iy,1,1)+dz*(curz(l+iz+ix &
	+iy,1,1) -curz(l-iz+ix+iy,1,1))-g)
	curz0=(f+dy*(g-f))*(.5)
	!   ---------

end subroutine interpolcurrent


!-------------------------------------------------------------------------------
! 						subroutine meanq_fld_cur			 
!											
! To compute derived fluid quantities, like density and density of ions	
! or 3 velocity of the fluid, or 3 velocity of ions
!  							
!-------------------------------------------------------------------------------

subroutine meanq_fld_cur(totname)

  implicit none
  
  ! local variables
  
  !save only density array
  character (len=5) totname
  real uloc, vloc, wloc, gamprt, addprtx, addprty
  integer ::count, uprank, dwnrank,nn0, i, j, k, lx1, lx2, ly1, ly2, &
       lz1, lz2, uptag, dwntag, comm,ierr,status(statsize)
  integer ::request(2), status1(statsize,2)
  real, allocatable :: buffer(:,:,:), bufferout(:,:,:)
  integer ::lfttag, rgttag, lftrank, rgtrank


  allocate(buffer(mx,my,2), bufferout(mx,my,2))
  
  if(debug .and. rank.eq.0) print *, rank,": ","in meanq", mx,my,mz
  
  uptag=100
  dwntag=200
  rgttag=uptag
  lfttag=dwntag
  
  comm=MPI_Comm_world
  
  maxhlf=maxptl/2
  curx=0.
  cury=0.
  curz=0.
  
  do nn0=1,maxptl
     if(nn0.gt.ions .and. nn0 .lt. maxhlf) cycle
     if(nn0 .gt. maxhlf+lecs) cycle
    gamprt=1./sqrt(1.+p(nn0)%u**2+p(nn0)%v**2+p(nn0)%w**2)

     addprtx=0.
     addprty=0.
     select case(totname)
        !density
     case('tdens')
        addprtx=p(nn0)%ch
     case('idens')
        if(nn0 .le. ions) then
           addprtx=p(nn0)%ch
        endif
     case('hdens')
     if(nn0 .ge. maxhlf+1 .and. p(nn0)%ind .gt. 0) then
        addprtx=p(nn0)%ch   
     endif
     case('ldens')
     if(nn0 .ge. maxhlf+1 .and. p(nn0)%ind .lt. 0) then
        addprtx=1.   
     endif
        !beam density
     case('btden')
        if(p(nn0)%ind .lt. 0) then
           addprtx=p(nn0)%ch
        endif
     case('biden')
        if(nn0 .le. ions .and. p(nn0)%ind .lt. 0) then
           addprtx=p(nn0)%ch
        endif
        !3-velocity
     case('tbetx')
        addprtx=p(nn0)%u*gamprt*p(nn0)%ch
        addprty=p(nn0)%ch
     case('ebetx')
     	if(nn0 .ge. maxhlf+1) then
     		addprtx=p(nn0)%u*gamprt*p(nn0)%ch
     		addprty=p(nn0)%ch
		endif        
     case('ibetx')
        if(nn0 .le. ions) then
           addprtx=p(nn0)%u*gamprt*p(nn0)%ch
           addprty=p(nn0)%ch        
        endif
     case('tbety')
        addprtx=p(nn0)%v*gamprt*p(nn0)%ch
        addprty=p(nn0)%ch
     case('ebety')
        if(nn0 .ge. maxhlf+1) then
           addprtx=p(nn0)%v*gamprt*p(nn0)%ch
           addprty=p(nn0)%ch        
        endif        
     case('ibety')
        if(nn0 .le. ions) then
           addprtx=p(nn0)%v*gamprt*p(nn0)%ch
           addprty=p(nn0)%ch        
        endif
     case('tbetz')
        addprtx=p(nn0)%w*gamprt*p(nn0)%ch
        addprty=p(nn0)%ch
     case('ebetz')
        if(nn0 .ge. maxhlf+1) then
           addprtx=p(nn0)%w*gamprt*p(nn0)%ch
           addprty=p(nn0)%ch        
        endif
     case('ibetz')
        if(nn0 .le. ions) then
           addprtx=p(nn0)%w*gamprt*p(nn0)%ch
           addprty=p(nn0)%ch        
        endif
        !4-velocity
     case('tmomx')
        addprtx=p(nn0)%u*p(nn0)%ch
        addprty=p(nn0)%ch
     case('imomx')
        if(nn0 .le. ions) then
           addprtx=p(nn0)%u*p(nn0)%ch
           addprty=p(nn0)%ch        
        endif
     case('tmomy')
        addprtx=p(nn0)%v*p(nn0)%ch
        addprty=p(nn0)%ch
     case('imomy')
        if(nn0 .le. ions) then
           addprtx=p(nn0)%v*p(nn0)%ch
           addprty=p(nn0)%ch        
        endif
     case('tmomz')
        addprtx=p(nn0)%w*p(nn0)%ch
        addprty=p(nn0)%ch
     case('imomz')
        if(nn0 .le. ions) then
           addprtx=p(nn0)%w*p(nn0)%ch
           addprty=p(nn0)%ch        
        endif
        !energy
     case('eener')
    	if(nn0 .ge. maxhlf+1) then
        addprtx=(1./gamprt-1.)*p(nn0)%ch
        addprty=p(nn0)%ch
        endif
     case('iener')
        if(nn0 .le. ions) then
           addprtx=(1./gamprt-1.)*p(nn0)%ch
           addprty=p(nn0)%ch        
        endif   
        !3 velocity square
     case('eetx2')
     	if(nn0 .ge. maxhlf+1) then
     		addprtx=(p(nn0)%u*gamprt)**2*p(nn0)%ch
     		addprty=p(nn0)%ch
		endif        
     case('ietx2')
        if(nn0 .le. ions) then
           addprtx=(p(nn0)%u*gamprt)**2*p(nn0)%ch
           addprty=p(nn0)%ch        
        endif
     case('eety2')
        if(nn0 .ge. maxhlf+1) then
           addprtx=(p(nn0)%v*gamprt)**2*p(nn0)%ch
           addprty=p(nn0)%ch        
        endif        
     case('iety2')
        if(nn0 .le. ions) then
           addprtx=(p(nn0)%v*gamprt)**2*p(nn0)%ch
           addprty=p(nn0)%ch        
        endif
     case('eetz2')
        if(nn0 .ge. maxhlf+1) then
           addprtx=(p(nn0)%w*gamprt)**2*p(nn0)%ch
           addprty=p(nn0)%ch        
        endif
     case('ietz2')
        if(nn0 .le. ions) then
           addprtx=(p(nn0)%w*gamprt)**2*p(nn0)%ch
           addprty=p(nn0)%ch        
        endif        
        
     end select

     if(nn0.le.ions .or. (nn0.gt.maxhlf .and. nn0.le.maxhlf+lecs)) then
        if(debug .and. ions .eq. 0) print *,rank,": ions=0",ions, lecs, nn0              
        
        i=p(nn0)%x
        j=p(nn0)%y
        k=p(nn0)%z
        
        lz1=max(k-idz,1)
        lz2=min(k+idz,mz)
        
#ifdef twoD
        lz1=1
        lz2=1
#endif
        
        ly1=max(j-idy,1)
        ly2=min(j+idy,my)
        
        lx1=max(i-idx,1)
        lx2=min(i+idx,mx)
        
        do k=lz1,lz2
           do j=ly1,ly2
              do i=lx1,lx2
              
                 curx(i,j,k)=curx(i,j,k)+addprtx
                 cury(i,j,k)=cury(i,j,k)+addprty

              enddo
           enddo
        enddo
        
     endif
  enddo                     !particle loop
  
  
	call exchange_current()
  !:::::::::::::::::::::::::::::::::::::::::::::::::::::::::
  !normalize density
  do k=1,mz
     do j=1,my
        do i=1,mx
           
           lz1=max(k-idz,1)
           lz2=min(k+idz,mz)
           
           ly1=max(j-idy,1)
           ly2=min(j+idy,my)
           
           lx1=max(i-idx,1)
           lx2=min(i+idx,mx)
           
#ifdef twoD
           lz1=1
           lz2=1
#endif
           curx(i,j,k)=curx(i,j,k)/((lx2-lx1+1)*(ly2-ly1+1)*(lz2-lz1+1))
           cury(i,j,k)=cury(i,j,k)/((lx2-lx1+1)*(ly2-ly1+1)*(lz2-lz1+1))
         
!           curx(i,j,k)=curx(i,j,k)/(idx*idy*idz)
!           cury(i,j,k)=cury(i,j,k)/(idx*idy*idz)
           
 
        enddo
     enddo
  enddo

  !to compute average 3-velocity, 4-velocity and energy
  if (totname .ne. 'tdens' .and. totname .ne. 'idens' .and. &
  	  totname .ne. 'hdens' .and. totname .ne. 'ldens') then 
     where(cury .ne. 0.) 
        curx=curx/cury
     elsewhere 
        curx=0.
     endwhere
  endif

  deallocate(buffer,bufferout)
  
end subroutine meanq_fld_cur


#ifdef twoD
end module m_output
#else
end module m_output_3d
#endif
