!
! domain module
!
! This module contains routines that alter the simulation domain, and 
! the routine that enlarges the simulation domain
!
!

#ifdef twoD 

module m_domain

	use m_system
	use m_aux
	use m_communications
	use m_fields
	use m_fieldboundaries
	use m_particles
	use m_globaldata
	use m_inputparser
	
#else

module m_domain_3d

	use m_system_3d
	use m_aux_3d
	use m_communications_3d
	use m_fields_3d
	use m_fieldboundaries_3d
	use m_particles_3d
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

	integer :: enlarge, shiftint, shiftstart
	character (len=5) :: rankchar
	character (len=34) :: fenlargeloc
	
	
!-------------------------------------------------------------------------------
!	INTERFACE DECLARATIONS
!-------------------------------------------------------------------------------
	
!-------------------------------------------------------------------------------
!	PUBLIC MODIFIERS
!-------------------------------------------------------------------------------

	! public functions
	
	public ::  enlarge_domain, enlarge_ydomain,  read_input_domain

	! public variables (used in outputs)

	public :: shiftint, shiftstart, rankchar
	
	! public variables (used in initialize)
	
	public :: fenlargeloc, enlarge
	

!-------------------------------------------------------------------------------
!	MODULE PROCEDURES AND FUNCTIONS
!-------------------------------------------------------------------------------

	contains



!-------------------------------------------------------------------------------
! 						subroutine read_input_domain					 
!																		
! Reads any variables related to (or needed by) this module
! 							
!-------------------------------------------------------------------------------

subroutine read_input_domain()

	implicit none

	call inputpar_geti_def("domain", "enlarge", 0, enlarge)	
	
end subroutine read_input_domain


!-------------------------------------------------------------------------------
! 						subroutine initialize_domain_default()					 
!																		
! Sets enlarge to 0 for the default problem
!							
!-------------------------------------------------------------------------------

!subroutine initialize_domain_default()
!
!	implicit none
!	
!	enlarge=0

!end subroutine initialize_domain_default



!-------------------------------------------------------------------------------
! 						subroutine enlarge_domain()					 
!																		
! Calls the enlarge_domain_with_disk subroutine, if needed
!							
!-------------------------------------------------------------------------------

subroutine enlarge_domain()

	implicit none
	
	if(enlarge .ne. 1) return 
	
#ifndef twoD 
	if(xinject2 .le. mx0-300 .or. periodicx .ne. 0) return 
#else
!	if(xinject2 .le. mx0-500 .or. periodicx .ne. 0) return 
	if(xinject2 .le. mx0-200 .or. periodicx .ne. 0) return 
#endif
	
	if(debug) print *, rank, ": b4 enlarge"
	call enlarge_domain_with_disk()	


end subroutine enlarge_domain



!-------------------------------------------------------------------------------
! 						subroutine enlarge_domain()					 
!																		
! Calls the enlarge_domain_with_disk subroutine, if needed
!							
!-------------------------------------------------------------------------------

subroutine enlarge_ydomain()

	implicit none
	
	if(enlarge .ne. 1) return 
	
#ifndef twoD 
	if(yinject2 .le. my0-20 .or. periodicy .ne. 0) return 
#else
!	if(yinject2 .le. my0-500 .or. periodicy .ne. 0) return 
	if(yinject2 .le. my0-20 .or. periodicy .ne. 0) return 
#endif
	
	if(debug) print *, rank, ": b4 enlarge"
	call enlarge_ydomain_with_disk()	


end subroutine enlarge_ydomain



!-------------------------------------------------------------------------------
! 						subroutine enlarge_domain_with_disk()					 
!													
! Enlarges the simulation domain as needed
!							
!-------------------------------------------------------------------------------

subroutine enlarge_domain_with_disk()

	implicit none
	external m_overload_mp_init_EMfields
	
	! local variables

	integer :: mxold, idummy, i , j, k, ierr, i1, i2
	character*20 testout
	
	! enlarge the size of the right end processors
	if(modulo(rank,sizex) .eq. sizex-1 ) then	
		mxold=mx	
		mx=mxold+300
		!mx=mxold+100

!	if(rank .eq. 0) print *, "outin-s done"
		iy=mx
		iz=iy*my
		lot=iz*mz
#ifdef twoD
		iz=0
		lot=mx*my
#endif
	
#ifndef twoD
!		buffsize = max(10*int(ppc0*c*max((mx-5)*(my-5),1*(mx-5)*(mz-5))) &
!				,10000)
		buffsize = max(10*int(ppc0*upsamp_e*c*max((mx-5)*(my-5),1*(mx-5)*(mz-5))) &
				,10000)	
				
#else
!		buffsize = max(3*int(ppc0*c*(1*(mx-5))),60000)
		buffsize = max(3*int(ppc0*upsamp_e*c*max((mx-5),(my-5))),60000)
		
!		if (splitparts) buffsize=buffsize*150 
#endif

		testout="err"//"."//trim(rankchar)
	
		print*, rank, "Enlarging domain!", mxold, mx
	
		!first save the arrays to disk as if I am doing restart.
!		call mpi_barrier(MPI_COMM_WORLD,ierr)
	
		open(unit=7,file=fenlargeloc, form='unformatted')
		rewind 7
		print *, rank,":writing enlarge file"
		write(7)(((bx(i,j,k),i=1,mxold-3),j=1,my),k=1,mz), &
		(((by(i,j,k),i=1,mxold-3),j=1,my),k=1,mz), &
		(((bz(i,j,k),i=1,mxold-3),j=1,my),k=1,mz), &
		(((ex(i,j,k),i=1,mxold-3),j=1,my),k=1,mz), &
		(((ey(i,j,k),i=1,mxold-3),j=1,my),k=1,mz), &
		(((ez(i,j,k),i=1,mxold-3),j=1,my),k=1,mz), &
		(((curx(i,j,k),i=1,mxold-3),j=1,my),k=1,mz), &
		(((cury(i,j,k),i=1,mxold-3),j=1,my),k=1,mz), &
		(((curz(i,j,k),i=1,mxold-3),j=1,my),k=1,mz)
		close(7)

		deallocate(bx,by,bz,ex,ey,ez,curx,cury,curz)

		deallocate(bufferin1,bufferin2)

		deallocate(bufferin1y,bufferin2y,bufferin1x,bufferin2x)

		deallocate(sendbufy,sendbufz)

		deallocate( temp, poutup, poutdwn, pinblw &
			, pinabv, poutlft, poutrgt,pinlft, pinrgt &
			, poutplus, poutminus,pinplus, pinminus, pall )

#ifdef filter2
		deallocate(xghost, yghost, zghost) !for filter2
! ! 	deallocate(temp_x, temp_y, temp_z) ! ghost|cur|ghost for filter2copy
	
#ifndef twoD
	call free_MPI_filter_datatypes()
	call create_MPI_filter_datatypes(3,mx-3,3,my-3,3,mz-3)
#endif
#endif!ifdef filter2

		allocate(bx(mx,my,mz),by(mx,my,mz),bz(mx,my,mz),ex(mx,my,mz),ey(mx &
			,my,mz),ez(mx,my,mz),curx(mx,my,mz),cury(mx,my,mz),curz(mx,my &
			,mz),bufferin1(mx,my,2),bufferin2(mx,my,2), &
			bufferin1y(mx,2,mz),bufferin2y(mx,2,mz))

	
		allocate(bufferin1x(2,my,mz),bufferin2x(2,my,mz))
		allocate(sendbufy(mx,2,mz),sendbufz(mx,my,2))
		allocate(temp(mx,my,mz))
#ifdef filter2
		allocate(xghost(2*ntimes,my,mz), yghost(mx,2*ntimes,mz), zghost(mx,my,2*ntimes)) !for filter2
#endif	
! ! 	allocate(temp_x(ntimes+(mx-5)+ntimes), temp_y(ntimes+(my-5)+ntimes), temp_z(ntimes+(mz-5)+ntimes)) ! ghost|cur|ghost for filter2copy

		allocate(poutup(buffsize),poutdwn(buffsize))
		allocate(pinabv(buffsize),pinblw(buffsize))
		allocate(poutrgt(buffsize),poutlft(buffsize))
		allocate(pinlft(buffsize),pinrgt(buffsize))
		allocate(poutplus(buffsize),poutminus(buffsize))
		allocate(pinminus(buffsize),pinplus(buffsize))
		allocate(pall(lot))
		
		pall=0
		curx=0.
		cury=0.
		curz=0.

!this initialization is specific to shocks

		bx=binit*cos(btheta)
		by=binit*sin(btheta)*sin(bphi)
		bz=binit*sin(btheta)*cos(bphi)				
		ex=0.			
		ey=(-beta)*bz
		ez=-(-beta)*by
			
!	do  k=1,mz
!		do  j=1,my
!			do  i=1,mx
!
!				bx(i,j,k)=Binit*cos(btheta) 
!				by(i,j,k)=Binit*sin(btheta)*sin(bphi)
!				bz(i,j,k)=Binit*sin(btheta)*cos(bphi)
!
!				ex(i,j,k)=0.
!				ey(i,j,k)=(-beta)*bz(i,j,k) 
!				ez(i,j,k)=-(-beta)*by(i,j,k)
!
!			enddo
!		enddo
!	enddo

		
		open(unit=7,file=fenlargeloc, form='unformatted')
		rewind 7
		print *, rank,":reading back enlarge file"
	
		read(7)(((bx(i,j,k),i=1,mxold-3),j=1,my),k=1,mz), &
			(((by(i,j,k),i=1,mxold-3),j=1,my),k=1,mz), &
			(((bz(i,j,k),i=1,mxold-3),j=1,my),k=1,mz), &
			(((ex(i,j,k),i=1,mxold-3),j=1,my),k=1,mz), &
			(((ey(i,j,k),i=1,mxold-3),j=1,my),k=1,mz), &
			(((ez(i,j,k),i=1,mxold-3),j=1,my),k=1,mz), &
			(((curx(i,j,k),i=1,mxold-3),j=1,my),k=1,mz), &
			(((cury(i,j,k),i=1,mxold-3),j=1,my),k=1,mz), &
			(((curz(i,j,k),i=1,mxold-3),j=1,my),k=1,mz)     
	
		close(7,status='delete')
		bufferin1=0.
		bufferin2=0.
		bufferin1y=0.
		bufferin2y=0.
		bufferin1x=0.
		bufferin2x=0.
		sendbufy=0.
		sendbufz=0.
		temp=0.

		
		endif	!if(modulo(rank,sizex) .eq. sizex-1 )

		call mpi_barrier(MPI_COMM_WORLD,ierr)
	
		call mpi_allgather(mx,1,mpi_integer,mxl,1,mpi_integer &
		  	   	  ,mpi_comm_world, ierr)

		i1=1
		i2=sizex
		mx0=sum(mxl(i1:i2)-5)+5
		
		i1=1
		i2=modulo(rank,sizex)+1
		mxcum=sum(mxl(i1:i2)-5)-(mxl(rank+1)-5)
	
		if(rank .eq. 0) print *, "enlarge done"
	
		
end subroutine enlarge_domain_with_disk


!-------------------------------------------------------------------------------
! 						subroutine enlarge_ydomain_with_disk()					 
!													
! Enlarges the simulation domain as needed
!							
!-------------------------------------------------------------------------------

subroutine enlarge_ydomain_with_disk()

	implicit none
	external m_overload_mp_init_EMfields
	
	! local variables

	integer :: myold, idummy, i , j, k, ierr, iglob
	integer :: j1, j2
	character*20 testout
	
	! enlarge the size of the right end processors
	if(modulo(rank,sizex*sizey)/sizex .eq. sizey-1 ) then	
		myold=my	
		my=myold+100
	
!	if(rank .eq. 0) print *, "outin-s done"
		iy=mx
		iz=iy*my
		lot=iz*mz
#ifdef twoD
		iz=0
		lot=mx*my
#endif
	
#ifndef twoD
!		buffsize = max(1*int(ppc0*c*max((mx-5)*(my-5),1*(mx-5)*(mz-5))) &
!				,10000)
		
		buffsize = max(10*int(ppc0*c*max((mx-5)*(my-5),1*(mx-5)*(mz-5))) &
				,10000)	
		
#else
!		buffsize = max(3*int(ppc0*c*(1*(mx-5))),60000)
		buffsize = max(3*int(ppc0*c*max((mx-5),(my-5))),60000)
		
!		if (splitparts) buffsize=buffsize*10 
#endif

		testout="err"//"."//trim(rankchar)
	
		print*, rank, "Enlarging y domain!", myold, my
	
		!first save the arrays to disk as if I am doing restart.
!		call mpi_barrier(MPI_COMM_WORLD,ierr)
	
		open(unit=7,file=fenlargeloc, form='unformatted')
		rewind 7
		print *, rank,":writing enlarge file"
		write(7)(((bx(i,j,k),i=1,mx),j=1,myold-3),k=1,mz), &
		(((by(i,j,k),i=1,mx),j=1,myold-3),k=1,mz), &
		(((bz(i,j,k),i=1,mx),j=1,myold-3),k=1,mz), &
		(((ex(i,j,k),i=1,mx),j=1,myold-3),k=1,mz), &
		(((ey(i,j,k),i=1,mx),j=1,myold-3),k=1,mz), &
		(((ez(i,j,k),i=1,mx),j=1,myold-3),k=1,mz), &
		(((curx(i,j,k),i=1,mx),j=1,myold-3),k=1,mz), &
		(((cury(i,j,k),i=1,mx),j=1,myold-3),k=1,mz), &
		(((curz(i,j,k),i=1,mx),j=1,myold-3),k=1,mz)
		close(7)
	
!		call mpi_barrier(MPI_COMM_WORLD,ierr)
	
		deallocate(bx,by,bz,ex,ey,ez,curx,cury,curz)

		deallocate(bufferin1,bufferin2)

		deallocate(bufferin1y,bufferin2y,bufferin1x,bufferin2x)

		deallocate(sendbufy,sendbufz)

		deallocate( temp, poutup, poutdwn, pinblw &
			, pinabv, poutlft, poutrgt,pinlft, pinrgt &
			, poutplus, poutminus,pinplus, pinminus, pall )
#ifdef filter2
		deallocate(xghost, yghost, zghost) !for filter2
! ! 	deallocate(temp_x, temp_y, temp_z) ! ghost|cur|ghost for filter2copy
	
!	call mpi_barrier(MPI_COMM_WORLD,ierr)

#ifndef twoD
	call free_MPI_filter_datatypes()
	call create_MPI_filter_datatypes(3,mx-3,3,my-3,3,mz-3)
#endif
#endif !ifdef filter2
		allocate(bx(mx,my,mz),by(mx,my,mz),bz(mx,my,mz),ex(mx,my,mz),ey(mx &
			,my,mz),ez(mx,my,mz),curx(mx,my,mz),cury(mx,my,mz),curz(mx,my &
			,mz),bufferin1(mx,my,2),bufferin2(mx,my,2), &
			bufferin1y(mx,2,mz),bufferin2y(mx,2,mz))
	
!	call mpi_barrier(MPI_COMM_WORLD,ierr)
	
		allocate(bufferin1x(2,my,mz),bufferin2x(2,my,mz))
		allocate(sendbufy(mx,2,mz),sendbufz(mx,my,2))
		allocate(temp(mx,my,mz))
#ifdef filter2
		allocate(xghost(2*ntimes,my,mz), yghost(mx,2*ntimes,mz), zghost(mx,my,2*ntimes)) !for filter2
#endif	
! ! 	allocate(temp_x(ntimes+(mx-5)+ntimes), temp_y(ntimes+(my-5)+ntimes), temp_z(ntimes+(mz-5)+ntimes)) ! ghost|cur|ghost for filter2copy

		allocate(poutup(buffsize),poutdwn(buffsize))
		allocate(pinabv(buffsize),pinblw(buffsize))
		allocate(poutrgt(buffsize),poutlft(buffsize))
		allocate(pinlft(buffsize),pinrgt(buffsize))
		allocate(poutplus(buffsize),poutminus(buffsize))
		allocate(pinminus(buffsize),pinplus(buffsize))
		allocate(pall(lot))

		pall=0
		curx=0.
		cury=0.
		curz=0.
	
	do  k=1,mz
		do  j=1,my
			do  i=1,mx
				iglob=i+mxcum			
				bx(i,j,k)=0.
				by(i,j,k)=0.
				bz(i,j,k)=0.

				ex(i,j,k)=0.
				ey(i,j,k)=0.
				ez(i,j,k)=0.
			enddo
		enddo
	enddo
	
		open(unit=7,file=fenlargeloc, form='unformatted')
		rewind 7
		print *, rank,":reading back y enlarge file"
	
		read(7)(((bx(i,j,k),i=1,mx),j=1,myold-3),k=1,mz), &
			(((by(i,j,k),i=1,mx),j=1,myold-3),k=1,mz), &
			(((bz(i,j,k),i=1,mx),j=1,myold-3),k=1,mz), &
			(((ex(i,j,k),i=1,mx),j=1,myold-3),k=1,mz), &
			(((ey(i,j,k),i=1,mx),j=1,myold-3),k=1,mz), &
			(((ez(i,j,k),i=1,mx),j=1,myold-3),k=1,mz), &
			(((curx(i,j,k),i=1,mx),j=1,myold-3),k=1,mz), &
			(((cury(i,j,k),i=1,mx),j=1,myold-3),k=1,mz), &
			(((curz(i,j,k),i=1,mx),j=1,myold-3),k=1,mz)     
	
		close(7,status='delete')
		bufferin1=0.
		bufferin2=0.
		bufferin1y=0.
		bufferin2y=0.
		bufferin1x=0.
		bufferin2x=0.
		sendbufy=0.
		sendbufz=0.
		temp=0.
	
		endif	!if(modulo(rank,sizex) .eq. sizex-1 )

		call mpi_barrier(MPI_COMM_WORLD,ierr)
	
		call mpi_allgather(my,1,mpi_integer,myl,1,mpi_integer &
		  	   	  ,mpi_comm_world, ierr)

		j1=1
		j2=sizex*(sizey-1)+1
		my0=sum(myl(j1:j2:sizex)-5)+5

		j1=1
		j2=modulo(rank,sizex*sizey)/sizex*sizex+1
		mycum=sum(myl(j1:j2:sizex)-5)-(myl(j2)-5)
		
		if(rank .eq. 0) print *, "enlarge y done"
	
	
		call mpi_barrier(MPI_COMM_WORLD,ierr)
		
end subroutine enlarge_ydomain_with_disk



#ifdef twoD
end module m_domain
#else
end module m_domain_3d
#endif
