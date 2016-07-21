
!
! User module
!
! This module contains functions that may be altered by a user of the code, and that are called 
! if caseinit variable is set to a number greater than 0. The functions that are going to be 
! called in such a case are: SetEMFieldsUser, ..., ...
!
! If a user wants to alter the functions that are called he/she may also alter the module m_overload
! which branches with the variable caseinit.


#ifdef twoD 

module m_dynamic_domain

	use m_globaldata
	use m_system
	use m_aux
	use m_communications
	use m_fields
	use m_fieldboundaries
	use m_particles
	use m_inputparser
	use m_fparser
	use m_domain
	
#else

module m_dynamic_domain_3d

	use m_globaldata_3d
	use m_system_3d
	use m_aux_3d
	use m_communications_3d
	use m_fields_3d
	use m_fieldboundaries_3d
	use m_particles_3d
	use m_inputparser_3d 
	use m_fparser_3d
	use m_domain_3d
	
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

	integer :: dynmx, dynmxstart, dynmxinterval, mxmin, &
			   dynmy, dynmystart, dynmyinterval, mymin, &
			   dynmz, dynmzstart, dynmzinterval, mzmin
			
		
!-------------------------------------------------------------------------------
!	INTERFACE DECLARATIONS
!-------------------------------------------------------------------------------
	
!-------------------------------------------------------------------------------
!	PUBLIC MODIFIERS
!-------------------------------------------------------------------------------

	public :: read_input_dynamic_grid
	public :: redist_x_domain, redist_y_domain, redist_z_domain
	
	public :: dynmx, dynmxstart, dynmxinterval, mxmin, &
			  dynmy, dynmystart, dynmyinterval, mymin, &
			  dynmz, dynmzstart, dynmzinterval, mzmin
!-------------------------------------------------------------------------------
!	MODULE PROCEDURES AND FUNCTIONS
!-------------------------------------------------------------------------------

	contains


!-------------------------------------------------------------------------------
! 						subroutine read_input_grid					 
!																		
! Reads any variables related to (or needed by) this module
!							
!-------------------------------------------------------------------------------

subroutine read_input_dynamic_grid()

	implicit none

	! local variables

	integer :: ldynmx, ldynmxstart, ldynmxinterval, lmxmin
	integer :: ldynmy, ldynmystart, ldynmyinterval, lmymin
	integer :: ldynmz, ldynmzstart, ldynmzinterval, lmzmin
	
	call inputpar_geti_def("dynamic_grid", "dynmx", 0, ldynmx)
	call inputpar_geti_def("dynamic_grid", "dynmxstart", 300, ldynmxstart)
	call inputpar_geti_def("dynamic_grid", "dynmxinterval", 400, ldynmxinterval)
	call inputpar_geti_def("dynamic_grid", "mxmin", 7, lmxmin)
	

	call inputpar_geti_def("dynamic_grid", "dynmy", 0, ldynmy)
	call inputpar_geti_def("dynamic_grid", "dynmystart", 200, ldynmystart)
	call inputpar_geti_def("dynamic_grid", "dynmyinterval", 400, ldynmyinterval)
	call inputpar_geti_def("dynamic_grid", "mymin", 7, lmymin)
	
	call inputpar_geti_def("dynamic_grid", "dynmz", 0, ldynmz)
	call inputpar_geti_def("dynamic_grid", "dynmzstart", 400, ldynmzstart)
	call inputpar_geti_def("dynamic_grid", "dynmzinterval", 400, ldynmzinterval)
	call inputpar_geti_def("dynamic_grid", "mzmin", 7, lmzmin)	

	
	dynmx=ldynmx
	dynmxstart=ldynmxstart
	dynmxinterval=ldynmxinterval
	mxmin=lmxmin
	
	dynmy=ldynmy
	dynmystart=ldynmystart
	dynmyinterval=ldynmyinterval
	mymin=lmymin
	
	dynmz=ldynmz
	dynmzstart=ldynmzstart
	dynmzinterval=ldynmzinterval	
	mzmin=lmzmin
		
end subroutine read_input_dynamic_grid


	
!-------------------------------------------------------------------------------
! 				subroutine redist_x_user()	
!										
! Re-distribute domain allocation (along the x direction)
! 
!-------------------------------------------------------------------------------

subroutine redist_x_domain()
	
	implicit none
	
	logical, parameter :: tmstop=.true.

	integer :: i, n1, n, j, k, i1, i2, j1, j2, k1, k2 ! dummy
	
	integer :: iL
	integer :: jL, jR
	
	logical :: in
	
	integer :: error , ierr
	integer :: recvrank, sendrank, buffsize1, buffsize0
	integer :: proc, plustag, minustag
	integer :: status(statsize)

	! number of particles escaping domain
	integer :: lenoutplus, lenoutminus, leninminus, leninplus
	integer :: lenioninplus0, lenioninminus0, lenlecinplus0, lenlecinminus0
	
	integer(8) :: leftend, rightend
	
	! variables for minimum mx constraint
	integer :: nmin, nmin0
	
	integer :: mxsend1, mxsend2, mxrecv1, mxrecv2
	integer, allocatable, dimension(:) :: mxsendL1, mxsendL2, mxrecvR1, mxrecvR2
	integer, allocatable, dimension(:) :: mxsendR1, mxsendR2, mxrecvL1, mxrecvL2
	integer, allocatable, dimension(:) :: countsendR, countrecvL, countsendL, countrecvR
	
	! new mx after domain redistribution
	integer :: mxnew, loc
	integer(8) :: mxnewloc, mxnewcum
	integer(8) totpartall, totpartcum, totpartcum_orig, totpartcum_corr, totpart, totpartmean
	integer :: nplus, nminus, nplus0, nminus0
	
	! arrays in the whole box
	integer, allocatable, dimension(:) :: mxnewl
	integer*8, allocatable, dimension(:) :: mxcuml, mxnewcuml
	integer*8, allocatable, dimension(:) :: totpartl, totpartcuml
	
	real :: cellfactor
	
	if(.not. ((dynmx.eq. 1) .and. (lap .ge. dynmxstart) .and. & 
		(modulo((lap-dynmxstart),dynmxinterval) .eq. 1))) return
		
		if(rank .eq. 0) print*, 'domain x-redistribution starts'

		allocate(mxnewl(size0),mxcuml(size0),mxnewcuml(size0))
		allocate(totpartl(size0),totpartcuml(size0))		
				
		totpart=ions+lecs		! total particles in a local cpu
				
		! the first rank in a row where my cpu belongs
		iL=(rank/sizex)*sizex
    	
     	! the first and last ranks in a column my cpu belongs
     	jL=modulo(rank,sizex)
     	jR=size0-sizex+modulo(rank,sizex)
  
     	! ex) sizex=4, sizey=3, sizez=2  --> size0=24
     	
     	!   jR  20   21   22   23      
        !      -------------------
	    !iL 8  | 8  | 9  | 10 | 11 | 
	    !      -------------------		rank/(sizex*sizey)=0
	    !   4  | 4  |  5 |  6 | 7  |
	    !      -------------------
	    !   0  | 0  |  1 |  2 | 3  |
	    !      --------------------> x
	    !	jL	0	  1    2   3
	    		
     	!   jR  20   21   22   23        
        !      ------------------
	    !iL	20 | 20 | 21 | 22 | 23 | 	rank/(sizex*sizey)=1
	    !      -------------------
	    !   16 | 16 | 17 | 18 | 19 |
	    !      -------------------
	    !   12 | 12 | 13 | 14 | 15 |
	    !      --------------------> x
	    !	jL	0	  1    2   3   		
	    		
		call mpi_allgather(mxcum,1,mpi_integer8,mxcuml,1,mpi_integer8 &
     	   	  ,mpi_comm_world, error)    	  
 	    			  	    	
		call mpi_allgather(totpart,1,mpi_integer8,totpartl,1,mpi_integer8 &
     	   	  ,mpi_comm_world, error)
 
 	    call mpi_allreduce(totpart,totpartall,1,mpi_integer8,mpi_sum &
			  ,mpi_comm_world,ierr)

		! mean particles/each x silce		
 	    totpartmean=totpartall/sizex
 	    
 	    if(rank.eq.0) print*, 'totpartmean=',totpartmean
 	    
		! cumulative total particle in each x slice
	    totpartcum=0
		do j=1, modulo(rank,sizex)
			i1=j
			!top rank on the top floor
			i2=size0-sizex+modulo(i1-1,sizex)+1
			totpartcum=totpartcum+sum(totpartl(i1:i2:sizex))
		enddo
		totpartcum_orig=totpartcum
		
!		if(rank/sizex .eq.0) then
!				print*, 'rank, totpartcum',rank,totpartcum
!		endif
 	    	
		totpartcum_corr=totpartall*mxcum/(mx0-5)
			
		! cellfactor is a correction factor to compansate computation time for cells
		cellfactor=0. 	! select it between [0,1]

		100 continue
		call MPI_BARRIER(MPI_COMM_WORLD,ierr)	
	
		totpartcum=totpartcum_orig-cellfactor*(totpartcum_orig-totpartcum_corr)
		  	
     	call mpi_allgather(totpartcum,1,mpi_integer8,totpartcuml,1,mpi_integer8 &
	       	   	  ,mpi_comm_world, error)    	  
	           	   	  
	    ! for example, consider domains in a row  	   	  
	   	!--------------------------------  
	    !totpartcum| 0  | 100 | 150 | 180|  cumulative particles before redistribution
	    !--------------------------------
	    !old rank  |  0 |  1  |  2  | 3  |  domains before redistribution
	    !--------------------------------
	    !new rank  |0| 1  |  2  |   3    |  domains after redistribution
	    !--------------------------------
	    !loc       |0| 0  |  1  |   2    |  the old dom     ain where the left boundary of a new domain is located
	    !--------------------------------	
	    !totpartcum|0| 50 | 100 |  150   |  cumulative particles after redistribution
	    !--------------------------------		
	    
	    ! find the old domain where the left boundary of a new domain is located
	    ! note: cpus in the same column do the same calculation
!	    call MPI_BARRIER(MPI_COMM_WORLD,ierr)	
	    		
	    do i=1, sizex-1
	    
	  	    k=iL+i
	    
			if((totpartmean*jL .ge. totpartcuml(k)) .and. &
	   			(totpartmean*jL .lt. totpartcuml(k+1))) then
     			loc=i-1
    			exit  			
			endif
     	enddo
		if(totpartmean*jL .ge. totpartcuml(iL+sizex)) then
			loc=sizex-1
		endif	
	   	   	  	    		    
	    ! left and right boundary of the particles
	    ! These can change depending on models
	    leftend=3
	    rightend=mx0-2
	    
	    !mxnewloc is the position of the left boundary of a new domain
	    k=iL+loc
	    if(modulo(rank,sizex) .eq. 0) then
	    	mxnewloc=0
	    else
	    	if(loc .eq. 0) then 			
	    		mxnewloc=int(1.0*(3+mxcuml(k+2)-leftend)/ &
	    				     (totpartcuml(k+2)-totpartcuml(k+1))* &
	    				     totpartmean*jL+leftend-3)      
	    	else if(loc .eq. sizex-1) then
				mxnewloc=int(1.0*(rightend-3-mxcuml(k+1))/(totpartall-totpartcuml(k+1))* &
        	 		(totpartmean*jL-totpartcuml(k+1))+mxcuml(k+1))
        	else				
				mxnewloc=int(1.0*(mxcuml(k+2)-mxcuml(k+1))/ &
							 (totpartcuml(k+2)-totpartcuml(k+1))* &
							 (totpartmean*jL-totpartcuml(k+1))+mxcuml(k+1))			
			endif
		endif
		
		!mxnewcuml is cumulative of mxnew
	    call mpi_allgather(mxnewloc,1,mpi_integer8,mxnewcuml,1,mpi_integer8 &
     	   	  ,mpi_comm_world, error)
     	   	  
     	!calcuate new mx
     	if(jL .eq. sizex-1) then
     		mxnew=(mx0-5)-mxnewcuml(iL+sizex)+5
     	else
     		mxnew=mxnewcuml(iL+jL+2)-mxnewcuml(iL+jL+1)+5	   	  
     	endif
									
		call mpi_allgather(mxnew,1,mpi_integer,mxnewl,1,mpi_integer &
    	   	  ,mpi_comm_world, error) 	   	  
	
    	! cumulative new mx
		mxnewcum=sum(mxnewl(iL+1:rank+1)-5)-(mxnewl(rank+1)-5)		
   		call mpi_allgather(mxnewcum,1,mpi_integer8,mxnewcuml,1,mpi_integer8 &
     	   	  ,mpi_comm_world, error) 	   	    
!		!*******************************************************	   	  

		if(cellfactor.ge.1) goto 200
    	!*******************************************************
    	! check if one domain < mxmin and re-arrange mx
    	!*******************************************************		
        if(mxnew .lt. mxmin) then
        	nmin=1
        else
        	nmin=0
        endif  	
	    ! number of cpus less than mxmin
	    call mpi_allreduce(nmin,nmin0,1,mpi_integer,mpi_sum &
            	,mpi_comm_world,ierr)

    	    	
        if(nmin0 .gt. 0) then
        	cellfactor=cellfactor+0.1
        	if(nmin.ne.0 .and. rank/sizex.eq.0) then
          		print*, '** detect too small mx=',mxnew,'at rank',rank 
                print*, '** redistribute with a relaxation factor =',cellfactor
            endif
    	endif
    			
    	if(nmin0 .gt. 0) goto 100
    		
    	200 continue
    	call MPI_BARRIER(MPI_COMM_WORLD,ierr)
    	
    	!******************************************************* 	
        ! Now we calculate nright and nleft. 
        ! These determine the miximum numbers of domain exhange.
	    ! nright: the difference between the farthest cpu sent to the right and the current cpu
	    ! nleft: the difference between the farthest cpu sent to the left from the current cpu
	
	    ! find the new domain where the right boundary of an old domain located
	    do i=1, sizex
	    
	    	k=iL+i
	    	if((mx-2+mxcum .gt. 3+mxnewcuml(k)) .and. &
	   			(mx-2+mxcum .le. mxnewl(k)-2+mxnewcuml(k))) then
     			loc=i-1
     			exit 
			endif 			
     	enddo
     	if(jL .eq. sizex-1) loc=sizex-1	
     	nplus=loc-jL
     	if(nplus .lt. 0) nplus=0
     	
     	call MPI_BARRIER(MPI_COMM_WORLD,ierr)	
     	
     	call mpi_allreduce(nplus,nplus0,1,mpi_integer,mpi_max &
			,mpi_comm_world,ierr)
     	nplus=nplus0
     	
     	! find the new domain where the left boundary of an old domain located
	    do i=1, sizex
	    	k=iL+i
			if((3+mxcum .ge. 3+mxnewcuml(k)) .and. &
	   			(3+mxcum .lt. mxnewl(k)-2+mxnewcuml(k))) then
     			loc=i-1
     			exit
     		endif
     	enddo
     	if(jL .eq. 0) loc=0
     	nminus=jL-loc
     	if(nminus .lt. 0) nminus=0
     	call mpi_allreduce(nminus,nminus0,1,mpi_integer,mpi_max &
			,mpi_comm_world,ierr)
     	nminus=nminus0
   	
     	if(rank.eq.0) then
     			print*, mxnewl(1:sizex)
     			print*, '***** nplus, nminus :', nplus,nminus  
     	endif
     	!*******************************************************
     
     	if(nplus.eq.0 .and. nminus.eq.0) then
     			
     		deallocate(mxnewl,mxcuml,mxnewcuml)
     		deallocate(totpartl,totpartcuml)	
     		
     		return
     		
     	endif		
     			
        allocate(mxsendR1(max(nplus,1)),mxsendR2(max(nplus,1)))
        allocate(mxrecvL1(max(nplus,1)),mxrecvL2(max(nplus,1)))
        allocate(mxsendL1(max(nminus,1)),mxsendL2(max(nminus,1)))
        allocate(mxrecvR1(max(nminus,1)),mxrecvR2(max(nminus,1)))
        allocate(countsendR(max(nplus,1)),countrecvL(max(nplus,1)))
        allocate(countsendL(max(nminus,1)),countrecvR(max(nminus,1)))
        
        countsendR=0
        countsendL=0
        countrecvR=0
        countrecvL=0
        
        
        deallocate(curx,cury,curz)
               
		allocate(curx(mxnew,my,mz),cury(mxnew,my,mz),curz(mxnew,my,mz))
			
		! initialization
		curx=0. !binit*cos(btheta)
		cury=0. !binit*sin(btheta)*sin(bphi)
		curz=0. !binit*sin(btheta)*cos(bphi)
					
		proc=modulo(rank,sizex)+1

		call MPI_BARRIER(MPI_COMM_WORLD,ierr)

		
		!***************************************************
		! send to the plus and receive from the minus
		!***************************************************
		do i=1, nplus

		  sendrank=(rank/sizex)*sizex + modulo(rank+i,sizex) 
		  recvrank=(rank/sizex)*sizex + modulo(rank-i,sizex) 	
		
    	  !**********************
    	  ! send to the plus
		  !**********************
    	  if(proc+i .gt. sizex) then
    	  		  
    	  		mxsendR1(i)=1
    	  		mxsendR2(i)=0
    	  else
    	  		   
    	  		if(3+mxnewcuml(proc+i) .gt. 3+mxcum) then	
!------------------			!------------------			!------------------		
!   | n   |					!   |   n    |			    !   | n   |	
!------------------   or    !------------------   or	!------------------
!  	    |n+i |				!     |n+i|					!    		 |n+i|	
!------------------			!------------------			!------------------    	    	
	    		 	if(3+mxnewcuml(proc+i) .le. mx-2+mxcum) then	
	    				if(mxnewl(proc+i)-2+mxnewcuml(proc+i) .gt. mx-2+mxcum) then
	    					mxsendR1(i)=3+mxnewcuml(proc+i)-mxcum	  	      
	    					mxsendR2(i)=mx-2
	    				else 
	    					mxsendR1(i)=3+mxnewcuml(proc+i)-mxcum	
	    					mxsendR2(i)=mxnewl(proc+i)-2+mxnewcuml(proc+i)-mxcum  	   		
	    				endif
	    			else
	    				mxsendR1(i)=1
	    				mxsendR2(i)=0
	    			endif
	    		else
!------------------			!------------------			!------------------	
!        |  n  |			!     |  n  |				!     | n  |	
!------------------   or	!------------------   or   	!------------------
!  |n+i|					!  |n+i |					!  |   n+i    |	
!------------------	    	!------------------			!------------------				    			
	    		 	if(mxnewl(proc+i)-2+mxnewcuml(proc+i) .le. 3+mxcum) then
	    		 		mxsendR1(i)=1
	    		 		mxsendR2(i)=0	
		    	
	    		 	else if(mxnewl(proc+i)-2+mxnewcuml(proc+i) .le. mx-2+mxcum) then		
	    				mxsendR1(i)=3
	    				mxsendR2(i)=mxnewl(proc+i)-2+mxnewcuml(proc+i)-mxcum
	    			else		    		
	    		 		mxsendR1(i)=3
	    		 		mxsendR2(i)=mx-2
	    		 	endif
	    		endif
	    	
	      endif !if(proc+i .gt. sizey)

	      countsendR(i)=(mxsendR2(i)-mxsendR1(i)+1)*my*mz
	      
		 !**********************
		 ! receive from the minus
		 !**********************
    	  if(proc-i .lt. 1) then
    	  		mxrecvL1(i)=1
    	  		mxrecvL2(i)=0
    	  else
    	  		  
    	  	  if(3+mxnewcum .gt. 3+mxcuml(proc-i)) then	
!------------------			!------------------			!------------------		
!   | n-i |					!   | n-i  |			 	!   | n-i |	
!------------------   or    !------------------   or	!------------------
!  	   | n   |				!     |n |					!    		 | n |	
!------------------			!------------------			!------------------    	    	
	    			if(3+mxnewcum .le. mxl(proc-i)-2+mxcuml(proc-i)) then	
	    				if(mxnew-2+mxnewcum .gt. mxl(proc-i)-2+mxcuml(proc-i)) then
	    					mxrecvL1(i)=3  	      
	    					mxrecvL2(i)=mxl(proc-i)-2+mxcuml(proc-i)-mxnewcum
	    				else
	    					mxrecvL1(i)=3
	    					mxrecvL2(i)=mxnew-2  	   		
	    				endif
	    			else
	    				mxrecvL1(i)=1
	    				mxrecvL2(i)=0
	    			endif
	    		else
!------------------			!------------------			!------------------	
!       | n-i |				!     | n-i  |				!     |n-i |	
!------------------   or	!------------------   or   	!------------------
!  | n |					!  | n   |					!  |    n     |	
!------------------	    	!------------------			!------------------				    			
	    			if(mxnew-2+mxnewcum .le. 3+mxcuml(proc-i)) then
	    				mxrecvL1(i)=1
	    				mxrecvL2(i)=0	
		    	
	    			else if(mxnew-2+mxnewcum .le. mxl(proc-i)-2+mxcuml(proc-i)) then		
	    				mxrecvL1(i)=3+mxcuml(proc-i)-mxnewcum
	    				mxrecvL2(i)=mxnew-2
	    			else		    		
	    				mxrecvL1(i)=3+mxcuml(proc-i)-mxnewcum
	    				mxrecvL2(i)=mxl(proc-i)-2+mxcuml(proc-i)-mxnewcum
	    			endif
	    		endif
	    	
	      endif !if(proc-i .lt. 1)
	      
	  
	      countrecvL(i)=(mxrecvL2(i)-mxrecvL1(i)+1)*my*mz
	     
				
	   enddo	!do i=1, nplus
	      
	   	!***************************************************
		! send to the minus and recieve from the plus
		!***************************************************
		do i=1, nminus
 	  
		  sendrank=(rank/sizex)*sizex + modulo(rank-i,sizex) 
		  recvrank=(rank/sizex)*sizex + modulo(rank+i,sizex) 	
		
    	  !**********************
    	  ! send to the minus
		  !**********************
    	  if(proc-i .lt. 1) then
    	  		  
    	  		mxsendL1(i)=1
    	  		mxsendL2(i)=0
    	  else
    	  		   
	    	if(3+mxnewcuml(proc-i) .gt. 3+mxcum) then	
!------------------			!------------------			!------------------		
!   | n   |					!   |   n    |			    !   | n   |	
!------------------   or    !------------------   or	!------------------
!  	    |n-i |				!     |n-i|					!    		 |n-i|	
!------------------			!------------------			!------------------    	    	
	    		if(3+mxnewcuml(proc-i) .le. mx-2+mxcum) then	
	    			if(mxnewl(proc-i)-2+mxnewcuml(proc-i) .gt. mx-2+mxcum) then
 	    				mxsendL1(i)=3+mxnewcuml(proc-i)-mxcum	  	      
	    				mxsendL2(i)=mx-2
	    			else 
	    				mxsendL1(i)=3+mxnewcuml(proc-i)-mxcum	
	    				mxsendL2(i)=mxnewl(proc-i)-2+mxnewcuml(proc-i)-mxcum  	   		
	    			endif
	    		else
		    		mxsendL1(i)=1
		    		mxsendL2(i)=0
		    	endif
	    	else
!------------------			!------------------			!------------------	
!        |  n  |			!     |  n  |				!     | n  |	
!------------------   or	!------------------   or   	!------------------
!  |n-i|					!  |n-i |					!  |   n-i    |	
!------------------	    	!------------------			!------------------				    			
	    		if(mxnewl(proc-i)-2+mxnewcuml(proc-i) .le. 3+mxcum) then
	    			mxsendL1(i)=1
	    			mxsendL2(i)=0	
		    	
	    		else if(mxnewl(proc-i)-2+mxnewcuml(proc-i) .le. mx-2+mxcum) then		
	    			mxsendL1(i)=3
	    			mxsendL2(i)=mxnewl(proc-i)-2+mxnewcuml(proc-i)-mxcum
	    		else		    		
	    			mxsendL1(i)=3
	    			mxsendL2(i)=mx-2
	    		endif
	    	endif
	    	
	      endif !if(proc+i .gt. sizey)

	      countsendL(i)=(mxsendL2(i)-mxsendL1(i)+1)*my*mz
	      
		!***************************************************
		! receive from the plus
		!***************************************************
    	  if(proc+i .gt. sizex) then
    	  		mxrecvR1(i)=1
    	  		mxrecvR2(i)=0
    	  else
    	  		  
    	  	if(3+mxnewcum .gt. 3+mxcuml(proc+i)) then	
!------------------			!------------------			!------------------		
!   | n+i |					!   | n+i  |			 	!   | n+i |	
!------------------   or    !------------------   or	!------------------
!  	   | n   |				!     |n |					!    		 | n |	
!------------------			!------------------			!------------------    	    	
	    		if(3+mxnewcum .le. mxl(proc+i)-2+mxcuml(proc+i)) then	
	    			if(mxnew-2+mxnewcum .gt. mxl(proc+i)-2+mxcuml(proc+i)) then
 	    				mxrecvR1(i)=3  	      
	    				mxrecvR2(i)=mxl(proc+i)-2+mxcuml(proc+i)-mxnewcum
	    			else
	    				mxrecvR1(i)=3
	    				mxrecvR2(i)=mxnew-2  	   		
	    			endif
	    		else
		    		mxrecvR1(i)=1
		    		mxrecvR2(i)=0
		    	endif
	    	else
!------------------			!------------------			!------------------	
!       | n+i |				!     | n+i  |				!     |n+i |	
!------------------   or	!------------------   or   	!------------------
!  | n |					!  | n   |					!  |    n     |	
!------------------	    	!------------------			!------------------				    			
	    		if(mxnew-2+mxnewcum .le. 3+mxcuml(proc+i)) then
	    			mxrecvR1(i)=1
	    			mxrecvR2(i)=0	
		    	
	    		else if(mxnew-2+mxnewcum .le. mxl(proc+i)-2+mxcuml(proc+i)) then		
	    			mxrecvR1(i)=3+mxcuml(proc+i)-mxnewcum
	    			mxrecvR2(i)=mxnew-2
	    		else		    		
	    			mxrecvR1(i)=3+mxcuml(proc+i)-mxnewcum
	    			mxrecvR2(i)=mxl(proc+i)-2+mxcuml(proc+i)-mxnewcum
	    		endif
	    	endif
	    	
	      endif !if(proc-i .lt. 1)
	      
	      countrecvR(i)=(mxrecvR2(i)-mxrecvR1(i)+1)*my*mz
	      
	      	
	   enddo	!do i=1, nleft  
	      
	   	   
		!***************************************************************		
		!send to its own processor
		!***************************************************************		
		if(3+mxnewcum .gt. 3+mxcum) then
!------------------			!------------------			!------------------		
!   |  n   |				!   |  n   |			 	!   |  n  |	
!------------------   or    !------------------   or	!------------------
!  	    | n  |				!     | n|					!    		 | n |	
!------------------			!------------------			!------------------
			if(3+mxnewcum .le. mx-2+mxcum) then
				if(mxnew-2+mxnewcum .gt. mx-2+mxcum) then
					mxsend1=3+mxnewcum-mxcum
					mxsend2=mx-2
				else
					mxsend1=3+mxnewcum-mxcum
					mxsend2=mxnew-2+mxnewcum-mxcum
				endif
			else
				mxsend1=1
				mxsend2=0
			endif
		else
!------------------			!------------------			!------------------	
!         | n  |			!     | n   |				!     |  n  |	
!------------------   or	!------------------   or   	!------------------
!  | n |					!  |  n  |					!  |    n      |	
!------------------	    	!------------------			!------------------				    			
	    	if(mxnew-2+mxnewcum .le. 3+mxcum) then
	    		mxsend1=1
	    		mxsend2=0	
	    			
	    	else if(mxnew-2+mxnewcum .le. mx-2+mxcum) then		
	    		mxsend1=3
	    		mxsend2=mxnew-2+mxnewcum-mxcum
	    	else		    		
	    		mxsend1=3
	    		mxsend2=mx-2
	    	endif
	    endif

	    
		!***************************************************************		
		!receive from its own processor
		!***************************************************************	    	  		  
    	if(3+mxnewcum .gt. 3+mxcum) then	
!------------------			!------------------			!------------------		
!   |   n |					!   |  n  |			 		!   |  n  |	
!------------------   or    !------------------   or	!------------------
!  	   | n  |				!     |n|					!    		 | n |	
!------------------			!------------------			!------------------    	    	
	    	if(3+mxnewcum .le. mx-2+mxcum) then	
	    		if(mxnew-2+mxnewcum .gt. mx-2+mxcum) then
 	    			mxrecv1=3  	      
	    			mxrecv2=mx-2+mxcum-mxnewcum
	    		else
	    			mxrecv1=3
	    			mxrecv2=mxnew-2  	   		
	    		endif
	    	else
		    	mxrecv1=1
		    	mxrecv2=0
		    endif
	    else
!------------------			!------------------			!------------------	
!       |  n  |				!     |  n  |				!     |  n  |	
!------------------   or	!------------------   or   	!------------------
!  | n|						!  | n   |					!  |     n    |	
!------------------	    	!------------------			!------------------				    			
	    	if(mxnew-2+mxnewcum .le. 3+mxcum) then
	    		mxrecv1=1
	    		mxrecv2=0	
		    	
	    	else if(mxnew-2+mxnewcum .le. mx-2+mxcum) then		
	    		mxrecv1=3+mxcum-mxnewcum
	    		mxrecv2=mxnew-2
	    	else		    		
	    		mxrecv1=3+mxcum-mxnewcum
	    		mxrecv2=mx-2+mxcum-mxnewcum
	    	endif
	    endif
	    	
	     
		!*******************************************************
		! Exchange B fields, send to plus and receive from minus
		!*******************************************************
		do i=1, nplus

		  sendrank=(rank/sizex)*sizex + modulo(rank+i,sizex) 
		  recvrank=(rank/sizex)*sizex + modulo(rank-i,sizex) 	
	    
		  plustag=100


		  call MPI_SendRecv(bx(mxsendR1(i):mxsendR2(i),:,:),countsendR(i),mpi_read,sendrank,plustag, &
				curx(mxrecvL1(i):mxrecvL2(i),:,:),countrecvL(i),mpi_read,recvrank,plustag, &
				MPI_Comm_world,status,ierr)
		  call MPI_SendRecv(by(mxsendR1(i):mxsendR2(i),:,:),countsendR(i),mpi_read,sendrank,plustag, &
				cury(mxrecvL1(i):mxrecvL2(i),:,:),countrecvL(i),mpi_read,recvrank,plustag, &
				MPI_Comm_world,status,ierr)
		  call MPI_SendRecv(bz(mxsendR1(i):mxsendR2(i),:,:),countsendR(i),mpi_read,sendrank,plustag, &
				curz(mxrecvL1(i):mxrecvL2(i),:,:),countrecvL(i),mpi_read,recvrank,plustag, &
				MPI_Comm_world,status,ierr)
				
		enddo

		!*******************************************************
		! Exchange B fields, send to minus and receive from plus
		!*******************************************************
		do i=1, nminus
 	  
		  sendrank=(rank/sizex)*sizex + modulo(rank-i,sizex) 
		  recvrank=(rank/sizex)*sizex + modulo(rank+i,sizex) 	
		  
		  minustag=200
		  
		  call MPI_SendRecv(bx(mxsendL1(i):mxsendL2(i),:,:),countsendL(i),mpi_read,sendrank,minustag, &
				curx(mxrecvR1(i):mxrecvR2(i),:,:),countrecvR(i),mpi_read,recvrank,minustag, &
				MPI_Comm_world,status,ierr)
		  call MPI_SendRecv(by(mxsendL1(i):mxsendL2(i),:,:),countsendL(i),mpi_read,sendrank,minustag, &
				cury(mxrecvR1(i):mxrecvR2(i),:,:),countrecvR(i),mpi_read,recvrank,minustag, &
				MPI_Comm_world,status,ierr)
		  call MPI_SendRecv(bz(mxsendL1(i):mxsendL2(i),:,:),countsendL(i),mpi_read,sendrank,minustag, &
				curz(mxrecvR1(i):mxrecvR2(i),:,:),countrecvR(i),mpi_read,recvrank,minustag, &
				MPI_Comm_world,status,ierr)	
		  
		enddo
		


	    
		!***************************************************
		! Exchange B fields, send to its own processor
		!***************************************************
		
	    if(mxrecv2 .ne. 0) then
	    	if(modulo(rank,sizex) .eq. sizex-1) then
	    		curx(mxrecv1:mxrecv2+2,:,:)=bx(mxsend1:mxsend2+2,:,:)
	    		cury(mxrecv1:mxrecv2+2,:,:)=by(mxsend1:mxsend2+2,:,:)
	    		curz(mxrecv1:mxrecv2+2,:,:)=bz(mxsend1:mxsend2+2,:,:)	
	    	else if(modulo(rank,sizex) .eq. 0) then
	    		curx(mxrecv1-2:mxrecv2,:,:)=bx(mxsend1-2:mxsend2,:,:)
	    		cury(mxrecv1-2:mxrecv2,:,:)=by(mxsend1-2:mxsend2,:,:)
	    		curz(mxrecv1-2:mxrecv2,:,:)=bz(mxsend1-2:mxsend2,:,:)
	    	else	    			
	    		curx(mxrecv1:mxrecv2,:,:)=bx(mxsend1:mxsend2,:,:)
	    		cury(mxrecv1:mxrecv2,:,:)=by(mxsend1:mxsend2,:,:)
	    		curz(mxrecv1:mxrecv2,:,:)=bz(mxsend1:mxsend2,:,:)
	    	endif
	    endif
		
		deallocate(bx,by,bz)
		
		allocate(bx(mxnew,my,mz), by(mxnew,my,mz), bz(mxnew,my,mz))
		
		bx=curx
		by=cury
		bz=curz
		
		curx=0.
		cury=0.
		curz=0.

		call MPI_BARRIER(MPI_COMM_WORLD,ierr)		

		!*******************************************************
		! Exchange E fields, send to plus and receive from minus
		!*******************************************************
		do i=1, nplus

		  sendrank=(rank/sizex)*sizex + modulo(rank+i,sizex) 
		  recvrank=(rank/sizex)*sizex + modulo(rank-i,sizex) 	
	    
		  plustag=100


		  call MPI_SendRecv(ex(mxsendR1(i):mxsendR2(i),:,:),countsendR(i),mpi_read,sendrank,plustag, &
				curx(mxrecvL1(i):mxrecvL2(i),:,:),countrecvL(i),mpi_read,recvrank,plustag, &
				MPI_Comm_world,status,ierr)
		  call MPI_SendRecv(ey(mxsendR1(i):mxsendR2(i),:,:),countsendR(i),mpi_read,sendrank,plustag, &
				cury(mxrecvL1(i):mxrecvL2(i),:,:),countrecvL(i),mpi_read,recvrank,plustag, &
				MPI_Comm_world,status,ierr)
		  call MPI_SendRecv(ez(mxsendR1(i):mxsendR2(i),:,:),countsendR(i),mpi_read,sendrank,plustag, &
				curz(mxrecvL1(i):mxrecvL2(i),:,:),countrecvL(i),mpi_read,recvrank,plustag, &
				MPI_Comm_world,status,ierr)
				
		enddo
		
		!*******************************************************
		! Exchange E fields, send to minus and receive from plus
		!*******************************************************
		do i=1, nminus
 	  
		  sendrank=(rank/sizex)*sizex + modulo(rank-i,sizex) 
		  recvrank=(rank/sizex)*sizex + modulo(rank+i,sizex) 	
		  
		  minustag=200
		  
		  call MPI_SendRecv(ex(mxsendL1(i):mxsendL2(i),:,:),countsendL(i),mpi_read,sendrank,minustag, &
				curx(mxrecvR1(i):mxrecvR2(i),:,:),countrecvR(i),mpi_read,recvrank,minustag, &
				MPI_Comm_world,status,ierr)
		  call MPI_SendRecv(ey(mxsendL1(i):mxsendL2(i),:,:),countsendL(i),mpi_read,sendrank,minustag, &
				cury(mxrecvR1(i):mxrecvR2(i),:,:),countrecvR(i),mpi_read,recvrank,minustag, &
				MPI_Comm_world,status,ierr)
		  call MPI_SendRecv(ez(mxsendL1(i):mxsendL2(i),:,:),countsendL(i),mpi_read,sendrank,minustag, &
				curz(mxrecvR1(i):mxrecvR2(i),:,:),countrecvR(i),mpi_read,recvrank,minustag, &
				MPI_Comm_world,status,ierr)		  
		  
		enddo
	    
		!***************************************************
		! Exchange E fields, send to its own processor
		!***************************************************
		
	    if(mxrecv2 .ne. 0) then
	    	if(modulo(rank,sizex) .eq. sizex-1) then
	    		curx(mxrecv1:mxrecv2+2,:,:)=ex(mxsend1:mxsend2+2,:,:)
	    		cury(mxrecv1:mxrecv2+2,:,:)=ey(mxsend1:mxsend2+2,:,:)
	    		curz(mxrecv1:mxrecv2+2,:,:)=ez(mxsend1:mxsend2+2,:,:)	
	    	else if(modulo(rank,sizex) .eq. 0) then
	    		curx(mxrecv1-2:mxrecv2,:,:)=ex(mxsend1-2:mxsend2,:,:)
	    		cury(mxrecv1-2:mxrecv2,:,:)=ey(mxsend1-2:mxsend2,:,:)
	    		curz(mxrecv1-2:mxrecv2,:,:)=ez(mxsend1-2:mxsend2,:,:)
	    	else	    			
	    		curx(mxrecv1:mxrecv2,:,:)=ex(mxsend1:mxsend2,:,:)
	    		cury(mxrecv1:mxrecv2,:,:)=ey(mxsend1:mxsend2,:,:)
	    		curz(mxrecv1:mxrecv2,:,:)=ez(mxsend1:mxsend2,:,:)
	    	endif
	    endif				
				
		deallocate(ex,ey,ez)
		
		allocate(ex(mxnew,my,mz), ey(mxnew,my,mz), ez(mxnew,my,mz))
		
		ex=curx
		ey=cury
		ez=curz
		
		curx=0.
		cury=0.
		curz=0.
		
		call MPI_BARRIER(MPI_COMM_WORLD,ierr)
		
	   	   
		! buffer size for particle exchange
		buffsize1=maxval(totpartl)*2
		buffsize0=maxval(totpartl)*4
			

		deallocate(poutminus,poutplus,pinminus,pinplus)	
		
		allocate(poutminus(buffsize1),poutplus(buffsize1))
		allocate(pinminus(buffsize0),pinplus(buffsize0))
		
	
		!****************************************************************
		! send ions to the plus and recive from the minus
		!****************************************************************	
		lenioninminus0=0
		
		do i=1, nplus

			LenIonOutPlus=0

			sendrank=(rank/sizex)*sizex + modulo(rank+i,sizex) 
			recvrank=(rank/sizex)*sizex + modulo(rank-i,sizex) 	
		
			if(countsendR(i) .ne. 0) then			
			n=1
			if (ions.gt.0) then
		52  	continue
				n1=n			
				if((p(n1)%x .ge. mxsendR1(i)) .and. (p(n1)%x .le. mxsendR2(i))) then
					in=.false.
					p(n1)%x=p(n1)%x+(mxcum-mxnewcuml(sendrank+1))
				else
					in=.true.
				endif
							
				if(in) go to 58    	
				LenIonOutPlus=LenIonOutPlus+1
				call copyprt(p(n1),poutplus(LenIonOutPlus))
				call copyprt(p(ions),p(n1))		
				ions=ions-1
				n=n-1
				
		58		n=n+1
	
				if(n.le.ions) go to 52
			endif !if (ions.gt.0) then
			
			endif !if(countsendR(i) .ne. 0) then
			
			plustag=100
			
			call MPI_SendRecv(LenIonOutPlus,1,mpi_integer,sendrank,plustag, &
						LenIonInMinus,1,mpi_integer,recvrank,plustag, &
						MPI_Comm_World,status,ierr)
								
			lenoutplus=max(LenIonOutPlus,1)
			leninminus=max(LenIonInMinus,1)
			
			if(lenioninminus0+leninminus .gt. buffsize0) then
					print *, "after x redistribution: BUFFERSIZE SMALLER THAN PARTICLE LIST"
					print *,rank,":","buffsize", buffsize0, "lenioninminus",lenioninminus0+leninminus
					stop
			endif
			
			call MPI_SendRecv(poutplus(1:lenoutplus),lenoutplus,particletype,sendrank,plustag, &
				pinminus(lenioninminus0+1:lenioninminus0+leninminus),leninminus,particletype,recvrank,plustag, &
				MPI_Comm_World,status,ierr)
						
			lenioninminus0=lenioninminus0+LenIonInMinus	
			
			if(debug) print*, 'rank, ions sent to the i''th right, ',rank, i
			
			call MPI_BARRIER(MPI_COMM_WORLD,ierr)
				
		enddo !do i=1, nplus	
		
		!****************************************************************
		! send ions to the minus and recive from the plus
		!****************************************************************	
		lenioninplus0=0
		
		do i=1, nminus

			LenIonOutMinus=0

			sendrank=(rank/sizex)*sizex + modulo(rank-i,sizex) 
			recvrank=(rank/sizex)*sizex + modulo(rank+i,sizex) 	
		
			if(countsendL(i) .ne. 0) then			
			n=1
			if (ions.gt.0) then
		62  	continue
				n1=n			
				if((p(n1)%x .ge. mxsendL1(i)) .and. (p(n1)%x .le. mxsendL2(i))) then
					in=.false.
					p(n1)%x=p(n1)%x+(mxcum-mxnewcuml(sendrank+1))
				else
					in=.true.
				endif
							
				if(in) go to 68    	
				LenIonOutMinus=LenIonOutMinus+1
				call copyprt(p(n1),poutminus(LenIonOutMinus))
				call copyprt(p(ions),p(n1))		
				ions=ions-1
				n=n-1
				
		68		n=n+1
	
				if(n.le.ions) go to 62
			endif !if (ions.gt.0) then
			
			endif !if(countsendR(i) .ne. 0) then
			
			minustag=200
			
			call MPI_SendRecv(LenIonOutMinus,1,mpi_integer,sendrank,minustag, &
						LenIonInPlus,1,mpi_integer,recvrank,minustag, &
						MPI_Comm_World,status,ierr)
								
			lenoutminus=max(LenIonOutMinus,1)
			leninplus=max(LenIonInPlus,1)
			
			if(lenioninplus0+leninplus .gt. buffsize0) then
					print *, "after x redistribution: BUFFERSIZE SMALLER THAN PARTICLE LIST"
					print *,rank,":","buffsize", buffsize0, "lenioninplus",lenioninplus0+leninplus
					stop
			endif			
			
			call MPI_SendRecv(poutminus(1:lenoutminus),lenoutminus,particletype,sendrank,minustag, &
				pinplus(lenioninplus0+1:lenioninplus0+leninplus),leninplus,particletype,recvrank,minustag, &
				MPI_Comm_World,status,ierr)
						
			lenioninplus0=lenioninplus0+LenIonInPlus	

			if(debug) print*, 'ions sent to the i''th minus',i
			
			call MPI_BARRIER(MPI_COMM_WORLD,ierr)
					
	enddo !do i=1, nleft
	
		
		! shift particle positons in a new domain
		p(1:ions)%x=p(1:ions)%x+(mxcum-mxnewcum)	
	
		! transfer particles in a new domain
		do n=1,lenioninplus0
			ions=ions+1
			call copyprt(pinplus(n),p(ions))
		enddo		
				
		do n=1,lenioninminus0
			ions=ions+1
			call copyprt(pinminus(n),p(ions))
		enddo

		
		!****************************************************************
		! send electrons to the plus and recive from the minus
		!****************************************************************	
		lenlecinminus0=0
		
		do i=1, nplus

			LenLecOutPlus=0

			sendrank=(rank/sizex)*sizex + modulo(rank+i,sizex) 
			recvrank=(rank/sizex)*sizex + modulo(rank-i,sizex) 	
		
			if(countsendR(i) .ne. 0) then			
			n=maxhlf+1
			if (lecs.gt.0) then
		53 	continue
				n1=n
				if((p(n1)%x .ge. mxsendR1(i)) .and. (p(n1)%x .le. mxsendR2(i))) then
					in=.false.
					p(n1)%x=p(n1)%x+(mxcum-mxnewcuml(sendrank+1))
				else
					in=.true.
				endif
							
				if(in) go to 59    	
				LenLecOutPlus=LenLecOutPlus+1
				call copyprt(p(n1),poutplus(LenLecOutPlus))
				call copyprt(p(maxhlf+lecs),p(n1))		
				lecs=lecs-1
				n=n-1
				
		59		n=n+1
	
				if(n.le.maxhlf+lecs) go to 53
			endif !if (lecs.gt.0) then
			
			endif !if(countsendR(i) .ne. 0) then
			
			plustag=100
						
			call MPI_SendRecv(LenLecOutPlus,1,mpi_integer,sendrank,plustag, &
						LenLecInMinus,1,mpi_integer,recvrank,plustag, &
						MPI_Comm_World,status,ierr)
								
			lenoutplus=max(LenLecOutPlus,1)
			leninminus=max(LenLecInMinus,1)
			
			if(lenlecinminus0+leninminus .gt. buffsize0) then
					print *, "after x redistribution: BUFFERSIZE SMALLER THAN PARTICLE LIST"
					print *,rank,":","buffsize", buffsize0, "lenlecinminus",lenioninminus0+leninminus
					stop
			endif
		
			
			call MPI_SendRecv(poutplus(1:lenoutplus),lenoutplus,particletype,sendrank,plustag, &
				pinminus(lenlecinminus0+1:lenlecinminus0+leninminus),leninminus,particletype,recvrank,plustag, &
				MPI_Comm_World,status,ierr)
						
			lenlecinminus0=lenlecinminus0+LenLecInMinus

			if(debug) print*, 'rank, electrons sent to the i''th plus',rank, i

			call MPI_BARRIER(MPI_COMM_WORLD,ierr)			
			
		enddo !do i=1, nplus		
		
		
		!****************************************************************
		! send electrons to the minus and recive from the plus
		!****************************************************************	
		lenlecinplus0=0
		
		do i=1, nminus

			LenLecOutMinus=0

			sendrank=(rank/sizex)*sizex + modulo(rank-i,sizex) 
			recvrank=(rank/sizex)*sizex + modulo(rank+i,sizex) 	
		
			if(countsendL(i) .ne. 0) then			
			n=maxhlf+1
			if (lecs.gt.0) then
		63 	continue
				n1=n
				if((p(n1)%x .ge. mxsendL1(i)) .and. (p(n1)%x .le. mxsendL2(i))) then
					in=.false.
					p(n1)%x=p(n1)%x+(mxcum-mxnewcuml(sendrank+1))
				else
					in=.true.
				endif
						
				if(in) go to 69    	
				LenLecOutMinus=LenLecOutMinus+1
				call copyprt(p(n1),poutminus(LenLecOutMinus))
				call copyprt(p(maxhlf+lecs),p(n1))		
				lecs=lecs-1
				n=n-1
				
		69		n=n+1
	
				if(n.le.maxhlf+lecs) go to 63
			endif !if (lecs.gt.0) then
			
			endif !if(countsendR(i) .ne. 0) then
			
			minustag=200
			
			call MPI_SendRecv(LenLecOutMinus,1,mpi_integer,sendrank,minustag, &
						LenLecInPlus,1,mpi_integer,recvrank,minustag, &
						MPI_Comm_World,status,ierr)
								
			lenoutminus=max(LenLecOutMinus,1)
			leninplus=max(LenLecInPlus,1)
			
			if(lenlecinplus0+leninplus .gt. buffsize0) then
					print *, "after x redistribution: BUFFERSIZE SMALLER THAN PARTICLE LIST"
					print *,rank,":","buffsize", buffsize0, "lenlecinplus",lenlecinplus0+leninplus
					stop
			endif			
			
			call MPI_SendRecv(poutminus(1:lenoutminus),lenoutminus,particletype,sendrank,minustag, &
				pinplus(lenlecinplus0+1:lenlecinplus0+leninplus),leninplus,particletype,recvrank,minustag, &
				MPI_Comm_World,status,ierr)
						
			lenlecinplus0=lenlecinplus0+LenLecInPlus				
			
			if(debug)  print*, 'rank, electrons sent to the i''th left', rank,i	
			
			call MPI_BARRIER(MPI_COMM_WORLD,ierr)
			
		enddo !do i=1, nleft	
		
		
		! shift particle positons in a new domain
		p(maxhlf+1:maxhlf+lecs)%x=p(maxhlf+1:maxhlf+lecs)%x+(mxcum-mxnewcum)
		! transfer particles in a new domain
				
		do n=1,lenlecinminus0
			lecs=lecs+1
			call copyprt(pinminus(n),p(maxhlf+lecs))
		enddo
			
		do n=1,lenlecinplus0
			lecs=lecs+1
			call copyprt(pinplus(n),p(maxhlf+lecs))
		enddo		
		
			
		mx=mxnew
		mxall=mx
			
		call mpi_allgather(mx,1,mpi_integer,mxl,1,mpi_integer &
		  	   	  ,mpi_comm_world, error)
			
		i1=1
		i2=modulo(rank,sizex)+1
		mxcum=sum(mxl(i1:i2)-5)-(mxl(i2)-5)
		
		ix=1
		iy=mx
		iz=iy*my
		lot=iz*mz
		
#ifdef twoD
		iz=0
		lot=mx*my
#endif			
				

		deallocate(mxsendR1,mxsendR2,mxrecvL1,mxrecvL2)
        deallocate(mxsendL1,mxsendL2,mxrecvR1,mxrecvR2)
        deallocate(countsendR,countrecvL,countsendL,countrecvR)
  			
		deallocate(mxnewl,mxcuml,mxnewcuml)
		deallocate(totpartl,totpartcuml)
	
		deallocate(bufferin1,bufferin2)
		deallocate(bufferin1y,bufferin2y,bufferin1x,bufferin2x)
		deallocate(sendbufy,sendbufz)
		deallocate(temp,poutup,poutdwn,pinblw &
			, pinabv,poutlft,poutrgt,pinlft, pinrgt &
			, poutplus,poutminus,pinplus,pinminus,pall )
#ifdef filter2
		deallocate(xghost,yghost,zghost) !for filter2	
		
#ifndef	twoD
	call free_MPI_filter_datatypes()
#endif
	
#ifndef twoD
	call create_MPI_filter_datatypes(3,mx-3,3,my-3,3,mz-3)
#endif
#endif !if filter2
!#else
!	call create_MPI_filter_datatypes(3,mx-3,3,my-3,1,1)
!#endif
		
		
#ifndef twoD
		buffsize1 = max(10*int(ppc0*c*max((mx-5)*(my-5),1*(mx-5)*(mz-5))) &
				,10000)					
#else
		buffsize1 = max(3*int(ppc0*c*max((mx-5),(my-5))),60000)
!		if (splitparts) buffsize1=buffsize1*10 
#endif
		
		allocate(bufferin1(mx,my,2),bufferin2(mx,my,2), &
			bufferin1y(mx,2,mz),bufferin2y(mx,2,mz))
	
		allocate(bufferin1x(2,my,mz),bufferin2x(2,my,mz))
		allocate(sendbufy(mx,2,mz),sendbufz(mx,my,2))
		allocate(temp(mx,my,mz))
#ifdef filter2
		allocate(xghost(2*ntimes,my,mz),yghost(mx,2*ntimes,mz),zghost(mx,my,2*ntimes)) !for filter2
#endif
		allocate(poutup(buffsize1),poutdwn(buffsize1))
		allocate(pinabv(buffsize1),pinblw(buffsize1))
		allocate(poutrgt(buffsize1),poutlft(buffsize1))
		allocate(pinlft(buffsize1),pinrgt(buffsize1))
		allocate(poutplus(buffsize1),poutminus(buffsize1))
		allocate(pinminus(buffsize1),pinplus(buffsize1))
		allocate(pall(lot))
		
		
		pall=0

		bufferin1=0.
		bufferin2=0.
		bufferin1y=0.
		bufferin2y=0.
		bufferin1x=0.
		bufferin2x=0.
		sendbufy=0.
		sendbufz=0.
		temp=0.
				
		call bc_b1()
		if(debug) print*, 'bc done', rank
		call bc_e1()
	
		if(rank .eq. 0) print *,"domain x-redistribution done  "
		
		call MPI_BARRIER(MPI_COMM_WORLD,ierr)
							
	
end subroutine redist_x_domain
	
!-------------------------------------------------------------------------------
! 				subroutine redist_y_user()	
!										
! Re-distribute domain allocation (along the y direction)
! 
!-------------------------------------------------------------------------------

subroutine redist_y_domain()
	
	implicit none
	
	logical, parameter :: tmstop=.true.

	integer :: i, n1, n, j, k, i1, i2, j1, j2, k1, k2 ! dummy
	
	integer :: iL, iR
	integer :: jL, jR
	
	logical :: in
	
	integer :: error , ierr
	integer :: recvrank, sendrank, buffsize1, buffsize0
	integer :: proc, rgttag, lfttag
	integer :: status(statsize)

	! number of particles escaping domain
	integer :: lenoutrgt, lenoutlft, leninlft, leninrgt
	integer :: lenioninrgt0, lenioninlft0, lenlecinrgt0, lenlecinlft0
	
	integer :: leftend, rightend
	
	! variables for minimum my constraint
	integer :: nmin, nmin0
	
	integer :: mysend1, mysend2, myrecv1, myrecv2
	integer, allocatable, dimension(:) :: mysendL1, mysendL2, myrecvR1, myrecvR2
	integer, allocatable, dimension(:) :: mysendR1, mysendR2, myrecvL1, myrecvL2
	integer, allocatable, dimension(:) :: countsendR, countrecvL, countsendL, countrecvR
	
	! new mx after domain redistribution
	integer :: mynew, loc
	integer(8) :: mynewloc, mynewcum
	integer(8) :: totpartall, totpartcum, totpartcum_orig, totpartcum_corr, totpart, totpartmean
	integer :: nright, nleft, nright0, nleft0
	
	! arrays in the whole box
	integer, allocatable, dimension(:) :: mynewl
	integer*8, allocatable, dimension(:) :: mycuml, mynewcuml
	integer*8, allocatable, dimension(:) :: totpartl, totpartcuml
	
	real :: cellfactor
	
	integer, allocatable, dimension(:) :: row
	
	if(.not. ((dynmy.eq. 1) .and. (lap .ge. dynmystart) .and. & 
		(modulo(lap-dynmystart,dynmyinterval) .eq. 1))) return
			
		if(rank .eq. 0) print*, 'domain y-redistribution starts'

!		call MPI_BARRIER(MPI_COMM_WORLD,ierr)
!		call timer(40)
		
		allocate(mynewl(size0),mycuml(size0),mynewcuml(size0))
		allocate(totpartl(size0),totpartcuml(size0))		
		allocate(row(sizey))
		
				
		totpart=ions+lecs		! total particles in alocal cpu
						
		! the first and last ranks in a row where my cpu belongs
		iL=(rank/sizex)*sizex
    	
     	! the first and last ranks in a column my cpu belongs
     	jL=modulo(rank,sizex)
     	jR=size0-sizex+modulo(rank,sizex)
  
     	! ex) sizex=4, sizey=3, sizez=2  --> size0=24
     	
     	!   jR  20   21   22   23      
        !      -------------------
	    !iL 8  | 8  | 9  | 10 | 11 | 
	    !      -------------------		rank/(sizex*sizey)=0
	    !   4  | 4  |  5 |  6 | 7  |
	    !      -------------------
	    !   0  | 0  |  1 |  2 | 3  |
	    !      --------------------> x
	    !	jL	0	  1    2   3
	    		
     	!   jR  20   21   22   23        
        !      ------------------
	    !iL	20 | 20 | 21 | 22 | 23 | 	rank/(sizex*sizey)=1
	    !      -------------------
	    !   16 | 16 | 17 | 18 | 19 |
	    !      -------------------
	    !   12 | 12 | 13 | 14 | 15 |
	    !      --------------------> x
	    !	jL	0	  1    2   3   		
	    			    		
	    do j=1, sizey
			row(j)=sizex*(j-1)+1
		enddo
						
		call mpi_allgather(mycum,1,mpi_integer8,mycuml,1,mpi_integer8 &
     	   	  ,mpi_comm_world, error)    	  
 	    			  	    	
		call mpi_allgather(totpart,1,mpi_integer8,totpartl,1,mpi_integer8 &
     	   	  ,mpi_comm_world, error)
 
 	    call mpi_allreduce(totpart,totpartall,1,mpi_integer8,mpi_sum &
			  ,mpi_comm_world,ierr)
			  
		! mean particles/each y slice		
 	    totpartmean=int(1.*totpartall/sizey)
 	    
 	    if(rank.eq.0) print*, 'totpartmean/row=',totpartmean
 	    	
		! cumulative total particle in each x slice
	    totpartcum=0
		do k=1, sizez
			i1=(k-1)*(sizex*sizey)+1
			if(modulo(rank,sizex*sizey)/sizex .ne. 0) then
				! left rank
				i2=modulo(rank/sizex - 1,sizey)*sizex+&
					rank/(sizex*sizey)*(sizex*sizey) &
					+ modulo(rank,sizex)
				i2=modulo(i2,sizex*sizey)/sizex*sizex+sizex
				!right-end rank on the k'th floor
				i2=i2+(k-1)*(sizex*sizey)
				totpartcum=totpartcum+sum(totpartl(i1:i2))
			endif
		enddo
		totpartcum_orig=totpartcum	
						
		totpartcum_corr=totpartall*mycum/(my0-5)

		! cellfactor is a correction factor to compansate computation time for cells
		cellfactor=0. 	! select it between [0,1]
	
		100 continue
		call MPI_BARRIER(MPI_COMM_WORLD,ierr)	
				
		totpartcum=totpartcum_orig-cellfactor*(totpartcum_orig-totpartcum_corr)
		  	     	   	  
    	call mpi_allgather(totpartcum,1,mpi_integer8,totpartcuml,1,mpi_integer8 &
	       	   	  ,mpi_comm_world, error)    
	       	   	    	   	  
 
	       	   		           	   	  
	    ! for example, consider domains in a column
	   	!--------------------------------  
	    !totpartcum| 0  | 100 | 150 | 180|  cumulative particles before redistribution
	    !--------------------------------
	    !rank/sizex|  0 |  1  |  2  | 3  |  domains before redistribution
	    !--------------------------------
	    !rank/sizex|0| 1  |  2  |   3    |  domains after redistribution
	    !--------------------------------
	    !loc       |0| 0  |  1  |   2    |  the old domain where the left boundary of a new domain is located
	    !--------------------------------	
	    !totpartcum|0| 50 | 100 |  150   |  cumulative particles after redistribution
	    !--------------------------------			    		
	    		
	 		
	    ! find the old domain where the left boundary of a new domain is located
	    ! note: cpus in the same row do the same calculation
	    		
	    j1=modulo(rank,sizex*sizey)/sizex  
	    do i=1, sizey-1
	      	    
			if(totpartmean*j1 .ge. totpartcuml(row(i)) .and. &
	   			totpartmean*j1 .lt. totpartcuml(row(i+1))) then
     			loc=i-1
    			exit  			
			endif
     	enddo
		if(totpartmean*j1 .ge. totpartcuml(row(sizey))) then
			loc=sizey-1
		endif
	   	   	  	    		    
	    ! left and right boundary of the particles
	    ! These can change depending on models
	    leftend=3
	    rightend=my0-2
	    
	    !mynewloc is the position of the left boundary of a new domain
	    if(modulo(rank,sizex*sizey)/sizex .eq. 0) then
	    	mynewloc=0
	    else
	    	if(loc .eq. 0) then 			
	    		mynewloc=int(1.0*(3+mycuml(row(loc+2))-leftend)/ &
	    				     (totpartcuml(row(loc+2))-totpartcuml(row(loc+1)))* &
	    				     totpartmean*j1+leftend-3)      
	    	else if(loc .eq. sizey-1) then
				mynewloc=int(1.0*(rightend-3-mycuml(row(loc+1)))/(totpartall-totpartcuml(row(loc+1)))* &
        	 		(totpartmean*j1-totpartcuml(row(loc+1)))+mycuml(row(loc+1)))
        	else				
				mynewloc=int(1.0*(mycuml(row(loc+2))-mycuml(row(loc+1)))/ &
							 (totpartcuml(row(loc+2))-totpartcuml(row(loc+1)))* &
							 (totpartmean*j1-totpartcuml(row(loc+1)))+mycuml(row(loc+1)))			
			endif
		endif
			    		 
		!mynewcuml is cumulative of mynew
	    call mpi_allgather(mynewloc,1,mpi_integer8,mynewcuml,1,mpi_integer8 &
     	   	  ,mpi_comm_world, error)
     	   	  	  
     	   	      	   	  
     	!calcuate new my
     	if(modulo(rank,sizex*sizey)/sizex .eq. sizey-1) then
     		mynew=(my0-5)-mynewcuml(row(sizey))+5
     	else
     		mynew=mynewcuml(row(modulo(rank,sizex*sizey)/sizex+2))-&
     					mynewcuml(row(modulo(rank,sizex*sizey)/sizex+1))+5	   	  
     	endif
									
		call mpi_allgather(mynew,1,mpi_integer,mynewl,1,mpi_integer &
    	   	  ,mpi_comm_world, error) 	   	  
	
    	! cumulative new my
    	j1=1
    	j2=modulo(rank,sizex*sizey)/sizex*sizex+1
    	mynewcum=sum(mynewl(j1:j2:sizex)-5)-(mynewl(j2)-5)		 			
!    	mynewcum=sum(mynewl(jL+1:modulo(rank,sizex*sizey)+1:sizex)-5)-&
!    				(mynewl(modulo(rank,sizex*sizey)+1)-5)		
   		call mpi_allgather(mynewcum,1,mpi_integer8,mynewcuml,1,mpi_integer8 &
     	   	  ,mpi_comm_world, error) 	   	    
!		!*******************************************************	   	  

		
		if(cellfactor.ge.1) goto 200
    	!*******************************************************
    	! check if one domain < mxmin and re-arrange mx
    	!*******************************************************		
        if(mynew .lt. mymin) then
        	nmin=1
        else
        	nmin=0
        endif  	
	    ! number of cpus less than mxmin
	    call mpi_allreduce(nmin,nmin0,1,mpi_integer,mpi_sum &
            	,mpi_comm_world,ierr)
      
        if(nmin0 .gt. 0) then
        	cellfactor=cellfactor+0.1
        	if(nmin.ne.0 .and. modulo(rank,sizex).eq.0) then
          		print*, '** detect too small my=',mynew,'at rank',rank 
                print*, '** redistribute with a relaxation factor =',cellfactor
            endif
    		goto 100
    	endif
    			
    		
    	200 continue
    	call MPI_BARRIER(MPI_COMM_WORLD,ierr)
		
        ! Now we calculate nright and nleft. 
        ! These determine the miximum numbers of domain exhange.
	    ! nright: the difference between the farthest cpu sent to the right and the current cpu
	    ! nleft: the difference between the farthest cpu sent to the left from the current cpu
	
	    ! find the new domain where the right boundary of an old domain located
	    do i=1, sizey
	    
	    	if((my-2+mycum .gt. 3+mynewcuml(row(i))) .and. &
	   			(my-2+mycum .le. mynewl(row(i))-2+mynewcuml(row(i)))) then
     			loc=i-1
     			exit 
			endif 			
     	enddo
     	if(modulo(rank,sizex*sizey)/sizex .eq. sizey-1) loc=sizey-1	
     	nright=loc-modulo(rank,sizex*sizey)/sizex
     	if(nright .lt. 0) nright=0
     	
     	call MPI_BARRIER(MPI_COMM_WORLD,ierr)	
     	
     	call mpi_allreduce(nright,nright0,1,mpi_integer,mpi_max &
			,mpi_comm_world,ierr)
     	nright=nright0
     	
     	! find the new domain where the left boundary of an old domain located
	    do i=1, sizey

			if((3+mycum .ge. 3+mynewcuml(row(i))) .and. &
	   			(3+mycum .lt. mynewl(row(i))-2+mynewcuml(row(i)))) then
     			loc=i-1
     			exit
     		endif
     	enddo
     	if(modulo(rank,sizex*sizey)/sizex .eq. 0) loc=0
     	nleft=modulo(rank,sizex*sizey)/sizex-loc
     	if(nleft .lt. 0) nleft=0
     	call mpi_allreduce(nleft,nleft0,1,mpi_integer,mpi_max &
			,mpi_comm_world,ierr)
     	nleft=nleft0
   	
     	if(rank.eq.0) then
     			print*, mynewl(1:(sizey-1)*sizex+1:sizex)
     			print*, '***** nright, nleft :', nright,nleft  
     	endif
     	!*******************************************************
     
     	if(nleft.eq.0 .and. nright.eq.0) then
     			
     		deallocate(mynewl,mycuml,mynewcuml)
     		deallocate(totpartl,totpartcuml)		
     		deallocate(row)
     		
     		return
     		
     	endif		
     			
        allocate(mysendR1(max(nright,1)),mysendR2(max(nright,1)))
        allocate(myrecvL1(max(nright,1)),myrecvL2(max(nright,1)))
        allocate(mysendL1(max(nleft,1)),mysendL2(max(nleft,1)))
        allocate(myrecvR1(max(nleft,1)),myrecvR2(max(nleft,1)))
        allocate(countsendR(max(nright,1)),countrecvL(max(nright,1)))
        allocate(countsendL(max(nleft,1)),countrecvR(max(nleft,1)))
        
        countsendR=0
        countsendL=0
        countrecvR=0
        countrecvL=0
               
        deallocate(curx,cury,curz)
               
		allocate(curx(mx,mynew,mz),cury(mx,mynew,mz),curz(mx,mynew,mz))
			
		! initialization
		curx=0. !binit*cos(btheta)
		cury=0. !binit*sin(btheta)*sin(bphi)
		curz=0. !binit*sin(btheta)*cos(bphi)
					
		proc=modulo(rank,sizex*sizey)/sizex+1
		
		call MPI_BARRIER(MPI_COMM_WORLD,ierr)
		
		
		!***************************************************
		! send to the right and receive from the left
		!***************************************************
		do i=1, nright

		  sendrank=modulo(rank/sizex + i,sizey)*sizex+rank/(sizex*sizey)*(sizex*sizey) &
				+ modulo(rank,sizex) 
		  recvrank=modulo(rank/sizex - i,sizey)*sizex+rank/(sizex*sizey)*(sizex*sizey) &
				+ modulo(rank,sizex) 
		
    	  !**********************
    	  ! send to the right
		  !**********************
    	  if(proc+i .gt. sizey) then
    	  		  
    	  		mysendR1(i)=1
    	  		mysendR2(i)=0
    	  else
    	  		   
    	  		if(3+mynewcuml(row(proc+i)) .gt. 3+mycum) then	
!------------------			!------------------			!------------------		
!   | n   |					!   |   n    |			    !   | n   |	
!------------------   or    !------------------   or	!------------------
!  	    |n+i |				!     |n+i|					!    		 |n+i|	
!------------------			!------------------			!------------------    	    	
	    		 	if(3+mynewcuml(row(proc+i)) .le. my-2+mycum) then	
	    				if(mynewl(row(proc+i))-2+mynewcuml(row(proc+i)) .gt. my-2+mycum) then
	    					mysendR1(i)=3+mynewcuml(row(proc+i))-mycum	  	      
	    					mysendR2(i)=my-2
	    				else 
	    					mysendR1(i)=3+mynewcuml(row(proc+i))-mycum	
	    					mysendR2(i)=mynewl(row(proc+i))-2+mynewcuml(row(proc+i))-mycum  	   		
	    				endif
	    			else
	    				mysendR1(i)=1
	    				mysendR2(i)=0
	    			endif
	    		else
!------------------			!------------------			!------------------	
!        |  n  |			!     |  n  |				!     | n  |	
!------------------   or	!------------------   or   	!------------------
!  |n+i|					!  |n+i |					!  |   n+i    |	
!------------------	    	!------------------			!------------------				    			
	    		 	if(mynewl(row(proc+i))-2+mynewcuml(row(proc+i)) .le. 3+mycum) then
	    		 		mysendR1(i)=1
	    		 		mysendR2(i)=0	
		    	
	    		 	else if(mynewl(row(proc+i))-2+mynewcuml(row(proc+i)) .le. my-2+mycum) then		
	    				mysendR1(i)=3
	    				mysendR2(i)=mynewl(row(proc+i))-2+mynewcuml(row(proc+i))-mycum
	    			else		    		
	    		 		mysendR1(i)=3
	    		 		mysendR2(i)=my-2
	    		 	endif
	    		endif
	    	
	      endif !if(proc+i .gt. sizey)

	      countsendR(i)=mx*(mysendR2(i)-mysendR1(i)+1)*mz
	      
		 !**********************
		 ! receive from the left
		 !**********************
    	  if(proc-i .lt. 1) then
    	  		myrecvL1(i)=1
    	  		myrecvL2(i)=0
    	  else
    	  		  
    	  	  if(3+mynewcum .gt. 3+mycuml(row(proc-i))) then	
!------------------			!------------------			!------------------		
!   | n-i |					!   | n-i  |			 	!   | n-i |	
!------------------   or    !------------------   or	!------------------
!  	   | n   |				!     |n |					!    		 | n |	
!------------------			!------------------			!------------------    	    	
	    			if(3+mynewcum .le. myl(row(proc-i))-2+mycuml(row(proc-i))) then	
	    				if(mynew-2+mynewcum .gt. myl(row(proc-i))-2+mycuml(row(proc-i))) then
	    					myrecvL1(i)=3  	      
	    					myrecvL2(i)=myl(row(proc-i))-2+mycuml(row(proc-i))-mynewcum
	    				else
	    					myrecvL1(i)=3
	    					myrecvL2(i)=mynew-2  	   		
	    				endif
	    			else
	    				myrecvL1(i)=1
	    				myrecvL2(i)=0
	    			endif
	    		else
!------------------			!------------------			!------------------	
!       | n-i |				!     | n-i  |				!     |n-i |	
!------------------   or	!------------------   or   	!------------------
!  | n |					!  | n   |					!  |    n     |	
!------------------	    	!------------------			!------------------				    			
	    			if(mynew-2+mynewcum .le. 3+mycuml(row(proc-i))) then
	    				myrecvL1(i)=1
	    				myrecvL2(i)=0	
		    	
	    			else if(mynew-2+mynewcum .le. myl(row(proc-i))-2+mycuml(row(proc-i))) then		
	    				myrecvL1(i)=3+mycuml(row(proc-i))-mynewcum
	    				myrecvL2(i)=mynew-2
	    			else		    		
	    				myrecvL1(i)=3+mycuml(row(proc-i))-mynewcum
	    				myrecvL2(i)=myl(row(proc-i))-2+mycuml(row(proc-i))-mynewcum
	    			endif
	    		endif
	    	
	      endif !if(proc-i .lt. 1)
	      
	  
	      countrecvL(i)=mx*(myrecvL2(i)-myrecvL1(i)+1)*mz
	     
				
	   enddo	!do i=1, nright
	      
	   	   
		!***************************************************
		! send to the left and recieve from the right
		!***************************************************
		do i=1, nleft
 	  
		  sendrank=modulo(rank/sizex - i,sizey)*sizex+rank/(sizex*sizey)*(sizex*sizey) &
				+ modulo(rank,sizex) 
		  recvrank=modulo(rank/sizex + i,sizey)*sizex+rank/(sizex*sizey)*(sizex*sizey) &
				+ modulo(rank,sizex) 
		
    	  !**********************
    	  ! send to the left
		  !**********************
    	  if(proc-i .lt. 1) then
    	  		  
    	  		mysendL1(i)=1
    	  		mysendL2(i)=0
    	  else
    	  		   
	    	if(3+mynewcuml(row(proc-i)) .gt. 3+mycum) then	
!------------------			!------------------			!------------------		
!   | n   |					!   |   n    |			    !   | n   |	
!------------------   or    !------------------   or	!------------------
!  	    |n-i |				!     |n-i|					!    		 |n-i|	
!------------------			!------------------			!------------------    	    	
	    		if(3+mynewcuml(row(proc-i)) .le. my-2+mycum) then	
	    			if(mynewl(row(proc-i))-2+mynewcuml(row(proc-i)) .gt. my-2+mycum) then
 	    				mysendL1(i)=3+mynewcuml(row(proc-i))-mycum	  	      
	    				mysendL2(i)=my-2
	    			else 
	    				mysendL1(i)=3+mynewcuml(row(proc-i))-mycum	
	    				mysendL2(i)=mynewl(row(proc-i))-2+mynewcuml(row(proc-i))-mycum  	   		
	    			endif
	    		else
		    		mysendL1(i)=1
		    		mysendL2(i)=0
		    	endif
	    	else
!------------------			!------------------			!------------------	
!        |  n  |			!     |  n  |				!     | n  |	
!------------------   or	!------------------   or   	!------------------
!  |n-i|					!  |n-i |					!  |   n-i    |	
!------------------	    	!------------------			!------------------				    			
	    		if(mynewl(row(proc-i))-2+mynewcuml(row(proc-i)) .le. 3+mycum) then
	    			mysendL1(i)=1
	    			mysendL2(i)=0	
		    	
	    		else if(mynewl(row(proc-i))-2+mynewcuml(row(proc-i)) .le. my-2+mycum) then		
	    			mysendL1(i)=3
	    			mysendL2(i)=mynewl(row(proc-i))-2+mynewcuml(row(proc-i))-mycum
	    		else		    		
	    			mysendL1(i)=3
	    			mysendL2(i)=my-2
	    		endif
	    	endif
	    	
	      endif

	      countsendL(i)=mx*(mysendL2(i)-mysendL1(i)+1)*mz
	      
		!***************************************************
		! receive from the right
		!***************************************************
    	  if(proc+i .gt. sizey) then
    	  		myrecvR1(i)=1
    	  		myrecvR2(i)=0
    	  else
    	  		  
    	  	if(3+mynewcum .gt. 3+mycuml(row(proc+i))) then	
!------------------			!------------------			!------------------		
!   | n+i |					!   | n+i  |			 	!   | n+i |	
!------------------   or    !------------------   or	!------------------
!  	   | n   |				!     |n |					!    		 | n |	
!------------------			!------------------			!------------------    	    	
	    		if(3+mynewcum .le. myl(row(proc+i))-2+mycuml(row(proc+i))) then	
	    			if(mynew-2+mynewcum .gt. myl(row(proc+i))-2+mycuml(row(proc+i))) then
 	    				myrecvR1(i)=3  	      
	    				myrecvR2(i)=myl(row(proc+i))-2+mycuml(row(proc+i))-mynewcum
	    			else
	    				myrecvR1(i)=3
	    				myrecvR2(i)=mynew-2  	   		
	    			endif
	    		else
		    		myrecvR1(i)=1
		    		myrecvR2(i)=0
		    	endif
	    	else
!------------------			!------------------			!------------------	
!       | n+i |				!     | n+i  |				!     |n+i |	
!------------------   or	!------------------   or   	!------------------
!  | n |					!  | n   |					!  |    n     |	
!------------------	    	!------------------			!------------------				    			
	    		if(mynew-2+mynewcum .le. 3+mycuml(row(proc+i))) then
	    			myrecvR1(i)=1
	    			myrecvR2(i)=0	
		    	
	    		else if(mynew-2+mynewcum .le. myl(row(proc+i))-2+mycuml(row(proc+i))) then		
	    			myrecvR1(i)=3+mycuml(row(proc+i))-mynewcum
	    			myrecvR2(i)=mynew-2
	    		else		    		
	    			myrecvR1(i)=3+mycuml(row(proc+i))-mynewcum
	    			myrecvR2(i)=myl(row(proc+i))-2+mycuml(row(proc+i))-mynewcum
	    		endif
	    	endif
	    	
	      endif !if(proc-i .lt. 1)
	      
	      countrecvR(i)=mx*(myrecvR2(i)-myrecvR1(i)+1)*mz
	      		
		
	   enddo	!do i=1, nleft  
	      
	   
		!***************************************************************		
		!send to its own processor
		!***************************************************************		
		if(3+mynewcum .gt. 3+mycum) then
!------------------			!------------------			!------------------		
!   |  n   |				!   |  n   |			 	!   |  n  |	
!------------------   or    !------------------   or	!------------------
!  	    | n  |				!     | n|					!    		 | n |	
!------------------			!------------------			!------------------
			if(3+mynewcum .le. my-2+mycum) then
				if(mynew-2+mynewcum .gt. my-2+mycum) then
					mysend1=3+mynewcum-mycum
					mysend2=my-2
				else
					mysend1=3+mynewcum-mycum
					mysend2=mynew-2+mynewcum-mycum
				endif
			else
				mysend1=1
				mysend2=0
			endif
		else
!------------------			!------------------			!------------------	
!         | n  |			!     | n   |				!     |  n  |	
!------------------   or	!------------------   or   	!------------------
!  | n |					!  |  n  |					!  |    n      |	
!------------------	    	!------------------			!------------------				    			
	    	if(mynew-2+mynewcum .le. 3+mycum) then
	    		mysend1=1
	    		mysend2=0	
	    			
	    	else if(mynew-2+mynewcum .le. my-2+mycum) then		
	    		mysend1=3
	    		mysend2=mynew-2+mynewcum-mycum
	    	else		    		
	    		mysend1=3
	    		mysend2=my-2
	    	endif
	    endif

	    
		!***************************************************************		
		!receive from its own processor
		!***************************************************************	    	  		  
    	if(3+mynewcum .gt. 3+mycum) then	
!------------------			!------------------			!------------------		
!   |   n |					!   |  n  |			 		!   |  n  |	
!------------------   or    !------------------   or	!------------------
!  	   | n  |				!     |n|					!    		 | n |	
!------------------			!------------------			!------------------    	    	
	    	if(3+mynewcum .le. my-2+mycum) then	
	    		if(mynew-2+mynewcum .gt. my-2+mycum) then
 	    			myrecv1=3  	      
	    			myrecv2=my-2+mycum-mynewcum
	    		else
	    			myrecv1=3
	    			myrecv2=mynew-2  	   		
	    		endif
	    	else
		    	myrecv1=1
		    	myrecv2=0
		    endif
	    else
!------------------			!------------------			!------------------	
!       |  n  |				!     |  n  |				!     |  n  |	
!------------------   or	!------------------   or   	!------------------
!  | n|						!  | n   |					!  |     n    |	
!------------------	    	!------------------			!------------------				    			
	    	if(mynew-2+mynewcum .le. 3+mycum) then
	    		myrecv1=1
	    		myrecv2=0	
		    	
	    	else if(mynew-2+mynewcum .le. my-2+mycum) then		
	    		myrecv1=3+mycum-mynewcum
	    		myrecv2=mynew-2
	    	else		    		
	    		myrecv1=3+mycum-mynewcum
	    		myrecv2=my-2+mycum-mynewcum
	    	endif
	    endif
	    	
	    
		!*******************************************************
		! Exchange B fields, send to plus and receive from minus
		!*******************************************************
		do i=1, nright

		  sendrank=modulo(rank/sizex + i,sizey)*sizex+rank/(sizex*sizey)*(sizex*sizey) &
				+ modulo(rank,sizex) 
		  recvrank=modulo(rank/sizex - i,sizey)*sizex+rank/(sizex*sizey)*(sizex*sizey) &
				+ modulo(rank,sizex) 
				
				
		  rgttag=100

		  call MPI_SendRecv(bx(:,mysendR1(i):mysendR2(i),:),countsendR(i),mpi_read,sendrank,rgttag, &
				curx(:,myrecvL1(i):myrecvL2(i),:),countrecvL(i),mpi_read,recvrank,rgttag, &
				MPI_Comm_world,status,ierr)
		  call MPI_SendRecv(by(:,mysendR1(i):mysendR2(i),:),countsendR(i),mpi_read,sendrank,rgttag, &
				cury(:,myrecvL1(i):myrecvL2(i),:),countrecvL(i),mpi_read,recvrank,rgttag, &
				MPI_Comm_world,status,ierr)
		  call MPI_SendRecv(bz(:,mysendR1(i):mysendR2(i),:),countsendR(i),mpi_read,sendrank,rgttag, &
				curz(:,myrecvL1(i):myrecvL2(i),:),countrecvL(i),mpi_read,recvrank,rgttag, &
				MPI_Comm_world,status,ierr)
				
		enddo
		
		!*******************************************************
		! Exchange B fields, send to minus and receive from plus
		!*******************************************************
		do i=1, nleft
 	  
		  sendrank=modulo(rank/sizex - i,sizey)*sizex+rank/(sizex*sizey)*(sizex*sizey) &
				+ modulo(rank,sizex) 
		  recvrank=modulo(rank/sizex + i,sizey)*sizex+rank/(sizex*sizey)*(sizex*sizey) &
				+ modulo(rank,sizex) 
			  
		  lfttag=200
		  
		  call MPI_SendRecv(bx(:,mysendL1(i):mysendL2(i),:),countsendL(i),mpi_read,sendrank,lfttag, &
				curx(:,myrecvR1(i):myrecvR2(i),:),countrecvR(i),mpi_read,recvrank,lfttag, &
				MPI_Comm_world,status,ierr)
		  call MPI_SendRecv(by(:,mysendL1(i):mysendL2(i),:),countsendL(i),mpi_read,sendrank,lfttag, &
				cury(:,myrecvR1(i):myrecvR2(i),:),countrecvR(i),mpi_read,recvrank,lfttag, &
				MPI_Comm_world,status,ierr)
		  call MPI_SendRecv(bz(:,mysendL1(i):mysendL2(i),:),countsendL(i),mpi_read,sendrank,lfttag, &
				curz(:,myrecvR1(i):myrecvR2(i),:),countrecvR(i),mpi_read,recvrank,lfttag, &
				MPI_Comm_world,status,ierr)		  
		  
		enddo
	    
		!***************************************************
		! Exchange B fields, send to its own processor
		!***************************************************
		
	    if(myrecv2 .ne. 0) then
	    	if(modulo(rank,sizex*sizey)/sizex .eq. sizey-1) then
	    		curx(:,myrecv1:myrecv2+2,:)=bx(:,mysend1:mysend2+2,:)
	    		cury(:,myrecv1:myrecv2+2,:)=by(:,mysend1:mysend2+2,:)
	    		curz(:,myrecv1:myrecv2+2,:)=bz(:,mysend1:mysend2+2,:)	
	    	else if(modulo(rank,sizex*sizey)/sizex .eq. 0) then
	    		curx(:,myrecv1-2:myrecv2,:)=bx(:,mysend1-2:mysend2,:)
	    		cury(:,myrecv1-2:myrecv2,:)=by(:,mysend1-2:mysend2,:)
	    		curz(:,myrecv1-2:myrecv2,:)=bz(:,mysend1-2:mysend2,:)
	    	else	    			
	    		curx(:,myrecv1:myrecv2,:)=bx(:,mysend1:mysend2,:)
	    		cury(:,myrecv1:myrecv2,:)=by(:,mysend1:mysend2,:)
	    		curz(:,myrecv1:myrecv2,:)=bz(:,mysend1:mysend2,:)
	    	endif
	    endif
		
		deallocate(bx,by,bz)
		
		allocate(bx(mx,mynew,mz), by(mx,mynew,mz), bz(mx,mynew,mz))
		
		bx=curx
		by=cury
		bz=curz
	    
		curx=0.
		cury=0.
		curz=0.
		
		call MPI_BARRIER(MPI_COMM_WORLD,ierr)		

		!*******************************************************
		! Exchange E fields, send to plus and receive from minus
		!*******************************************************
		do i=1, nright

		  sendrank=modulo(rank/sizex + i,sizey)*sizex+rank/(sizex*sizey)*(sizex*sizey) &
				+ modulo(rank,sizex) 
		  recvrank=modulo(rank/sizex - i,sizey)*sizex+rank/(sizex*sizey)*(sizex*sizey) &
				+ modulo(rank,sizex) 
				
				
		  rgttag=100

		  call MPI_SendRecv(ex(:,mysendR1(i):mysendR2(i),:),countsendR(i),mpi_read,sendrank,rgttag, &
				curx(:,myrecvL1(i):myrecvL2(i),:),countrecvL(i),mpi_read,recvrank,rgttag, &
				MPI_Comm_world,status,ierr)
		  call MPI_SendRecv(ey(:,mysendR1(i):mysendR2(i),:),countsendR(i),mpi_read,sendrank,rgttag, &
				cury(:,myrecvL1(i):myrecvL2(i),:),countrecvL(i),mpi_read,recvrank,rgttag, &
				MPI_Comm_world,status,ierr)
		  call MPI_SendRecv(ez(:,mysendR1(i):mysendR2(i),:),countsendR(i),mpi_read,sendrank,rgttag, &
				curz(:,myrecvL1(i):myrecvL2(i),:),countrecvL(i),mpi_read,recvrank,rgttag, &
				MPI_Comm_world,status,ierr)
				
		enddo
		
		!*******************************************************
		! Exchange B fields, send to minus and receive from plus
		!*******************************************************
		do i=1, nleft
 	  
		  sendrank=modulo(rank/sizex - i,sizey)*sizex+rank/(sizex*sizey)*(sizex*sizey) &
				+ modulo(rank,sizex) 
		  recvrank=modulo(rank/sizex + i,sizey)*sizex+rank/(sizex*sizey)*(sizex*sizey) &
				+ modulo(rank,sizex) 
			  
		  lfttag=200
		  
		  call MPI_SendRecv(ex(:,mysendL1(i):mysendL2(i),:),countsendL(i),mpi_read,sendrank,lfttag, &
				curx(:,myrecvR1(i):myrecvR2(i),:),countrecvR(i),mpi_read,recvrank,lfttag, &
				MPI_Comm_world,status,ierr)
		  call MPI_SendRecv(ey(:,mysendL1(i):mysendL2(i),:),countsendL(i),mpi_read,sendrank,lfttag, &
				cury(:,myrecvR1(i):myrecvR2(i),:),countrecvR(i),mpi_read,recvrank,lfttag, &
				MPI_Comm_world,status,ierr)
		  call MPI_SendRecv(ez(:,mysendL1(i):mysendL2(i),:),countsendL(i),mpi_read,sendrank,lfttag, &
				curz(:,myrecvR1(i):myrecvR2(i),:),countrecvR(i),mpi_read,recvrank,lfttag, &
				MPI_Comm_world,status,ierr)		  
		  
		enddo
	    
		!***************************************************
		! Exchange B fields, send to its own processor
		!***************************************************
		
	    if(myrecv2 .ne. 0) then
	    	if(modulo(rank,sizex*sizey)/sizex .eq. sizey-1) then
	    		curx(:,myrecv1:myrecv2+2,:)=ex(:,mysend1:mysend2+2,:)
	    		cury(:,myrecv1:myrecv2+2,:)=ey(:,mysend1:mysend2+2,:)
	    		curz(:,myrecv1:myrecv2+2,:)=ez(:,mysend1:mysend2+2,:)	
	    	else if(modulo(rank,sizex*sizey)/sizex .eq. 0) then
	    		curx(:,myrecv1-2:myrecv2,:)=ex(:,mysend1-2:mysend2,:)
	    		cury(:,myrecv1-2:myrecv2,:)=ey(:,mysend1-2:mysend2,:)
	    		curz(:,myrecv1-2:myrecv2,:)=ez(:,mysend1-2:mysend2,:)
	    	else	    			
	    		curx(:,myrecv1:myrecv2,:)=ex(:,mysend1:mysend2,:)
	    		cury(:,myrecv1:myrecv2,:)=ey(:,mysend1:mysend2,:)
	    		curz(:,myrecv1:myrecv2,:)=ez(:,mysend1:mysend2,:)
	    	endif
	    endif
		
		deallocate(ex,ey,ez)
		
		allocate(ex(mx,mynew,mz), ey(mx,mynew,mz), ez(mx,mynew,mz))
		
		ex=curx
		ey=cury
		ez=curz
	    
		curx=0.
		cury=0.
		curz=0.
	    
		call MPI_BARRIER(MPI_COMM_WORLD,ierr)
	   	   
		! buffer size for particle exchange
	    ! pinlft and pinrgt have larger buffer size than poutrgt and poulgt
	    ! becuause one processor can receive particles from several cpus
		buffsize1=maxval(totpartl)*2
		buffsize0=maxval(totpartl)*4
			
		deallocate(poutrgt,poutlft,pinlft,pinrgt)	
		
		allocate(poutrgt(buffsize1),poutlft(buffsize1))
		allocate(pinlft(buffsize0),pinrgt(buffsize0))
	
		!****************************************************************
		! send ions to the right and recive from the left
		!****************************************************************	
		lenioninlft0=0
		
		do i=1, nright

			LenIonOutRgt=0

			sendrank=modulo(rank/sizex + i,sizey)*sizex+rank/(sizex*sizey)*(sizex*sizey) &
				+ modulo(rank,sizex) 
			recvrank=modulo(rank/sizex - i,sizey)*sizex+rank/(sizex*sizey)*(sizex*sizey) &
				+ modulo(rank,sizex) 

			if(countsendR(i) .ne. 0) then			
			n=1
			if (ions.gt.0) then
		52  	continue
				n1=n			
				if((p(n1)%y .ge. mysendR1(i)) .and. (p(n1)%y .le. mysendR2(i))) then
					in=.false.
					p(n1)%y=p(n1)%y+(mycum-mynewcuml(sendrank+1))
				else
					in=.true.
				endif
							
				if(in) go to 58    	
				LenIonOutRgt=LenIonOutRgt+1
				call copyprt(p(n1),poutrgt(LenIonOutRgt))
				call copyprt(p(ions),p(n1))		
				ions=ions-1
				n=n-1
				
		58		n=n+1
	
				if(n.le.ions) go to 52
			endif !if (ions.gt.0) then
			
			endif !if(countsendR(i) .ne. 0) then
			
			rgttag=100
			
			call MPI_SendRecv(LenIonOutRgt,1,mpi_integer,sendrank,rgttag, &
						LenIonInLft,1,mpi_integer,recvrank,rgttag, &
						MPI_Comm_World,status,ierr)
								
			lenoutrgt=max(LenIonOutRgt,1)
			leninlft=max(LenIonInLft,1)
			
			if(lenioninlft0+leninlft .gt. buffsize0) then
					print *, "after y redistribution: BUFFERSIZE SMALLER THAN PARTICLE LIST"
					print *,rank,":","buffsize", buffsize0, "lenioninlft",lenioninlft0+leninlft
					stop
			endif
			
			call MPI_SendRecv(poutrgt(1:lenoutrgt),lenoutrgt,particletype,sendrank,rgttag, &
				pinlft(lenioninlft0+1:lenioninlft0+leninlft),leninlft,particletype,recvrank,rgttag, &
				MPI_Comm_World,status,ierr)
						
			lenioninlft0=lenioninlft0+LenIonInLft	
			
			if(debug) print*, 'rank, ions sent to the i''th right, ',rank, i
			
			call MPI_BARRIER(MPI_COMM_WORLD,ierr)

				
		enddo !do i=1, nright	
		
		!****************************************************************
		! send ions to the left and recive from the right
		!****************************************************************	
		lenioninrgt0=0
		
		do i=1, nleft

			LenIonOutLft=0

			sendrank=modulo(rank/sizex - i,sizey)*sizex+rank/(sizex*sizey)*(sizex*sizey) &
				+ modulo(rank,sizex) 
			recvrank=modulo(rank/sizex + i,sizey)*sizex+rank/(sizex*sizey)*(sizex*sizey) &
				+ modulo(rank,sizex) 
		
			if(countsendL(i) .ne. 0) then			
			n=1
			if (ions.gt.0) then
		62  	continue
				n1=n			
				if((p(n1)%y .ge. mysendL1(i)) .and. (p(n1)%y .le. mysendL2(i))) then
					in=.false.
					p(n1)%y=p(n1)%y+(mycum-mynewcuml(sendrank+1))
				else
					in=.true.
				endif
							
				if(in) go to 68    	
				LenIonOutLft=LenIonOutLft+1
				call copyprt(p(n1),poutlft(LenIonOutLft))
				call copyprt(p(ions),p(n1))		
				ions=ions-1
				n=n-1
				
		68		n=n+1
	
				if(n.le.ions) go to 62
			endif !if (ions.gt.0) then
			
			endif !if(countsendR(i) .ne. 0) then
			
			lfttag=200
			
			call MPI_SendRecv(LenIonOutLft,1,mpi_integer,sendrank,lfttag, &
						LenIonInRgt,1,mpi_integer,recvrank,lfttag, &
						MPI_Comm_World,status,ierr)
								
			lenoutlft=max(LenIonOutLft,1)
			leninrgt=max(LenIonInRgt,1)
			
			if(lenioninrgt0+leninrgt .gt. buffsize0) then
					print *, "after y redistribution: BUFFERSIZE SMALLER THAN PARTICLE LIST"
					print *,rank,":","buffsize", buffsize0, "lenioninrgt",lenioninrgt0+leninrgt
					stop
			endif			
			
			call MPI_SendRecv(poutlft(1:lenoutlft),lenoutlft,particletype,sendrank,lfttag, &
				pinrgt(lenioninrgt0+1:lenioninrgt0+leninrgt),leninrgt,particletype,recvrank,lfttag, &
				MPI_Comm_World,status,ierr)
						
			lenioninrgt0=lenioninrgt0+LenIonInRgt	

			if(debug) print*, 'ions sent to the i''th left',i
			
			call MPI_BARRIER(MPI_COMM_WORLD,ierr)

					
	enddo !do i=1, nleft	
		
		! shift particle positons in a new domain
		p(1:ions)%y=p(1:ions)%y+(mycum-mynewcum)	
	
		! transfer particles in a new domain
		do n=1,lenioninrgt0
			ions=ions+1
			call copyprt(pinrgt(n),p(ions))
		enddo		
				
		do n=1,lenioninlft0
			ions=ions+1
			call copyprt(pinlft(n),p(ions))
		enddo

		
		!****************************************************************
		! send electrons to the right and recive from the left
		!****************************************************************	
		lenlecinlft0=0
		
		do i=1, nright

			LenLecOutRgt=0

			sendrank=modulo(rank/sizex + i,sizey)*sizex+rank/(sizex*sizey)*(sizex*sizey) &
				+ modulo(rank,sizex) 
			recvrank=modulo(rank/sizex - i,sizey)*sizex+rank/(sizex*sizey)*(sizex*sizey) &
				+ modulo(rank,sizex) 

			if(countsendR(i) .ne. 0) then			
			n=maxhlf+1
			if (lecs.gt.0) then
		53 	continue
				n1=n
				if((p(n1)%y .ge. mysendR1(i)) .and. (p(n1)%y .le. mysendR2(i))) then
					in=.false.
					p(n1)%y=p(n1)%y+(mycum-mynewcuml(sendrank+1))
				else
					in=.true.
				endif
							
				if(in) go to 59    	
				LenLecOutRgt=LenLecOutRgt+1
				call copyprt(p(n1),poutrgt(LenLecOutRgt))
				call copyprt(p(maxhlf+lecs),p(n1))		
				lecs=lecs-1
				n=n-1
				
		59		n=n+1
	
				if(n.le.maxhlf+lecs) go to 53
			endif !if (lecs.gt.0) then
			
			endif !if(countsendR(i) .ne. 0) then
			
			rgttag=100
						
			call MPI_SendRecv(LenLecOutRgt,1,mpi_integer,sendrank,rgttag, &
						LenLecInLft,1,mpi_integer,recvrank,rgttag, &
						MPI_Comm_World,status,ierr)
								
			lenoutrgt=max(LenLecOutRgt,1)
			leninlft=max(LenLecInLft,1)
			
			if(lenlecinlft0+leninlft .gt. buffsize0) then
					print *, "after x redistribution: BUFFERSIZE SMALLER THAN PARTICLE LIST"
					print *,rank,":","buffsize", buffsize0, "lenlecinlft",lenioninlft0+leninlft
					stop
			endif
		
			
			call MPI_SendRecv(poutrgt(1:lenoutrgt),lenoutrgt,particletype,sendrank,rgttag, &
				pinlft(lenlecinlft0+1:lenlecinlft0+leninlft),leninlft,particletype,recvrank,rgttag, &
				MPI_Comm_World,status,ierr)
						
			lenlecinlft0=lenlecinlft0+LenLecInLft

			if(debug) print*, 'rank, electrons sent to the i''th right',rank, i

			call MPI_BARRIER(MPI_COMM_WORLD,ierr)
			
			
		enddo !do i=1, nright		
		
		
		!****************************************************************
		! send electrons to the left and recive from the right
		!****************************************************************	
		lenlecinrgt0=0
		
		do i=1, nleft

			LenLecOutLft=0
                                            
			sendrank=modulo(rank/sizex - i,sizey)*sizex+rank/(sizex*sizey)*(sizex*sizey) &
				+ modulo(rank,sizex) 
			recvrank=modulo(rank/sizex + i,sizey)*sizex+rank/(sizex*sizey)*(sizex*sizey) &
				+ modulo(rank,sizex) 
		
			if(countsendL(i) .ne. 0) then			
			n=maxhlf+1
			if (lecs.gt.0) then
		63 	continue
				n1=n
				if((p(n1)%y .ge. mysendL1(i)) .and. (p(n1)%y .le. mysendL2(i))) then
					in=.false.
					p(n1)%y=p(n1)%y+(mycum-mynewcuml(sendrank+1))
				else
					in=.true.
				endif
						
				if(in) go to 69    	
				LenLecOutLft=LenLecOutLft+1
				call copyprt(p(n1),poutlft(LenLecOutLft))
				call copyprt(p(maxhlf+lecs),p(n1))		
				lecs=lecs-1
				n=n-1
				
		69		n=n+1
	
				if(n.le.maxhlf+lecs) go to 63
			endif !if (lecs.gt.0) then
			
			endif !if(countsendR(i) .ne. 0) then
			
			lfttag=200
			
			call MPI_SendRecv(LenLecOutLft,1,mpi_integer,sendrank,lfttag, &
						LenLecInRgt,1,mpi_integer,recvrank,lfttag, &
						MPI_Comm_World,status,ierr)
								
			lenoutlft=max(LenLecOutLft,1)
			leninrgt=max(LenLecInRgt,1)
			
			if(lenlecinrgt0+leninrgt .gt. buffsize0) then
					print *, "after x redistribution: BUFFERSIZE SMALLER THAN PARTICLE LIST"
					print *,rank,":","buffsize", buffsize0, "lenlecinrgt",lenlecinrgt0+leninrgt
					stop
			endif			
			
			call MPI_SendRecv(poutlft(1:lenoutlft),lenoutlft,particletype,sendrank,lfttag, &
				pinrgt(lenlecinrgt0+1:lenlecinrgt0+leninrgt),leninrgt,particletype,recvrank,lfttag, &
				MPI_Comm_World,status,ierr)
						
			lenlecinrgt0=lenlecinrgt0+LenLecInRgt				
			
			if(debug)  print*, 'rank, electrons sent to the i''th left', rank,i	
			
			call MPI_BARRIER(MPI_COMM_WORLD,ierr)

			
		enddo !do i=1, nleft	
		
		call MPI_BARRIER(MPI_COMM_WORLD,ierr)
		
		! shift particle positons in a new domain
		p(maxhlf+1:maxhlf+lecs)%y=p(maxhlf+1:maxhlf+lecs)%y+(mycum-mynewcum)
		! transfer particles in a new domain
				
		do n=1,lenlecinlft0
			lecs=lecs+1
			call copyprt(pinlft(n),p(maxhlf+lecs))
		enddo
			
		do n=1,lenlecinrgt0
			lecs=lecs+1
			call copyprt(pinrgt(n),p(maxhlf+lecs))
		enddo		
		
			
		my=mynew
		myall=my
			
		call mpi_allgather(my,1,mpi_integer,myl,1,mpi_integer &
		  	   	  ,mpi_comm_world, error)
			
		j1=1
		j2=modulo(rank,sizex*sizey)/sizex*sizex+1
		mycum=sum(myl(j1:j2:sizex)-5)-(myl(j2)-5)
	
		ix=1
		iy=mx
		iz=iy*my
		lot=iz*mz
		
#ifdef twoD
		iz=0
		lot=mx*my
#endif			
				 
     	deallocate(mysendR1,mysendR2,myrecvL1,myrecvL2)
        deallocate(mysendL1,mysendL2,myrecvR1,myrecvR2)
        deallocate(countsendR,countrecvL,countsendL,countrecvR)
  			
		deallocate(mynewl,mycuml,mynewcuml)
		deallocate(totpartl,totpartcuml)
		deallocate(row)	
		
		deallocate(bufferin1,bufferin2)
		deallocate(bufferin1y,bufferin2y,bufferin1x,bufferin2x)
		deallocate(sendbufy,sendbufz)
		deallocate(temp,poutup,poutdwn,pinblw &
			, pinabv,poutlft,poutrgt,pinlft, pinrgt &
			, poutplus,poutminus,pinplus,pinminus,pall )
#ifdef filter2
		deallocate(xghost,yghost,zghost) !for filter2	
		
#ifndef	twoD
	call free_MPI_filter_datatypes()
#endif
	
#ifndef twoD
	call create_MPI_filter_datatypes(3,mx-3,3,my-3,3,mz-3)
#endif
#endif !if filter2
!#else
!	call create_MPI_filter_datatypes(3,mx-3,3,my-3,1,1)
!#endif
		
		
#ifndef twoD
		buffsize1 = max(10*int(ppc0*c*max((mx-5)*(my-5),1*(mx-5)*(mz-5))) &
				,10000)					
#else
		buffsize1 = max(3*int(ppc0*c*max((mx-5),(my-5))),60000)
!		if (splitparts) buffsize1=buffsize1*10
#endif
		
		allocate(bufferin1(mx,my,2),bufferin2(mx,my,2), &
			bufferin1y(mx,2,mz),bufferin2y(mx,2,mz))
	
		allocate(bufferin1x(2,my,mz),bufferin2x(2,my,mz))
		allocate(sendbufy(mx,2,mz),sendbufz(mx,my,2))
		allocate(temp(mx,my,mz))
#ifdef filter2
		allocate(xghost(2*ntimes,my,mz), yghost(mx,2*ntimes,mz), zghost(mx,my,2*ntimes)) !for filter2
#endif
		allocate(poutup(buffsize1),poutdwn(buffsize1))
		allocate(pinabv(buffsize1),pinblw(buffsize1))
		allocate(poutrgt(buffsize1),poutlft(buffsize1))
		allocate(pinlft(buffsize1),pinrgt(buffsize1))
		allocate(poutplus(buffsize1),poutminus(buffsize1))
		allocate(pinminus(buffsize1),pinplus(buffsize1))
		allocate(pall(lot))
		
		
		pall=0

		bufferin1=0.
		bufferin2=0.
		bufferin1y=0.
		bufferin2y=0.
		bufferin1x=0.
		bufferin2x=0.
		sendbufy=0.
		sendbufz=0.
		temp=0.

		
		call bc_b1()
		if(debug) print*, 'bc done', rank
		call bc_e1()

		if(rank .eq. 0) print *,"domain redistribution done  "			
		
		call MPI_BARRIER(MPI_COMM_WORLD,ierr)
						
	
end subroutine redist_y_domain  	



!-------------------------------------------------------------------------------
! 				subroutine redist_z_user()	
!										
! Re-distribute domain allocation (along the y direction)
! 
!-------------------------------------------------------------------------------

subroutine redist_z_domain()
	
	implicit none
	
	logical, parameter :: tmstop=.true.

	integer :: i, n1, n, j, k, i1, i2, j1, j2, k1, k2 ! dummy
	
	integer :: iL, iR
	integer :: jL, jR
	
	logical :: in
	
	integer :: error , ierr
	integer :: recvrank, sendrank, buffsize1, buffsize0
	integer :: proc, uptag, dwntag
	integer :: status(statsize)

	! number of particles escaping domain
	integer :: lenoutup, lenoutdwn, leninblw, leninabv
	integer :: lenioninabv0, lenioninblw0, lenlecinabv0, lenlecinblw0
	
	integer :: leftend, rightend
	
	! variables for minimum mz constraint
	integer :: nmin, nmin0
	
	integer :: mzsend1, mzsend2, mzrecv1, mzrecv2
	integer, allocatable, dimension(:) :: mzsendL1, mzsendL2, mzrecvR1, mzrecvR2
	integer, allocatable, dimension(:) :: mzsendR1, mzsendR2, mzrecvL1, mzrecvL2
	integer, allocatable, dimension(:) :: countsendR, countrecvL, countsendL, countrecvR
	
	! new mx after domain redistribution
	integer :: mznew, loc
	integer(8) :: mznewloc, mznewcum
	integer(8) totpartall, totpartcum, totpartcum_orig, totpartcum_corr, totpart, totpartmean
	integer :: nup, ndwn, nup0, ndwn0
	
	! arrays in the whole box
	integer, allocatable, dimension(:) :: mznewl
	integer*8, allocatable, dimension(:) :: mzcuml, mznewcuml
	integer*8, allocatable, dimension(:) :: totpartl, totpartcuml
	
	real :: cellfactor
	
	integer, allocatable, dimension(:) :: row
		
	if(.not. ((dynmz.eq. 1) .and. (lap .ge. dynmzstart) .and. & 
		(modulo(lap-dynmzstart,dynmzinterval) .eq. 1))) return
			
		if(rank .eq. 0) print*, 'domain z-redistribution starts'

!		call MPI_BARRIER(MPI_COMM_WORLD,ierr)
!		call timer(40)
		
		allocate(mznewl(size0),mzcuml(size0),mznewcuml(size0))
		allocate(totpartl(size0),totpartcuml(size0))		
		allocate(row(sizez))
		
				
		totpart=ions+lecs		! total particles in alocal cpu
						
		! the first and last ranks in a row where mz cpu belongs
		iL=(rank/sizex)*sizex
    	
     	! the first and last ranks in a column mz cpu belongs
     	jL=modulo(rank,sizex)
     	jR=size0-sizex+modulo(rank,sizex)
  
     	! ex) sizex=4, sizey=3, sizez=2  --> size0=24
     	
     	!   jR  20   21   22   23      
        !      -------------------
	    !iL 8  | 8  | 9  | 10 | 11 | 
	    !      -------------------		rank/(sizex*sizey)=0
	    !   4  | 4  |  5 |  6 | 7  |
	    !      -------------------
	    !   0  | 0  |  1 |  2 | 3  |
	    !      --------------------> x
	    !	jL	0	  1    2   3
	    		
     	!   jR  20   21   22   23        
        !      ------------------
	    !iL	20 | 20 | 21 | 22 | 23 | 	rank/(sizex*sizey)=1
	    !      -------------------
	    !   16 | 16 | 17 | 18 | 19 |
	    !      -------------------
	    !   12 | 12 | 13 | 14 | 15 |
	    !      --------------------> x
	    !	jL	0	  1    2   3   		
	    			    		
	    do k=1, sizez
			row(k)=(sizex*sizey)*(k-1)+1
		enddo
						
		call mpi_allgather(mzcum,1,mpi_integer8,mzcuml,1,mpi_integer8 &
     	   	  ,mpi_comm_world, error)    	  
 	    			  	    	
		call mpi_allgather(totpart,1,mpi_integer8,totpartl,1,mpi_integer8 &
     	   	  ,mpi_comm_world, error)
 
 	    call mpi_allreduce(totpart,totpartall,1,mpi_integer8,mpi_sum &
			  ,mpi_comm_world,ierr)

		! mean particles/each z slice		
 	    totpartmean=totpartall/sizez
 	    
    
 	    if(rank.eq.0) print*, 'totpartmean/row=',totpartmean
 	 	    
		! cumulative total particle in each z slice
	    totpartcum=0
	    i2=rank/(sizex*sizey)*(sizex*sizey)
	    totpartcum=sum(totpartl(1:i2))

		totpartcum_orig=totpartcum	
						
		totpartcum_corr=totpartall*mzcum/(mz0-5)
		

		! cellfactor is a correction factor to compansate computation time for cells
		cellfactor=0. 	! select it between [0,1]
	
		100 continue
		call MPI_BARRIER(MPI_COMM_WORLD,ierr)	
				
		totpartcum=totpartcum_orig-cellfactor*(totpartcum_orig-totpartcum_corr)
		  	     	   	  
    	call mpi_allgather(totpartcum,1,mpi_integer8,totpartcuml,1,mpi_integer8 &
	       	   	  ,mpi_comm_world, error)    
	           	   	  
	    ! for example, consider domains in a column
	   	!--------------------------------  
	    !totpartcum| 0  | 100 | 150 | 180|  cumulative particles before redistribution
	    !--------------------------------
	    !rank/sizex|  0 |  1  |  2  | 3  |  domains before redistribution
	    !--------------------------------
	    !rank/sizex|0| 1  |  2  |   3    |  domains after redistribution
	    !--------------------------------
	    !loc       |0| 0  |  1  |   2    |  the old domain where the left boundary of a new domain is located
	    !--------------------------------	
	    !totpartcum|0| 50 | 100 |  150   |  cumulative particles after redistribution
	    !--------------------------------			    		
	    		
	 		
	    ! find the old domain where the left boundary of a new domain is located
	    ! note: cpus in the same row do the same calculation
	    		
	    k1=rank/(sizex*sizey)
	    do i=1, sizez-1
	      	    
			if(totpartmean*k1 .ge. totpartcuml(row(i)) .and. &
	   			totpartmean*k1 .lt. totpartcuml(row(i+1))) then
     			loc=i-1
    			exit  			
			endif
     	enddo
		if(totpartmean*k1 .ge. totpartcuml(row(sizez))) then
			loc=sizez-1
		endif
	   	   	  	    		    
	    ! left and right boundary of the particles
	    ! These can change depending on models
	    leftend=3
	    rightend=mz0-2
	    
	    !mznewloc is the position of the left boundary of a new domain
	    if(rank/(sizex*sizey) .eq. 0) then
	    	mznewloc=0
	    else
	    	if(loc .eq. 0) then 			
	    		mznewloc=int(1.0*(3+mzcuml(row(loc+2))-leftend)/ &
	    				     (totpartcuml(row(loc+2))-totpartcuml(row(loc+1)))* &
	    				     totpartmean*k1+leftend-3)      
	    	else if(loc .eq. sizez-1) then
				mznewloc=int(1.0*(rightend-3-mzcuml(row(loc+1)))/(totpartall-totpartcuml(row(loc+1)))* &
        	 		(totpartmean*k1-totpartcuml(row(loc+1)))+mzcuml(row(loc+1)))
        	else				
				mznewloc=int(1.0*(mzcuml(row(loc+2))-mzcuml(row(loc+1)))/ &
							 (totpartcuml(row(loc+2))-totpartcuml(row(loc+1)))* &
							 (totpartmean*k1-totpartcuml(row(loc+1)))+mzcuml(row(loc+1)))			
			endif
		endif
			    		 
		!mznewcuml is cumulative of mznew
	    call mpi_allgather(mznewloc,1,mpi_integer8,mznewcuml,1,mpi_integer8 &
     	   	  ,mpi_comm_world, error)
     	   	      	   	  
     	!calcuate new mz
     	if(rank/(sizex*sizey) .eq. sizez-1) then
     		mznew=(mz0-5)-mznewcuml(row(sizez))+5
     	else
     		mznew=mznewcuml(row(rank/(sizex*sizey)+2))-&
     					mznewcuml(row(rank/(sizex*sizey)+1))+5	   	  
     	endif
									
		call mpi_allgather(mznew,1,mpi_integer,mznewl,1,mpi_integer &
    	   	  ,mpi_comm_world, error) 	   	  
	
    	! cumulative new mz
    	k1=1
    	k2=rank/(sizex*sizey)*(sizex*sizey)+1
    	mznewcum=sum(mznewl(k1:k2:(sizex*sizey))-5)-(mznewl(k2)-5)		
 
  		call mpi_allgather(mznewcum,1,mpi_integer8,mznewcuml,1,mpi_integer8 &
     	   	  ,mpi_comm_world, error) 	   	    
!		!*******************************************************	   	  

		if(cellfactor.ge.1) goto 200
    	!*******************************************************
    	! check if one domain < mxmin and re-arrange mx
    	!*******************************************************		
        if(mznew .lt. mzmin) then
        	nmin=1
        else
        	nmin=0
        endif  	
	    ! number of cpus less than mxmin
	    call mpi_allreduce(nmin,nmin0,1,mpi_integer,mpi_sum &
            	,mpi_comm_world,ierr)
      
        if(nmin0 .gt. 0) then
        	cellfactor=cellfactor+0.1
        	if(nmin.ne.0 .and. modulo(rank,sizex).eq.0) then
          		print*, '** detect too small mz=',mznew,'at rank',rank 
                print*, '** redistribute with a relaxation factor =',cellfactor
            endif
    		goto 100
    	endif
    			
    		
    	200 continue
    	call MPI_BARRIER(MPI_COMM_WORLD,ierr)
		
        ! Now we calculate nright and nleft. 
        ! These determine the miximum numbers of domain exhange.
	    ! nright: the difference between the farthest cpu sent to the right and the current cpu
	    ! nleft: the difference between the farthest cpu sent to the left from the current cpu
	
	    ! find the new domain where the right boundary of an old domain located
	    do i=1, sizez
	    
	    	if((mz-2+mzcum .gt. 3+mznewcuml(row(i))) .and. &
	   			(mz-2+mzcum .le. mznewl(row(i))-2+mznewcuml(row(i)))) then
     			loc=i-1
     			exit 
			endif 			
     	enddo
     	if(rank/(sizex*sizey) .eq. sizez-1) loc=sizez-1	
     	nup=loc-rank/(sizex*sizey)
     	if(nup .lt. 0) nup=0
     	
     	call mpi_allreduce(nup,nup0,1,mpi_integer,mpi_max &
			,mpi_comm_world,ierr)
     	nup=nup0
     	
     	! find the new domain where the left boundary of an old domain located
	    do i=1, sizez

			if((3+mzcum .ge. 3+mznewcuml(row(i))) .and. &
	   			(3+mzcum .lt. mznewl(row(i))-2+mznewcuml(row(i)))) then
     			loc=i-1
     			exit
     		endif
     	enddo
     	if(rank/(sizex*sizey) .eq. 0) loc=0
     	ndwn=rank/(sizex*sizey)-loc
     	if(ndwn .lt. 0) ndwn=0
     	call mpi_allreduce(ndwn,ndwn0,1,mpi_integer,mpi_max &
			,mpi_comm_world,ierr)
     	ndwn=ndwn0
   	
     	if(rank.eq.0) then
     			print*, mznewl(1:(sizez-1)*(sizex*sizey)+1:(sizex*sizey))
     			print*, '***** nup, ndwn :', nup,ndwn 
     	endif
     	!*******************************************************
     
     	if(nup.eq.0 .and. ndwn.eq.0) then
     			
     		deallocate(mznewl,mzcuml,mznewcuml)
     		deallocate(totpartl,totpartcuml)		
     		deallocate(row)
     		
     		return
     		
     	endif		
     			
     			
        allocate(mzsendR1(max(nup,1)),mzsendR2(max(nup,1)))
        allocate(mzrecvL1(max(nup,1)),mzrecvL2(max(nup,1)))
        allocate(mzsendL1(max(ndwn,1)),mzsendL2(max(ndwn,1)))
        allocate(mzrecvR1(max(ndwn,1)),mzrecvR2(max(ndwn,1)))
        allocate(countsendR(max(nup,1)),countrecvL(max(nup,1)))
        allocate(countsendL(max(ndwn,1)),countrecvR(max(ndwn,1)))
        
        countsendR=0
        countsendL=0
        countrecvR=0
        countrecvL=0
               
        deallocate(curx,cury,curz)
               
		allocate(curx(mx,my,mznew),cury(mx,my,mznew),curz(mx,my,mznew))
			
		! initialization
		curx=0. !binit*cos(btheta)
		cury=0. !binit*sin(btheta)*sin(bphi)
		curz=0. !binit*sin(btheta)*cos(bphi)
	
		proc=rank/(sizex*sizey)+1
		
		call MPI_BARRIER(MPI_COMM_WORLD,ierr)
		
		
		!***************************************************
		! send to the up and receive from the down
		!***************************************************
		do i=1, nup

		  sendrank=modulo((rank/(sizex*sizey) + i),sizez)*(sizex*sizey) + & 
		 		  modulo(rank,sizex*sizey)
		  recvrank=modulo((rank/(sizex*sizey) - i),sizez)*(sizex*sizey) + & 
			      modulo(rank,sizex*sizey)
		
    	  !**********************
    	  ! send to the up
		  !**********************
    	  if(proc+i .gt. sizez) then
    	  		  
    	  		mzsendR1(i)=1
    	  		mzsendR2(i)=0
    	  else
    	  		   
    	  		if(3+mznewcuml(row(proc+i)) .gt. 3+mzcum) then	
!------------------			!------------------			!------------------		
!   | n   |					!   |   n    |			    !   | n   |	
!------------------   or    !------------------   or	!------------------
!  	    |n+i |				!     |n+i|					!    		 |n+i|	
!------------------			!------------------			!------------------    	    	
	    		 	if(3+mznewcuml(row(proc+i)) .le. mz-2+mzcum) then	
	    				if(mznewl(row(proc+i))-2+mznewcuml(row(proc+i)) .gt. mz-2+mzcum) then
	    					mzsendR1(i)=3+mznewcuml(row(proc+i))-mzcum	  	      
	    					mzsendR2(i)=mz-2
	    				else 
	    					mzsendR1(i)=3+mznewcuml(row(proc+i))-mzcum	
	    					mzsendR2(i)=mznewl(row(proc+i))-2+mznewcuml(row(proc+i))-mzcum  	   		
	    				endif
	    			else
	    				mzsendR1(i)=1
	    				mzsendR2(i)=0
	    			endif
	    		else
!------------------			!------------------			!------------------	
!        |  n  |			!     |  n  |				!     | n  |	
!------------------   or	!------------------   or   	!------------------
!  |n+i|					!  |n+i |					!  |   n+i    |	
!------------------	    	!------------------			!------------------				    			
	    		 	if(mznewl(row(proc+i))-2+mznewcuml(row(proc+i)) .le. 3+mzcum) then
	    		 		mzsendR1(i)=1
	    		 		mzsendR2(i)=0	
		    	
	    		 	else if(mznewl(row(proc+i))-2+mznewcuml(row(proc+i)) .le. mz-2+mzcum) then		
	    				mzsendR1(i)=3
	    				mzsendR2(i)=mznewl(row(proc+i))-2+mznewcuml(row(proc+i))-mzcum
	    			else		    		
	    		 		mzsendR1(i)=3
	    		 		mzsendR2(i)=mz-2
	    		 	endif
	    		endif
	    	
	      endif !if(proc+i .gt. sizey)

	      countsendR(i)=mx*my*(mzsendR2(i)-mzsendR1(i)+1)
	      
		 !**********************
		 ! receive from the down
		 !**********************
    	  if(proc-i .lt. 1) then
    	  		mzrecvL1(i)=1
    	  		mzrecvL2(i)=0
    	  else
    	  		  
    	  	  if(3+mznewcum .gt. 3+mzcuml(row(proc-i))) then	
!------------------			!------------------			!------------------		
!   | n-i |					!   | n-i  |			 	!   | n-i |	
!------------------   or    !------------------   or	!------------------
!  	   | n   |				!     |n |					!    		 | n |	
!------------------			!------------------			!------------------    	    	
	    			if(3+mznewcum .le. mzl(row(proc-i))-2+mzcuml(row(proc-i))) then	
	    				if(mznew-2+mznewcum .gt. mzl(row(proc-i))-2+mzcuml(row(proc-i))) then
	    					mzrecvL1(i)=3  	      
	    					mzrecvL2(i)=mzl(row(proc-i))-2+mzcuml(row(proc-i))-mznewcum
	    				else
	    					mzrecvL1(i)=3
	    					mzrecvL2(i)=mznew-2  	   		
	    				endif
	    			else
	    				mzrecvL1(i)=1
	    				mzrecvL2(i)=0
	    			endif
	    		else
!------------------			!------------------			!------------------	
!       | n-i |				!     | n-i  |				!     |n-i |	
!------------------   or	!------------------   or   	!------------------
!  | n |					!  | n   |					!  |    n     |	
!------------------	    	!------------------			!------------------				    			
	    			if(mznew-2+mznewcum .le. 3+mzcuml(row(proc-i))) then
	    				mzrecvL1(i)=1
	    				mzrecvL2(i)=0	
		    	
	    			else if(mznew-2+mznewcum .le. mzl(row(proc-i))-2+mzcuml(row(proc-i))) then		
	    				mzrecvL1(i)=3+mzcuml(row(proc-i))-mznewcum
	    				mzrecvL2(i)=mznew-2
	    			else		    		
	    				mzrecvL1(i)=3+mzcuml(row(proc-i))-mznewcum
	    				mzrecvL2(i)=mzl(row(proc-i))-2+mzcuml(row(proc-i))-mznewcum
	    			endif
	    		endif
	    	
	      endif !if(proc-i .lt. 1)
	      
	  
	      countrecvL(i)=mx*my*(mzrecvL2(i)-mzrecvL1(i)+1)
	     
					
	   enddo	!do i=1, nright
	      	   	   
		!***************************************************
		! send to the down and recieve from the up
		!***************************************************
		do i=1, ndwn
 	  
		  sendrank=modulo((rank/(sizex*sizey) - i),sizez)*(sizex*sizey) + & 
			      modulo(rank,sizex*sizey)
		  recvrank=modulo((rank/(sizex*sizey) + i),sizez)*(sizex*sizey) + & 
		 		  modulo(rank,sizex*sizey) 
		
    	  !**********************
    	  ! send to the down
		  !**********************
    	  if(proc-i .lt. 1) then
    	  		  
    	  		mzsendL1(i)=1
    	  		mzsendL2(i)=0
    	  else
    	  		   
	    	if(3+mznewcuml(row(proc-i)) .gt. 3+mzcum) then	
!------------------			!------------------			!------------------		
!   | n   |					!   |   n    |			    !   | n   |	
!------------------   or    !------------------   or	!------------------
!  	    |n-i |				!     |n-i|					!    		 |n-i|	
!------------------			!------------------			!------------------    	    	
	    		if(3+mznewcuml(row(proc-i)) .le. mz-2+mzcum) then	
	    			if(mznewl(row(proc-i))-2+mznewcuml(row(proc-i)) .gt. mz-2+mzcum) then
 	    				mzsendL1(i)=3+mznewcuml(row(proc-i))-mzcum	  	      
	    				mzsendL2(i)=mz-2
	    			else 
	    				mzsendL1(i)=3+mznewcuml(row(proc-i))-mzcum	
	    				mzsendL2(i)=mznewl(row(proc-i))-2+mznewcuml(row(proc-i))-mzcum  	   		
	    			endif
	    		else
		    		mzsendL1(i)=1
		    		mzsendL2(i)=0
		    	endif
	    	else
!------------------			!------------------			!------------------	
!        |  n  |			!     |  n  |				!     | n  |	
!------------------   or	!------------------   or   	!------------------
!  |n-i|					!  |n-i |					!  |   n-i    |	
!------------------	    	!------------------			!------------------				    			
	    		if(mznewl(row(proc-i))-2+mznewcuml(row(proc-i)) .le. 3+mzcum) then
	    			mzsendL1(i)=1
	    			mzsendL2(i)=0	
		    	
	    		else if(mznewl(row(proc-i))-2+mznewcuml(row(proc-i)) .le. mz-2+mzcum) then		
	    			mzsendL1(i)=3
	    			mzsendL2(i)=mznewl(row(proc-i))-2+mznewcuml(row(proc-i))-mzcum
	    		else		    		
	    			mzsendL1(i)=3
	    			mzsendL2(i)=mz-2
	    		endif
	    	endif
	    	
	      endif

	      countsendL(i)=mx*my*(mzsendL2(i)-mzsendL1(i)+1)
	      
		!***************************************************
		! receive from the up
		!***************************************************
    	  if(proc+i .gt. sizez) then
    	  		mzrecvR1(i)=1
    	  		mzrecvR2(i)=0
    	  else
    	  		  
    	  	if(3+mznewcum .gt. 3+mzcuml(row(proc+i))) then	
!------------------			!------------------			!------------------		
!   | n+i |					!   | n+i  |			 	!   | n+i |	
!------------------   or    !------------------   or	!------------------
!  	   | n   |				!     |n |					!    		 | n |	
!------------------			!------------------			!------------------    	    	
	    		if(3+mznewcum .le. mzl(row(proc+i))-2+mzcuml(row(proc+i))) then	
	    			if(mznew-2+mznewcum .gt. mzl(row(proc+i))-2+mzcuml(row(proc+i))) then
 	    				mzrecvR1(i)=3  	      
	    				mzrecvR2(i)=mzl(row(proc+i))-2+mzcuml(row(proc+i))-mznewcum
	    			else
	    				mzrecvR1(i)=3
	    				mzrecvR2(i)=mznew-2  	   		
	    			endif
	    		else
		    		mzrecvR1(i)=1
		    		mzrecvR2(i)=0
		    	endif
	    	else
!------------------			!------------------			!------------------	
!       | n+i |				!     | n+i  |				!     |n+i |	
!------------------   or	!------------------   or   	!------------------
!  | n |					!  | n   |					!  |    n     |	
!------------------	    	!------------------			!------------------				    			
	    		if(mznew-2+mznewcum .le. 3+mzcuml(row(proc+i))) then
	    			mzrecvR1(i)=1
	    			mzrecvR2(i)=0	
		    	
	    		else if(mznew-2+mznewcum .le. mzl(row(proc+i))-2+mzcuml(row(proc+i))) then		
	    			mzrecvR1(i)=3+mzcuml(row(proc+i))-mznewcum
	    			mzrecvR2(i)=mznew-2
	    		else		    		
	    			mzrecvR1(i)=3+mzcuml(row(proc+i))-mznewcum
	    			mzrecvR2(i)=mzl(row(proc+i))-2+mzcuml(row(proc+i))-mznewcum
	    		endif
	    	endif
	    	
	      endif !if(proc-i .lt. 1)
	      
	      countrecvR(i)=mx*my*(mzrecvR2(i)-mzrecvR1(i)+1)
	      
		
	   enddo	!do i=1, nleft  
	      
	   
		!***************************************************************		
		!send to its own processor
		!***************************************************************		
		if(3+mznewcum .gt. 3+mzcum) then
!------------------			!------------------			!------------------		
!   |  n   |				!   |  n   |			 	!   |  n  |	
!------------------   or    !------------------   or	!------------------
!  	    | n  |				!     | n|					!    		 | n |	
!------------------			!------------------			!------------------
			if(3+mznewcum .le. mz-2+mzcum) then
				if(mznew-2+mznewcum .gt. mz-2+mzcum) then
					mzsend1=3+mznewcum-mzcum
					mzsend2=mz-2
				else
					mzsend1=3+mznewcum-mzcum
					mzsend2=mznew-2+mznewcum-mzcum
				endif
			else
				mzsend1=1
				mzsend2=0
			endif
		else
!------------------			!------------------			!------------------	
!         | n  |			!     | n   |				!     |  n  |	
!------------------   or	!------------------   or   	!------------------
!  | n |					!  |  n  |					!  |    n      |	
!------------------	    	!------------------			!------------------				    			
	    	if(mznew-2+mznewcum .le. 3+mzcum) then
	    		mzsend1=1
	    		mzsend2=0	
	    			
	    	else if(mznew-2+mznewcum .le. mz-2+mzcum) then		
	    		mzsend1=3
	    		mzsend2=mznew-2+mznewcum-mzcum
	    	else		    		
	    		mzsend1=3
	    		mzsend2=mz-2
	    	endif
	    endif

	    
		!***************************************************************		
		!receive from its own processor
		!***************************************************************	    	  		  
    	if(3+mznewcum .gt. 3+mzcum) then	
!------------------			!------------------			!------------------		
!   |   n |					!   |  n  |			 		!   |  n  |	
!------------------   or    !------------------   or	!------------------
!  	   | n  |				!     |n|					!    		 | n |	
!------------------			!------------------			!------------------    	    	
	    	if(3+mznewcum .le. mz-2+mzcum) then	
	    		if(mznew-2+mznewcum .gt. mz-2+mzcum) then
 	    			mzrecv1=3  	      
	    			mzrecv2=mz-2+mzcum-mznewcum
	    		else
	    			mzrecv1=3
	    			mzrecv2=mznew-2  	   		
	    		endif
	    	else
		    	mzrecv1=1
		    	mzrecv2=0
		    endif
	    else
!------------------			!------------------			!------------------	
!       |  n  |				!     |  n  |				!     |  n  |	
!------------------   or	!------------------   or   	!------------------
!  | n|						!  | n   |					!  |     n    |	
!------------------	    	!------------------			!------------------				    			
	    	if(mznew-2+mznewcum .le. 3+mzcum) then
	    		mzrecv1=1
	    		mzrecv2=0	
		    	
	    	else if(mznew-2+mznewcum .le. mz-2+mzcum) then		
	    		mzrecv1=3+mzcum-mznewcum
	    		mzrecv2=mznew-2
	    	else		    		
	    		mzrecv1=3+mzcum-mznewcum
	    		mzrecv2=mz-2+mzcum-mznewcum
	    	endif
	    endif
	    	
	    
		!*******************************************************
		! Exchange B fields, send up and receive below
		!*******************************************************
		do i=1, nup

		  sendrank=modulo((rank/(sizex*sizey) + i),sizez)*(sizex*sizey) + & 
		 		  modulo(rank,sizex*sizey)
		  recvrank=modulo((rank/(sizex*sizey) - i),sizez)*(sizex*sizey) + & 
			      modulo(rank,sizex*sizey)
				
				
		  uptag=100

		  call MPI_SendRecv(bx(:,:,mzsendR1(i):mzsendR2(i)),countsendR(i),mpi_read,sendrank,uptag, &
				curx(:,:,mzrecvL1(i):mzrecvL2(i)),countrecvL(i),mpi_read,recvrank,uptag, &
				MPI_Comm_world,status,ierr)
		  call MPI_SendRecv(by(:,:,mzsendR1(i):mzsendR2(i)),countsendR(i),mpi_read,sendrank,uptag, &
				cury(:,:,mzrecvL1(i):mzrecvL2(i)),countrecvL(i),mpi_read,recvrank,uptag, &
				MPI_Comm_world,status,ierr)
		  call MPI_SendRecv(bz(:,:,mzsendR1(i):mzsendR2(i)),countsendR(i),mpi_read,sendrank,uptag, &
				curz(:,:,mzrecvL1(i):mzrecvL2(i)),countrecvL(i),mpi_read,recvrank,uptag, &
				MPI_Comm_world,status,ierr)
				
		enddo
		
		!*******************************************************
		! Exchange B fields, send down and receive above
		!*******************************************************
		do i=1, ndwn
 	  
		  sendrank=modulo((rank/(sizex*sizey) - i),sizez)*(sizex*sizey) + & 
			      modulo(rank,sizex*sizey)
		  recvrank=modulo((rank/(sizex*sizey) + i),sizez)*(sizex*sizey) + & 
		 		  modulo(rank,sizex*sizey) 
			  
		  dwntag=200
		  
		  call MPI_SendRecv(bx(:,:,mzsendL1(i):mzsendL2(i)),countsendL(i),mpi_read,sendrank,dwntag, &
				curx(:,:,mzrecvR1(i):mzrecvR2(i)),countrecvR(i),mpi_read,recvrank,dwntag, &
				MPI_Comm_world,status,ierr)
		  call MPI_SendRecv(by(:,:,mzsendL1(i):mzsendL2(i)),countsendL(i),mpi_read,sendrank,dwntag, &
				cury(:,:,mzrecvR1(i):mzrecvR2(i)),countrecvR(i),mpi_read,recvrank,dwntag, &
				MPI_Comm_world,status,ierr)
		  call MPI_SendRecv(bz(:,:,mzsendL1(i):mzsendL2(i)),countsendL(i),mpi_read,sendrank,dwntag, &
				curz(:,:,mzrecvR1(i):mzrecvR2(i)),countrecvR(i),mpi_read,recvrank,dwntag, &
				MPI_Comm_world,status,ierr)		  
		  
		enddo
	    
		!***************************************************
		! Exchange B fields, send to its own processor
		!***************************************************
		
	    if(mzrecv2 .ne. 0) then
	    	if(rank/(sizex*sizey) .eq. sizez-1) then
	    		curx(:,:,mzrecv1:mzrecv2+2)=bx(:,:,mzsend1:mzsend2+2)
	    		cury(:,:,mzrecv1:mzrecv2+2)=by(:,:,mzsend1:mzsend2+2)
	    		curz(:,:,mzrecv1:mzrecv2+2)=bz(:,:,mzsend1:mzsend2+2)	
	    	else if(rank/(sizex*sizey) .eq. 0) then
	    		curx(:,:,mzrecv1-2:mzrecv2)=bx(:,:,mzsend1-2:mzsend2)
	    		cury(:,:,mzrecv1-2:mzrecv2)=by(:,:,mzsend1-2:mzsend2)
	    		curz(:,:,mzrecv1-2:mzrecv2)=bz(:,:,mzsend1-2:mzsend2)
	    	else	    			
	    		curx(:,:,mzrecv1:mzrecv2)=bx(:,:,mzsend1:mzsend2)
	    		cury(:,:,mzrecv1:mzrecv2)=by(:,:,mzsend1:mzsend2)
	    		curz(:,:,mzrecv1:mzrecv2)=bz(:,:,mzsend1:mzsend2)
	    	endif
	    endif
		
		deallocate(bx,by,bz)
		
		allocate(bx(mx,my,mznew), by(mx,my,mznew), bz(mx,my,mznew))
		
		bx=curx
		by=cury
		bz=curz
	    	
		curx=0.
		cury=0.
		curz=0.
	    
		call MPI_BARRIER(MPI_COMM_WORLD,ierr)
			    
		
		!*******************************************************
		! Exchange E fields, send up and receive below
		!*******************************************************
		do i=1, nup

		  sendrank=modulo((rank/(sizex*sizey) + i),sizez)*(sizex*sizey) + & 
		 		  modulo(rank,sizex*sizey)
		  recvrank=modulo((rank/(sizex*sizey) - i),sizez)*(sizex*sizey) + & 
			      modulo(rank,sizex*sizey)
				
				
		  uptag=100

		  call MPI_SendRecv(ex(:,:,mzsendR1(i):mzsendR2(i)),countsendR(i),mpi_read,sendrank,uptag, &
				curx(:,:,mzrecvL1(i):mzrecvL2(i)),countrecvL(i),mpi_read,recvrank,uptag, &
				MPI_Comm_world,status,ierr)
		  call MPI_SendRecv(ey(:,:,mzsendR1(i):mzsendR2(i)),countsendR(i),mpi_read,sendrank,uptag, &
				cury(:,:,mzrecvL1(i):mzrecvL2(i)),countrecvL(i),mpi_read,recvrank,uptag, &
				MPI_Comm_world,status,ierr)
		  call MPI_SendRecv(ez(:,:,mzsendR1(i):mzsendR2(i)),countsendR(i),mpi_read,sendrank,uptag, &
				curz(:,:,mzrecvL1(i):mzrecvL2(i)),countrecvL(i),mpi_read,recvrank,uptag, &
				MPI_Comm_world,status,ierr)
				
		enddo
		
		!*******************************************************
		! Exchange E fields, send down and receive above
		!*******************************************************
		do i=1, ndwn
 	  
		  sendrank=modulo((rank/(sizex*sizey) - i),sizez)*(sizex*sizey) + & 
			      modulo(rank,sizex*sizey)
		  recvrank=modulo((rank/(sizex*sizey) + i),sizez)*(sizex*sizey) + & 
		 		  modulo(rank,sizex*sizey) 
			  
		  dwntag=200
		  
		  call MPI_SendRecv(ex(:,:,mzsendL1(i):mzsendL2(i)),countsendL(i),mpi_read,sendrank,dwntag, &
				curx(:,:,mzrecvR1(i):mzrecvR2(i)),countrecvR(i),mpi_read,recvrank,dwntag, &
				MPI_Comm_world,status,ierr)
		  call MPI_SendRecv(ey(:,:,mzsendL1(i):mzsendL2(i)),countsendL(i),mpi_read,sendrank,dwntag, &
				cury(:,:,mzrecvR1(i):mzrecvR2(i)),countrecvR(i),mpi_read,recvrank,dwntag, &
				MPI_Comm_world,status,ierr)
		  call MPI_SendRecv(ez(:,:,mzsendL1(i):mzsendL2(i)),countsendL(i),mpi_read,sendrank,dwntag, &
				curz(:,:,mzrecvR1(i):mzrecvR2(i)),countrecvR(i),mpi_read,recvrank,dwntag, &
				MPI_Comm_world,status,ierr)		  
		  
		enddo
	    
		!***************************************************
		! Exchange E fields, send to its own processor
		!***************************************************
		
	    if(mzrecv2 .ne. 0) then
	    	if(rank/(sizex*sizey) .eq. sizez-1) then
	    		curx(:,:,mzrecv1:mzrecv2+2)=ex(:,:,mzsend1:mzsend2+2)
	    		cury(:,:,mzrecv1:mzrecv2+2)=ey(:,:,mzsend1:mzsend2+2)
	    		curz(:,:,mzrecv1:mzrecv2+2)=ez(:,:,mzsend1:mzsend2+2)	
	    	else if(rank/(sizex*sizey) .eq. 0) then
	    		curx(:,:,mzrecv1-2:mzrecv2)=ex(:,:,mzsend1-2:mzsend2)
	    		cury(:,:,mzrecv1-2:mzrecv2)=ey(:,:,mzsend1-2:mzsend2)
	    		curz(:,:,mzrecv1-2:mzrecv2)=ez(:,:,mzsend1-2:mzsend2)
	    	else	    			
	    		curx(:,:,mzrecv1:mzrecv2)=ex(:,:,mzsend1:mzsend2)
	    		cury(:,:,mzrecv1:mzrecv2)=ey(:,:,mzsend1:mzsend2)
	    		curz(:,:,mzrecv1:mzrecv2)=ez(:,:,mzsend1:mzsend2)
	    	endif
	    endif
		
		deallocate(ex,ey,ez)
		
		allocate(ex(mx,my,mznew), ey(mx,my,mznew), ez(mx,my,mznew))
		
		ex=curx
		ey=cury
		ez=curz
		
		curx=0.
		cury=0.
		curz=0.

		call MPI_BARRIER(MPI_COMM_WORLD,ierr)		
			   	   
		! buffer size for particle exchange
	    ! pinlft and pinrgt have larger buffer size than poutrgt and poulgt
	    ! becuause one processor can receive particles from several cpus
		buffsize1=maxval(totpartl)*2
		buffsize0=maxval(totpartl)*4
			
		deallocate(poutup,poutdwn,pinblw,pinabv)	
		
		allocate(poutup(buffsize1),poutdwn(buffsize1))
		allocate(pinblw(buffsize0),pinabv(buffsize0))
	
		
!		do n=maxhlf+1,maxhlf+lecs
!				n1=n
!				if(p(n1)%x<0 .or. p(n1)%y<0 .or. p(n1)%z<0 .or. &
!				isnan(p(n1)%x*p(n)%y*p(n1)%z) .or. p(n1)%x>mx-1.or.p(n1)%y>my-1.or. &
!				p(n1)%z>mznew-1) then !.or.zt>mz-1) then
!				print*, 'found a ghost particle before redistribution....'
!				print*, "prb in lec dep, ind, proc, lap, rank: "&
!				,p(n1)%ind," ",p(n1)%proc," ",&
!				lap," ",rank," x=",p(n1)%x," y=",p(n1)%y," z=",p(n1)%z,&
!					"mx,y,z ",mx,my,mz
!               endif
!		enddo		
				
		!****************************************************************
		! send ions to the up and recive from the down
		!****************************************************************	
		lenioninblw0=0
		
		do i=1, nup

			LenIonOutUp=0

			sendrank=modulo((rank/(sizex*sizey) + i),sizez)*(sizex*sizey) + & 
		 		  modulo(rank,sizex*sizey) 
			recvrank=modulo((rank/(sizex*sizey) - i),sizez)*(sizex*sizey) + & 
			      modulo(rank,sizex*sizey) 

			if(countsendR(i) .ne. 0) then			
			n=1
			if (ions.gt.0) then
		52  	continue
				n1=n			
				if((p(n1)%z .ge. mzsendR1(i)) .and. (p(n1)%z .le. mzsendR2(i))) then
					in=.false.
					p(n1)%z=p(n1)%z+(mzcum-mznewcuml(sendrank+1))
				else
					in=.true.
				endif
							
				if(in) go to 58    	
				LenIonOutUp=LenIonOutUp+1
				call copyprt(p(n1),poutup(LenIonOutUp))
				call copyprt(p(ions),p(n1))		
				ions=ions-1
				n=n-1
				
		58		n=n+1
	
				if(n.le.ions) go to 52
			endif !if (ions.gt.0) then
			
			endif !if(countsendR(i) .ne. 0) then
			
			uptag=100
			
			call MPI_SendRecv(LenIonOutUp,1,mpi_integer,sendrank,uptag, &
						LenIonInBlw,1,mpi_integer,recvrank,uptag, &
						MPI_Comm_World,status,ierr)
								
			lenoutup=max(LenIonOutUp,1)
			leninblw=max(LenIonInBlw,1)
			
			if(lenioninblw0+leninblw .gt. buffsize0) then
					print *, "after z redistribution: BUFFERSIZE SMALLER THAN PARTICLE LIST"
					print *,rank,":","buffsize", buffsize1, "lenioninblw",lenioninblw0+leninblw
					stop
			endif
			
			call MPI_SendRecv(poutup(1:lenoutup),lenoutup,particletype,sendrank,uptag, &
				pinblw(lenioninblw0+1:lenioninblw0+leninblw),leninblw,particletype,recvrank,uptag, &
				MPI_Comm_World,status,ierr)
						
			lenioninblw0=lenioninblw0+LenIonInBlw	
			
			if(debug) print*, 'rank, ions sent to the i''th right, ',rank, i
			
			call MPI_BARRIER(MPI_COMM_WORLD,ierr)

				
		enddo !do i=1, nright	
		
		!****************************************************************
		! send ions to the down and recive from the up
		!****************************************************************	
		lenioninabv0=0
		
		do i=1, ndwn

			LenIonOutDwn=0

			sendrank=modulo(rank/(sizex*sizey) - i,sizez)*(sizex*sizey) + & 
			      modulo(rank,sizex*sizey)
			recvrank=modulo(rank/(sizex*sizey) + i,sizez)*(sizex*sizey) + & 
		 		  modulo(rank,sizex*sizey) 
		
			if(countsendL(i) .ne. 0) then			
			n=1
			if (ions.gt.0) then
		62  	continue
				n1=n			
				if((p(n1)%z .ge. mzsendL1(i)) .and. (p(n1)%z .le. mzsendL2(i))) then
					in=.false.
					p(n1)%z=p(n1)%z+(mzcum-mznewcuml(sendrank+1))
				else
					in=.true.
				endif
							
				if(in) go to 68    	
				LenIonOutDwn=LenIonOutDwn+1
				call copyprt(p(n1),poutdwn(LenIonOutDwn))
				call copyprt(p(ions),p(n1))		
				ions=ions-1
				n=n-1
				
		68		n=n+1
	
				if(n.le.ions) go to 62
			endif !if (ions.gt.0) then
			
			endif !if(countsendR(i) .ne. 0) then
			
			dwntag=200
			
			call MPI_SendRecv(LenIonOutDwn,1,mpi_integer,sendrank,dwntag, &
						LenIonInAbv,1,mpi_integer,recvrank,dwntag, &
						MPI_Comm_World,status,ierr)
								
			lenoutdwn=max(LenIonOutDwn,1)
			leninabv=max(LenIonInAbv,1)
			
			if(lenioninabv0+leninabv .gt. buffsize0) then
					print *, "after z redistribution: BUFFERSIZE SMALLER THAN PARTICLE LIST"
					print *,rank,":","buffsize", buffsize0, "lenioninabv",lenioninabv0+leninabv
					stop
			endif			
			
			call MPI_SendRecv(poutdwn(1:lenoutdwn),lenoutdwn,particletype,sendrank,dwntag, &
				pinabv(lenioninabv0+1:lenioninabv0+leninabv),leninabv,particletype,recvrank,dwntag, &
				MPI_Comm_World,status,ierr)
						
			lenioninabv0=lenioninabv0+LenIonInAbv

			if(debug) print*, 'ions sent to the i''th left',i
			
			call MPI_BARRIER(MPI_COMM_WORLD,ierr)

					
	enddo !do i=1, nleft	
		
		! shift particle positons in a new domain
		p(1:ions)%z=p(1:ions)%z+(mzcum-mznewcum)	
	
		! transfer particles in a new domain
		do n=1,lenioninabv0
			ions=ions+1
			call copyprt(pinabv(n),p(ions))
		enddo		
				
		do n=1,lenioninblw0
			ions=ions+1
			call copyprt(pinblw(n),p(ions))
		enddo

		
		!****************************************************************
		! send electrons to the up and recive from the down
		!****************************************************************	
		lenlecinblw0=0
		
		do i=1, nup

			LenLecOutUp=0

			sendrank=modulo((rank/(sizex*sizey) + i),sizez)*(sizex*sizey) + & 
		 		  modulo(rank,sizex*sizey) 
			recvrank=modulo((rank/(sizex*sizey) - i),sizez)*(sizex*sizey) + & 
			      modulo(rank,sizex*sizey) 
			
			if(countsendR(i) .ne. 0) then			
			n=maxhlf+1
			if (lecs.gt.0) then
		53 	continue
				n1=n
				if((p(n1)%z .ge. mzsendR1(i)) .and. (p(n1)%z .le. mzsendR2(i))) then
					in=.false.
					p(n1)%z=p(n1)%z+(mzcum-mznewcuml(sendrank+1))
				else
					in=.true.
				endif
							
				if(in) go to 59    	
				LenLecOutUp=LenLecOutUp+1
				call copyprt(p(n1),poutup(LenLecOutUp))
				call copyprt(p(maxhlf+lecs),p(n1))		
				lecs=lecs-1
				n=n-1
				
		59		n=n+1
	
				if(n.le.maxhlf+lecs) go to 53
			endif !if (lecs.gt.0) then
			
			endif !if(countsendR(i) .ne. 0) then
			
			uptag=100
						
			call MPI_SendRecv(LenLecOutUp,1,mpi_integer,sendrank,uptag, &
						LenLecInBlw,1,mpi_integer,recvrank,uptag, &
						MPI_Comm_World,status,ierr)
								
			lenoutup=max(LenLecOutUp,1)
			leninblw=max(LenLecInBlw,1)
			
			if(lenlecinblw0+leninblw .gt. buffsize0) then
					print *, "after z redistribution: BUFFERSIZE SMALLER THAN PARTICLE LIST"
					print *,rank,":","buffsize", buffsize0, "lenlecinblw",lenioninblw0+leninblw
					stop
			endif
		
			
			call MPI_SendRecv(poutup(1:lenoutup),lenoutup,particletype,sendrank,uptag, &
				pinblw(lenlecinblw0+1:lenlecinblw0+leninblw),leninblw,particletype,recvrank,uptag, &
				MPI_Comm_World,status,ierr)
						
			lenlecinblw0=lenlecinblw0+LenLecInBlw

			if(debug) print*, 'rank, electrons sent to the i''th right',rank, i

			call MPI_BARRIER(MPI_COMM_WORLD,ierr)
			
			
		enddo !do i=1, nup		
		
		
		!****************************************************************
		! send electrons to the down and recive from the up
		!****************************************************************	
		lenlecinabv0=0
		
		do i=1, ndwn

			LenLecOutDwn=0
                                            
			sendrank=modulo(rank/(sizex*sizey) - i,sizez)*(sizex*sizey) + & 
			      modulo(rank,sizex*sizey)
			recvrank=modulo(rank/(sizex*sizey) + i,sizez)*(sizex*sizey) + & 
		 		  modulo(rank,sizex*sizey) 
		
			if(countsendL(i) .ne. 0) then			
			n=maxhlf+1
			if (lecs.gt.0) then
		63 	continue
				n1=n
				if((p(n1)%z .ge. mzsendL1(i)) .and. (p(n1)%z .le. mzsendL2(i))) then
					in=.false.
					p(n1)%z=p(n1)%z+(mzcum-mznewcuml(sendrank+1))
				else
					in=.true.
				endif
						
				if(in) go to 69    	
				LenLecOutDwn=LenLecOutDwn+1
				call copyprt(p(n1),poutdwn(LenLecOutDwn))
				call copyprt(p(maxhlf+lecs),p(n1))		
				lecs=lecs-1
				n=n-1
				
		69		n=n+1
	
				if(n.le.maxhlf+lecs) go to 63
			endif !if (lecs.gt.0) then
			
			endif !if(countsendR(i) .ne. 0) then
			
			dwntag=200
			
			call MPI_SendRecv(LenLecOutDwn,1,mpi_integer,sendrank,dwntag, &
						LenLecInAbv,1,mpi_integer,recvrank,dwntag, &
						MPI_Comm_World,status,ierr)
								
			lenoutdwn=max(LenLecOutDwn,1)
			leninabv=max(LenLecInAbv,1)
			
			if(lenlecinabv0+leninabv .gt. buffsize0) then
					print *, "After z redistribution : BUFFERSIZE SMALLER THAN PARTICLE LIST"
					print *,rank,":","buffsize", buffsize1, "lenlecinabv",lenlecinabv0+leninabv
					stop
			endif			
			
			call MPI_SendRecv(poutdwn(1:lenoutdwn),lenoutdwn,particletype,sendrank,dwntag, &
				pinabv(lenlecinabv0+1:lenlecinabv0+leninabv),leninabv,particletype,recvrank,dwntag, &
				MPI_Comm_World,status,ierr)
						
			lenlecinabv0=lenlecinabv0+LenLecInAbv			
			
			if(debug)  print*, 'rank, electrons sent to the i''th down', rank,i	
			
			call MPI_BARRIER(MPI_COMM_WORLD,ierr)
	
		enddo !do i=1, ndwn	
		
		! shift particle positons in a new domain
!		p(maxhlf+1:maxhlf+lecs)%z=p(maxhlf+1:maxhlf+lecs)%z+(mzcum-mznewcum)
		! transfer particles in a new domain

		p(maxhlf+1:maxhlf+lecs)%z=p(maxhlf+1:maxhlf+lecs)%z+(mzcum-mznewcum)
								
		do n=1,lenlecinblw0
			lecs=lecs+1
			call copyprt(pinblw(n),p(maxhlf+lecs))
		enddo
			
		do n=1,lenlecinabv0
			lecs=lecs+1
			call copyprt(pinabv(n),p(maxhlf+lecs))
		enddo	
				
		mz=mznew
		mzall=mz
			
		call mpi_allgather(mz,1,mpi_integer,mzl,1,mpi_integer &
		  	   	  ,mpi_comm_world, error)
			
		k1=1
		k2=rank/(sizex*sizey)*(sizex*sizey)+1
		mzcum=sum(mzl(k1:k2:(sizex*sizey))-5)-(mzl(k2)-5)
	
		ix=1
		iy=mx
		iz=iy*my
		lot=iz*mz
		
#ifdef twoD
		iz=0
		lot=mx*my
#endif			
				
		totpart=ions+lecs		! total particles in a local cpu
		call mpi_allgather(totpart,1,mpi_integer8,totpartl,1,mpi_integer8 &
     	   	  ,mpi_comm_world, error)
 
     	deallocate(mzsendR1,mzsendR2,mzrecvL1,mzrecvL2)
        deallocate(mzsendL1,mzsendL2,mzrecvR1,mzrecvR2)
        deallocate(countsendR,countrecvL,countsendL,countrecvR)
  			
		deallocate(mznewl,mzcuml,mznewcuml)
		deallocate(totpartl,totpartcuml)
		deallocate(row)	
		
		deallocate(bufferin1,bufferin2)
		deallocate(bufferin1y,bufferin2y,bufferin1x,bufferin2x)
		deallocate(sendbufy,sendbufz)
		deallocate(temp,poutup,poutdwn,pinblw &
			, pinabv,poutlft,poutrgt,pinlft, pinrgt &
			, poutplus,poutminus,pinplus,pinminus,pall )
#ifdef filter2
		deallocate(xghost,yghost,zghost) !for filter2	
	
#ifndef twoD		
	call free_MPI_filter_datatypes()
#endif

#ifndef twoD
	call create_MPI_filter_datatypes(3,mx-3,3,my-3,3,mz-3)
#endif
!#else
!	call create_MPI_filter_datatypes(3,mx-3,3,my-3,1,1)
!#endif
#endif !if filter2
				
#ifndef twoD
		buffsize1 = max(10*int(ppc0*c*max(1*(mx-5)*(my-5),1*(mx-5)*(mz-5))) &
				,10000)					
#else
		buffsize1 = max(3*int(ppc0*c*max((mx-5),(my-5))),60000)
#endif
		allocate(bufferin1(mx,my,2),bufferin2(mx,my,2), &
			bufferin1y(mx,2,mz),bufferin2y(mx,2,mz))
	
		allocate(bufferin1x(2,my,mz),bufferin2x(2,my,mz))
		allocate(sendbufy(mx,2,mz),sendbufz(mx,my,2))
		allocate(temp(mx,my,mz))
#ifdef filter2
		allocate(xghost(2*ntimes,my,mz), yghost(mx,2*ntimes,mz), zghost(mx,my,2*ntimes)) !for filter2
#endif
		allocate(poutup(buffsize1),poutdwn(buffsize1))
		allocate(pinabv(buffsize1),pinblw(buffsize1))
		allocate(poutrgt(buffsize1),poutlft(buffsize1))
		allocate(pinlft(buffsize1),pinrgt(buffsize1))
		allocate(poutplus(buffsize1),poutminus(buffsize1))
		allocate(pinminus(buffsize1),pinplus(buffsize1))
		allocate(pall(lot))
		
		
		pall=0

		bufferin1=0.
		bufferin2=0.
		bufferin1y=0.
		bufferin2y=0.
		bufferin1x=0.
		bufferin2x=0.
		sendbufy=0.
		sendbufz=0.
		temp=0.
		
		call bc_b1()
		if(debug) print*, 'bc done', rank
		call bc_e1()
				

		if(rank .eq. 0) print *,"domain z redistribution done  "			
		
		call MPI_BARRIER(MPI_COMM_WORLD,ierr)
				
	
end subroutine redist_z_domain


#ifdef twoD
end module m_dynamic_domain
#else
end module m_dynamic_domain_3d
#endif

