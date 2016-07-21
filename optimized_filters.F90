!-------------------------------------------------------------------------------
! 				subroutine apply_filter2_opt()		
!										
! Improved optimization using local temporary variables and eliminating 
! temporary arrays. Only works if all dimensions >= ntimes
!
!-------------------------------------------------------------------------------

subroutine apply_filter2_opt()
	implicit none
	
	! local variables
	
	real wtm1,wt,wtp1,winv
! 	real, allocatable, dimension(:,:,:) :: xghost, yghost, zghost
	integer :: istr, ifin, jstr, jfin, kstr, kfin,istr_filter, ifin_filter
	integer :: i, j, k, n, ierr
	

	! istr,ifin are the boundaries for periodicity
	! istr_filter,ifin_filter are the filtering boundaries
	! i
	istr=3
	ifin=mx-3
	istr_filter=istr
	ifin_filter=ifin
!	if(wall) then 
!		istr_filter=3		!right going injection
!		ifin_filter=min(mx-3,int(xinject2)+10)
!	endif

! 	print *,"Wall: ",wall
! 	print *,"istr:",istr,"; ifin",ifin

! 	print *,"mx",mx,"my",my,"mz",mz

	! j
	jstr=3
	jfin=my-3

	! k
#ifndef twoD
	kstr=3
	kfin=mz-3 !2,mz-1 !3,mz-3
#else
	kstr=1
	kfin=1
#endif

! 	allocate(xghost(2*ntimes,my,mz), yghost(mx,2*ntimes,mz), zghost(mx,my,2*ntimes))
! 	xghost=0
! 	yghost=0
! 	zghost=0
	
	! MPI Datatypes
! 	if(lap.eq.lapst) then
! 		call create_MPI_filter_datatypes(istr,ifin,jstr,jfin,kstr,kfin)
! 	endif
	call timer(20)
	! filter curx in x direction
	!!! Populate xghost with curx data !!!
	if(periodicx.eq.1) then
		!call deep_copylayrx(curx,xghost,istr,ifin,jstr,jfin,kstr,kfin)
		call deep_copy_layrx1(curx,xghost,istr,ifin,jstr,jfin,kstr,kfin)
	else
		!call nonper_copylayrx(curx,xghost,istr,ifin,jstr,jfin,kstr,kfin)
		call deep_copy_layrx2(curx,xghost,istr,ifin,jstr,jfin,kstr,kfin)
	endif
	call timer(20,tmstop=.true.)

	call timer(21)
	call filter_x(curx,xghost,istr_filter,ifin_filter,jstr,jfin,kstr,kfin)
	call timer(21,tmstop=.true.)

	call timer(22)
	!filter curx in y direction
	!!! Populate yghost with curx data !!!
	if(periodicy.eq.1) then
		call deep_copy_layry1(curx,yghost,istr,ifin,jstr,jfin,kstr,kfin)
	else
		call deep_copy_layry2(curx,yghost,istr,ifin,jstr,jfin,kstr,kfin)
	endif
	call timer(22,tmstop=.true.)

	call timer(23)
	call filter_y(curx,yghost,istr_filter,ifin_filter,jstr,jfin,kstr,kfin)
	call timer(23,tmstop=.true.)

	call timer(24)
	!filter curx in z direction
#ifndef twoD
	!!! Populate zghost with curx data !!!
	if(periodicz.eq.1) then
		call deep_copy_layrz1(curx,zghost,istr,ifin,jstr,jfin,kstr,kfin)
	else
		!transfer z information among processors, not periodic between 0,n-1
		call deep_copy_layrz2(curx,zghost,istr,ifin,jstr,jfin,kstr,kfin)
	endif
	call timer(24,tmstop=.true.)

	call timer(25)
	call filter_z(curx,zghost,istr_filter,ifin_filter,jstr,jfin,kstr,kfin)
	call timer(25,tmstop=.true.)
#endif
	
	call timer(26)
	! filter cury in x direction
	!!! Populate xghost with cury data !!!
	if(periodicx.eq.1) then
		call deep_copylayrx(cury,xghost,istr,ifin,jstr,jfin,kstr,kfin)
	else
		call nonper_copylayrx(cury,xghost,istr,ifin,jstr,jfin,kstr,kfin)
	endif
	call timer(26,tmstop=.true.)
		
	call timer(27)
	call filter_x(cury,xghost,istr_filter,ifin_filter,jstr,jfin,kstr,kfin)
	call timer(27,tmstop=.true.)
	
	call timer(28)
	!filter cury in y direction
	!!! Populate yghost with cury data !!!
	if(periodicy.eq.1) then
		call deep_copy_layry1(cury,yghost,istr,ifin,jstr,jfin,kstr,kfin)
	else
		call deep_copy_layry2(cury,yghost,istr,ifin,jstr,jfin,kstr,kfin)
	endif
	call timer(28,tmstop=.true.)
	
	call timer(29)
	call filter_y(cury,yghost,istr_filter,ifin_filter,jstr,jfin,kstr,kfin)
	call timer(29,tmstop=.true.)
	
	call timer(30)
	!filter cury in z direction
#ifndef twoD
	!!! Populate zghost with cury data !!!
	if(periodicz.eq.1) then
		call deep_copy_layrz1(cury,zghost,istr,ifin,jstr,jfin,kstr,kfin)
	else
		!transfer z information among processors, not periodic between two ends of the grid
		call deep_copy_layrz2(cury,zghost,istr,ifin,jstr,jfin,kstr,kfin)
	endif
	call timer(30,tmstop=.true.)

	call timer(31)
	call filter_z(cury,zghost,istr_filter,ifin_filter,jstr,jfin,kstr,kfin)
	call timer(31,tmstop=.true.)
#endif
	
	call timer(32)
	! filter curz in x direction
	!!! Populate xghost with curz data !!!
	if(periodicx.eq.1) then
		call deep_copylayrx(curz,xghost,istr,ifin,jstr,jfin,kstr,kfin)
	else
		call nonper_copylayrx(curz,xghost,istr,ifin,jstr,jfin,kstr,kfin)
	endif
	call timer(32,tmstop=.true.)
	
	call timer(33)
	call filter_x(curz,xghost,istr_filter,ifin_filter,jstr,jfin,kstr,kfin)
	call timer(33,tmstop=.true.)
	
	call timer(34)
	!filter curz in y direction
	!!! Populate yghost with curz data !!!
	if(periodicy.eq.1) then
		call deep_copy_layry1(curz,yghost,istr,ifin,jstr,jfin,kstr,kfin)
	else
		call deep_copy_layry2(curz,yghost,istr,ifin,jstr,jfin,kstr,kfin)
	endif
	call timer(34,tmstop=.true.)

	call timer(35)
	call filter_y(curz,yghost,istr_filter,ifin_filter,jstr,jfin,kstr,kfin)
	call timer(35,tmstop=.true.)
	
	call timer(36)
	!filter curz in z direction
#ifndef twoD
	!!! Populate zghost with curz data !!!
	if(periodicz.eq.1) then
		call deep_copy_layrz1(curz,zghost,istr,ifin,jstr,jfin,kstr,kfin)
	else
		!transfer z information among processors, not periodic between two ends of the grid
		call deep_copy_layrz2(curz,zghost,istr,ifin,jstr,jfin,kstr,kfin)
	endif
	call timer(36,tmstop=.true.)

	call timer(37)
	call filter_z(curz,zghost,istr_filter,ifin_filter,jstr,jfin,kstr,kfin)
	call timer(37,tmstop=.true.)

#endif

! 	deallocate(xghost)
! 	deallocate(yghost)
! 	deallocate(zghost)

! 	print *,"Ending filter"

! 	do i=20,37
! 		print *,i,real(tm(i),4)
! 	enddo
	
! 	if(rank == 0) then
! 		print *,"curx"
! 		print *,"filter_x:",real(tm(21),4)
! 		print *,"filter_y:",real(tm(23),4)
! 		print *,"filter_z:",real(tm(25),4)
! 	
! 		print *,"cury"
! 		print *,"filter_x:",real(tm(27),4)
! 		print *,"filter_y:",real(tm(29),4)
! 		print *,"filter_z:",real(tm(31),4)
! 	
! 		print *,"curz"
! 		print *,"filter_x:",real(tm(33),4)
! 		print *,"filter_y:",real(tm(35),4)
! 		print *,"filter_z:",real(tm(37),4)
! 
! ! 		print *,"nocopy_filter2:",real(tm(21)+tm(23)+tm(25)+tm(27)+tm(29)+tm(31)+tm(33)+tm(35)+tm(37),4)
! 	endif


end subroutine apply_filter2_opt

!-------------------------------------------------------------------------------
! 				subroutine apply_filter2_opt()		
!										
! Same as apply_filter2_opt(), except copies before filtering
!
!-------------------------------------------------------------------------------
subroutine apply_filter2_copy_opt()
	implicit none
	
	! local variables
	
	real wtm1,wt,wtp1,winv
! 	real, allocatable, dimension(:,:,:) :: xghost, yghost, zghost
	integer :: istr, ifin, jstr, jfin, kstr, kfin,istr_filter, ifin_filter
	integer :: i, j, k, n, ierr
	
	! istr,ifin are the boundaries for periodicity
	! istr_filter,ifin_filter are the filtering boundaries
	! i
	istr=3
	ifin=mx-3
	istr_filter=istr
	ifin_filter=ifin
! ! ! 	if(wall) then 
! ! ! 	   istr_filter=2		!right going injection
! ! ! 	   ifin_filter=min(mx-1,int(xinject2)+10)
! ! ! 	endif

! 	print *,"Wall: ",wall
! 	print *,"istr:",istr,"; ifin",ifin

! 	print *,"mx",mx,"my",my,"mz",mz

	! j
	jstr=3
	jfin=my-3

	! k
#ifndef twoD
	kstr=3
	kfin=mz-3 !2,mz-1 !3,mz-3
#else
	kstr=1
	kfin=1
#endif

! 	allocate(xghost(2*ntimes,my,mz), yghost(mx,2*ntimes,mz), zghost(mx,my,2*ntimes))
! 	xghost=0
! 	yghost=0
! 	zghost=0
	
	! MPI Datatypes
! 	if(lap.eq.lapst) then
! 		call create_MPI_filter_datatypes(istr,ifin,jstr,jfin,kstr,kfin)
! 	endif
	call timer(20)
	! filter curx in x direction
	!!! Populate xghost with curx data !!!
	if(periodicx.eq.1) then
		call deep_copylayrx(curx,xghost,istr,ifin,jstr,jfin,kstr,kfin)
	else
		call nonper_copylayrx(curx,xghost,istr,ifin,jstr,jfin,kstr,kfin)
	endif
	call timer(20,tmstop=.true.)

	call timer(121)
	call copy_filter_x(curx,xghost,istr_filter,ifin_filter,jstr,jfin,kstr,kfin)
	call timer(121,tmstop=.true.)

	call timer(22)
	!filter curx in y direction
	!!! Populate yghost with curx data !!!
	if(periodicy.eq.1) then
		call deep_copy_layry1(curx,yghost,istr,ifin,jstr,jfin,kstr,kfin)
	else
		print *, "Radiation BC not implemented in apply_filter"
	endif
	call timer(22,tmstop=.true.)

	call timer(123)
	call copy_filter_y(curx,yghost,istr_filter,ifin_filter,jstr,jfin,kstr,kfin)
	call timer(123,tmstop=.true.)

	call timer(24)
	!filter curx in z direction
#ifndef twoD
	!!! Populate zghost with curx data !!!
	if(periodicz.eq.1) then
		call deep_copy_layrz1(curx,zghost,istr,ifin,jstr,jfin,kstr,kfin)
	else
		!transfer z information among processors, not periodic between 0,n-1
		call deep_copy_layrz2(curx,zghost,istr,ifin,jstr,jfin,kstr,kfin)
	endif
	call timer(24,tmstop=.true.)

	call timer(125)
	call copy_filter_z(curx,zghost,istr_filter,ifin_filter,jstr,jfin,kstr,kfin)
	call timer(125,tmstop=.true.)
#endif
	
	call timer(26)
	! filter cury in x direction
	!!! Populate xghost with cury data !!!
	if(periodicx.eq.1) then
		call deep_copylayrx(cury,xghost,istr,ifin,jstr,jfin,kstr,kfin)
	else
		call nonper_copylayrx(cury,xghost,istr,ifin,jstr,jfin,kstr,kfin)
	endif
	call timer(26,tmstop=.true.)
		
	call timer(127)
	call copy_filter_x(cury,xghost,istr_filter,ifin_filter,jstr,jfin,kstr,kfin)
	call timer(127,tmstop=.true.)
	
	call timer(28)
	!filter cury in y direction
	!!! Populate yghost with cury data !!!
	if(periodicy.eq.1) then
		call deep_copy_layry1(cury,yghost,istr,ifin,jstr,jfin,kstr,kfin)
	else
		print *, "Radiation BC not implemented in apply_filter"
	endif
	call timer(28,tmstop=.true.)
	
	call timer(129)
	call copy_filter_y(cury,yghost,istr_filter,ifin_filter,jstr,jfin,kstr,kfin)
	call timer(129,tmstop=.true.)
	
	call timer(30)
	!filter cury in z direction
#ifndef twoD
	!!! Populate zghost with cury data !!!
	if(periodicz.eq.1) then
		call deep_copy_layrz1(cury,zghost,istr,ifin,jstr,jfin,kstr,kfin)
	else
		!transfer z information among processors, not periodic between 0,n-1
		call deep_copy_layrz2(cury,zghost,istr,ifin,jstr,jfin,kstr,kfin)
	endif
	call timer(30,tmstop=.true.)

	call timer(131)
	call copy_filter_z(cury,zghost,istr_filter,ifin_filter,jstr,jfin,kstr,kfin)
	call timer(131,tmstop=.true.)
#endif
	
	call timer(32)
	! filter curz in x direction
	!!! Populate xghost with curz data !!!
	if(periodicx.eq.1) then
		call deep_copylayrx(curz,xghost,istr,ifin,jstr,jfin,kstr,kfin)
	else
		call nonper_copylayrx(curz,xghost,istr,ifin,jstr,jfin,kstr,kfin)
	endif
	call timer(32,tmstop=.true.)
	
	call timer(133)
	call copy_filter_x(curz,xghost,istr_filter,ifin_filter,jstr,jfin,kstr,kfin)
	call timer(133,tmstop=.true.)
	
	call timer(34)
	!filter curz in y direction
	!!! Populate yghost with curz data !!!
	if(periodicy.eq.1) then
		call deep_copy_layry1(curz,yghost,istr,ifin,jstr,jfin,kstr,kfin)
	else
		print *, "Radiation BC not implemented in apply_filter"
	endif
	call timer(34,tmstop=.true.)

	call timer(135)
	call copy_filter_y(curz,yghost,istr_filter,ifin_filter,jstr,jfin,kstr,kfin)
	call timer(135,tmstop=.true.)
	
	call timer(36)
	!filter curz in z direction
#ifndef twoD
	!!! Populate zghost with curz data !!!
	if(periodicz.eq.1) then
		call deep_copy_layrz1(curz,zghost,istr,ifin,jstr,jfin,kstr,kfin)
	else
		!transfer z information among processors, not periodic between 0,n-1
		call deep_copy_layrz2(curz,zghost,istr,ifin,jstr,jfin,kstr,kfin)
	endif
	call timer(36,tmstop=.true.)

	call timer(137)
	call copy_filter_z(curz,zghost,istr_filter,ifin_filter,jstr,jfin,kstr,kfin)
	call timer(137,tmstop=.true.)

#endif

! 	deallocate(xghost)
! 	deallocate(yghost)
! 	deallocate(zghost)

! 	print *,"Ending filter"

! 	do i=20,37
! 		print *,i,real(tm(i),4)
! 	enddo
	
! 	if(rank == 0) then
! 		print *,"curx"
! 		print *,"copy_filter_x:",real(tm(121),4)
! 		print *,"copy_filter_y:",real(tm(123),4)
! 		print *,"copy_filter_z:",real(tm(125),4)
! 	
! 		print *,"cury"
! 		print *,"copy_filter_x:",real(tm(127),4)
! 		print *,"copy_filter_y:",real(tm(129),4)
! 		print *,"copy_filter_z:",real(tm(131),4)
! 	
! 		print *,"curz"
! 		print *,"copy_filter_x:",real(tm(133),4)
! 		print *,"copy_filter_y:",real(tm(135),4)
! 		print *,"copy_filter_z:",real(tm(137),4)
! 
! ! 		print *,"copy_filter2:",real(tm(121)+tm(123)+tm(125)+tm(127)+tm(129)+tm(131)+tm(133)+tm(135)+tm(137),4)
! 	endif

end subroutine apply_filter2_copy_opt

!-------------------------------------------------------------------------------
! 				subroutine filter_x()		
!										
! 	Does "ntimes" passes in the x direction on the specified array, using the 
!	specified ghost cells
!
!  Assumes at least 3 real cells.  Assumes nothing about ntimes
!-------------------------------------------------------------------------------
subroutine filter_x(cur, ghost, istr, ifin, jstr, jfin, kstr, kfin)
	implicit none
	real, dimension(:,:,:) :: cur, ghost
	real :: temp1, temp2
	real :: wtm1,wt,wtp1
	integer :: i,j,k,n
	integer :: istr, ifin, jstr, jfin, kstr, kfin
	logical :: even, ntimes_even

	wtm1 = .25
	wt = .5
	wtp1 = .25

	ntimes_even = .false.
	if(mod(ntimes,2) .eq. 0) then
		ntimes_even = .true.
	endif

	even = .false.
	if(mod(ifin-istr+1,2) .eq. 0) then
		even = .true.
	endif

	do k=kstr,kfin
		do j=jstr,jfin
			n=1
			do while(n .le. ntimes)
				! Filter the ghost cells preceeding the cur array
				temp2 = ghost(1,j,k)
				do i=2,ntimes-2,2
					temp1=wtm1*ghost(i-1,j,k)+wt*ghost(i,j,k)+wtp1*ghost(i+1,j,k)
					ghost(i-1,j,k)=temp2
					temp2=wtm1*ghost(i,j,k)+wt*ghost(i+1,j,k)+wtp1*ghost(i+2,j,k)
					ghost(i,j,k)=temp1
				enddo
				if(ntimes .gt. 1) then
					if(.not.ntimes_even) then
						! i = ntimes-1
						temp1=wtm1*ghost(i-1,j,k)+wt*ghost(i,j,k)+wtp1*ghost(i+1,j,k)
						ghost(i-1,j,k)=temp2
						temp2=temp1
					endif
					i = ntimes
					temp1=wtm1*ghost(i-1,j,k)+wt*ghost(i,j,k)+wtp1*cur(istr,j,k)
					ghost(i-1,j,k)=temp2
				else
					i=1
					temp1=temp2
				endif
				temp2=wtm1*ghost(i,j,k)+wt*cur(istr,j,k)+wtp1*cur(istr+1,j,k)
				ghost(i,j,k)=temp1
				! Filter the cur array
				do i=istr+1,ifin-2,2
					temp1=wtm1*cur(i-1,j,k)+wt*cur(i,j,k)+wtp1*cur(i+1,j,k)
					cur(i-1,j,k)=temp2
					temp2=wtm1*cur(i,j,k)+wt*cur(i+1,j,k)+wtp1*cur(i+2,j,k)
					cur(i,j,k)=temp1
				enddo
				if(.not. even) then ! i = ifin-1
					temp1=wtm1*cur(i-1,j,k)+wt*cur(i,j,k)+wtp1*cur(i+1,j,k)
					cur(i-1,j,k)=temp2
				else
					temp1 = temp2
				endif
				i=ifin
				if(ifin .ne. mx-3) then
					temp2=wtm1*cur(i-1,j,k)+wt*cur(i,j,k)+wtp1*cur(i+1,j,k)
					cur(i-1,j,k)=temp1
					cur(i,j,k)=temp2
				else
					temp2=wtm1*cur(i-1,j,k)+wt*cur(i,j,k)+wtp1*ghost(ntimes+1,j,k)
					cur(i-1,j,k)=temp1
					if(ntimes .gt. 1) then
						temp1=wtm1*cur(i,j,k)+wt*ghost(ntimes+1,j,k)+wtp1*ghost(ntimes+2,j,k)
					endif
					cur(i,j,k)=temp2
					! Filter the ghost cells following the cur array
					do i=ntimes+2,2*ntimes-2,2
						temp2=wtm1*ghost(i-1,j,k)+wt*ghost(i,j,k)+wtp1*ghost(i+1,j,k)
						ghost(i-1,j,k)=temp1
						temp1=wtm1*ghost(i,j,k)+wt*ghost(i+1,j,k)+wtp1*ghost(i+2,j,k)
						ghost(i,j,k)=temp2
					enddo
					if(ntimes .gt. 1) then
						if(.not.ntimes_even) then
							! i = 2*ntimes-1
							temp2=wtm1*ghost(i-1,j,k)+wt*ghost(i,j,k)+wtp1*ghost(i+1,j,k)
							ghost(i-1,j,k)=temp1
							temp1=temp2
						endif
						i = 2*ntimes
						ghost(i-1,j,k)=temp1
					endif
				endif
				n=n+1
			enddo
		enddo
	enddo

end subroutine filter_x

! filter_x, with copy->filter->writeback 
subroutine copy_filter_x(cur, ghost, istr, ifin, jstr, jfin, kstr, kfin)
	implicit none
	real, dimension(:,:,:) :: cur, ghost
! 	real, allocatable, dimension(:) :: temp_x
	real :: temp1, temp2
	real :: wtm1,wt,wtp1
	integer :: i,j,k,n
	integer :: istr, ifin, jstr, jfin, kstr, kfin
	integer :: mx_real
	logical :: even
	real :: total1, total2, total3

	total1 = 0.
	total2 = 0.
	total3 = 0.

	wtm1 = .25
	wt = .5
	wtp1 = .25

	mx_real = ifin-istr+1
! 	allocate(temp_x(ntimes+mx_real+ntimes)) ! ghost|cur|ghost

	even = .false.
	if(mod(ntimes+mx_real+ntimes,2) .eq. 0) then
		even = .true.
	endif

	do k=kstr,kfin
		do j=jstr,jfin
			call timer(38)
			! Copy to temp array
			temp_x(:ntimes) = ghost(:ntimes,j,k)
			temp_x(ntimes+1:ntimes+mx_real) = cur(istr:ifin,j,k)
			temp_x(ntimes+mx_real+1:) = ghost(ntimes+1:,j,k)
			call timer(38,tmstop=.true.)
			total1 = total1 + tm(38)

			call timer(39)
			! Filter
			n=1
			do while(n .le. ntimes)
				temp2=temp_x(1)
				do i=2,ntimes+mx_real+ntimes-2,2
					temp1=wtm1*temp_x(i-1)+wt*temp_x(i)+wtp1*temp_x(i+1)
					temp_x(i-1)=temp2
					temp2=wtm1*temp_x(i)+wt*temp_x(i+1)+wtp1*temp_x(i+2)
					temp_x(i)=temp1
				enddo
				if(.not.even) then
					! i = ntimes+mx_real+ntimes-1
					temp1=wtm1*temp_x(i-1)+wt*temp_x(i)+wtp1*temp_x(i+1)
					temp_x(i)=temp1
				endif
				temp_x(i-1)=temp2
				n=n+1
			enddo
			call timer(39, tmstop=.true.)
			total2 = total2 + tm(39)

			call timer(40)
			! Write out
			cur(istr:ifin,j,k) = temp_x(ntimes+1:ntimes+mx_real)
			call timer(40, tmstop=.true.)
			total3 = total3 + tm(40)
		enddo
	enddo

! 	if(rank .eq. 0) print *,"filterx totals",total1,total2,total3

! 	deallocate(temp_x)

end subroutine copy_filter_x

!-------------------------------------------------------------------------------
! 				subroutine filter_y()		
!										
! 	Does "ntimes" passes in the y direction on the specified array, using the 
!	specified ghost cells
!	Copies the cells to be filtered to temporary memory to avoid
!	cross-memory-grain access
!
!-------------------------------------------------------------------------------
subroutine filter_y(cur, ghost, istr, ifin, jstr, jfin, kstr, kfin)
	implicit none
	real, dimension(:,:,:) :: cur, ghost
	real :: temp1, temp2
	real :: wtm1,wt,wtp1
	integer :: i,j,k,n
	integer :: istr, ifin, jstr, jfin, kstr, kfin
	logical :: even, ntimes_even

	wtm1 = .25
	wt = .5
	wtp1 = .25

	ntimes_even = .false.
	if(mod(ntimes,2) .eq. 0) then
		ntimes_even = .true.
	endif

	even = .false.
	if(mod(jfin-jstr+1,2) .eq. 0) then
		even = .true.
	endif

	do k=kstr,kfin
		do i=istr,ifin
			n=1
			do while(n .le. ntimes)
				! Filter the ghost cells preceeding the cur array
				temp2 = ghost(i,1,k)
				do j=2,ntimes-2,2
					temp1=wtm1*ghost(i,j-1,k)+wt*ghost(i,j,k)+wtp1*ghost(i,j+1,k)
					ghost(i,j-1,k)=temp2
					temp2=wtm1*ghost(i,j,k)+wt*ghost(i,j+1,k)+wtp1*ghost(i,j+2,k)
					ghost(i,j,k)=temp1
				enddo
				if(ntimes .gt. 1) then
					if(.not.ntimes_even) then
						! i = ntimes-1
						temp1=wtm1*ghost(i,j-1,k)+wt*ghost(i,j,k)+wtp1*ghost(i,j+1,k)
						ghost(i,j-1,k)=temp2
						temp2=temp1
					endif
					j = ntimes
					temp1=wtm1*ghost(i,j-1,k)+wt*ghost(i,j,k)+wtp1*cur(i,jstr,k)
					ghost(i,j-1,k)=temp2
				else
					j=1
					temp1=temp2
				endif
				temp2=wtm1*ghost(i,j,k)+wt*cur(i,jstr,k)+wtp1*cur(i,jstr+1,k)
				ghost(i,j,k)=temp1
				! Filter the cur array
				do j=jstr+1,jfin-2,2
					temp1=wtm1*cur(i,j-1,k)+wt*cur(i,j,k)+wtp1*cur(i,j+1,k)
					cur(i,j-1,k)=temp2
					temp2=wtm1*cur(i,j,k)+wt*cur(i,j+1,k)+wtp1*cur(i,j+2,k)
					cur(i,j,k)=temp1
				enddo
				if(.not. even) then ! i = ifin-1
					temp1=wtm1*cur(i,j-1,k)+wt*cur(i,j,k)+wtp1*cur(i,j+1,k)
					cur(i,j-1,k)=temp2
				else
					temp1 = temp2
				endif
				j=jfin
				temp2=wtm1*cur(i,j-1,k)+wt*cur(i,j,k)+wtp1*ghost(i,ntimes+1,k)
				cur(i,j-1,k)=temp1
				if(ntimes .gt. 1) then
					temp1=wtm1*cur(i,j,k)+wt*ghost(i,ntimes+1,k)+wtp1*ghost(i,ntimes+2,k)
				endif
				cur(i,j,k)=temp2
				! Filter the ghost cells following the cur array
				do j=ntimes+2,2*ntimes-2,2
					temp2=wtm1*ghost(i,j-1,k)+wt*ghost(i,j,k)+wtp1*ghost(i,j+1,k)
					ghost(i,j-1,k)=temp1
					temp1=wtm1*ghost(i,j,k)+wt*ghost(i,j+1,k)+wtp1*ghost(i,j+2,k)
					ghost(i,j,k)=temp2
				enddo
				if(ntimes .gt. 1) then
					if(.not.ntimes_even) then
						! i = 2*ntimes-1
						temp2=wtm1*ghost(i,j-1,k)+wt*ghost(i,j,k)+wtp1*ghost(i,j+1,k)
						ghost(i,j-1,k)=temp1
						temp1=temp2
					endif
					j = 2*ntimes
					ghost(i,j-1,k)=temp1
				endif
				n=n+1
			enddo
		enddo
	enddo

end subroutine filter_y

subroutine copy_filter_y(cur, ghost, istr, ifin, jstr, jfin, kstr, kfin)
	implicit none
	real, dimension(:,:,:) :: cur, ghost
! 	real, allocatable, dimension(:) :: temp_y
	real :: temp1, temp2
	real :: wtm1,wt,wtp1
	integer :: i,j,k,n
	integer :: istr, ifin, jstr, jfin, kstr, kfin
	integer :: my_real
	logical :: even
	real :: total1, total2, total3

	total1 = 0.
	total2 = 0.
	total3 = 0.

	wtm1 = .25
	wt = .5
	wtp1 = .25

	my_real = jfin-jstr+1
! 	allocate(temp_y(ntimes+my_real+ntimes)) ! ghost|cur|ghost

	even = .false.
	if(mod(ntimes+my_real+ntimes,2) .eq. 0) then
		even = .true.
	endif

	do k=kstr,kfin
		do i=istr,ifin
			call timer(38)
			! Copy to temp array
			temp_y(:ntimes) = ghost(i,:ntimes,k)
			temp_y(ntimes+1:ntimes+my_real) = cur(i,jstr:jfin,k)
			temp_y(ntimes+my_real+1:) = ghost(i,ntimes+1:,k)
			call timer(38,tmstop=.true.)
			total1 = total1 + tm(38)

			call timer(39)
			! Filter
			n=1
			do while(n .le. ntimes)
				temp2=temp_y(1)
				do j=2,ntimes+my_real+ntimes-2,2
					temp1=wtm1*temp_y(j-1)+wt*temp_y(j)+wtp1*temp_y(j+1)
					temp_y(j-1)=temp2
					temp2=wtm1*temp_y(j)+wt*temp_y(j+1)+wtp1*temp_y(j+2)
					temp_y(j)=temp1
				enddo
				if(.not.even) then
					! j = ntimes+my_real+ntimes-1
					temp1=wtm1*temp_y(j-1)+wt*temp_y(j)+wtp1*temp_y(j+1)
					temp_y(j)=temp1
				endif
				temp_y(j-1)=temp2
				n=n+1
			enddo
			call timer(39, tmstop=.true.)
			total2 = total2 + tm(39)

			call timer(40)
			! Write out
			cur(i,jstr:jfin,k) = temp_y(ntimes+1:ntimes+my_real)
			call timer(40, tmstop=.true.)
			total3 = total3 + tm(40)
		enddo
	enddo

! 	if(rank .eq. 0) print *,"filtery totals",total1,total2,total3

! 	deallocate(temp_y)

end subroutine copy_filter_y

!-------------------------------------------------------------------------------
! 				subroutine filter_z()		
!										
! 	Does "ntimes" passes in the z direction on the specified array, using the 
!	specified ghost cells
!	Copies the cells to be filtered to temporary memory to avoid
!	cross-memory-grain access
!
!-------------------------------------------------------------------------------
subroutine filter_z(cur, ghost, istr, ifin, jstr, jfin, kstr, kfin)
	implicit none
	real, dimension(:,:,:) :: cur, ghost
	real :: temp1, temp2
	real :: wtm1,wt,wtp1
	integer :: i,j,k,n
	integer :: istr, ifin, jstr, jfin, kstr, kfin
	logical :: even, ntimes_even

	wtm1 = .25
	wt = .5
	wtp1 = .25

	ntimes_even = .false.
	if(mod(ntimes,2) .eq. 0) then
		ntimes_even = .true.
	endif

	even = .false.
	if(mod(kfin-kstr+1,2) .eq. 0) then
		even = .true.
	endif

	do j=jstr,jfin
		do i=istr,ifin
			n=1
			do while(n .le. ntimes)
				! Filter the ghost cells preceeding the cur array
				temp2 = ghost(i,j,1)
				do k=2,ntimes-2,2
					temp1=wtm1*ghost(i,j,k-1)+wt*ghost(i,j,k)+wtp1*ghost(i,j,k+1)
					ghost(i,j,k-1)=temp2
					temp2=wtm1*ghost(i,j,k)+wt*ghost(i,j,k+1)+wtp1*ghost(i,j,k+2)
					ghost(i,j,k)=temp1
				enddo
				if(ntimes .gt. 1) then
					if(.not.ntimes_even) then
						! i = ntimes-1
						temp1=wtm1*ghost(i,j,k-1)+wt*ghost(i,j,k)+wtp1*ghost(i,j,k+1)
						ghost(i,j,k-1)=temp2
						temp2=temp1
					endif
					k = ntimes
					temp1=wtm1*ghost(i,j,k-1)+wt*ghost(i,j,k)+wtp1*cur(i,j,kstr)
					ghost(i,j,k-1)=temp2
				else
					k=1
					temp1=temp2
				endif
				temp2=wtm1*ghost(i,j,k)+wt*cur(i,j,kstr)+wtp1*cur(i,j,kstr+1)
				ghost(i,j,k)=temp1
				! Filter the cur array
				do k=kstr+1,kfin-2,2
					temp1=wtm1*cur(i,j,k-1)+wt*cur(i,j,k)+wtp1*cur(i,j,k+1)
					cur(i,j,k-1)=temp2
					temp2=wtm1*cur(i,j,k)+wt*cur(i,j,k+1)+wtp1*cur(i,j,k+2)
					cur(i,j,k)=temp1
				enddo
				if(.not. even) then ! i = ifin-1
					temp1=wtm1*cur(i,j,k-1)+wt*cur(i,j,k)+wtp1*cur(i,j,k+1)
					cur(i,j,k-1)=temp2
				else
					temp1 = temp2
				endif
				k=kfin
				temp2=wtm1*cur(i,j,k-1)+wt*cur(i,j,k)+wtp1*ghost(i,j,ntimes+1)
				cur(i,j,k-1)=temp1
				if(ntimes .gt. 1) then
					temp1=wtm1*cur(i,j,k)+wt*ghost(i,j,ntimes+1)+wtp1*ghost(i,j,ntimes+2)
				endif
				cur(i,j,k)=temp2
				! Filter the ghost cells following the cur array
				do k=ntimes+2,2*ntimes-2,2
					temp2=wtm1*ghost(i,j,k-1)+wt*ghost(i,j,k)+wtp1*ghost(i,j,k+1)
					ghost(i,j,k-1)=temp1
					temp1=wtm1*ghost(i,j,k)+wt*ghost(i,j,k+1)+wtp1*ghost(i,j,k+2)
					ghost(i,j,k)=temp2
				enddo
				if(ntimes .gt. 1) then
					if(.not.ntimes_even) then
						! i = 2*ntimes-1
						temp2=wtm1*ghost(i,j,k-1)+wt*ghost(i,j,k)+wtp1*ghost(i,j,k+1)
						ghost(i,j,k-1)=temp1
						temp1=temp2
					endif
					k = 2*ntimes
					ghost(i,j,k-1)=temp1
				endif
				n=n+1
			enddo
		enddo
	enddo

end subroutine filter_z

subroutine copy_filter_z(cur, ghost, istr, ifin, jstr, jfin, kstr, kfin)
	implicit none
	real, dimension(:,:,:) :: cur, ghost
! 	real, allocatable, dimension(:) :: temp_z
	real :: temp1, temp2
	real :: wtm1,wt,wtp1
	integer :: i,j,k,n
	integer :: istr, ifin, jstr, jfin, kstr, kfin
	integer :: mz_real
	logical :: even

	wtm1 = .25
	wt = .5
	wtp1 = .25

	mz_real = kfin-kstr+1
! 	allocate(temp_z(ntimes+mz_real+ntimes)) ! ghost|cur|ghost

	even = .false.
	if(mod(ntimes+mz_real+ntimes,2) .eq. 0) then
		even = .true.
	endif

	do j=jstr,jfin
		do i=istr,ifin
			! Copy to temp array
			temp_z(:ntimes) = ghost(i,j,:ntimes)
			temp_z(ntimes+1:ntimes+mz_real) = cur(i,j,kstr:kfin)
			temp_z(ntimes+mz_real+1:) = ghost(i,j,ntimes+1:)
			! Filter
			n=1
			do while(n .le. ntimes)
				temp2=temp_z(1)
				do k=2,ntimes+mz_real+ntimes,2
					temp1=wtm1*temp_z(k-1)+wt*temp_z(k)+wtp1*temp_z(k+1)
					temp_z(k-1)=temp2
					temp2=wtm1*temp_z(k)+wt*temp_z(k+1)+wtp1*temp_z(k+2)
					temp_z(k)=temp1
				enddo
				if(.not.even) then
					! j = ntimes+my_real+ntimes-1
					temp1=wtm1*temp_z(k-1)+wt*temp_z(k)+wtp1*temp_z(k+1)
					temp_z(k)=temp1
				endif
				temp_z(k-1)=temp2
				n=n+1
			enddo
			! Write out
			cur(i,j,kstr:kfin) = temp_z(ntimes+1:ntimes+mz_real)
		enddo
	enddo

! 	deallocate(temp_z)

end subroutine copy_filter_z

subroutine apply_filter3_opt()
	implicit none
	
	call chunk_filter(curx)
	call chunk_filter(cury)
	call chunk_filter(curz)
	
end subroutine apply_filter3_opt

subroutine chunk_filter(cur)
	implicit none

	real, allocatable, dimension(:,:,:) :: cur

	! local variables
	
	real wtm1,wt,wtp1,winv
! 	real, allocatable, dimension(:,:,:) :: xghost, yghost, zghost
	integer :: istr, ifin, jstr, jfin, kstr, kfin,istr_filter, ifin_filter
	integer :: i, j, k, n, ierr
	integer :: num_x_chunks, num_y_chunks, num_z_chunks
	integer :: last_x_size, last_y_size, last_z_size

! 	x_chunk_len=(mx-5)/128
! 	y_chunk_len=(my-5)/4
! 	z_chunk_len=(mz-5)

	! istr,ifin are the boundaries for periodicity
	! istr_filter,ifin_filter are the filtering boundaries
	! i
	istr=3
	ifin=mx-3
	istr_filter=istr
	ifin_filter=ifin
! ! ! 	if(wall) then 
! ! ! 	   istr_filter=2		!right going injection
! ! ! 	   ifin_filter=min(mx-1,int(xinject2)+10)
! ! ! 	endif

! 	print *,"Wall: ",wall
! 	print *,"istr:",istr,"; ifin",ifin

! 	print *,"mx",mx,"my",my,"mz",mz

	! j
	jstr=3
	jfin=my-3

	! k
#ifndef twoD
	kstr=3
	kfin=mz-3 !2,mz-1 !3,mz-3
#else
	kstr=1
	kfin=1
#endif

! 	allocate(xghost(2*ntimes,my,mz), yghost(mx,2*ntimes,mz), zghost(mx,my,2*ntimes))
! 	xghost=0
! 	yghost=0
! 	zghost=0

	! Calculate number of chunks
	num_x_chunks = (mx-5)/x_chunk_len
	num_y_chunks = (my-5)/y_chunk_len
	num_z_chunks = (mz-5)/z_chunk_len

	last_x_size = (mx-5) - num_x_chunks*x_chunk_len
	last_y_size = (my-5) - num_y_chunks*y_chunk_len
	last_z_size = (mz-5) - num_z_chunks*z_chunk_len

	! Periodic boundary conditions
	!!! Populate xghost with curx data !!!
	if(periodicx.eq.1) then
		call deep_copylayrx(cur,xghost,istr,ifin,jstr,jfin,kstr,kfin)
	else
		call nonper_copylayrx(cur,xghost,istr,ifin,jstr,jfin,kstr,kfin)
	endif
	!!! Populate yghost with curx data !!!
	if(periodicy.eq.1) then
		call deep_copy_layry1(cur,yghost,istr,ifin,jstr,jfin,kstr,kfin)
	else
		print *, "Radiation BC not implemented in apply_filter"
	endif
#ifndef twoD
	!!! Populate zghost with curx data !!!
	if(periodicz.eq.1) then
		call deep_copy_layrz1(cur,zghost,istr,ifin,jstr,jfin,kstr,kfin)
	else
		!transfer z information among processors, not periodic between 0,n-1
		call deep_copy_layrz2(cur,zghost,istr,ifin,jstr,jfin,kstr,kfin)
	endif
#endif

	! Filter the chunks
	do k=1,num_z_chunks
		do j=1,num_y_chunks
			do i=1,num_x_chunks
				! Copy boundary conditions (need to add corners)
				! Lower x

				if(i .eq. 1) then
				   chunk(:ntimes,ntimes+1:ntimes &
					+y_chunk_len,ntimes+1:ntimes+z_chunk_len) &
					=xghost(:ntimes,3+(j-1)*y_chunk_len:3+ &
					j*y_chunk_len-1,3+(k-1)*z_chunk_len:3+ &
					k*z_chunk_len-1)
				else
				   chunk(:ntimes,ntimes+1:ntimes &
					+y_chunk_len,ntimes+1:ntimes+z_chunk_len) &
					=cur(3+(i-1)*x_chunk_len-ntimes:3+ &
					(i-1)*x_chunk_len-1,3+(j-1)*y_chunk_len:3+ &
					j*y_chunk_len-1,3+(k-1)*z_chunk_len:3+ &
					k*z_chunk_len-1)
				endif

				! Upper x
				if(i .eq. num_x_chunks) then
					chunk(ntimes+x_chunk_len+1:,ntimes+1:ntimes+y_chunk_len,ntimes+1:ntimes+z_chunk_len) =xghost(ntimes+1:,3+(j-1)*y_chunk_len:3+j*y_chunk_len-1,3+(k-1)*z_chunk_len:3+k*z_chunk_len-1)
				else
					chunk(ntimes+x_chunk_len+1:,ntimes+1:ntimes+y_chunk_len,ntimes+1:ntimes+z_chunk_len) =cur(3+i*x_chunk_len:3+i*x_chunk_len+ntimes-1,3+(j-1)*y_chunk_len:3+j*y_chunk_len-1,3+(k-1)*z_chunk_len:3+k*z_chunk_len-1)
				endif

				! Lower y
				if(j .eq. 1) then
					chunk(ntimes+1:ntimes+x_chunk_len,:ntimes,ntimes+1:ntimes+z_chunk_len) =yghost(3+(i-1)*x_chunk_len:3+i*x_chunk_len-1,:ntimes,3+(k-1)*z_chunk_len:3+k*z_chunk_len-1)
				else
					chunk(ntimes+1:ntimes+x_chunk_len,:ntimes,ntimes+1:ntimes+z_chunk_len) =cur(3+(i-1)*x_chunk_len:3+i*x_chunk_len-1,3+(j-1)*y_chunk_len-ntimes:3+(j-1)*y_chunk_len-1,3+(k-1)*z_chunk_len:3+k*z_chunk_len-1)
				endif

				! Upper y
				if(j .eq. num_y_chunks) then
					chunk(ntimes+1:ntimes+x_chunk_len,ntimes+y_chunk_len+1:,ntimes+1:ntimes+z_chunk_len) =yghost(3+(i-1)*x_chunk_len:3+i*x_chunk_len-1,ntimes+1:,3+(k-1)*z_chunk_len:3+k*z_chunk_len-1)
				else
					chunk(ntimes+1:ntimes+x_chunk_len,ntimes+y_chunk_len+1:,ntimes+1:ntimes+z_chunk_len) =cur(3+(i-1)*x_chunk_len:3+i*x_chunk_len-1,3+j*y_chunk_len:3+j*y_chunk_len+ntimes-1,3+(k-1)*z_chunk_len:3+k*z_chunk_len-1)
				endif

#ifndef twoD
				! Lower z
				if(k .eq. 1) then
					chunk(ntimes+1:ntimes+x_chunk_len,ntimes+1:ntimes+y_chunk_len,:ntimes) =zghost(3+(i-1)*x_chunk_len:3+i*x_chunk_len-1,3+(j-1)*y_chunk_len:3+j*y_chunk_len-1,:ntimes)
				else
					chunk(ntimes+1:ntimes+x_chunk_len,ntimes+1:ntimes+y_chunk_len,:ntimes) =cur(3+(i-1)*x_chunk_len:3+i*x_chunk_len-1,3+(j-1)*y_chunk_len:3+j*y_chunk_len-1,3+(k-1)*z_chunk_len-ntimes:3+(k-1)*z_chunk_len-1)
				endif

				! Upper z
				if(k .eq. num_z_chunks) then
					chunk(ntimes+1:ntimes+x_chunk_len,ntimes+1:ntimes+y_chunk_len,ntimes+z_chunk_len+1:) =zghost(3+(i-1)*x_chunk_len:3+i*x_chunk_len-1,3+(j-1)*y_chunk_len:3+j*y_chunk_len-1,ntimes+1:)
				else
					chunk(ntimes+1:ntimes+x_chunk_len,ntimes+1:ntimes+y_chunk_len,ntimes+z_chunk_len+1:) =cur(3+(i-1)*x_chunk_len:3+i*x_chunk_len-1,3+(j-1)*y_chunk_len:3+j*y_chunk_len-1,3+k*z_chunk_len:3+k*z_chunk_len+ntimes-1)
				endif
#endif 

				! Copy this chunk to the chunk array
				chunk(ntimes+1:ntimes+x_chunk_len,ntimes+1:ntimes+y_chunk_len,ntimes+1:ntimes+z_chunk_len) =cur(3+(i-1)*x_chunk_len:3+i*x_chunk_len-1,3+(j-1)*y_chunk_len:3+j*y_chunk_len-1,3+(k-1)*z_chunk_len:3+k*z_chunk_len-1)

				! Filter the chunk
				call chunk_filter_x(chunk,2,ntimes+x_chunk_len+ntimes-1,1,ntimes+y_chunk_len+ntimes,1,ntimes+z_chunk_len+ntimes)

			
				call chunk_filter_y(chunk,1,ntimes+x_chunk_len+ntimes,2,ntimes+y_chunk_len+ntimes-1,1,ntimes+z_chunk_len+ntimes)
#ifndef twoD

				call chunk_filter_z(chunk,1,ntimes+x_chunk_len+ntimes,1,ntimes+y_chunk_len+ntimes,2,ntimes+z_chunk_len+ntimes-1)
#endif

				! Writeback the chunk
				temp(3+(i-1)*x_chunk_len:3+i*x_chunk_len-1,3+(j-1)*y_chunk_len:3+j*y_chunk_len-1,3+(k-1)*z_chunk_len:3+k*z_chunk_len-1) =chunk(ntimes+1:ntimes+x_chunk_len,ntimes+1:ntimes+y_chunk_len,ntimes+1:ntimes+z_chunk_len)
			enddo
		enddo
	enddo

	cur(3:mx-3,3:my-3,3:mz-3)=temp(3:mx-3,3:my-3,3:mz-3)

end subroutine chunk_filter

! filter_x 
subroutine chunk_filter_x(cur, istr, ifin, jstr, jfin, kstr, kfin)
	implicit none
	real, dimension(:,:,:) :: cur
! 	real, allocatable, dimension(:) :: temp_x
	real :: temp1, temp2
	real :: wtm1,wt,wtp1
	integer :: i,j,k,n
	integer :: istr, ifin, jstr, jfin, kstr, kfin
	integer :: mx_real
	logical :: even
	real :: total1, total2, total3

	total1 = 0.
	total2 = 0.
	total3 = 0.

	wtm1 = .25
	wt = .5
	wtp1 = .25

	mx_real = ifin-istr+1
! 	allocate(temp_x(ntimes+mx_real+ntimes)) ! ghost|cur|ghost

	even = .false.
	if(mod(mx_real,2) .eq. 0) then
		even = .true.
	endif

	do k=kstr,kfin
		do j=jstr,jfin
			!call timer(39)
			! Filter
			n=1
			do while(n .le. ntimes)
				temp2=cur(istr-1,j,k)
				do i=istr,ifin-1,2
					temp1=wtm1*cur(i-1,j,k)+wt*cur(i,j,k)+wtp1*cur(i+1,j,k)
					cur(i-1,j,k)=temp2
					temp2=wtm1*cur(i,j,k)+wt*cur(i+1,j,k)+wtp1*cur(i+2,j,k)
					cur(i,j,k)=temp1
				enddo
				if(.not.even) then
					! i = ifin
					temp1=wtm1*cur(i-1,j,k)+wt*cur(i,j,k)+wtp1*cur(i+1,j,k)
					cur(i,j,k)=temp1
				endif
				cur(i-1,j,k)=temp2
				n=n+1
			enddo
			!call timer(39, tmstop=.true.)
			!total2 = total2 + tm(39)
		enddo
	enddo

! 	if(rank .eq. 0) print *,"filterx totals",total1,total2,total3

! 	deallocate(temp_x)

end subroutine chunk_filter_x

! filter_x, with copy->filter->writeback 
subroutine chunk_filter_y(cur, istr, ifin, jstr, jfin, kstr, kfin)
	implicit none
	real, dimension(:,:,:) :: cur
! 	real, allocatable, dimension(:) :: temp_x
	real :: temp1, temp2
	real :: wtm1,wt,wtp1
	integer :: i,j,k,n
	integer :: istr, ifin, jstr, jfin, kstr, kfin
	integer :: my_real
	logical :: even
	real :: total1, total2, total3

	total1 = 0.
	total2 = 0.
	total3 = 0.

	wtm1 = .25
	wt = .5
	wtp1 = .25

	my_real = jfin-jstr+1
! 	allocate(temp_x(ntimes+mx_real+ntimes)) ! ghost|cur|ghost

	even = .false.
	if(mod(my_real,2) .eq. 0) then
		even = .true.
	endif

	do k=kstr,kfin
		do i=istr,ifin
			!call timer(39)
			! Filter
			n=1
			do while(n .le. ntimes)
				temp2=cur(i,jstr-1,k)
				do j=jstr,jfin-1,2
					temp1=wtm1*cur(i,j-1,k)+wt*cur(i,j,k)+wtp1*cur(i,j+1,k)
					cur(i,j-1,k)=temp2
					temp2=wtm1*cur(i,j,k)+wt*cur(i,j+1,k)+wtp1*cur(i,j+2,k)
					cur(i,j,k)=temp1
				enddo
				if(.not.even) then
					! i = ifin
					temp1=wtm1*cur(i,j-1,k)+wt*cur(i,j,k)+wtp1*cur(i,j+1,k)
					cur(i,j,k)=temp1
				endif
				cur(i,j-1,k)=temp2
				n=n+1
			enddo
			!call timer(39, tmstop=.true.)
			!total2 = total2 + tm(39)
		enddo
	enddo

! 	if(rank .eq. 0) print *,"filterx totals",total1,total2,total3

! 	deallocate(temp_x)

end subroutine chunk_filter_y

! filter_x, with copy->filter->writeback 
subroutine chunk_filter_z(cur, istr, ifin, jstr, jfin, kstr, kfin)
	implicit none
	real, dimension(:,:,:) :: cur
! 	real, allocatable, dimension(:) :: temp_x
	real :: temp1, temp2
	real :: wtm1,wt,wtp1
	integer :: i,j,k,n
	integer :: istr, ifin, jstr, jfin, kstr, kfin
	integer :: mz_real
	logical :: even
	real :: total1, total2, total3

	total1 = 0.
	total2 = 0.
	total3 = 0.

	wtm1 = .25
	wt = .5
	wtp1 = .25

	mz_real = kfin-kstr+1
! 	allocate(temp_x(ntimes+mx_real+ntimes)) ! ghost|cur|ghost

	even = .false.
	if(mod(mz_real,2) .eq. 0) then
		even = .true.
	endif

	do j=jstr,jfin
		do i=istr,ifin
			!call timer(39)
			! Filter
			n=1
			do while(n .le. ntimes)
				temp2=cur(i,j,kstr-1)
				do k=kstr,kfin-1,2
					temp1=wtm1*cur(i,j,k-1)+wt*cur(i,j,k)+wtp1*cur(i,j,k+1)
					cur(i,j,k-1)=temp2
					temp2=wtm1*cur(i,j,k)+wt*cur(i,j,k+1)+wtp1*cur(i,j,k+2)
					cur(i,j,k)=temp1
				enddo
				if(.not.even) then
					! i = ifin
					temp1=wtm1*cur(i,j,k-1)+wt*cur(i,j,k)+wtp1*cur(i,j,k+1)
					cur(i,j,k)=temp1
				endif
				cur(i,j,k-1)=temp2
				n=n+1
			enddo
			!call timer(39, tmstop=.true.)
			!total2 = total2 + tm(39)
		enddo
	enddo

! 	if(rank .eq. 0) print *,"filterx totals",total1,total2,total3

! 	deallocate(temp_x)

end subroutine chunk_filter_z

!-------------------------------------------------------------------------------
! 						subroutine deep_copylayrx					 
!																		
! Copies ntimes layers from one side of the grid to the ghost grid (x direction)
!							
!-------------------------------------------------------------------------------

subroutine deep_copylayrx(b,ghost,ix,fx,iy,fy,iz,fz)

	implicit none
	
	! dummy variables

	integer, intent(in) :: ix, fx, iy, fy, iz, fz
	real(sprec), dimension(:,:,:), intent(in) :: b
	real(sprec), dimension(:,:,:), intent(out) :: ghost
	
	! local variables 
	
	integer :: i, j, k
	
! 	do k=iz,fz
! 		do j=iy,fy
! 			do i=1,ntimes
! 				ghost(i,j,k)=b(fx-ntimes+i,j,k)
! 			enddo
! 		enddo
! 	enddo
	ghost(:ntimes,iy:fy,iz:fz)=b(fx-ntimes+1:fx,iy:fy,iz:fz)
	
! 	do k=iz,fz
! 		do j=iy,fy
! 			do i=ntimes+1,2*ntimes
! 				ghost(i,j,k)=b(ix+i-ntimes-1,j,k)
! 			enddo
! 		enddo
! 	enddo
	ghost(ntimes+1:,iy:fy,iz:fz)=b(ix:ix+ntimes-1,iy:fy,iz:fz)

end subroutine deep_copylayrx

!-------------------------------------------------------------------------------
! 						subroutine nonper_copylayrx					 
!																		
! Simply copies the boundary values from b into the adjacent ghost grid
!							
!-------------------------------------------------------------------------------

subroutine nonper_copylayrx(b,ghost,ix,fx,iy,fy,iz,fz)

	implicit none
	
	! dummy variables

	integer, intent(in) :: ix, fx, iy, fy, iz, fz
	real(sprec), dimension(:,:,:), intent(in) :: b
	real(sprec), dimension(:,:,:), intent(out) :: ghost
	
	! local variables 
	
	integer :: i, j, k
	
	do k=iz,fz
		do j=iy,fy
			ghost(:ntimes,j,k)=b(ix,j,k)
		enddo
	enddo
	
	do k=iz,fz
		do j=iy,fy
			ghost(ntimes+1:2*ntimes,j,k)=b(fx,j,k)
		enddo
	enddo
	
end subroutine nonper_copylayrx

!-------------------------------------------------------------------------------
! 						subroutine deep_copy_layry1					 
!																		
! Copies ntimes layers from one side of the grid to the ghost grid across procs (y direction)
! 							
!-------------------------------------------------------------------------------

subroutine deep_copy_layrx1(b,ghost,ix,fx,iy,fy,iz,fz) 

	implicit none
	
	! dummy variables

	integer, intent(in) :: ix, fx, iy, fy, iz, fz
	real(sprec), dimension(:,:,:), intent(in) :: b
	real(sprec), dimension(:,:,:), intent(out) :: ghost
	
	! local variables
	
	integer :: plusrank, minusrank, iperiodic,comm,count,plustag,minustag,ierr,status(statsize)
	integer :: request(2), status1(statsize,2)
	integer :: i, j, k
! ! 	real(sprec), allocatable, dimension(:,:,:) :: bufferout

! ! 	allocate(bufferout(mx,ntimes,mz))

! 	print*,ix,fx,iy,fy,iz,fz

	plustag=100
	minustag=200
	comm=MPI_Comm_world
	
	plusrank=(rank/sizex)*sizex + modulo(rank+1,sizex)
	minusrank=(rank/sizex)*sizex + modulo(rank-1,sizex)		
		
! ! 	count=mx*ntimes*mz
	count=1

	!send up and recv from below
	call MPI_SendRecv(b,count,xhighcap,plusrank,plustag, &
	ghost,count,xghostlowcap,minusrank,plustag, &
	comm,status,ierr)

! ! ! 	ghost(ix:fx,:ntimes,iz:fz)=b(ix:fx,fy-ntimes+1:fy,iz:fz)

	!send dwn and recv from above
	call MPI_SendRecv(b,count,xlowcap,minusrank,minustag, &
	ghost,count,xghosthighcap,plusrank,minustag, &
	comm,status,ierr)

! ! ! 	ghost(ix:fx,ntimes+1:,iz:fz)=b(ix:fx,iy:iy+ntimes-1,iz:fz)

! ! 	call MPI_SendRecv(bufferout,count,mpi_read,lftrank,lfttag, &
! ! 	ghost(:,ntimes+1:,:),count,mpi_read,rgtrank,lfttag, &
! ! 	comm,status,ierr)

! ! 	!----------------- x field
! ! 
! ! 	!send up and recv from below
! ! 
! ! 	bufferout(:,:,:) = b(:,fz-ntimes+1:fz,:)
! ! 	
! ! 	call MPI_SendRecv(bufferout,count,mpi_read,rgtrank,rgttag, &
! ! 	ghost(:,:ntimes,:),count,mpi_read,lftrank,rgttag, &
! ! 	comm,status,ierr)
! ! 
! ! ! 	ghost(:,:ntimes,:) = bufferin(:,:,:)
! ! 
! ! 	call MPI_SendRecv(b(:,fz-ntimes+1:fz,:),count,mpi_read,rgtrank,rgttag, &
! ! 	ghost(:,:ntimes,:),count,mpi_read,lftrank,rgttag, &
! ! 	comm,status,ierr)
! ! 	
! ! 	!send dwn and recv from above
! ! 
! ! 	bufferout(:,:,:) = b(:,iz:iz+ntimes-1,:)
! ! 
! ! 	call MPI_SendRecv(bufferout,count,mpi_read,lftrank,lfttag, &
! ! 	ghost(:,ntimes+1:,:),count,mpi_read,rgtrank,lfttag, &
! ! 	comm,status,ierr)
! ! 
! ! ! 	ghost(:,ntimes+1:,:) = bufferin(:,:,:)
! ! 
! ! 	call MPI_SendRecv(b(:,iz:iz+ntimes-1,:),count,mpi_read,lftrank,lfttag, &
! ! 	ghost(:,ntimes+1:,:),count,mpi_read,rgtrank,lfttag, &
! ! 	comm,status,ierr)
! ! 
! ! 	deallocate(bufferout)

end subroutine deep_copy_layrx1

!-------------------------------------------------------------------------------
! 						subroutine deep_copy_layry2
!								
! Copies ntimes layers from one side of the grid to the ghost grid across procs (y direction)
! 							
!-------------------------------------------------------------------------------

subroutine deep_copy_layrx2(b,ghost,ix,fx,iy,fy,iz,fz) 

	implicit none
	
	! dummy variables

	integer, intent(in) :: ix, fx, iy, fy, iz, fz
	real(sprec), dimension(:,:,:), intent(in) :: b
	real(sprec), dimension(:,:,:), intent(out) :: ghost
	
	! local variables
	
	integer :: plusrank, minusrank, iperiodic,comm,count,plustag,minustag,ierr,status(statsize)
	integer :: request(2), status1(statsize,2)
	integer :: i, j, k
! ! 	real(sprec), allocatable, dimension(:,:,:) :: bufferout

! ! 	allocate(bufferout(mx,ntimes,mz))


	plustag=100
	minustag=200
	comm=MPI_Comm_world
	
	plusrank=(rank/sizex)*sizex + modulo(rank+1,sizex)
	minusrank=(rank/sizex)*sizex + modulo(rank-1,sizex)		
				

! ! 	count=mx*ntimes*mz
	count=1

	!send right and recv from left
	call MPI_SendRecv(b,count,xhighcap,plusrank,plustag, &
	ghost,count,xghostlowcap,minusrank,plustag, &
	comm,status,ierr)

!check if this processor is in the low-x range, and should not receive from the other side of the grid
	if(modulo(rank,sizex)+1 .eq. 1) then 
		do k=iz,fz
		   do j=iy,fy			
		   	ghost(:ntimes,j,k) = b(ix,j,k)
		   enddo
		 enddo		

	endif			

! ! ! 	ghost(ix:fx,:ntimes,iz:fz)=b(ix:fx,fy-ntimes+1:fy,iz:fz)

	!send left and recv from right
	call MPI_SendRecv(b,count,xlowcap,minusrank,minustag, &
	ghost,count,xghosthighcap,minusrank,minustag, &
	comm,status,ierr)

!check if this processor is in the high-y range, and should not receive from the other side of the grid
	if(modulo(rank,sizex)+1 .eq. sizex) then 
		do k=iz,fz
		   do j=iy,fy			
		   	ghost(ntimes+1:,j,k) = b(fx,j,k)
		   enddo
		 enddo		
	endif			

! ! ! 	ghost(ix:fx,ntimes+1:,iz:fz)=b(ix:fx,iy:iy+ntimes-1,iz:fz)

end subroutine deep_copy_layrx2



!-------------------------------------------------------------------------------
! 						subroutine deep_copy_layry1					 
!																		
! Copies ntimes layers from one side of the grid to the ghost grid across procs (y direction)
! 							
!-------------------------------------------------------------------------------

subroutine deep_copy_layry1(b,ghost,ix,fx,iy,fy,iz,fz) 

	implicit none
	
	! dummy variables

	integer, intent(in) :: ix, fx, iy, fy, iz, fz
	real(sprec), dimension(:,:,:), intent(in) :: b
	real(sprec), dimension(:,:,:), intent(out) :: ghost
	
	! local variables
	
	integer :: rgtrank, lftrank, iperiodic,comm,count,rgttag,lfttag,ierr,status(statsize)
	integer :: request(2), status1(statsize,2)
	integer :: i, j, k
! ! 	real(sprec), allocatable, dimension(:,:,:) :: bufferout

! ! 	allocate(bufferout(mx,ntimes,mz))

! 	print*,ix,fx,iy,fy,iz,fz

	rgttag=100
	lfttag=200
	comm=MPI_Comm_world
	
!	rgtrank=modulo((rank/sizex + 1),sizey)*sizex + modulo(rank,sizex) 
!	lftrank=modulo((rank/sizex - 1),sizey)*sizex + modulo(rank,sizex)

	rgtrank=modulo(rank/sizex + 1,sizey)*sizex+rank/(sizex*sizey)*(sizex*sizey) &
				+ modulo(rank,sizex) 
	lftrank=modulo(rank/sizex - 1,sizey)*sizex+rank/(sizex*sizey)*(sizex*sizey) &
				+ modulo(rank,sizex)		
		
! ! 	count=mx*ntimes*mz
	count=1

	!send up and recv from below
	call MPI_SendRecv(b,count,yhighcap,rgtrank,rgttag, &
	ghost,count,yghostlowcap,lftrank,rgttag, &
	comm,status,ierr)

! ! ! 	ghost(ix:fx,:ntimes,iz:fz)=b(ix:fx,fy-ntimes+1:fy,iz:fz)

	!send dwn and recv from above
	call MPI_SendRecv(b,count,ylowcap,lftrank,lfttag, &
	ghost,count,yghosthighcap,rgtrank,lfttag, &
	comm,status,ierr)

! ! ! 	ghost(ix:fx,ntimes+1:,iz:fz)=b(ix:fx,iy:iy+ntimes-1,iz:fz)

! ! 	call MPI_SendRecv(bufferout,count,mpi_read,lftrank,lfttag, &
! ! 	ghost(:,ntimes+1:,:),count,mpi_read,rgtrank,lfttag, &
! ! 	comm,status,ierr)

! ! 	!----------------- x field
! ! 
! ! 	!send up and recv from below
! ! 
! ! 	bufferout(:,:,:) = b(:,fz-ntimes+1:fz,:)
! ! 	
! ! 	call MPI_SendRecv(bufferout,count,mpi_read,rgtrank,rgttag, &
! ! 	ghost(:,:ntimes,:),count,mpi_read,lftrank,rgttag, &
! ! 	comm,status,ierr)
! ! 
! ! ! 	ghost(:,:ntimes,:) = bufferin(:,:,:)
! ! 
! ! 	call MPI_SendRecv(b(:,fz-ntimes+1:fz,:),count,mpi_read,rgtrank,rgttag, &
! ! 	ghost(:,:ntimes,:),count,mpi_read,lftrank,rgttag, &
! ! 	comm,status,ierr)
! ! 	
! ! 	!send dwn and recv from above
! ! 
! ! 	bufferout(:,:,:) = b(:,iz:iz+ntimes-1,:)
! ! 
! ! 	call MPI_SendRecv(bufferout,count,mpi_read,lftrank,lfttag, &
! ! 	ghost(:,ntimes+1:,:),count,mpi_read,rgtrank,lfttag, &
! ! 	comm,status,ierr)
! ! 
! ! ! 	ghost(:,ntimes+1:,:) = bufferin(:,:,:)
! ! 
! ! 	call MPI_SendRecv(b(:,iz:iz+ntimes-1,:),count,mpi_read,lftrank,lfttag, &
! ! 	ghost(:,ntimes+1:,:),count,mpi_read,rgtrank,lfttag, &
! ! 	comm,status,ierr)
! ! 
! ! 	deallocate(bufferout)

end subroutine deep_copy_layry1

!-------------------------------------------------------------------------------
! 						subroutine deep_copy_layry2
!								
! Copies ntimes layers from one side of the grid to the ghost grid across procs (y direction)
! 							
!-------------------------------------------------------------------------------

subroutine deep_copy_layry2(b,ghost,ix,fx,iy,fy,iz,fz) 

	implicit none
	
	! dummy variables

	integer, intent(in) :: ix, fx, iy, fy, iz, fz
	real(sprec), dimension(:,:,:), intent(in) :: b
	real(sprec), dimension(:,:,:), intent(out) :: ghost
	
	! local variables
	
	integer :: rgtrank, lftrank, iperiodic,comm,count,rgttag,lfttag,ierr,status(statsize)
	integer :: request(2), status1(statsize,2)
	integer :: i, j, k
! ! 	real(sprec), allocatable, dimension(:,:,:) :: bufferout

! ! 	allocate(bufferout(mx,ntimes,mz))


	rgttag=100
	lfttag=200
	comm=MPI_Comm_world
	
!	rgtrank=modulo((rank/sizex + 1),sizey)*sizex + modulo(rank,sizex) 
!	lftrank=modulo((rank/sizex - 1),sizey)*sizex + modulo(rank,sizex)
		
	rgtrank=modulo(rank/sizex + 1,sizey)*sizex+rank/(sizex*sizey)*(sizex*sizey) &
				+ modulo(rank,sizex) 
	lftrank=modulo(rank/sizex - 1,sizey)*sizex+rank/(sizex*sizey)*(sizex*sizey) &
				+ modulo(rank,sizex)		
		

! ! 	count=mx*ntimes*mz
	count=1

	!send right and recv from left
	call MPI_SendRecv(b,count,yhighcap,rgtrank,rgttag, &
	ghost,count,yghostlowcap,lftrank,rgttag, &
	comm,status,ierr)

!check if this processor is in the low-y range, and should not receive from the other side of the grid
	if(modulo(rank,sizex*sizey)/sizex+1 .eq. 1) then 
		do k=iz,fz
		   do i=ix,fx			
		   	ghost(i,:ntimes,k) = b(i,iy,k)
		   enddo
		 enddo		

	endif			

! ! ! 	ghost(ix:fx,:ntimes,iz:fz)=b(ix:fx,fy-ntimes+1:fy,iz:fz)

	!send left and recv from right
	call MPI_SendRecv(b,count,ylowcap,lftrank,lfttag, &
	ghost,count,yghosthighcap,rgtrank,lfttag, &
	comm,status,ierr)

!check if this processor is in the high-y range, and should not receive from the other side of the grid
	if(modulo(rank,sizex*sizey)/sizex+1 .eq. sizey) then 
		do k=iz,fz
		   do i=ix,fx			
		   	ghost(i,ntimes+1:,k) = b(i,fy,k)
		   enddo
		 enddo		
	endif			

! ! ! 	ghost(ix:fx,ntimes+1:,iz:fz)=b(ix:fx,iy:iy+ntimes-1,iz:fz)

end subroutine deep_copy_layry2


!-------------------------------------------------------------------------------
! 						subroutine deep_copy_layrz1			 
!													
! Copies a layer from one side of the grid to the other across procs (z direction)
! (using sendrecv)				
!-------------------------------------------------------------------------------

subroutine deep_copy_layrz1(b,ghost,ix,fx,iy,fy,iz,fz)

	implicit none
	
	! dummy variables

	integer, intent(in) :: ix, fx, iy, fy, iz, fz
	real(sprec), dimension(:,:,:), intent(in) :: b
	real(sprec), dimension(:,:,:), intent(out) :: ghost

	! local variables
	
	integer ::uprank, dwnrank, iperiodic,comm,count,uptag,dwntag,ierr,status(statsize)
	integer ::request(2), status1(statsize,2)
	integer :: i, j, k
! ! 	real(sprec), allocatable, dimension(:,:,:) :: bufferout, bufferin

! ! 	allocate(bufferout(mx,my,ntimes), bufferin(mx,my,ntimes))

	uptag=100
	dwntag=200
	comm=MPI_Comm_world
	
#ifndef twoD 
	
	!if twoD don't exchange anything
	
	uprank=modulo(rank/(sizex*sizey) + 1,sizez)*(sizex*sizey) + & 
		 		  modulo(rank,sizex*sizey) 
	dwnrank=modulo(rank/(sizex*sizey) - 1,sizez)*(sizex*sizey) + & 
			      modulo(rank,sizex*sizey)		
		
! ! 	count=mx*my*ntimes
	count=1

	!send up and recv from below
	call MPI_SendRecv(b,count,zhighcap,uprank,uptag, &
	ghost,count,zghostlowcap,dwnrank,uptag, &
	comm,status,ierr)

! ! ! 	ghost(ix:fx,iy:fy,:ntimes)=b(ix:fx,iy:fy,fz-ntimes+1:fz)

	!send dwn and recv from above
	call MPI_SendRecv(b,count,zlowcap,dwnrank,dwntag, &
	ghost,count,zghosthighcap,uprank,dwntag, &
	comm,status,ierr)


! ! ! 	ghost(ix:fx,iy:fy,ntimes+1:)=b(ix:fx,iy:fy,iz:iz+ntimes-1)

! ! 	!----------------- x field
! ! 	
! ! 	!send up and recv from below
! ! 
! ! 	bufferout(:,:,:) = b(:,:,fz-ntimes+1:fz)
! ! 
! ! 	call MPI_SendRecv(bufferout,count,mpi_read,uprank,uptag, &
! ! 	bufferin,count,mpi_read,dwnrank,uptag, &
! ! 	comm,status,ierr)
! ! 
! ! 	ghost(:,:,:ntimes) = bufferin(:,:,:)
! ! 	
! ! 	!send dwn and recv from above
! ! 
! ! 	bufferout(:,:,:) = b(:,:,iz:iz+ntimes-1)
! ! 	
! ! 	call MPI_SendRecv(bufferout,count,mpi_read,dwnrank,dwntag, &
! ! 	bufferin,count,mpi_read,uprank,dwntag, &
! ! 	comm,status,ierr)
! ! 
! ! 	ghost(:,:,ntimes+1:) = bufferin(:,:,:)
		
#else
	ghost(1,1,1)=0. ! Prevents the compiler from complaining
#endif

! 	deallocate(bufferout, bufferin)

end subroutine deep_copy_layrz1


!-------------------------------------------------------------------------------
! 						subroutine copy_layrz2					 
!																		
! Copies a layer from one side of the grid to the other across procs (z direction)
! (using sendrecv)											
!-------------------------------------------------------------------------------
subroutine deep_copy_layrz2(b,ghost,ix,fx,iy,fy,iz,fz) 

	implicit none
	
	! dummy variables

	integer, intent(in) :: ix, fx, iy, fy, iz, fz
	real(sprec), dimension(:,:,:), intent(in) :: b
	real(sprec), dimension(:,:,:), intent(out) :: ghost
	integer :: i, j, k

	! local variables
	
	integer ::uprank, dwnrank, iperiodic,comm,count,uptag,dwntag,ierr,status(statsize)
	integer ::request(2), status1(statsize,2)
	real(sprec), allocatable, dimension(:,:,:) :: bufferin
	
	allocate(bufferin(mx,my,2*ntimes))


	uptag=100
	dwntag=200
	comm=MPI_Comm_world
	
#ifndef twoD 
	
	!if twoD don't exchange anything

	uprank=modulo(rank/(sizex*sizey) + 1,sizez)*(sizex*sizey) + & 
		 		  modulo(rank,sizex*sizey) 
	dwnrank=modulo(rank/(sizex*sizey) - 1,sizez)*(sizex*sizey) + & 
			      modulo(rank,sizex*sizey)
	
! ! 	count=mx*my*ntimes
	count=1
	
	!send up and recv from below
	
	call MPI_SendRecv(b,count,zhighcap,uprank,uptag, &
	bufferin,count,zghostlowcap,dwnrank,uptag, &
	comm,status,ierr)
	
	if(rank/(sizex*sizey).ge.1) then  !these processors are not in the bottom row in z
		ghost(:,:,:ntimes) = bufferin(:,:,:ntimes)
	else  !these are in the bottom, just copy the first element into the ghost zones
		do j=iy,fy
		   do i=ix,fx			
		   	ghost(i,j,:ntimes) = b(i,j,iz)
		   enddo
		 enddo		   
	endif	

	!send dwn and recv from above
	
	call MPI_SendRecv(b,count,zlowcap,dwnrank,dwntag, &
	bufferin,count,zghosthighcap,uprank,dwntag, &
	comm,status,ierr)

	if(rank/(sizex*sizey)+1 .ne. sizez ) then !these processors are not in the top row
		ghost(:,:,ntimes+1:) = bufferin(:,:,ntimes+1:)
	else !these are in the top row
		do j=iy,fy
		   do i=ix,fx			
		   	ghost(i,j,ntimes+1:) = b(i,j,fz)
		   enddo
		 enddo	
	endif

#else
	ghost(1,1,1)=0. ! Prevents the compiler from complaining	
#endif

	deallocate(bufferin)

end subroutine deep_copy_layrz2
