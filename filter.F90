!-------------------------------------------------------------------------------
! 				subroutine apply_filter1_opt()		
!										
! 
!
!-------------------------------------------------------------------------------

subroutine apply_filter1_opt()

	implicit none
	
	! local variables
	
!	real wtm1,wt,wtp1,winv
	integer :: istr, ifin, i, j, k, n

 real wtl,wtr,wtu,wtb,wtlt,wtrt,wtlb,wtrb,wt, winv
#ifndef twoD
 real wtlz,wtrz,wtuz,wtbz,wtltz,wtrtz,wtlbz,wtrbz, wtz
#endif

! apply digital filter to currents to smooth them 
! 1 2 1 in every direction
! done as convolution of 3 one-dimensional sweeps
! filtering is repeated ntimes	
	
	n=1
	
	istr=2  +1
	ifin=mx-1 -2

#ifdef twoD
   winv=1./16.
   wtl=2.*winv
   wt=4.*winv
   wtr=2.*winv
   wtu=2.*winv
   wtb=2.*winv
   wtlt=1.*winv
   wtrt=1.*winv
   wtlb=1.*winv
   wtrb=1.*winv

#else
   winv=1./64.
   wtl=4.*winv
   wt=8.*winv
   wtr=4.*winv
   wtu=4.*winv
   wtb=4.*winv
   wtlt=2.*winv
   wtrt=2.*winv
   wtlb=2.*winv
   wtrb=2.*winv

   wtlz=2.*winv
   wtz=4.*winv
   wtrz=2.*winv
   wtuz=2.*winv
   wtbz=2.*winv
   wtltz=1.*winv
   wtrtz=1.*winv
   wtlbz=1.*winv
   wtrbz=1.*winv

#endif
	do while(n .le. ntimes) 
		temp=0. !a little risky, but we probably don't need to clean temp
		
		!filter in x direction
		if(periodicx.eq.1) then 
!			call copy_layrx1(curx,cury,curz,mx,my,mz,2,mx-3,mx-2,3)
			call copy_layrx1_opt(curx,cury,curz,mx,my,mz,2,mx-3,mx-2,3)
		else
!			call copy_layrx2(curx,cury,curz,mx,my,mz,2,mx-3,mx-2,3)
			call copy_layrx2_opt(curx,cury,curz,mx,my,mz,2,mx-3,mx-2,3)
		endif

		if(periodicy.eq.1) then 
!			call copy_layry1(curx,cury,curz,mx,my,mz,2,my-3,my-2,3)
			call copy_layry1_opt(curx,cury,curz,mx,my,mz,2,my-3,my-2,3)
		else
!			call copy_layry2(curx,cury,curz,mx,my,mz,2,my-3,my-2,3)
			call copy_layry2_opt(curx,cury,curz,mx,my,mz,2,my-3,my-2,3)
		endif
#ifndef twoD
			if(periodicz.eq.1) then
				!call copy_layrz1(curx,cury,curz,mx,my,mz,2,mz-3,mz-2,3)
				 call copy_layrz1_opt(curx,cury,curz,mx,my,mz,2,mz-3,mz-2,3)		
			else
				!transfer z information among processors, not periodic between 0,n-1
				!call copy_layrz2(curx,cury,curz,mx,my,mz,2,mz-3,mz-2,3)
				call copy_layrz2_opt(curx,cury,curz,mx,my,mz,2,mz-3,mz-2,3)
			endif
#endif

#ifndef twoD
		do k=3,mz-3 
#else
		do k=1,1
#endif
		do j=3,my-3 
			do i=istr,ifin 
				temp(i,j,k)=wtl*curx(i-1,j,k)+wt*curx(i,j,k)+wtr*curx(i+1,j,k)+ &
                                            wtb*curx(i,j-1,k)+wtu*curx(i,j+1,k)+ &
                                            wtlt*curx(i-1,j+1,k)+wtrt*curx(i+1,j+1,k)+ &
                                            wtlb*curx(i-1,j-1,k)+wtrb*curx(i+1,j-1,k)
#ifndef twoD
              temp(i,j,k)=temp(i,j,k)+wtlz*curx(i-1,j,k-1)+wtz*curx(i,j,k-1)+wtrz*curx(i+1,j,k-1)+ &
                                            wtbz*curx(i,j-1,k-1)+wtuz*curx(i,j+1,k-1)+ &
                                            wtltz*curx(i-1,j+1,k-1)+wtrtz*curx(i+1,j+1,k-1)+ &
                                            wtlbz*curx(i-1,j-1,k-1)+wtrbz*curx(i+1,j-1,k-1)+ &
                                      wtlz*curx(i-1,j,k+1)+wtz*curx(i,j,k+1)+wtrz*curx(i+1,j,k+1)+ &
                                            wtbz*curx(i,j-1,k+1)+wtuz*curx(i,j+1,k+1)+ &
                                            wtltz*curx(i-1,j+1,k+1)+wtrtz*curx(i+1,j+1,k+1)+ &
                                            wtlbz*curx(i-1,j-1,k+1)+wtrbz*curx(i+1,j-1,k+1)
#endif
			enddo
		enddo
		enddo


#ifndef twoD
		do k=3,mz-3
#else
		do k=1,1
#endif
			do j=3,my-3 
				do i=istr,ifin 
					curx(i,j,k)=temp(i,j,k)
				enddo
			enddo
		enddo


#ifndef twoD
		do k=3,mz-3 
#else
		do k=1,1
#endif
		do j=3,my-3 
			do i=istr,ifin 
				temp(i,j,k)=wtl*cury(i-1,j,k)+wt*cury(i,j,k)+wtr*cury(i+1,j,k)+ &
                                            wtb*cury(i,j-1,k)+wtu*cury(i,j+1,k)+ &
                                            wtlt*cury(i-1,j+1,k)+wtrt*cury(i+1,j+1,k)+ &
                                            wtlb*cury(i-1,j-1,k)+wtrb*cury(i+1,j-1,k)

#ifndef twoD
              temp(i,j,k)=temp(i,j,k)+wtlz*cury(i-1,j,k-1)+wtz*cury(i,j,k-1)+wtrz*cury(i+1,j,k-1)+ &
                                            wtbz*cury(i,j-1,k-1)+wtuz*cury(i,j+1,k-1)+ &
                                            wtltz*cury(i-1,j+1,k-1)+wtrtz*cury(i+1,j+1,k-1)+ &
                                            wtlbz*cury(i-1,j-1,k-1)+wtrbz*cury(i+1,j-1,k-1)+ &
                                      wtlz*cury(i-1,j,k+1)+wtz*cury(i,j,k+1)+wtrz*cury(i+1,j,k+1)+ &
                                            wtbz*cury(i,j-1,k+1)+wtuz*cury(i,j+1,k+1)+ &
                                            wtltz*cury(i-1,j+1,k+1)+wtrtz*cury(i+1,j+1,k+1)+ &
                                            wtlbz*cury(i-1,j-1,k+1)+wtrbz*cury(i+1,j-1,k+1)
#endif

			enddo
		enddo
		enddo

#ifndef twoD
		do k=3,mz-3
#else
		do k=1,1
#endif
			do j=3,my-3 
				do i=istr,ifin 
					cury(i,j,k)=temp(i,j,k)
				enddo
			enddo
		enddo



#ifndef twoD
		do k=3,mz-3 
#else
		do k=1,1
#endif
		do j=3,my-3 
			do i=istr,ifin 
				temp(i,j,k)=wtl*curz(i-1,j,k)+wt*curz(i,j,k)+wtr*curz(i+1,j,k)+ &
                                            wtb*curz(i,j-1,k)+wtu*curz(i,j+1,k)+ &
                                            wtlt*curz(i-1,j+1,k)+wtrt*curz(i+1,j+1,k)+ &
                                            wtlb*curz(i-1,j-1,k)+wtrb*curz(i+1,j-1,k)
#ifndef twoD
              temp(i,j,k)=temp(i,j,k)+wtlz*curz(i-1,j,k-1)+wtz*curz(i,j,k-1)+wtrz*curz(i+1,j,k-1)+ &
                                            wtbz*curz(i,j-1,k-1)+wtuz*curz(i,j+1,k-1)+ &
                                            wtltz*curz(i-1,j+1,k-1)+wtrtz*curz(i+1,j+1,k-1)+ &
                                            wtlbz*curz(i-1,j-1,k-1)+wtrbz*curz(i+1,j-1,k-1)+ &
                                      wtlz*curz(i-1,j,k+1)+wtz*curz(i,j,k+1)+wtrz*curz(i+1,j,k+1)+ &
                                            wtbz*curz(i,j-1,k+1)+wtuz*curz(i,j+1,k+1)+ &
                                            wtltz*curz(i-1,j+1,k+1)+wtrtz*curz(i+1,j+1,k+1)+ &
                                            wtlbz*curz(i-1,j-1,k+1)+wtrbz*curz(i+1,j-1,k+1)
#endif
			enddo
		enddo
		enddo

#ifndef twoD
		do k=3,mz-3
#else
		do k=1,1
#endif
			do j=3,my-3 
				do i=istr,ifin 
					curz(i,j,k)=temp(i,j,k)
				enddo
			enddo
		enddo
 
		n=n+1
	enddo

end subroutine apply_filter1_opt
