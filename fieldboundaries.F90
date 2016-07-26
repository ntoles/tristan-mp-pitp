!
! Field boundaries module
!
! This module contains the boundary conditions for the fields, including communication 
! (MPI related) functions
!
!

#ifdef twoD 

module m_fieldboundaries

	use m_system
	use m_aux
	use m_communications
	use m_fields
	use m_globaldata
	use m_inputparser
	
#else

module m_fieldboundaries_3d

	use m_system_3d
	use m_aux_3d
	use m_communications_3d
	use m_fields_3d
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
!-------------------------------------------------------------------------------
!	INTERFACE DECLARATIONS
!-------------------------------------------------------------------------------
	
!-------------------------------------------------------------------------------
!	PUBLIC MODIFIERS
!-------------------------------------------------------------------------------

	! public functions

	public :: exchange_current, pre_bc_b, pre_bc_e, copylayrx, post_bc_b, &
			  post_bc_e, bc_b1, bc_b2, bc_e2, bc_e1, apply_filter1_opt, &
			  read_input_boundaries, conduct_bc
#ifdef filter2
        public :: apply_filter2_opt, apply_filter2_copy_opt, apply_filter3_opt
#endif	
!-------------------------------------------------------------------------------
!	MODULE PROCEDURES AND FUNCTIONS
!-------------------------------------------------------------------------------

	contains



!-------------------------------------------------------------------------------
! 					subroutine read_input_boundaries					 
!																		
! Reads any variables related to (or needed by) this module
! 							
!-------------------------------------------------------------------------------

subroutine read_input_boundaries()

	implicit none


	call inputpar_geti_def("boundaries", "periodicx", 1, periodicx)
	call inputpar_geti_def("boundaries", "periodicy", 1, periodicy)
	call inputpar_geti_def("boundaries", "periodicz", 1, periodicz)
	
	radiationx=1-periodicx
	radiationy=1-periodicy
	radiationz=1-periodicz
#ifdef twoD
	if(radiationy .eq. 1) radiationz=1 !this is to trigger radiaiton BCs in 2D
#endif

	if(periodicx .eq. 0) then 
	x1in=3 !set the location of planes where particles are removed from simulation, perpendicular to x. 
	x2in=mx0-2
	! these are the default locations, which can be overwritten in the user module 
        ! (init_particle_distribution_user)
	endif

end subroutine read_input_boundaries



!-------------------------------------------------------------------------------
! 						subroutine pre_bc_b					 
!																		
! 
! 							
!-------------------------------------------------------------------------------

subroutine pre_bc_b()

	implicit none
	

	if(radiationx.eq.1 .and. radiationy .eq. 1 .and. radiationz.eq. 1) then
!	if(radiationx.eq.1 .and. radiationy .eq. 1 .and. radiationz.eq. 1) then
!		if(size0.gt.1) then 
!			if(rank.eq.0)print *, rank,": fix preledge b"! 
!		else
!	   print *, "b4 preledge" !hack
#ifndef twoD
			call preledge(by,bz,bx,ey,ez,ex,iy,iz,ix,my,mz,mx,1,c)
!			print *, "b4 preledge 2" !hack
			call preledge(bz,bx,by,ez,ex,ey,iz,ix,iy,mz,mx,my,1,c)
!			print *, "b4 preledge 3" !hack
			call preledge(bx,by,bz,ex,ey,ez,ix,iy,iz,mx,my,mz,1,c)
!		endif
!			print *, "aft preledge b4 bc_b1" !hack
#endif
			call bc_b1()
	endif

end subroutine pre_bc_b



!-------------------------------------------------------------------------------
! 						subroutine pre_bc_e					 
!																		
!
! 							
!-------------------------------------------------------------------------------

subroutine pre_bc_e()

	implicit none

	!radiation BCs are not implemented yet
	if(radiationx.eq.1 .and. radiationy .eq. 1 ) then !.and. radiationz.eq. 1) then
!	if(radiationx.eq.1 .and. radiationy .eq. 1 .and. radiationz.eq. 1) then

!		if(size0.gt.1) then 
!			if(rank .eq. 0) print *, rank,": fix preledge e"! 
!		else
#ifndef twoD
			call preledge(ey,ez,ex,by,bz,bx,-iy,-iz,-ix,my,mz,mx,lot,c)
			call preledge(ez,ex,ey,bz,bx,by,-iz,-ix,-iy,mz,mx,my,lot,c)
			call preledge(ex,ey,ez,bx,by,bz,-ix,-iy,-iz,mx,my,mz,lot,c)
#endif
!		endif

				call bc_e1()
	endif

end subroutine pre_bc_e



!-------------------------------------------------------------------------------
! 						subroutine bc_b1					 
!																		
!	Boundary conditions for the magnetic field
!							
!-------------------------------------------------------------------------------

subroutine bc_b1()

	implicit none
	

	if(sizex .ne. 1) then
		if(periodicx.eq.1) then
			call copy_layrx1_opt(bx,by,bz,mx,my,mz,2,mx-3,mx-2,3)
			if( highorder) then
				call copy_layrx1_opt(bx,by,bz,mx,my,mz,1,mx-4,mx-1,4)
			endif
		else
			call copy_layrx2_opt(bx,by,bz,mx,my,mz,2,mx-3,mx-2,3)
			if( highorder) then
				call copy_layrx2_opt(bx,by,bz,mx,my,mz,1,mx-4,mx-1,4)
			endif
		endif
	else
		if(periodicx.eq.1) then
			call copylayrx(bx,by,bz,mx,my,mz,2,mx-3,mx-2,3)
			if( highorder) then
				call copylayrx(bx,by,bz,mx,my,mz,1,mx-4,mx-1,4)
			endif
		endif
	endif	
			
	if(sizey .ne. 1) then
		if(periodicy.eq.1) then
			call copy_layry1_opt(bx,by,bz,mx,my,mz,2,my-3,my-2,3)
			if( highorder) then
				call copy_layry1_opt(bx,by,bz,mx,my,mz,1,my-4,my-1,4)
			endif
		else
			call copy_layry2_opt(bx,by,bz,mx,my,mz,2,my-3,my-2,3)
			if( highorder) then
				call copy_layry2_opt(bx,by,bz,mx,my,mz,1,my-4,my-1,4)
			endif
		endif
	else
		if(periodicy.eq.1) then
			call copylayry(bx,by,bz,mx,my,mz,2,my-3,my-2,3)
			if( highorder) then
				call copylayry(bx,by,bz,mx,my,mz,1,my-4,my-1,4)
			endif
		endif
	endif				
					
	if(periodicz.eq.1) then
		call copy_layrz1_opt(bx,by,bz,mx,my,mz,2,mz-3,mz-2,3)
		if( highorder) then
			call copy_layrz1_opt(bx,by,bz,mx,my,mz,1,mz-4,mz-1,4)
		endif
	else
		!transfer z information among processors, not periodic between 0,n-1
		call copy_layrz2_opt(bx,by,bz,mx,my,mz,2,mz-3,mz-2,3)
		if( highorder) then
			call copy_layrz2_opt(bx,by,bz,mx,my,mz,1,mz-4,mz-1,4)
		endif
	endif

end subroutine bc_b1



!-------------------------------------------------------------------------------
! 						subroutine bc_b2					 
!																		
! Boundary conditions for the magnetic field
! 							
!-------------------------------------------------------------------------------

subroutine bc_b2()

	implicit none
	

	if( (radiationx .eq. 1)) then 
		call surface(by,bz,bx,ey,ez,ex,iy,iz,ix,my,mz,mx,1,c) 
	endif

	if(radiationy .eq. 1) then 
		call surface(bz,bx,by,ez,ex,ey,iz,ix,iy,mz,mx,my,1,c) 
	endif

#ifndef twoD
	if(radiationz .eq. 1) then 
		call surface(bx,by,bz,ex,ey,ez,ix,iy,iz,mx,my,mz,1,c)
	endif
#endif

	call bc_b1()

end subroutine bc_b2



!-------------------------------------------------------------------------------
! 						subroutine bc_e1					 
!																		
! Boundary conditions for the electric field
!							
!-------------------------------------------------------------------------------

subroutine bc_e1()

	implicit none

	if(sizex .ne. 1) then
		if(periodicx.eq.1) then
			call copy_layrx1_opt(ex,ey,ez,mx,my,mz,2,mx-3,mx-2,3)
			if( highorder) then
				call copy_layrx1_opt(ex,ey,ez,mx,my,mz,1,mx-4,mx-1,4)
			endif
		else
			call copy_layrx2_opt(ex,ey,ez,mx,my,mz,2,mx-3,mx-2,3)
			if( highorder) then
				call copy_layrx2_opt(ex,ey,ez,mx,my,mz,1,mx-4,mx-1,4)
			endif
		endif
	else
		if(periodicx.eq.1) then
			call copylayrx(ex,ey,ez,mx,my,mz,2,mx-3,mx-2,3)
			if( highorder) then
				call copylayrx(ex,ey,ez,mx,my,mz,1,mx-4,mx-1,4)
			endif
		endif
	endif	
	
	if(sizey .ne. 1) then
		if(periodicy.eq.1) then
			call copy_layry1_opt(ex,ey,ez,mx,my,mz,2,my-3,my-2,3)
			if( highorder) then
				call copy_layry1_opt(ex,ey,ez,mx,my,mz,1,my-4,my-1,4)
			endif
		else
			call copy_layry2_opt(ex,ey,ez,mx,my,mz,2,my-3,my-2,3)
			if( highorder) then
				call copy_layry2_opt(ex,ey,ez,mx,my,mz,1,my-4,my-1,4)
			endif
		endif
	else	
		if(periodicy.eq.1) then
			call copylayry(ex,ey,ez,mx,my,mz,2,my-3,my-2,3)
			if( highorder) then
				call copylayry(ex,ey,ez,mx,my,mz,1,my-4,my-1,4)
			endif
		endif
	endif			
			
	if(periodicz.eq.1) then
		call copy_layrz1_opt(ex,ey,ez,mx,my,mz,2,mz-3,mz-2,3)
		if( highorder) then
			call copy_layrz1_opt(ex,ey,ez,mx,my,mz,1,mz-4,mz-1,4)
		endif
	else
		!transfer z information among processors, not periodic between 0,n-1
		call copy_layrz2_opt(ex,ey,ez,mx,my,mz,2,mz-3,mz-2,3)
		if( highorder) then
			call copy_layrz2_opt(ex,ey,ez,mx,my,mz,1,mz-4,mz-1,4)
		endif
	endif

end subroutine bc_e1



!-------------------------------------------------------------------------------
! 						subroutine bc_e2					 
!																		
! Boundary conditions for the electric field
!							
!-------------------------------------------------------------------------------

subroutine bc_e2()

	implicit none
	
	integer :: i, j , k
	

	if((radiationx .eq. 1)) then 
		call surface(ey,ez,ex,by,bz,bx,-iy,-iz,-ix,my,mz,mx,lot,c)
	endif

	if(radiationy .eq. 1) then 
		call surface(ez,ex,ey,bz,bx,by,-iz,-ix,-iy,mz,mx,my,lot,c)
	endif

#ifndef twoD
	if(radiationz .eq. 1) then 
		call surface(ex,ey,ez,bx,by,bz,-ix,-iy,-iz,mx,my,mz,lot,c)
	endif
#endif

	call bc_e1()

end subroutine bc_e2



!-------------------------------------------------------------------------------
! 						subroutine post_bc_b					 
!																		
! 
!  							
!-------------------------------------------------------------------------------

subroutine post_bc_b()

	implicit none
	

	if(radiationx .eq. 1 .and. radiationy .eq. 1 .and. radiationz .eq.1) then 
#ifndef twoD
			call postedge(by,bz,bx,ey,ez,ex,iy,iz,ix,my,mz,mx,1,c)
			call postedge(bz,bx,by,ez,ex,ey,iz,ix,iy,mz,mx,my,1,c)
			call postedge(bx,by,bz,ex,ey,ez,ix,iy,iz,mx,my,mz,1,c)
#endif

			call bc_b1()
	endif

end subroutine post_bc_b



!-------------------------------------------------------------------------------
! 						subroutine post_bc_e					 
!																		
! 
! 							
!-------------------------------------------------------------------------------

subroutine post_bc_e()

	implicit none
	

	if(radiationx .eq. 1 .and. radiationy .eq. 1 .and. radiationz .eq.1) then 
!	if(radiationx .eq. 1 .and. radiationy .eq. 1  .and. radiationz .eq.1) then
!		if(size0.gt.1) then 	 
!			if(rank.eq.0)print *, "fix postedge e"! 
!		else							!call postedge
#ifndef twoD
			call postedge(ey,ez,ex,by,bz,bx,-iy,-iz,-ix,my,mz,mx,lot,c)
			call postedge(ez,ex,ey,bz,bx,by,-iz,-ix,-iy,mz,mx,my,lot,c)
			call postedge(ex,ey,ez,bx,by,bz,-ix,-iy,-iz,mx,my,mz,lot,c)	
#endif		
			call bc_e1()
!		endif
	endif

end subroutine post_bc_e



!-------------------------------------------------------------------------------
! 						subroutine surface					 
!																		
! 
! 							
!-------------------------------------------------------------------------------

subroutine surface(bx,by,bz,ex,ey,ez,ix,iy,iz,mx,my,mz,m00,c)

	implicit none

	! dummy variables

	real(sprec) :: bx(1), by(1), bz(1), ex(1), ey(1), ez(1)
!	real(sprec), dimension(:), intent(inout) :: bx, by, bz, ex, ey, ez
	integer, intent(in) :: ix, iy, iz, mx, my, mz, m00
	real(sprec), intent(in) :: c
	
	! local variables

	integer :: m, n 	
	real(sprec) :: rs, s, os
	
	rs=2.*c/(1.+c)
	s=.4142136
	os=.5*(1.-s)*rs

#ifndef twoD
		do m=m00+iz*(mz-1),m00+iz*(mz-1)+iy*(my-2),iy
			do n=m,m+ix*(mx-2),ix
				bz(n)=bz(n)+.5*c*(ex(n+iy)-ex(n)-ey(n+ix)+ey(n))
			enddo
			do n=m+ix,m+ix*(mx-2),ix
				bx(n)=bx(n)+rs*(bx(n-iz)-bx(n)+s*(bz(n)-bz(n-ix)))-os*( &
				ez(n+iy)-ez(n))-(os-c)*(ez(n+iy-iz)-ez(n-iz)) &
				-c*(ey(n)-ey(n-iz))
			enddo
		enddo
		
		do m=m00+iz*(mz-1),m00+iz*(mz-1)+ix*(mx-2),ix
			do n=m+iy,m+iy*(my-2),iy
				by(n)=by(n)+rs*(by(n-iz)-by(n)+s*(bz(n)-bz(n-iy)))+os*( &
				ez(n+ix)-ez(n))+(os-c)*(ez(n+ix-iz)-ez(n-iz))+c*(ex(n)-ex(n-iz))
			enddo
			do n=m,m+iy*(my-2),iy
				bz(n)=bz(n)+.5*c*(ex(n+iy)-ex(n)-ey(n+ix)+ey(n))
			enddo
		enddo
#else

		if(ix .eq. 0) then 
		do m=m00+iz*(mz-1),m00+iz*(mz-1)+iy*(my-2),iy
			 n=m !,m+ix*(mx-2),ix
				bz(n)=bz(n)+.5*c*(ex(n+iy)-ex(n)-ey(n+ix)+ey(n))
			
!			do n=m+ix,m+ix*(mx-2),ix
				bx(n)=bx(n)+rs*(bx(n-iz)-bx(n)+s*(bz(n)-bz(n-ix)))-os*( &
				ez(n+iy)-ez(n))-(os-c)*(ez(n+iy-iz)-ez(n-iz)) &
				-c*(ey(n)-ey(n-iz))
!			enddo
		enddo
		
!		do m=m00+iz*(mz-1),m00+iz*(mz-1)+ix*(mx-2),ix
		   m=m00+iz*(mz-1)
			do n=m+iy,m+iy*(my-2),iy
				by(n)=by(n)+rs*(by(n-iz)-by(n)+s*(bz(n)-bz(n-iy)))+os*( &
				ez(n+ix)-ez(n))+(os-c)*(ez(n+ix-iz)-ez(n-iz))+c*(ex(n)-ex(n-iz))
			enddo
			do n=m,m+iy*(my-2),iy
				bz(n)=bz(n)+.5*c*(ex(n+iy)-ex(n)-ey(n+ix)+ey(n))
			enddo
!		enddo
		endif




		if(iy .eq. 0) then 
		   m=m00+iz*(mz-1)
	
		   do n=m,m+ix*(mx-2),ix
		      bz(n)=bz(n)+.5*c*(ex(n+iy)-ex(n)-ey(n+ix)+ey(n))
		   enddo
		   do n=m+ix,m+ix*(mx-2),ix
		      bx(n)=bx(n)+rs*(bx(n-iz)-bx(n)+s*(bz(n)-bz(n-ix)))-os*( &
		      ez(n+iy)-ez(n))-(os-c)*(ez(n+iy-iz)-ez(n-iz)) &
		      -c*(ey(n)-ey(n-iz))
		   enddo
		   do m=m00+iz*(mz-1),m00+iz*(mz-1)+ix*(mx-2),ix
		      n=m
		      by(n)=by(n)+rs*(by(n-iz)-by(n)+s*(bz(n)-bz(n-iy)))+os*( &
		      ez(n+ix)-ez(n))+(os-c)*(ez(n+ix-iz)-ez(n-iz))+c*(ex(n)-ex(n-iz))
		      bz(n)=bz(n)+.5*c*(ex(n+iy)-ex(n)-ey(n+ix)+ey(n))
		   enddo
		endif
		if(iz .eq. 0) then 
		do m=m00+iz*(mz-1),m00+iz*(mz-1)+iy*(my-2),iy
			do n=m,m+ix*(mx-2),ix
				bz(n)=bz(n)+.5*c*(ex(n+iy)-ex(n)-ey(n+ix)+ey(n))
			enddo
			do n=m+ix,m+ix*(mx-2),ix
				bx(n)=bx(n)+rs*(bx(n-iz)-bx(n)+s*(bz(n)-bz(n-ix)))-os*( &
				ez(n+iy)-ez(n))-(os-c)*(ez(n+iy-iz)-ez(n-iz)) &
				-c*(ey(n)-ey(n-iz))
			enddo
		enddo
		
		do m=m00+iz*(mz-1),m00+iz*(mz-1)+ix*(mx-2),ix
			do n=m+iy,m+iy*(my-2),iy
				by(n)=by(n)+rs*(by(n-iz)-by(n)+s*(bz(n)-bz(n-iy)))+os*( &
				ez(n+ix)-ez(n))+(os-c)*(ez(n+ix-iz)-ez(n-iz))+c*(ex(n)-ex(n-iz))
			enddo
			do n=m,m+iy*(my-2),iy
				bz(n)=bz(n)+.5*c*(ex(n+iy)-ex(n)-ey(n+ix)+ey(n))
			enddo
		enddo
		endif
#endif


end subroutine surface

subroutine conduct_bc()
	if(radiationx.eq.1 .and. radiationy .eq. 1 .and. radiationz.eq. 1) then
!	   call conduct(ix,iy,iz,mx,my,mz,ey,ez,.25) !should split this by sides
!	   call conduct(iy,iz,ix,my,mz,mx,ez,ex,.25)
!	   call conduct(iz,ix,iy,mz,mx,my,ex,ey,.25)
	endif
end subroutine conduct_bc

subroutine conduct(ix,iy,iz,mx,my,mz,ey,ez,s)
	implicit none
!       dimension ey(1),ez(1)
	real(sprec) :: ey(1),ez(1), s, temp
	integer, intent(in) :: ix, iy, iz, mx, my, mz
	integer i,j,k
	
!  field arrays are treated as one-dimensional arrays in this subroutine
!  Current density components jy=s*ey, jz=s*ez are spread laterally in the
!  proportion .25,.5,.25 required by our smoothing convention for all field
!  sources.
	do  i=1+2*ix,1+ix*(mx-3),ix*(mx-5)
	   do  j=i+2*iy,i+iy*(my-4),iy
	      do  k=j+2*iz,j+iz*(mz-4),iz
		 temp=.25*s*ey(k)
		 ey(k-ix)=ey(k-ix)-temp
		 ey(k+ix)=ey(k+ix)-temp
		 ey(k)=ey(k)-2.*temp
		 temp=.25*s*ez(k)
		 ez(k-ix)=ez(k-ix)-temp
		 ez(k+ix)=ez(k+ix)-temp
		 ez(k)=ez(k)-2.*temp
	      enddo
	   enddo
	enddo
end subroutine conduct


subroutine surface_par(bx,by,bz,ex,ey,ez,ix,iy,iz,mx,my,mz,m00,c)

	implicit none

	! dummy variables

	real(sprec) :: bx(1), by(1), bz(1), ex(1), ey(1), ez(1)
!	real(sprec), dimension(:), intent(inout) :: bx, by, bz, ex, ey, ez
	integer, intent(in) :: ix, iy, iz, mx, my, mz, m00
	real(sprec), intent(in) :: c
	
	! local variables

	integer :: m, n 	
	real(sprec) :: rs, s, os
	
	rs=2.*c/(1.+c)
	s=.4142136
	os=.5*(1.-s)*rs

#ifndef twoD
		do m=m00+iz*(mz-1),m00+iz*(mz-1)+iy*(my-2),iy
			do n=m,m+ix*(mx-2),ix
				bz(n)=bz(n)+.5*c*(ex(n+iy)-ex(n)-ey(n+ix)+ey(n))
			enddo
			do n=m+ix,m+ix*(mx-2),ix
				bx(n)=bx(n)+rs*(bx(n-iz)-bx(n)+s*(bz(n)-bz(n-ix)))-os*( &
				ez(n+iy)-ez(n))-(os-c)*(ez(n+iy-iz)-ez(n-iz)) &
				-c*(ey(n)-ey(n-iz))
			enddo
		enddo
		
		do m=m00+iz*(mz-1),m00+iz*(mz-1)+ix*(mx-2),ix
			do n=m+iy,m+iy*(my-2),iy
				by(n)=by(n)+rs*(by(n-iz)-by(n)+s*(bz(n)-bz(n-iy)))+os*( &
				ez(n+ix)-ez(n))+(os-c)*(ez(n+ix-iz)-ez(n-iz))+c*(ex(n)-ex(n-iz))
			enddo
			do n=m,m+iy*(my-2),iy
				bz(n)=bz(n)+.5*c*(ex(n+iy)-ex(n)-ey(n+ix)+ey(n))
			enddo
		enddo
#else
		m=m00+iz*(mz-1)
	
		do n=m,m+ix*(mx-2),ix
			bz(n)=bz(n)+.5*c*(ex(n+iy)-ex(n)-ey(n+ix)+ey(n))
		enddo
		do n=m+ix,m+ix*(mx-2),ix
			bx(n)=bx(n)+rs*(bx(n-iz)-bx(n)+s*(bz(n)-bz(n-ix)))-os*( &
			ez(n+iy)-ez(n))-(os-c)*(ez(n+iy-iz)-ez(n-iz)) &
			-c*(ey(n)-ey(n-iz))
		enddo
		do m=m00+iz*(mz-1),m00+iz*(mz-1)+ix*(mx-2),ix
			n=m
		by(n)=by(n)+rs*(by(n-iz)-by(n)+s*(bz(n)-bz(n-iy)))+os*( &
				ez(n+ix)-ez(n))+(os-c)*(ez(n+ix-iz)-ez(n-iz))+c*(ex(n)-ex(n-iz))
			bz(n)=bz(n)+.5*c*(ex(n+iy)-ex(n)-ey(n+ix)+ey(n))
		enddo
#endif


end subroutine surface_par

!-------------------------------------------------------------------------------
! 						subroutine copylayrx					 
!																		
! Copies a layer from one side of the grid to the other (x direction)
!							
!-------------------------------------------------------------------------------

subroutine copylayrx(bx,by,bz,mx,my,mz,lt,ls,nt,ns)

	implicit none
	
	! dummy variables

	integer, intent(in) :: mx,my,mz,lt,ls,nt,ns
	real(sprec), dimension(:,:,:), intent(inout) :: bx, by, bz
	
	! local variables 
	
	integer :: j, k
	
	
	do k=1,mz
		do j=1,my
			bx(lt,j,k)=bx(ls,j,k)
			by(lt,j,k)=by(ls,j,k)
			bz(lt,j,k)=bz(ls,j,k)
		enddo
	enddo
	
	do k=1,mz
		do j=1,my
			bx(nt,j,k)=bx(ns,j,k)
			by(nt,j,k)=by(ns,j,k)
			bz(nt,j,k)=bz(ns,j,k)
		enddo
	enddo

end subroutine copylayrx



!-------------------------------------------------------------------------------
! 						subroutine copy_layry					 
!																		
! Copies a layer from one side of the grid to the other (y direction)
! 							
!-------------------------------------------------------------------------------

subroutine copylayry(bx,by,bz,mx,my,mz,lt,ls,nt,ns)

	implicit none

	! dummy variables

	integer, intent(in) :: mx,my,mz,lt,ls,nt,ns
	real(sprec), dimension(:,:,:), intent(inout) :: bx, by, bz
	
	! local variables 
	
	integer :: i, k
	
	
	do k=1,mz
		do i=1,mx
			bx(i,lt,k)=bx(i,ls,k)
			by(i,lt,k)=by(i,ls,k)
			bz(i,lt,k)=bz(i,ls,k)
		enddo
	enddo

	do k=1,mz
		do i=1,mx
			bx(i,nt,k)=bx(i,ns,k)
			by(i,nt,k)=by(i,ns,k)
			bz(i,nt,k)=bz(i,ns,k)
		enddo
	enddo

end subroutine copylayry



!-------------------------------------------------------------------------------
! 						subroutine copy_layrx1					 
!																		
! Copies a layer from one side of the grid to the other across procs (x direction)
! 							
!-------------------------------------------------------------------------------

subroutine copy_layrx1(bx,by,bz,mx,my,mz,lt,ls,nt,ns)

	implicit none
	
	! dummy variables

	integer, intent(in) ::mx,my,mz,lt,ls,nt,ns
	real(sprec), dimension(:,:,:), intent(inout) :: bx, by, bz
	
	! local variables
	
	integer :: plusrank, minusrank, iperiodic,comm,count,plustag,minustag,ierr,status(statsize)
	integer :: request(2), status1(statsize,2)


	plustag=100
	minustag=200
	comm=MPI_Comm_world
	
	plusrank=(rank/sizex)*sizex + modulo(rank+1,sizex)
	minusrank=(rank/sizex)*sizex + modulo(rank-1,sizex)

	count=my*mz

	!----------------- x field

	!send up and recv from below
	
	call MPI_SendRecv(bx(ls,:,:),count,mpi_read,plusrank,plustag, &
	bx(lt,:,:),count,mpi_read,minusrank,plustag, &
	comm,status,ierr)
	
	!send dwn and recv from above
	
	call MPI_SendRecv(bx(ns,:,:),count,mpi_read,minusrank,minustag, &
	bx(nt,:,:),count,mpi_read,plusrank,minustag, &
	comm,status,ierr)
	
	!----------------- y field
	
	!send up and recv from below
	
	call MPI_SendRecv(by(ls,:,:),count,mpi_read,plusrank,plustag, &
	by(lt,:,:),count,mpi_read,minusrank,plustag, &
	comm,status,ierr)
	
	!send dwn and recv from above
	
	call MPI_SendRecv(by(ns,:,:),count,mpi_read,minusrank,minustag, &
	by(nt,:,:),count,mpi_read,plusrank,minustag, &
	comm,status,ierr)
	
	!----------------- z field
	
	!send up and recv from below
	
	call MPI_SendRecv(bz(ls,:,:),count,mpi_read,plusrank,plustag, &
	bz(lt,:,:),count,mpi_read,minusrank,plustag, &
	comm,status,ierr)
	
	!send dwn and recv from above
	
	call MPI_SendRecv(bz(ns,:,:),count,mpi_read,minusrank,minustag, &
	bz(nt,:,:),count,mpi_read,plusrank,minustag, &
	comm,status,ierr)

end subroutine copy_layrx1


!-------------------------------------------------------------------------------
! 						subroutine copy_layrx2		 
!																		
! Copies a layer from one side of the grid to the other across procs (y direction)
! 							
!-------------------------------------------------------------------------------

subroutine copy_layrx2(bx,by,bz,mx,my,mz,lt,ls,nt,ns)

	implicit none
	
	! dummy variables

	integer, intent(in) ::mx,my,mz,lt,ls,nt,ns
	real(sprec), dimension(:,:,:), intent(inout) :: bx, by, bz
	
	! local variables
	
	integer :: plusrank, minusrank, iperiodic,comm,count,plustag,minustag,ierr,status(statsize)
	integer :: request(2), status1(statsize,2)			
			
	real bufferin(my,mz)

	plustag=100
	minustag=200
	comm=MPI_Comm_world
	
	plusrank=(rank/sizex)*sizex + modulo(rank+1,sizex)
	minusrank=(rank/sizex)*sizex + modulo(rank-1,sizex)

	count=my*mz

	!----------------- x field

	!send up and recv from below
	
	call MPI_SendRecv(bx(ls,:,:),count,mpi_read,plusrank,plustag, &
	bufferin,count,mpi_read,minusrank,plustag, &
	comm,status,ierr)

	if(modulo(rank,sizex)+1 .ne. 1)bx(lt,:,:)=bufferin
	
	!send dwn and recv from above
	
	call MPI_SendRecv(bx(ns,:,:),count,mpi_read,minusrank,minustag, &
	bufferin,count,mpi_read,plusrank,minustag, &
	comm,status,ierr)

	if (modulo(rank,sizex)+1 .ne. sizex) bx(nt,:,:)=bufferin

	!----------------- y field
	
	!send up and recv from below
	
	call MPI_SendRecv(by(ls,:,:),count,mpi_read,plusrank,plustag, &
	bufferin,count,mpi_read,minusrank,plustag, &
	comm,status,ierr)

	if(modulo(rank,sizex)+1 .ne. 1)by(lt,:,:)=bufferin
	
	!send dwn and recv from above
	
	call MPI_SendRecv(by(ns,:,:),count,mpi_read,minusrank,minustag, &
	bufferin,count,mpi_read,plusrank,minustag, &
	comm,status,ierr)

	if (modulo(rank,sizex)+1 .ne. sizex) by(nt,:,:)=bufferin
	
	!----------------- z field
	
	!send up and recv from below
	
	call MPI_SendRecv(bz(ls,:,:),count,mpi_read,plusrank,plustag, &
	bufferin,count,mpi_read,minusrank,plustag, &
	comm,status,ierr)

	if(modulo(rank,sizex)+1 .ne. 1)bz(lt,:,:)=bufferin
	
	!send dwn and recv from above
	
	call MPI_SendRecv(bz(ns,:,:),count,mpi_read,minusrank,minustag, &
	bufferin,count,mpi_read,plusrank,minustag, &
	comm,status,ierr)

	if (modulo(rank,sizex)+1 .ne. sizex) bz(nt,:,:)=bufferin

end subroutine copy_layrx2





!-------------------------------------------------------------------------------
! 						subroutine copy_layry1					 
!																		
! Copies a layer from one side of the grid to the other across procs (y direction)
! 							
!-------------------------------------------------------------------------------

subroutine copy_layry1(bx,by,bz,mx,my,mz,lt,ls,nt,ns)

	implicit none
	
	! dummy variables

	integer, intent(in) ::mx,my,mz,lt,ls,nt,ns
	real(sprec), dimension(:,:,:), intent(inout) :: bx, by, bz
	
	! local variables
	
	integer :: rgtrank, lftrank, iperiodic,comm,count,rgttag,lfttag,ierr,status(statsize)
	integer :: request(2), status1(statsize,2)


	rgttag=100
	lfttag=200
	comm=MPI_Comm_world
	
!	rgtrank=modulo((rank/sizex + 1),sizey)*sizex + modulo(rank,sizex) 
!	lftrank=modulo((rank/sizex - 1),sizey)*sizex + modulo(rank,sizex)

	rgtrank=modulo(rank/sizex + 1,sizey)*sizex+rank/(sizex*sizey)*(sizex*sizey) &
				+ modulo(rank,sizex) 
	lftrank=modulo(rank/sizex - 1,sizey)*sizex+rank/(sizex*sizey)*(sizex*sizey) &
				+ modulo(rank,sizex) 
				
    count=mx*mz

	!----------------- x field

	!send up and recv from below
	
	call MPI_SendRecv(bx(:,ls,:),count,mpi_read,rgtrank,rgttag, &
	bx(:,lt,:),count,mpi_read,lftrank,rgttag, &
	comm,status,ierr)
	
	!send dwn and recv from above
	
	call MPI_SendRecv(bx(:,ns,:),count,mpi_read,lftrank,lfttag, &
	bx(:,nt,:),count,mpi_read,rgtrank,lfttag, &
	comm,status,ierr)
	
	!----------------- y field
	
	!send up and recv from below
	
	call MPI_SendRecv(by(:,ls,:),count,mpi_read,rgtrank,rgttag, &
	by(:,lt,:),count,mpi_read,lftrank,rgttag, &
	comm,status,ierr)
	
	!send dwn and recv from above
	
	call MPI_SendRecv(by(:,ns,:),count,mpi_read,lftrank,lfttag, &
	by(:,nt,:),count,mpi_read,rgtrank,lfttag, &
	comm,status,ierr)
	
	!----------------- z field
	
	!send up and recv from below
	
	call MPI_SendRecv(bz(:,ls,:),count,mpi_read,rgtrank,rgttag, &
	bz(:,lt,:),count,mpi_read,lftrank,rgttag, &
	comm,status,ierr)
	
	!send dwn and recv from above
	
	call MPI_SendRecv(bz(:,ns,:),count,mpi_read,lftrank,lfttag, &
	bz(:,nt,:),count,mpi_read,rgtrank,lfttag, &
	comm,status,ierr)

end subroutine copy_layry1

!!!need to implement copy_layry2 as copy_layrz2 
! if (modulo(rank,sizey)+1 .eq. 1) 
! if (modulo(rank,sizey)+1 .eq. sizey)

!-------------------------------------------------------------------------------
! 						subroutine copy_layry2		 
!																		
! Copies a layer from one side of the grid to the other across procs (y direction)
! 							
!-------------------------------------------------------------------------------

subroutine copy_layry2(bx,by,bz,mx,my,mz,lt,ls,nt,ns)

	implicit none
	
	! dummy variables

	integer, intent(in) ::mx,my,mz,lt,ls,nt,ns
	real(sprec), dimension(:,:,:), intent(inout) :: bx, by, bz
	
	! local variables
	
	integer :: rgtrank, lftrank, iperiodic,comm,count,rgttag,lfttag,ierr,status(statsize)
	integer :: request(2), status1(statsize,2)
	
	real bufferin(mx,mz)

	rgttag=100
	lfttag=200
	comm=MPI_Comm_world
	
!	rgtrank=modulo((rank/sizex + 1),sizey)*sizex + modulo(rank,sizex) 
!	lftrank=modulo((rank/sizex - 1),sizey)*sizex + modulo(rank,sizex)

	rgtrank=modulo(rank/sizex + 1,sizey)*sizex+rank/(sizex*sizey)*(sizex*sizey) &
				+ modulo(rank,sizex) 
	lftrank=modulo(rank/sizex - 1,sizey)*sizex+rank/(sizex*sizey)*(sizex*sizey) &
				+ modulo(rank,sizex) 
		
	count=mx*mz

	!----------------- x field

	!send up and recv from below
	
	call MPI_SendRecv(bx(:,ls,:),count,mpi_read,rgtrank,rgttag, &
	bufferin,count,mpi_read,lftrank,rgttag, &
	comm,status,ierr)

!	if(rank/sizex+1 .ne. 1) bx(:,lt,:)=bufferin
	if(modulo(rank,sizex*sizey)/sizex .ne. 0) bx(:,lt,:)=bufferin	
		
	
	!send dwn and recv from above
	
	call MPI_SendRecv(bx(:,ns,:),count,mpi_read,lftrank,lfttag, &
	bufferin,count,mpi_read,rgtrank,lfttag, &
	comm,status,ierr)

!	if(rank/sizex+1 .ne. sizey) bx(:,nt,:)=bufferin
	if(modulo(rank,sizex*sizey)/sizex .ne. sizey-1) bx(:,nt,:)=bufferin	

	!----------------- y field
	
	!send up and recv from below
	
	call MPI_SendRecv(by(:,ls,:),count,mpi_read,rgtrank,rgttag, &
	bufferin,count,mpi_read,lftrank,rgttag, &
	comm,status,ierr)

!	if(rank/sizex+1 .ne. 1)by(:,lt,:)=bufferin
	if(modulo(rank,sizex*sizey)/sizex .ne. 0) by(:,lt,:)=bufferin
		
	!send dwn and recv from above
	
	call MPI_SendRecv(by(:,ns,:),count,mpi_read,lftrank,lfttag, &
	bufferin,count,mpi_read,rgtrank,lfttag, &
	comm,status,ierr)

!	if(rank/sizex+1 .ne. sizey) by(:,nt,:)=bufferin
	if(modulo(rank,sizex*sizey)/sizex .ne. sizey-1) by(:,nt,:)=bufferin		
	
	!----------------- z field
	
	!send up and recv from below
	
	call MPI_SendRecv(bz(:,ls,:),count,mpi_read,rgtrank,rgttag, &
	bufferin,count,mpi_read,lftrank,rgttag, &
	comm,status,ierr)

!	if(rank/sizex+1 .ne. 1)bz(:,lt,:)=bufferin
	if(modulo(rank,sizex*sizey)/sizex .ne. 0) bz(:,lt,:)=bufferin	
	
	!send dwn and recv from above
	
	call MPI_SendRecv(bz(:,ns,:),count,mpi_read,lftrank,lfttag, &
	bufferin,count,mpi_read,rgtrank,lfttag, &
	comm,status,ierr)

!	if(rank/sizex+1 .ne. sizey) bz(:,nt,:)=bufferin
	if(modulo(rank,sizex*sizey)/sizex .ne. sizey-1) bz(:,nt,:)=bufferin	

end subroutine copy_layry2


!-------------------------------------------------------------------------------
! 						subroutine copy_layrx1_opt		 
!																		
! Copies a layer from one side of the grid to the other across procs (y direction)
! 							
!-------------------------------------------------------------------------------

subroutine copy_layrx1_opt(bx,by,bz,mx,my,mz,lt,ls,nt,ns)

	implicit none
	
	! dummy variables

	integer, intent(in) ::mx,my,mz,lt,ls,nt,ns
	real(sprec), dimension(:,:,:), intent(inout) :: bx, by, bz
	
	! local variables
	
	integer :: plusrank, minusrank, iperiodic,comm,count,plustag,minustag,ierr,status(statsize)
	integer :: request(2), status1(statsize,2)
	real bufferin(my,mz,3), bufferout(my,mz,3)

	plustag=100
	minustag=200
	comm=MPI_Comm_world
	
	plusrank=(rank/sizex)*sizex + modulo(rank+1,sizex)
	minusrank=(rank/sizex)*sizex + modulo(rank-1,sizex)
	
	count=my*mz*3

	!send up and recv from below
	
	bufferout(:,:,1)=bx(ls,:,:)
	bufferout(:,:,2)=by(ls,:,:)
	bufferout(:,:,3)=bz(ls,:,:)

	call MPI_SendRecv(bufferout,count,mpi_read,plusrank,plustag, &
	bufferin,count,mpi_read,minusrank,plustag, &
	comm,status,ierr)

!	if(modulo(rank,sizey)+1 .ne. 1) then 
!	if(rank/sizex+1 .ne. 1) then
	   bx(lt,:,:)=bufferin(:,:,1)
	   by(lt,:,:)=bufferin(:,:,2)
	   bz(lt,:,:)=bufferin(:,:,3)
!	endif

	!send dwn and recv from above

	bufferout(:,:,1)=bx(ns,:,:)
	bufferout(:,:,2)=by(ns,:,:)
	bufferout(:,:,3)=bz(ns,:,:)

	call MPI_SendRecv(bufferout,count,mpi_read,minusrank,minustag, &
	bufferin,count,mpi_read,plusrank,minustag, &
	comm,status,ierr)

!	if (modulo(rank,sizey)+1 .ne. sizey) then
!	if(rank/sizex+1 .ne. sizey) then
	   bx(nt,:,:)=bufferin(:,:,1)
	   by(nt,:,:)=bufferin(:,:,2)
	   bz(nt,:,:)=bufferin(:,:,3)
!	endif

end subroutine copy_layrx1_opt


!-------------------------------------------------------------------------------
! 						subroutine copy_layrx2_opt		 
!																		
! Copies a layer from one side of the grid to the other across procs (y direction)
! 							
!-------------------------------------------------------------------------------

subroutine copy_layrx2_opt(bx,by,bz,mx,my,mz,lt,ls,nt,ns)

	implicit none
	
	! dummy variables

	integer, intent(in) ::mx,my,mz,lt,ls,nt,ns
	real(sprec), dimension(:,:,:), intent(inout) :: bx, by, bz
	
	! local variables
	
	integer :: plusrank, minusrank, iperiodic,comm,count,plustag,minustag,ierr,status(statsize)
	integer :: request(2), status1(statsize,2)
	real bufferin(my,mz,3), bufferout(my,mz,3)

	plustag=100
	minustag=200
	comm=MPI_Comm_world
	
	plusrank=(rank/sizex)*sizex + modulo(rank+1,sizex)
	minusrank=(rank/sizex)*sizex + modulo(rank-1,sizex)
	
	count=my*mz*3

	!send up and recv from below
	
	bufferout(:,:,1)=bx(ls,:,:)
	bufferout(:,:,2)=by(ls,:,:)
	bufferout(:,:,3)=bz(ls,:,:)

	call MPI_SendRecv(bufferout,count,mpi_read,plusrank,plustag, &
	bufferin,count,mpi_read,minusrank,plustag, &
	comm,status,ierr)

	if(modulo(rank,sizex)+1 .ne. 1) then
	   bx(lt,:,:)=bufferin(:,:,1)
	   by(lt,:,:)=bufferin(:,:,2)
	   bz(lt,:,:)=bufferin(:,:,3)
	endif

	!send dwn and recv from above

	bufferout(:,:,1)=bx(ns,:,:)
	bufferout(:,:,2)=by(ns,:,:)
	bufferout(:,:,3)=bz(ns,:,:)

	call MPI_SendRecv(bufferout,count,mpi_read,minusrank,minustag, &
	bufferin,count,mpi_read,plusrank,minustag, &
	comm,status,ierr)

	if(modulo(rank,sizex)+1 .ne. sizex) then
	   bx(nt,:,:)=bufferin(:,:,1)
	   by(nt,:,:)=bufferin(:,:,2)
	   bz(nt,:,:)=bufferin(:,:,3)
	endif

end subroutine copy_layrx2_opt



!-------------------------------------------------------------------------------
! 						subroutine copy_layry1_opt		 
!																		
! Copies a layer from one side of the grid to the other across procs (y direction)
! 							
!-------------------------------------------------------------------------------

subroutine copy_layry1_opt(bx,by,bz,mx,my,mz,lt,ls,nt,ns)

	implicit none
	
	! dummy variables

	integer, intent(in) ::mx,my,mz,lt,ls,nt,ns
	real(sprec), dimension(:,:,:), intent(inout) :: bx, by, bz
	
	! local variables
	
	integer :: rgtrank, lftrank, iperiodic,comm,count,rgttag,lfttag,ierr,status(statsize)
	integer :: request(2), status1(statsize,2)
	real bufferin(mx,mz,3), bufferout(mx,mz,3)

	rgttag=100
	lfttag=200
	comm=MPI_Comm_world
	
!	rgtrank=(rank/sizey)*sizey + modulo(rank+1,sizey) !modulo(rank+1,size0)
!	lftrank=(rank/sizey)*sizey + modulo(rank-1,sizey) !modulo(rank-1,size0)
	
	rgtrank=modulo(rank/sizex + 1,sizey)*sizex+rank/(sizex*sizey)*(sizex*sizey) &
				+ modulo(rank,sizex) 
	lftrank=modulo(rank/sizex - 1,sizey)*sizex+rank/(sizex*sizey)*(sizex*sizey) &
				+ modulo(rank,sizex) 
	
	count=mx*mz*3

	!send up and recv from below
	
	bufferout(:,:,1)=bx(:,ls,:)
	bufferout(:,:,2)=by(:,ls,:)
	bufferout(:,:,3)=bz(:,ls,:)

	call MPI_SendRecv(bufferout,count,mpi_read,rgtrank,rgttag, &
	bufferin,count,mpi_read,lftrank,rgttag, &
	comm,status,ierr)

!	if(modulo(rank,sizey)+1 .ne. 1) then 
!	if(rank/sizex+1 .ne. 1) then
	   bx(:,lt,:)=bufferin(:,:,1)
	   by(:,lt,:)=bufferin(:,:,2)
	   bz(:,lt,:)=bufferin(:,:,3)
!	endif

	!send dwn and recv from above

	bufferout(:,:,1)=bx(:,ns,:)
	bufferout(:,:,2)=by(:,ns,:)
	bufferout(:,:,3)=bz(:,ns,:)

	call MPI_SendRecv(bufferout,count,mpi_read,lftrank,lfttag, &
	bufferin,count,mpi_read,rgtrank,lfttag, &
	comm,status,ierr)

!	if (modulo(rank,sizey)+1 .ne. sizey) then
!	if(rank/sizex+1 .ne. sizey) then
	   bx(:,nt,:)=bufferin(:,:,1)
	   by(:,nt,:)=bufferin(:,:,2)
	   bz(:,nt,:)=bufferin(:,:,3)
!	endif

end subroutine copy_layry1_opt


!-------------------------------------------------------------------------------
! 						subroutine copy_layry2_opt		 
!																		
! Copies a layer from one side of the grid to the other across procs (y direction)
! 							
!-------------------------------------------------------------------------------

subroutine copy_layry2_opt(bx,by,bz,mx,my,mz,lt,ls,nt,ns)

	implicit none
	
	! dummy variables

	integer, intent(in) ::mx,my,mz,lt,ls,nt,ns
	real(sprec), dimension(:,:,:), intent(inout) :: bx, by, bz
	
	! local variables
	
	integer :: rgtrank, lftrank, iperiodic,comm,count,rgttag,lfttag,ierr,status(statsize)
	integer :: request(2), status1(statsize,2)
	real bufferin(mx,mz,3), bufferout(mx,mz,3)

	rgttag=100
	lfttag=200
	comm=MPI_Comm_world
	
!	rgtrank=(rank/sizey)*sizey + modulo(rank+1,sizey) !modulo(rank+1,size0)
!	lftrank=(rank/sizey)*sizey + modulo(rank-1,sizey) !modulo(rank-1,size0)
	
!	rgtrank=modulo((rank/sizex + 1),sizey)*sizex + modulo(rank,sizex) 
!	lftrank=modulo((rank/sizex - 1),sizey)*sizex + modulo(rank,sizex)

!	rgtrank=modulo((rank/sizex + 1),sizey)*sizex + modulo(rank,sizex) 
!	lftrank=modulo((rank/sizex - 1),sizey)*sizex + modulo(rank,sizex)

	rgtrank=modulo(rank/sizex + 1,sizey)*sizex+rank/(sizex*sizey)*(sizex*sizey) &
				+ modulo(rank,sizex) 
	lftrank=modulo(rank/sizex - 1,sizey)*sizex+rank/(sizex*sizey)*(sizex*sizey) &
				+ modulo(rank,sizex) 
		
	count=mx*mz*3

	!send up and recv from below
	
	bufferout(:,:,1)=bx(:,ls,:)
	bufferout(:,:,2)=by(:,ls,:)
	bufferout(:,:,3)=bz(:,ls,:) 

	call MPI_SendRecv(bufferout,count,mpi_read,rgtrank,rgttag, &
	bufferin,count,mpi_read,lftrank,rgttag, &
	comm,status,ierr)

!	if(modulo(rank,sizey)+1 .ne. 1) then 
!	if(rank/sizex+1 .ne. 1) then
	if(modulo(rank,sizex*sizey)/sizex .ne. 0) then
	   bx(:,lt,:)=bufferin(:,:,1)
	   by(:,lt,:)=bufferin(:,:,2)
	   bz(:,lt,:)=bufferin(:,:,3)
	endif

	!send dwn and recv from above

	bufferout(:,:,1)=bx(:,ns,:)
	bufferout(:,:,2)=by(:,ns,:)
	bufferout(:,:,3)=bz(:,ns,:)

	call MPI_SendRecv(bufferout,count,mpi_read,lftrank,lfttag, &
	bufferin,count,mpi_read,rgtrank,lfttag, &
	comm,status,ierr)

!	if (modulo(rank,sizey)+1 .ne. sizey) then
!	if(rank/sizex+1 .ne. sizey) then
	if(modulo(rank,sizex*sizey)/sizex .ne. sizey-1) then
	   bx(:,nt,:)=bufferin(:,:,1)
	   by(:,nt,:)=bufferin(:,:,2)
	   bz(:,nt,:)=bufferin(:,:,3)
	endif

end subroutine copy_layry2_opt


!-------------------------------------------------------------------------------
! 						subroutine copy_layrz1					 
!																		
! Copies a layer from one side of the grid to the other across procs (z direction)
! (using sendrecv)				
!-------------------------------------------------------------------------------

subroutine copy_layrz1(bx,by,bz,mx,my,mz,lt,ls,nt,ns)

	implicit none
	
	! dummy variables

	real(sprec), dimension(:,:,:), intent(inout) :: bx, by, bz
	integer, intent(in) ::mx,my,mz,lt,ls,nt,ns

	! local variables
	
	integer ::uprank, dwnrank, iperiodic,comm,count,uptag,dwntag,ierr,status(statsize)
	integer ::request(2), status1(statsize,2)


	uptag=100
	dwntag=200
	comm=MPI_Comm_world
	
#ifndef twoD 
	
		!if twoD don't exchange anything

		uprank=modulo(rank/(sizex*sizey) + 1,sizez)*(sizex*sizey) + & 
		 		  modulo(rank,sizex*sizey) 
		dwnrank=modulo(rank/(sizex*sizey) - 1,sizez)*(sizex*sizey) + & 
			      modulo(rank,sizex*sizey)
		
		count=mx*my

		!----------------- x field
		
		!send up and recv from below
		
		call MPI_SendRecv(bx(:,:,ls),count,mpi_read,uprank,uptag, &
		bx(:,:,lt),   count,mpi_read,dwnrank,uptag, &
		comm,status,ierr)
		
		!send dwn and recv from above
		
		call MPI_SendRecv(bx(:,:,ns),count,mpi_read,dwnrank,dwntag, &
		bx(:,:,nt),count,mpi_read,uprank,dwntag, &
		comm,status,ierr)
		
		!----------------- y field
		
		!send up and recv from below
		
		call MPI_SendRecv(by(:,:,ls),count,mpi_read,uprank,uptag, &
		by(:,:,lt),   count,mpi_read,dwnrank,uptag, &
		comm,status,ierr)
		
		!send dwn and recv from above
		
		call MPI_SendRecv(by(:,:,ns),count,mpi_read,dwnrank,dwntag, &
		by(:,:,nt),count,mpi_read,uprank,dwntag, &
		comm,status,ierr)
		
		!----------------- z field
		
		!send up and recv from below
		
		call MPI_SendRecv(bz(:,:,ls),count,mpi_read,uprank,uptag, &
		bz(:,:,lt),   count,mpi_read,dwnrank,uptag, &
		comm,status,ierr)
		
		!send dwn and recv from above
		
		call MPI_SendRecv(bz(:,:,ns),count,mpi_read,dwnrank,dwntag, &
		bz(:,:,nt),count,mpi_read,uprank,dwntag, &
		comm,status,ierr)
		
#endif

end subroutine copy_layrz1



!-------------------------------------------------------------------------------
! 						subroutine copy_layrz2					 
!																		
! Copies a layer from one side of the grid to the other across procs (z direction)
! (using sendrecv)											
!-------------------------------------------------------------------------------

subroutine copy_layrz2(bx,by,bz,mx,my,mz,lt,ls,nt,ns)

	implicit none
	
	! dummy variables

	real(sprec), dimension(:,:,:), intent(inout) :: bx, by, bz
	integer, intent(in) ::mx,my,mz,lt,ls,nt,ns
	
	! local variables
	
	integer ::uprank, dwnrank, iperiodic,comm,count,uptag,dwntag,ierr,status(statsize)
	integer ::request(2), status1(statsize,2)
	real bufferin(mx,my)!,bufferout(mx,my)

#ifndef twoD 
	uptag=100
	dwntag=200
	comm=MPI_Comm_world
	
	uprank=modulo(rank/(sizex*sizey) + 1,sizez)*(sizex*sizey) + & 
		 		  modulo(rank,sizex*sizey) 
	dwnrank=modulo(rank/(sizex*sizey) - 1,sizez)*(sizex*sizey) + & 
			      modulo(rank,sizex*sizey)
	 
	
	count=mx*my

	!----------------- x field
	
	!send up and recv from below
	
	call MPI_SendRecv(bx(:,:,ls),count,mpi_read,uprank,uptag, &
	bufferin,count,mpi_read,dwnrank,uptag, &
	comm,status,ierr)
	
	if(rank/(sizex*sizey) .ne. 0) bx(:,:,lt)=bufferin
	
	!send dwn and recv from above
	
	call MPI_SendRecv(bx(:,:,ns),count,mpi_read,dwnrank,dwntag, &
	bufferin,count,mpi_read,uprank,dwntag, &
	comm,status,ierr)
	
	if(rank/(sizex*sizey) .ne. sizez-1) bx(:,:,nt)=bufferin
	
	!----------------- y field
	
	!send up and recv from below
	
	call MPI_SendRecv(by(:,:,mz-3),count,mpi_read,uprank,uptag, &
	bufferin,   count,mpi_read,dwnrank,uptag, &
	comm,status,ierr)
	
	if(rank/(sizex*sizey) .ne. 0) by(:,:,2)=bufferin
	
	!send dwn and recv from above
	
	call MPI_SendRecv(by(:,:,ns),count,mpi_read,dwnrank,dwntag, &
	bufferin,count,mpi_read,uprank,dwntag, &
	comm,status,ierr)
	
	if(rank/(sizex*sizey) .ne. sizez-1) by(:,:,nt)=bufferin
	
	!----------------- z field
	
	!send up and recv from below
	
	call MPI_SendRecv(bz(:,:,ls),count,mpi_read,uprank,uptag, &
	bufferin,   count,mpi_read,dwnrank,uptag, &
	comm,status,ierr)
	
	if(rank/(sizex*sizey) .ne. 0) bz(:,:,lt)=bufferin
	
	!send dwn and recv from above
	
	call MPI_SendRecv(bz(:,:,ns),count,mpi_read,dwnrank,dwntag, &
	bufferin,count,mpi_read,uprank,dwntag, &
	comm,status,ierr)

	if(rank/(sizex*sizey) .ne. sizez-1) bz(:,:,nt)=bufferin
#endif
end subroutine copy_layrz2



!-------------------------------------------------------------------------------
! 						subroutine copy_layry1_opt		 
!																		
! Copies a layer from one side of the grid to the other across procs (y direction)
! 							
!-------------------------------------------------------------------------------

subroutine copy_layrz1_opt(bx,by,bz,mx,my,mz,lt,ls,nt,ns)

	implicit none
	
	! dummy variables

	integer, intent(in) ::mx,my,mz,lt,ls,nt,ns
	real(sprec), dimension(:,:,:), intent(inout) :: bx, by, bz
	
	! local variables
	
	integer ::uprank, dwnrank, iperiodic,comm,count,uptag,dwntag,ierr,status(statsize)
	integer ::request(2), status1(statsize,2)
	real bufferin(mx,my,3), bufferout(mx,my,3)

	uptag=100
	dwntag=200
	comm=MPI_Comm_world
#ifndef twoD 	
	
!	rgtrank=(rank/sizey)*sizey + modulo(rank+1,sizey) !modulo(rank+1,size0)
!	lftrank=(rank/sizey)*sizey + modulo(rank-1,sizey) !modulo(rank-1,size0)
	
	uprank=modulo((rank/(sizex*sizey) + 1),sizez)*(sizex*sizey) + & 
		 		  modulo(rank,sizex*sizey) 
	dwnrank=modulo((rank/(sizex*sizey) - 1),sizez)*(sizex*sizey) + & 
			      modulo(rank,sizex*sizey) 
	
	count=mx*my*3

	!send up and recv from below
	
	bufferout(:,:,1)=bx(:,:,ls)
	bufferout(:,:,2)=by(:,:,ls)
	bufferout(:,:,3)=bz(:,:,ls)

	call MPI_SendRecv(bufferout,count,mpi_read,uprank,uptag, &
	bufferin,count,mpi_read,dwnrank,uptag, &
	comm,status,ierr)

!	if(modulo(rank,sizey)+1 .ne. 1) then 
!	if(rank/sizex+1 .ne. 1) then
	   bx(:,:,lt)=bufferin(:,:,1)
	   by(:,:,lt)=bufferin(:,:,2)
	   bz(:,:,lt)=bufferin(:,:,3)
!	endif

	!send dwn and recv from above

	bufferout(:,:,1)=bx(:,:,ns)
	bufferout(:,:,2)=by(:,:,ns)
	bufferout(:,:,3)=bz(:,:,ns)

	call MPI_SendRecv(bufferout,count,mpi_read,dwnrank,dwntag, &
	bufferin,count,mpi_read,uprank,dwntag, &
	comm,status,ierr)

!	if (modulo(rank,sizey)+1 .ne. sizey) then
!	if(rank/sizex+1 .ne. sizey) then
	   bx(:,:,nt)=bufferin(:,:,1)
	   by(:,:,nt)=bufferin(:,:,2)
	   bz(:,:,nt)=bufferin(:,:,3)
!	endif
#endif
end subroutine copy_layrz1_opt


!-------------------------------------------------------------------------------
! 						subroutine copy_layry2_opt		 
!																		
! Copies a layer from one side of the grid to the other across procs (y direction)
! 							
!-------------------------------------------------------------------------------

subroutine copy_layrz2_opt(bx,by,bz,mx,my,mz,lt,ls,nt,ns)

	implicit none
	
	! dummy variables

	integer, intent(in) ::mx,my,mz,lt,ls,nt,ns
	real(sprec), dimension(:,:,:), intent(inout) :: bx, by, bz
	
	! local variables
	
	integer ::uprank, dwnrank, iperiodic,comm,count,uptag,dwntag,ierr,status(statsize)
	integer ::request(2), status1(statsize,2)
	real bufferin(mx,my,3), bufferout(mx,my,3)

	uptag=100
	dwntag=200
	comm=MPI_Comm_world

#ifndef twoD 	
!	rgtrank=(rank/sizey)*sizey + modulo(rank+1,sizey) !modulo(rank+1,size0)
!	lftrank=(rank/sizey)*sizey + modulo(rank-1,sizey) !modulo(rank-1,size0)
	
	
	uprank=modulo((rank/(sizex*sizey) + 1),sizez)*(sizex*sizey) + & 
		 		  modulo(rank,sizex*sizey) 
	dwnrank=modulo((rank/(sizex*sizey) - 1),sizez)*(sizex*sizey) + & 
			      modulo(rank,sizex*sizey) 
	
	count=mx*my*3


	!send up and recv from below
	
	bufferout(:,:,1)=bx(:,:,ls)
	bufferout(:,:,2)=by(:,:,ls)
	bufferout(:,:,3)=bz(:,:,ls)

	call MPI_SendRecv(bufferout,count,mpi_read,uprank,uptag, &
	bufferin,count,mpi_read,dwnrank,uptag, &
	comm,status,ierr)

!	if(modulo(rank,sizey)+1 .ne. 1) then 
	if(rank/(sizex*sizey) .ne. 0) then
	   bx(:,:,lt)=bufferin(:,:,1)
	   by(:,:,lt)=bufferin(:,:,2)
	   bz(:,:,lt)=bufferin(:,:,3)
	endif

	!send dwn and recv from above

	bufferout(:,:,1)=bx(:,:,ns)
	bufferout(:,:,2)=by(:,:,ns)
	bufferout(:,:,3)=bz(:,:,ns)

	call MPI_SendRecv(bufferout,count,mpi_read,dwnrank,dwntag, &
	bufferin,count,mpi_read,uprank,dwntag, &
	comm,status,ierr)

!	if (modulo(rank,sizey)+1 .ne. sizey) then
	if(rank/(sizex*sizey) .ne. sizez-1) then
	   bx(:,:,nt)=bufferin(:,:,1)
	   by(:,:,nt)=bufferin(:,:,2)
	   bz(:,:,nt)=bufferin(:,:,3)
	endif
#endif
end subroutine copy_layrz2_opt




!-------------------------------------------------------------------------------
! 						subroutine exchange_current					 
!																		
! 
!  							
!-------------------------------------------------------------------------------

subroutine exchange_current()

	implicit none

	! local variables
	
	integer :: uprank, dwnrank, comm,count,uptag,dwntag,ierr, &
			   status(statsize), rgtrank, lftrank, rgttag, lfttag, &
			   plusrank, minusrank, plustag, minustag
	integer :: request(2), status1(statsize,2)
	logical iperiodic
	integer :: i,j,k, k1,k2, j1,j2, i1, i2

	uptag=100
	dwntag=200
	rgttag=uptag
	lfttag=dwntag
	plustag=uptag
	minustag=dwntag
	
	comm=MPI_Comm_world
	
	
	if(sizex .eq. 1) then ! no communication is required
			
	!add periodic from x

	 if(periodicx .eq. 1) then
	 		 
		bufferin1x(1:2,:,:)=curx(mx-2:mx-1,:,:)
		bufferin2x(1:2,:,:)=curx(1:2,:,:)
		curx(3:4,:,:)=curx(3:4,:,:)+bufferin1x(1:2,:,:) 
		curx(mx-4:mx-3,:,:)=curx(mx-4:mx-3,:,:)+bufferin2x(1:2,:,:) 
		
		bufferin1x(1:2,:,:)=cury(mx-2:mx-1,:,:)
		bufferin2x(1:2,:,:)=cury(1:2,:,:)
		cury(3:4,:,:)=cury(3:4,:,:)+bufferin1x(1:2,:,:) 
		cury(mx-4:mx-3,:,:)=cury(mx-4:mx-3,:,:)+bufferin2x(1:2,:,:) 
		
		bufferin1x(1:2,:,:)=curz(mx-2:mx-1,:,:) 
		bufferin2x(1:2,:,:)=curz(1:2,:,:)
		curz(3:4,:,:)=curz(3:4,:,:)+bufferin1x(1:2,:,:)
		curz(mx-4:mx-3,:,:)=curz(mx-4:mx-3,:,:)+bufferin2x(1:2,:,:)

	 endif			
	
	else	! if sizex > 1, comminucation is required
	
	!add current from x

		plusrank=(rank/sizex)*sizex + modulo(rank+1,sizex)
		minusrank=(rank/sizex)*sizex + modulo(rank-1,sizex)

		k1=1 !3
		k2=mz !mz-3

		i1=1 !3
		i2=mx !mx-3

		j1=1 !3 
		j2=my !my-3 
	
#ifdef twoD
		count=(j2-j1+1)*2
		k1=1
		k2=1
#else
		count=(j2-j1+1)*(k2-k1+1)*2
#endif

		call MPI_SendRecv(curx(mx-2:mx-1,j1:j2,k1:k2),count,mpi_read,plusrank &
		,plustag,bufferin1x(1:2,j1:j2,k1:k2), count,mpi_read,minusrank,plustag, &
		comm,status,ierr)

		iperiodic=.true.
		if(modulo(rank,sizex).eq.0 .and. periodicx .ne. 1) iperiodic=.false.
		if(iperiodic) curx(3:4,j1:j2,k1:k2)=curx(3:4,j1:j2,k1:k2)+bufferin1x(1:2,j1:j2,k1:k2)
			
		!send lft and recv from rgt
	
		call MPI_SendRecv(curx(1:2,j1:j2,k1:k2),count,mpi_read,minusrank,minustag, &
		bufferin2x(1:2,j1:j2,k1:k2),count,mpi_read,plusrank,minustag, comm &
		,status ,ierr)
	
		iperiodic=.true.
		if(modulo(rank,sizex).eq.sizex-1 .and. periodicx .ne. 1) iperiodic=.false.
		if(iperiodic) curx(mx-4:mx-3,j1:j2,k1:k2)=curx(mx-4:mx-3,j1:j2,k1:k2)+bufferin2x(1:2,j1:j2,k1:k2)
	
	
		!----------------- y current
	
		!send rgt and recv from lft

		if(debug) print *, "rankminus=",rank     

		call MPI_SendRecv(cury(mx-2:mx-1,j1:j2,k1:k2),count,mpi_read,plusrank &
		,plustag,bufferin1x(1:2,j1:j2,k1:k2), count,mpi_read,minusrank,plustag, &
		comm,status,ierr)

		iperiodic=.true.
		if(modulo(rank,sizex).eq.0 .and. periodicx .ne. 1) iperiodic=.false.
		if(iperiodic) cury(3:4,j1:j2,k1:k2)=cury(3:4,j1:j2,k1:k2)+bufferin1x(1:2,j1:j2,k1:k2)
			
		!send lft and recv from rgt
	
		call MPI_SendRecv(cury(1:2,j1:j2,k1:k2),count,mpi_read,minusrank,minustag, &
		bufferin2x(1:2,j1:j2,k1:k2),count,mpi_read,plusrank,minustag, comm &
		,status ,ierr)
	
		iperiodic=.true.
		if(modulo(rank,sizex).eq.sizex-1 .and. periodicx .ne. 1) iperiodic=.false.
		if(iperiodic) cury(mx-4:mx-3,j1:j2,k1:k2)=cury(mx-4:mx-3,j1:j2,k1:k2)+bufferin2x(1:2,j1:j2,k1:k2)
	
	
		!----------------- z current

		if(debug) print *, "rank plus=",rank     

		call MPI_SendRecv(curz(mx-2:mx-1,j1:j2,k1:k2),count,mpi_read,plusrank &
		,plustag,bufferin1x(1:2,j1:j2,k1:k2), count,mpi_read,minusrank,plustag, &
		comm,status,ierr)

		iperiodic=.true.
		if(modulo(rank,sizex).eq.0 .and. periodicx .ne. 1) iperiodic=.false.
		if(iperiodic) curz(3:4,j1:j2,k1:k2)=curz(3:4,j1:j2,k1:k2)+bufferin1x(1:2,j1:j2,k1:k2)
			
		!send lft and recv from rgt
	
		call MPI_SendRecv(curz(1:2,j1:j2,k1:k2),count,mpi_read,minusrank,minustag, &
		bufferin2x(1:2,j1:j2,k1:k2),count,mpi_read,plusrank,minustag, comm &
		,status ,ierr)
	
		iperiodic=.true.
		if(modulo(rank,sizex).eq.sizex-1 .and. periodicx .ne. 1) iperiodic=.false.
		if(iperiodic) curz(mx-4:mx-3,j1:j2,k1:k2)=curz(mx-4:mx-3,j1:j2,k1:k2)+bufferin2x(1:2,j1:j2,k1:k2)
	
	endif ! if sizex .ne. 1
	
	
	if(sizey .eq. 1) then ! no communication is required
			
	!add periodic from x

	 if(periodicy .eq. 1) then
!	 		 
		bufferin1y(:,1:2,:)=curx(:,my-2:my-1,:)
		bufferin2y(:,1:2,:)=curx(:,1:2,:)
		curx(:,3:4,:)=curx(:,3:4,:)+bufferin1y(:,1:2,:) 
		curx(:,my-4:my-3,:)=curx(:,my-4:my-3,:)+bufferin2y(:,1:2,:) 
		
		bufferin1y(:,1:2,:)=cury(:,my-2:my-1,:)
		bufferin2y(:,1:2,:)=cury(:,1:2,:)
		cury(:,3:4,:)=cury(:,3:4,:)+bufferin1y(:,1:2,:) 
		cury(:,my-4:my-3,:)=cury(:,my-4:my-3,:)+bufferin2y(:,1:2,:) 
		
		bufferin1y(:,1:2,:)=curz(:,my-2:my-1,:) 
		bufferin2y(:,1:2,:)=curz(:,1:2,:)
		curz(:,3:4,:)=curz(:,3:4,:)+bufferin1y(:,1:2,:)
		curz(:,my-4:my-3,:)=curz(:,my-4:my-3,:)+bufferin2y(:,1:2,:)

	 endif		
	
	else	! communication is required
			
	!add current from y

!	rgtrank=modulo((rank/sizex + 1),sizey)*sizex + modulo(rank,sizex) 
!	lftrank=modulo((rank/sizex - 1),sizey)*sizex + modulo(rank,sizex)

		rgtrank=modulo(rank/sizex + 1,sizey)*sizex+rank/(sizex*sizey)*(sizex*sizey) &
				+ modulo(rank,sizex) 
		lftrank=modulo(rank/sizex - 1,sizey)*sizex+rank/(sizex*sizey)*(sizex*sizey) &
				+ modulo(rank,sizex) 		
		
		k1=1 !3
		k2=mz !mz-3

		i1=1 !3
		i2=mx !mx-3

		j1=1 !3 
		j2=my !my-3 

!hack
!	k1=1
!	k2=mz
!	j1=1
!	j2=my
!	i1=1
!	i2=mx
!end hack

		
#ifdef twoD
		count=(i2-i1+1)*2
		k1=1
		k2=1
#else
		count=(i2-i1+1)*(k2-k1+1)*2
#endif

		call MPI_SendRecv(curx(i1:i2,my-2:my-1,k1:k2),count,mpi_read,rgtrank &
		,rgttag,bufferin1y(i1:i2,1:2,k1:k2), count,mpi_read,lftrank,rgttag, &
		comm,status,ierr)

		iperiodic=.true.
!		if(rank/sizex.eq.0 .and. periodicy .ne. 1) iperiodic=.false.
		if(modulo(rank,sizex*sizey)/sizex .eq. 0 .and. periodicy .ne. 1) iperiodic=.false.		
		if(iperiodic) curx(i1:i2,3:4,k1:k2)=curx(i1:i2,3:4,k1:k2)+bufferin1y(i1:i2,1:2,k1:k2)
			
		!send lft and recv from rgt
	
		call MPI_SendRecv(curx(i1:i2,1:2,k1:k2),count,mpi_read,lftrank,lfttag, &
		bufferin2y(i1:i2,1:2,k1:k2),count,mpi_read,rgtrank,lfttag, comm &
		,status ,ierr)
	
		iperiodic=.true.
!		if(rank/sizex.eq.sizey-1 .and. periodicy .ne. 1) iperiodic=.false.
		if(modulo(rank,sizex*sizey)/sizex .eq. sizey-1 .and. periodicy .ne. 1) iperiodic=.false.
		if(iperiodic) curx(i1:i2,my-4:my-3,k1:k2)=curx(i1:i2,my-4:my-3,k1:k2)+bufferin2y(i1:i2,1:2,k1:k2)
	
		!----------------- y current
	
		!send rgt and recv from lft

		if(debug) print *, "rankleft=",rank     

		call MPI_SendRecv(cury(i1:i2,my-2:my-1,k1:k2),count,mpi_read,rgtrank &
		,rgttag,bufferin1y(i1:i2,1:2,k1:k2), count,mpi_read,lftrank,rgttag, &
		comm,status,ierr)
		
		iperiodic=.true.
!		if(rank/sizex.eq.0 .and. periodicy .ne. 1) iperiodic=.false.
		if(modulo(rank,sizex*sizey)/sizex .eq. 0 .and. periodicy .ne. 1) iperiodic=.false.		
		if(iperiodic) cury(i1:i2,3:4,k1:k2)=cury(i1:i2,3:4,k1:k2)+bufferin1y(i1:i2,1:2,k1:k2)

	
		!send lft and recv from rgt
	
		call MPI_SendRecv(cury(i1:i2,1:2,k1:k2),count,mpi_read,lftrank,lfttag, &
		bufferin2y(i1:i2,1:2,k1:k2),count,mpi_read,rgtrank,lfttag, comm,status &
		,ierr)
	
		iperiodic=.true.
!		if(rank/sizex.eq.sizey-1 .and. periodicy .ne. 1) iperiodic=.false.
		if(modulo(rank,sizex*sizey)/sizex .eq. sizey-1 .and. periodicy .ne. 1) iperiodic=.false.	
		if(iperiodic) cury(i1:i2,my-4:my-3,k1:k2)=cury(i1:i2,my-4:my-3,k1:k2)+bufferin2y(i1:i2,1:2,k1:k2)
	
		!----------------- z current

		if(debug) print *, "rank rgt=",rank     

		call MPI_SendRecv(curz(i1:i2,my-2:my-1,k1:k2),count,mpi_read,rgtrank &
		,rgttag,bufferin1y(i1:i2,1:2,k1:k2),count,mpi_read,lftrank,rgttag &
		,comm,status,ierr)

		iperiodic=.true.
!		if(rank/sizex.eq.0 .and. periodicy .ne. 1) iperiodic=.false.
		if(modulo(rank,sizex*sizey)/sizex .eq. 0 .and. periodicy .ne. 1) iperiodic=.false.	
		if(iperiodic) curz(i1:i2,3:4,k1:k2)=curz(i1:i2,3:4,k1:k2)+bufferin1y(i1:i2,1:2,k1:k2)

		!send lft and recv from rgt

		if(debug) print *, "zlleftrgt=",rank     

		call MPI_SendRecv(curz(i1:i2,1:2,k1:k2),count,mpi_read,lftrank,lfttag, &
		bufferin2y(i1:i2,1:2,k1:k2),count,mpi_read,rgtrank,lfttag, &
		comm,status,ierr)

		iperiodic=.true.
!		if(rank/sizex.eq.sizey-1 .and. periodicy .ne. 1) iperiodic=.false.
		if(modulo(rank,sizex*sizey)/sizex .eq. sizey-1 .and. periodicy .ne. 1) iperiodic=.false.		
		if(iperiodic) curz(i1:i2,my-4:my-3,k1:k2)=curz(i1:i2,my-4:my-3,k1:k2)+bufferin2y(i1:i2,1:2,k1:k2)
	
	endif	! if sizey .ne. 1
		
#ifndef twoD
		
		!add current from z

		uprank=modulo(rank/(sizex*sizey) + 1,sizez)*(sizex*sizey) + & 
		 		  modulo(rank,sizex*sizey) 
		dwnrank=modulo(rank/(sizex*sizey) - 1,sizez)*(sizex*sizey) + & 
			      modulo(rank,sizex*sizey)
	 
		k1=1 !3
		k2=mz !mz-3

		i1=1 !3
		i2=mx !mx-3

		j1=1 !3 
		j2=my !my-3 
		
		count=(i2-i1+1)*(j2-j1+1)*2
		
		!----------------- x current
		
		!send up and recv from below
		
		call MPI_SendRecv(curx(i1:i2,j1:j2,mz-2:mz-1),count,mpi_read,uprank,uptag, &
		bufferin1(i1:i2,j1:j2,1:2),count,mpi_read,dwnrank,uptag, comm &
		,status,ierr)

		iperiodic=.true.
		if(rank/(sizex*sizey) .eq. 0 .and. periodicz .ne. 1) iperiodic=.false.
		if(iperiodic)curx(i1:i2,j1:j2,3:4)=curx(i1:i2,j1:j2,3:4)+bufferin1(i1:i2,j1:j2,1:2)	

		!send dwn and recv from above
		
		call MPI_SendRecv(curx(i1:i2,j1:j2,1:2),count,mpi_read,dwnrank,dwntag, &
		bufferin2(i1:i2,j1:j2,1:2),count,mpi_read,uprank,dwntag, comm,status &
		,ierr)

		iperiodic=.true.
		if(rank/(sizex*sizey) .eq. sizez-1 .and. periodicz .ne. 1) iperiodic=.false.
		if(iperiodic) curx(i1:i2,j1:j2,mz-4:mz-3)=curx(i1:i2,j1:j2,mz-4:mz-3)+bufferin2(i1:i2,j1:j2,1:2)
		
		!----------------- y current
		
		!send up and recv from below
		
		call MPI_SendRecv(cury(i1:i2,j1:j2,mz-2:mz-1),count,mpi_read,uprank,uptag, &
		bufferin1(i1:i2,j1:j2,1:2), count,mpi_read,dwnrank,uptag, comm,status &
		,ierr)

		iperiodic=.true.
		if(rank/(sizex*sizey) .eq. 0 .and. periodicz .ne. 1) iperiodic=.false.
		if(iperiodic)cury(i1:i2,j1:j2,3:4)=cury(i1:i2,j1:j2,3:4)+bufferin1(i1:i2,j1:j2,1:2)	
		
		!send dwn and recv from above
		
		call MPI_SendRecv(cury(i1:i2,j1:j2,1:2),count,mpi_read,dwnrank,dwntag, &
		bufferin2(i1:i2,j1:j2,1:2),count,mpi_read,uprank,dwntag, comm,status &
		,ierr)

		iperiodic=.true.
		if(rank/(sizex*sizey) .eq. sizez-1 .and. periodicz .ne. 1) iperiodic=.false.
		if(iperiodic) cury(i1:i2,j1:j2,mz-4:mz-3)=cury(i1:i2,j1:j2,mz-4:mz-3)+bufferin2(i1:i2,j1:j2,1:2)
				
		!----------------- z current
		
		!send up and recv from below

		call MPI_SendRecv(curz(i1:i2,j1:j2,mz-2:mz-1),count,mpi_read,uprank,uptag, &
		bufferin1(i1:i2,j1:j2,1:2),count,mpi_read,dwnrank,uptag, &
		comm,status,ierr)

		iperiodic=.true.
		if(rank/(sizex*sizey) .eq. 0 .and. periodicz .ne. 1) iperiodic=.false.
		if(iperiodic)curz(i1:i2,j1:j2,3:4)=curz(i1:i2,j1:j2,3:4)+bufferin1(i1:i2,j1:j2,1:2)	
		
		!send dwn and recv from above
		
		call MPI_SendRecv(curz(i1:i2,j1:j2,1:2),count,mpi_read,dwnrank,dwntag, &
		bufferin2(i1:i2,j1:j2,1:2),count,mpi_read,uprank,dwntag, &
		comm,status,ierr)

		iperiodic=.true.
		if(rank/(sizex*sizey) .eq. sizez-1 .and. periodicz .ne. 1) iperiodic=.false.
		if(iperiodic) curz(i1:i2,j1:j2,mz-4:mz-3)=curz(i1:i2,j1:j2,mz-4:mz-3)+bufferin2(i1:i2,j1:j2,1:2)

#endif 

end subroutine exchange_current



!-------------------------------------------------------------------------------
! 						subroutine preledge()					 
!																		
! 
!
!-------------------------------------------------------------------------------

subroutine preledge(bx,by,bz,ex,ey,ez,ix,iy,iz,mx,my,mz,m,c)

!	implicit none

	! dummy variables

	real :: bx(1), by(1), bz(1), ex(1), ey(1), ez(1)
!	real(sprec), dimension(:), intent(inout) :: bx, by, bz, ex, ey, ez
	integer, intent(in) :: ix,iy,iz,mx,my,mz,m
	real(sprec), intent(in) :: c
	
	! local variables
	
	real(sprec) :: s, t, r, p, q, temp
	integer :: n
	
	s=.4142136
	t=c/(2.*c+1.+s)

#ifndef twoD
	do n=m+iy*(my-1)+iz*(mz-1)+ix,m+iy*(my-1)+iz*(mz-1)+ix*(mx-2),ix
		bx(n)=bx(n-iy-iz)+(1.-4.*t)*bx(n)+(1.-2.*t)*(bx(n-iy)+bx(n-iz)) &
		 +s*t*(by(n)-by(n-ix)+by(n-iz)-by(n-ix-iz) &
		 +bz(n)-bz(n-ix)+bz(n-iy)-bz(n-ix-iy))
	enddo
 
	r=4./(2.*c+2.+s)

	do n=m+iz*(mz-1),m+iz*(mz-1)+iy*(my-2),iy
		bx(n)=(1.-c*r)*(bx(n)+bx(n+ix))+bx(n-iz)+bx(n+ix-iz)-r*(bz(n) &
		 +c*((1.-s)*(ex(n+iy)-ex(n))+(1.+s)*.25*(ez(n+iy) &
		 -ez(n)+ez(n+ix+iy)-ez(n+ix)+ez(n+iy-iz)-ez(n-iz) &
		 +ez(n+ix+iy-iz)-ez(n+ix-iz))))
	enddo

	do n=m+iy*(my-1),m+iy*(my-1)+iz*(mz-2),iz
		bx(n)=(1.-c*r)*(bx(n)+bx(n+ix))+bx(n-iy)+bx(n+ix-iy) &
		 -r*(by(n)-c*((1.-s)*(ex(n+iz)-ex(n)) &
		  +(1.+s)*.25*(ey(n+iz)-ey(n)+ey(n+ix+iz) &
		 -ey(n+ix)+ey(n+iz-iy)-ey(n-iy)+ey(n+ix+iz-iy) &
		 -ey(n+ix-iy))))
	enddo
	p=(1.+c)*2./(1.+2.*c*(1.+c*s))
	q=c*s*2./(1.+2.*c*(1.+c*s))
	
	do n=m+iy*(my-1)+iz*(mz-1),m+ix*(mx-2)+iy*(my-1)+iz*(mz-1),ix
		temp=bz(n)-.5*c*(1.-s)*(ey(n+ix)-ey(n)+ey(n+ix-iy)-ey(n-iy))
		bz(n)=bz(n-iy)-bz(n)+p*temp+q*by(n)
		by(n)=by(n-iz)-by(n)+p*by(n)+q*temp
	enddo
#else

	if(ix .eq. 0) then 
!	do n=m+iy*(my-1)+iz*(mz-1)+ix,m+iy*(my-1)+iz*(mz-1)+ix*(mx-2),ix
	 n=m+iy*(my-1)+iz*(mz-1)+ix
		bx(n)=bx(n-iy-iz)+(1.-4.*t)*bx(n)+(1.-2.*t)*(bx(n-iy)+bx(n-iz)) &
		 +s*t*(by(n)-by(n-ix)+by(n-iz)-by(n-ix-iz) &
		 +bz(n)-bz(n-ix)+bz(n-iy)-bz(n-ix-iy))
!	enddo
 
	r=4./(2.*c+2.+s)

	do n=m+iz*(mz-1),m+iz*(mz-1)+iy*(my-2),iy
		bx(n)=(1.-c*r)*(bx(n)+bx(n+ix))+bx(n-iz)+bx(n+ix-iz)-r*(bz(n) &
		 +c*((1.-s)*(ex(n+iy)-ex(n))+(1.+s)*.25*(ez(n+iy) &
		 -ez(n)+ez(n+ix+iy)-ez(n+ix)+ez(n+iy-iz)-ez(n-iz) &
		 +ez(n+ix+iy-iz)-ez(n+ix-iz))))
	enddo

	do n=m+iy*(my-1),m+iy*(my-1)+iz*(mz-2),iz
		bx(n)=(1.-c*r)*(bx(n)+bx(n+ix))+bx(n-iy)+bx(n+ix-iy) &
		 -r*(by(n)-c*((1.-s)*(ex(n+iz)-ex(n)) &
		  +(1.+s)*.25*(ey(n+iz)-ey(n)+ey(n+ix+iz) &
		 -ey(n+ix)+ey(n+iz-iy)-ey(n-iy)+ey(n+ix+iz-iy) &
		 -ey(n+ix-iy))))
	enddo
	p=(1.+c)*2./(1.+2.*c*(1.+c*s))
	q=c*s*2./(1.+2.*c*(1.+c*s))

!	do n=m+iy*(my-1)+iz*(mz-1),m+ix*(mx-2)+iy*(my-1)+iz*(mz-1),ix
	 n=m+iy*(my-1)+iz*(mz-1)
		temp=bz(n)-.5*c*(1.-s)*(ey(n+ix)-ey(n)+ey(n+ix-iy)-ey(n-iy))
		bz(n)=bz(n-iy)-bz(n)+p*temp+q*by(n)
		by(n)=by(n-iz)-by(n)+p*by(n)+q*temp
!	enddo

	endif
	if(iy .eq. 0) then 
	do n=m+iy*(my-1)+iz*(mz-1)+ix,m+iy*(my-1)+iz*(mz-1)+ix*(mx-2),ix
		bx(n)=bx(n-iy-iz)+(1.-4.*t)*bx(n)+(1.-2.*t)*(bx(n-iy)+bx(n-iz)) &
		 +s*t*(by(n)-by(n-ix)+by(n-iz)-by(n-ix-iz) &
		 +bz(n)-bz(n-ix)+bz(n-iy)-bz(n-ix-iy))
	enddo
 
	r=4./(2.*c+2.+s)
!	do n=m+iz*(mz-1),m+iz*(mz-1)+iy*(my-2),iy
	 n=m+iz*(mz-1)
		bx(n)=(1.-c*r)*(bx(n)+bx(n+ix))+bx(n-iz)+bx(n+ix-iz)-r*(bz(n) &
		 +c*((1.-s)*(ex(n+iy)-ex(n))+(1.+s)*.25*(ez(n+iy) &
		 -ez(n)+ez(n+ix+iy)-ez(n+ix)+ez(n+iy-iz)-ez(n-iz) &
		 +ez(n+ix+iy-iz)-ez(n+ix-iz))))
!	enddo

	do n=m+iy*(my-1),m+iy*(my-1)+iz*(mz-2),iz
		bx(n)=(1.-c*r)*(bx(n)+bx(n+ix))+bx(n-iy)+bx(n+ix-iy) &
		 -r*(by(n)-c*((1.-s)*(ex(n+iz)-ex(n)) &
		  +(1.+s)*.25*(ey(n+iz)-ey(n)+ey(n+ix+iz) &
		 -ey(n+ix)+ey(n+iz-iy)-ey(n-iy)+ey(n+ix+iz-iy) &
		 -ey(n+ix-iy))))
	enddo
	p=(1.+c)*2./(1.+2.*c*(1.+c*s))
	q=c*s*2./(1.+2.*c*(1.+c*s))
	
	do n=m+iy*(my-1)+iz*(mz-1),m+ix*(mx-2)+iy*(my-1)+iz*(mz-1),ix
		temp=bz(n)-.5*c*(1.-s)*(ey(n+ix)-ey(n)+ey(n+ix-iy)-ey(n-iy))
		bz(n)=bz(n-iy)-bz(n)+p*temp+q*by(n)
		by(n)=by(n-iz)-by(n)+p*by(n)+q*temp
	enddo

	endif
	if(iz .eq. 0) then 

	do n=m+iy*(my-1)+iz*(mz-1)+ix,m+iy*(my-1)+iz*(mz-1)+ix*(mx-2),ix
		bx(n)=bx(n-iy-iz)+(1.-4.*t)*bx(n)+(1.-2.*t)*(bx(n-iy)+bx(n-iz)) &
		 +s*t*(by(n)-by(n-ix)+by(n-iz)-by(n-ix-iz) &
		 +bz(n)-bz(n-ix)+bz(n-iy)-bz(n-ix-iy))
	enddo
 
	r=4./(2.*c+2.+s)

	do n=m+iz*(mz-1),m+iz*(mz-1)+iy*(my-2),iy
		bx(n)=(1.-c*r)*(bx(n)+bx(n+ix))+bx(n-iz)+bx(n+ix-iz)-r*(bz(n) &
		 +c*((1.-s)*(ex(n+iy)-ex(n))+(1.+s)*.25*(ez(n+iy) &
		 -ez(n)+ez(n+ix+iy)-ez(n+ix)+ez(n+iy-iz)-ez(n-iz) &
		 +ez(n+ix+iy-iz)-ez(n+ix-iz))))
	enddo

!	do n=m+iy*(my-1),m+iy*(my-1)+iz*(mz-2),iz
	 n=m+iy*(my-1)
		bx(n)=(1.-c*r)*(bx(n)+bx(n+ix))+bx(n-iy)+bx(n+ix-iy) &
		 -r*(by(n)-c*((1.-s)*(ex(n+iz)-ex(n)) &
		  +(1.+s)*.25*(ey(n+iz)-ey(n)+ey(n+ix+iz) &
		 -ey(n+ix)+ey(n+iz-iy)-ey(n-iy)+ey(n+ix+iz-iy) &
		 -ey(n+ix-iy))))
!	enddo

	p=(1.+c)*2./(1.+2.*c*(1.+c*s))
	q=c*s*2./(1.+2.*c*(1.+c*s))
	
	do n=m+iy*(my-1)+iz*(mz-1),m+ix*(mx-2)+iy*(my-1)+iz*(mz-1),ix
		temp=bz(n)-.5*c*(1.-s)*(ey(n+ix)-ey(n)+ey(n+ix-iy)-ey(n-iy))
		bz(n)=bz(n-iy)-bz(n)+p*temp+q*by(n)
		by(n)=by(n-iz)-by(n)+p*by(n)+q*temp
	enddo

	endif
#endif



end subroutine preledge



!-------------------------------------------------------------------------------
! 						subroutine postedge()					 
!																		
! 
!
!-------------------------------------------------------------------------------

subroutine postedge(bx,by,bz,ex,ey,ez,ix,iy,iz,mx,my,mz,m,c)

!	implicit none

	! dummy variables

	real :: bx(1), by(1), bz(1), ex(1), ey(1), ez(1)
!	real(sprec), dimension(:), intent(inout) :: bx, by, bz, ex, ey, ez
	integer, intent(in) :: ix,iy,iz,mx,my,mz,m
	real(sprec), intent(in) :: c

	! local variables
	
	real(sprec) :: s, p, r, t, q, temp
	integer :: n

	s=.4142136
	p=(1.+c)*2./(1.+2.*c*(1.+c*s))
	q=c*s*2./(1.+2.*c*(1.+c*s))

#ifndef twoD
	
	do n=m+iy*(my-1)+iz*(mz-1),m+ix*(mx-2)+iy*(my-1)+iz*(mz-1),ix
		temp=by(n-iz)-.5*c*(1.-s)*(ez(n+ix)-ez(n)+ez(n+ix-iz)-ez(n-iz))
		bz(n)=bz(n-iy)+bz(n)-q*temp-p*bz(n-iy)
		by(n)=by(n-iz)+by(n)-q*bz(n-iy)-p*temp
	enddo
	
	t=c/(2.*c+1.+s)
	
	do n=m+iy*(my-1)+iz*(mz-1)+ix,m+iy*(my-1)+iz*(mz-1)+ix*(mx-2),ix
		bx(n)=bx(n)-(1.-4.*t)*bx(n-iy-iz)-(1.-2.*t)*(bx(n-iy)+bx(n-iz)) &
		+s*t*(by(n)-by(n-ix)+by(n-iz)-by(n-ix-iz) &
		+bz(n)-bz(n-ix)+bz(n-iy)-bz(n-ix-iy))
	enddo

	r=4./(2.*c+2.+s)
	
	do n=m+iz*(mz-1),m+iz*(mz-1)+iy*(my-2),iy
		bx(n)=bx(n)-bx(n+ix)-(1.-c*r)*(bx(n-iz)+bx(n+ix-iz))+r*bz(n)
	enddo
	
	do n=m+iy*(my-1),m+iy*(my-1)+iz*(mz-2),iz
		bx(n)=bx(n)-bx(n+ix)-(1.-c*r)*(bx(n-iy)+bx(n+ix-iy))+r*by(n)
	enddo

#else
	if(ix .eq. 0) then 
!	do n=m+iy*(my-1)+iz*(mz-1),m+ix*(mx-2)+iy*(my-1)+iz*(mz-1),ix
	n=m+iy*(my-1)+iz*(mz-1)
		temp=by(n-iz)-.5*c*(1.-s)*(ez(n+ix)-ez(n)+ez(n+ix-iz)-ez(n-iz))
		bz(n)=bz(n-iy)+bz(n)-q*temp-p*bz(n-iy)
		by(n)=by(n-iz)+by(n)-q*bz(n-iy)-p*temp
!	enddo
	
	t=c/(2.*c+1.+s)
	
!	do n=m+iy*(my-1)+iz*(mz-1)+ix,m+iy*(my-1)+iz*(mz-1)+ix*(mx-2),ix
	 n=m+iy*(my-1)+iz*(mz-1)+ix
		bx(n)=bx(n)-(1.-4.*t)*bx(n-iy-iz)-(1.-2.*t)*(bx(n-iy)+bx(n-iz)) &
		+s*t*(by(n)-by(n-ix)+by(n-iz)-by(n-ix-iz) &
		+bz(n)-bz(n-ix)+bz(n-iy)-bz(n-ix-iy))
!	enddo

	r=4./(2.*c+2.+s)
	
	do n=m+iz*(mz-1),m+iz*(mz-1)+iy*(my-2),iy
		bx(n)=bx(n)-bx(n+ix)-(1.-c*r)*(bx(n-iz)+bx(n+ix-iz))+r*bz(n)
	enddo
	
	do n=m+iy*(my-1),m+iy*(my-1)+iz*(mz-2),iz
		bx(n)=bx(n)-bx(n+ix)-(1.-c*r)*(bx(n-iy)+bx(n+ix-iy))+r*by(n)
	enddo

	endif

	if(iy .eq. 0) then 
	do n=m+iy*(my-1)+iz*(mz-1),m+ix*(mx-2)+iy*(my-1)+iz*(mz-1),ix
		temp=by(n-iz)-.5*c*(1.-s)*(ez(n+ix)-ez(n)+ez(n+ix-iz)-ez(n-iz))
		bz(n)=bz(n-iy)+bz(n)-q*temp-p*bz(n-iy)
		by(n)=by(n-iz)+by(n)-q*bz(n-iy)-p*temp
	enddo
	
	t=c/(2.*c+1.+s)
	
	do n=m+iy*(my-1)+iz*(mz-1)+ix,m+iy*(my-1)+iz*(mz-1)+ix*(mx-2),ix
		bx(n)=bx(n)-(1.-4.*t)*bx(n-iy-iz)-(1.-2.*t)*(bx(n-iy)+bx(n-iz)) &
		+s*t*(by(n)-by(n-ix)+by(n-iz)-by(n-ix-iz) &
		+bz(n)-bz(n-ix)+bz(n-iy)-bz(n-ix-iy))
	enddo

	r=4./(2.*c+2.+s)
	
!	do n=m+iz*(mz-1),m+iz*(mz-1)+iy*(my-2),iy
	 n=m+iz*(mz-1)
		bx(n)=bx(n)-bx(n+ix)-(1.-c*r)*(bx(n-iz)+bx(n+ix-iz))+r*bz(n)
!	enddo
	
	do n=m+iy*(my-1),m+iy*(my-1)+iz*(mz-2),iz
		bx(n)=bx(n)-bx(n+ix)-(1.-c*r)*(bx(n-iy)+bx(n+ix-iy))+r*by(n)
	enddo

	endif

	if(iz .eq. 0) then 
	do n=m+iy*(my-1)+iz*(mz-1),m+ix*(mx-2)+iy*(my-1)+iz*(mz-1),ix
		temp=by(n-iz)-.5*c*(1.-s)*(ez(n+ix)-ez(n)+ez(n+ix-iz)-ez(n-iz))
		bz(n)=bz(n-iy)+bz(n)-q*temp-p*bz(n-iy)
		by(n)=by(n-iz)+by(n)-q*bz(n-iy)-p*temp
	enddo
	
	t=c/(2.*c+1.+s)
	
	do n=m+iy*(my-1)+iz*(mz-1)+ix,m+iy*(my-1)+iz*(mz-1)+ix*(mx-2),ix
		bx(n)=bx(n)-(1.-4.*t)*bx(n-iy-iz)-(1.-2.*t)*(bx(n-iy)+bx(n-iz)) &
		+s*t*(by(n)-by(n-ix)+by(n-iz)-by(n-ix-iz) &
		+bz(n)-bz(n-ix)+bz(n-iy)-bz(n-ix-iy))
	enddo

	r=4./(2.*c+2.+s)
	
	do n=m+iz*(mz-1),m+iz*(mz-1)+iy*(my-2),iy
		bx(n)=bx(n)-bx(n+ix)-(1.-c*r)*(bx(n-iz)+bx(n+ix-iz))+r*bz(n)
	enddo
	
!	do n=m+iy*(my-1),m+iy*(my-1)+iz*(mz-2),iz
	 n=m+iy*(my-1)
		bx(n)=bx(n)-bx(n+ix)-(1.-c*r)*(bx(n-iy)+bx(n+ix-iy))+r*by(n)
!	enddo

	endif

#endif

end subroutine postedge



#ifdef filter1  !use direct averaging filter if filter1 is undefined; otherwise -- repeated passes through grid
!-------------------------------------------------------------------------------
! 				subroutine apply_filter1_opt()		
!										
! 
!
!-------------------------------------------------------------------------------

subroutine apply_filter1_opt()

	implicit none   
	
	! local variables
	
	real wtm1,wt,wtp1,winv
	integer :: istr, ifin, i, j, k, n

! apply digital filter to currents to smooth them 
! 1 2 1 in every direction
! done as convolution of 3 one-dimensional sweeps
! filtering is repeated ntimes	
	
	n=1
	
	istr=2  +1
	ifin=mx-1 -2


!	if(wall) then 
!	   istr=2		!right going injection
!	   ifin=mx-1
!	   ifin=min(mx-1,int(xinject2)+10)
!	endif

! ! 	if(rank .eq. 0) print *,"xinject2",xinject2

		wtm1=1.
		wt=2.
		wtp1=1.
		winv=1./4.**3
  
  winv=1./4. !if we do winv multiplication right away
#ifdef twoD
 winv=1./16.
 winv=1./4.
#endif

   wtm1=wtm1*winv
   wt=wt*winv
   wtp1=wtp1*winv

	do while(n .le. ntimes) 
		temp=0.
		
		!filter in x direction

		if(periodicx.eq.1) then 
!			call copy_layrx1(curx,cury,curz,mx,my,mz,2,mx-3,mx-2,3)
			call copy_layrx1_opt(curx,cury,curz,mx,my,mz,2,mx-3,mx-2,3)
		else
!			call copy_layrx2(curx,cury,curz,mx,my,mz,2,mx-3,mx-2,3)
			call copy_layrx2_opt(curx,cury,curz,mx,my,mz,2,mx-3,mx-2,3)
		endif


#ifndef twoD
		do k=3,mz-3 
#else
		do k=1,1
#endif
			do j=3,my-3 
				do i=istr,ifin 
					temp(i,j,k)=wtm1*curx(i-1,j,k)+wt*curx(i,j,k)+wtp1*curx(i+1,j,k)
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
					temp(i,j,k)=wtm1*cury(i-1,j,k)+wt*cury(i,j,k)+wtp1*cury(i+1,j,k)
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
					temp(i,j,k)=wtm1*curz(i-1,j,k)+wt*curz(i,j,k)+wtp1*curz(i+1,j,k)
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
		
		!filter in y direction
		
		if(periodicy.eq.1) then 
!			call copy_layry1(curx,cury,curz,mx,my,mz,2,my-3,my-2,3)
			call copy_layry1_opt(curx,cury,curz,mx,my,mz,2,my-3,my-2,3)
		else
!			call copy_layry2(curx,cury,curz,mx,my,mz,2,my-3,my-2,3)
			call copy_layry2_opt(curx,cury,curz,mx,my,mz,2,my-3,my-2,3)
		endif
		
#ifndef twoD
		do k=3,mz-3 
#else
		do k=1,1
#endif
			do j=3,my-3
				do i=istr,ifin 
					temp(i,j,k)=wtm1*curx(i,j-1,k)+wt*curx(i,j,k)+wtp1*curx(i,j+1,k)
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
					temp(i,j,k)=wtm1*cury(i,j-1,k)+wt*cury(i,j,k)+wtp1*cury(i,j+1,k)
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
					temp(i,j,k)=wtm1*curz(i,j-1,k)+wt*curz(i,j,k)+wtp1*curz(i,j+1,k)
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
		
		!filter in z direction
#ifndef twoD
			if(periodicz.eq.1) then
				!call copy_layrz1(curx,cury,curz,mx,my,mz,2,mz-3,mz-2,3)
				 call copy_layrz1_opt(curx,cury,curz,mx,my,mz,2,my-3,my-2,3)		
			else
				!transfer z information among processors, not periodic between 0,n-1
				!call copy_layrz2(curx,cury,curz,mx,my,mz,2,mz-3,mz-2,3)
				call copy_layrz2_opt(curx,cury,curz,mx,my,mz,2,mz-3,mz-2,3)
			endif
			
			do k=3,mz-3
				do j=3,my-3 
					do i=istr,ifin
						temp(i,j,k)=wtm1*curx(i,j,k-1)+wt*curx(i,j,k)+wtp1*curx(i,j,k+1)
					enddo
				enddo
			enddo
			
			do k=3,mz-3 
				do j=3,my-3 
					do i=istr,ifin 
						curx(i,j,k)=temp(i,j,k) !*winv
					enddo
				enddo
			enddo
			
			do k=3,mz-3
				do j=3,my-3 
					do i=istr,ifin 
						temp(i,j,k)=wtm1*cury(i,j,k-1)+wt*cury(i,j,k)+wtp1*cury(i,j,k+1)
					enddo
				enddo
			enddo
			
			do k=3,mz-3 
				do j=3,my-3 
					do i=istr,ifin 
						cury(i,j,k)=temp(i,j,k) !*winv
					enddo
				enddo
			enddo
			
			do k=3,mz-3
				do j=3,my-3 
					do i=istr,ifin
						temp(i,j,k)=wtm1*curz(i,j,k-1)+wt*curz(i,j,k)+wtp1*curz(i,j,k+1)
					enddo
				enddo
			enddo
			
			do k=3,mz-3 
				do j=3,my-3 
					do i=istr,ifin 
						curz(i,j,k)=temp(i,j,k) !*winv
					enddo
				enddo
			enddo
			
#else
!#define dontskipme
#ifdef dontskipme  !don't need this if winv is absorbed into the weights
			!twod: skip z averaging, use new winv
			do k=1,1
				do j=3,my-3 
					do i=istr,ifin 
						curx(i,j,k)=curx(i,j,k) !*winv
					enddo
				enddo
			enddo

			do k=1,1
				do j=3,my-3 
					do i=istr,ifin 
						cury(i,j,k)=cury(i,j,k) !*winv
					enddo
				enddo
			enddo

			do k=1,1
				do j=3,my-3 
					do i=istr,ifin 
						curz(i,j,k)=curz(i,j,k) !*winv
					enddo
				enddo
			enddo
#endif !dontskipme			
#endif
		n=n+1

!if(rank .eq. 10) then 
!			do k=1,1
!				do j=3,my-3 
!					do i=100,200 
!						print *, n, "curx",i,j,k, curx(i,j,k)
!					enddo
!				enddo
!			enddo
!endif
	enddo

end subroutine apply_filter1_opt
#else !ifdef filter 1
#include "filter.F90" !real apply_filter1_opt
#endif


#ifdef filter2
#include "optimized_filters.F90"
#endif !ifdef filter2

#ifdef twoD
end module m_fieldboundaries
#else
end module m_fieldboundaries_3d
#endif
