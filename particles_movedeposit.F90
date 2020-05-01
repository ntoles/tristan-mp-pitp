

!
! Particle module
!
! Includes the particle data structures and the routines for depositing
! and moving the particles
!
!

#ifdef twoD 

module m_particles_movedeposit

	use m_globaldata
	use m_system
	use m_aux
	use m_communications
	use m_fields
	use m_inputparser
	use m_fparser
	use m_particles
        use m_user
#else

module m_particles_movedeposit_3d

	use m_globaldata_3d
	use m_system_3d
	use m_aux_3d
	use m_communications_3d
	use m_fields_3d
	use m_inputparser_3d
	use m_fparser_3d
        use m_particles_3d
        use m_user_3d
#endif


	
	implicit none
        
!-------------------------------------------------------------------------------
!	INTERFACE DECLARATIONS
!-------------------------------------------------------------------------------
	
!-------------------------------------------------------------------------------
!	PUBLIC MODIFIERS
!-------------------------------------------------------------------------------

	! public functions
	
	public :: deposit_particles, move_particles

        contains
!-------------------------------------------------------------------------------
! 						subroutine move_particles()					 
!																		
! Pushes the particle's velocity and positions in time
!
!-------------------------------------------------------------------------------

subroutine move_particles()

	implicit none
	
	integer ierr

		call mover(1,ions,qmi)			! ions				
		call mover(1+maxhlf,lecs+maxhlf,qme)	! electrons
	
	if(debug) print *,rank,": done mover", "step=", lap

end subroutine move_particles



!-------------------------------------------------------------------------------
! 						subroutine mover()					 
! 
!
!-------------------------------------------------------------------------------

subroutine mover(n1,n2,qm)

	implicit none

	! dummy variables
	
	integer :: n2, n1
	real(sprec) :: qm

	! local variables

	integer :: npr, l, n, i, j, k
	real dx, dy, dz, f, g, ex0, ey0, ez0, bx0, by0, bz0
	real u0,v0,w0,u1,v1,w1, cinv, g1, corrqm
	real qm0

#ifdef vay
	real ustar,sig,tx,ty,tz,vx0,vy0,vz0
#endif
	real bx_ext, by_ext, bz_ext, ex_ext, ey_ext, ez_ext

	cinv=1./c	
	qm0=qm

	  do n=n1,n2

	      npr=n
	      
	!##########################################################################
	!### 	tri-linear field interpolation (original version)
	!##########################################################################
			
		i=aint(p(npr)%x)
		dx=p(npr)%x-i
		j=p(npr)%y
		dy=p(npr)%y-j
		k=p(npr)%z
		dz=p(npr)%z-k

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
		
		ex0=(f+dz*(g-f))*(.25*qm)
		
		!   ------------
		f=ey(l,1,1)+ey(l-iy,1,1)+dy*(ey(l+iy,1,1)-ey(l-iy,1,1))
		f=f+dz*(ey(l+iz,1,1)+ey(l-iy+iz,1,1)+dy*(ey(l+iy+iz,1,1)-ey(l-iy &
		+iz,1,1))-f)
		g=ey(l+ix,1,1)+ey(l-iy+ix,1,1)+dy*(ey(l+iy+ix,1,1)-ey(l-iy+ix,1 &
		,1))
		g=g+dz* &
		(ey(l+iz+ix,1,1)+ey(l-iy+iz+ix,1,1)+dy*(ey(l+iy+iz+ix,1,1) &
		-ey(l-iy+iz+ix,1,1))-g)
		
		ey0=(f+dx*(g-f))*(.25*qm)
		
		!   ------------
		f=ez(l,1,1)+ez(l-iz,1,1)+dz*(ez(l+iz,1,1)-ez(l-iz,1,1))
		f=f+dx*(ez(l+ix,1,1)+ez(l-iz+ix,1,1)+dz*(ez(l+iz+ix,1,1)-ez(l-iz &
		+ix,1,1))-f)
		g=ez(l+iy,1,1)+ez(l-iz+iy,1,1)+dz*(ez(l+iz+iy,1,1)-ez(l-iz+iy,1 &
		,1))
		g=g+dx* &
		(ez(l+ix+iy,1,1)+ez(l-iz+ix+iy,1,1)+dz*(ez(l+iz+ix+iy,1,1) &
		-ez(l-iz+ix+iy,1,1))-g)
		
		ez0=(f+dy*(g-f))*(.25*qm)
		
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
		
		bx0=(f+dx*(g-f))*(.125*qm*cinv)
		
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
		
		by0=(f+dy*(g-f))*(.125*qm*cinv)
		
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
		
		bz0=(f+dz*(g-f))*(.125*qm*cinv)
		
			
                if(external_fields) then
                call get_external_fields(real(p(npr)%x,sprec),real(p(npr)%y,sprec),&
                p(npr)%z,ex_ext,ey_ext,ez_ext,bx_ext,by_ext,bz_ext,qm)
                   
                bx0=bx0+bx_ext*0.5*qm*cinv
                by0=by0+by_ext*0.5*qm*cinv
                bz0=bz0+bz_ext*0.5*qm*cinv
                ex0=ex0+ex_ext*0.5*qm
                ey0=ey0+ey_ext*0.5*qm
                ez0=ez0+ez_ext*0.5*qm

		endif

		!   First half electric acceleration, with relativity's gamma:		
#ifdef vay
			!Use Vay 2008 particle mover
			g=1./sqrt(1.+p(npr)%u**2+p(npr)%v**2+p(npr)%w**2) !reciprocal of the Lorentz factor
			vx0=c*p(npr)%u*g !3-velocity of the particle
			vy0=c*p(npr)%v*g
			vz0=c*p(npr)%w*g
			
			u1=c*p(npr)%u+2.*ex0+vy0*bz0-vz0*by0 !u-prime, taking into account 
                                                             !that cinv is already incorporated within B
			v1=c*p(npr)%v+2.*ey0+vz0*bx0-vx0*bz0
			w1=c*p(npr)%w+2.*ez0+vx0*by0-vy0*bx0

			!Lorentz factor for uprime			

			ustar=cinv*(u1*bx0+v1*by0+w1*bz0)
			sig=cinv*cinv*(c**2+u1**2+v1**2+w1**2)-(bx0**2+by0**2+bz0**2)
			g=1./sqrt(0.5*(sig+sqrt(sig**2+4.*(bx0**2+by0**2+bz0**2+ustar**2))))
			tx=bx0*g
			ty=by0*g
			tz=bz0*g
			f=1./(1.+tx**2+ty**2+tz**2)
			
			u0=f*(u1+(u1*tx+v1*ty+w1*tz)*tx+v1*tz-w1*ty)
			v0=f*(v1+(u1*tx+v1*ty+w1*tz)*ty+w1*tx-u1*tz)
			w0=f*(w1+(u1*tx+v1*ty+w1*tz)*tz+u1*ty-v1*tx)				
#else
			!Use Boris algorithm
			!   First half electric acceleration, with Lorentz gamma:
			u0=c*p(npr)%u+ex0
			v0=c*p(npr)%v+ey0
			w0=c*p(npr)%w+ez0
			!   First half magnetic rotation, with Lorentz gamma:
			g=c/sqrt(c**2+u0**2+v0**2+w0**2)
			bx0=g*bx0
			by0=g*by0
			bz0=g*bz0
			
			f=2./(1.+bx0*bx0+by0*by0+bz0*bz0)
			u1=(u0+v0*bz0-w0*by0)*f
			v1=(v0+w0*bx0-u0*bz0)*f
			w1=(w0+u0*by0-v0*bx0)*f
			!   Second half mag. rot'n &   el. acc'n:
			u0=u0+v1*bz0-w1*by0+ex0 
			v0=v0+w1*bx0-u1*bz0+ey0 
			w0=w0+u1*by0-v1*bx0+ez0 
#endif
			!   Get normalized 4-velocity:
			p(npr)%u=u0*cinv
			p(npr)%v=v0*cinv
			p(npr)%w=w0*cinv
		
		!   Position advance:
			g=c/sqrt(c**2+u0**2+v0**2+w0**2)

			p(npr)%x=p(npr)%x + p(npr)%u*g*c
			p(npr)%y=p(npr)%y + p(npr)%v*g*c 
			p(npr)%z=p(npr)%z + p(npr)%w*g*c

	enddo

end subroutine mover

!-------------------------------------------------------------------------------
! 						subroutine deposit_particles()					 
!																		
! Deposits the particles on the grid
!
!-------------------------------------------------------------------------------
subroutine deposit_particles()

	implicit none
	
	! local variables
	
	integer ::ind1,n1,n
	logical in
	real(sprec) :: x0, y0, perx, pery
	real z0,perz,invgam
	
	integer ierr, ions0, lecs0
	
	integer i1, j1, k1

	real gc_x, gc_y, gc_z !hack, vars for printing
	
	real xt,yt,zt !hack variables for printing 	
        real midx, midy, midz, maxx, maxy, maxz, minx, miny, minz

	real(dprec) :: x1,x2, y1, y2
	real(sprec) :: z1,z2

        real x1sp, x2sp, y1sp, y2sp
	
	! relay point (xr,yr,zr)
!	real(dprec) :: xr, yr
        real :: xr, yr
	real zr
	! charge fluxes in 1st order scheme
	real Fx1, Fx2, Fy1, Fy2, Fz1, Fz2
	! shape functions
	real Wx1, Wx2, Wy1, Wy2, Wz1, Wz2, onemWx1, onemWx2, onemWy1, onemWy2, onemWz1, onemWz2
	integer :: i2, j2, k2  
       
        integer :: i1p1, i2p1, j1p1, j2p1
        
#ifndef twoD
        integer :: k1p1, k2p1
#endif
 
 integer :: nlev, levels, m

 levels=8

	LenIonOutUp=0
	LenIonOutDwn=0
	LenLecOutUp=0
	LenLecOutDwn=0
	LenIonInBlw=0
	LenIonInAbv=0
	LenLecInBlw=0
	LenLecInAbv=0
	
	LenIonOutLft=0
	LenIonOutRgt=0
	LenLecOutLft=0
	LenLecOutRgt=0
	LenIonInLft=0
	LenIonInRgt=0
	LenLecInLft=0
	LenLecInRgt=0
	
	LenIonOutMinus=0
	LenIonOutPlus=0
	LenLecOutMinus=0
	LenLecOutPlus=0
	LenIonInMinus=0
	LenIonInPlus=0
	LenLecInMinus=0
	LenLecInPlus=0	
		
	nionout=0
	nlecout=0	
	
	!particles are in local coordinates
	n=1

        maxx=mx-2.
        minx=3.
        maxy=my-2.
        miny=3.

        midx=.5*(maxx-minx)
        midy=.5*(maxy-miny)

#ifndef twoD
        maxz=mz-2.
        minz=3.
#else
        minz=3.
        maxz= 6.-2.
#endif
        midz = .5*(maxz-minz)

	if (ions.gt.0) then

#define cleancur

#ifdef cleancur	
    do n=1,ions
	n1=n
#endif
             invgam=1./sqrt(1+p(n1)%u**2+p(n1)%v**2+p(n1)%w**2)
		
		x0=p(n1)%x-p(n1)%u*invgam*c
		y0=p(n1)%y-p(n1)%v*invgam*c
		z0=p(n1)%z-p(n1)%w*invgam*c
		
		q=p(n1)%ch*qi ! real(splitratio)**(1.-real(p(n1)%splitlev))*qi 

!		call zigzag(p(n1)%x,p(n1)%y,p(n1)%z,x0,y0,z0,in)
  
	! zigzag method (1st order version) (Ref: Umeda et al. 2003)

        x1sp=x0
        x2sp=p(n1)%x
        y1sp=y0
        y2sp=p(n1)%y

        z1=z0
        z2=p(n1)%z

        i1=int(x1sp)
        i2=int(x2sp)
        j1=int(y1sp)
        j2=int(y2sp)
        k1=int(z1)
        k2=int(z2)


	xr=min(real(min(i1,i2)+1),max(real(max(i1,i2)),.5*(x1sp+x2sp)))
	yr=min(real(min(j1,j2)+1),max(real(max(j1,j2)),.5*(y1sp+y2sp)))
	zr=min(real(min(k1,k2)+1),max(real(max(k1,k2)),.5*(z1+z2)))
	
#ifdef twoD
		k1=1
		k2=1
#endif

	!-q to include -j in the Ampere's equation, to be consistent
	
	Fx1=-q*(xr-x1sp)
	Fy1=-q*(yr-y1sp)
	Fz1=-q*(zr-z1)
	
	Wx1=.5*(x1sp+xr)-i1
	Wy1=.5*(y1sp+yr)-j1

#ifndef twoD
        Wz1=.5*(z1+zr)-k1
#endif

	Wx2=.5*(x2sp+xr)-i2
	Wy2=.5*(y2sp+yr)-j2

#ifndef twoD
		Wz2=.5*(z2+zr)-k2
#endif
	
	Fx2=-q*(x2sp-xr)
	Fy2=-q*(y2sp-yr)
	Fz2=-q*(z2-zr)
	
#ifdef twoD
		Wz1=0
		Wz2=0
#endif

  onemWx1=1.-Wx1
  onemWx2=1.-Wx2
  onemWy1=1.-Wy1
  onemWy2=1.-Wy2
  onemWz1=1.-Wz1
  onemWz2=1.-Wz2

  i1p1=i1+1
  i2p1=i2+1
  j1p1=j1+1
  j2p1=j2+1
#ifndef twoD
  k1p1=k1+1
  k2p1=k2+1
#endif  

#ifdef cleancur

  curx(i1,j1,k1)=curx(i1,j1,k1)+Fx1 * onemWy1 * onemWz1
  curx(i1,j1p1,k1)=curx(i1,j1p1,k1)+Fx1 * Wy1 * onemWz1 
#ifndef twoD
  curx(i1,j1,  k1p1)= curx(i1,j1,  k1p1)+Fx1 * onemWy1 * Wz1
  curx(i1,j1p1,k1p1)= curx(i1,j1p1,k1p1)+Fx1 *  Wy1    * Wz1
#endif

  curx(i2,j2,k2)=curx(i2,j2,k2)+Fx2 * onemWy2 * onemWz2 
  curx(i2,j2p1,k2)=curx(i2,j2p1,k2)+Fx2 * Wy2 * onemWz2 
#ifndef twoD
  curx(i2,j2,  k2p1)= curx(i2,j2,  k2p1)+Fx2 * onemWy2* Wz2
  curx(i2,j2p1,k2p1)= curx(i2,j2p1,k2p1)+Fx2 *  Wy2    * Wz2
#endif

  cury(i1,j1,k1)=cury(i1,j1,k1)+Fy1 * onemWx1 * onemWz1 
  cury(i1p1,j1,k1)=cury(i1p1,j1,k1)+Fy1 * Wx1 * onemWz1  
#ifndef twoD
  cury(i1  ,j1,k1p1)= cury(i1  ,j1,k1p1)+Fy1 * onemWx1 * Wz1
  cury(i1p1,j1,k1p1)= cury(i1p1,j1,k1p1)+Fy1 *  Wx1    * Wz1
#endif
  cury(i2,j2,k2)=cury(i2,j2,k2)+Fy2 * onemWx2 * onemWz2 
  cury(i2p1,j2,k2)=cury(i2p1,j2,k2)+Fy2 *  Wx2 * onemWz2
#ifndef twoD
  cury(i2  ,j2,k2p1)= cury(i2  ,j2,k2p1)+Fy2 * onemWx2 * Wz2
  cury(i2p1,j2,k2p1)= cury(i2p1,j2,k2p1)+Fy2 *  Wx2    * Wz2
#endif

  curz(i1,j1,k1)=curz(i1,j1,k1)+ Fz1 * onemWx1 * onemWy1
  curz(i1p1,j1,k1)=curz(i1p1,j1,k1)+Fz1 * Wx1 * onemWy1
  curz(i1,j1p1,k1)=curz(i1,j1p1,k1)+Fz1 * onemWx1 * Wy1
  curz(i1p1,j1p1,k1)=curz(i1p1,j1p1,k1)+Fz1 * Wx1 * Wy1 

  curz(i2,j2,k2)=curz(i2,j2,k2)+ Fz2 * onemWx2 * onemWy2 
  curz(i2p1,j2,k2)=curz(i2p1,j2,k2)+ Fz2 * Wx2 * onemWy2 
  curz(i2,j2p1,k2)=curz(i2,j2p1,k2)+ Fz2 * onemWx2 * Wy2
  curz(i2p1,j2p1,k2)=curz(i2p1,j2p1,k2)+ Fz2 * Wx2 * Wy2 

#endif !cleancur
        enddo !m,n=1,ions

#define dontSKIP
#ifdef dontSKIP
        pind(1:ions)=0

        do n=1,ions
        in=.true.
	n1=n	
  perx=0.
  pery=0.
  perz=0.

  if(p(n1)%x .lt. minx .or. p(n1)%x .gt. maxx) then 
 		perx=sign(midx,p(n1)%x-minx)+sign(midx,p(n1)%x-maxx)    
  endif
  if(p(n1)%y .lt. miny .or. p(n1)%y .gt. maxy) then 
 		pery=sign(midy,p(n1)%y-miny)+sign(midy,p(n1)%y-maxy)    
  endif
  if(p(n1)%z .lt. minz .or. p(n1)%z .gt. maxz) then 
		perz=sign(midz,p(n1)%z-minz)+sign(midz,p(n1)%z-maxz)
  endif

		if(periodicx.eq.0) then      
			in=(p(n1)%x+mxcum .gt. x1in) .and. (p(n1)%x+mxcum .lt. x2in)
		endif
		
		if(periodicy.eq.0 .and. in) then 
			in=(p(n1)%y+mycum .gt. y1in) .and. (p(n1)%y+mycum .lt. y2in)
		endif
#ifndef twoD						
		if(periodicz.eq.0 .and. in) then      
			in=(p(n1)%z+mzcum .gt. z1in) .and. (p(n1)%z+mzcum .lt. z2in)
		endif		
#endif

		! to discard a particle
		if(.not. in) then
			perx=0
			pery=0
			perz=0
		endif		
		
		! to send to another processor
		if(perx .ne. 0 .and. in .and. sizex .ne. 1) then
			in=.false.	
			if(sizey .ne. 1) pery=0
			if(sizez .ne. 1) perz=0
		endif
				
		! assume non-uniform mx
		if(perx .lt. 0 .and. sizex .ne. 1) then		
			! minus rank (-x direction)		
			i1=(rank/sizex)*sizex + modulo(rank-1,sizex) 
			perx=-(mxl(i1+1)-5.)				
		endif
		
		p(n1)%x=p(n1)%x-perx		
				
		! to send to another processor
		if(pery .ne. 0 .and. in .and. sizey .ne. 1) then
			in=.false.
			if(sizex .ne. 1) perx=0
			if(sizez .ne. 1) perz=0
		endif
		
		! assume non-uniform my
		if(pery .lt. 0 .and. sizey .ne. 1) then		
			! left rank	(-y direction)	
!			j1=modulo((rank/sizex - 1),sizey)*sizex + modulo(rank,sizex)
			j1=modulo(rank/sizex - 1,sizey)*sizex+rank/(sizex*sizey)*(sizex*sizey) &
				+ modulo(rank,sizex) 
			pery=-(myl(j1+1)-5.)				
		endif
		
		p(n1)%y=p(n1)%y-pery
		
#ifndef twoD
		! to send to another processor
		if(perz .ne. 0 .and. in ) then
			in=.false.
			if(sizex .ne. 1) perx=0
			if(sizey .ne. 1) pery=0
		endif	
		! assume non-uniform mz
		if(perz .lt. 0) then		
			! down rank	(-z direction)	
			k1=modulo(rank/(sizex*sizey) - 1,sizez)*(sizex*sizey) + & 
			      modulo(rank,sizex*sizey) 
			perz=-(mzl(k1+1)-5.)				
		endif
#endif
		
		p(n1)%z=p(n1)%z-perz		
		
		
		if(in) go to 58
	
		nionout=nionout+1
			
		!check if the paricle is going to another proc
#ifndef	twoD	
		if(perz .lt. 0) then
			!this particle is going to the processor below
		      LenIonOutDwn=LenIonOutDwn+1
!		      poutdwn(LenIonOutDwn)=p(n1)           
		      call copyprt(p(n1), poutdwn(LenIonOutDwn))
		endif
	
		if(perz .gt. 0) then
			!this particle is going to the processor above
			LenIonOutUp=LenIonOutUp+1
	!		poutup(LenIonOutUp)=p(n1)
			call copyprt(p(n1),poutup(LenIonOutUp))

		endif
#endif
		
		if(sizey .ne. 1) then
		if(pery .lt. 0) then
			!this particle is going to the processor on the left
			LenIonOutLft=LenIonOutLft+1
!			poutlft(LenIonOutLft)=p(n1)
			call copyprt(p(n1),poutlft(LenIonOutLft))
		endif

		if(pery .gt. 0) then
			!this particle is going to the processor on the right
			LenIonOutRgt=LenIonOutRgt+1
			call copyprt(p(n1),poutrgt(LenIonOutRgt))
		endif
		endif !if(sizey .ne. 1)
		
		if(sizex .ne. 1) then
		if(perx .lt. 0) then
			!this particle is going to the processor on the left
			LenIonOutMinus=LenIonOutMinus+1
			call copyprt(p(n1),poutminus(LenIonOutMinus))
		endif

		if(perx .gt. 0) then
			!this particle is going to the processor on the right
			LenIonOutPlus=LenIonOutPlus+1
			call copyprt(p(n1),poutplus(LenIonOutPlus))
		endif
		endif !if(sizex .ne. 1) then
  
  pind(n1)=1
!  print *, rank, "pind=1", perx, pery, perz, n1
58 continue	

  enddo !n=1,ions
  
!now clean particles
  ions0=ions
  do n=1,ions0
     if(pind(n) .ne. 0) then
        do while (pind(n) .ne. 0)
!           print *, "discarding", pind(n), n, p(n)%x, p(n)%y, p(n)%z, ions, ions0, lap, rank
           call copyprt(p(ions),p(n))
           pind(n)=pind(ions)
           pind(ions)=0
           ions=ions-1
        enddo
     endif
  enddo
#endif !dontSKIP

	endif !if ions .gt. 0
		
	if(debug) print *, rank,":","deposited ions","step=",lap 
	
	if(lecs.gt.0) then
!		53      continue

        do n=1,lecs
        in=.true.

		n1=n+maxhlf !indp(n)

		invgam=1./sqrt(1+p(n1)%u**2+p(n1)%v**2+p(n1)%w**2)
		
		x0=p(n1)%x-p(n1)%u*invgam*c
		y0=p(n1)%y-p(n1)%v*invgam*c
		z0=p(n1)%z-p(n1)%w*invgam*c
		
		q=p(n1)%ch*qe 
		
!		call zigzag(p(n1)%x,p(n1)%y,p(n1)%z,x0,y0,z0,in)

!  goto 107
        x1sp=x0
        x2sp=p(n1)%x
        y1sp=y0
        y2sp=p(n1)%y

        z1=z0
        z2=p(n1)%z

        i1=int(x1sp)
        i2=int(x2sp)
        j1=int(y1sp)
        j2=int(y2sp)
        k1=int(z1)
        k2=int(z2)


	xr=min(real(min(i1,i2)+1),max(real(max(i1,i2)),.5*(x1sp+x2sp)))
	yr=min(real(min(j1,j2)+1),max(real(max(j1,j2)),.5*(y1sp+y2sp)))
	zr=min(real(min(k1,k2)+1),max(real(max(k1,k2)),.5*(z1+z2)))
	
#ifdef twoD
		k1=1
		k2=1
#endif

	!-q to include -j in the Ampere's equation, to be consistent
	
	Fx1=-q*(xr-x1sp)
	Fy1=-q*(yr-y1sp)
	Fz1=-q*(zr-z1)
	
	Wx1=.5*(x1sp+xr)-i1
	Wy1=.5*(y1sp+yr)-j1

#ifndef twoD
        Wz1=.5*(z1+zr)-k1
#endif

	Wx2=.5*(x2sp+xr)-i2
	Wy2=.5*(y2sp+yr)-j2

#ifndef twoD
		Wz2=.5*(z2+zr)-k2
#endif
	
	Fx2=-q*(x2sp-xr)
	Fy2=-q*(y2sp-yr)
	Fz2=-q*(z2-zr)
	
#ifdef twoD
		Wz1=0
		Wz2=0
#endif

  onemWx1=1.-Wx1
  onemWx2=1.-Wx2
  onemWy1=1.-Wy1
  onemWy2=1.-Wy2
  onemWz1=1.-Wz1
  onemWz2=1.-Wz2

  i1p1=i1+1
  i2p1=i2+1
  j1p1=j1+1
  j2p1=j2+1
#ifndef twoD
  k1p1=k1+1
  k2p1=k2+1
#endif  

#ifdef cleancur
  
!  j(i1+j1*iy+k1*iz)=  j(i1+j1*iy+k1*iz)+Fx1 * onemWy1 * onemWz1

!if(1<0) then 
  curx(i1,j1,k1)=curx(i1,j1,k1)+Fx1 * onemWy1 * onemWz1
  curx(i1,j1p1,k1)=curx(i1,j1p1,k1)+Fx1 * Wy1 * onemWz1 
#ifndef twoD
  curx(i1,j1,  k1p1)= curx(i1,j1,  k1p1)+Fx1 * onemWy1 * Wz1
  curx(i1,j1p1,k1p1)= curx(i1,j1p1,k1p1)+Fx1 *  Wy1    * Wz1
#endif

  curx(i2,j2,k2)=curx(i2,j2,k2)+Fx2 * onemWy2 * onemWz2 
  curx(i2,j2p1,k2)=curx(i2,j2p1,k2)+Fx2 * Wy2 * onemWz2 
#ifndef twoD
  curx(i2,j2,  k2p1)= curx(i2,j2,  k2p1)+Fx2 * onemWy2* Wz2
  curx(i2,j2p1,k2p1)= curx(i2,j2p1,k2p1)+Fx2 *  Wy2    * Wz2
#endif

  cury(i1,j1,k1)=cury(i1,j1,k1)+Fy1 * onemWx1 * onemWz1 
  cury(i1p1,j1,k1)=cury(i1p1,j1,k1)+Fy1 * Wx1 * onemWz1  
#ifndef twoD
  cury(i1  ,j1,k1p1)= cury(i1  ,j1,k1p1)+Fy1 * onemWx1 * Wz1
  cury(i1p1,j1,k1p1)= cury(i1p1,j1,k1p1)+Fy1 *  Wx1    * Wz1
#endif
  cury(i2,j2,k2)=cury(i2,j2,k2)+Fy2 * onemWx2 * onemWz2 
  cury(i2p1,j2,k2)=cury(i2p1,j2,k2)+Fy2 *  Wx2    * onemWz2
#ifndef twoD
  cury(i2  ,j2,k2p1)= cury(i2  ,j2,k2p1)+Fy2 * onemWx2 * Wz2
  cury(i2p1,j2,k2p1)= cury(i2p1,j2,k2p1)+Fy2 *  Wx2    * Wz2
#endif
  curz(i1,j1,k1)=curz(i1,j1,k1)+ Fz1 * onemWx1 * onemWy1
  curz(i1p1,j1,k1)=curz(i1p1,j1,k1)+Fz1 * Wx1 * onemWy1
  curz(i1,j1p1,k1)=curz(i1,j1p1,k1)+Fz1 * onemWx1 * Wy1
  curz(i1p1,j1p1,k1)=curz(i1p1,j1p1,k1)+Fz1 * Wx1 * Wy1 
  curz(i2,j2,k2)=curz(i2,j2,k2)+ Fz2 * onemWx2 * onemWy2 
  curz(i2p1,j2,k2)=curz(i2p1,j2,k2)+ Fz2 * Wx2 * onemWy2 
  curz(i2,j2p1,k2)=curz(i2,j2p1,k2)+ Fz2 * onemWx2 * Wy2
  curz(i2p1,j2p1,k2)=curz(i2p1,j2p1,k2)+ Fz2 * Wx2 * Wy2 
!endif

#endif !cleancur

107 continue
  enddo !n=1,lecs



  pind(1:lecs)=0

  do n=1,lecs
     n1=n+maxhlf
     in=.true.
     perx=0.
     pery=0.
     perz=0.

     if(p(n1)%x .lt. minx .or. p(n1)%x .gt. maxx) then 
 		perx=sign(midx,p(n1)%x-minx)+sign(midx,p(n1)%x-maxx)    
     endif

     if(p(n1)%y .lt. miny .or. p(n1)%y .gt. maxy) then 
 		pery=sign(midy,p(n1)%y-miny)+sign(midy,p(n1)%y-maxy)    
     endif

     if(p(n1)%z .lt. minz .or. p(n1)%z .gt. maxz) then 
 		perz=sign(midz,p(n1)%z-minz)+sign(midz,p(n1)%z-maxz)    
     endif
		
		if(periodicx.eq.0) then      
			in=(p(n1)%x+mxcum .gt. x1in) .and. (p(n1)%x+mxcum .lt. x2in)
		endif
		
		if(periodicy.eq.0 .and. in) then 
			in=(p(n1)%y+mycum .gt. y1in) .and. (p(n1)%y+mycum .lt. y2in)
		endif
		
#ifndef twoD			
		if(periodicz.eq.0 .and. in) then      
			in=(p(n1)%z+mzcum .gt. z1in) .and. (p(n1)%z+mzcum .lt. z2in)
		endif		
#endif
  !to discard a particle
		if(.not. in) then
			perx=0
			pery=0
			perz=0
		endif		
		
		! to send to another processor
		if(perx.ne.0 .and. in .and. sizex.ne.1) then
			in=.false.	
			if(sizey .ne. 1) pery=0
			if(sizez .ne. 1) perz=0
		endif
		
		! assume non-uniform mx in general
		if(perx.lt.0 .and. sizex.ne.1) then		
			! minus rank (-x direction)		
			i1=(rank/sizex)*sizex + modulo(rank-1,sizex) 
			perx=-(mxl(i1+1)-5.)				
		endif
		
		p(n1)%x=p(n1)%x-perx		
		
		! to send to another processor
		if(pery.ne.0 .and. in .and. sizey.ne.1) then
			in=.false.	
			if(sizex .ne. 1) perx=0
			if(sizez .ne. 1) perz=0
		endif
		
		! assume non-uniform my
		if(pery.lt.0 .and. sizey.ne.1) then		
			! left rank	(-y direction)	
			j1=modulo(rank/sizex - 1,sizey)*sizex+rank/(sizex*sizey)*(sizex*sizey) &
				+ modulo(rank,sizex) 
			pery=-(myl(j1+1)-5.)				
		endif
		
		p(n1)%y=p(n1)%y-pery
		
#ifndef twoD
		! to send to another processor
		if(perz .ne. 0 .and. in ) then
			in=.false.	
			if(sizex .ne. 1) perx=0
			if(sizey .ne. 1) pery=0
		endif	
		! assume non-uniform mz
		if(perz .lt. 0) then		
			! down rank	(-z direction)	
			k1=modulo(rank/(sizex*sizey) - 1,sizez)*(sizex*sizey) + & 
			      modulo(rank,sizex*sizey) 
			perz=-(mzl(k1+1)-5.)				
		endif
#endif

		p(n1)%z=p(n1)%z-perz		
			
		if(in) go to 59
		nlecout=nlecout+1
		
#ifndef twoD		
		if(perz .lt. 0) then
			!this particle is going to the processor below
			LenLecOutDwn=LenLecOutDwn+1
			call copyprt(p(n1),poutdwn(LenIonOutDwn+LenLecOutDwn))
		endif

		if(perz .gt. 0) then
			!this particle is going to the processor above
			LenLecOutUp=LenLecOutUp+1
			call copyprt(p(n1),poutup(LenIonOutUp+LenLecOutUp))
                endif
#endif

		if(sizey .ne. 1) then
                if(pery .lt. 0) then
			!this particle is going to the processor on the left
			LenLecOutLft=LenLecOutLft+1
			call copyprt(p(n1),poutlft(LenIonOutLft+LenLecOutLft))
		endif

		if(pery .gt. 0) then
			!this particle is going to the processor on the right
			LenLecOutRgt=LenLecOutRgt+1
			call copyprt(p(n1),poutrgt(LenIonOutRgt+LenLecOutRgt))
                endif
		endif  !if(sizey .ne. 1)
		
		if(sizex .ne. 1) then
		if(perx .lt. 0) then
			!this particle is going to the processor on the left
			LenLecOutMinus=LenLecOutMinus+1
			call copyprt(p(n1),poutminus(LenIonOutMinus+LenLecOutMinus))
		endif

		if(perx .gt. 0) then
			!this particle is going to the processor on the right
			LenLecOutPlus=LenLecOutPlus+1
			call copyprt(p(n1),poutplus(LenIonOutPlus+LenLecOutPlus))
		endif
		endif	!if(sizex .ne. 1)
		
  pind(n)=1 !mark this particle for deletion; use n, not n1
59 continue
  
enddo !n=1,lecs
		if(debug) print *, rank,":","deposited lecs","step=",lap 
!now clean particles
  lecs0=lecs
  do n=1,lecs0
     if(pind(n) .ne. 0) then
        do while (pind(n) .ne. 0) 
           call copyprt(p(maxhlf+lecs),p(maxhlf+n))
           pind(n)=pind(lecs)
           pind(lecs)=0
           lecs=lecs-1
        enddo
     endif
  enddo		
	endif !if lecs .gt. 0
	
	nioneject=LenIonOutUp+LenIonOutDwn+LenIonOutLft+LenIonOutRgt+&
			  LenIonOutMinus+LenIonOutPlus
	
	nleceject=LenLecOutUp+LenLecOutDwn+LenLecOutLft+LenLecOutRgt+&
			  LenLecOutMinus+LenLecOutPlus
	
			  
	
end subroutine deposit_particles

#ifdef twoD
end module m_particles_movedeposit
#else
end module m_particles_movedeposit_3d
#endif
