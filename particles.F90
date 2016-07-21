!
! Particle module
!
! Includes the particle data structures and the routines for depositing
! and moving the particles
!
!

#ifdef twoD

module m_particles

	use m_globaldata
	use m_system
	use m_aux
	use m_communications
	use m_fields
	use m_inputparser
	use m_fparser
	
#else

module m_particles_3d

	use m_globaldata_3d
	use m_system_3d
	use m_aux_3d
	use m_communications_3d
	use m_fields_3d
	use m_inputparser_3d
	use m_fparser_3d

#endif


	
	implicit none
	
	private

!-------------------------------------------------------------------------------
!	PARAMETERS
!-------------------------------------------------------------------------------

	integer, parameter :: pdf_sz=1000
	
!-------------------------------------------------------------------------------
!	TYPE DEFINITIONS
!-------------------------------------------------------------------------------

	type :: particle
		sequence
		real(sprec) :: x, y, z, u, v, w, ch
		integer (kind=4) :: ind, proc, splitlev !need dum to keep size even # of bytes
	end type particle

	
	type :: prtnum
		sequence
		integer (kind=4) :: ind, proc
	end type prtnum
	
!-------------------------------------------------------------------------------
!	VARIABLES
!-------------------------------------------------------------------------------

	real(sprec) :: gamma0, ppc0, delgam, npnb, betshock, dummyvar, sigma, mi, me, qmi, &
				   qme, leftwall, acool, qi, qe, q, c_omp, movwinoffset
				   
    real(dprec) :: pi
    real(sprec) :: walloc
    
    integer :: upsamp_e, upsamp_i

	
	integer :: totalpartnum, shockloc, leftclean, nionout, nlecout, splitnmax,  &
			   splitratio, nioneject, nleceject, buffsize, movingshock, pcosthmult
			   
	integer, allocatable, dimension(:) :: pall 
	
	integer(8) :: maxptl0
	integer :: ions, lecs, maxhlf, lapreorder, maxptl, LenIonOutUp,LenIonOutDwn,LenLecOutUp, &
               LenLecOutDwn,LenIonInBlw,LenIonInAbv,LenLecInBlw, LenLecInAbv, LenIonOutLft, &
               LenIonOutRgt,LenLecOutLft, LenLecOutRgt, LenIonInLft, LenIonInRgt,LenLecInLft,LenLecInRgt, &
               LenIonOutPlus, LenIonOutMinus, LenIonInPlus, LenIonInMinus, &
               LenLecOutPlus, LenLecOutMinus, LenLecInPlus, LenLecInMinus, &
	       outcorner
      
	integer :: receivedions, receivedlecs, injectedions, injectedlecs

	! these variables for writing test particles in output.F90
	logical prt_first_lec, prt_first_ion
	
	
	real(dprec) :: time_cut     
	
	type(particle),allocatable :: p(:)
	type(particle),allocatable :: tempp(:) 
	type(particle),allocatable :: poutup(:), pinabv(:),poutdwn(:),pinblw(:) &
	,poutlft(:),poutrgt(:),pinlft(:),pinrgt(:) &
	,poutminus(:),poutplus(:),pinminus(:),pinplus(:)
	
 	integer, allocatable :: pind(:) !for optimized BC
	
	! Variables related to particle spliting
!	real, allocatable, dimension(:) :: split_E_ions, split_E_lecs
!	real splitoffs          ! maximal offset of child particles [in cells]
!	integer :: split_auto      ! if =1 then split level thresholds automatically adjust  
!	integer :: split_lap_start, split_lap_rate
	
!	integer ::split_prev_ions ! number of unsplit ions in last cycle split
!	integer ::split_prev_lecs ! number of unsplit ions in last cycle split
!	real split_frac_ions    ! fraction of split ions w.r.t total
!	real split_frac_lecs    ! fraction of split lecs w.r.t total
!	integer ::split_gamma_res ! resolution of energy histogram
!	real split_lgamma_m1_min ! min(log(gamma)-1) of histogram
!	real split_lgamma_m1_max ! max(log(gamma)-1) of histogram
!	real split_min_ratio    ! minimal ratio betwen consecutive split thresholds. 
!	real(dprec) :: split_timer
	
	
	logical :: external_fields, user_part_bcs

!-------------------------------------------------------------------------------
!	INTERFACE DECLARATIONS
!-------------------------------------------------------------------------------
	
!-------------------------------------------------------------------------------
!	PUBLIC MODIFIERS
!-------------------------------------------------------------------------------

	! public functions
	
        public :: reorder_particles, exchange_particles, zigzag, &
			  inject_others, allocate_particles, maxwell_dist,  powerlaw3d, &
			  read_input_particles, check_overflow, check_overflow_num, init_maxw_table, &
			  inject_from_wall, inject_plasma_region, copyprt 
!        public ::  split_from_wall, split_plasma_region, init_split_parts, split_E_ions, split_E_lecs, 


	! public variables 

	public :: ppc0,upsamp_e,upsamp_i, buffsize, poutup, pinabv, poutdwn, pinblw, &
			  poutlft, poutrgt, pinlft, pinrgt, &
			  poutminus, poutplus, pinminus, pinplus, &
			  pall, ions, lecs, particle, p, maxhlf, maxptl0, &           ! used in domain
			  leftwall, walloc, &							   ! used in restart
			  sigma, gamma0, delgam, mi, me, dummyvar, nionout, nlecout, & ! used in outputs
			  injectedions, injectedlecs, receivedions, receivedlecs, nioneject, nleceject, &      ! used in outputs
			  splitratio, qi, qe, q, prtnum, c_omp, maxptl, movingshock, &						   ! used in outputs
			  time_cut, prt_first_lec, prt_first_ion, betshock, leftclean, tempp, pind, &	!pind for optim		
! used in initialize
!			  gamma_table, pdf_table, gamma_table1, pdf_table1, gamma_table_cold, pdf_table_cold, &! used in initialize
!			  splitnmax, splitoffs, split_auto, split_lap_start, split_lap_rate, split_prev_ions, &! used in initialize
!			  split_prev_lecs, split_frac_ions, split_frac_lecs, split_gamma_res, &				   ! used in initialize
!			  split_lgamma_m1_min, split_lgamma_m1_max, split_min_ratio, split_timer, &! used in initialize
			  pdf_sz, npnb,qme, qmi, lapreorder, &											   				   ! used in initialize
			  totalpartnum, pcosthmult, LenIonOutUp, &		   ! used in the m_user_* module
			  LenIonOutDwn,LenLecOutUp, LenLecOutDwn,LenIonInBlw,LenIonInAbv,LenLecInBlw, LenLecInAbv, &        ! used in the m_user_* module
			  LenIonOutLft, LenIonOutRgt,LenLecOutLft, LenLecOutRgt, LenIonInLft, &   ! used in the m_user_* module
			  LenIonInRgt,LenLecInLft,LenLecInRgt, &
			  LenIonOutPlus, LenIonOutMinus, LenIonInPlus, LenIonInMinus, &
              LenLecOutPlus, LenLecOutMinus, LenLecInPlus, LenLecInMinus, &
			  shockloc, outcorner, & 										   ! used in the m_user_* module
			  movwinoffset, pi, external_fields, user_part_bcs
		
			  
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

subroutine read_input_particles()

	implicit none
	
	real(sprec) :: lmaxptl0
	real(sprec) :: omp

	call inputpar_getd_def("particles", "sigma", 0._sprec, sigma)
	call inputpar_getd_def("particles", "maxptl0", 1.e7_sprec, lmaxptl0)
	call inputpar_getd_def("particles", "ppc0", 1._sprec, ppc0)

	call inputpar_getd_def("particles", "delgam", 1.e-8_sprec, delgam)
	call inputpar_getd_def("particles", "me", 1._sprec, me)
	call inputpar_getd_def("particles", "mi", 20._sprec, mi)

	call inputpar_getd_def("particles", "gamma0", 1._sprec, gamma0)

	call inputpar_getd_def("particles", "c_omp", 0._sprec, c_omp)
	
	call inputpar_geti_def("problem", "upsamp_e", 1, upsamp_e)
	call inputpar_geti_def("problem", "upsamp_i", 1, upsamp_i)		
	
	
        !non-rel drift velocities can be specified with gamma<1 in input file. gamma = beta in that case. 
	if(gamma0 .lt. 1) gamma0=sqrt(1./(1.-gamma0**2)) 

	maxptl0=lmaxptl0
!	maxptl0=lmaxptl0*sizez	!weak scaling hack

	!plasma reaction
	omp=c/c_omp
	        
 	beta=sqrt(1.-1./gamma0**2) !drift velocity
	
	!note that gamma0 is used in the definition of omega_p to get charge. 
	qe=-(omp**2*gamma0)/((ppc0*.5)*(1+me/mi)) 
	!no sqrt because qe/me = 1 
	
	qi=-qe 
	
	me=me*abs(qi)
	mi=mi*abs(qi)
	
	qme=qe/me	! equal to -1/me in input
	qmi=qi/mi	! equal to  1/mi in input
 
       	Binit=sqrt(gamma0*ppc0*.5*c**2*(me*(1+me/mi))*sigma)

        pcosthmult=1 !by default, Maxwellians will be initialized in 3D (nonzero spread in 3 directions)

!	qi=0.
!	qe=0.

end subroutine read_input_particles

!-------------------------------------
!       Subroutine copyprt
!       copies one particle type to another explicitly, without relying on F90 intrinsics
!
!-------------------------------------
	subroutine copyprt(inp,outp)
	type(particle), intent(in)::inp
	type(particle), intent(out)::outp
	outp%x=inp%x
	outp%y=inp%y
	outp%z=inp%z
	outp%u=inp%u
	outp%v=inp%v
	outp%w=inp%w
	outp%ch=inp%ch
	outp%ind=inp%ind
	outp%proc=inp%proc
	outp%splitlev=inp%splitlev
	end subroutine copyprt


!-------------------------------------------------------------------------------
! 						subroutine allocate_particles				 
!																		
! Allocates the particle module's variables, some initialization of default vars
! 							
!-------------------------------------------------------------------------------

subroutine allocate_particles()

	implicit none
	
	! local variables
		
	pi=3.1415927
	
!	Binit_default=sqrt((gamma0-1)*.5*ppc0*c**2*(mi+me)*sigma) ! needed in fields module	
	
	!effectively me and mi = abs(qe) and qi, for use in computing energy density  

	leftwall=15.   !location of the left wall to reflect the stream
	movwinoffset = 0.

	!particle buffer size
	
	maxptl=maxptl0/size0
	maxhlf=maxptl/2
		
	! buffer size for particle communication

#ifndef twoD
	buffsize = max(10*int(ppc0*max(upsamp_e,upsamp_i)*c*max((mx-5)*(my-5),1*(mx-5)*(mz-5))),10000)
#else
!	buffsize = max(3*int(ppc0*c*(1*(mx-5))),60000) !5 for splitting particles
	buffsize = max(3*int(ppc0*max(upsamp_e,upsamp_i)*c*max((mx-5),(my-5))),60000)
!        if (splitparts) buffsize=buffsize*150 
#endif
	
	! allocate particle buffers
	
	allocate(p(maxptl),tempp(maxhlf))
        allocate(pind(maxhlf)) !for optimization
	p%x=0
	p%y=0
	p%z=0
	p%u=0
	p%v=0
	p%w=0
	p%ch=1
	p%splitlev=1
	
	allocate(poutup(buffsize),poutdwn(buffsize),pinblw(buffsize),pinabv(buffsize) &
	,poutlft(buffsize),poutrgt(buffsize),pinlft(buffsize),pinrgt(buffsize) &
	,poutminus(buffsize),poutplus(buffsize),pinminus(buffsize),pinplus(buffsize))
	
	! lot = mx*my for 2D ; total cells/processor
	allocate(pall(lot)) ! this is used in the particle module to reorder particles
 
        splitratio=10. !legacy splitting -- should not be used 
	splitnmax=15            ! max number of levels when splitting is used
		
	x1in=3.
	x2in=mx0-2.
	y1in=3.
	y2in=my0-2.
	z1in=3.
	z2in=mz0-2.

        totalpartnum=0
	
!	allocate(split_E_ions(splitnmax))
!	allocate(split_E_lecs(splitnmax))	
	
end subroutine allocate_particles




!-------------------------------------------------------------------------------
! 						subroutine check_overflow()
!														
! 
!
!-------------------------------------------------------------------------------
subroutine check_overflow()

	if(ions .gt. maxhlf .or. lecs .gt. maxhlf) then
		print *,rank,": Particle OVERFLOW: ions,lecs,maxhlf",ions &
		,lecs,maxhlf, "Increase the maxptl variable in input file and restart"
		stop
	endif

end subroutine check_overflow

!-------------------------------------------------------------------------------
! 						subroutine check_overflow_num()
!										
!
!-------------------------------------------------------------------------------
subroutine check_overflow_num(num)
	real num
	if(num .gt. maxhlf) then 
	   print *, rank, ": OVERFLOW: number of injected &
	   particles will exceed the maximum. Increase maxptl in &
	   input file and restart" , num, maxhlf
	   stop
	endif
end subroutine check_overflow_num

!-------------------------------------------------------------------------------
! 						subroutine reorder_particles()					 
!															
! 
!
!-------------------------------------------------------------------------------

subroutine reorder_particles()

	implicit none
	
	if(lap .eq. lapreorder+1) return
	
	if(modulo(lap,10) .ne. 0) return
	

	call reorder_particles_()	
	lapreorder=lap


end subroutine reorder_particles



!-------------------------------------------------------------------------------
! 						subroutine reorder_particles_()					 
!			
! 
!
!-------------------------------------------------------------------------------

subroutine reorder_particles_()

	implicit none
	
	! local variables
	
	integer ::n, celli, k,j, i

	if(Rank.eq.0)print *, "reordering particles"    
	pall=0
	!count the number of particles in each cell
 
	do n=1,ions
		celli=int(p(n)%x)+int(p(n)%y)*iy+int(p(n)%z)*iz !ix=1
		pall(celli)=pall(celli)+1
	enddo
	
	!convert Pall into an allocation
	!      if(rank.eq.0)print *,"done count"
	
	k=1
	do i=1,lot
		j=pall(i) 
		pall(i)=k
		k=k+j
	enddo
	
	if(debug) print *,rank,": done realloc"
	if(debug) print *, rank,":","ions=",ions

	do n=1,ions
		celli=int(p(n)%x)+int(p(n)%y)*iy+int(p(n)%z)*iz !ix=1
		j=pall(celli)
		pall(celli)=pall(celli)+1
		tempp(j)=p(n)
	enddo
	
	do n=1,ions
		p(n)=tempp(n)
	enddo
	
	!now for lecs
	!            if(Rank .eq. 0) print *, "now lecs"

	pall=0

	!count the number of particles in each cell

	if(debug) print *, rank,":","maxhlf",maxhlf,"lecs",lecs

	do n=maxhlf+1,maxhlf+lecs
		celli=int(p(n)%x)+int(p(n)%y)*iy+int(p(n)%z)*iz !ix=1
		pall(celli)=pall(celli)+1
	enddo
	
	!convert Pall into an allocation
	
	k=1
	
	do i=1,lot
		j=pall(i) 
		pall(i)=k
		k=k+j
	enddo
	
	!            if(Rank .eq. 0) print *,"done realloc"
	
	if(debug) print *,rank,": done realloc"

	do n=maxhlf+1,maxhlf+lecs
		celli=int(p(n)%x)+int(p(n)%y)*iy+int(p(n)%z)*iz !ix=1
		j=pall(celli)
		pall(celli)=pall(celli)+1
		tempp(j)=p(n)
	enddo
	
	do n=1,lecs
		p(maxhlf+n)=tempp(n)
	enddo
	
	if(debug) print *,rank,": done reorder"

end subroutine reorder_particles_



!-------------------------------------------------------------------------------
! 						subroutine index_particles()					 
!				
! Indexing particles, currently turned off
!
!-------------------------------------------------------------------------------

subroutine index_particles()

	implicit none
	
	! local variables

	integer ::n
	
	
!	if(rank.eq.0)print *,rank,":","indexing particles...",ions,lecs
!	
!	do n=1,ions
!		icell(n)=int(p(n)%x)+int(p(n)%y)*iy+int(p(n)%z)*iz !ix=1
!	enddo
!	
!	do n=maxhlf+1,maxhlf+lecs
!		icell(n)=int(p(n)%x)*ix +int(p(n)%y)*iy+int(p(n)%z)*iz
!	enddo
!	
!	
!	do n=1,ions
!		indprevrs(indp(n))=n
!	enddo
!	
!	do n=maxhlf+1,maxhlf+lecs
!		indprevrs(indp(n))=n
!	enddo
!	
!	if(rank.eq.0)print *,rank,":", "done indexing" 
	
end subroutine index_particles


!-------------------------------------------------------------------------------
! 				subroutine zigzag()					 
!																		
! Charge (current) deposition on the grid
!
!
!-------------------------------------------------------------------------------

subroutine zigzag(x2,y2,z2,x1,y1,z1,in)

	implicit none
	
	! dummy variables
	
	logical :: in
	real(sprec) :: x1,x2, y1, y2
	real(sprec) :: z1,z2
	
	! relay point (xr,yr,zr)
	real(sprec) :: xr, yr
	real zr
	! charge fluxes in 1st order scheme
	real Fx1, Fx2, Fy1, Fy2, Fz1, Fz2
	! shape functions
	real Wx1, Wx2, Wy1, Wy2, Wz1, Wz2
	integer ::i1, i2, j1, j2, k1, k2  

	! dummy variables
	integer :: s1, s2, s3
	! charge fluxes in 2nd order scheme
	real, dimension(2) :: flx, fly, flz
	real xp, yp, zp	
	! shape functions in 2nd order scheme
	real, dimension(3) :: Sx, Sy, Sz
	
	! zigzag method (1st order version) (Ref: Umeda et al. 2003)
	
	i1=aint(x1)
	i2=aint(x2)
	j1=aint(y1)
	j2=aint(y2)
	k1=aint(z1)
	k2=aint(z2)
		
		xr=min(real(min(i1,i2)+1),max(real(max(i1,i2)),.5*(x1+x2)))
		yr=min(real(min(j1,j2)+1),max(real(max(j1,j2)),.5*(y1+y2)))
		zr=min(real(min(k1,k2)+1),max(real(max(k1,k2)),.5*(z1+z2)))
	
#ifdef twoD
		k1=1
		k2=1
#endif

	!-q to include -j in the Ampere's equation, to be consistent
	
	Fx1=-q*(xr-x1)
	Fy1=-q*(yr-y1)
	Fz1=-q*(zr-z1)
	
	Wx1=.5*(x1+xr)-i1
	Wy1=.5*(y1+yr)-j1

#ifndef twoD
		Wz1=.5*(z1+zr)-k1
#endif

	Wx2=.5*(x2+xr)-i2
	Wy2=.5*(y2+yr)-j2

#ifndef twoD
		Wz2=.5*(z2+zr)-k2
#endif
	
	Fx2=-q*(x2-xr)
	Fy2=-q*(y2-yr)
	Fz2=-q*(z2-zr)
	
#ifdef twoD
		Wz1=0
		Wz2=0
#endif

	curx(i1,j1,  k1  )= curx(i1,j1,  k1  )+Fx1 * (1.-Wy1)*(1.-Wz1)
	curx(i1,j1+1,k1  )= curx(i1,j1+1,k1  )+Fx1 *  Wy1    *(1.-Wz1)

#ifndef twoD
		curx(i1,j1,  k1+1)= curx(i1,j1,  k1+1)+Fx1 * (1-Wy1) * Wz1
		curx(i1,j1+1,k1+1)= curx(i1,j1+1,k1+1)+Fx1 *  Wy1    * Wz1
#endif
	
	curx(i2,j2,  k2  )= curx(i2,j2,  k2  )+Fx2 * (1.-Wy2)*(1.-Wz2)
	curx(i2,j2+1,k2  )= curx(i2,j2+1,k2  )+Fx2 *  Wy2    *(1.-Wz2)
#ifndef twoD
		curx(i2,j2,  k2+1)= curx(i2,j2,  k2+1)+Fx2 * (1.-Wy2)* Wz2
		curx(i2,j2+1,k2+1)= curx(i2,j2+1,k2+1)+Fx2 *  Wy2    * Wz2
#endif

	cury(i1  ,j1,k1  )= cury(i1  ,j1,k1  )+Fy1 * (1.-Wx1)*(1.-Wz1)
	cury(i1+1,j1,k1  )= cury(i1+1,j1,k1  )+Fy1 *  Wx1    *(1.-Wz1) 
#ifndef twoD
		cury(i1  ,j1,k1+1)= cury(i1  ,j1,k1+1)+Fy1 * (1.-Wx1)* Wz1
		cury(i1+1,j1,k1+1)= cury(i1+1,j1,k1+1)+Fy1 *  Wx1    * Wz1
#endif

	cury(i2  ,j2,k2  )= cury(i2  ,j2,k2  )+Fy2 * (1.-Wx2)*(1.-Wz2)
	cury(i2+1,j2,k2  )= cury(i2+1,j2,k2  )+Fy2 *  Wx2    *(1.-Wz2) 
#ifndef twoD
		cury(i2  ,j2,k2+1)= cury(i2  ,j2,k2+1)+Fy2 * (1.-Wx2)* Wz2
		cury(i2+1,j2,k2+1)= cury(i2+1,j2,k2+1)+Fy2 *  Wx2    * Wz2
#endif

	curz(i1  ,j1  ,k1)= curz(i1  ,j1  ,k1)+Fz1 * (1.-Wx1)*(1.-Wy1)
	curz(i1+1,j1  ,k1)= curz(i1+1,j1  ,k1)+Fz1 *  Wx1    *(1.-Wy1)
	curz(i1  ,j1+1,k1)= curz(i1  ,j1+1,k1)+Fz1 * (1.-Wx1)* Wy1
	curz(i1+1,j1+1,k1)= curz(i1+1,j1+1,k1)+Fz1 *  Wx1    * Wy1
	
	curz(i2  ,j2  ,k2)= curz(i2  ,j2  ,k2)+Fz2 * (1.-Wx2)*(1.-Wy2)
	curz(i2+1,j2  ,k2)= curz(i2+1,j2  ,k2)+Fz2 *  Wx2    *(1.-Wy2)
	curz(i2  ,j2+1,k2)= curz(i2  ,j2+1,k2)+Fz2 * (1.-Wx2)* Wy2
	curz(i2+1,j2+1,k2)= curz(i2+1,j2+1,k2)+Fz2 *  Wx2    * Wy2
		

	in=.true.
		
end subroutine zigzag

!-------------------------------------------------------------------------------
! 						subroutine inject_others()					 
!																		
! Injects particles that cross from other processors
!
!-------------------------------------------------------------------------------
subroutine inject_others()

	implicit none

	real(dprec) :: time01, time02, time03, time04
	real perz, pery
	integer ::n, j1, k1

!	print *, rank,": in InjO, b4", "lap",lap, LenIonOutUp, LenIonOutDwn, LenLecOutUp, LenLecOutDwn

	LenIonOutPlus=0
	LenLecOutPlus=0
	LenIonOutMinus=0
	LenLecOutMinus=0	
			
	LenIonOutRgt=0
	LenLecOutRgt=0
	LenIonOutLft=0
	LenLecOutLft=0

	LenIonOutUp=0 !resetting these in anticipation of needing to send more
	LenIonOutDwn=0 !particles in z direction to accommodate corner crossing particles
	LenLecOutUp=0
	LenLecOutDwn=0
	
	receivedions=0
	receivedlecs=0
		
#ifndef twoD	
	do n=1,LenIonInAbv
		ions=ions+1
!		p(ions)=pinabv(n)
		call copyprt(pinabv(n),p(ions))
	enddo
	
	do n=1,LenLecInAbv
		lecs=lecs+1
!		p(maxhlf+lecs)=pinabv(LenIonInAbv+n)
		call copyprt(pinabv(LenIonInAbv+n),p(maxhlf+lecs))
	enddo

	!those that entered from below

	do n=1,LenIonInBlw
		ions=ions+1
!		p(ions)=pinblw(n)
		call copyprt(pinblw(n),p(ions))
	enddo
	
	do n=1,LenLecInBlw
		lecs=lecs+1
!		p(maxhlf+lecs)=pinblw(LenIonInBlw+n)
		call copyprt(pinblw(LenIonInBlw+n),p(maxhlf+lecs))
	enddo
#endif
	

	if(sizey .ne. 1) then
			
	!now inject particles that entered from left and right in y		
	do n=1,LenIonInRgt
		ions=ions+1
!		p(ions)=pinrgt(n)
		call copyprt(pinrgt(n),p(ions))

#ifndef twoD
!#ifdef CORNER


!hack:  to take care of particles that go through corners
		perz=sign(.5*(mz-5.),p(ions)%z-3.)+sign(.5*(mz-5.),p(ions)%z-mz +2.)
!hack: need to add special rules for perz in case last mz is not equal to mzall
!	if(perz .ne. 0) then 
!	   print *, rank," :", "perz ne 0, ir", lap,n,ions, perz,p(ions)%x, p(ions)%y, p(ions)%z,p(ions)%z-mz +2.,p(ions)%z-3.
!	endif

		if(perz .lt. 0) then		
			! down rank	(-z direction)	
			k1=modulo(rank/(sizex*sizey) - 1,sizez)*(sizex*sizey) + & 
			      modulo(rank,sizex*sizey) 
			perz=-(mzl(k1+1)-5.)				
		endif		
		
        p(ions)%z=p(ions)%z-perz   
		if(perz .lt. 0) then
		!this particle is going to the processor below
		      LenIonOutDwn=LenIonOutDwn+1
		     ! poutdwn(LenIonOutDwn)=p(ions)           
		      call copyprt(p(ions),poutdwn(LenIonOutDwn))
!		   if(periodicz .eq. 0 .and. & 
!                      ((p(ions)%z +perz+ rank/sizey*(mzall-5).lt.z1in) &
!                       .or.(p(ions)%z+perz + rank/sizey*(mzall-5).gt.z2in)) ) then
!			 LenIonOutDwn=LenIonOutDwn-1  !let the particle escape 
!		   endif
		endif
	
		if(perz .gt. 0) then
			!this particle is going to the processor above
			LenIonOutUp=LenIonOutUp+1
		   !	poutup(LenIonOutUp)=p(ions)
			call copyprt(p(ions),poutup(LenIonOutUp))
!		   if ( periodicz .eq. 0 .and. & 
!                       ((p(ions)%z+perz + rank/sizey*(mzall-5).lt.z1in) &
!                     .or.(p(ions)%z+perz + rank/sizey*(mzall-5).gt.z2in)) ) then
!			 LenIonOutUp=LenIonOutUp-1  !let the particle escape 
!		   endif
		endif
	if(perz .ne. 0) then	
		ions=ions-1	
	endif	
 !end hacking
!end ifdef CORNER
!#endif
#endif
	enddo  !do n=1,LenIonInRgt
	
	!those that entered from left

	do n=1,LenIonInLft
		ions=ions+1
!		p(ions)=pinlft(n)
		call copyprt(pinlft(n),p(ions))

#ifndef twoD
!#ifdef CORNER
!hacking:  to take care of particles that go through corners
	perz=sign(.5*(mz-5.),p(ions)%z-3.)+sign(.5*(mz-5.),p(ions)%z-mz +2.)
!hack: need to add special rules for perz in case last mz is not equal to mzall
!	if(perz .ne. 0) then 
!	   print *, rank," :", "perz ne 0, il", lap,n,ions, perz,p(ions)%x, p(ions)%y, p(ions)%z,p(ions)%ind, p(ions)%proc
!	endif		
		if(perz .lt. 0) then		
			! down rank	(-z direction)	
			k1=modulo(rank/(sizex*sizey) - 1,sizez)*(sizex*sizey) + & 
			modulo(rank,sizex*sizey) 
			perz=-(mzl(k1+1)-5.)				
		endif
        p(ions)%z=p(ions)%z-perz   
		if(perz .lt. 0) then
		!this particle is going to the processor below
		      LenIonOutDwn=LenIonOutDwn+1
!		      poutdwn(LenIonOutDwn)=p(ions)           
		      call copyprt(p(ions),poutdwn(LenIonOutDwn))

!		   if ( periodicz .eq. 0 .and. & 
!                        ((p(ions)%z +perz+ rank/sizey*(mzall-5).lt.z1in) &
!                         .or.(p(ions)%z+perz + rank/sizey*(mzall-5).gt.z2in)) ) then
!			 LenIonOutDwn=LenIonOutDwn-1  !let the particle escape 
!		   endif
		endif
	
		if(perz .gt. 0) then
			!this particle is going to the processor above
			LenIonOutUp=LenIonOutUp+1
		!	poutup(LenIonOutUp)=p(ions)
			call copyprt(p(ions),poutup(LenIonOutUp))
			
!		   if ( periodicz .eq. 0 .and. & 
!                       ((p(ions)%z+perz + rank/sizey*(mzall-5).lt.z1in) &
!                     .or.(p(ions)%z+perz + rank/sizey*(mzall-5).gt.z2in)) ) then
!			 LenIonOutUp=LenIonOutUp-1  !let the particle escape 
!		   endif
		endif
	if(perz .ne. 0) then	
		ions=ions-1	
	endif	
 !end hacking
!end ifdef CORNER
!#endif 
#endif
	enddo	!do n=1,LenIonInLft

	do n=1,LenLecInRgt

		lecs=lecs+1
!		p(maxhlf+lecs)=pinrgt(LenIonInRgt+n)
		call copyprt(pinrgt(LenIonInRgt+n),p(maxhlf+lecs))


#ifndef twoD
!#ifdef CORNER
!hacking:  to take care of particles that go through corners
	perz=sign(.5*(mz-5.),p(maxhlf+lecs)%z-3.)+sign(.5*(mz-5.),p(maxhlf+lecs)%z- mz +2.)
!hack: need to add special rules for perz in case last mz is not equal to mzall
 
		if(perz .lt. 0) then		
			! down rank	(-z direction)	
			k1=modulo(rank/(sizex*sizey) - 1,sizez)*(sizex*sizey) + & 
			      modulo(rank,sizex*sizey) 
			perz=-(mzl(k1+1)-5.)				
		endif
		
		p(maxhlf+lecs)%z=p(maxhlf+lecs)%z-perz   

!	if(perz .ne. 0) then 
!	   print *, rank," :", "perz ne 0, er", lap,n,lecs, perz,p(maxhlf+lecs)%x, p(maxhlf+lecs)%y, p(maxhlf+lecs)%z,p(maxhlf+lecs)%ind, p(maxhlf+lecs)%proc
!	endif

		if(perz .lt. 0) then
		!this particle is going to the processor below
		      LenLecOutDwn=LenLecOutDwn+1
	!	      poutdwn(LenIonOutDwn+LenLecOutDwn)=p(maxhlf+lecs)           
		     call copyprt(p(maxhlf+lecs),poutdwn(LenIonOutDwn+LenLecOutDwn))

!		   if ( periodicz .eq. 0 .and. & 
!                       ((p(maxhlf+lecs)%z +perz+ rank/sizey*(mzall-5).lt.z1in) &
!                         .or.(p(maxhlf+lecs)%z+perz + rank/sizey*(mzall-5).gt.z2in)) ) then
!			 LenLecOutDwn=LenLecOutDwn-1  !let the particle escape 
!		   endif
		endif
	
		if(perz .gt. 0) then
			!this particle is going to the processor above
			LenLecOutUp=LenLecOutUp+1
!			poutup(LenIonOutUp+LenLecOutUp)=p(maxhlf+lecs)
			call copyprt(p(maxhlf+lecs),poutup(LenIonOutUp+LenLecOutUp))
			
!		   if ( periodicz .eq. 0 .and. & 
!                        ((p(maxhlf+lecs)%z+perz + rank/sizey*(mzall-5).lt.z1in) &
!                     .or.(p(maxhlf+lecs)%z+perz+ rank/sizey*(mzall-5).gt.z2in)) ) then
!			 LenLecOutUp=LenLecOutUp-1  !let the particle escape 
!		   endif
		endif
	if(perz .ne. 0) then	
		lecs=lecs-1	
	endif	
!end ifdef CORNER
!#endif 
#endif
	enddo !do n=1,LenLecInRgt
	
	do n=1,LenLecInLft
		lecs=lecs+1
!		p(maxhlf+lecs)=pinlft(LenIonInLft+n)
		call copyprt(pinlft(LenIonInLft+n),p(maxhlf+lecs))
#ifndef twoD
!#ifdef CORNER
!hacking:  to take care of particles that go through corners
	perz=sign(.5*(mz-5.),p(maxhlf+lecs)%z-3.)+sign(.5*(mz-5.),p(maxhlf+lecs)%z- mz +2.)
!hack: need to add special rules for perz in case last mz is not equal to mzall

		if(perz .lt. 0) then		
			! down rank	(-z direction)	
			k1=modulo(rank/(sizex*sizey) - 1,sizez)*(sizex*sizey) + & 
			      modulo(rank,sizex*sizey) 
			perz=-(mzl(k1+1)-5.)				
		endif		
		
!	if(perz .ne. 0) then 
!	   print *, rank," :", "perz ne 0, el", lap,n,lecs, perz,p(maxhlf+lecs)%x, p(maxhlf+lecs)%y, p(maxhlf+lecs)%z,p(maxhlf+lecs)%ind, p(maxhlf+lecs)%proc
!	endif

        p(maxhlf+lecs)%z=p(maxhlf+lecs)%z-perz   
		if(perz .lt. 0) then
		!this particle is going to the processor below
		      LenLecOutDwn=LenLecOutDwn+1
!		      poutdwn(LenIonOutDwn+LenLecOutDwn)=p(maxhlf+lecs)           
		      call copyprt(p(maxhlf+lecs),poutdwn(LenIonOutDwn+LenLecOutDwn))

!		   if(periodicz .eq. 0 .and. & 
!                       ((p(maxhlf+lecs)%z +perz+ rank/sizey*(mzall-5).lt.z1in) &
!                       .or.(p(maxhlf+lecs)%z+perz + rank/sizey*(mzall-5).gt.z2in)) ) then
!			 LenLecOutDwn=LenLecOutDwn-1  !let the particle escape 
!		   endif
		endif
	
		if(perz .gt. 0) then
			!this particle is going to the processor above
			LenLecOutUp=LenLecOutUp+1
!			poutup(LenIonOutUp+LenLecOutUp)=p(maxhlf+lecs)
			call copyprt(p(maxhlf+lecs),poutup(LenIonOutUp+LenLecOutUp))

!		   if(periodicz .eq. 0 .and. & 
!                        ((p(maxhlf+lecs)%z+perz + rank/sizey*(mzall-5).lt.z1in) &
!                     .or.(p(maxhlf+lecs)%z+perz+ rank/sizey*(mzall-5).gt.z2in)) ) then
!			 LenLecOutUp=LenLecOutUp-1  !let the particle escape 
!		   endif
		endif
	if(perz .ne. 0) then	
		lecs=lecs-1	
	endif	
!end ifdef CORNER
!#endif 
#endif
	enddo  !do n=1,LenLecInLft
	
	endif	!if(sizey .ne. 1)
	
	if(sizex .ne. 1) then

	do n=1,LenIonInPlus
		ions=ions+1
		call copyprt(pinplus(n),p(ions))
		
!#ifdef CORNER
         if(sizey .ne. 1) then 
		pery=sign(.5*(my-5.),p(ions)%y-3.)+sign(.5*(my-5.),p(ions)%y-my+2.)
		if(pery .lt. 0) then		
			! left rank	(-y direction)	
			j1=modulo(rank/sizex - 1,sizey)*sizex+rank/(sizex*sizey)*(sizex*sizey) &
				+ modulo(rank,sizex) 
			pery=-(myl(j1+1)-5.)				
		endif
		p(ions)%y=p(ions)%y-pery
		if(pery .lt. 0) then
			!this particle is going to the processor on the left
			LenIonOutLft=LenIonOutLft+1
!			poutlft(LenIonOutLft)=p(n1)
			call copyprt(p(ions),poutlft(LenIonOutLft))
!		   if(periodicy .eq. 0 .and. & 
!                      ((p(n1)%y+pery + mycum .lt. y1in) &
!                         .or.(p(n1)%y+pery+ mycum .gt. y2in)) ) then
!			 LenIonOutLft=LenIonOutLft-1  !let the particle escape 
!		   endif
		endif
		if(pery .gt. 0) then
			!this particle is going to the processor on the right
			LenIonOutRgt=LenIonOutRgt+1
!			poutrgt(LenIonOutRgt)=p(n1)
			call copyprt(p(ions),poutrgt(LenIonOutRgt))
!		   if (periodicy .eq. 0 .and. & 
!                     ((p(n1)%y+pery + mycum .lt. y1in) &
!                  .or.(p(n1)%y+pery+ mycum .gt. y2in)) ) then
!			 LenIonOutRgt=LenIonOutRgt-1  !let the particle escape 
!		   endif
		endif
	if(pery .ne. 0) then	
		ions=ions-1	
	endif	
!end hacking
!end ifdef CORNER
!#endif 
        endif !if sizey ne 1
	enddo !do n=1,LenIonInPlus

	do n=1,LenIonInMinus
		ions=ions+1
		call copyprt(pinminus(n),p(ions))
		
!#ifdef CORNER
           if(sizey .ne. 1) then 
		pery=sign(.5*(my-5.),p(ions)%y-3.)+sign(.5*(my-5.),p(ions)%y-my+2.)
		if(pery .lt. 0) then		
			! left rank	(-y direction)	
			j1=modulo(rank/sizex - 1,sizey)*sizex+rank/(sizex*sizey)*(sizex*sizey) &
				+ modulo(rank,sizex) 
			pery=-(myl(j1+1)-5.)				
		endif
		p(ions)%y=p(ions)%y-pery
		if(pery .lt. 0) then
			!this particle is going to the processor on the left
			LenIonOutLft=LenIonOutLft+1
!			poutlft(LenIonOutLft)=p(n1)
			call copyprt(p(ions),poutlft(LenIonOutLft))
!		   if(periodicy .eq. 0 .and. & 
!                      ((p(n1)%y+pery + mycum .lt. y1in) &
!                         .or.(p(n1)%y+pery+ mycum .gt. y2in)) ) then
!			 LenIonOutLft=LenIonOutLft-1  !let the particle escape 
!		   endif
		endif

		if(pery .gt. 0) then
			!this particle is going to the processor on the right
			LenIonOutRgt=LenIonOutRgt+1
!			poutrgt(LenIonOutRgt)=p(n1)
			call copyprt(p(ions),poutrgt(LenIonOutRgt))
!		   if (periodicy .eq. 0 .and. & 
!                     ((p(n1)%y+pery + mycum .lt. y1in) &
!                  .or.(p(n1)%y+pery+ mycum .gt. y2in)) ) then
!			 LenIonOutRgt=LenIonOutRgt-1  !let the particle escape 
!		   endif
		endif
	if(pery .ne. 0) then	
		ions=ions-1	
	endif	
 !end hacking
 !end ifdef CORNER
!#endif 		
        endif !if sizey ne 1		
	enddo	!do n=1,LenIonInMinus
	
	do n=1,LenLecInPlus
		lecs=lecs+1
		call copyprt(pinplus(LenIonInPlus+n),p(maxhlf+lecs))
		
!#ifdef CORNER
         if(sizey .ne. 1) then      
		pery=sign(.5*(my-5.),p(maxhlf+lecs)%y-3.)+sign(.5*(my-5.),p(maxhlf+lecs)%y-my+2.)
		if(pery .lt. 0) then		
			! left rank	(-y direction)	
			j1=modulo(rank/sizex - 1,sizey)*sizex+rank/(sizex*sizey)*(sizex*sizey) &
				+ modulo(rank,sizex) 
			pery=-(myl(j1+1)-5.)				
		endif
		p(maxhlf+lecs)%y=p(maxhlf+lecs)%y-pery
		if(pery .lt. 0) then
			!this particle is going to the processor on the left
			LenLecOutLft=LenLecOutLft+1
		!	poutlft(LenIonOutLft+LenLecOutLft)=p(n1)
			call copyprt(p(maxhlf+lecs),poutlft(LenIonOutLft+LenLecOutLft))
!		   if(periodicy .eq. 0 .and. & 
!                   ((p(n1)%y + pery + mycum .lt. y1in) &
!                         .or.(p(n1)%y+ pery + mycum .gt. y2in)) ) then
!                   LenLecOutLft=LenLecOutLft-1  !let the particle escape 
!		   endif
		endif
		if(pery .gt. 0) then
			!this particle is going to the processor on the right
			LenLecOutRgt=LenLecOutRgt+1
		!	poutrgt(LenIonOutRgt+LenLecOutRgt)=p(n1)
			call copyprt(p(maxhlf+lecs),poutrgt(LenIonOutRgt+LenLecOutRgt))
			
!		   if(periodicy .eq. 0 .and. & 
!                    ((p(n1)%y + pery + mycum .lt. y1in) &
!                     .or.(p(n1)%y+ pery+ mycum .gt. y2in)) ) then
!			 LenLecOutRgt=LenLecOutRgt-1  !let the particle escape 
!		   endif
		endif
	if(pery .ne. 0) then	
		lecs=lecs-1	
	endif	
!end hacking
!end ifdef CORNER
!#endif 
        endif !if sizey ne 1		
	enddo !do n=1,LenLecInPlus
	
	do n=1,LenLecInMinus
		lecs=lecs+1
		call copyprt(pinminus(LenIonInMinus+n),p(maxhlf+lecs))
		
!#ifdef CORNER
          if(sizey .ne. 1) then 
		pery=sign(.5*(my-5.),p(maxhlf+lecs)%y-3.)+sign(.5*(my-5.),p(maxhlf+lecs)%y-my+2.)
		if(pery .lt. 0) then		
			! left rank	(-y direction)	
			j1=modulo(rank/sizex - 1,sizey)*sizex+rank/(sizex*sizey)*(sizex*sizey) &
				+ modulo(rank,sizex) 
			pery=-(myl(j1+1)-5.)				
		endif
		p(maxhlf+lecs)%y=p(maxhlf+lecs)%y-pery
		if(pery .lt. 0) then
			!this particle is going to the processor on the left
			LenLecOutLft=LenLecOutLft+1
		!	poutlft(LenIonOutLft+LenLecOutLft)=p(n1)
			call copyprt(p(maxhlf+lecs),poutlft(LenIonOutLft+LenLecOutLft))
!		   if(periodicy .eq. 0 .and. & 
!                   ((p(n1)%y + pery + mycum .lt. y1in) &
!                         .or.(p(n1)%y+ pery + mycum .gt. y2in)) ) then
!                   LenLecOutLft=LenLecOutLft-1  !let the particle escape 
!		   endif
		endif
		if(pery .gt. 0) then
			!this particle is going to the processor on the right
			LenLecOutRgt=LenLecOutRgt+1
		!	poutrgt(LenIonOutRgt+LenLecOutRgt)=p(n1)
			call copyprt(p(maxhlf+lecs),poutrgt(LenIonOutRgt+LenLecOutRgt))
			
!		   if(periodicy .eq. 0 .and. & 
!                    ((p(n1)%y + pery + mycum .lt. y1in) &
!                     .or.(p(n1)%y+ pery+ mycum .gt. y2in)) ) then
!			 LenLecOutRgt=LenLecOutRgt-1  !let the particle escape 
!		   endif
		endif
	if(pery .ne. 0) then	
		lecs=lecs-1	
	endif	
!end hacking
!end ifdef CORNER
!#endif 
         endif !if sizey ne 1		
	enddo	!do n=1,LenLecInMinus

	endif !if(sizex .ne. 1)
	
	if(rank.eq.0) time03=mpi_wtime()
	
	receivedions=lenioninabv+lenioninblw+lenioninrgt+lenioninlft+&
				 lenioninplus+lenioninminus
	receivedlecs=lenlecinabv+lenlecinblw+lenlecinrgt+lenlecinlft+&
				 lenlecinplus+lenlecinminus

!	print *, rank,": in InjO, a4", "lap",lap,  LenIonOutUp, LenIonOutDwn, LenLecOutUp, LenLecOutDwn	

end subroutine inject_others





!-------------------------------------------------------------------------------
! 						subroutine exchange_particles()					 
!																		
! Exchanges particles between diferent processors
!
!-------------------------------------------------------------------------------

subroutine exchange_particles()

	implicit none
	
	! local variables
	
	integer ::arrsizeup(2),arrsizedwn(2), lenup, lendwn, ierr
	integer ::arrsizeblw(2),arrsizeabv(2),status(statsize)
	integer ::uprank, dwnrank,uptag, dwntag, lenabv, lenblw
	
	integer ::arrsizeoutrgt(2),arrsizeoutlft(2), lenoutrgt, lenoutlft
	integer ::arrsizeinlft(2),arrsizeinrgt(2)
	integer ::rgtrank, lftrank,rgttag, lfttag, leninrgt, leninlft
	
	integer ::arrsizeoutplus(2),arrsizeoutminus(2), lenoutplus, lenoutminus
	integer ::arrsizeinminus(2),arrsizeinplus(2)
	integer ::plusrank, minusrank, plustag, minustag, leninplus, leninminus
	
	arrsizeup(1)=LenIonOutUp
	arrsizeup(2)=LenLecOutUp
	lenup=max((arrsizeup(1)+arrsizeup(2)),1)
	
	arrsizedwn(1)=LenIonOutDwn
	arrsizedwn(2)=LenLecOutDwn
	lendwn=max((arrsizedwn(1)+arrsizedwn(2)),1)

#ifndef twoD
	uprank=modulo((rank/(sizex*sizey) + 1),sizez)*(sizex*sizey) + & 
		   modulo(rank,sizex*sizey) 
	dwnrank=modulo((rank/(sizex*sizey) - 1),sizez)*(sizex*sizey) + & 
		   modulo(rank,sizex*sizey) 
	! in 2D, uprank=modulo(rank,sizey), dwnrank=modulo(rank,sizey)
	uptag=100
	dwntag=200

	!send the info about the size of the array
	!send up and recv from below

#ifdef MPI      
		call MPI_SendRecv(arrsizeup,2,mpi_integer,uprank,uptag, &
		arrsizeblw,2,mpi_integer,dwnrank,uptag, &
		MPI_Comm_World,status,ierr)
#else
		call MPI_SendRecv(arrsizeup,2,mpi_integer,uprank,uptag, &
		arrsizeblw,2,mpi_integer,dwnrank,uptag, &
		MPI_Comm_World,status,ierr)
#endif
	
	!send dwn and recv from above
#ifdef MPI
		call MPI_SendRecv(arrsizedwn,2,mpi_integer,dwnrank,dwntag, &
		arrsizeabv,2,mpi_integer,uprank,dwntag, &
		MPI_Comm_World,status,ierr)
#else
		call MPI_SendRecv(arrsizedwn,2,mpi_integer,dwnrank,dwntag, &
		arrsizeabv,2,mpi_integer,uprank,dwntag, &
		MPI_Comm_World,status,ierr)
#endif

	LenIonInBlw=arrsizeblw(1)
	LenLecInBlw=arrsizeblw(2)
	
	LenIonInAbv=arrsizeabv(1)
	LenLecInAbv=arrsizeabv(2)
	
	lenblw=max((LenIonInBlw+LenLecInBlw),1)
	lenabv=max((LenIonInAbv+LenLecInAbv),1)
	
	if(lenblw .gt. buffsize .or. lenabv .gt. buffsize) then
		print*, "ERROR: BUFFERSIZE SMALLER THAN PARTICLE LIST"
		print*,rank,":","buffsize", buffsize, "lenblw",lenblw,"lenabv" &
		,lenabv
		stop
	endif

	!send and receive particle arrays
	!send up and recv from below
#ifdef MPI
		call MPI_SendRecv(poutup,lenup,particletype,uprank,uptag, &
		pinblw,lenblw,particletype,dwnrank,uptag, &
		MPI_Comm_World,status,ierr)
#else
		pinblw=poutup	
#endif

	!send dwn and recv from above
	!      print *,"lenup", lendwn, lenabv
	
#ifdef MPI
		call MPI_SendRecv(poutdwn,lendwn,particletype,dwnrank,dwntag, &
		pinabv,lenabv,particletype,uprank,dwntag, &
		MPI_Comm_World,status,ierr)
#else
		pinabv=poutdwn
#endif
!end ifdef twoD	
#endif 

	!----------------------------------------------------------
	! now send and receive particles left and right along the y direction
	!----------------------------------------------------------
	
	arrsizeoutrgt(1)=LenIonOutRgt
	arrsizeoutrgt(2)=LenLecOutRgt
	lenoutrgt=max((arrsizeoutrgt(1)+arrsizeoutrgt(2)),1)
	
	arrsizeoutlft(1)=LenIonOutLft
	arrsizeoutlft(2)=LenLecOutLft
	lenoutlft=max((arrsizeoutlft(1)+arrsizeoutlft(2)),1)

	if(sizey .ne. 1) then
	
	rgtrank=modulo(rank/sizex + 1,sizey)*sizex+rank/(sizex*sizey)*(sizex*sizey) &
				+ modulo(rank,sizex) 
	lftrank=modulo(rank/sizex - 1,sizey)*sizex+rank/(sizex*sizey)*(sizex*sizey) &
				+ modulo(rank,sizex) 

	rgttag=100
	lfttag=200
	
	!send the info about the size of the array
	!send up and recv from below
	if(debug)print *,rank,": outrgt",arrsizeoutrgt

#ifdef MPI
		call MPI_SendRecv(arrsizeoutrgt,2,mpi_integer,rgtrank,rgttag, &
		arrsizeinlft,2,mpi_integer,lftrank,rgttag, &
		MPI_Comm_World,status,ierr)
#else
		call MPI_SendRecv(arrsizeoutrgt,2,mpi_integer,rgtrank,rgttag, &
		arrsizeinlft,2,mpi_integer,lftrank,rgttag, &
		MPI_Comm_World,status,ierr)
#endif

	!send dwn and recv from above
	if(debug)print *,rank,": outlft",arrsizeoutlft

#ifdef MPI
		call MPI_SendRecv(arrsizeoutlft,2,mpi_integer,lftrank,lfttag, &
		arrsizeinrgt,2,mpi_integer,rgtrank,lfttag, &
		MPI_Comm_World,status,ierr)
#else
		call MPI_SendRecv(arrsizeoutlft,2,mpi_integer,lftrank,lfttag, &
		arrsizeinrgt,2,mpi_integer,rgtrank,lfttag, &
		MPI_Comm_World,status,ierr)
#endif
	
	LenIonInLft=arrsizeinlft(1)
	LenLecInLft=arrsizeinlft(2)
	
	LenIonInRgt=arrsizeinrgt(1)
	LenLecInRgt=arrsizeinrgt(2)
	
	leninlft=max((LenIonInLft+LenLecInLft),1)
	leninrgt=max((LenIonInRgt+LenLecInRgt),1)
	
	if(leninlft .gt. buffsize .or. leninrgt .gt. buffsize) then
		print *, "ERROR: BUFFERSIZE SMALLER THAN PARTICLE LIST"
		print *,rank,":","buffsize", buffsize, "leninlft",leninlft &
		,"leninrgt",leninrgt
		stop
	endif

	!send and receive particle arrays
	!send up and recv from below
	
	if(debug)print *,rank,"lenoutrgt",lenoutrgt,"leninlft",leninlft
#ifdef MPI
		call MPI_SendRecv(poutrgt,lenoutrgt,particletype,rgtrank,rgttag, &
		pinlft,leninlft,particletype,lftrank,rgttag, &
		MPI_Comm_World,status,ierr)
#else
		pinlft=poutrgt
#endif
	
	
#ifdef MPI   
		call MPI_SendRecv(poutlft,lenoutlft,particletype,lftrank,lfttag, &
		pinrgt,leninrgt,particletype,rgtrank,lfttag, &
		MPI_Comm_World,status,ierr)
#else
		pinrgt=poutlft
#endif

	endif ! if(sizey .ne. 1)

	!----------------------------------------------------------
	! now send and receive particles plus and minus (along the x direction)
	!----------------------------------------------------------
	
	arrsizeoutplus(1)=LenIonOutPlus
	arrsizeoutplus(2)=LenLecOutPlus
	lenoutplus=max((arrsizeoutplus(1)+arrsizeoutplus(2)),1)
	
	arrsizeoutminus(1)=LenIonOutMinus
	arrsizeoutminus(2)=LenLecOutMinus
	lenoutminus=max((arrsizeoutminus(1)+arrsizeoutminus(2)),1)
	
	if(sizex .ne. 1) then
	
	plusrank=(rank/sizex)*sizex + modulo(rank+1,sizex) 
	minusrank=(rank/sizex)*sizex + modulo(rank-1,sizex) 
	plustag=100
	minustag=200
	
	!send the info about the size of the array
	!send up and recv from below
	if(debug)print *,rank,": outplus",arrsizeoutplus

	call MPI_SendRecv(arrsizeoutplus,2,mpi_integer,plusrank,plustag, &
		arrsizeinminus,2,mpi_integer,minusrank,plustag, &
		MPI_Comm_World,status,ierr)

	!send dwn and recv from above
	if(debug)print *,rank,": outminus",arrsizeoutminus


	call MPI_SendRecv(arrsizeoutminus,2,mpi_integer,minusrank,minustag, &
		arrsizeinplus,2,mpi_integer,plusrank,minustag, &
		MPI_Comm_World,status,ierr)

	
	LenIonInMinus=arrsizeinminus(1)
	LenLecInMinus=arrsizeinminus(2)
	
	LenIonInPlus=arrsizeinplus(1)
	LenLecInPlus=arrsizeinplus(2)
	
	leninminus=max((LenIonInMinus+LenLecInMinus),1)
	leninplus=max((LenIonInPlus+LenLecInPlus),1)
	
	if(leninminus .gt. buffsize .or. leninplus .gt. buffsize) then
		print *, "ERROR: BUFFERSIZE SMALLER THAN PARTICLE LIST"
		print *,rank,":","buffsize", buffsize, "leninminus",leninminus &
		,"leninplus",leninplus
		stop
	endif

	!send and receive particle arrays
	!send up and recv from below
	
	if(debug)print *,rank,"lenoutplus",lenoutplus,"leninminus",leninminus

	call MPI_SendRecv(poutplus,lenoutplus,particletype,plusrank,plustag, &
		pinminus,leninminus,particletype,minusrank,plustag, &
		MPI_Comm_World,status,ierr)
   
	call MPI_SendRecv(poutminus,lenoutminus,particletype,minusrank,minustag, &
		pinplus,leninplus,particletype,plusrank,minustag, &
		MPI_Comm_World,status,ierr)

	endif ! if(sizex .ne.1)
	
end subroutine exchange_particles


!-------------------------------------------------------------------------------
! 						subroutine init_maxw_table() 
!																		
! Initialize probability distribution tables for quick maxwellian generation
!
!-------------------------------------------------------------------------------

subroutine init_maxw_table(delgam, gamma_table, pdf_table)
	real maxg, delgam, gamma_table(:), pdf_table(:), func(pdf_sz)
	integer i

	!initialize gamma_table and pdf_table
	!tables are integrated PDFs for initialization of maxwellians for electrons and ions
	!delgam in the input file sets the temperature in delta gamma, here it is converted 
	!between species by multiplying by mass.
	
	maxg=(delgam)*20+1. !20 times the temperature is usually enough
	
	do i=1,pdf_sz
		gamma_table(i)=(maxg-1.)/(pdf_sz-1)*(i-1)!+1. remove 1. from definition to avoid underflows. 
	enddo
	
#ifndef twoD
! a safer way of writing gamma*sqrt(gamma^2-1)*exp(-(gamma-1)/delgam) -- no underflows
	func=(gamma_table+1.)*sqrt(gamma_table*(gamma_table+2.))*exp(-(gamma_table)/delgam) !3D distribution
#else
!	if(pcosthmult .eq. 0) func=gamma_table*exp(-(gamma_table-1)/delgam)
 	if(pcosthmult .eq. 0) func=(gamma_table+1.)*exp(-(gamma_table)/delgam)
! if magnetized shocks, initialize a full 3d distribution
!	if(pcosthmult .eq. 1) func=gamma_table*sqrt(gamma_table**2-1)*exp(-(gamma_table-1)/delgam) 	
	if(pcosthmult .eq. 1) func=(gamma_table+1.)*sqrt(gamma_table*(gamma_table+2.))*exp(-(gamma_table)/delgam) 	
#endif
	
	pdf_table(1)=0.

	do i=2,pdf_sz
		pdf_table(i)=sum(func(1:i))
	enddo
	
	!normalize pdf_table
	
	pdf_table=pdf_table/pdf_table(pdf_sz)
	

end subroutine init_maxw_table

!-------------------------------------------------------------------------------
! 						subroutine maxwell_dist()					 
!																		
! Maxwell distribution for the particles (similar to above)
!
!-------------------------------------------------------------------------------

subroutine maxwell_dist(gamma0, cd, dseed, u,v,w,gamma_table,pdf_table,pdf_sz)
	
	implicit none
	
	! dummy variables
	
	integer :: pdf_sz
	real(sprec), dimension(:), intent(in) :: gamma_table, pdf_table
	real(sprec) ::  cd, gamma0,  u,v,w
	real(dprec) :: dseed

	! local variables
	
	integer :: i
	logical flag
	real v0, rannum, gam
	real pcosth, pphi, psinth, v0t, ut1, vt1, wt1, ptx, pty, ptz
	real px0, px1, py0, py1, pz0, pz1, gam1, ek 
	
	
	v0=cd*sqrt(1.-1./gamma0**2)
	
	rannum=random(dseed)
	if(rannum .eq. 1.0) rannum=random(dseed)
	
	i=1
	flag=.true.
	gam=0. !gam = gamma-1., only kinetic energy to avoid underflows.
	!choose gamma from the table
	do while (flag)
		
		if(i.eq.pdf_sz) then
			gam=gamma_table(pdf_sz)
			flag=.false.
		endif         
		if(rannum .ge. pdf_table(i) .and. rannum .lt. pdf_table(i+1)) then
			
			gam=gamma_table(i)+(gamma_table(i+1)-gamma_table(i)) &
			/(pdf_table(i+1)-pdf_table(i))*(rannum-pdf_table(i)) 
			
			flag=.false.
		endif
		i=i+1
		
	enddo
	
	!choose phase space angles
#ifdef twoD
		pcosth=(2*random(dseed)-1)*pcosthmult
		if (sigma .ne. 0.) pcosth=2*random(dseed)-1
#else 
		pcosth=2*random(dseed)-1
#endif
	
	pphi=random(dseed)*2*pi
	psinth=sqrt(1-pcosth**2)
	
!	v0t=cd*sqrt((gam-1)*(gam+1))/gam !*sqrt(1.-1./gam**2)
	v0t=cd*sqrt((gam)*(gam+2.))/(1.+gam) !*sqrt(1.-1./gam**2)
	
	ut1=v0t*psinth*cos(pphi)
	vt1=v0t*psinth*sin(pphi)
	wt1=v0t*pcosth
	
	ptx=(1.+gam)*ut1
	pty=(1.+gam)*vt1
	ptz=(1.+gam)*wt1
	
	px1=(ptx+v0*(gam+1.))*gamma0 !this will include the sign of gamma0
 	py1=pty
	pz1=ptz
	
	gam1=sqrt(1+(px1**2+py1**2+pz1**2)/cd**2)
	u=px1/cd                 
	v=py1/cd                  
	w=pz1/cd                  
	
	
end subroutine maxwell_dist



!-------------------------------------------------------------------------------
! 						subroutine powerlaw3d()					 
!																		
! Calculates a power law for the Bell instability (CR) problem 
!
!-------------------------------------------------------------------------------

subroutine powerlaw3d(gammashock, pmin, cd, dseed, u,v,w)
	
	implicit none
	
	! dummy variables

	real(sprec) :: gammashock, pmin, cd, u,v,w
	real(dprec), intent(in) :: dseed

	! local variables
	
	logical flag
	real v0, rannum, gam, potencia, p	
	real pphi, mu, v0t, ut1, vt1, wt1, ptx, pty, ptz
	real px1, py1, pz1 
	
	v0=cd*sqrt(1.-1./gammashock**2)
	
	rannum=random(dseed)
	flag=.true.
	do while (flag) !putting an arbitrary cut off to the momentum of CRs
		flag=.false.
		if(rannum .gt. .9999) then
			rannum=random(dseed)
			flag=.true.
		endif
	enddo
	
	p = pmin!/(1.-rannum)
	gam = sqrt(1.+p**2)
	
	!choose phase space angles
	!      pcosth=2*random(dseed)-1
	pphi=random(dseed)*2*pi
	mu=1.!*random(dseed)! - 1.
	
	v0t=cd*sqrt((gam-1)*(gam+1))/gam !*sqrt(1.-1./gam**2)
	
	ut1=v0t*mu
	vt1=v0t*sqrt(1.-mu**2)*sin(pphi)
	wt1=v0t*sqrt(1.-mu**2)*cos(pphi)
	
	ptx=gam*ut1
	pty=gam*vt1
	ptz=gam*wt1

	px1=(ptx+v0*gam)*gammashock
	py1=pty
	pz1=ptz
	
	u=px1/cd !/gam1
	v=py1/cd !/gam1
	w=pz1/cd !/gam1

end subroutine powerlaw3d 



!-------------------------------------------------------------------------------
! 				subroutine inject_from_wall()	
!															
! Injects plasma from a wall set by global coordinates. Wall is perpendicular to x axis.
! Plasma can be a Maxwellian
! with individual ion and electron temperatures (set as spread in gamma factor,
! the distribution in momentum space is exp(-gamma/delgame), exp(-gamma/delgami) )
! 
!-------------------------------------------------------------------------------

subroutine inject_from_wall(x1,x2,y1,y2,z1,z2,ppc,gamma_drift,delgam_i,delgam_e,&
	wall_speed,weight,upsamp_e_in,upsamp_i_in)
	implicit none
	real(sprec):: x1,x2,xt,y1,y2,z1,z2, x1new,x2new
	real(sprec):: y1new,y2new, z1new, z2new

	real ppc,delgam_i,delgam_e,beta_inj,gamma_drift,wall_speed,beta_wall,weight
	integer direction
	integer, optional :: upsamp_e_in, upsamp_i_in
        integer :: upsamp_e, upsamp_i

	if(.not. present(upsamp_e_in)) then 
	   upsamp_e = 1
	else
	   upsamp_e = upsamp_e_in
	endif
	if(.not. present(upsamp_i_in)) then
	   upsamp_i = 1
	else
	   upsamp_i = upsamp_i_in
	endif

	if(abs(wall_speed) .ge. 1) then 
	   beta_wall=sign(sqrt(1-1/wall_speed**2),wall_speed) !if wall_speed > 1, it is interpreted as gamma of wall
	else
	   beta_wall=wall_speed !otherwise, it's wall speed normalized by c
	endif
	
	if(abs(gamma_drift) .ge. 1) then 
	   beta_inj=sign(sqrt(1-1./gamma_drift**2),gamma_drift)
	else
	   beta_inj=gamma_drift
	endif
	
	if(x1 .ne. x2 .and. y1 .ne. y2) then
!#ifndef twoD
	    if(z1 .ne. z2) then
!#endif
	       print *,rank,": in injection from wall: wall is not along x,y,or z"
	       stop
!#ifndef twoD
	    endif
!#endif
	endif

!determine which wall is the inection wall and set coordinates of a small 
!volume that the plasma emitted from the injector will occupy in one step, 
!also including the possibility of wall motion.
!With magnetic field, this way may be too aggressive, may have to inject really
! on a line. This works fine without external fields. 

	if(x1 .eq. x2) then 
	   x1new=x1
	   x2new=x1+(beta_inj-beta_wall)*c
	   direction=1

	   if(x2new .lt. x1new) then !plasma_region expects x2>x1
	      xt=x1new
	      x1new=x2new
	      x2new=xt
	   endif   
	else
	   x1new=x1
	   x2new=x2
	endif

	if(y1 .eq. y2) then 
	   y1new=y1
	   y2new=y1+(beta_inj-beta_wall)*c
	   direction=2

	   if(y2new .lt. y1new) then 
	      xt=y1new
	      y1new=y2new
	      y2new=xt
	   endif   
	else
	   y1new=y1
	   y2new=y2
	endif

	if(z1 .eq. z2) then 
	   z1new=z1
	   z2new=z1+(beta_inj-beta_wall)*c
	   direction=3

	   if(z2new .lt. z1new) then 
	      xt=z1new
	      z1new=z2new
	      z2new=xt
	   endif   
	else
	   z1new=z1
	   z2new=z2
	endif

	call inject_plasma_region(x1new,x2new,y1new,y2new,z1new,z2new,ppc,gamma_drift,delgam_i,delgam_e,&
           weight,direction,upsamp_e,upsamp_i)

end subroutine inject_from_wall

!-------------------------------------------------------------------------------
! 				subroutine inject_plasma_region()	
!															
! Injects plasma in a region set by global coordinates. Plasma can be a Maxwellian
! with individual ion and electron temperatures (set as spread in gamma factor,
! the distribution in momentum space is exp(-gamma/delgame), exp(-gamma/delgami) )
!
!-------------------------------------------------------------------------------

subroutine inject_plasma_region(x1,x2,y1,y2,z1,z2,ppc,gamma_drift_in, &
          delgam_i,delgam_e,weight,direction,upsamp_e_in,upsamp_i_in)
	implicit none

	real(sprec):: x1,x2,y1,y2,z1,z2 !to accurately calculate differences between x-s, which may be >64k, use doubles
	real ppc,delgam_i,delgam_e,gamma_drift,gamma_drift_in,weight
	logical use_density_profile
	integer direction, ions0, ions1, n1

	real(sprec), dimension(3) :: values

	real dx, dy, dz, len !hack
	
!tables for pre-computed maxwellian distributions. Can add more here
	real, dimension(pdf_sz) :: pdf_table_i, gamma_table_i, pdf_table_e, gamma_table_e

	integer :: i, n, iglob_min, iglob_max, jglob_min, jglob_max
	real :: minz, maxz, delta_z, numps, tmp
	real(sprec) :: minx, maxx, miny, maxy, inject_minx, inject_maxx, inject_miny, inject_maxy, delta_y, delta_x

#ifndef twoD
	integer :: kglob_min, kglob_max
	real :: inject_minz, inject_maxz
#endif

	integer downsmp, subsamp, ionsdwn, ionindx !hack for downsmp
	integer ne, ni
	real corrfact
	integer, optional, intent(in) :: upsamp_e_in, upsamp_i_in
	integer :: upsamp_e,upsamp_i

	real gammaprt, gammaperp, betaxprt, corrweight, beta_drift

	if(.not. present(upsamp_e_in)) then 
	   upsamp_e = 1
	else
	   upsamp_e = upsamp_e_in
	endif
	if(.not. present(upsamp_i_in)) then
	   upsamp_i = 1
	else
	   upsamp_i = upsamp_i_in
	endif

	if(upsamp_e < 1) then 
	   print *,"upsamp_e<1. Reduce ppc. Stopping!"
	   stop
	endif
	
	 if(abs(gamma_drift_in) .lt. 1) then
	    gamma_drift=sign(sqrt(1./(1.-gamma_drift_in**2)),gamma_drift_in) !gamma_drift_in can be used as velocity/c, to simplify non-relativistic drift assignment.
	 else
	    gamma_drift=gamma_drift_in
	 endif
	 

! init_maxw assumes dimensionality of the Maxwellian based on the dimensionality of simulation
! and whether B field is non-zero (then you can have 3D maxwellians even in 2D). 
! choice of dimensionality should be made user assigned in a more transparent way.
! Additional tables can be defined in particles.F90

	call init_maxw_table(delgam_i, gamma_table_i, pdf_table_i)
	call init_maxw_table(delgam_e, gamma_table_e, pdf_table_e)

!	call init_maxw_table(1e-4*mi/me, gamma_table_cold, pdf_table_cold)


	iglob_min=3+mxcum
	iglob_max=(mx-2)+mxcum		
		
	minx=3.
	maxx=3.
	
	if(x2 < x1) then 
	   print *, rank, ": Inject_Plasma_Region. x2 < x1. x1=",x1," x2=",x2 
	   stop
	endif
	
	inject_minx=x1
	inject_maxx=x2
	
	if(inject_minx<iglob_min) minx=3.
	if(inject_maxx<iglob_min) maxx=3.
	if(inject_minx>=iglob_max) minx=mx-2
	if(inject_maxx>=iglob_max) maxx=mx-2
	
	if(inject_minx>=iglob_min .and. inject_minx<iglob_max) then
	   minx=inject_minx-1.*mxcum
	endif
	if(inject_maxx>=iglob_min .and. inject_maxx<iglob_max) then
	   maxx=inject_maxx-1.*mxcum
	endif
	
	delta_x=(maxx-minx)
	
	jglob_min=3+mycum !global extent of the y boundaries on this processor
	jglob_max=(my-2)+mycum
	
	miny=3.
	maxy=3.
	
	if(y2 < y1) then 
	   print *, rank, ": Inject_Plasma_Region. y2 < y1. y1=",y1," y2=",y2 
	   stop
	endif
	
	inject_miny=y1
	inject_maxy=y2
	
	! inject_miny and inject_maxy are global coordinates
	! jglob_min and jglob_max are global coordinates
	! minx and miny are local coordinates

	if(inject_miny<jglob_min) miny=3.
	if(inject_maxy<jglob_min) maxy=3.
	if(inject_miny>=jglob_max) miny=my-2
	if(inject_maxy>=jglob_max) maxy=my-2

	if(inject_miny>=jglob_min .and. inject_miny<jglob_max) then
	   miny=inject_miny-mycum
	endif
	if(inject_maxy>=jglob_min .and. inject_maxy<jglob_max) then
	   maxy=inject_maxy-mycum
	endif

	delta_y=(maxy-miny)

	if(delta_y < 0) print *,"deltay", miny, maxy, inject_miny, inject_maxy

	if (maxy==0 .or. miny==0) then !injection region is outside this CPU's domain
		delta_y=0 
		maxy=0
		miny=0
	endif
	
	 
#ifdef twoD
	maxz=4.
	minz=3.
	delta_z=(maxz-minz)
#endif
	
#ifndef twoD

	kglob_min=3+mzcum !global extent of the y boundaries on this processor
	kglob_max=(mz-2)+mzcum
	
	minz=3.
	maxz=3.
	
	if(z2 < z1) then 
	   print *, rank, ": In inject_Plasma_Region. z2 < z1. z1=",z1," z2=",z2 
	   stop
	endif
	
	inject_minz=z1
	inject_maxz=z2

	if(inject_minz<kglob_min) minz=3.
	if(inject_maxz<kglob_min) maxz=3.
	if(inject_minz>=kglob_max) minz=mz-2.
	if(inject_maxz>=kglob_max) maxz=mz-2.

	if(inject_minz>=kglob_min .and. inject_minz<kglob_max) then
!	   minz=inject_minz-rank/(sizex*sizey)*(mzall-5)
	   minz=inject_minz-mzcum
	endif

	if(inject_maxz>=kglob_min .and. inject_maxz<kglob_max) then
!	   maxz=inject_maxz-rank/(sizex*sizey)*(mzall-5)
	   maxz=inject_maxz-mzcum
	endif

	delta_z=(maxz-minz)

	if (maxz==0 .or. minz==0) then !injection region is outside this CPU's domain
		delta_z=0 
		maxz=0
		minz=0
	endif

#endif

!now ready to inject the plasma region

	n=0
	
	!compute number to be injected from Poisson statistics
	
	numps=(.5*ppc)*delta_z*delta_y*delta_x

	if(numps < 10) then 
	   if(numps .ne. 0.) then 
	      numps=poisson(numps)
	   else
	      numps=0. !avoid unlikely event of poisson returning nonzero number for numps=0
	   endif
	else 
	   numps=ceiling(numps)
	endif

	ions0=ions

	call check_overflow_num(numps)	

	do while(n < int(numps))
		
		n=n+1
		ions=ions+1

		!  Add some random spread:
		
		p(ions)%x=minx+delta_x * random(dseed) 
		p(ions)%y=miny+delta_y * random(dseed) 
		p(ions)%z=minz+delta_z * random(dseed) 
		
		call maxwell_dist(gamma_drift,c,dseed,p(ions)%u,p(ions)%v,p(ions &
		)%w,gamma_table_i, pdf_table_i,pdf_sz ) ! , 0., 1.)	

		p(ions)%ch=weight 
	
		if(direction.eq.2)then
		   tmp=p(ions)%u
		   p(ions)%u=p(ions)%v
		   p(ions)%v=tmp
		endif

		if(direction.eq.3)then
		   tmp=p(ions)%u
		   p(ions)%u=p(ions)%w
		   p(ions)%w=tmp
		endif

		beta_drift=sign(sqrt(1.-1./gamma_drift**2),gamma_drift)

		gammaprt=sqrt(1.+p(ions)%u**2+p(ions)%v**2+p(ions)%w**2)
		gammaperp=1.+p(ions)%v**2+p(ions)%w**2
		betaxprt=p(ions)%u/gammaprt
		corrweight=gammaprt**2*(1.+beta_drift*betaxprt)/ &
		(gammaprt**2+gammaperp*gamma_drift**2-gammaperp)
		p(ions)%ch=p(ions)%ch*corrweight


		totalpartnum =totalpartnum+1
		p(ions)%ind=totalpartnum
		p(ions)%proc=rank
		p(ions)%splitlev=1

		if (use_density_profile) then
			values(1)=(p(ions)%x-3.)/c_omp
			values(2)=((p(ions)%y+modulo(rank,sizey)*(myall-5)) -3.)/c_omp
			values(3)=((p(ions)%z+(rank/sizey)*(mzall-5))-3.)/c_omp
			p(ions)%ch=evalf(1,real(values,8))
		endif


		!  Place electrons in the same locations as ions for zero charge
		!  density (consistent with zero or uniform inital electric fields):



		ne=0
		do while (ne < upsamp_e)
		   ne=ne+1
		   lecs=lecs+1

		   p(maxhlf+lecs)%x=p(ions)%x 
		   p(maxhlf+lecs)%y=p(ions)%y
		   p(maxhlf+lecs)%z=p(ions)%z

		   p(maxhlf+lecs)%ch=weight/upsamp_e

		if (use_density_profile) then
			values(1)=(p(ions)%x-3.)/c_omp
			values(2)=((p(ions)%y+modulo(rank,sizey)*(myall-5)) -3.)/c_omp
			values(3)=((p(ions)%z+(rank/sizey)*(mzall-5))-3.)/c_omp
			p(maxhlf+lecs)%ch=evalf(1,real(values,8))
		endif
	
		   call maxwell_dist(gamma_drift,c,dseed,p(lecs+maxhlf)%u,p(lecs+maxhlf)%v,p(lecs &
		   +maxhlf)%w,gamma_table_e, pdf_table_e,pdf_sz) !, 0., 1.)	

		   if(direction.eq.2)then
		      tmp=p(maxhlf+lecs)%u
		      p(maxhlf+lecs)%u=p(maxhlf+lecs)%v
		      p(maxhlf+lecs)%v=tmp
		   endif

		   if(direction.eq.3)then
		      tmp=p(maxhlf+lecs)%u
		      p(maxhlf+lecs)%u=p(maxhlf+lecs)%w
		      p(maxhlf+lecs)%w=tmp
		   endif

		beta_drift=sign(sqrt(1.-1./gamma_drift**2),gamma_drift)

		gammaprt=sqrt(1.+p(lecs+maxhlf)%u**2+p(lecs+maxhlf)%v**2+p(lecs+maxhlf)%w**2)
		gammaperp=1.+p(lecs+maxhlf)%v**2+p(lecs+maxhlf)%w**2
		betaxprt=p(lecs+maxhlf)%u/gammaprt
		corrweight=gammaprt**2*(1.+beta_drift*betaxprt)/ &
		(gammaprt**2+gammaperp*gamma_drift**2-gammaperp)
		p(maxhlf+lecs)%ch=p(maxhlf+lecs)%ch*corrweight

		   totalpartnum =totalpartnum+1
		   p(maxhlf+lecs)%ind=totalpartnum
		   p(maxhlf+lecs)%proc=rank
		   p(maxhlf+lecs)%splitlev=1

		enddo

	     enddo

		injectedions=injectedions+n
		injectedlecs=injectedlecs+ne

		ions1=ions

!now go through all ions I just created and split them
		if(upsamp_i .ne. 1) print *, "upsamp_i", upsamp_i
		if(upsamp_i > 1) then 
		   n=ions0
		   do while(n <= ions1) ! go through all injected ions
		      n1=n
		      if( p(n1)%ch .gt. .7 ) then !check if the ion is up for splitting; all should be here
			 ni=0
			 do while (ni < upsamp_i)
			    ni=ni+1
			    ions=ions+1
			    call copyprt(p(n1), p(ions))
			    p(ions)%ch=p(ions)%ch/upsamp_i

			    call maxwell_dist(gamma_drift,c,dseed,p(ions)%u,p(ions)%v,p(ions &
			    )%w,gamma_table_i, pdf_table_i,pdf_sz ) ! , 0., 1.)	
			    
			    if(direction.eq.2)then
			       tmp=p(ions)%u
			       p(ions)%u=p(ions)%v
			       p(ions)%v=tmp
			    endif

			    if(direction.eq.3)then
			       tmp=p(ions)%u
			       p(ions)%u=p(ions)%w
			       p(ions)%w=tmp
			    endif

			    beta_drift=sign(sqrt(1.-1./gamma_drift**2),gamma_drift) 

			    gammaprt=sqrt(1.+p(ions)%u**2+p(ions)%v**2+p(ions)%w**2) 
			    gammaperp=1.+p(ions)%v**2+p(ions)%w**2 !this assumes drift in x only
			    betaxprt=p(ions)%u/gammaprt
			    corrweight=gammaprt**2*(1.+beta_drift*betaxprt)/ &
			    (gammaprt**2+gammaperp*gamma_drift**2-gammaperp)
			    p(ions)%ch=p(ions)%ch*corrweight


			    totalpartnum =totalpartnum+1
			    p(ions)%ind=totalpartnum
			    p(ions)%proc=rank
			    p(ions)%splitlev=1

			 enddo
			    
!remove the original high weight ion
			 call copyprt(p(ions),p(n1))
			 ions=ions-1
			 n=n-1
!			 totalpartnum=totalpartnum -1 
		      endif	!if ch > .7
		      
		      n=n+1
		   enddo	!while n < ions
		   
		endif		!if upsamp_i > 1

		
end subroutine inject_plasma_region



#ifdef twoD
end module m_particles
#else
end module m_particles_3d
#endif
