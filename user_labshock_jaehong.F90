
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

module m_user

	use m_globaldata
	use m_system
	use m_aux
	use m_communications
	use m_fields
	use m_particles
	use m_inputparser
	use m_fparser
	use m_domain
	
#else

module m_user_3d

	use m_globaldata_3d
	use m_system_3d
	use m_aux_3d
	use m_communications_3d
	use m_fields_3d
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

	real(sprec) :: temperature_ratio, sigma_ext, bz_ext0, bext0, gamma_piston
	character(len=256) :: density_profile
 
	real :: betainj,left_wall_speedup !, beta_splitinj, split_domain
	integer :: betainj_interval, rightclean, rightwall !split_e, split_injector,
	real :: Tp_ratio, Te_ratio, dens_ratio
	integer :: ywidth, zwidth, rightclean_interval, rightclean_start
	logical :: background
	
!-------------------------------------------------------------------------------
!	INTERFACE DECLARATIONS
!-------------------------------------------------------------------------------
	
!-------------------------------------------------------------------------------
!	PUBLIC MODIFIERS
!-------------------------------------------------------------------------------

	public :: init_EMfields_user, init_particle_distribution_user, &
	inject_particles_user, read_input_user, field_bc_user, get_external_fields, &
	particle_bc_user, shift_domain_user


	public :: betainj !, beta_splitinj, split_domain
!	public :: ywidth, zwidth

!-------------------------------------------------------------------------------
!	MODULE PROCEDURES AND FUNCTIONS
!-------------------------------------------------------------------------------

	contains



!-------------------------------------------------------------------------------
! 						subroutine read_input_shock		
!									
! Reads any variables related to (or needed by) this module
! 							
!-------------------------------------------------------------------------------

subroutine read_input_user()

	implicit none
	integer lextflds, luserpartbcs
	integer :: lwall
	integer :: lbackground

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!CHANGE THIS NAME IF YOU ARE CREATING A NEW USER FILE
!This helps to identify which user file is being compiled through Makefile. 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	if(rank.eq.0)print *, "Using user file user_piston_shock.F90"

	call inputpar_getd_def("problem", "temperature_ratio", 1._sprec, Temperature_ratio)
	call inputpar_gets_def("problem", "density_profile", "0", density_profile)
	call inputpar_geti_def("problem", "leftclean", 0, leftclean)

	call inputpar_getd_def("problem","sigma_ext",0._sprec,sigma_ext)

	call inputpar_geti_def("problem","external_fields",0,lextflds)

	if(lextflds .eq. 1) then 
	   external_fields =.true.
	else
	   external_fields =.false.
	endif

!	if(external_fields) bz_ext0 = sqrt((gamma0-1)*.5*ppc0*c**2*(mi+me)*sigma_ext)
	if(external_fields) bext0 = sqrt((gamma0-1)*.5*ppc0*c**2*(mi+me)*sigma_ext)

	call inputpar_geti_def("problem","user_part_bcs",0,luserpartbcs)

	if(luserpartbcs .eq. 1) then 
	   user_part_bcs = .true.
	else
	   user_part_bcs = .false.
	endif

	call inputpar_geti_def("problem", "caseinit", 0, caseinit)

	call inputpar_geti_def("problem", "wall", 0, lwall)
	
	call inputpar_getd_def("problem","wallgam",1._sprec,wallgam)

	if(wallgam .eq. 0.) wallgam=1.	
	if(wallgam<1) wallgam=sign(1./sqrt(1-wallgam**2),wallgam)

        call inputpar_getd_def("problem", "left_wall_speedup",1._sprec,left_wall_speedup) !extra speed multiplier with which left wall will effectively move; wallgam will set the speed for reflecting particles off the wall; this is to eliminate downstream in upstream frame

	call inputpar_getd_def("problem", "betainj", .99999, betainj)
	call inputpar_geti_def("problem", "betainj_interval",1000, betainj_interval)

!	call inputpar_geti_def("problem", "split_injector",0, split_injector)
!	call inputpar_getd_def("problem", "beta_splitinj", 0., beta_splitinj) !move the wall to be splitdomain behind the splitter
!	call inputpar_getd_def("problem", "split_domain", 0., split_domain) !move the wall to be splitdomain behind the splitter
!	call inputpar_geti_def("particles", "split_e", 1, split_e) !move the wall to be splitdomain behind the splitter	

	call inputpar_getd_def("problem", "gamma_piston",0., gamma_piston) 
	call inputpar_geti_def("problem", "rightclean",0, rightclean) 
	call inputpar_geti_def("problem", "rightclean_interval",300, rightclean_interval) 
	call inputpar_geti_def("problem", "rightclean_start",500, rightclean_start) 
	call inputpar_geti_def("problem", "rightwall",0, rightwall) 

	call inputpar_getd_def("problem", "Tp_ratio",1., Tp_ratio) 
	call inputpar_getd_def("problem", "Te_ratio",1., Te_ratio) 
	call inputpar_getd_def("problem", "dens_ratio",1., dens_ratio) 
	call inputpar_geti_def("problem", "background",1, lbackground) 
	
	call inputpar_geti_def("problem", "ywidth",10000, ywidth) 
	call inputpar_geti_def("problem", "zwidth",10000, zwidth) 


	if (lwall==1) then
		wall=.true.
	else
		wall=.false.
	endif
	
	if(wall) user_part_bcs=.true.
	
	
	if(lbackground==1) then
		background=.true.
	else
		background=.false.
	endif

end subroutine read_input_user

!-------------------------------------------------------------
!     Compute external fields to be added to the mover. 
!     These fields do not evolve via Maxwell Eqs, but can depend on time
!-------------------------------------------------------------
	subroutine get_external_fields(x,y,z,ex_ext, ey_ext, ez_ext, bx_ext,by_ext,bz_ext, qm)
	
	real,intent(inout):: bx_ext,by_ext,bz_ext, ex_ext, ey_ext, ez_ext
	real, intent(in):: x,y,z
	real, optional :: qm
	
	ex_ext=0.
	ey_ext=0.
	ez_ext=0.
	bx_ext=0.
	by_ext=0.
	bz_ext=0.


	
	end subroutine get_external_fields


!-------------------------------------------------------------------------------
! 						subroutine parse_density_profile_function()
!												
! Parses the mathematical function that defines the density profile, defined in
! the input file as density_profile
!-------------------------------------------------------------------------------

subroutine parse_density_profile_function(use_density_profile)

	implicit none
	
	! dummy variables
	
	logical, intent(out) :: use_density_profile
	
	! local variables
	
	character(len=1), dimension(3) :: vars=(/'x','y','z'/)
	logical, save :: initialized=.false.

	use_density_profile=.true.	
	
	if (density_profile=="0") then
		use_density_profile=.false.
		return
	endif
	
	if (.not. initialized) then	
		call initf(10)
		call parsef (1, density_profile, vars)
		initialized=.true.
	endif
	
end subroutine parse_density_profile_function



!-------------------------------------------------------------------------------
! 						subroutine init_EMfields_shock		 
!												
! Sets the electromagnetic fields of any specific user purpose
!							
!-------------------------------------------------------------------------------

subroutine init_EMfields_user()
	
	! local variables
	
	integer :: i, j, k, iglob_min, iglob_max, i1,i2
	
	real xmin, xmax, betawall
	
	!determine initial magnetic field based on magnetization sigma which 
        !is magnetic energy density/ kinetic energy density
	!this definition works even for nonrelativistic flows. 
	
	btheta=btheta/180.*pi
	bphi=bphi/180.*pi
	
    binit=sqrt(ppc0*.5*c**2*me*sigma)

	!initialize B field to be set by Binit and the inclination angle -- used for shocks
	do  k=1,mz
		do  j=1,my
			do  i=1,mx
							
				bx(i,j,k)=binit*cos(btheta) 
				by(i,j,k)=binit*sin(btheta)*sin(bphi)
				bz(i,j,k)=binit*sin(btheta)*cos(bphi)

				ex(i,j,k)=0.
				ey(i,j,k)=(-beta)*bz(i,j,k) 
				ez(i,j,k)=-(-beta)*by(i,j,k)

			enddo
		enddo
	enddo
	
	
end subroutine init_EMfields_user


!-------------------------------------------------------------------------------
! 						subroutine init_particle_distribution_shock()	
!											
! Sets the particle distrubtion for a user defined case
!
!-------------------------------------------------------------------------------

subroutine init_particle_distribution_user()

	implicit none

	! local variables
	
	real(sprec), dimension(pdf_sz) :: func
	integer :: i, n, direction, ierr
	real gamma_drift, gamma_drift_in, delgam_i, delgam_e
	real, dimension(pdf_sz) :: pdf_table_i, pdf_table_e, gamma_table_i, gamma_table_e	
	real  ppc, weight
	logical :: use_density_profile
!	integer :: upsamp_e_local, upsamp_i_local
	real :: x1,x2,y1,y2,z1,z2
	character (len=20) :: dummychar
	integer ions0, lecs0, dions, dlecs
	
	call parse_density_profile_function(use_density_profile)
	
	call init_split_parts() !set split points for particle splits, not used unless splitpart is set to true in input
	
	pcosthmult=1 !if 0 and 2D run, the Maxwellian distribution corresponding 
                     !to temperature is initialized in 2D, 1 for 3D.  
	             !when sigma > 0, it is set to 3D automatically
	             		 
	! -------- Set particle boundaries ------------
	!set initial injection points for shocks 	
	leftwall=20.
	xinject=leftwall
	xinject2=mx0-50.
        walloc=leftwall
	! ------------------------------------------
	
	! piston width
!	ywidth=min(5.*c_omp*sqrt(mi/me),my0-5.)	! piston width in units of c/ompi
!	zwidth=min(ywidth,mz0-5.)
		
	totalpartnum=0 !for purpose of keeping track of the total number of particles injected on this cpu

	!jaehongp; 9/2013	
	call read_restart_totpartnum_info()
		
	ions0=0
    lecs0=0
	
	if(background) then
	!**********************************
	!background plasma
	!**********************************
	x1=xinject!+50
	x2=xinject2
	y1=3. 
	y2=my0-2.
	z1=3.
	z2=mz0-2. !if 2D, it will reset automatically

    !**** proton + electron
    ppc= ppc0 
    weight=1.0
    
!    beta = 0. !this will ensure that E fields in the expanding box are not initialized with 
              !the assumption of a drift

    gamma_drift= -gamma0 !0. ! negative gamma_drift will send the plasma in the negative direction
    delgam_i=delgam
    delgam_e=delgam*mi/me*Temperature_ratio

    direction=1 !drift along x
    call inject_plasma_region(x1,x2,y1,y2,z1,z2,ppc,&
    gamma_drift,delgam_i,delgam_e,weight,&
         use_density_profile,direction)     
                       
      
   ! index background plasma as negative signs to distingush from the injected plasma
    p(1:ions)%ind=-p(1:ions)%ind
    p(maxhlf:maxhlf+lecs)%ind=-p(maxhlf:maxhlf+lecs)%ind
    
    endif	!if(background)
    
               
	x1in=3 !set the location of planes where particles are removed from simulation, perpendicular to x. 
	x2in=mx0-2
	
	if(irestart .eq. 1) then 
	   open(unit=20,file="./restart/adjust_param.txt",action="read")
	   read (20,*) dummychar,betainj
!	   read (20,*) dummychar,beta_splitinj
!	   read (20,*) dummychar,split_domain
	endif

	if(rank .eq. 0) then
	    open(unit=10,file="./output/adjust_param.txt",action="write",status="replace")
	    write (10,*) "betainj=",betainj
            write (10,*) "left_wall_speedup=",left_wall_speedup
!	    write (10,*) "beta_splitinj=",beta_splitinj
!	    write (10,*) "split_domain=",split_domain

	    close (10)
	endif
	
	call MPI_BARRIER(MPI_COMM_WORLD,ierr)
	
	call check_overflow()
	call reorder_particles()
	
end subroutine init_particle_distribution_user


!-------------------------------------------------------------------------------
! jaehongp; 9/2013
! This subroutine is designed to prevent newly genarted particles
! when restart begins from being indexed by 1,2,3,4....again.!
!	subroutine read_restart_Ncell_info					 
!-------------------------------------------------------------------------------

subroutine read_restart_totpartnum_info()

	implicit none
	
	! local variables
	
	integer :: dummy
	character (len=70) :: frestartprt
	character (len=10) :: jobname
	logical :: exst

	if (irestart .ne. 1) return
	
	if(size0 .le. 1000) then 
		write(rankchar,"(i3.3)") rank
	else
		write(rankchar,"(i4.4)") rank
	endif
	
	!use jobname from the file created by submit script
	inquire(file="jobname",exist=exst)

	if(exst) then 
		open(unit=20,file="jobname",form='formatted')
		read(20,*) jobname
	else
		jobname="default"
	endif
	
	frestartprt="restart/restprtl."//trim(jobname)//"."// &
	trim(rankchar)//".d"
	
	open(unit=30,file=frestartprt,form='unformatted')	
	rewind(30)
	read(30) dummy,dummy,dummy,dummy,totalpartnum
	close(30)
	
	print*,"totalpartnum :", totalpartnum, "at rank :",rank
	
   	if(mod(totalpartnum,2) .eq. 1) totalpartnum=totalpartnum+1

end subroutine read_restart_totpartnum_info



!-------------------------------------------------------------------------------
! 				subroutine inject_particles_shock()					 
!										
! Injects new particles in the simulation on every step. To skip injection, set ppc=0 below
!
!-------------------------------------------------------------------------------

subroutine inject_particles_user()

	implicit none
	real :: x1,x2,y1,y2,z1,z2
	real delgam_i, delgam_e, injector_speed, gamma_drift, gamma_drift_in
	real ppc, weight
	real, dimension(pdf_sz) :: pdf_table_i, pdf_table_e, gamma_table_i, gamma_table_e	
	logical use_density_profile
	character (len=20) :: dummychar
	integer direction, ierr
	integer ions0, lecs0, dions, dlecs
	real wall_speed
	integer upsamp_e, upsamp_i
	
	real injectvol, num, num1, xglob, yglob, zglob
	real xmin, xmax, ymin, ymax, zmin, zmax

	integer  upsamp_e_local, upsamp_i_local, merge_e, merge_i !merge not working yet
	integer n, n1 !for check
	
!	rightclean=1 --> rightclean in assigned in input file in this new verson
	upsamp_e=1
	upsamp_i=1
	
	use_density_profile=.false. !no profile possible in injection now

	injectedions=0 !set counter to 0 here, so that all possible other injections on this step add up
	injectedlecs=0

	if(rightclean .eq. 1) then
!	   if(betainj .ne. 0) then 
	     xinject2=xinject2+c*.99999 !*0 - beta*c !hack
!	   endif !if betainj was set to 0, we hit the wall, no need to move the injector
	else
	   xinject2=xinject2+c*betainj!*0 -beta*c !hack 	!for when rightclean=0 and the injector moves slowly
	endif

	if(xinject2 .gt. mx0-15) then !stop expansion of the injector has reached the end of the domain
	   xinject2=mx0-15.
	   betainj=0. !injector hit the wall, stop moving the injector
	endif

!	   betainj=0. !injector hit the wall, stop moving the injector

	
	if(rightclean .eq. 1) then 
	   injector_speed = .99999
	else
	   injector_speed = betainj
	endif

 if(lap .gt. 600) then 
!    betainj=0.
!    rightclean=0.
 endif
	ions0=ions
	lecs0=lecs 
	  	 			 		
  	!*******************************
	!******** piston plasma	
	!*******************************           
	x1=xinject
	x2=x1
	y1=max(3.,(my0-5)*.5+3.-.5*ywidth)
	y2=min(my0-2.,(my0-5)*.5+3.+.5*ywidth) 
	z1=max(3.,(mz0-5)*.5+3.-.5*zwidth)
	z2=min(mz0-2.,(mz0-5)*.5+3.+.5*zwidth)

          
    !**** proton  + electron	
    ppc=ppc0*dens_ratio
    weight=1.0
	   	
	gamma_drift=gamma_piston
	wall_speed=0. !speed of the injecting wall
	delgam_i=delgam*Tp_ratio
	delgam_e=delgam*mi/me*Te_ratio
		
	call inject_from_wall(x1,x2,y1,y2,z1,z2,ppc,gamma_drift,&
	  delgam_i,delgam_e,wall_speed,weight, &
         use_density_profile,upsamp_e,upsamp_i)
           
   
	if(background) then	
		
    ions0=ions
    lecs0=lecs	             
          
  	!*******************************
	!******** background plasma	
	!*******************************           
	x1=xinject2 
	x2=x1 !x2=x1 if the injection is parallel to x

	y1=3. !global coordinates
	y2=my0-2.  
	z1= 3.
	z2= mz0-2. !if 2D, it will reset automatically	
        
    !**** proton  + electron	
    ppc=ppc0
    weight=1.0	
    
!    if(lap .gt. 600) ppc =0 !hack
 
    gamma_drift=-gamma0 !0. ! negative gamma_drift will send the plasma in the negative direction
    delgam_i=delgam
    delgam_e=delgam*mi/me*Temperature_ratio

    call inject_from_wall(x1,x2,y1,y2,z1,z2,ppc,gamma_drift,delgam_i, &
       	  delgam_e,injector_speed,weight, &
       	  use_density_profile,upsamp_e,upsamp_i)
      	  
    p(ions0+1:ions)%ind=-p(ions0+1:ions)%ind
    p(maxhlf+lecs0+1:maxhlf+lecs)%ind=-p(maxhlf+lecs0+1:maxhlf+lecs)%ind
   
    endif !if(background)
    
    x1in=3.
    if(rightclean .eq. 1) then
  		x2in=mx0-2.	              
    else
    	x2in=xinject2  !for cases when there is no right wall and we want CRs to leave through the right injector

    	if(rightwall >0) x2in=xinject2+10 !not sure about +10, but setting to xinject2 was killing subsonic electrons, 
    					!not allowing them to scatter off the right wall. 
    endif
  		    
	! *** rightclean
	if(rightclean .eq. 1 .and. lap .gt. rightclean_start .and. modulo(lap, rightclean_interval) .eq. 0) then
		xinject2=1.*int(xinject2-rightclean_interval*c*(.99999 - betainj))
		x2in=xinject2
!remove particles beyond x2in
  if(1>0) then 
	n=1
	do while (n <= ions)
 	   if(p(n)%x + mxcum > x2in) then 
	      call copyprt(p(ions),p(n))
	      ions=ions-1
	      n=n-1
	   endif
	   n=n+1
	enddo

	n=1
	do while (n <= lecs)
	   if(p(n+maxhlf)%x + mxcum > x2in) then 
	      call copyprt(p(maxhlf+lecs),p(maxhlf+n))
	      lecs=lecs-1
	      n=n-1
	   endif
	   n=n+1
	enddo
endif !if 1<0

	endif	

	if(modulo(lap,betainj_interval) .eq. 0) then
      
	    open(unit=20,file="./output/adjust_param.txt",action="read")
	    read (20,*) dummychar,betainj
            read (20,*) dummychar,left_wall_speedup
!	    read (20,*) dummychar,beta_splitinj
!	    read (20,*) dummychar,split_domain

	    if(rank .eq. 0) print *, "read parameters from file", betainj !, beta_splitinj, split_domain, rank

	    close (20)	
	    call MPI_BARRIER(MPI_COMM_WORLD,ierr)	
	    
	    if(rank .eq. 0) then !save into restart directory
	       open(unit=10,file="./restart/adjust_param.txt",action="write",status="replace")
	       write (10,*) "betainj=",betainj
               write (10,*) "left_wall_speedup=",left_wall_speedup
!	       write (10,*) "beta_splitinj=",beta_splitinj
!	       write (10,*) "split_domain=",split_domain
	       
	       close (10)
	    endif
	 endif
	
!	  call MPI_BARRIER(MPI_COMM_WORLD,ierr)


end subroutine inject_particles_user



!-------------------------------------------------------------------------------
! 				subroutine field_bc_shock()	
!										
! Applies boundary conditions specific to user problem. 
! 
!-------------------------------------------------------------------------------

	subroutine field_bc_user()
	implicit none
	integer i,j,k, j1, j2, jglob_min, jglob_max, jglob
	integer i1, i2, iglob_min, iglob_max, iglob
	integer k1, k2, kglob_min, kglob_max, kglob	
	integer ierr
	real xmin, xmax, ymin, ymax, zmin, zmax

!reset fields on the right end of the grid, where the plasma is injected
!make it drifting fields, even though there is no plasma there

	ymin=(my0-5)*.5+3.-ywidth*.5-5
	ymax=(my0-5)*.5+3.+ywidth*.5+5 
	zmin=(mz0-5)*.5+3.-zwidth*.5-5
	zmax=(mz0-5)*.5+3.+zwidth*.5+5
	
	jglob_min=3+mycum
	jglob_max=(my-2)+mycum
	if(ymin<jglob_min) j1=1
	if(ymax<jglob_min) j2=1
	if(ymin>=jglob_max) j1=my
	if(ymax>=jglob_max) j2=my
	if(ymin>=jglob_min .and. ymin<jglob_max) then
		j1=int(ymin)-mycum
	endif
	if(ymax>=jglob_min .and. ymax<jglob_max) then
		j2=int(ymax)-mycum
	endif

#ifndef twoD	
	kglob_min=3+mzcum
	kglob_max=(mz-2)+mzcum
	if(zmin<kglob_min) k1=1
	if(zmax<kglob_min) k2=1
	if(zmin>=kglob_max) k1=mz
	if(zmax>=kglob_max) k2=mz
	if(zmin>=kglob_min .and. zmin<kglob_max) then
		k1=int(zmin)-mzcum
	endif
	if(zmax>=kglob_min .and. zmax<kglob_max) then
		k2=int(zmax)-mzcum
	endif	
#endif	

	xmin=1
	xmax=xinject-5.
	iglob_min=3+mxcum
	iglob_max=(mx-2)+mxcum
	if(xmin<iglob_min) i1=1
	if(xmax<iglob_min) i2=1
	if(xmin>=iglob_max) i1=mx
	if(xmax>=iglob_max) i2=mx
	if(xmin>=iglob_min .and. xmin<iglob_max) then
		i1=int(xmin)-mxcum
	endif
	if(xmax>=iglob_min .and. xmax<iglob_max) then
		i2=int(xmax)-mxcum
	endif	
	
		! unmagnetized piston
#ifdef twoD
		if((i1 .ne. i2) .and. (j1 .ne. j2)) then
			ey(i1:i2,j1:j2,:)=0.
			ez(i1:i2,j1:j2,:)=0.	
		endif
#else
		if((i1 .ne. i2) .and. (j1 .ne. j2) .and. (k1 .ne. k2)) then
			ey(i1:i2,j1:j2,k1:k2)=0.
			ez(i1:i2,j1:j2,k1:k2)=0.							
		endif	
#endif

	xmin=1.
	xmax=5.
	iglob_min=3+mxcum
	iglob_max=(mx-2)+mxcum
	if(xmin<iglob_min) i1=1
	if(xmax<iglob_min) i2=1
	if(xmin>=iglob_max) i1=mx
	if(xmax>=iglob_max) i2=mx
	if(xmin>=iglob_min .and. xmin<iglob_max) then
		i1=int(xmin)-mxcum
	endif
	if(xmax>=iglob_min .and. xmax<iglob_max) then
		i2=int(xmax)-mxcum
	endif	
		
	if(i1 .ne. i2) then
!		ex(i1:i2,:,:)=0.
		ey(i1:i2,:,:)=0.
		ez(i1:i2,:,:)=0.	
	endif
	
	xmin=mx0-10.
	xmax=mx0
	iglob_min=3+mxcum
	iglob_max=(mx-2)+mxcum
	if(xmin<iglob_min) i1=1
	if(xmax<iglob_min) i2=1
	if(xmin>=iglob_max) i1=mx
	if(xmax>=iglob_max) i2=mx
	if(xmin>=iglob_min .and. xmin<iglob_max) then
		i1=int(xmin)-mxcum
	endif
	if(xmax>=iglob_min .and. xmax<iglob_max) then
		i2=int(xmax)-mxcum
	endif	
		
	if(i1 .ne. i2) then
		bx(i1:i2,:,:)=binit*cos(btheta) 
		by(i1:i2,:,:)=binit*sin(btheta)*sin(bphi)
		bz(i1:i2,:,:)=binit*sin(btheta)*cos(bphi)
                ex(i1:i2,:,:)=0.				
		ey(i1:i2,:,:)=(-beta)*bz(i1:i2,:,:) 
		ez(i1:i2,:,:)=-(-beta)*by(i1:i2,:,:)    

	endif

#ifdef null
	xmin=xinject2-1000.
	xmax=xinject2+100.
	iglob_min=3+mxcum
	iglob_max=(mx-2)+mxcum
	if(xmin<iglob_min) i1=1
	if(xmax<iglob_min) i2=1
	if(xmin>=iglob_max) i1=mx
	if(xmax>=iglob_max) i2=mx
	if(xmin>=iglob_min .and. xmin<iglob_max) then
		i1=int(xmin)-mxcum
	endif
	if(xmax>=iglob_min .and. xmax<iglob_max) then
		i2=int(xmax)-mxcum
	endif	
		
	if(i1 .ne. i2) then
               ! ex(i1:i2,:,:)=0
                do i=i1,i2
!                 ex(i,:,:)=ex(i,:,:)/(1.+(i+mxcum-xmin)/(1000.))  
!                   ex(i,:,:)=0.
                enddo
	endif
#endif


!	rightclean =1
			
	! rightclean
	if(rightclean .eq. 1 .and. lap .gt. rightclean_start .and. modulo(lap, rightclean_interval) .eq. 1 ) then !hack
			! global variables
			xmin=xinject2+1 !xinject2+1 !+100 !+100 !hack
			xmax=mx0-3
			iglob_min=3+mxcum
			iglob_max=(mx-2)+mxcum
			if(xmin<iglob_min) i1=3
			if(xmax<iglob_min) i2=3
			if(xmin>=iglob_max) i1=mx-2
			if(xmax>=iglob_max) i2=mx-2
			if(xmin>=iglob_min .and. xmin<iglob_max) then
				i1=xmin-mxcum
			endif
			if(xmax>=iglob_min .and. xmax<iglob_max) then
				i2=xmax-mxcum
			endif			
			if(i1 .ne. i2) then
				bx(i1:i2,:,:)=binit*cos(btheta) 
				by(i1:i2,:,:)=binit*sin(btheta)*sin(bphi)
				bz(i1:i2,:,:)=binit*sin(btheta)*cos(bphi)
				ex(i1:i2,:,:)=0.                         
				ey(i1:i2,:,:)=(-beta)*bz(i1:i2,:,:) 
				ez(i1:i2,:,:)=-(-beta)*by(i1:i2,:,:)    
			endif
	endif	! rightclean
			
	

!	call MPI_BARRIER(MPI_COMM_WORLD,ierr)


	end subroutine field_bc_user
	

!-------------------------------------------------------------------------------
! 				subroutine shift_domain_user
!										
! shift fields and particles backward
!-------------------------------------------------------------------------------
subroutine shift_domain_user
	
	implicit none

	integer :: i, n1, n, i1, i2 ! dummy 
	
	integer :: iL, iR	! the first and last ranks in a row
	
	logical :: in
	
	integer :: error , ierr
	
	real, allocatable, dimension(:,:,:) :: tempbx, tempby, tempbz, tempex, tempey,tempez

	integer mxold, mxave, shiftlength
	
	integer buffsize1
 
    shiftlength=walloc-leftwall	! length shifted backward
    if(mxl(1)-shiftlength .lt. leftwall+3 .or. shiftlength .eq. 0) return

	! This condition can be changed.
 if(walloc .gt. 4*leftwall .or. walloc .gt. mxl(1)-10. ) then
    
    if(rank.eq.0)print*, 'domain shift starts', walloc, 4*leftwall, mxl(1)-10.
   
    
    if(modulo(rank,sizex) .eq. 0) then
       
       if(rank.eq.0)print*, 'old walloc=',walloc, 'new walloc=',walloc-shiftlength,"shiftlength", shiftlength
					
       mxold=mx
       mx=mxold-shiftlength
       ix=1
       iy=mx
       iz=iy*my
       lot=mx*my
       
#ifdef twoD
       iz=0
       lot=mx*my
#endif	
       
       i1=1
       i2=i1+shiftlength
       
       if(rank.eq.0)print*, 'mxold, mx=',mxold, mx
     		
       allocate(tempbx(mx,my,mz),tempby(mx,my,mz),tempbz(mx,my,mz))
       allocate(tempex(mx,my,mz),tempey(mx,my,mz),tempez(mx,my,mz))
       
       tempbx=0. !binit*cos(btheta)
       tempby=0. !binit*sin(btheta)*sin(bphi)
       tempbz=0. !binit*sin(btheta)*cos(bphi)				
       tempex=0. !-(-beta)*tempbz				
       tempey=0.
       tempez=0. !-(beta)*tempbx
       
       tempbx(i1:mx,:,:)=bx(i2:mxold,:,:)
       tempby(i1:mx,:,:)=by(i2:mxold,:,:)
       tempbz(i1:mx,:,:)=bz(i2:mxold,:,:)
       tempex(i1:mx,:,:)=ex(i2:mxold,:,:)
       tempey(i1:mx,:,:)=ey(i2:mxold,:,:)
       tempez(i1:mx,:,:)=ez(i2:mxold,:,:)
   
       i1=int(walloc-shiftlength)
       tempex(1:i1-5,:,:)=0.
       tempey(1:i1-5,:,:)=0.
       tempez(1:i1-5,:,:)=0.
       
			
       p(1:ions)%x=p(1:ions)%x-1.*shiftlength
       p(maxhlf+1:maxhlf+lecs)%x=p(maxhlf+1:maxhlf+lecs)%x-1.*shiftlength

		
       deallocate(bx,by,bz,ex,ey,ez,curx,cury,curz)
       deallocate(bufferin1,bufferin2)
       deallocate(bufferin1y,bufferin2y)
       deallocate(sendbufy,sendbufz)
       deallocate(temp)
#ifdef filter2
       deallocate(yghost, zghost) !for filter2
#endif
       deallocate(poutminus, poutplus,pinminus,pinplus)
       deallocate(poutlft, poutrgt,pinlft,pinrgt)
       deallocate(pall)
			
		
       buffsize1 = max(3*int(ppc0*c*max((mx-5),(my-5))),60000)
       
       allocate(bx(mx,my,mz),by(mx,my,mz),bz(mx,my,mz))
       allocate(ex(mx,my,mz),ey(mx,my,mz),ez(mx,my,mz))
       allocate(curx(mx,my,mz),cury(mx,my,mz),curz(mx,my,mz))
       allocate(bufferin1(mx,my,2),bufferin2(mx,my,2))		
       allocate(bufferin1y(mx,2,mz),bufferin2y(mx,2,mz))
       allocate(sendbufy(mx,2,mz),sendbufz(mx,my,2))
       allocate(temp(mx,my,mz))
#ifdef filter2
       allocate(yghost(mx,2*ntimes,mz))
       allocate(zghost(mx,my,2*ntimes)) !for filter2	
#endif
       allocate(poutminus(buffsize1),poutplus(buffsize1))		
       allocate(pinminus(buffsize1),pinplus(buffsize1))
       allocate(poutrgt(buffsize1),poutlft(buffsize1))		
       allocate(pinrgt(buffsize1),pinlft(buffsize1))
       allocate(pall(lot))
       
       pall=0
       curx=0.
       cury=0.
       curz=0.

       bufferin1=0.
       bufferin2=0.
       bufferin1y=0.
       bufferin2y=0.
       sendbufy=0.
       sendbufz=0.		
       temp=0.
       
       bx=tempbx
       by=tempby
       bz=tempbz
       ex=tempex
       ey=tempey
       ez=tempez
       
       deallocate(tempbx,tempby,tempbz,tempex,tempey,tempez)
	
    endif !if(modulo(rank,sizex) .eq. 0)
     	
    
    call MPI_BARRIER(MPI_COMM_WORLD,ierr)	
     	
    ! shift global variables
    walloc=walloc-1.*shiftlength
    xinject2=xinject2-1.*shiftlength
    xinject3=xinject3-1.*shiftlength	
    x2in=x2in-1.*shiftlength

    call mpi_allgather(mx,1,mpi_integer,mxl,1,mpi_integer &
         ,mpi_comm_world, error)
			
    iL=1
    iR=modulo(rank,sizex)+1
     	
    mxcum=sum(mxl(iL:iR)-5)-(mxl(iR)-5)
    mx0=sum(mxl(1:sizex)-5)+5		
		
    call MPI_BARRIER(MPI_COMM_WORLD,ierr)
    if(rank .eq. 0) print *,"domain shift done  "
 endif
end subroutine shift_domain_user

!------------------------------------------------------------------------------

	subroutine particle_bc_user()
	implicit none
	real invgam, gammawall, betawall, gamma, z0
	real(sprec) :: uout, vout, wout
	real zcolis,q0,tfrac
	real :: x0, xcolis, y0, ycolis, walloc0, wallocrgt
	integer n1,i0,i1, iter, i00
	logical in
	integer reflect, thermal
	logical rightcleanoff
	integer ierr

	
	if(wall) then
				
	   gammawall=wallgam
	   betawall=sqrt(1.-1/gammawall**2)
		
	   walloc=walloc+betawall*c 
	   if(lap .gt. 20000 .and. modulo(lap, 200) .eq. 0 .and. betawall .ne. 0.) walloc=walloc +200*((left_wall_speedup-1)*betawall*c) ! + .014*c !.0125*c

!	   walloc=leftwall
		
!	   if(lap .gt. 20000 .and. split_injector .eq. 1 .and. split_domain > 0 ) then 
!	      walloc = max(walloc, real(xinject3 - split_domain,8)) 
!	      if(walloc .gt. mxl(1)-9) walloc = mxl(1)-9
!	   endif

	   if(movwin) walloc=walloc-movwinoffset

	   do iter=1,2
	      if(iter.eq.1) then 
		 i0=1
		 i1=ions
		 q0=qi
		 i00=0
	      else
		 i0=maxhlf+1
		 i1=maxhlf+lecs
		 q0=qe
		 i00=maxhlf
	      endif
		 
!	   do n1=i0,i1	
	      n1=i0
	     do while(n1 <= i1) 
	    
	   if(p(n1)%x+mxcum .lt. walloc) then 
	      gamma=sqrt(1+(p(n1)%u**2+p(n1)%v**2+p(n1)%w**2))
	      
	      !this algorithm ignores change in y and z coordinates
	      !during the scattering. Including it can result in rare
	      !conditions where a particle gets stuck in the ghost zones. 
	      !This can be improved. 

	      !unwind x location of particle
	      x0=p(n1)%x-p(n1)%u/gamma*c
	      y0=p(n1)%y !-p(n1)%v/gamma*c
	      z0=p(n1)%z !-p(n1)%w/gamma*c

	      !unwind wall location
	      walloc0=walloc-betawall*c-mxcum
	      
	      !where did they meet?
	      tfrac=abs((x0-walloc0)/(betawall*c-p(n1)%u/gamma*c))

	      if(tfrac .lt. 1) then 
		 xcolis=x0+p(n1)%u/gamma*c*tfrac
		 ycolis=y0	!+p(n1)%v/gamma*c*tfrac
		 zcolis=z0	!+p(n1)%w/gamma*c*tfrac

	      !deposit current upto intersection
		 q=p(n1)%ch*q0 !real(splitratio)**(1.-real(p(n1)%splitlev))*q0
		 call zigzag(xcolis,ycolis,zcolis,x0,y0,z0,in)

	      !reset particle momentum, getting a kick from the wall
		 p(n1)%u=gammawall**2*gamma*(2*betawall - p(n1)%u/gamma*(1 + betawall**2))
		 gamma=sqrt(1+(p(n1)%u**2+p(n1)%v**2+p(n1)%w**2))
	      
		 tfrac=min(abs((p(n1)%x-xcolis)/max(abs(p(n1)%x-x0),1e-9)),1.)
              !move particle from the wall position with the new velocity
 
		 p(n1)%x = xcolis + p(n1)%u/gamma*c * tfrac
		 p(n1)%y = ycolis !+ p(n1)%v/gamma*c * tfrac
		 p(n1)%z = zcolis !+ p(n1)%w/gamma*c * tfrac
	     
!now clean up the piece of trajectory behind the wall, that deposit_particles will be adding when it 
!unwinds the position of the particle by the full timestep. 
	      
		 q=-q
		 call zigzag(xcolis,ycolis,zcolis,p(n1)%x-p(n1)%u/gamma*c, & 
		 p(n1)%y-p(n1)%v/gamma*c,p(n1)%z-p(n1)%w/gamma*c,in)
		      
	      else		!if tfrac < 1? 
!       remove the particle
		 call copyprt(p(i1),p(n1))
		 n1=n1-1
		 i1=i1-1
	      endif		!if tfrac < 1
	   endif
	     n1=n1+1
	enddo			! while n1 (particles)

	if(iter .eq. 1)  ions=i1-i00 
	if(iter .eq. 2)  lecs=i1-i00
	enddo			! Do iter  (species)

	endif	! if(wall)
    	

 !=============================================================
!add a second wall at the injector
!	
! algorithm: right wall is moving with background flow starting from xinject2 on every step
! in the process it reflects whatever crossed it, imparting the momentum from reflection 
! off moving wall. 

	if(wall .and. rightwall > 0 ) then 

	   gammawall=-gamma0 !1.
	   betawall= -sqrt(1.-1/gammawall**2) !0.  !minus sign to denote left-moving nature of the wall
	    

!	   wallocrgt=xinject2  +betawall*c !hack  !+ betawall*c*lap 
	   wallocrgt=xinject2  !+ betainj*c !+betainj*c
	   if(movwin) wallocrgt=wallocrgt-movwinoffset


	   do iter=1,2
	      if(iter.eq.1) then 
		 i0=1
		 i1=ions
		 q0=qi
	      else
		 i0=maxhlf+1
		 i1=maxhlf+lecs
		 q0=qe
	      endif
		 
	   do n1=i0,i1	

	   if(p(n1)%x+mxcum .gt. wallocrgt .and. p(n1)%x+mxcum .lt. x2in) then ! .and. p(n1)%ind .gt. 0*totalpartnum - 1e4) then 
	     ! print *, "rank=", rank, n1, p(n1)%y+mycum, walloc, p(n1)%y, walloc-mycum

!	      if( p(n1)%x - (wallocrgt-mxcum) .gt. 1) then 
!		 if( xinject2-mxcum .gt. 3 .and. xinject2-mxcum .lt. mx-2) print *, "xinject2 is on rank", rank
!		 print *, "distance larger than 1 cell behind right wall", p(n1)%x+mxcum, "xinject2=", xinject2, wallocrgt-mxcum, rank, mxcum, n1, mx, x2in
!	      endif

	      gamma=sqrt(1+(p(n1)%u**2+p(n1)%v**2+p(n1)%w**2))
	     
	      !this algorithm ignores change in y and z coordinates
	      !during the scattering. Including it can result in rare
	      !conditions where a particle gets stuck in the ghost zones. 
	      !This can be improved. 

	      !unwind y location of particle
	      x0=p(n1)%x-p(n1)%u/gamma*c
	      y0=p(n1)%y !-p(n1)%v/gamma*c
	      z0=p(n1)%z !-p(n1)%w/gamma*c

	      !unwind wall location
	      walloc0=wallocrgt-betawall*c-mxcum
	      
	      !where did they meet?
	      tfrac=abs((x0-walloc0)/max( abs(betawall*c-p(n1)%u/gamma*c),1e-9))

!	      tfrac=min(abs((x0-walloc0)/max(abs(betawall*c-p(n1)%u/gamma*c),1e-9)),1.)

	      xcolis=x0+p(n1)%u/gamma*c*tfrac !+mycum
	      ycolis=y0 !+p(n1)%v/gamma*c*tfrac
	      zcolis=z0 !+p(n1)%w/gamma*c*tfrac

	      !deposit current upto intersection
	      q=p(n1)%ch*q0 !real(splitratio)**(1.-real(p(n1)%splitlev))*q0
	      call zigzag(xcolis,ycolis,zcolis,x0,y0,z0,in)

	      !reset particle momentum, getting a kick from the wall
	      p(n1)%u=gammawall**2*gamma*(2*betawall - p(n1)%u/gamma*(1 + betawall**2))
	      gamma=sqrt(1+(p(n1)%u**2+p(n1)%v**2+p(n1)%w**2))
	      
	      tfrac=min(abs((p(n1)%x-xcolis)/max(abs(p(n1)%x-x0),1e-9)),1.)
              !move particle from the wall position with the new velocity
 
	      p(n1)%x = xcolis + p(n1)%u/gamma*c * tfrac
	      p(n1)%y = ycolis !+ p(n1)%v/gamma*c * tfrac
	      p(n1)%z = zcolis !+ p(n1)%w/gamma*c * tfrac
	     
!now clean up the piece of trajectory behind the wall, that deposit_particles will be adding when it 
!unwinds the position of the particle by the full timestep. 
	      
	      q=-q
	      call zigzag(xcolis,ycolis,zcolis,p(n1)%x-p(n1)%u/gamma*c, & 
	              p(n1)%y-p(n1)%v/gamma*c,p(n1)%z-p(n1)%w/gamma*c,in)
	      
	     endif
	   
	   
	enddo	! Do n1 (particles)
	enddo	! Do iter  (species)

	endif !if wall
	       		
        		
	
	end subroutine particle_bc_user


#ifdef twoD
end module m_user
#else
end module m_user_3d
#endif

