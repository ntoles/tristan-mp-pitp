!
! Main loop file
!
! Contains the main loop, which calls all other functions and procedures as
! the simulation progresses
!
!

#ifdef twoD 

module m_tris

	use m_globaldata
	use m_system
	use m_aux
	use m_communications
	use m_output
	use m_particles
	use m_fields
	use m_fieldboundaries
	use m_domain
	use m_overload
	use m_particles_movedeposit
	use m_dynamic_domain
#else

module m_tris_3d

	use m_globaldata_3d
	use m_system_3d
	use m_aux_3d
	use m_communications_3d
	use m_output_3d
	use m_particles_3d
	use m_fields_3d
	use m_fieldboundaries_3d
	use m_domain_3d
	use m_overload_3d
	use m_particles_movedeposit_3d
	use m_dynamic_domain_3d
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

	public :: mainloop, finalize

!-------------------------------------------------------------------------------
!	MODULE PROCEDURES AND FUNCTIONS
!-------------------------------------------------------------------------------

	contains



!-------------------------------------------------------------------------------
! 						subroutine mainloop					 
!																		
! Performs the main loop of tristan, calling all relevant functions/proceedures
! 							
!-------------------------------------------------------------------------------

subroutine mainloop()

	implicit none
	
	logical, parameter :: tmstop=.true.
	integer :: ierr, n
	real (dprec) :: temp

		integer n1
		real xg,yg,zg
		
	
	do lap=lapst,last !every time step is called "lap"
	
	    call timer(1)
		call timer(2)

		! Advance the Magnetic field half time step	
!		print *, "prebc" !hack

		call pre_bc_b()
!		print *, "advance b half" !hack		
		call advance_Bhalfstep()
				
!		print *, "bc_b1" !hack
		call bc_b1() !periodic bcs, no surface
			
		if(debug) print *,rank,":advanced b", "step=", lap
		
		call timer(2,tmstop=tmstop)
		call timer(3)

		!can calcualte energy diagnostics here -- see energy.f
		
		! Advance particle velocities
		temp=mpi_wtime()

		call move_particles()
		
		call timer(3,tmstop=tmstop)
		call timer(4)

		call timer(14) !hack
		call advance_Bhalfstep()
				
		if(debug) print *,rank,"b4 bc_b2", "step=", lap
		call timer(14,tmstop=tmstop) !hack
		call timer(15) !hack
		call bc_b2() !with radiation bc

		if(debug) print *,rank,"aft bc_b2", "step=", lap
		call timer(15,tmstop=tmstop) !hack
		call timer(16) !hack
		call field_bc_user()
		
		call timer(16,tmstop=tmstop) !hack

		call timer(24) !hack

		call post_bc_b()
		if(debug) print *,rank,": advanced b2", "step=", lap
		call pre_bc_e()
			
		call advance_Efield()
							
		call field_bc_user()
		
		call timer(24,tmstop=tmstop) !hack
		call timer(25) !hack
		call bc_e2() !with radiation bc

		call post_bc_e()

		call field_bc_user()

		if(debug) print *,rank,": advanced e", "step=", lap,ions,lecs
		call timer(25,tmstop=tmstop) !hack

		call timer(4,tmstop=tmstop)

		call timer(5)

		call reset_currents() !set arrays for current to 0.

                if(user_part_bcs) call particle_bc_user()  !e.g., reflecting walls
		call deposit_particles()
		
		call timer(5,tmstop=tmstop)
		
		if(debug) print *,rank,": deposit", "step=", lap

		call timer(6)

		call exchange_particles()
							
		call timer(6,tmstop=tmstop)

		call timer(7)
		
		if(debug) then
			print *,rank,": exchange p", "step=", lap
			call mpi_barrier(mpi_comm_world,ierr)
			print *,rank,"b4 exchange" ,"step=",lap
		endif
		call exchange_current() !include periodic bc on current
		
		if(debug) print *,rank,": exchcur", "step=", lap
		
		call timer(7,tmstop=tmstop)
		call timer(8)

		! ! Uncomment one of the following filters.
		! ! Check domain.F90 and fields.F90 to uncomment array allocations if using a filter other than 1 or 2.

		! Original. Possibly fastest at 1 or 2 cores, or when using a moving wall
		!
#ifdef twoD
                call apply_filter1_opt()
#else
		! Reordered. Usually fastest at more cores.  Only works if all dimensions >= ntimes
#ifdef filter2
		if(ntimes .le. my .or. ntimes .le. mz) then 
		call apply_filter2_opt()  !comment out to hack no filter
!!		call apply_filter1_opt()
		else
                call apply_filter1_opt()
	     endif
#else
                call apply_filter1_opt()
!#ifdef filter2
#endif 
		  
#endif
		! Reordered with copy. Probably never the fastest
		!call apply_filter2_copy_opt()

		! Chunk filtering. Warning: boundary conditions not fully implemented
		!call apply_filter3_opt()

		call timer(8,tmstop=tmstop)
		
		if(debug) print *,rank,": filter", "step=", lap
		
		call timer(9)

		call add_current()
		
		call conduct_bc()

		call bc_e1() 

		call field_bc_user()

		call timer(9,tmstop=tmstop)
		
		if(debug) print *,rank,": injecting", " step=", lap

		call timer(10)
		
		call timer(20)
		call inject_others() !inject particles that enter from other processors
		call timer(20, tmstop=tmstop)

		if(debug) print *,rank,": aft inject_others", " step=", lap

		call timer(21)
!		call mpi_allreduce((LenIonOutDwn+LenIonOutUp+LenLecOutDwn+LenLecOutUp),outcorner, &
!                 1,mpi_integer,mpi_max,mpi_comm_world,ierr)
		  call timer(21,tmstop=tmstop)
		  
		  outcorner=1
		  if(outcorner .ne. 0) then
		   if(debug) print *, rank,": b4 second exchange", " step=", lap
		   call timer(22)	   
			call exchange_particles() !to fix the corner going particles, moving in z 
			call timer(22, tmstop=tmstop)
		
			call timer(23)
		   if(debug) print *, rank,": b4 second inject", " step=", lap
			call inject_others()
			call timer(23,tmstop=tmstop)
		endif
		call timer(24)

		call timer(24,tmstop=tmstop)
		
		call timer(19)
		call inject_particles() 
		call timer(19,tmstop=tmstop)
		if(debug) print *,rank,": aft inject_particles", " step=", lap
		  
		  call timer(25)
		call check_overflow()
		call timer(25,tmstop=tmstop)

		if(debug) print *, "aft inject", "step=",lap, ions, lecs, receivedlecs
		
		call timer(10,tmstop=tmstop)
			
		call Diagnostics()		
		
		call timer(11)
		
		call enlarge_domain()
		
		call timer(11,tmstop=tmstop)
				
		call timer(12)
		
		call reorder_particles()
		   
		call timer(12,tmstop=tmstop)
			
		call timer(31)
				
		call redist_x_domain() 
		
		call timer(31,tmstop=tmstop)
		
		call timer(32)
		
		call redist_y_domain() 
		
		call shift_domain() 		
		
		call timer(32,tmstop=tmstop)

#ifndef twoD		
		call timer(33)
		
		call redist_z_domain() 
		
		call timer(33,tmstop=tmstop)		
#endif		
		
!		call split_parts()
		
		if(debug) print *,rank,": endofloop", "step=", lap

		!if(rank.eq.0) time1=mpi_wtime()
		
		call timer(1,tmstop=tmstop)
		
		call print_timers()
		
		
		call pause_simulation()
		
		call mpi_barrier(MPI_COMM_WORLD,ierr)
		
	enddo

end subroutine mainloop



!-------------------------------------------------------------------------------
! 						subroutine pause_simulation					 
!																		
! Pauses the simulation if a "pause" file exists
! 							
!-------------------------------------------------------------------------------

subroutine pause_simulation()

	implicit none

	logical :: exst2
	
	if(rank .eq. 0) then 
		inquire(file="pause",exist=exst2)
		do while(exst2)
			call sleep(20)
			print *, "rank 0: waiting"
			inquire(file="pause",exist=exst2)
		enddo
	endif


end subroutine pause_simulation



!-------------------------------------------------------------------------------
! 						subroutine energy()					 
!																		
! Calculates the energy in the simulation, for energy conservation purposes
! (should not be in this module, has to be re-writen to fit another module)							
!-------------------------------------------------------------------------------

subroutine energy()

	implicit none

	! local variables
	
	real enarred(3),enarr(3)
	real een,ben,ken,uxtot,uytot,uztot,utot,gamtot,vav,gav
	integer :: elap, ierr
	character fnamen*20
	logical ext
	integer ::statusfile
	integer ::i,j,k,i1,i2
	
	
	! to be checked for electron-ion shocks and with receding injector; also, not sure if using the total energy instead of the kinetic energy is making any difference
	
	elap=5 
	
	if(periodicx.eq.1) then
		i1=3
		i2=mx-3
	else 
		if(wall) then
			i1=leftwall
			i2=min(mx-3,int(xinject2))
		endif
	endif
	
	if(rank.eq.0 .and. lap .lt. lapst+elap) then !first time
		write(fnamen,"(a6)") "energy"
		if(irestart.eq.0) then 
			open(unit=21,file=fnamen)
			close(21,status='delete') 
			open(unit=21,file=fnamen)
		endif
	
		if(irestart.eq.1) then 
			inquire(file=fnamen,exist=ext)
			statusfile=0
			if(ext) statusfile=1
			if(statusfile .eq. 0) then 
				open(unit=21,file=fnamen)
			else
				open(unit=21,file=fnamen,status='old',position='APPEND')
			endif
		endif
	
	endif
	
	
	!      goto 809                  !skip energy determination
	if (modulo(lap,elap) .eq. 0) then
		een=0
		ben=0
#ifdef twoD        
		do k=1,1
#else
		do k=3,mz-3
#endif 
			do j=3,my-3
				do i=i1,i2
					een=een+(ex(i,j,k)**2+ey(i,j,k)**2+ez(i,j,k)**2)
					ben=ben+(bx(i,j,k)**2+by(i,j,k)**2+bz(i,j,k)**2)
				enddo
			enddo
		enddo
		
		een=een
		ben=ben
		
		!     kinetic energy
		ken=sum((sqrt(1.+p(1:ions)%u**2+p(1:ions)%v**2+p(1:ions) &
		%w**2)-1))*c**2*mi+sum((sqrt(1.+p(maxhlf+1:maxhlf+lecs)%u**2 &
		+p(maxhlf+1:maxhlf+lecs)%v**2+p(maxhlf+1:maxhlf+lecs)%w**2 &
		)-1))*c**2*me
				
		enarr(1)=een
		enarr(2)=ben
		enarr(3)=ken
		enarred=0.
		call MPI_allreduce(enarr,enarred,3,mpi_real,mpi_sum &
		,MPI_Comm_World,ierr)
		
		if(rank.eq.0) then 
			een=enarred(1)
			ben=enarred(2)
			ken=enarred(3)
			
#ifdef twoD
				een=.5*enarred(1)/((i2-i1+1)*(my0-5)) 
				ben=.5*enarred(2)/((i2-i1+1)*(my0-5))
				ken=enarred(3)/((i2-i1+1)*(my0-5))
#else
				een=.5*enarred(1)/((i2-i1+1)*(my0-5)*(mz0-5))
				ben=.5*enarred(2)/((i2-i1+1)*(my0-5)*(mz0-5))
				ken=enarred(3)/((i2-i1+1)*(my0-5)*(mz0-5))
#endif
			
			
			write(21,fmt="(i7,' ',4(E15.6,' '))")lap &
			,een,ben,ken,een+ben+ken
			
			print *, "energy_params",lap,i1,i2,my0,mz0
			
		endif                     !if(rank.eq.0)
	endif                     !if (modulo(lap,elap) .eq. 0)
	
	809  continue

end subroutine energy



!-------------------------------------------------------------------------------
! 						subroutine finalize					 
!																		
! Terminates MPI library
! 							
!-------------------------------------------------------------------------------

subroutine finalize()

	implicit none
	
	call MPI_Finalize(MPI_Comm_World)

end subroutine finalize

#ifdef twoD
end module m_tris
#else
end module m_tris_3d
#endif
