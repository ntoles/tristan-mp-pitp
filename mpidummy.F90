!
! Dummy MPI module
!
! This module defines the MPI functions / procedures used in the code, only 
! to be used when no MPI library is available. 
! 
!

#ifdef twoD 

module m_mpidummy

	use m_globaldata
	
#else

module m_mpidummy_3d

	use m_globaldata_3d
	
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

	interface MPI_SendRecv
		module procedure MPI_SendRecv_int
		module procedure MPI_SendRecv_sprec
	end interface MPI_SendRecv
	
	interface MPI_ISend
		module procedure MPI_iSend_int
		module procedure MPI_iSend_dprec
		module procedure MPI_iSend_sprec
	end interface MPI_ISend

	interface MPI_IRecv
		module procedure MPI_IRecv_int
		module procedure MPI_IRecv_dprec
		module procedure MPI_IRecv_sprec
	end interface MPI_IRecv
	
	interface MPI_Send
		module procedure MPI_Send_int
		module procedure MPI_Send_dprec
		module procedure MPI_Send_sprec
	end interface MPI_Send

	interface MPI_Recv
		module procedure MPI_Recv_int
		module procedure MPI_Recv_dprec
		module procedure MPI_Recv_sprec
	end interface MPI_Recv

	
!-------------------------------------------------------------------------------
!	PUBLIC MODIFIERS
!-------------------------------------------------------------------------------

!-------------------------------------------------------------------------------
!	MODULE PROCEDURES AND FUNCTIONS
!-------------------------------------------------------------------------------

	contains



!-------------------------------------------------------------------------------
! Dummy procedure
!-------------------------------------------------------------------------------

subroutine MPI_Init(ierr)
	
	implicit none
	
	! dummy variables
	
	integer :: ierr
	
end



!-------------------------------------------------------------------------------
! Dummy procedure
!-------------------------------------------------------------------------------

subroutine MPI_Comm_rank(MPI_Comm_world, rank, ierr)

	implicit none
	
	! dummy variables

	integer :: rank, MPI_Comm_world, ierr


	rank=0

end



!-------------------------------------------------------------------------------
! Dummy procedure
!-------------------------------------------------------------------------------

subroutine MPI_Comm_size(MPI_Comm_world, size, ierr)

	implicit none
	
	! dummy variables

	integer :: size, MPI_Comm_world, ierr

	
	size=1

end



!-------------------------------------------------------------------------------
! Dummy procedure
!-------------------------------------------------------------------------------

subroutine MPI_SendRecv_sprec(arrsend, lensend, mpi_size, sendrank &
,sendtag,arrrecv, lenrecv, mpi_size1, recvrank, recvtag, mpi_comm_world, status, ierr)

	implicit none

	! dummy variables
	
	integer :: lensend, sendrank, sendtag, lenrecv, recvrank, recvtag, &
			   mpi_comm_world, ierr, status, mpi_size,mpi_size1
	real(sprec) :: arrsend(lensend), arrrecv(lenrecv)

	
	arrrecv=arrsend

end subroutine MPI_SendRecv_sprec



!-------------------------------------------------------------------------------
! Dummy procedure
!-------------------------------------------------------------------------------

subroutine MPI_SendRecv_int(arrsend, lensend, mpi_size, sendrank &
,sendtag,arrrecv, lenrecv, mpi_size1, recvrank, recvtag, mpi_comm_world, status, ierr)

	implicit none
	
	! dummy variables
	
	integer :: lensend, sendrank, sendtag, lenrecv, recvrank, recvtag, &
			   mpi_comm_world, ierr, status, mpi_size,mpi_size1
	integer :: arrsend(lensend), arrrecv(lenrecv)


	arrrecv=arrsend

end subroutine MPI_SendRecv_int



!-------------------------------------------------------------------------------
! Dummy procedure
!-------------------------------------------------------------------------------

subroutine MPI_iSend_int(arrsend,lensend,mpi_size,sendrank, sendtag, &
comm, req, ierr)

	implicit none
	
	integer, dimension(:), intent(in) :: arrsend
	integer :: lensend, mpi_size,sendrank, sendtag, comm, req, ierr

end subroutine MPI_iSend_int



!-------------------------------------------------------------------------------
! Dummy procedure
!-------------------------------------------------------------------------------

subroutine MPI_iSend_sprec(arrsend,lensend,mpi_size,sendrank, sendtag, &
comm, req, ierr)

	implicit none
	
	real(sprec), dimension(:), intent(in) :: arrsend
	integer :: lensend, mpi_size,sendrank, sendtag, comm, req, ierr

end subroutine MPI_iSend_sprec



!-------------------------------------------------------------------------------
! Dummy procedure
!-------------------------------------------------------------------------------

subroutine MPI_iSend_dprec(arrsend,lensend,mpi_size,sendrank, sendtag, &
comm, req, ierr)

	implicit none
	
	real(dprec), dimension(:), intent(in) :: arrsend
	integer :: lensend, mpi_size,sendrank, sendtag, comm, req, ierr

end subroutine MPI_iSend_dprec



!-------------------------------------------------------------------------------
! Dummy procedure
!-------------------------------------------------------------------------------

subroutine MPI_iRecv_int(arrrecv,lensend,mpi_size,sendrank,sendtag, comm, req, ierr)

	implicit none
	
	! dummy variables
	
	integer, dimension(:), intent(out) :: arrrecv 
	integer :: lensend,mpi_size,sendrank,sendtag, comm, req, ierr


end subroutine MPI_iRecv_int



!-------------------------------------------------------------------------------
! Dummy procedure
!-------------------------------------------------------------------------------

subroutine MPI_iRecv_sprec(arrrecv,lensend,mpi_size,sendrank,sendtag, comm, req, ierr)

	implicit none
	
	! dummy variables
	
	real(sprec), dimension(:), intent(out) :: arrrecv 
	integer :: lensend,mpi_size,sendrank,sendtag, comm, req, ierr


end subroutine MPI_iRecv_sprec



!-------------------------------------------------------------------------------
! Dummy procedure
!-------------------------------------------------------------------------------

subroutine MPI_iRecv_dprec(arrrecv,lensend,mpi_size,sendrank,sendtag, comm, req, ierr)

	implicit none
	
	! dummy variables
	
	real(dprec), dimension(:), intent(out) :: arrrecv 
	integer :: lensend,mpi_size,sendrank,sendtag, comm, req, ierr


end subroutine MPI_iRecv_dprec



!-------------------------------------------------------------------------------
! Dummy function
!-------------------------------------------------------------------------------

real function MPI_wtime()
	
	implicit none

	mpi_wtime=0.

end function mpi_wtime



!-------------------------------------------------------------------------------
! Dummy procedure
!-------------------------------------------------------------------------------

subroutine MPI_Recv_int(arrrecv,lensend,mpi_size,sendrank,sendtag, comm, status, ierr)

	implicit none
	
	! dummy variables
	
	integer, dimension(:), intent(in) :: arrrecv
	integer :: lensend,mpi_size,sendrank,sendtag, comm, status, ierr

end subroutine MPI_Recv_int



!-------------------------------------------------------------------------------
! Dummy procedure
!-------------------------------------------------------------------------------

subroutine MPI_Recv_sprec(arrrecv,lensend,mpi_size,sendrank,sendtag, comm, status, ierr)

	implicit none
	
	! dummy variables
	
	real(sprec), dimension(:), intent(in) :: arrrecv
	integer :: lensend,mpi_size,sendrank,sendtag, comm, status, ierr

end subroutine MPI_Recv_sprec



!-------------------------------------------------------------------------------
! Dummy procedure
!-------------------------------------------------------------------------------

subroutine MPI_Recv_dprec(arrrecv,lensend,mpi_size,sendrank,sendtag, comm, status, ierr)

	implicit none
	
	! dummy variables
	
	real(dprec), dimension(:), intent(in) :: arrrecv
	integer :: lensend,mpi_size,sendrank,sendtag, comm, status, ierr

end subroutine MPI_Recv_dprec



!-------------------------------------------------------------------------------
! Dummy procedure
!-------------------------------------------------------------------------------

subroutine MPI_Send_int(arrsend,lensend,mpi_size,sendrank,sendtag, comm, status, ierr)

	implicit none
	
	! dummy variables
	
	integer, dimension(:), intent(in) :: arrsend
	integer :: lensend,mpi_size,sendrank,sendtag, comm, status, ierr

end subroutine MPI_Send_int



!-------------------------------------------------------------------------------
! Dummy procedure
!-------------------------------------------------------------------------------

subroutine MPI_Send_sprec(arrsend,lensend,mpi_size,sendrank,sendtag, comm, status, ierr)

	implicit none
	
	! dummy variables
	
	real(sprec), dimension(:), intent(in) :: arrsend
	integer :: lensend,mpi_size,sendrank,sendtag, comm, status, ierr

end subroutine MPI_Send_sprec



!-------------------------------------------------------------------------------
! Dummy procedure
!-------------------------------------------------------------------------------

subroutine MPI_Send_dprec(arrsend,lensend,mpi_size,sendrank,sendtag, comm, status, ierr)

	implicit none
	
	! dummy variables
	
	real(dprec), dimension(:), intent(in) :: arrsend
	integer :: lensend,mpi_size,sendrank,sendtag, comm, status, ierr

end subroutine MPI_Send_dprec



!-------------------------------------------------------------------------------
! Dummy procedure
!-------------------------------------------------------------------------------

subroutine MPI_reduce(ia,ib,icount, mpi_size,mpi_sum,inode,mcomm,ierr) 

	implicit none
	
	! dummy variables

	integer :: ia(icount), ib(icount)
	integer :: icount, mpi_size,mpi_sum,inode,mcomm,ierr


	ib=ia

end subroutine MPI_reduce



!-------------------------------------------------------------------------------
! Dummy procedure
!-------------------------------------------------------------------------------

subroutine MPI_ALLREDUCE(curx,curx1,all, mpi_real,mpi_sum,MPI_Comm_world,ierr)

	implicit none
	
	! dummy variables

	integer ::all,i,mpi_real,mpi_sum,MPI_Comm_world,ierr
	real curx(all), curx1(all)


	do i=1,all
		curx1(i)=curx(i)
	enddo

end subroutine MPI_ALLREDUCE



!-------------------------------------------------------------------------------
! Dummy procedure
!-------------------------------------------------------------------------------

subroutine mpi_wait(request, status, ierr)

	implicit none
	
	! dummy variables
	
	integer :: request, status, ierr

end subroutine mpi_wait



!-------------------------------------------------------------------------------
! Dummy procedure
!-------------------------------------------------------------------------------

subroutine mpi_waitall(request, status, ierr)

	implicit none
	
	! dummy variables
	
	integer :: request, status, ierr

end subroutine mpi_waitall



!-------------------------------------------------------------------------------
! Dummy procedure
!-------------------------------------------------------------------------------

subroutine  MPI_barrier(MPI_Comm_world, ierr)

	implicit none
	
	! dummy variables
	
	integer :: (MPI_Comm_world, ierr

end



!-------------------------------------------------------------------------------
! Dummy procedure
!-------------------------------------------------------------------------------

subroutine  MPI_finalize(MPI_Comm_world, ierr)

	implicit none
	
	! dummy variables
	
	integer :: MPI_Comm_world, ierr

end



!-------------------------------------------------------------------------------
! Dummy procedure
!-------------------------------------------------------------------------------

subroutine  MPI_type_commit(dummy, ierr)

	implicit none 
	
	! dummy variables
	
	integer :: dummy, ierr

end



!-------------------------------------------------------------------------------
! Dummy procedure
!-------------------------------------------------------------------------------

subroutine  MPI_type_extent(MPI_Comm_world,dumm,ierr)

	implicit none
	
	! dummy variables
	
	integer :: MPI_Comm_world,dumm,ierr

end



!-------------------------------------------------------------------------------
! Dummy procedure
!-------------------------------------------------------------------------------

subroutine  MPI_type_struct(MPI_Comm_world, ierr)

	implicit none
	
	! dummy variables
	
	integer :: MPI_Comm_world, ierr

end



!-------------------------------------------------------------------------------
! Dummy procedure
!-------------------------------------------------------------------------------

subroutine  MPI_gather_sprec(x1,len,n,x2,len1,rank,mpi_comm_world,ierr)

	implicit none
	
	! dummy variables
	
	real(sprec), dimension(:) :: x2, x1
	integer :: rank,MPI_Comm_world, ierr, x1, len, n, x2, len1

	x2=x1

end subroutine MPI_gather_sprec



!-------------------------------------------------------------------------------
! Dummy procedure
!-------------------------------------------------------------------------------

subroutine  MPI_gather_dprec(x1,len,n,x2,len1,rank,mpi_comm_world,ierr)

	implicit none
	
	! dummy variables
	
	real(dprec), dimension(:) :: x2, x1
	integer :: rank,MPI_Comm_world, ierr, x1, len, n, x2, len1

	x2=x1

end subroutine MPI_gather_dprec



!-------------------------------------------------------------------------------
! Dummy procedure
!-------------------------------------------------------------------------------

subroutine  MPI_gather_int(x1,len,n,x2,len1,rank,mpi_comm_world,ierr)

	implicit none
	
	! dummy variables
	
	integer, dimension(:) :: x2, x1
	integer :: rank,MPI_Comm_world, ierr, x1, len, n, x2, len1

	x2=x1

end subroutine MPI_gather_int


#ifdef twoD
end module m_mpidummy
#else
end module m_mpidummy_3d
#endif