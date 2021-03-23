MODULE COO_routines
	
	! Subroutines included
	!	matvectcoo - matrix-vector multiplication in COO storage format
	!	MPIvars_init - Initiatialization required for matvectcoo
	
	USE kinds ! describes precision
	USE mp_module ! MPI routines
	
	IMPLICIT NONE
	
	REAL(KIND = RP) :: my_Total_mul_time = 0.D0, my_Total_comm_mkl = 0.D0
	INTEGER(KIND = IP) :: my_range_n, my_range_t
	INTEGER(KIND = IP), ALLOCATABLE :: 	scounts_n(:), scounts_t(:), &
										sdispls_n(:), sdispls_t(:), &
										recvcounts_n(:), recvdispls_n(:), &
										recvcounts_t(:), recvdispls_t(:)
	LOGICAL :: calc_scounts_n, calc_scounts_t
	 
	CONTAINS
	
	SUBROUTINE MPIvars_init( )
	
		IMPLICIT NONE
	
		ALLOCATE(scounts_n(world_size), sdispls_n(world_size), &
				 recvcounts_n(world_size), recvdispls_n(world_size))
		ALLOCATE(scounts_t(world_size), sdispls_t(world_size), &
				 recvcounts_t(world_size), recvdispls_t(world_size))
	
	END SUBROUTINE MPIvars_init
	
	SUBROUTINE matvectcoo(Amat, RowInd, ColInd, nnz, x, y, nrows, &
							everyones_rows, transa)

!		---------------------------------------------------
!		Subroutine discription
!		---------------------------------------------------
!		Matrix*vector in COO format
!		A * x = y


!		The subroutine is designed to have optimal storage using MPI
!		x passed to this routine (and the one stored on each processor) is 
! 		my_x = x(index_start:(index_start+index_range-1))
!		where,
!				index_start = minval(ColInd)
!				index_range = maxval(ColInd) - minval(ColInd) + 1 

! 		However, if you wish to pass the entire x vector
!				assign 
!				index_start = 1
!				index_range = size(x)

!		----------------------------------------------------
!		INPUT DATA DISCRIPTION:
!		----------------------------------------------------

!		If data is stored on each processor
!				nrows - number of rows stored on my_id
!				(Amat(nnz), RowInd(nnz), ColInd(nnz)) - sparse-COO disption of A matrix
!				x(index_range) - input vector
!				nnz(1) - number of non-zeros in Amat (length of Amat)
!				index_start(1) - minimum index for x vector
!				index_range(1) - length of x 

!				everyones_rows(world_size) - number of rows stored on every processor.

!		---------------------------------------------------
!		OUTPUT DATA DISCRIPTION:
!		---------------------------------------------------
!				y(:) = A(minval(RowInd):maxval(RowInd), :)*x(:)

		
		IMPLICIT NONE
		INCLUDE 'mpif.h'
		INCLUDE 'mkl.fi'
		INTEGER(KIND = IP) :: 	nrows
		INTEGER(KIND = IP) :: 	scounts(world_size), sdispls(world_size), &
								recvcounts(world_size), recvdispls(world_size)
		INTEGER(KIND = IP) :: 	my_maxCol, my_minCol, my_range
		INTEGER(KIND = IP) ::	index_start, index_range,  RowStart, RowFinal, &
								Row, Col, start, final, my_start, &
								minCol(world_size), maxCol(world_size), nnz
		INTEGER :: n, i
		INTEGER(KIND=IP), DIMENSION(nnz) :: RowInd, ColInd
		COMPLEX(KIND=CP), DIMENSION(nnz) :: Amat
		COMPLEX(KIND=CP), DIMENSION(nrows) :: x
		COMPLEX(KIND=CP), DIMENSION(nrows) :: y
		INTEGER(KIND = IP), DIMENSION(world_size), INTENT(IN) :: everyones_rows
		INTEGER(KIND = IP), DIMENSION(world_size+1) :: cumulative_rows
		INTEGER :: i1, ws
		COMPLEX(KIND = CP), ALLOCATABLE :: temp(:)
		REAL(KIND =RP) :: finish, tstart, start_mkl, finish_mkl
		CHARACTER(len= 1) :: transa
		
		IF (.not. ALLOCATED(scounts_n)) THEN
			CALL MPIvars_init ( )
		ENDIF
		
		
		cumulative_rows(1)=1
		DO ws=2, world_size+1
			cumulative_rows(ws) = cumulative_rows(ws-1) + everyones_rows(ws-1)
		ENDDO
		
		!----------------------------------------------------------------------------------------
		! Redistribute b and x to correspoinding processors
		
		! Multiplication require x(minval(RowInd):max(RowInd)) be on my processor
		! This snippet calculates the _counts and _displs required for that process
		!
		! This only needs to be done once, the first time the code is called
		! Unfortunately it only works if there is one matrix around.
		IF ((calc_scounts_n) .or. (calc_scounts_t)) THEN
			my_maxCol = maxval(ColInd)
			my_minCol = minval(ColInd)			
			
			my_range = my_maxCol - my_minCol+1
			
			CALL MPI_ALLGATHER(my_minCol, 1, mp_int, minCol, 1, mp_int, comm, ierr)
			CALL MPI_ALLGATHER(my_maxCol, 1, mp_int, maxCol, 1, mp_int, comm, ierr)
			
			my_start = cumulative_rows(my_id+1)-1
			scounts(:) = 0
			sdispls(:) = 0
			! This loop is of O(N)
			DO i1 = 1, nrows
				DO ws = 1, world_size
					IF (((i1+my_start).ge.minCol(ws)) .AND. ((i1+my_start).le.maxCol(ws))) THEN
						IF ((i1+my_start).eq.minCol(ws)) sdispls(ws) = i1-1
						scounts(ws) = scounts(ws)+1
					ENDIF
				ENDDO
			ENDDO

			CALL MPI_ALLTOALL(scounts, 1, mp_int, recvcounts, 1, mp_int, comm, ierr)

			CALL calculate_displs(recvcounts, recvdispls) !O(np)
			
			IF (transa.eq.'N') THEN
				calc_scounts_n = .false. ! Set this to true whenever counts 
										 ! dipls are required to be calculated
				my_range_n = my_range
				scounts_n = scounts
				sdispls_n = sdispls
				recvcounts_n = recvcounts
				recvdispls_n = recvdispls
			ELSEIF (transa.eq.'T') THEN
				calc_scounts_t = .false. ! Set this to true whenever counts 
										 ! dipls are required to be calculated
				my_range_t = my_range
				scounts_t = scounts
				sdispls_t = sdispls
				recvcounts_t = recvcounts
				recvdispls_t = recvdispls
			ENDIF
		ENDIF 
		!----------------------------------------------------------------------------------------

		IF (transa.eq.'N') THEN
			my_range = my_range_n
			scounts = scounts_n
			sdispls = sdispls_n
			recvcounts = recvcounts_n
			recvdispls = recvdispls_n
		ELSEIF (transa.eq.'T') THEN
			my_range = my_range_t
			scounts = scounts_t
			sdispls = sdispls_t
			recvcounts = recvcounts_t
			recvdispls = recvdispls_t
		ENDIF
		
		my_maxCol = maxval(ColInd)
		my_minCol = minval(ColInd)
		
		ALLOCATE(temp(my_range))
		
	
		IF (world_size.eq.1) THEN
			index_start = 1
		ENDIF
		
		IF (world_size.eq.1) THEN
			temp(:) = x(:)
		ELSE
			start_mkl = dsecnd()
			CALL MPI_ALLTOALLV(	x(1), int(scounts, KIND = 4), int(sdispls, KIND = 4), &
								mp_complex, temp(1), int(recvcounts, KIND = 4), &
								int(recvdispls, KIND = 4), mp_complex, comm, ierr)
			finish_mkl = dsecnd()
		ENDIF
			
		start = minval(RowInd)
		final = maxval(RowInd)
				
		y(:) = 0.D0
		
		! Multiply
		
		CALL cpu_time(tstart)
		DO n = 1, nnz
			Row = RowInd(n) - start + 1
			Col = ColInd(n) - my_minCol + 1
			y(Row) = y(Row) + Amat(n)*temp(Col) ! this is not thread safe and cannot be accelerated.
		ENDDO
		CALL cpu_time(finish)
		
		my_Total_comm_mkl = my_Total_comm_mkl + finish_mkl - start_mkl
		my_Total_mul_time = my_Total_mul_time + finish - tstart
		
		RowStart = start - index_start
		RowFinal = RowStart + final - start
		
		DEALLOCATE(temp)
		
	END SUBROUTINE matvectcoo
	
END MODULE COO_routines
