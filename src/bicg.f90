MODULE bicg

	!	DO NOT CHANGE THESE VARIABLES INSIDE YOUR SUBROUTINES OR PROGRAMS!!! ---------------------
	!	These are defined in mpi_init subroutine of the mp_module
	USE mp_module ! Defined mp-global variables inside
	! --------------------------------------------------------------------------------------- 
	CONTAINS
	
	SUBROUTINE mpibicgcoo (	my_Amat, my_RowInd, my_ColInd, my_x, my_b, my_nnz, &
							my_nrows, nrows, everyones_rows, tol, maxit, resvec, flag, &
							iterpause, restartfile)

!		----------------------------------------------------
!		SUBROUTINE DISCRIPTION
!		----------------------------------------------------
!
!		Parallelized version of Bi-Conjugate Gradient method to solve
!		System of linear equations A*x = b for structurally symmetric matrices

!		Conventions :-
!			Variables with prefix my_* are local variables, and need 
!			to be calculated a priori
!			Variables with postfix *_t are variables for transpose matrix

!		----------------------------------------------------
!		INPUT DATA DISCRIPTION
!		----------------------------------------------------

!		(my_Amat, my_RowInd, my_ColInd) is the sparse COO representation of local A-matrix 
!			stored on local machine
!		my_x - Local x vector (For zero-vector initial guess my_x(1:my_nrows) = 0.D0) 
!		my_b - Local b vector
!		my_nnz - number of elements in my_Amat (Number of non-zeros in full - Amatrix)
!		my_nrows - Number of rows of Amat contained on local machine
!		nrows - Total number of rows (hence, my_nrows \approx nrows/np)
!			where, np is the number of processors
!		everyones_rows - rows contained on every processor. Can be calculated using 
!			MPI_ALLGATHER procedure on my_nrows
!		cumulative_rows - cumulative representation on everyones_rows
!		tol - tolerance for BICG
!		maxit - maximum number of iterations for BICG

!		----------------------------------------------------
!		OUTPUT DATA DISCRIPTION
!		----------------------------------------------------

!		my_x - Local solution x
!		resvec - vector of residuals calculated at every iteration
!		flag - flag = 0 Biconjugate gradient converged
!			   flag = 1 BICG did not converge
!			   flag = 4 One of the scalar quantities calculated during bicg became too 
!						small or too large to continue computing.
!						(TIP: Non-dimensionalizing the LA problem helps)

!		-----------------------------------------------------

		USE kinds
		USE essentials
		USE COO_routines
		USE blas95
		USE f95_precision
		IMPLICIT NONE
		INCLUDE 'mpif.h' ! MPI header file
		INCLUDE 'mkl.fi' ! MKL header file
		!---------------------------------
		! MPI variables
		INTEGER(KIND = IP), DIMENSION(world_size) :: scounts, sdispls, &
													 recvcounts, recvdispls
		!------------------------------
		! Matrix/vector distribution variables
		INTEGER(KIND=IP) :: ierr
		INTEGER(KIND=IP) :: my_nnz, my_nrows, nrows
		INTEGER(KIND=IP) :: counter, ws, i1 ! counters
		INTEGER(KIND=IP), DIMENSION(my_nnz) :: my_colInd, my_RowInd
		INTEGER(KIND=IP) :: my_start, my_size
		INTEGER(KIND=IP), DIMENSION(world_size) :: minCol, maxCol 
		COMPLEX(KIND=CP), DIMENSION(my_nnz) :: my_Amat
		INTEGER(KIND=IP), ALLOCATABLE, DIMENSION(:) :: my_RowInd_t, my_ColInd_t
		COMPLEX(KIND=CP), ALLOCATABLE, DIMENSION(:) :: my_Amat_t
		INTEGER(KIND=IP), DIMENSION(world_size) :: everyones_rows
		INTEGER(KIND=IP), DIMENSION(world_size+1) :: cumulative_rows
		COMPLEX(KIND=CP), DIMENSION(my_nrows) :: my_x, my_b
		!-------------------------------------------
		! bicg variables 
		REAL(KIND=RP) :: tol
		INTEGER(KIND=4) :: maxit
		REAL(KIND=RP) :: my_n2b, n2b, tolb, my_normr, normr
		INTEGER(KIND=IP) :: iter
		INTEGER(KIND=4)  ::	flag
		INTEGER(KIND=IP) :: ii ! counter
		LOGICAL :: call_terminate
		REAL(KIND=RP) :: relres 
		COMPLEX(KIND=CP), DIMENSION(my_nrows) :: my_r, my_rt, my_z, my_zt, my_p, my_q, my_pt, my_qt, my_y, my_yt
		REAL(KIND=4), DIMENSION(maxit+1) :: resvec
		COMPLEX(KIND=CP) :: beta, alpha, my_ptq, ptq, my_rho, rho, rho1
		!-------------------------------------------
		! Timing variables
		REAL(KIND=RP) :: ci, cf, c1, c2, cTranspose_Time
		REAL(KIND=RP) :: cTotal_Time, rate, my_cTranspose_Time, my_cTotal_Time
		INTEGER(KIND=IP) :: cr, cm, s1, s2, si , sf, sTranspose_Time, sTotal_Time, communication_time 
		LOGICAL :: count_time
		REAL(KIND=RP) :: finish, start, start_mkl, finish_mkl
		!-------------------------------------------
		! resatart variables
		INTEGER :: iterpause
		CHARACTER(len = 256) :: restartfile

		resvec(:) = 0.D0
		! Default values
		call_terminate =.false.
		count_time = .false. ! set this to false to set timing off
		
		! Assuming that the data is distributed to their respective nodes, start from Here!!
		
		CALL check_mpinit( ) ! check if MPI has been called
		
		! calculate the transpose and distribute it =============================================
		! Use this block of code to calculate the transpose of a matrix and
		! redistribute it to its children
		
		! Amat is rearranged such that Amat_t is in continuous memory.
		! i.e. vector Amat is transformed such that it is stored as if it
		! were in CSC format 
		
		cumulative_rows(1)=1
		DO ws=2, world_size+1
			cumulative_rows(ws) = cumulative_rows(ws-1) + everyones_rows(ws-1)
		ENDDO
		
		CALL cpu_time(start)
		CALL heapsort(my_ColInd, my_RowInd, my_Amat, my_nnz) ! O(n*log(n))
		CALL cpu_time(finish)
						
		counter = 0
		ws = 1
		i1 = 1
		scounts(:) = 0
		
		CALL MPI_BARRIER(comm, ierr)
		DO WHILE (i1.le.my_nnz) !O(n)
			IF ((my_ColInd(i1).ge.cumulative_rows(ws)) .AND. &
				(my_ColInd(i1).lt.cumulative_rows(ws+1))) THEN
				counter = counter+1
				scounts(ws) = counter
				i1 = i1 + 1
			ELSE
				counter = 0
				ws = ws + 1
			ENDIF
		ENDDO
		
		start_mkl = dsecnd()
		CALL MPI_ALLTOALL(scounts, 1, mp_int, recvcounts, 1, mp_int, comm, ierr)
		finish_mkl = dsecnd()
		
		my_Total_comm_mkl = my_Total_comm_mkl + finish_mkl - start_mkl


		CALL calculate_displs(scounts, sdispls) !O(np)
		CALL calculate_displs(recvcounts, recvdispls) !O(np)
		
		my_size = sum(recvcounts)

		ALLOCATE(my_Amat_t(my_size), my_RowInd_t(my_size), my_ColInd_t(my_size))
		
		! Distribute Amat and corresponding Row and Column indices
		start_mkl = dsecnd()
		CALL MPI_ALLTOALLV (my_Amat, int(scounts, KIND = 4), int(sdispls, KIND = 4), &
							mp_complex, my_Amat_t, int(recvcounts, KIND = 4), int(recvdispls, KIND = 4), &
							mp_complex, comm, ierr)
		CALL MPI_ALLTOALLV (my_ColInd, int(scounts, KIND = 4), int(sdispls, KIND = 4), &
							mp_int, my_RowInd_t, int(recvcounts, KIND = 4), int(recvdispls, KIND = 4), &
							mp_int, comm, ierr)
		CALL MPI_ALLTOALLV (my_RowInd, int(scounts, KIND = 4), int(sdispls, KIND = 4), &
							mp_int, my_ColInd_t, int(recvcounts, KIND = 4), int(recvdispls, KIND = 4), &
							mp_int, comm, ierr)
		finish_mkl = dsecnd()
		
		my_Total_comm_mkl = my_Total_comm_mkl + finish_mkl - start_mkl
		
		my_Amat_t = conjg(my_Amat_t)
		CALL heapsort(my_RowInd, my_ColInd, my_Amat, my_nnz) ! O(n*log(n))
		
		
		
		CALL heapsort(my_RowInd_t, my_ColInd_t, my_Amat_t, my_size) ! O(n*log(n))
				
		calc_scounts_n = .true.
		calc_scounts_t = .true.
		
	
		! BEGIN BI-CONJUGATE GRADIENT
		! The code from here on closely represents the BICG matlab routine. Although, this version
		! utilizes MPI-capabilities.
		! Ref: Templates for the Solution of Linear Systems: Building Blocks for Iterative Methods by Berrett et. al.
		
		! Checking vaidity of tolerance values
		IF (root_node) THEN
			IF (tol.lt.eps) THEN
				print *, 'WARNING: Tolerance is below machine error'
				tol=eps
			ELSEIF (tol.ge.1) THEN
				print *, 'WARNING: Tolerance is too large'
				tol= 1-eps
			ENDIF
		ENDIF
		
		CALL MPI_BCAST(tol, 1, mp_real, root_process, comm, ierr)

		my_n2b=nrm2(my_b)
		my_n2b = my_n2b*my_n2b

		start_mkl = dsecnd()
		CALL MPI_ALLREDUCE(my_n2b, n2b, 1, mp_real, mp_sumr, comm, ierr)
		finish_mkl = dsecnd()
		
		my_Total_comm_mkl = my_Total_comm_mkl + finish_mkl - start_mkl
		
		n2b = sqrt(n2b)
		
		ii = 0
		
		IF (root_node) THEN
			IF (n2b.eq.0) THEN
				my_x(:)= 0.D0
				flag=0
				relres=0.D0
				iter=0
				resvec(1)=0.D0
				print *, 'NORM2(b)=0, hence the solution is a trivial solution'
				call_terminate=.true.
			ENDIF
		ENDIF
		
		CALL MPI_BCAST(call_terminate, 1, mp_logical, root_process, comm, ierr)
		
		IF (call_terminate) THEN
			CALL terminate(flag, call_terminate, normr, n2b, ii)
			RETURN
		ENDIF
		
		! Set up for the method
		flag= 1
		tolb=tol*n2b ! normalized tolerance
		rho=1
		CALL matvectcoo(my_Amat, my_RowInd, my_ColInd, my_nnz, &
							my_x, my_r, my_nrows, everyones_rows, 'N')
!		CALL matvectcoo(my_Amat, my_RowInd, my_ColInd, temp, my_nnz, my_minCol, my_range, my_r, my_nrows)

		my_r= my_b - my_r
		my_normr= nrm2(my_r)
		my_normr = my_normr*my_normr
		
		start_mkl = dsecnd()
		CALL MPI_ALLREDUCE(my_normr, normr, 1, mp_real, mp_sumr, comm, ierr)
		finish_mkl = dsecnd()
		
		my_Total_comm_mkl = my_Total_comm_mkl + finish_mkl - start_mkl
		
		normr = sqrt(normr) ! Calculated error
				
		IF (normr .le. tolb)THEN	! Initial guess was good enough
			IF (root_node) THEN
				flag=0
				relres= normr/n2b
				iter=0
				resvec(1)=normr
				call_terminate=.true.
			ENDIF
		ENDIF
		
!		IF (ISNAN(normr)) THEN
!			flag = 2
!			call_terminate = .true.
!		ENDIF

		CALL MPI_BCAST(call_terminate, 1, mp_logical, root_process, comm, ierr)
		
		IF (call_terminate) THEN
			CALL terminate(flag, call_terminate, normr, n2b, ii)
			RETURN
		ENDIF
		
		my_rt= my_r ! shadow residual
		resvec(1)= normr


		DO ii=1,maxit
			my_y=my_r
			my_z=my_y
			my_yt=my_rt
			my_zt=my_yt
			rho1=rho
			my_rho = dotc(my_rt, my_z)
			
			start_mkl = dsecnd()
			CALL MPI_ALLREDUCE(my_rho, rho, 1, mp_complex, mp_sumc, comm, ierr) 
			finish_mkl = dsecnd()
		
			my_Total_comm_mkl = my_Total_comm_mkl + finish_mkl - start_mkl
			
			
			IF (rho.eq.0) THEN! need to add the infinity portion here
				IF (root_node) THEN
					flag=4
					call_terminate=.true.
				ENDIF
				EXIT
			ENDIF
			IF (ii.eq.1) THEN
				my_p=my_z
				my_pt=my_zt
			ELSE
				beta= rho/rho1 			! rho is gathered from all children
				IF (beta.eq.0) THEN		! need to add the infinity portion here
					IF (root_node) THEN
						flag=4
						call_terminate=.true.
					ENDIF
					EXIT
				ENDIF
				my_p= my_z+beta*my_p
				my_pt= my_zt+ conjg(beta)*my_pt
			ENDIF ! ii IF
			
!			IF (count_time .AND. root_process) CALL system_clock(s1)	!Timing	
!			CALL MPI_ALLTOALLV(my_p, scounts, sdispls, mp_complex, temp, recvcounts, &
!							   recvdispls, mp_complex, comm, ierr)
!			IF (count_time .AND. root_process) CALL system_clock(s2)  	!Timing
!			IF (count_time .AND. root_process) communication_time = communication_time + s2-s1
			
			CALL matvectcoo(my_Amat, my_RowInd, my_ColInd, my_nnz, &
								my_p, my_q, my_nrows, everyones_rows, 'N')
!			CALL matvectcoo(my_Amat, my_RowInd, my_ColInd, temp, my_nnz, my_minCol, my_range, my_q, my_nrows)
			
			
!			IF (count_time .AND. root_process) CALL system_clock(s1)	!Timing	
!			CALL MPI_ALLTOALLV(my_pt, scounts, sdispls, mp_complex, temp, recvcounts, &
!							   recvdispls, mp_complex, comm, ierr)
!			IF (count_time .AND. root_process) CALL system_clock(s2)  	!Timing
!			IF (count_time .AND. root_process) communication_time = communication_time + s2-s1
			
			
			CALL matvectcoo(my_Amat_t, my_RowInd_t, my_ColInd_t, my_size, &
								my_pt, my_qt, my_nrows, everyones_rows, 'T')
!			CALL matvectcoo(my_Amat_t, my_RowInd_t, my_ColInd_t, temp, my_nnz, my_minCol, my_range, my_qt, my_nrows)
			
			my_ptq=dotc(my_pt, my_q)
			
			start_mkl = dsecnd()
			CALL MPI_ALLREDUCE(my_ptq, ptq, 1, mp_complex, mp_sumc, comm, ierr)
			finish_mkl = dsecnd()
		
			my_Total_comm_mkl = my_Total_comm_mkl + finish_mkl - start_mkl
					
			IF (ptq.eq.0) THEN! need to add the infinity portion here
				IF (root_node) THEN
					flag=4
					call_terminate=.true.
				ENDIF
				EXIT
			ELSE
				alpha=rho/ptq
			ENDIF
			
			my_x= my_x+alpha*my_p
			my_r= my_r-alpha*my_q
			my_rt= my_rt - conjg(alpha)*my_qt ! Need to check this with MATLAB code
			my_normr= nrm2(my_r)
			my_normr = my_normr*my_normr
			
			
			start_mkl = dsecnd()			
			CALL MPI_ALLREDUCE(my_normr, normr, 1, mp_real, mp_sumr, comm, ierr)
			finish_mkl = dsecnd()
		
			my_Total_comm_mkl = my_Total_comm_mkl + finish_mkl - start_mkl
			
			normr = sqrt(normr)
			resvec(ii+1) = normr

			IF (root_node) THEN
				IF (MOD(ii,100).eq.0.0) THEN
					WRITE(stdout, '(a)') '	'
					WRITE(stdout, '(a, I6, a, E10.3)') 'Iter #', ii, ' relres = ', normr/tolb
				ENDIF
			ENDIF

			! IF (normr.gt.HUGE(normr)) THEN
				! IF (root_node) THEN
					! flag = 4
					! call_terminate = .true.
				! ENDIF
				! EXIT
			! ENDIF
			
			! IF (normr.ne.normr) THEN ! checking for NaN values
				! IF (root_node) THEN
					! flag = 4
					! call_terminate = .true.
				! ENDIF
				! EXIT
			! ENDIF

			! check for convergence
			IF (normr.le.tolb) THEN
				IF (my_id.eq.root_process) THEN
					flag=0
					call_terminate=.true.
				ENDIF
				EXIT
			ENDIF
						
		ENDDO

		CALL MPI_BCAST(call_terminate, 1, mp_logical, root_process, comm, ierr)
		
		IF ((.NOT. call_terminate) .AND. (flag.ne.1) .AND. (my_id.eq.root_process)) THEN
			print *, 'ERROR: flag==0 and call_terminate==false'
		ENDIF
				
		CALL terminate(flag, call_terminate, normr, n2b, ii)
		
		!----------------!
			
10	FORMAT(E12.4)
11 	FORMAT(i5)

	END SUBROUTINE mpibicgcoo
	
	SUBROUTINE mpibicgstabcoo (	my_Amat, my_RowInd, my_ColInd, my_x, my_b, my_nnz, &
							my_nrows, nrows, everyones_rows, tol, maxit, resvec, flag, &
							iterpause, restartfile)

!		----------------------------------------------------
!		SUBROUTINE DISCRIPTION
!		----------------------------------------------------
!
!		Parallelized version of Bi-Conjugate Gradient method to solve
!		System of linear equations A*x = b for structurally symmetric matrices

!		Conventions :-
!			Variables with prefix my_* are local variables, and need 
!			to be calculated a priori
!			Variables with postfix *_t are variables for transpose matrix

!		----------------------------------------------------
!		INPUT DATA DISCRIPTION
!		----------------------------------------------------

!		(my_Amat, my_RowInd, my_ColInd) is the sparse COO representation of local A-matrix 
!			stored on local machine
!		my_x - Local x vector (For zero-vector initial guess my_x(1:my_nrows) = 0.D0) 
!		my_b - Local b vector
!		my_nnz - number of elements in my_Amat (Number of non-zeros in full - Amatrix)
!		my_nrows - Number of rows of Amat contained on local machine
!		nrows - Total number of rows (hence, my_nrows \approx nrows/np)
!			where, np is the number of processors
!		everyones_rows - rows contained on every processor. Can be calculated using 
!			MPI_ALLGATHER procedure on my_nrows
!		cumulative_rows - cumulative representation on everyones_rows
!		tol - tolerance for BICG
!		maxit - maximum number of iterations for BICG

!		----------------------------------------------------
!		OUTPUT DATA DISCRIPTION
!		----------------------------------------------------

!		my_x - Local solution x
!		resvec - vector of residuals calculated at every iteration
!		flag - flag = 0 Biconjugate gradient converged
!			   flag = 1 BICG did not converge
!			   flag = 4 One of the scalar quantities calculated during bicg became too 
!						small or too large to continue computing.
!						(TIP: Non-dimensionalizing the LA problem helps)

!		-----------------------------------------------------

		USE kinds
		USE essentials
		USE COO_routines
		USE blas95
		USE f95_precision
		IMPLICIT NONE
		INCLUDE 'mpif.h' ! MPI header file
		INCLUDE 'mkl.fi' ! MKL header file
		!---------------------------------
		! MPI variables
		INTEGER(KIND = IP), DIMENSION(world_size) :: scounts, sdispls, &
													 recvcounts, recvdispls
		!------------------------------
		! Matrix/vector distribution variables
		INTEGER(KIND=IP) :: ierr
		INTEGER(KIND=IP) :: my_nnz, my_nrows, nrows
		INTEGER(KIND=IP) :: counter, ws, i1 ! counters
		INTEGER(KIND=IP), DIMENSION(my_nnz) :: my_colInd, my_RowInd
		INTEGER(KIND=IP) :: my_start, my_size
		INTEGER(KIND=IP), DIMENSION(world_size) :: minCol, maxCol 
		COMPLEX(KIND=CP), DIMENSION(my_nnz) :: my_Amat
		INTEGER(KIND=IP), ALLOCATABLE, DIMENSION(:) :: my_RowInd_t, my_ColInd_t
		COMPLEX(KIND=CP), ALLOCATABLE, DIMENSION(:) :: my_Amat_t
		INTEGER(KIND=IP), DIMENSION(world_size) :: everyones_rows
		INTEGER(KIND=IP), DIMENSION(world_size+1) :: cumulative_rows
		COMPLEX(KIND=CP), DIMENSION(my_nrows) :: my_x, my_b
		!-------------------------------------------
		! bicg variables 
		REAL(KIND=RP) :: tol
		INTEGER(KIND=4) :: maxit
		REAL(KIND=RP) :: my_n2b, n2b, tolb, my_normr, normr
		INTEGER(KIND=IP) :: iter
		INTEGER(KIND=4)  ::	flag
		INTEGER(KIND=IP) :: ii ! counter
		LOGICAL :: call_terminate
		REAL(KIND=RP) :: relres 
		COMPLEX(KIND=CP), DIMENSION(my_nrows) :: my_r, my_rt, my_p, my_xh, my_s, my_sh, &
												 my_ph, my_t, my_v
		REAL(KIND=4), DIMENSION(2*maxit+1) :: resvec
		COMPLEX(KIND=CP) :: beta, alpha, my_rtv, rtv, my_rho, rho, rho1, omega, my_tts, &
							my_ttt, tts, ttt
		!-------------------------------------------
		! Timing variables
		REAL(KIND=RP) :: ci, cf, c1, c2, cTranspose_Time
		REAL(KIND=RP) :: cTotal_Time, rate, my_cTranspose_Time, my_cTotal_Time
		INTEGER(KIND=IP) :: cr, cm, s1, s2, si , sf, sTranspose_Time, sTotal_Time, communication_time 
		LOGICAL :: count_time
		REAL(KIND=RP) :: finish, start, start_mkl, finish_mkl
		!--------------------------------------------
		! restart variables
		! resatart variables
		INTEGER :: iterpause
		CHARACTER(len = 256) :: restartfile


		resvec(:) = 0.D0
		! Default values
		call_terminate =.false.
		count_time = .false. ! set this to false to set timing off
		
		! Assuming that the data is distributed to their respective nodes, start from Here!!
		
		CALL check_mpinit( ) ! check if MPI has been called
		
		calc_scounts_n = .true.
		calc_scounts_t = .false.
		
	
		! BEGIN BI-CONJUGATE GRADIENT
		! The code from here on closely represents the BICG matlab routine. Although, this version
		! utilizes MPI-capabilities.
		! Ref: Templates for the Solution of Linear Systems: Building Blocks for Iterative Methods by Berrett et. al.
		
		! Checking vaidity of tolerance values
		IF (root_node) THEN
			IF (tol.lt.eps) THEN
				print *, 'WARNING: Tolerance is below machine error'
				tol=eps
			ELSEIF (tol.ge.1) THEN
				print *, 'WARNING: Tolerance is too large'
				tol= 1-eps
			ENDIF
		ENDIF
		
		CALL MPI_BCAST(tol, 1, mp_real, root_process, comm, ierr)

		my_n2b=nrm2(my_b)
		my_n2b = my_n2b*my_n2b

		start_mkl = dsecnd()
		CALL MPI_ALLREDUCE(my_n2b, n2b, 1, mp_real, mp_sumr, comm, ierr)
		finish_mkl = dsecnd()
		
		my_Total_comm_mkl = my_Total_comm_mkl + finish_mkl - start_mkl
		
		n2b = sqrt(n2b)
		
		ii = 0
		
		IF (root_node) THEN
			IF (n2b.eq.0) THEN
				my_x(:)= 0.D0
				flag=0
				relres=0.D0
				iter=0
				resvec(1)=0.D0
				print *, 'NORM2(b)=0, hence the solution is a trivial solution'
				call_terminate=.true.
			ENDIF
		ENDIF
		
		CALL MPI_BCAST(call_terminate, 1, mp_logical, root_process, comm, ierr)
		
		IF (call_terminate) THEN
			CALL terminate(flag, call_terminate, normr, n2b, ii)
			RETURN
		ENDIF
		
		! Set up for the method
		flag= 1
		tolb=tol*n2b ! normalized tolerance
		rho=1
		omega = 1
		alpha = 0
		CALL matvectcoo(my_Amat, my_RowInd, my_ColInd, my_nnz, &
							my_x, my_r, my_nrows, everyones_rows, 'N')
!		CALL matvectcoo(my_Amat, my_RowInd, my_ColInd, temp, my_nnz, my_minCol, my_range, my_r, my_nrows)

		my_r= my_b - my_r
		my_normr= nrm2(my_r)
		my_normr = my_normr*my_normr
		
		start_mkl = dsecnd()
		CALL MPI_ALLREDUCE(my_normr, normr, 1, mp_real, mp_sumr, comm, ierr)
		finish_mkl = dsecnd()
		
		my_Total_comm_mkl = my_Total_comm_mkl + finish_mkl - start_mkl
		
		normr = sqrt(normr) ! Calculated error
				
		IF (normr .le. tolb)THEN	! Initial guess was good enough
			IF (root_node) THEN
				flag=0
				relres= normr/n2b
				iter=0
				resvec(1)=normr
				call_terminate=.true.
			ENDIF
		ENDIF
		
!		IF (ISNAN(normr)) THEN
!			flag = 2
!			call_terminate = .true.
!		ENDIF

		CALL MPI_BCAST(call_terminate, 1, mp_logical, root_process, comm, ierr)
		
		IF (call_terminate) THEN
			CALL terminate(flag, call_terminate, normr, n2b, ii)
			RETURN
		ENDIF
		
		my_rt= my_r ! shadow residual
		resvec(1)= normr
		


		DO ii=1,maxit
			rho1=rho
			my_rho = dotc(my_rt, my_r)
			
			start_mkl = dsecnd()
			CALL MPI_ALLREDUCE(my_rho, rho, 1, mp_complex, mp_sumc, comm, ierr) 
			finish_mkl = dsecnd()
		
			my_Total_comm_mkl = my_Total_comm_mkl + finish_mkl - start_mkl
			
			
			IF (rho.eq.0) THEN! need to add the infinity portion here
				IF (root_node) THEN
					flag=4
					call_terminate=.true.
				ENDIF
				EXIT
			ENDIF
			IF (ii.eq.1) THEN
				my_p=my_r
			ELSE
				beta= (rho/rho1)*(alpha/omega) 			! rho is gathered from all children
				IF (beta.eq.0) THEN		! need to add the infinity portion here
					IF (root_node) THEN
						flag=4
						call_terminate=.true.
					ENDIF
					EXIT
				ENDIF
				my_p= my_r+beta*(my_p - omega*my_v)
			ENDIF ! ii IF
			my_ph = my_p
			
			CALL matvectcoo(my_Amat, my_RowInd, my_ColInd, my_nnz, &
							my_ph, my_v, my_nrows, everyones_rows, 'N')
			
			my_rtv=dotc(my_rt, my_v)
			
			start_mkl = dsecnd()
			CALL MPI_ALLREDUCE(my_rtv, rtv, 1, mp_complex, mp_sumc, comm, ierr)
			finish_mkl = dsecnd()
		
			my_Total_comm_mkl = my_Total_comm_mkl + finish_mkl - start_mkl
					
			IF (rtv.eq.0) THEN! need to add the infinity portion here
				IF (root_node) THEN
					flag=4
					call_terminate=.true.
				ENDIF
				EXIT
			ENDIF
			alpha = rho/rtv
			
			my_xh= my_x+alpha*my_ph
			my_s= my_r-alpha*my_v
			my_normr= nrm2(my_s)
			my_normr = my_normr*my_normr
			
!			Write down the solution to the file every iterpause iterations
			IF (MOD(ii,iterpause).eq.0.0) THEN
				CALL MPI_BARRIER(comm, ierr)
				OPEN(unit = 639, file = restartfile, form = 'unformatted')
				WRITE(639) my_x
				CLOSE(unit = 639)
			ENDIF
			
			start_mkl = dsecnd()			
			CALL MPI_ALLREDUCE(my_normr, normr, 1, mp_real, mp_sumr, comm, ierr)
			finish_mkl = dsecnd()
		
			my_Total_comm_mkl = my_Total_comm_mkl + finish_mkl - start_mkl
			
			normr = sqrt(normr)
			resvec(2*ii) = normr

			IF (root_node) THEN
				IF (MOD(ii,100).eq.0.0) THEN
					WRITE(stdout, '(a)') '	'
					WRITE(stdout, '(a, I6, a, E10.3)') 'Iter #', ii, ' relres = ', normr/tolb
				ENDIF
			ENDIF

			! check for convergence
			IF (normr.le.tolb) THEN
				my_x = my_x + alpha*my_p
				CALL MPI_BARRIER(comm, ierr)
				OPEN(unit = 639, file = restartfile, form = 'unformatted')
				WRITE(639) my_x
				CLOSE(unit = 639)
				IF (my_id.eq.root_process) THEN
					flag=0
					call_terminate=.true.
				ENDIF
				EXIT
			ENDIF
			
			my_sh = my_s
			CALL matvectcoo(my_Amat, my_RowInd, my_ColInd, my_nnz, &
							my_sh, my_t, my_nrows, everyones_rows, 'N')
			
			my_tts = dotc(my_t, my_s)
			
			start_mkl = dsecnd()
			CALL MPI_ALLREDUCE(my_tts, tts, 1, mp_complex, mp_sumc, comm, ierr) 
			finish_mkl = dsecnd()
		
			my_Total_comm_mkl = my_Total_comm_mkl + finish_mkl - start_mkl
			
			my_ttt = dotc(my_t, my_t)
			
			start_mkl = dsecnd()
			CALL MPI_ALLREDUCE(my_ttt, ttt, 1, mp_complex, mp_sumc, comm, ierr) 
			finish_mkl = dsecnd()
		
			my_Total_comm_mkl = my_Total_comm_mkl + finish_mkl - start_mkl
			
			omega = tts/ttt
			
			my_x = my_xh + omega*my_sh
			my_r = my_s - omega*my_t
			
			my_normr = nrm2(my_r)
			my_normr = my_normr*my_normr
			
			start_mkl = dsecnd()			
			CALL MPI_ALLREDUCE(my_normr, normr, 1, mp_real, mp_sumr, comm, ierr)
			finish_mkl = dsecnd()
		
			my_Total_comm_mkl = my_Total_comm_mkl + finish_mkl - start_mkl
			
			normr = sqrt(normr)
			resvec(2*ii+1) = normr
			
			IF (normr.le.tolb) THEN
				my_x = my_x + alpha*my_p
				CALL MPI_BARRIER(comm, ierr)
				OPEN(unit = 639, file = restartfile, form = 'unformatted')
				WRITE(639) my_x
				CLOSE(unit = 639)
				IF (my_id.eq.root_process) THEN
					flag=0
					call_terminate=.true.
				ENDIF
				EXIT
			ENDIF
			
		ENDDO

		CALL MPI_BCAST(call_terminate, 1, mp_logical, root_process, comm, ierr)
		
		IF ((.NOT. call_terminate) .AND. (flag.ne.1) .AND. (my_id.eq.root_process)) THEN
			print *, 'ERROR: flag==0 and call_terminate==false'
		ENDIF
				
		CALL terminate(flag, call_terminate, normr, n2b, ii)
		
		!----------------!
			
10	FORMAT(E12.4)
11 	FORMAT(i5)

	END SUBROUTINE mpibicgstabcoo
	
	SUBROUTINE terminate(flag, call_terminate, normr, n2b, iter)
	
		USE kinds, ONLY: RP
		IMPLICIT NONE
		INTEGER(KIND=4)		:: flag
		LOGICAL					:: call_terminate
		REAL(KIND=RP) 			:: normr, n2b
		INTEGER(KIND=IP)		:: iter
	
		IF (io_node) THEN
			IF ((flag.eq.0) .AND. (call_terminate)) THEN
				print *, 'Biconjugate Gradient converged with relative residual = ', normr/n2b, 'on iteration number ', iter
			ENDIF
		ENDIF
	
		IF ((.NOT. call_terminate) .AND. (flag.eq.0)) THEN
			IF (my_id.eq.root_process) THEN 
				print *, 'Biconjugate Gradient did not converge (relres = ', normr/n2b, ')'
				print *, 'Writing the current solution for restart'
			ENDIF
			! write the my_x's to a file for restarts
		ENDIF
	
		IF (io_node) THEN
			IF (flag.eq.1) THEN
				write(stdout, '(a, F13.7)') 'Biconjugate gradient did not converge with a final relative residual', &
										normr/n2b
			ENDIF
			
			IF (flag.eq.4) THEN
				WRITE(stdout, '(a)') 'Error blew up'
			ENDIF
		ENDIF
	END SUBROUTINE

end module bicg
