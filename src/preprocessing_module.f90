MODULE preprocessing_module
	
	USE kinds
	USE mp_module
	IMPLICIT NONE
	CONTAINS
	
!	====================================================================================================

	SUBROUTINE get_natsc(at, atc, nat, natsc, crystal_coordinates)
	! determine the number of atoms in the supercell, natsc, using information about the primitive
		
		USE essentials
		IMPLICIT NONE
		
		REAL(KIND = RP)		::	Vsc, V, at(3,3), atsc(3,3)
		REAL				::	atc(3,3)
		INTEGER				::	natsc, nsc, nat
		INTEGER				::	i
		LOGICAL				::	crystal_coordinates
		
		
		IF (crystal_coordinates) THEN
			atc = at 	! If working in crystal coordinates then supercell 
							! is the primitive cell
		ELSE
			atc(:,:) = 0.D0
			DO i = 1, 3
				atc(i,i) = 1.0_RP
			ENDDO
		ENDIF
		
		atsc = atc
		
		CALL cell_volume(at, 1.D0, V)
		
		CALL cell_volume(atsc, 1.D0, Vsc)
		
				
		IF ((Vsc/V-NINT(Vsc/V)).gt.1.0E-4) THEN
			WRITE(*, *) 'ERROR : Volume of super cell ~= n * Volume of primitive lattice'
			STOP
		ENDIF
		

		
		nsc = NINT(Vsc/V)
		
		natsc = nat*nsc
		
	END SUBROUTINE get_natsc
	
!	====================================================================================================
	
	SUBROUTINE get_disp(asr, na_ifc, fd, has_zstar, q, alat1, alat2, at1, at2, &
						nat, ntyp, nr1, nr2, nr3, ibrav, mode, ityp1, ityp2, epsil, zeu1, &
						zeu2, vg, omega1, omega2, frc1, frc2, w2, f_of_q1, f_of_q2, amass1, amass2, &
						tau1, tau2, z)

		USE dispersion
		USE rigid
		USE sum_rule
		IMPLICIT NONE
		
		CHARACTER(len = 256)					::	asr
		LOGICAL									:: 	na_ifc, fd, has_zstar(2)
		REAL(KIND = RP)							:: 	qq, qhat(3), q(3), alat1,alat2, at1(3,3), at2(3,3)
		INTEGER									::	nat(2), nr1, nr2, nr3, ibrav(2), mode, ntyp(2), errore
		INTEGER									::	ityp1(nat(1)), ityp2(nat(2))
		REAL(KIND = RP)							::	epsil(3,3,2), zeu1(3,3,nat(1)), zeu2(3,3,nat(2))
		REAL(KIND = RP)							::	vg(3)
		REAL(KIND = RP)							:: 	omega1, omega2
		REAL(KIND = RP)							::	frc1(nr1,nr2,nr3,3,3,nat(1),nat(1)), &
													frc2(nr1,nr2,nr3,3,3,nat(2),nat(2))
		REAL(KIND = RP)							::	w2(3*nat(1)), f_of_q1(3,3,nat(1),nat(1)), &
													f_of_q2(3,3,nat(2),nat(2)), tau1(3, nat(1)), tau2(3, nat(2)), &
													amass1(ntyp(1)), amass2(ntyp(2))
		REAL(KIND = RP), ALLOCATABLE 			:: 	m_loc(:,:)
		COMPLEX(KIND = CP)						::	z(3*nat(1), 3*nat(1))
		
		!get_qpoint: see next subroutine below; if qpoint has negative group velocity, it pick -q instead.
		call get_qpoint(asr, na_ifc, fd, has_zstar, q, alat1, alat2, at1, at2, &
							nat, ntyp, nr1, nr2, nr3, ibrav, mode, ityp1, ityp2, epsil, zeu1, &
							zeu2, vg, omega1, omega2, frc1, frc2, w2, f_of_q1, f_of_q2, amass1, amass2, &
							tau1, tau2)
		
		!disp: in dispersion.f90 module					
		call disp(	frc1, f_of_q1, tau1, zeu1, m_loc, nr1, nr2, nr3, epsil(:,:,1), nat(1), &
					ibrav(1), alat1, at1, ntyp(1), ityp1, amass1, omega1, &
					has_zstar(1), na_ifc, fd, asr, q, w2, z)
					
		IF (asr /= 'no') THEN
				CALL set_asr (asr, nr1, nr2, nr3, frc2, zeu2, &
					nat(2), ibrav(2), tau2)
		END IF
		
	END SUBROUTINE

!	===============================================================================================================
	
	SUBROUTINE get_qpoint(asr, na_ifc, fd, has_zstar, q, alat1, alat2, at1, at2, &
							nat, ntyp, nr1, nr2, nr3, ibrav, mode, ityp1, ityp2, epsil, zeu1, &
							zeu2, vg, omega1, omega2, frc1, frc2, w2, f_of_q1, f_of_q2, amass1, amass2, &
							tau1, tau2)

! 	In the case that the group velocity is moving away from the Transmitted PML the incident
!	phonon mode needs to be chosen appropriately i.e. if vg(q)[3] < 0 then set q(3) = -q(3)
!	Assuming that the PML is always oriented normal to z-dimention
	
		USE dispersion
		USE rigid
		
		IMPLICIT NONE
		
		CHARACTER(len = 256)					::	asr
		LOGICAL									:: 	na_ifc, fd, has_zstar(2)
		REAL(KIND = RP)							:: 	qq, qhat(3), q(3), alat1,alat2, at1(3,3), at2(3,3)
		INTEGER									::	nat(2), nr1, nr2, nr3, ibrav(2), mode, ntyp(2), errore
		INTEGER									::	ityp1(nat(1)), ityp2(nat(2))
		REAL(KIND = RP)							::	epsil(3,3,2), zeu1(3,3,nat(1)), zeu2(3,3,nat(2))
		REAL(KIND = RP)							::	vg(3)
		REAL(KIND = RP)							:: 	omega1, omega2
		REAL(KIND = RP)							::	frc1(nr1,nr2,nr3,3,3,nat(1),nat(1)), &
													frc2(nr1,nr2,nr3,3,3,nat(2),nat(2))
		REAL(KIND = RP)							::	w2(3*nat(1)), f_of_q1(3,3,nat(1),nat(1)), &
													f_of_q2(3,3,nat(2),nat(2)), tau1(3, nat(1)), tau2(3, nat(2)), &
													amass1(ntyp(1)), amass2(ntyp(2))
		REAL(KIND = RP), ALLOCATABLE 			:: 	m_loc(:,:)
		
		IF (root_node) THEN
		
			IF (na_ifc) THEN
				  qq=sqrt(q(1)**2+q(2)**2+q(3)**3)
				  if(qq == 0.0) qq=1.0
				  qhat(1)=q(1)/qq
				  qhat(2)=q(2)/qq
				  qhat(3)=q(3)/qq
		
				CALL nonanal_ifc(nat(1), nat(1), ityp1, epsil(:,:,1), qhat, zeu1, &
								 omega1, nr1,nr2,nr3,f_of_q1 )
				CALL nonanal_ifc(nat(2), nat(2), ityp2, epsil(:,:,2), qhat, zeu2, &
								 omega2, nr1,nr2,nr3,f_of_q2 )
			ENDIF
			
			CALL Group_velocity	(	frc1, f_of_q1, tau1, zeu1, m_loc, nr1, nr2, nr3, epsil(:,:,1), nat(1), &
					ibrav(1), alat1, at1, ntyp(1), ityp1, amass1, omega1, &
					has_zstar(1), na_ifc, fd, asr, q, vg, mode)
			
			IF (vg(3).lt.0.0) THEN
				q(3) = -q(3)
				vg(3) = -vg(3)
				
				IF (na_ifc) THEN
					  qq=sqrt(q(1)**2+q(2)**2+q(3)**3)
					  if(qq == 0.0) qq=1.0
					  qhat(1)=q(1)/qq
					  qhat(2)=q(2)/qq
					  qhat(3)=q(3)/qq
			
					CALL nonanal_ifc(nat(1), nat(1), ityp1, epsil(:,:,1), qhat, zeu1, &
									 omega1, nr1,nr2,nr3,f_of_q1 )
					CALL nonanal_ifc(nat(2), nat(2), ityp2, epsil(:,:,2), qhat, zeu2, &
									 omega2, nr1,nr2,nr3,f_of_q2 )
				ENDIF
			ENDIF
			
			IF (vg(3).eq.0) THEN
				WRITE(stdout, '(a)') 'Group velocity = 0'
				CALL MPI_ABORT(comm, errore, ierr)
			ENDIF
			
		ENDIF
		
		CALL MPI_BCAST(q, 3, mp_real, root_process, comm, ierr)
		
		CALL MPI_BCAST(f_of_q1, 9*nat(1)*nat(1), mp_real, root_process, comm, ierr)
		CALL MPI_BCAST(f_of_q2, 9*nat(1)*nat(1), mp_real, root_process, comm, ierr)
		
	END SUBROUTINE
	
!	=========================================================================================
	
	SUBROUTINE get_supercell(at1, tau1, ityp1, nat, z, r_cell, na_vec, ib_vec1, &
								ib_vec2, natc, tauc, itypc, atc, zc, crystal_coordinates, &
								nr1, nr2, nr3)
	
		USE essentials
		IMPLICIT NONE
		
		LOGICAL								:: 	crystal_coordinates
		INTEGER								::	i
		REAL(KIND = RP)						::	at1(3,3)
		INTEGER								::	nat(2)
		COMPLEX(KIND = CP)					::	z(3*nat(1))
		INTEGER								::	ityp1(nat(1))
		REAL(KIND = RP) 					:: 	atsc(3,3)
		REAL								::	atc(3,3)
		INTEGER 							:: 	natc
		REAL(KIND = RP)				 		:: 	tausc(3, natc), r_cell(3, natc), tauc(3,natc)
		INTEGER				 				:: 	itypc(natc)
		COMPLEX(KIND = CP)				 	:: 	zsc(3*natc), zc(3*natc)
		INTEGER				 				:: 	ib_vec1(nat(1), -2*nr1:2*nr1, -2*nr2:2*nr2, -2*nr3:2*nr3, natc), &
												ib_vec2(nat(2), -2*nr1:2*nr1, -2*nr2:2*nr2, -2*nr3:2*nr3, natc), &
												na_vec(natc)
		REAL(KIND = RP)						::	tau1(3, nat(1))
		INTEGER								::	nr1, nr2, nr3
	
		atsc = atc
		CALL Supercell(	atsc, at1, tau1, ityp1, nat(1), tauc, itypc, &
						natc, z, zc, r_cell, na_vec, ib_vec1, nr1, nr2, nr3 )
		!	___c variables listed below are current variables and are chosen based on 
		!	the choice of coordinate system you want to work in.
		
		ib_vec2 = ib_vec1
	
		IF (real(zc(1)).LT.0.0) THEN
			zc = CMPLX(-1.0, 0.0, KIND = CP)*zc(:)
		ENDIF
		
	END SUBROUTINE

!	=========================================================================================
	
	SUBROUTINE gen_TD(domain_file, mass_file, my_amass_TD, my_ityp_TD, PD, TD, periodic, ntyp, &
						amass1,	amass2, natc, itypc, mass_input, LPML, nr3, my_TD3, TD3_start)

		USE essentials

		IMPLICIT NONE
		INCLUDE 'mpif.h' ! MPI header file
		character(len = 256)			::	domain_file, mass_file, domain_prefix_id, mass_prefix_id, &
											format_string
		REAL							::	PD(3), TD(3)
		INTEGER, AlLOCATABLE			::	my_ityp_PD(:,:,:,:), my_ityp_TD(:,:,:,:)
		INTEGER							::	nr3, n1, n2, n3, natc, na, ntyp(2), my_n3, isub(4), &
											nentries
		logical							::	mass_input, periodic
		real(kind = RP), ALLOCATABLE	::	my_amass_PD(:,:,:,:), my_amass_TD(:,:,:,:)
		real(kind = RP)					::	amass1(ntyp(1)), amass2(ntyp(2))
		integer							::	LPML, summation
		INTEGER							::	itypc(natc)
		INTEGER(KIND = IP)				::	recv_values(world_size), displs(world_size)
		INTEGER							::	my_PD3, PD3_start, my_TD3, TD3_start, &
											everyones_PD3_start(world_size), rank, location, testing, &
											get_buffer(world_size, 2), put_buffer(world_size, 2)
		INTEGER							:: 	status(MPI_STATUS_SIZE), fh
		INTEGER(KIND = MPI_OFFSET_KIND)	:: 	offset
		INTEGER							:: 	type_size, counter
		INTEGER(KIND = IP)				::	my_natoms
		INTEGER							::	yplane, i
		INTEGER					 		:: 	win
		INTEGER							::	recvdispls(world_size), recvcounts(world_size), &
											sdispls(world_size), scounts(world_size)
		
		
!		There is a bug in the code where the MPI_READ only works if the number of processors set
!		during this run is the same as what was set when generating the domain and mass files.
!		I've not been able to locate the source of error and spent ton of time on it. 
!		Rohit

		CALL get_nr3(int(PD(3)), my_PD3, PD3_start)
		ALLOCATE(my_ityp_PD(natc,int(PD(1)),int(PD(2)),my_PD3), &
					my_amass_PD(natc,int(PD(1)),int(PD(2)),my_PD3))
		
		my_natoms = size(my_ityp_PD)
		
		IF (my_id .lt. 10) then
			format_string = "(a, a, I1, a)"
		ELSEIF (my_id .lt. 100) THEN
			format_string = "(a, a, I2, a)"
		ELSE 
			format_string = "(a, a, I3, a)"
		ENDIF

		write(domain_prefix_id, format_string) trim(domain_file), '_', my_id, '.dat'
		write(mass_prefix_id, format_string) trim(mass_file), '_', my_id, '.dat'
		
		open (unit  = 1, file = domain_prefix_id, form = 'unformatted')
		read (unit = 1) my_ityp_PD
		close(unit = 1)
		
		IF (mass_input) THEN
			open (unit  = 1, file = mass_prefix_id, form = 'unformatted')
			read (unit = 1) my_amass_PD
			close(unit = 1)
		ENDIF
							
		IF (io_node) WRITE(stdout, "(a,i5)") 'TD(3)=',int(TD(3))
		CALL get_nr3(int(TD(3)), my_TD3, TD3_start)
				
		CALL MPI_ALLGATHER(PD3_start, 1, MPI_INT, everyones_PD3_start, 1, MPI_INT, comm, ierr)
		
		ALLOCATE(my_ityp_TD(natc, int(TD(1)), int(TD(2)), (-2*nr3+1):(my_TD3 + 2*nr3)) , &
				 my_amass_TD(natc, int(TD(1)), int(TD(2)), (-2*nr3+1):(my_TD3 + 2*nr3)))
		
		get_buffer(:, 1) = huge(1)
		get_buffer(:, 2) = -1
		put_buffer(:, 1) = huge(1)
		put_buffer(:, 2) = -1
!		summation = 0
		
		IF (periodic) THEN !Um, what does it do if it's NOT periodic?  Looks like it's broken for that case
			DO n1 = 1, TD(1)
				DO n2 = 1, TD(2)
					DO my_n3 = -2*nr3+1, my_TD3 + 2*nr3
						n3 = TD3_start + my_n3
						IF ((n1.gt.(TD(1)/2.D0-PD(1)/2.D0)) .and. &
							(n1.le.(TD(1)/2.D0+PD(1)/2.D0)) .and. &
							(n2.gt.(TD(2)/2.D0-PD(2)/2.D0)) .and. &
							(n2.le.(TD(2)/2.D0+PD(2)/2.D0)) .and. &
							(n3.gt.(TD(3)/2.D0-PD(3)/2.D0)) .and. &
							(n3.le.(TD(3)/2.D0+PD(3)/2.D0))) THEN
!							If the atom is inside the Primary domain
							CALL get_rank(n3 - LPML, rank, location, everyones_PD3_start)
							location = location * int(PD(1)) * int(PD(2)) * natc
							IF (get_buffer(rank, 1) .gt. location) &
											get_buffer(rank, 1) = location
							IF (get_buffer(rank, 2) .lt. (location + int(PD(1)) * int(PD(2)) * natc)) &
											get_buffer(rank, 2) = location + int(PD(1)) * int(PD(2)) * natc
							
							IF (put_buffer(rank, 1) .gt. ((my_n3 - 1 + 2*nr3)* int(TD(1)) * int(TD(2)) *natc)) &
												put_buffer(rank, 1) = ((my_n3-1 + 2*nr3)* int(TD(1)) * int(TD(2)) *natc)
							IF (put_buffer(rank, 2) .lt. ((my_n3 + 2*nr3)* int(TD(1)) * int(TD(2)) *natc)) &
												put_buffer(rank, 2) = ((my_n3 + 2*nr3)* int(TD(1)) * int(TD(2)) *natc)

						ELSEIF (n3.lt.TD(3)/2) THEN
							my_ityp_TD(:, n1, n2, my_n3) = 1
							DO na = 1, natc
								my_amass_TD(na, n1, n2, my_n3) = amass1(itypc(na))
							END DO
						ELSEIF (n3.ge.TD(3)/2) THEN
							my_ityp_TD(:, n1, n2, my_n3) = 2
							DO na = 1, natc
								my_amass_TD(na, n1, n2, my_n3) = amass2(itypc(na))
							END DO
						ENDIF
					ENDDO
				ENDDO
			ENDDO
		ENDIF
		

!		WRITE (stdout, *) '--------------------------------------------------------'
!		WRITE (stdout, *) 'writing get buffer for ', my_id
!		WRITE (stdout, *) '--------------------------------------------------------'
		
!		DO i = 1, world_size
!			print *, get_buffer(i, 1), get_buffer(i, 2)
!			print *, put_buffer(i, 1), put_buffer(i, 2)
!		END DO


		DO rank = 1, world_size
			IF ((put_buffer(rank, 1) .eq. huge(1)) .or. &
				(put_buffer(rank, 2) .eq. -1)) THEN
				get_buffer(rank, 1) = 0
				get_buffer(rank, 2) = 0
				recvcounts(rank) = 0
				recvdispls(rank) = 0
			ELSE
				get_buffer(rank, 2) = get_buffer(rank, 2) - get_buffer(rank, 1)
				recvdispls(rank) = put_buffer(rank, 1)
				recvcounts(rank) = put_buffer(rank, 2) - put_buffer(rank, 1)
			END IF			
		END DO
		
!!		WRITE (stdout, *) '--------------------------------------------------------'
!!		WRITE (stdout, *) 'writing get buffer for ', my_id
!!		WRITE (stdout, *) '--------------------------------------------------------'
		
!!		DO i = 1, world_size
!!			print *, get_buffer(i, 1), get_buffer(i, 2)
!!			print *, recvdispls(i), recvcounts(i)
!!		END DO
		
		CALL MPI_ALLTOALL(get_buffer(:,1), 1, MPI_INT, sdispls, 1, MPI_INT, comm, ierr)
		CALL MPI_ALLTOALL(get_buffer(:,2), 1, MPI_INT, scounts, 1, MPI_INT, comm, ierr)
		
!		WRITE (stdout, *) '--------------------------------------------------------'
!		WRITE (stdout, *) 'writing get buffer for ', my_id + 1
!		WRITE (stdout, *) '--------------------------------------------------------'
		
!		DO i = 1, world_size
!			print *, 'send', i, sdispls(i), scounts(i)
!			print *, 'recv', i, recvdispls(i), recvcounts(i)
!		END DO
		
		CALL MPI_ALLTOALLV(	my_ityp_PD, scounts, sdispls, MPI_INT, &
							my_ityp_TD, recvcounts, recvdispls, MPI_INT, comm, ierr)
							
		IF (mass_input) THEN
			CALL MPI_ALLTOALLV( my_amass_PD, scounts, sdispls, mp_real, &
							my_amass_TD, recvcounts, recvdispls, mp_real, comm, ierr)
		ELSE
			DO n1 = 1, TD(1)
				DO n2 = 1, TD(2)
					DO my_n3 = -2*nr3+1, my_TD3 + 2*nr3
						DO na = 1, natc
							IF (my_ityp_TD(na, n1, n2, my_n3) .eq. 1) THEN
								my_amass_TD(na, n1, n2, my_n3) = amass1(itypc(na))
							ELSE IF (my_ityp_TD(na, n1, n2, my_n3) .eq. 2) THEN
								my_amass_TD(na, n1, n2, my_n3) = amass2(itypc(na))
							ENDIF
						END DO
					END DO
				END DO
			END DO
		ENDIF
!!		print *, size(my_ityp_TD) - summation
		
		yplane = int(TD(2)/2.D0)
		
		WRITE (stdout, *) '--------------------------------------------------------'
		WRITE (stdout, *) 'Cross-sectional image of PD (check if this is correct)', my_id
		WRITE (stdout, *) '--------------------------------------------------------'

			!
			! I only want to show what atoms are in the central domain...not which ones that connect with!
!!					DO n1 = 1, TD(1)
!!		DO n3 = -2*nr3 + 1, my_TD3 + 2*nr3
!!				IF (n3.eq.my_TD3 + 2*nr3) THEN
!!					write(stdout, fmt = '(I2)') my_ityp_TD(1, n1,yplane,n3)
!!				ELSE
!!					write(stdout, fmt = '(I2)', advance = 'no') &
!!											my_ityp_TD(1, n1,yplane,n3)
!!				ENDIF
!!			ENDDO
!!		ENDDO
		
		DO n1 = 1, TD(1)
			!
			! I only want to show what atoms are in the central domain...not which ones that connect with!
			DO n3 = 1, my_TD3
				IF (n3.eq.my_TD3) THEN
					! only advance to next line after writing all other n3 entries
					write(stdout, fmt = '(I2)') my_ityp_TD(1, n1,yplane,n3)
				ELSE
					write(stdout, fmt = '(I2)', advance = 'no') &
											my_ityp_TD(1, n1,yplane,n3)
				ENDIF
			ENDDO
		ENDDO

!		print *, 'Hello'
		
!		RETURN
				
	END SUBROUTINE

!	==========================================================================================
	
	SUBROUTINE print_domain(PD, TD, ityp_PD, ityp_TD)

		USE mp_module
		USE kinds
		use, intrinsic 				:: 	iso_fortran_env, only : stdin=>input_unit, &
												  stdout=>output_unit, &
												  stderr=>error_unit
		IMPLICIT NONE
		INTEGER						:: 	yplane, n1, n2, n3
		REAL						::	PD(3), TD(3)
		INTEGER						::	ityp_PD(int(PD(1)),int(PD(2)),int(PD(3)))
		INTEGER						::	ityp_TD(int(TD(1)),int(TD(2)),int(TD(3)))			
	
		yplane = int(PD(2)/2.D0)
		
		IF (io_node) THEN
			WRITE (stdout, *) '--------------------------------------------------------'
			WRITE (stdout, *) 'Cross-sectional image of PD (check if this is correct)'
			WRITE (stdout, *) '--------------------------------------------------------'
			DO n1 = 1, PD(1)
				DO n3 = 1, PD(3)
					IF (n3.eq.PD(3)) THEN
						write(stdout, fmt = '(I2)') ityp_PD(n1,yplane,n3)
					ELSE
						write(stdout, fmt = '(I2)', advance = 'no') &
												ityp_PD(n1,yplane,n3)
					ENDIF
				ENDDO
			ENDDO
		ENDIF
		 
		yplane = int(TD(2)/2.D0)
		 
		IF (io_node) THEN
			WRITE (stdout, *) '	' 
			WRITE (stdout, *) 'Cross-sectional image of TD (check if this is correct)'
			WRITE (stdout, *) '--------------------------------------------------------'
			DO n1 = 1, TD(1)
				DO n2 = 1, TD(3)
					IF (n2.eq.TD(3)) THEN
						write(stdout, fmt = '(I2)') ityp_TD(n1,yplane,n2)
					ELSE
						write(stdout, fmt = '(I2)', advance = 'no') &
													ityp_TD(n1,yplane,n2)
					ENDIF
				ENDDO
			ENDDO
		ENDIF
		
		CALL MPI_ABORT(comm, ierr)
		
	END SUBROUTINE
	
	!	====================================================================================================
	
	SUBROUTINE get_nr3(nr3, my_nr3, nr3_start)
	
		IMPLICIT NONE
		INCLUDE 'mpif.h' ! MPI header file
		
		INTEGER				:: nr3, my_nr3, everyones_nr3(world_size), nr3_start, rem, i, cumulative
		
	
		my_nr3 = nr3/world_size
		
		rem = MOD(nr3,world_size)
			
		IF (my_id.gt.(world_size-rem-1)) THEN
			my_nr3 = my_nr3+1
		ENDIF
		WRITE(stdout, '(a11, i3, a15, i5)') 'Processor ',my_id,' says my_nr3 = ', my_nr3
		
		
		CALL MPI_ALLGATHER(my_nr3, 1, MPI_INT, everyones_nr3, 1, MPI_INT, comm, ierr)
		
		cumulative = 0
		DO i = 0, world_size-1
			IF (i .eq. my_id) THEN
				nr3_start = cumulative
			END IF
			cumulative = cumulative + everyones_nr3(i+1)
		END DO
		
		IF (cumulative .ne. nr3) THEN
			IF (io_node) WRITE(stdout, '(a)') 'Problems with splitting up the domain'
			CALL MPI_ABORT(comm, ierr)
		END IF
		
	END SUBROUTINE
	
	!	====================================================================================================
	
	SUBROUTINE get_rank(n3, rank, location, everyones_PD3_start)
	! what is the goal of this subroutine?  what is 'location'?  location of what?
	
		IMPLICIT NONE
		INTEGER			:: n3, rank, location, everyones_PD3_start(world_size), i
		
		DO i = 1, world_size-1
			IF (n3 .gt. everyones_PD3_start(i+1)) THEN
				CYCLE
			ELSE
				rank = i
				location = n3 - everyones_PD3_start(i) - 1
				RETURN
			END IF
		END DO
		
		rank = world_size
		location = n3 - everyones_PD3_start(world_size) - 1
	END SUBROUTINE
	
END MODULE preprocessing_module



