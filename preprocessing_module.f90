MODULE preprocessing_module
	
	USE kinds
	USE mp_module
	IMPLICIT NONE
	CONTAINS
	
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
		
		call get_qpoint(asr, na_ifc, fd, has_zstar, q, alat1, alat2, at1, at2, &
							nat, ntyp, nr1, nr2, nr3, ibrav, mode, ityp1, ityp2, epsil, zeu1, &
							zeu2, vg, omega1, omega2, frc1, frc2, w2, f_of_q1, f_of_q2, amass1, amass2, &
							tau1, tau2)
							
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
		REAL(KIND = RP) 					:: 	at_conv(3,3)
		REAL								::	atc(3,3)
		INTEGER 							:: 	natsc, natc
		REAL(KIND = RP), ALLOCATABLE 		:: 	tausc(:,:), r_cell(:,:), tauc(:,:)
		INTEGER, ALLOCATABLE 				:: 	itypsc(:), itypc(:)
		COMPLEX(KIND = CP), ALLOCATABLE 	:: 	zsc(:), zc(:)
		INTEGER, ALLOCATABLE 				:: 	ib_vec1(:,:,:,:,:), ib_vec2(:,:,:,:,:), &
												na_vec(:)
		REAL(KIND = RP)						::	tau1(3, nat(1))
		INTEGER								::	nr1, nr2, nr3
		
		IF (crystal_coordinates) THEN
			at_conv = at1 	! If working in crystal coordinates then supercell 
							! is the primitive cell
		ELSE
			at_conv(:,:) = 0.D0
			DO i = 1, 3
				at_conv(i,i) = 1.0_RP
			ENDDO
		ENDIF
		
	
		CALL Supercell(	at_conv, at1, tau1, ityp1, nat(1), tausc, itypsc, &
						natsc, z, zsc, r_cell, na_vec, ib_vec1, nr1, nr2, nr3 )
		!	___c variables listed below are current variables and are chosen based on 
		!	the choice of coordinate system you want to work in.
		natc = natsc
		if (.not. allocated(tauc)) ALLOCATE (tauc(3,natc), itypc(natc), zc(3*natc))
		if (.not. allocated(ib_vec2)) ALLOCATE(ib_vec2(nat(2), -2*nr1:2*nr1, &
												-2*nr2:2*nr2, -2*nr3:2*nr3, natc))
		ib_vec2 = ib_vec1
		tauc = tausc
		itypc = itypsc
		atc = at_conv
		zc = zsc
	
		IF (real(zc(1)).LT.0.0) THEN
			zc = CMPLX(-1.0, 0.0, KIND = CP)*zc(:)
		ENDIF
		
	END SUBROUTINE

!	=========================================================================================
	
	SUBROUTINE gen_TD(domain_file, mass_file, amass_TD, ityp_TD, PD, TD, periodic, ntyp, &
						amass1,	amass2, natc, itypc, mass_input, LPML)


		IMPLICIT NONE
		INCLUDE 'mpif.h' ! MPI header file
		character(len = 256)			::	domain_file, mass_file
		REAL							::	PD(3), TD(3)
		INTEGER, ALLOCATABLE			:: 	ityp_PD(:,:,:,:), ityp_TD(:,:,:,:)
		INTEGER							::	nr1, nr2, nr3, n1, n2, n3, natc, na, ntyp(2)
		logical							::	mass_input, periodic
		real(kind = RP), ALLOCATABLE	::	amass_PD(:,:,:,:), amass_TD(:,:,:,:)
		real(kind = RP)					::	amass1(ntyp(1)), amass2(ntyp(2))
		integer							::	LPML
		INTEGER							::	itypc(natc)
				
		ALLOCATE(ityp_PD(int(PD(1)), int(PD(2)), int(PD(3)), natc))
		if (mass_input) ALLOCATE (amass_PD(int(PD(1)), int(PD(2)), int(PD(3)), natc))
		
		IF (io_node) THEN
			open (unit  = 1, file = domain_file, form = 'unformatted')
			read (unit = 1) ityp_PD
			close(unit = 1)
		ENDIF
		
		
		CALL MPI_BCAST(ityp_PD, size(ityp_PD), MPI_INT, root_process, comm, ierr)
		
		print *, mass_input
	
		IF (mass_input) THEN
			IF (io_node) THEN
				open (unit  = 1, file = mass_file, form = 'unformatted')
				read (unit = 1) amass_PD
				close(unit = 1)
			ENDIF
			CALL MPI_BCAST(amass_PD, size(amass_PD), mp_real, root_process, comm, ierr)	
		ENDIF
		
		
	!	----------------------------------------------------------------------------------------
	!	Generating the Total Domain
	!	This section adds a layer of PML of length LPML on the edges of the 
	!	primary domain.
	!	The Total domain defines the type of atom for every atom inside the simulation
	!	cell
	
		ALLOCATE(ityp_TD(int(TD(1)), int(TD(2)), int(TD(3)), natc))
		ALLOCATE(amass_TD(int(TD(1)), int(TD(2)), int(TD(3)), natc))
	
		print *, int(TD(1)), int(TD(2)), int(TD(3)), natc
	
		IF (periodic) THEN
			DO n1 = 1, TD(1)
				DO n2 = 1, TD(2)
					DO n3 = 1, TD(3)
						DO na = 1, natc
							IF ((n1.gt.(TD(1)/2.D0-PD(1)/2.D0)) .and. &
								(n1.le.(TD(1)/2.D0+PD(1)/2.D0)) .and. &
								(n2.gt.(TD(2)/2.D0-PD(2)/2.D0)) .and. &
								(n2.le.(TD(2)/2.D0+PD(2)/2.D0)) .and. &
								(n3.gt.(TD(3)/2.D0-PD(3)/2.D0)) .and. &
								(n3.le.(TD(3)/2.D0+PD(3)/2.D0))) THEN
	!							If the atom is inside the Primary domain
								ityp_TD(n1,n2,n3,na)= ityp_PD(n1, n2, n3-LPML, na)
								IF (mass_input) THEN
									amass_TD(n1,n2,n3,na) = amass_PD(n1,n2,n3-LPML, na)
								ELSE
									IF (ityp_TD(n1,n2,n3,na).eq.1) THEN
										amass_TD(n1,n2,n3,na) = amass1(itypc(na))
									ELSEIF (ityp_TD(n1,n2,n3,na).eq.2) THEN
										amass_TD(n1,n2,n3,na) = amass1(itypc(na))
									ENDIF
								ENDIF
							ELSEIF (n3.lt.TD(3)/2) THEN
								ityp_TD(n1,n2,n3,na) = 1
								amass_TD(n1,n2,n3,na) = amass1(itypc(na))
							ELSEIF (n3.ge.TD(3)/2) THEN
								ityp_TD(n1,n2,n3,na) = 2
								amass_TD(n1,n2,n3,na) = amass2(itypc(na))
							ENDIF
						ENDDO
					ENDDO
				ENDDO
			ENDDO
		ELSE
			DO n1 = 1, TD(1)
				DO n2 = 1, TD(2)
					DO n3 = 1, TD(3)
						DO na = 1, natc
							IF ((n1.gt.(TD(1)/2.D0-PD(1)/2.D0)) .and. &
								(n1.lt.(TD(1)/2.D0+PD(1)/2.D0)) .and. &
								(n2.gt.(TD(2)/2.D0-PD(2)/2.D0)) .and. &
								(n2.lt.(TD(2)/2.D0+PD(2)/2.D0)) .and. &
								(n3.gt.(TD(3)/2.D0-PD(3)/2.D0)) .and. &
								(n3.lt.(TD(3)/2.D0+PD(3)/2.D0))) THEN
	!							If the atom is inside the Primary domain
								ityp_TD(n1,n2,n3,na)= &
									ityp_PD(n1-LPML, n2-LPML, n3-LPML, na)
								IF (mass_input) THEN
									amass_TD(n1,n2,n3,na) = amass_PD(n1-LPML,n2-LPML,&
																		n3-LPML, na)
								ELSE
									IF (ityp_TD(n1,n2,n3,na).eq.1) THEN
										amass_TD(n1,n2,n3,na) = amass1(itypc(na))
									ELSEIF (ityp_TD(n1,n2,n3,na).eq.2) THEN
										amass_TD(n1,n2,n3,na) = amass2(itypc(na))
									ENDIF
								ENDIF
							ELSE
								ityp_TD(n1,n2,n3,na) = 1
								amass_TD(n1,n2,n3,na) = amass1(itypc(na))
							ENDIF
						ENDDO
					ENDDO
				ENDDO
			ENDDO
		ENDIF
		
		call print_domain(PD, TD, ityp_PD(:,:,:,1), ityp_TD(:,:,:,2))
		
		
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
	END SUBROUTINE
	
END MODULE preprocessing_module
