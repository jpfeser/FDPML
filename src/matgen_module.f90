MODULE matgen_module

	CONTAINS
	
!	=============================================================================================

	SUBROUTINE gen_uinc(natc, nSub, my_natoms, my_nrows, zclist, wavetype, my_uinc, atoms_start, &
						atc, r_cell, LPML, qlist, nq)
		
		USE kinds
		USE essentials
		USE constants
		USE mp_module
		IMPLICIT NONE
		
		INTEGER								::	counter, p, i, natc, LPML, nq
		INTEGER(KIND = IP)					::	my_natoms, my_nrows
		COMPLEX(KIND = CP), ALLOCATABLE		::	my_uinc(:)
		INTEGER(KIND = IP)					::	atom_tuple(4)
		COMPLEX(KIND = CP)					::	zc(3*natc), zclist(3*natc,nq)
		REAL								::	arg
		CHARACTER(len = 256)				::	wavetype
		INTEGER(KIND = IP)					::	nSub(4)
		INTEGER(KIND = IP)					::	atoms_start(world_size), n, n1, n2, n3, na
		REAL(KIND = RP)						::	r(3), q(3), qlist(3, nq)
		REAL(KIND = RP)						::	r_cell(3, natc)
		REAL								::	atc(3,3)
		INTEGER								::	qpoint
		
		
		
		
		ALLOCATE(my_uinc(my_nrows))
		
		my_uinc(:) = 0.0_RP
		
		DO qpoint = 1, nq
			q = qlist(:,qpoint)
			zc = zclist(:, qpoint)
			counter = 0
			DO p = 1,my_natoms
				n = p + atoms_start(my_id+1)
				atom_tuple = ind2sub(n, nSub) 	!	convert a 1-D indexial scheme to a 3-D.
												!	(Find exact location of nth atom in 3-D
												!	space) 
				na = atom_tuple(1)
				n1 = atom_tuple(2)
				n2 = atom_tuple(3)
				n3 = atom_tuple(4)
				! Location from origin
				r = (n1-1)*atc(:,1) + (n2-1)*atc(:,2) + (n3-1)*atc(:,3) 
				r = r + r_cell(:,na)
				arg = tpi*(q(1)*r(1) + q(2)*r(2) + q(3)*r(3))
				DO i= 1, 3
				counter = counter + 1 
				! wave with amplitude = 1
					IF (wavetype .eq. 'full') THEN
						my_uinc(counter)= my_uinc(counter) + zc(3*na-(3-i))* &
											CMPLX(COS(arg),-SIN(arg),KIND=CP) / nq**2
					ELSEIF (wavetype .eq. 'half') THEN
						IF (n3.lt.LPML) THEN
							my_uinc(counter)= my_uinc(counter) + zc(3*na-(3-i))* &
												CMPLX(COS(arg),-SIN(arg),KIND=CP) / nq**2
						ELSE
							my_uinc(counter) = my_uinc(counter) + CMPLX(0.0_RP, 0.0_RP, KIND = CP)
						ENDIF
					ENDIF
				ENDDO
			ENDDO
		ENDDO
		
	END SUBROUTINE

!	=========================================================================================

	SUBROUTINE gen_sig(sig, my_natoms, sigmamax, LPML, periodic, atc, atom_tuple, &
						TD, PD, atoms_start, my_ityp_TD, tauc, natc, my_TD3, nr3, TD3_start)
!	generate the damping coefficients for every atom inside the simulation domain						

	
		USE kinds
		USE mp_module
		USE essentials
		IMPLICIT NONE
		INTEGER								::	counter, p, i1, i2, i3, ia, ityp, natc, my_TD3, nr3, &
												TD3_start, i3max
		INTEGER(KIND = IP)					::	my_natoms
		REAL(KIND = RP)						::	sigmamax
		LOGICAL								::	periodic
		INTEGER(KIND = IP)					::	nSub(4), n
		REAL								::	atc(3,3)
		REAL(KIND = RP), ALLOCATABLE		::	sig(:,:)
		INTEGER(KIND = IP)					::	atom_tuple(4), atoms_start(world_size)
		REAL(KIND = RP)						::	r(3)
		REAL								::	TD(3), PD(3)
		INTEGER								::	my_ityp_TD(natc, int(TD(1)), &
															int(TD(2)), (-2*nr3+1):(my_TD3 + 2*nr3)), &
															LPML
		REAL(KIND = RP)						::	tauc(3,natc)
		
		
		
		
		nSub = (/ natc, int(TD(1)), int(TD(2)), int(TD(3)) /)
		
		ALLOCATE(sig(my_natoms, 3))
		sig(:,:) = 0.0_RP
		
		DO p = 1, my_natoms
			n = p + atoms_start(my_id+1)
			atom_tuple = ind2sub(n, nSub) !need to figure out the global coordinates
			ia = atom_tuple(1) 
			i1 = atom_tuple(2)
			i2 = atom_tuple(3)
			i3 = atom_tuple(4) - TD3_start
			
			ityp = my_ityp_TD(ia, i1, i2, i3)
			r = (i1-1)*atc(:,1) + (i2-1)*atc(:,2) + (i3 + TD3_start -1)*atc(:,3) + tauc(:,ia)
			IF (r(1).lt.((TD(1)-1.D0)/2.D0-(PD(1)-1.D0)/2.D0)) THEN
		        sig(p,1)= abs(r(1)-((TD(1)-1.D0)/2.D0-(PD(1)-1.D0)/2.D0))**2
			ELSEIF (r(1).gt.((TD(1)-1.D0)/2.D0+(PD(1)-1.D0)/2.D0)) THEN
				sig(p,1)= abs(r(1)-((TD(1)-1.D0)/2.D0+(PD(1)-1.D0)/2.D0))**2
			ENDIF
			IF (r(2).lt.((TD(2)-1.D0)/2.D0-(PD(2)-1.D0)/2.D0)) THEN
				sig(p,2)= abs(r(2)-((TD(2)-1.D0)/2.D0-(PD(2)-1.D0)/2.D0))**2
			ELSEIF (r(2).gt.((TD(2)-1.D0)/2.D0+(PD(2)-1.D0)/2.D0)) THEN
				sig(p,2)= abs(r(2)-((TD(2)-1.D0)/2.D0+(PD(2)-1.D0)/2.D0))**2
			ENDIF
			IF (r(3).lt.((TD(3)-1.D0)/2.D0-(PD(3)-1.D0)/2.D0)) THEN
				sig(p,3)= abs(r(3)-((TD(3)-1.D0)/2.D0-(PD(3)-1.D0)/2.D0))**2
			ELSEIF (r(3).gt.((TD(3)-1.D0)/2.D0+(PD(3)-1.D0)/2.D0)) THEN
				sig(p,3)= abs(r(3)-((TD(3)-1.D0)/2.D0+(PD(3)-1.D0)/2.D0))**2
			ENDIF
		ENDDO
		
		sig = sigmamax*sig/(LPML**2)
		
		IF (periodic) THEN
			sig(:,1) = sig(:,3)
			sig(:,2) = sig(:,3)
		ELSE
			sig(:,1) = sig(:,1) + sig(:,2) + sig(:,3)
			sig(:,2) = sig(:,1)
			sig(:,3) = sig(:,1)
		ENDIF
	END SUBROUTINE gen_sig

!	========================================================================================

	SUBROUTINE gen_Amat(nr1, nr2, nr3, nat, at1, at2, atc, nrwsx, TD, ityp_TD, na_vec, natc, ib_vec1, &
						ib_vec2, r_cell, tau1, tau2, sig, my_nrows, my_natoms, atoms_start, amass_TD, &
						frc1, frc2, f_of_q1, f_of_q2, w2, mode, fd, ilist, jlist, Alist, counter1, nSub, &
						borderlogic, counter2, wscache, atws, iSub, r_ws, my_nnz)
	
		USE kinds
		USE mp_module
		USE essentials
		use, intrinsic :: iso_fortran_env, only : stdin=>input_unit, &
											  stdout=>output_unit, &
											  stderr=>error_unit
		USE ws
		USE constants
		IMPLICIT NONE
		INTEGER									::	nr1, nr2, nr3, nat(2), natc, na, m1, &
													m2, m3
		REAL(KIND = RP), ALLOCATABLE			::	wscache(:,:,:,:,:)
		REAL(KIND = RP)							:: 	atws(3,3), at1(3,3), at2(3,3), weight
		REAL									::	atc(3,3)
		INTEGER									::	nrwsx
		INTEGER									::	nrws
		REAL(RP) 								::	rws(0:3,nrwsx)
		INTEGER									::	nb, n1, n2, n3, i, p, ipol, jpol, ib
		INTEGER(KIND=IP)						::	n, atom_tuple(4), i1, i2, i3, ia, &
													nSub(4), Location, iSub(4)
		REAL									::	TD(3)
		INTEGER									::	ityp_TD(int(TD(1)), int(TD(2)),&
														int(TD(3)), natc), ityp
		INTEGER									::	na_vec(natc), &
													ib_vec1(nat(1), -2*nr1:2*nr1, &
														-2*nr2:2*nr2, -2*nr3:2*nr3, natc), &
													ib_vec2(nat(2), -2*nr1:2*nr1, &
														-2*nr2:2*nr2, -2*nr3:2*nr3, natc)
		REAL(KIND = RP)							::	r(3), r_ws(3)
		LOGICAL									:: 	periodic
		INTEGER(KIND = IP)						:: 	counter1, counter2, counter
		REAL(KIND = RP)							::	r_cell(3, natc), tau1(3, nat(1)), &
													tau2(3, nat(2))
		INTEGER(KIND = IP)						::	my_natoms, my_nrows, atoms_start(world_size), &
													my_nnz
		REAL(KIND = RP)							::	sig(my_natoms,3)
		INTEGER(KIND = IP), ALLOCATABLE			::	borderlogic(:)
		REAL(KIND = RP)							::	mass1, mass2, amass_TD(int(TD(1)), &
														int(TD(2)), int(TD(3)), natc), &
													frc1(nr1,nr2,nr3,3,3,nat(1),nat(1)), &
													frc2(nr1,nr2,nr3,3,3,nat(2),nat(2)), &
													IFC(3,3), f_of_q1(3,3,nat(1),nat(2)), &
													f_of_q2(3,3,nat(2), nat(2))
		REAL(KIND = RP)							::	w2(3*nat(1))
		INTEGER									:: 	mode
		INTEGER(KIND = IP), ALLOCATABLE			::	ilist(:), jlist(:)
		COMPLEX(KIND = CP), ALLOCATABLE			::	Alist(:)
		LOGICAL									::	fd
		
													
		
		

		
		
		
		
		
		!	=======================================================================================
	
	!	Set up A-matrix
	
		IF (io_node) write (stdout, *) '	'
		IF (io_node) write (stdout, *) '-----------------------------------------------------------------------'
		IF (io_node) write (stdout, *) 'Setting up A matrix'
	
	
	!	----------------------------------------------------------------------------------------
	!	Copied as is from qunatum espresso
	!	Please refer FAQ page on qunatum-espresso (phonon section, question 7.7)
	!	link: http://www.quantum-espresso.org/resources/faq/phonons
	
		print *, nrwsx
	
		ALLOCATE(wscache(-2*nr3:2*nr3, -2*nr2:2*nr2, -2*nr1:2*nr1, nat(1), nat(1)))
	
		atws(:,1) = at1(:,1)*DBLE(nr1)
		atws(:,2) = at1(:,2)*DBLE(nr2)
		atws(:,3) = at1(:,3)*DBLE(nr3)
		! initialize WS r-vectors
		CALL wsinit(rws,nrwsx,nrws,atws)
	
		counter = 0
		DO na=1, nat(1)
			DO nb=1, nat(1)
			!
				DO n1=-2*nr1,2*nr1
					DO n2=-2*nr2,2*nr2
						DO n3=-2*nr3,2*nr3
							DO i=1, 3
								r(i) = n1*at1(i,1)+n2*at1(i,2)+n3*at1(i,3)
								r_ws(i) = r(i) + tau1(i,na)-tau1(i,nb)
								if (fd) r_ws(i) = r(i) + tau1(i,nb)-tau1(i,na)
							END DO
							wscache(n3,n2,n1,nb,na) = wsweight(r_ws,rws,nrws)
							IF (wscache(n3,n2,n1,nb,na).gt.0.D0) THEN
								counter = counter + 1 
							ENDIF	
						ENDDO
					ENDDO
				ENDDO
			ENDDO
		ENDDO
		
		CALL MPI_BCAST(wscache, ((4*nr1+1)*(4*nr2+1)*(4*nr3+1)*nat(1)*nat(2)), mp_real, root_process, comm, ierr)
	!	----------------------------------------------------------------------------------------
	
	!	Counter loop to find how much memory is needed for Alist, ilist and jlist
		counter1 = 0	!	counter for ilist, jlist and Alist
		counter2 = 0	!	counter for atoms connected ouside the domain
	!	CALL cpu_time(start)
		DO p = 1, my_natoms
			n = p + atoms_start(my_id+1)
			atom_tuple = ind2sub(n, nSub)
			ia = atom_tuple(1)
			i1 = atom_tuple(2)
			i2 = atom_tuple(3)
			i3 = atom_tuple(4)
			ityp = ityp_TD(i1,i2,i3,ia)
			na = na_vec(ia)
			IF (ityp.eq.1) THEN
				DO n1 = -2*nr1,2*nr1
					DO n2 = -2*nr2,2*nr2
						DO n3 = -2*nr3,2*nr3
							DO nb = 1, nat(1)
								weight = wscache(n3,n2,n1,nb,na)
								IF (weight.gt.0) THEN
	!								Please refer QE link provide above wsweights to 
	!								understand m1, m2, m3
									m1 = MOD(n1+1,nr1)
									IF(m1.LE.0) m1=m1+nr1
									m2 = MOD(n2+1,nr2)
									IF(m2.LE.0) m2=m2+nr2
									m3 = MOD(n3+1,nr3)
									IF(m3.LE.0) m3=m3+nr3
									r = (i1-1)*atc(:,1) + (i2-1)*atc(:,2) + &
										(i3-1)*atc(:,3)
									r = r + r_cell(:,ia) + n1*at1(:,1) + n2*at1(:,2) + &
										n3*at1(:,3)+ tau1(:,nb)
									r = floor(r)
									ib = ib_vec1(nb,n3,n2,n1,ia) !	Defined in supercell 
																 !	routine
									IF (ib.eq.-1) THEN
										WRITE (stdout, '(a)') 'ERROR : wrong ib'
										STOP
									ENDIF
	!								Apply periodic boundaries
									IF (periodic) THEN
										r(1)= MOD(r(1),(TD(1)-1))
		                                IF (r(1).lt.0) THEN
											r(1)=r(1)+TD(1)-1
										ENDIF
										r(2)= MOD(r(2),(TD(2)-1))
										IF (r(2).lt.0) THEN
											r(2)=r(2)+TD(2)-1
										ENDIF
									ENDIF
	!								If the connecting atom is outside TD, then skip loop								
									IF ((r(1).lt.0) .or. (r(2).lt.0) .or. &
										(r(3).lt.0) .or. (r(1).gt.(TD(1)-1)) .or. &
										(r(2).gt.(TD(2)-1)) .or. (r(3).gt.(TD(3)-1))) THEN
										DO i = 1, 3
											counter2 = counter2+1
										ENDDO
										CYCLE
									ENDIF
									DO ipol = 1,3
										DO jpol = 1,3
											counter1 = counter1 + 1
											IF (counter1.eq.huge(counter1)) THEN
												WRITE(stdout, '(a, a, I4)') &
												'Size of Amat is larger',& 
												'than that is allowed by int kind', IP
												STOP
											ENDIF
										ENDDO
									ENDDO
								ENDIF
							ENDDO
						ENDDO
					ENDDO
				ENDDO
			ELSEIF (ityp.eq.2) THEN
	!		This loop is identical to loop above ityp == 1. It only accounts 
	!		for a change in lattice cell type for material 2
				DO n1 = -2*nr1,2*nr1
					DO n2 = -2*nr2,2*nr2
						DO n3 = -2*nr3,2*nr3
							DO nb = 1, nat(2)
								weight = wscache(n3,n2,n1,nb,na)
								IF (weight.gt.0) THEN
									m1 = MOD(n1+1,nr1)
									IF(m1.LE.0) m1=m1+nr1
									m2 = MOD(n2+1,nr2)
									IF(m2.LE.0) m2=m2+nr2
									m3 = MOD(n3+1,nr3)
									IF(m3.LE.0) m3=m3+nr3
									r = (i1-1)*atc(:,1) + (i2-1)*atc(:,2) + &
										(i3-1)*atc(:,3)
									r = r + r_cell(:,ia) + n1*at2(:,1) + n2*at2(:,2) + &
										n3*at2(:,3)+ tau2(:,nb)
									r = floor(r)
									ib = ib_vec2(nb,n3,n2,n1,ia)
									IF (ib.eq.-1) THEN
										WRITE (stdout, '(a)') 'ERROR : wrong ib'
										STOP
									ENDIF
									IF (periodic) THEN
										r(1)= MOD(r(1),(TD(1)-1))
		                                IF (r(1).lt.0) THEN
											r(1)=r(1)+TD(1)-1
										ENDIF
										r(2)= MOD(r(2),(TD(2)-1))
										IF (r(2).lt.0) THEN
											r(2)=r(2)+TD(2)-1
										ENDIF
									ENDIF
									IF ((r(1).lt.0) .or. (r(2).lt.0) .or. &
										(r(3).lt.0) .or. (r(1).gt.(TD(1)-1)) .or. &
										(r(2).gt.(TD(2)-1)) .or. (r(3).gt.(TD(3)-1))) THEN
										DO i = 1, 3
											counter2 = counter2+1
										ENDDO
										CYCLE
									ENDIF
									DO ipol = 1,3
										DO jpol = 1,3
											counter1 = counter1 + 1
											IF (counter1.eq.huge(counter1)) THEN
												WRITE(stdout, '(a, a, I4)') &
												'Size of Amat is larger',& 
												'than that is allowed by int kind', IP
												STOP
											ENDIF
										ENDDO
									ENDDO
								ENDIF
							ENDDO
						ENDDO
					ENDDO
				ENDDO
			ENDIF	
		ENDDO
		
		DO p = 1, my_natoms
			DO ipol = 1,3
				IF (sig(p,ipol).eq.0.D0) THEN
					CYCLE
				ELSE
					counter1 = counter1 + 1
					IF (counter1.eq.huge(counter1)) THEN
						WRITE(stdout, '(a, a, I4)') &
							'Size of Amat is larger',& 
							'than that is allowed by int kind', IP
						STOP
					ENDIF
				ENDIF
			ENDDO
		ENDDO
	!	----------------------------------------------------------------------------------------
	
		my_nnz = counter1+my_nrows
		ALLOCATE(	ilist(counter1+my_nrows), jlist(counter1+my_nrows), &
					Alist(counter1+my_nrows)	) ! my_nrows accounts for -m*w**2 
												  ! subtraction from every row of A-mat
												  ! Amat = -(m*w**2+K) equn A3 on 
												  !	PHYSICAL REVIEW B 95, 125434 (2017)
	
		WRITE (stdout, '(a, I)') 'nnz = ', counter1+my_nrows
	
		ALLOCATE(borderlogic(counter2))		  ! Defines every atom connected outside
												  ! TD
		
		counter1 = 0
		counter2 = 0
		
!		my_BW = 0
	
		DO p = 1, my_natoms
			n = p + atoms_start(my_id+1)
			atom_tuple = ind2sub(n, nSub)
			ia = atom_tuple(1)
			i1 = atom_tuple(2)
			i2 = atom_tuple(3)
			i3 = atom_tuple(4)
			ityp = ityp_TD(i1,i2,i3,ia)
			na = na_vec(ia)
			IF (ityp.eq.1) THEN
				DO n1 = -2*nr1,2*nr1
					DO n2 = -2*nr2,2*nr2
						DO n3 = -2*nr3,2*nr3
							DO nb = 1, nat(1)
								weight = wscache(n3,n2,n1,nb,na)
								IF (weight.gt.0) THEN
									m1 = MOD(n1+1,nr1)
									IF(m1.LE.0) m1=m1+nr1
									m2 = MOD(n2+1,nr2)
									IF(m2.LE.0) m2=m2+nr2
									m3 = MOD(n3+1,nr3)
									IF(m3.LE.0) m3=m3+nr3
									r = (i1-1)*atc(:,1) + (i2-1)*atc(:,2) + &
										(i3-1)*atc(:,3)
									r = r + r_cell(:,ia) + n1*at1(:,1) + n2*at1(:,2) + &
										n3*at1(:,3)+ tau1(:,nb)
									r = floor(r)
									ib = ib_vec1(nb,n3,n2,n1,ia)
									IF (periodic) THEN
										r(1)= MOD(r(1),(TD(1)-1))
		                                IF (r(1).lt.0) THEN
											r(1)=r(1)+TD(1)-1
										ENDIF
										r(2)= MOD(r(2),(TD(2)-1))
										IF (r(2).lt.0) THEN
											r(2)=r(2)+TD(2)-1
										ENDIF
									ENDIF
									IF ((r(1).lt.0) .or. (r(2).lt.0) .or. &
										(r(3).lt.0) .or. (r(1).gt.(TD(1)-1)) .or. &
										(r(2).gt.(TD(2)-1)) .or. (r(3).gt.(TD(3)-1))) THEN
	!								IF the atom is connected outside the domain then
	!								record it.
										DO i = 1,3
											counter2 = counter2 + 1
											borderlogic(counter2) = 3*p - (3-i)
										ENDDO
										CYCLE
									ENDIF
									iSub = (/ ib, int(r(1))+1, int(r(2))+1, &
											int(r(3))+1 /)
											
	!								Location defines the location in uinc for 
	!								(n1, n2, n3, nb) atom. 
	!								Location essentially defines the column number in A
	
									Location = sub2ind(iSub, nSub) ! transform 3-D space 
																   ! to 1-D
	
									mass1 = amass_TD(i1,i2,i3,ia)
									mass2 = amass_TD(int(r(1))+1,int(r(2))+1,&
													 int(r(3))+1,ib)
	!								mass1 = amass1(ityp2(na))
	!								IF (ityp_TD(int(r(1))+1, int(r(2))+1, &
	!											int(r(3))+1, ib).eq.2) THEN
	!									mass2 = amass2(ityp2(nb))
	!								ELSEIF (ityp_TD(int(r(1))+1, int(r(2))+1, &
	!												int(r(3))+1, ib).eq.1) THEN
	!									mass2 = amass1(ityp1(nb))
	!								ENDIF
	!								Adding non-analytic part of FCs 
									IFC = (frc1(m1,m2,m3,:,:,na,nb)+f_of_q1(:,:,na,nb))
									DO ipol = 1,3
										DO jpol = 1,3
											counter1 = counter1 + 1
											ilist(counter1) = 3*n-(3-ipol)
		                                    jlist(counter1) = 3*Location-(3-jpol)
											Alist(counter1) = IFC(ipol,jpol)* &
															  weight/(amu_ry* &
															  sqrt(mass1*mass2)* &
															  w2(mode))
!!											my_BW = maximum(BW, KIABS((ilist(counter1)-jlist(counter1))))				  
											IF ((sqrt(mass1*mass2).eq.0.D0) .or. &
												(w2(mode).eq.0.D0)) THEN
												print *, 'HELLO'
											ENDIF
										ENDDO
									ENDDO
								ENDIF
							ENDDO
						ENDDO
					ENDDO
				ENDDO
			ELSEIF (ityp.eq.2) THEN
	!		Identical loop as above but accouting for lattice type of material 2
				DO n1 = -2*nr1,2*nr1
					DO n2 = -2*nr2,2*nr2
						DO n3 = -2*nr3,2*nr3
							DO nb = 1, nat(2)
								weight = wscache(n3,n2,n1,nb,na)
								IF (weight.gt.0) THEN
									m1 = MOD(n1+1,nr1)
									IF(m1.LE.0) m1=m1+nr1
									m2 = MOD(n2+1,nr2)
									IF(m2.LE.0) m2=m2+nr2
									m3 = MOD(n3+1,nr3)
									IF(m3.LE.0) m3=m3+nr3
									r = (i1-1)*atc(:,1) + (i2-1)*atc(:,2) + &
										(i3-1)*atc(:,3)
									r = r + r_cell(:,ia) + n1*at2(:,1) + n2*at2(:,2) + &
										n3*at2(:,3)+ tau2(:,nb)
									r = floor(r)
									ib = ib_vec2(nb,n3,n2,n1,ia)
									IF (periodic) THEN
										r(1)= MOD(r(1),TD(1)-1)
		                                IF (r(1)<0) THEN
											r(1)=r(1)+TD(1)-1
										ENDIF
										r(2)= MOD(r(2),TD(2)-1)
										IF (r(2)<0) THEN
											r(2)=r(2)+TD(2)-1
										ENDIF
									ENDIF
									IF ((r(1).lt.0) .or. (r(2).lt.0) .or. &
										(r(3).lt.0) .or. (r(1).gt.(TD(1)-1)) .or. &
										(r(2).gt.(TD(2)-1)) .or. (r(3).gt.(TD(3)-1))) THEN
										DO i = 1,3
											counter2 = counter2 + 1
											borderlogic(counter2) = 3*p - (3-i)
										ENDDO
										CYCLE
									ENDIF
									isub = (/ ib, int(r(1))+1, int(r(2))+1, &
												int(r(3))+1 /)
									Location = sub2ind(iSub, nSub)
									mass1 = amass_TD(i1,i2,i3,ia)
									mass2 = amass_TD(int(r(1))+1,int(r(2))+1,&
													 int(r(3))+1,ib)
	!								mass1 = amass2(ityp2(na))
	!								IF (ityp_TD(int(r(1))+1, int(r(2))+1, &
	!											int(r(3))+1, ib).eq.2) THEN
	!									mass2 = amass2(ityp2(nb))
	!								ELSEIF (ityp_TD(int(r(1))+1, int(r(2))+1, &
	!												int(r(3))+1, ib).eq.1) THEN
	!									mass2 = amass1(ityp1(nb))
	!								ENDIF
									IFC = (frc2(m1,m2,m3,:,:,na,nb)+f_of_q2(:,:,na,nb))
									DO ipol = 1,3
										DO jpol = 1,3
											counter1 = counter1 + 1
											ilist(counter1) = 3*n-(3-ipol)
		                                    jlist(counter1) = 3*Location-(3-jpol)
		                                    
!		                                    my_BW = maximum(BW, KIABS((ilist(counter1)-jlist(counter1))))
											Alist(counter1) = IFC(ipol,jpol)* &
															  weight/(amu_ry* &
															  sqrt(mass1*mass2)* &
															  w2(mode))
															  
											IF ((sqrt(mass1*mass2).eq.0.D0) .or. &
												(w2(mode).eq.0.D0)) THEN
												print *, 'HELLO'
											ENDIF
										ENDDO
									ENDDO
								ENDIF
							ENDDO
						ENDDO
					ENDDO
				ENDDO
			ENDIF	
		ENDDO
		
		
	!	CALL MPI_REDUCE(my_BW, BW, 1, mp_int, mp_maxi, root_node, comm, ierr)
		
	!	IF (io_node) THEN
	!		WRITE (stdout, '(a, I)') 'Bandwidth(A)  = ', BW
	!	END IF
	!	Adding -m*w**2 part after normalizing by m*w**2. So adding -1 :-P
		DO p = 1, my_natoms
			n = p + atoms_start(my_id+1)
			DO ipol = 1,3
				counter1 = counter1+1
				ilist(counter1) = 3*n-(3-ipol)
				jlist(counter1) = 3*n-(3-ipol)
				Alist(counter1) = CMPLX(-1.0, 0.0, KIND = CP)
			ENDDO
		ENDDO
		
		Alist(:) = CMPLX(-1.0, 0.0, kind = CP)*Alist(:)
			
		IF (io_node) WRITE (stdout, *) '	'
		IF (io_node) WRITE (stdout, *) 'DONE'
		
		print *, 'counter2', counter2
	END SUBROUTINE gen_Amat
	
!	========================================================================================
	
END MODULE matgen_module
