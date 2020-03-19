MODULE matgen_module

	CONTAINS
	
!	=============================================================================================

	SUBROUTINE gen_uinc(natc, nSub, my_natoms, my_nrows, zc, wavetype, my_uinc, atoms_start, &
						atc, r_cell, LPML, q)
		
		USE kinds
		USE essentials
		USE constants
		USE mp_module
		IMPLICIT NONE
		
		INTEGER								::	counter, p, i, natc, LPML
		INTEGER(KIND = IP)					::	my_natoms, my_nrows
		COMPLEX(KIND = CP), ALLOCATABLE		::	my_uinc(:)
		INTEGER(KIND = IP)					::	atom_tuple(4)
		COMPLEX(KIND = CP)					::	zc(3*natc)
		REAL								::	arg, atc(3,3)
		CHARACTER(len = 256)				::	wavetype
		INTEGER(KIND = IP)					::	nSub(4)
		INTEGER(KIND = IP)					::	atoms_start(world_size), n, n1, n2, n3, na
		REAL(KIND = RP)						::	r(3), q(3) 
		REAL(KIND = RP)						::	r_cell(3, natc)
		
		
		counter = 0
		ALLOCATE(my_uinc(my_nrows))
	
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
					my_uinc(counter)= zc(3*na-(3-i))* &
										CMPLX(COS(arg),-SIN(arg),KIND=CP)
				ELSEIF (wavetype .eq. 'half') THEN
					IF (n3.lt.LPML) THEN
						my_uinc(counter)= zc(3*na-(3-i))* &
											CMPLX(COS(arg),-SIN(arg),KIND=CP)
					ELSE
						my_uinc(counter) = CMPLX(0.0_RP, 0.0_RP, KIND = CP)
					ENDIF
				ENDIF
			ENDDO
		ENDDO
		
	END SUBROUTINE

!	=========================================================================================

	SUBROUTINE gen_sig(sig, my_natoms, sigmamax, LPML, periodic, nSub, atc, atom_tuple, &
						TD, PD, atoms_start, ityp_TD, tauc, natc)
	
		USE kinds
		USE mp_module
		USE essentials
		IMPLICIT NONE
		INTEGER								::	counter, p, i1, i2, i3, ia, ityp, natc
		INTEGER(KIND = IP)					::	my_natoms
		REAL(KIND = RP)						::	sigmamax
		LOGICAL								::	periodic
		INTEGER(KIND = IP)					::	nSub(4), n
		REAL								::	atc(3,3)
		REAL(KIND = RP), ALLOCATABLE		::	sig(:,:)
		INTEGER(KIND = IP)					::	atom_tuple(4), atoms_start(world_size)
		REAL(KIND = RP)						::	r(3)
		REAL								::	TD(3), PD(3)
		INTEGER								::	ityp_TD(int(TD(1)), int(TD(2)),&
														int(TD(3)), natc), LPML
		REAL(KIND = RP)						::	tauc(3,natc)
		
		
		
		
		
		ALLOCATE(sig(my_natoms, 3))
		sig(:,:) = 0.0_RP
		
		DO p = 1, my_natoms
			n = p + atoms_start(my_id+1)
			atom_tuple = ind2sub(n, nSub)
			ia = atom_tuple(1)
			i1 = atom_tuple(2)
			i2 = atom_tuple(3)
			i3 = atom_tuple(4)
			ityp = ityp_TD(i1,i2,i3,ia)
			r = (i1-1)*atc(:,1) + (i2-1)*atc(:,2) + (i3-1)*atc(:,3) + tauc(:,ia)
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
	
END MODULE matgen_module
