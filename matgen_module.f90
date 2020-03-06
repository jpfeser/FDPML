MODULE matgen_module

	CONTAINS
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
END MODULE matgen_module
