PROGRAM calculate_error

!	This program calculates the error between 2 sets of solutions created 
!	using FDPML. If the error is less than the tol then this program returns 
!	0 to calling script else returns a 1

	USE lapack95
	USE blas95
	USE f95_precision
	USE kinds
	use, intrinsic :: iso_fortran_env, only : stdin=>input_unit, &
											  stdout=>output_unit, &
											  stderr=>error_unit
	IMPLICIT NONE
	INCLUDE 'mkl.fi' ! MKL header file
	
	CHARACTER(len = 256)					::	filename, ref_filename, tmp_dir, flfrc, &
												format_string
	REAL(KIND = RP)							::	PD, TD
	INTEGER(KIND = IP)						::	LPML, natoms, nrows, nrowsm1
	LOGICAL									::	periodic, crystal_coordinates
	REAL(KIND = RP)							::	at, atc
	INTEGER									::	nat, natc
	INTEGER									::	nprocs
	INTEGER(KIND = IP), ALLOCATABLE			::	everyones_rows
	REAL(KIND = RP), ALLOCATABLE			::	uscat_calc(:), uscat_ref(:)
	REAL(KIND = RP)							::	error, tol
	
	tol = 1e-6
	
	NAMELIST /filenames/ flfrc, tmp_dir, ref_filename &
			 /system/ PD, LPML, periodic, crystal_coordinates &
			 /simulation/ nprocs
	
	ALLOCATE(everyones_rows(nprocs))
			 
	CALL readfc( flfrc1, at )
			 
	CALL get_natsc(at, atc, nat, natc, crystal_coordinates)
	
	natoms = int(TD(1)*TD(2)*TD(3)*natc)
	
	tol = tol*natoms
	
	CALL get_everyones_rows(nprocs, natoms, everyones_rows)
	
	IF (periodic) THEN
		TD = PD + (/ 0, 0, 2*LPML/)
	ELSE
		TD = PD + (/ 2*LPML, 2*LPML, 2*LPML/)
	ENDIF
	
	nrowsm1 = 1
	DO p = 1, nprocs
		IF (p .lt. 10) then
			format_string = "(a, a, a, I1, a)"
		ELSEIF (p.lt. 100) THEN
			format_string = "(a, a, a, I2, a)"
		ELSE 
			format_string = "(a, a, a, I3, a)"
		ENDIF
		nrows = everyones_rows(p)
		WRITE(filename, format) trim(tmp_dir), '/', 'uscat_', my_id, '.save'
		open (unit  = 639, file = filename, form = 'unformatted')
		read (unit = 639) uscat_calc(nrowsm1:nrows)
		close(unit = 639)
		nrowsm1 = nrows+1
	ENDIF
	
	open (unit  = 639, file = ref_filename, form = 'unformatted')
	read (unit = 639) uscat_ref
	close(unit = 639)
	
	CALL clean_uscat(uscat_calc, TD, PD, periodic)
	CALL clean_uscat(uscat_ref, TD, PD, periodic)
	
	error = nrm2((uscat_calc - uscat_ref))
	
	IF (error.le.tol) THEN
		STOP 0
	ELSE
		STOP 1
	END IF
	
END PROGRAM calculate_error

SUBROUTINE get_everyones_rows(nprocs, natoms, everyones_rows)
						
!	Calculate the number of atoms(rows) assigned to every processor

		USE kinds 
		IMPLICIT NONE
		INTEGER			 					:: 	p, nprocs
		INTEGER(KIND = IP)					::	natoms, rem
		INTEGER(KIND = IP)				 	:: 	everyones_rows(nprocs)
		
		DO p = 1, nprocs
			everyones_rows(p) = natoms/nprocs
			rem = MOD(natoms,nprocs)
		
			IF (p.gt.(nrpocs-rem)) THEN
					everyones_rows(p) = everyones_rows(p)+1
			ENDIF
		END DO
		
		everyones_rows = 3*everyones_rows
		
END SUBROUTINE

SUBROUTINE clean_uscat(uscat, TD, PD, periodic)

!	This subroutine implies the following rules to the scattered wave
!	For all p:
!		IF ( uscat(p) lies inside PD ) THEN
!			uscat(p) = abs(uscat(p))
!		ELSE
!			uscat(p) = p
!		END
!	END
	
	USE kinds
	USE essentials
	IMPLICIT NONE
	INCLUDE 'mpif.h' ! MPI header file
	REAL							::	PD(3), TD(3)
	INTEGER							::	n1, n2, n3, natc, na
	logical							::	periodic
	REAL(KIND = RP)					::	uscat

	nSub = (/natc, int(TD(1)), int(TD(2)), int(TD(3))/)

	IF (periodic) THEN
		DO n1 = 1, TD(1)
			DO n2 = 1, TD(2)
				DO n3 = 1, TD(3)
					DO na = 1, natc
						iSub = (/ na, n1, n2, n3/)
						iG = sub2ind(iSub, nSub)
						IF ((n1.gt.(TD(1)/2.D0-PD(1)/2.D0)) .and. &
							(n1.le.(TD(1)/2.D0+PD(1)/2.D0)) .and. &
							(n2.gt.(TD(2)/2.D0-PD(2)/2.D0)) .and. &
							(n2.le.(TD(2)/2.D0+PD(2)/2.D0)) .and. &
							(n3.gt.(TD(3)/2.D0-PD(3)/2.D0)) .and. &
							(n3.le.(TD(3)/2.D0+PD(3)/2.D0))) THEN
							uscat(iG) = abs(uscat(iG))
						ELSE (n3.lt.TD(3)/2) THEN
							uscat(iG) = 0.0_RP
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
							uscat(iG) = abs(uscat(iG))
						ELSE
							uscat(iG) = 0.0_RP
						ENDIF
					ENDDO
				ENDDO
			ENDDO
		ENDDO
	ENDIF		
		
END SUBROUTINE

!-----------------------------------------------------------------------
SUBROUTINE readfc ( flfrc, at )
!-----------------------------------------------------------------------
  
	USE mp_module
	USE constants,  ONLY : amu_ry
	USE essentials
	!
	IMPLICIT NONE
	! I/O variable
	INCLUDE 'mpif.h'
	CHARACTER(LEN=256) :: flfrc
	INTEGER :: ibrav, nr1,nr2,nr3,nat, ntyp
	REAL(KIND = RP) :: alat, at(3,3), epsil(3,3)
	LOGICAL :: has_zstar
	! local variables
	INTEGER :: i, j, na, nb, m1,m2,m3
	INTEGER :: ibid, jbid, nabid, nbbid, m1bid,m2bid,m3bid
	REAL(KIND = RP), ALLOCATABLE :: amass(:)
	REAL(KIND = RP) :: amass_from_file, omega
	INTEGER :: nt
	REAL(KIND = RP) :: celldm(6)
	CHARACTER(LEN=3) :: atm
	REAL(KIND = RP), ALLOCATABLE :: frc(:,:,:,:,:,:,:), tau(:,:),  zeu(:,:,:), m_loc(:,:)
	! frc : interatomic force constants in real space
	! tau_blk : atomic positions for the original cell
	! zeu : effective charges for the original cell
	! m_loc: the magnetic moments of each atom
	INTEGER, ALLOCATABLE  :: ityp(:)
	! ityp_blk : atomic types for each atom of the original cell
	!
	
	! Default values
	celldm(:) = 0.D0
	!
	!
	IF (io_node) OPEN (unit=1,file=flfrc,status='old',form='formatted')
	!
	!  read cell data
	!
	IF (io_node)THEN
	 READ(1,*) ntyp,nat,ibrav,(celldm(i),i=1,6)
	 if (ibrav==0) then
		read(1,*) ((at(i,j),i=1,3),j=1,3)
	 end if
	ENDIF
	CALL MPI_BCAST(ntyp, 1, MPI_INT, root_process, comm, ierr)
	CALL MPI_BCAST(nat, 1, MPI_INT, root_process, comm, ierr)
	CALL MPI_BCAST(ibrav, 1, MPI_INT, root_process, comm, ierr)
	CALL MPI_BCAST(celldm, 6, mp_real, root_process, comm, ierr)
	
	IF (ibrav==0) THEN
	 CALL MPI_BCAST(at, 9, mp_real, root_process, comm, ierr)
	ENDIF
	
	ALLOCATE(amass(ntyp))
	
	amass(:)=0.0_RP
	!
	CALL latgen(ibrav,celldm,at(1,1),at(1,2),at(1,3),omega)
	alat = celldm(1)
	at = at / alat !  bring at in units of alat
	CALL MPI_BCAST(at, 9, mp_real, root_process, comm, ierr)
	CALL cell_volume(at, alat, omega)
	CALL MPI_BCAST(omega, 1, mp_real, root_process, comm, ierr)
	!
	!  read atomic types, positions and masses
	!
	DO nt = 1,ntyp
	 IF (io_node) READ(1,*) i,atm,amass_from_file
	
	 CALL MPI_BCAST(i, 1, MPI_INT, root_process, comm, ierr)
	 CALL MPI_BCAST(atm, ntyp, MPI_CHARACTER, root_process, comm, ierr)
	 CALL MPI_BCAST(amass_from_file, 1, mp_real, root_process, comm, ierr)
	 IF (i.NE.nt) THEN
		write(*,*) ' ERROR'
		write(*,*) ' Error reading mass from file '
		write(*,*) ' ABORTING....'
		STOP
	 ENDIF
	 IF (amass(nt).EQ.0.d0) THEN
		amass(nt) = amass_from_file/amu_ry
	 ELSE
		WRITE(*,*) 'for atomic type',nt,' mass from file not used'
	 END IF
	END DO
	
	CALL MPI_BCAST(amass, ntyp, mp_real, root_process, comm, ierr)
	!
	ALLOCATE (tau(3,nat), ityp(nat), zeu(3,3,nat))
	!
	DO na=1,nat
	 IF (io_node) READ(1,*) i,ityp(na),(tau(j,na),j=1,3)
	 CALL MPI_BCAST(i, 1, MPI_INT, root_process, comm, ierr)
	 IF (i.NE.na) THEN
		write(*,*) ' ERROR'
		write(*,*) ' Error reading ityp from file'
		write(*,*) ' ABORTING....'
		STOP
	 ENDIF
	END DO
	CALL MPI_BCAST(ityp, nat, MPI_INT, root_process, comm, ierr)
	CALL MPI_BCAST(tau, 3*nat, mp_real, root_process, comm, ierr)
	!
	!  read macroscopic variable
	!
	IF (io_node) READ (1,*) has_zstar
	CALL MPI_BCAST(has_zstar, 1, mp_logical, root_process, comm, ierr)
	IF (has_zstar) THEN
	 IF (io_node) READ(1,*) ((epsil(i,j),j=1,3),i=1,3)
	 CALL MPI_BCAST(epsil, 9, mp_real, root_process, comm, ierr)
	 IF (io_node) THEN
		DO na=1,nat
		   READ(1,*)
		   READ(1,*) ((zeu(i,j,na),j=1,3),i=1,3)
		END DO
	 ENDIF
	  CALL MPI_BCAST(zeu, 9, mp_real, root_process, comm, ierr)
	ELSE
	 zeu  (:,:,:) = 0.d0
	 epsil(:,:) = 0.d0
	END IF
	!
	IF (io_node) READ (1,*) nr1,nr2,nr3
	CALL MPI_BCAST(nr1, 1, MPI_INT, root_process, comm, ierr)
	CALL MPI_BCAST(nr2, 1, MPI_INT, root_process, comm, ierr)
	CALL MPI_BCAST(nr3, 1, MPI_INT, root_process, comm, ierr)
	!
	!  read real-space interatomic force constants
	!
	ALLOCATE ( frc(nr1,nr2,nr3,3,3,nat,nat) )
	frc(:,:,:,:,:,:,:) = 0.d0
	DO i=1,3
	 DO j=1,3
		DO na=1,nat
		   DO nb=1,nat
			  IF (io_node) READ (1,*) ibid, jbid, nabid, nbbid
			  CALL MPI_BCAST(ibid, 1, MPI_INT, root_process, comm, ierr)
			  CALL MPI_BCAST(jbid, 1, MPI_INT, root_process, comm, ierr)
			  CALL MPI_BCAST(nabid, 1, MPI_INT, root_process, comm, ierr)
			  CALL MPI_BCAST(nbbid, 1, MPI_INT, root_process, comm, ierr)
			  IF(i .NE.ibid  .OR. j .NE.jbid .OR.                   &
				 na.NE.nabid .OR. nb.NE.nbbid) THEN
				 write(*,*) ' ERROR'
				 write(*,*) ' Error in reading force constants from frc file'
				 write(*,*) ' ABORTING....'
				 STOP
			  ENDIF
			  IF (io_node) READ (1,*) (((m1bid, m2bid, m3bid,        &
						  frc(m1,m2,m3,i,j,na,nb),                  &
						   m1=1,nr1),m2=1,nr2),m3=1,nr3)
			   
			  CALL MPI_BCAST(frc(:,:,:,i,j,na,nb), (nr1*nr2*nr3), mp_real, root_process, comm, ierr)
		   END DO
		END DO
	 END DO
	END DO
	!
	IF (io_node) CLOSE(unit=1)
	!
	RETURN
END SUBROUTINE readfc

!-------------------------------------------------------------------------
subroutine latgen(ibrav,celldm,a1,a2,a3,omega)
!-----------------------------------------------------------------------
!     sets up the crystallographic vectors a1, a2, and a3.
!
!     ibrav is the structure index:
!       1  cubic P (sc)                8  orthorhombic P
!       2  cubic F (fcc)               9  1-face (C) centered orthorhombic
!       3  cubic I (bcc)              10  all face centered orthorhombic
!       4  hexagonal and trigonal P   11  body centered orthorhombic
!       5  trigonal R, 3-fold axis c  12  monoclinic P (unique axis: c)
!       6  tetragonal P (st)          13  one face (base) centered monoclinic
!       7  tetragonal I (bct)         14  triclinic P
!     Also accepted:
!       0  "free" structure          -12  monoclinic P (unique axis: b)
!      -5  trigonal R, threefold axis along (111) 
!      -9  alternate description for base centered orthorhombic
!     -13  one face (base) centered monoclinic (unique axis: b)
!      91  1-face (A) centered orthorombic
!
!     celldm are parameters which fix the shape of the unit cell
!     omega is the unit-cell volume
!
!     NOTA BENE: all axis sets are right-handed
!     Boxes for US PPs do not work properly with left-handed axis
!
	implicit none
	integer, intent(in) :: ibrav
	real(RP), intent(inout) :: celldm(6)
	real(RP), intent(inout) :: a1(3), a2(3), a3(3)
	real(RP), intent(out) :: omega
	!
	real(RP), parameter:: sr2 = 1.414213562373d0, &
						sr3 = 1.732050807569d0
	integer :: i,j,k,l,iperm,ir
	real(RP) :: term, cbya, s, term1, term2, singam, sen
	!
	!  user-supplied lattice vectors
	!
	if (ibrav == 0) then
	 if (SQRT( a1(1)**2 + a1(2)**2 + a1(3)**2 ) == 0 )  THEN
		 write(*,*) ' Error in input data at lattice generation. Check input.'
		 STOP
	 endif
	 if (SQRT( a2(1)**2 + a2(2)**2 + a2(3)**2 ) == 0 )  then
		 write(*,*) ' Error in input data at lattice generation. Check input.'
		 stop
	 endif
	 if (SQRT( a3(1)**2 + a3(2)**2 + a3(3)**2 ) == 0 )  then
		 write(*,*) ' Error in input data at lattice generation. Check input.'
		 stop
	 endif
	
	 if ( celldm(1) /= 0.D0 ) then
	 !
	 ! ... input at are in units of alat => convert them to a.u.
	 !
		 a1(:) = a1(:) * celldm(1)
		 a2(:) = a2(:) * celldm(1)
		 a3(:) = a3(:) * celldm(1)
	 else
	 !
	 ! ... input at are in atomic units: define celldm(1) from a1
	 !
		 celldm(1) = SQRT( a1(1)**2 + a1(2)**2 + a1(3)**2 )
	 end if
	 !
	else
	 a1(:) = 0.d0
	 a2(:) = 0.d0
	 a3(:) = 0.d0
	end if
	!
	if (celldm (1) <= 0.d0) then
	 write(*,*) ' Error in input data at lattice generation. Check input.'
	 stop
	endif
	!
	!  index of bravais lattice supplied
	!
	if (ibrav == 1) then
	 !
	 !     simple cubic lattice
	 !
	 a1(1)=celldm(1)
	 a2(2)=celldm(1)
	 a3(3)=celldm(1)
	 !
	else if (ibrav == 2) then
	 !
	 !     fcc lattice
	 !
	 term=celldm(1)/2.d0
	 a1(1)=-term
	 a1(3)=term
	 a2(2)=term
	 a2(3)=term
	 a3(1)=-term
	 a3(2)=term
	 !
	else if (ibrav == 3) then
	 !
	 !     bcc lattice
	 !
	 term=celldm(1)/2.d0
	 do ir=1,3
		a1(ir)=term
		a2(ir)=term
		a3(ir)=term
	 end do
	 a2(1)=-term
	 a3(1)=-term
	 a3(2)=-term
	 !
	else if (ibrav == 4) then
	 !
	 !     hexagonal lattice
	 !
	 if (celldm (3) <= 0.d0) then 
		write(*,*) ' Error in input data at lattice generation. Check input.'
		stop
	 endif
	 !
	 cbya=celldm(3)
	 a1(1)=celldm(1)
	 a2(1)=-celldm(1)/2.d0
	 a2(2)=celldm(1)*sr3/2.d0
	 a3(3)=celldm(1)*cbya
	 !
	else if (ABS(ibrav) == 5) then
	 !
	 !     trigonal lattice
	 !
	 if (celldm (4) <= -0.5_RP .or. celldm (4) >= 1.0_RP) then
		write(*,*) ' Error in input data at lattice generation. Check input.'
		stop
	 endif
	 !
	 term1=sqrt(1.0_RP + 2.0_RP*celldm(4))
	 term2=sqrt(1.0_RP - celldm(4))
	 !
	 IF ( ibrav == 5) THEN
		!     threefold axis along c (001)
		a2(2)=sr2*celldm(1)*term2/sr3
		a2(3)=celldm(1)*term1/sr3
		a1(1)=celldm(1)*term2/sr2
		a1(2)=-a1(1)/sr3
		a1(3)= a2(3)
		a3(1)=-a1(1)
		a3(2)= a1(2)
		a3(3)= a2(3)
	 ELSE IF ( ibrav == -5) THEN
		!     threefold axis along (111)
		! Notice that in the cubic limit (alpha=90, celldm(4)=0, term1=term2=1)
		! does not yield the x,y,z axis, but an equivalent rotated triplet:
		!    a/3 (-1,2,2), a/3 (2,-1,2), a/3 (2,2,-1)
		! If you prefer the x,y,z axis as cubic limit, you should modify the
		! definitions of a1(1) and a1(2) as follows:'
		!    a1(1) = celldm(1)*(term1+2.0_dp*term2)/3.0_dp
		!    a1(2) = celldm(1)*(term1-term2)/3.0_dp
		! (info by G. Pizzi and A. Cepellotti)
		!
		a1(1) = celldm(1)*(term1-2.0_RP*term2)/3.0_RP
		a1(2) = celldm(1)*(term1+term2)/3.0_RP
		a1(3) = a1(2)
		a2(1) = a1(3)
		a2(2) = a1(1)
		a2(3) = a1(2)
		a3(1) = a1(2)
		a3(2) = a1(3)
		a3(3) = a1(1)
	 END IF
	else if (ibrav == 6) then
	 !
	 !     tetragonal lattice
	 !
	 if (celldm (3) <= 0.d0) then
		write(*,*) ' Error in input data at lattice generation. Check input.'
		stop
	 endif
	 !
	 cbya=celldm(3)
	 a1(1)=celldm(1)
	 a2(2)=celldm(1)
	 a3(3)=celldm(1)*cbya
	 !
	else if (ibrav == 7) then
	 !
	 !     body centered tetragonal lattice
	 !
	 if (celldm (3) <= 0.d0) then
		write(*,*) ' Error in input data at lattice generation. Check input.'
		stop
	 endif
	 !
	 cbya=celldm(3)
	 a2(1)=celldm(1)/2.d0
	 a2(2)=a2(1)
	 a2(3)=cbya*celldm(1)/2.d0
	 a1(1)= a2(1)
	 a1(2)=-a2(1)
	 a1(3)= a2(3)
	 a3(1)=-a2(1)
	 a3(2)=-a2(1)
	 a3(3)= a2(3)
	 !
	else if (ibrav == 8) then
	 !
	 !     Simple orthorhombic lattice
	 !
	 if (celldm (2) <= 0.d0) then
		write(*,*) ' Error in input data at lattice generation. Check input.'
		stop
	 endif
	 if (celldm (3) <= 0.d0) then
		write(*,*) ' Error in input data at lattice generation. Check input.'
		stop
	 endif
	 !
	 a1(1)=celldm(1)
	 a2(2)=celldm(1)*celldm(2)
	 a3(3)=celldm(1)*celldm(3)
	 !
	else if ( ABS(ibrav) == 9) then
	 !
	 !     One face (base) centered orthorhombic lattice  (C type)
	 !
	 if (celldm (2) <= 0.d0) then
		write(*,*) ' Error in input data at lattice generation. Check input.'
		stop
	 endif
	 if (celldm (3) <= 0.d0) then
		write(*,*) ' Error in input data at lattice generation. Check input.'
		stop
	 endif
	 !
	 IF ( ibrav == 9 ) THEN
		!   old PWscf description
		a1(1) = 0.5d0 * celldm(1)
		a1(2) = a1(1) * celldm(2)
		a2(1) = - a1(1)
		a2(2) = a1(2)
	 ELSE
		!   alternate description
		a1(1) = 0.5d0 * celldm(1)
		a1(2) =-a1(1) * celldm(2)
		a2(1) = a1(1)
		a2(2) =-a1(2)
	 END IF
	 a3(3) = celldm(1) * celldm(3)
	 !
	else if ( ibrav == 91 ) then
	 !
	 !     One face (base) centered orthorhombic lattice  (A type)
	 !
	 if (celldm (2) <= 0.d0) then
		write(*,*) ' Error in input data at lattice generation. Check input.'
		stop
	 endif
	 if (celldm (3) <= 0.d0) then
		write(*,*) ' Error in input data at lattice generation. Check input.'
		stop
	 endif
	 !
	 a1(1) = celldm(1)
	 a2(2) = celldm(1) * celldm(2) * 0.5_RP
	 a2(3) = - celldm(1) * celldm(3) * 0.5_RP
	 a3(2) = a2(2)
	 a3(3) = - a2(3)
	 !
	else if (ibrav == 10) then
	 !
	 !     All face centered orthorhombic lattice
	 !
	 if (celldm (2) <= 0.d0) then
		write(*,*) ' Error in input data at lattice generation. Check input.'
		stop
	 endif
	 if (celldm (3) <= 0.d0) then
		write(*,*) ' Error in input data at lattice generation. Check input.'
		stop
	 endif
	 !
	 a2(1) = 0.5d0 * celldm(1)
	 a2(2) = a2(1) * celldm(2)
	 a1(1) = a2(1)
	 a1(3) = a2(1) * celldm(3)
	 a3(2) = a2(1) * celldm(2)
	 a3(3) = a1(3)
	 !
	else if (ibrav == 11) then
	 !
	 !     Body centered orthorhombic lattice
	 !
	 if (celldm (2) <= 0.d0) then
		write(*,*) ' Error in input data at lattice generation. Check input.'
		stop
	 endif
	 if (celldm (3) <= 0.d0) then
		write(*,*) ' Error in input data at lattice generation. Check input.'
		stop
	 endif
	 !
	 a1(1) = 0.5d0 * celldm(1)
	 a1(2) = a1(1) * celldm(2)
	 a1(3) = a1(1) * celldm(3)
	 a2(1) = - a1(1)
	 a2(2) = a1(2)
	 a2(3) = a1(3)
	 a3(1) = - a1(1)
	 a3(2) = - a1(2)
	 a3(3) = a1(3)
	 !
	else if (ibrav == 12) then
	 !
	 !     Simple monoclinic lattice, unique (i.e. orthogonal to a) axis: c
	 !
	 if (celldm (2) <= 0.d0) then
		write(*,*) ' Error in input data at lattice generation. Check input.'
		stop
	 endif
	 if (celldm (3) <= 0.d0) then
		write(*,*) ' Error in input data at lattice generation. Check input.'
		stop
	 endif
	 if (abs(celldm(4))>=1.d0) then
		write(*,*) ' Error in input data at lattice generation. Check input.'
		stop
	 endif
	 !
	 sen=sqrt(1.d0-celldm(4)**2)
	 a1(1)=celldm(1)
	 a2(1)=celldm(1)*celldm(2)*celldm(4)
	 a2(2)=celldm(1)*celldm(2)*sen
	 a3(3)=celldm(1)*celldm(3)
	 !
	else if (ibrav ==-12) then
	 !
	 !     Simple monoclinic lattice, unique axis: b (more common)
	 !
	 if (celldm (2) <= 0.d0) then
		write(*,*) ' Error in input data at lattice generation. Check input.'
		stop
	 endif
	 if (celldm (3) <= 0.d0) then
		write(*,*) ' Error in input data at lattice generation. Check input.'
		stop
	 endif
	 if (abs(celldm(5))>=1.d0) then
		write(*,*) ' Error in input data at lattice generation. Check input.'
		stop
	 endif
	 !
	 sen=sqrt(1.d0-celldm(5)**2)
	 a1(1)=celldm(1)
	 a2(2)=celldm(1)*celldm(2)
	 a3(1)=celldm(1)*celldm(3)*celldm(5)
	 a3(3)=celldm(1)*celldm(3)*sen
	 !
	else if (ibrav == 13) then
	 !
	 !     One face centered monoclinic lattice unique axis c
	 !
	 if (celldm (2) <= 0.d0) then
		write(*,*) ' Error in input data at lattice generation. Check input.'
		stop
	 endif
	 if (celldm (3) <= 0.d0) then
		write(*,*) ' Error in input data at lattice generation. Check input.'
		stop
	 endif
	 if (abs(celldm(4))>=1.d0) then
		write(*,*) ' Error in input data at lattice generation. Check input.'
		stop
	 endif
	 !
	 sen = sqrt( 1.d0 - celldm(4) ** 2 )
	 a1(1) = 0.5d0 * celldm(1) 
	 a1(3) =-a1(1) * celldm(3)
	 a2(1) = celldm(1) * celldm(2) * celldm(4)
	 a2(2) = celldm(1) * celldm(2) * sen
	 a3(1) = a1(1)
	 a3(3) =-a1(3)
	else if (ibrav == -13) then
	 !
	 !     One face centered monoclinic lattice unique axis b
	 !
	 if (celldm (2) <= 0.d0) then
		write(*,*) ' Error in input data at lattice generation. Check input.'
		stop
	 endif
	 if (celldm (3) <= 0.d0) then
		write(*,*) ' Error in input data at lattice generation. Check input.'
		stop
	 endif
	 if (abs(celldm(5))>=1.d0) then
		write(*,*) ' Error in input data at lattice generation. Check input.'
		stop
	 endif
	 !
	 sen = sqrt( 1.d0 - celldm(5) ** 2 )
	 a1(1) = 0.5d0 * celldm(1) 
	 a1(2) =-a1(1) * celldm(2)
	 a2(1) = a1(1)
	 a2(2) =-a1(2)
	 a3(1) = celldm(1) * celldm(3) * celldm(5)
	 a3(3) = celldm(1) * celldm(3) * sen
	 !
	else if (ibrav == 14) then
	 !
	 !     Triclinic lattice
	 !
	 if (celldm (2) <= 0.d0) then
		write(*,*) ' Error in input data at lattice generation. Check input.'
		stop
	 endif
	 if (celldm (3) <= 0.d0) then
		write(*,*) ' Error in input data at lattice generation. Check input.'
		stop
	 endif
	 if (abs(celldm(4))>=1.d0) then
		write(*,*) ' Error in input data at lattice generation. Check input.'
		stop
	 endif
	 if (abs(celldm(5))>=1.d0) then
		write(*,*) ' Error in input data at lattice generation. Check input.'
		stop
	 endif
	 if (abs(celldm(6))>=1.d0) then
		write(*,*) ' Error in input data at lattice generation. Check input.'
		stop
	 endif
	 !
	 singam=sqrt(1.d0-celldm(6)**2)
	 term= (1.d0+2.d0*celldm(4)*celldm(5)*celldm(6)             &
		  -celldm(4)**2-celldm(5)**2-celldm(6)**2)
	 if (term < 0.d0) then
		write(*,*) ' Error in input data at lattice generation. Check input.'
		stop
	 endif
	 term= sqrt(term/(1.d0-celldm(6)**2))
	 a1(1)=celldm(1)
	 a2(1)=celldm(1)*celldm(2)*celldm(6)
	 a2(2)=celldm(1)*celldm(2)*singam
	 a3(1)=celldm(1)*celldm(3)*celldm(5)
	 a3(2)=celldm(1)*celldm(3)*(celldm(4)-celldm(5)*celldm(6))/singam
	 a3(3)=celldm(1)*celldm(3)*term
	 !
	else
	 !
	 write(*,*) ' Error lattice type non-existent'
	 stop
	 !
	end if
	!
	!  calculate unit-cell volume omega
	!
	omega=0.d0
	s=1.d0
	i=1
	j=2
	k=3
	!
	101 do iperm=1,3
	 omega=omega+s*a1(i)*a2(j)*a3(k)
	 l=i
	 i=j
	 j=k
	 k=l
	end do
	!
	i=2
	j=1
	k=3
	s=-s
	if(s < 0.d0) go to 101
	omega=abs(omega)
	return
!
end subroutine latgen
