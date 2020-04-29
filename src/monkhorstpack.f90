MODULE monkhorstpack

  USE kinds,      ONLY : RP
  use, intrinsic :: iso_fortran_env, only : stdin=>input_unit, &
											 stdout=>output_unit, &
											 stderr=>error_unit
    
  INTEGER :: &
       s(3,3,48),            &! symmetry matrices, in crystal axis
       invs(48),             &! index of inverse operation: S^{-1}_i=S(invs(i))
       ftau(3,48),           &! fractional translations, in FFT coordinates
       nrot,                 &! number of bravais lattice symmetries
       nsym = 1,             &! total number of crystal symmetries
       nsym_ns = 0,          &! nonsymmorphic (fractional translation) symms
       nsym_na = 0            ! excluded nonsymmorphic symmetries because
                              ! fract. transl. is noncommensurate with FFT grid
  REAL (RP) :: &
       ft (3,48),            &! fractional translations, in crystal axis
       sr (3,3,48),          &! symmetry matrices, in cartesian axis
       accep = 1.0d-5         ! initial value of the acceptance threshold
                              ! for position comparison by eqvect in checksym
  !
  CHARACTER(len=45) ::  sname(48)   ! name of the symmetries
  INTEGER :: &
       t_rev(48) = 0          ! time reversal flag, for noncolinear magnetism
  INTEGER, ALLOCATABLE :: &
       irt(:,:)               ! symmetric atom for each atom and sym.op.
  LOGICAL :: &
       time_reversal=.true., &! if .TRUE. the system has time reversal symmetry
       invsym,               &! if .TRUE. the system has inversion symmetry
       nofrac= .false.,      &! if .TRUE. fract. translations are not allowed
       allfrac= .false.,     &! if .TRUE. all fractionary translations allowed,
                              ! even those not commensurate with FFT grid
       nosym = .false.,      &! if .TRUE. no symmetry is used
       nosym_evc = .false.,  &! if .TRUE. symmetry is used only to symmetrize
                              ! k points
       no_t_rev=.false.       ! if .TRUE. remove the symmetries that
                              ! require time reversal
  REAL(RP),TARGET :: &
       d1(3,3,48),           &! matrices for rotating spherical
       d2(5,5,48),           &! harmonics (d1 for l=1, ...)
       d3(7,7,48)             !
  !
  REAL(RP) :: at(3,3) = RESHAPE( (/ 0.0_RP /), (/ 3, 3 /), (/ 0.0_RP /) )
  REAL(RP) :: bg(3,3) = RESHAPE( (/ 0.0_RP /), (/ 3, 3 /), (/ 0.0_RP /) )

  CONTAINS
  
  SUBROUTINE gen_qpoints (ibrav, at_, bg_, nat, tau, ityp, nk1, nk2, nk3, &
     ntetra, nqx, nq, q, tetra)
	  !-----------------------------------------------------------------------
	  !
	  IMPLICIT NONE
	  ! input
	  INTEGER :: ibrav, nat, nk1, nk2, nk3, ntetra, ityp(*)
	  REAL(RP) :: at_(3,3), bg_(3,3), tau(3,nat)
	  ! output
	  INTEGER :: nqx, nq, tetra(4,ntetra)
	  REAL(RP) :: q(3,nqx)
	  ! local
	  REAL(RP) :: xqq(3), wk(nqx), mdum(3,nat)
	  LOGICAL :: magnetic_sym=.FALSE., skip_equivalence=.false.
	  !
	  time_reversal = .false.
	  t_rev(:) = 0
	  xqq (:) =0.d0
	  at = at_
	  bg = bg_
	  CALL set_sym_bl ( )
	  !
	  CALL kpoint_grid ( nrot, time_reversal, skip_equivalence, s, t_rev, bg, nqx, &
	                           0,0,0, nk1,nk2,nk3, nq, q, wk)
	  !
!	  CALL find_sym ( nat, tau, ityp, 6, 6, 6, .not.time_reversal, mdum )
	  !
!	  CALL irreducible_BZ (nrot, s, nsym, time_reversal, magnetic_sym, &
!	                       at, bg, nqx, nq, q, wk, t_rev)
	  !
!	  IF (ntetra /= 6 * nk1 * nk2 * nk3) &
!	       CALL errore ('gen_qpoints','inconsistent ntetra',1)
	  !
!	  CALL tetrahedra (nsym, s, time_reversal, t_rev, at, bg, nqx, 0, 0, 0, &
!	       nk1, nk2, nk3, nq, q, wk, ntetra, tetra)
	  !
	  RETURN
  END SUBROUTINE gen_qpoints
  
  !-----------------------------------------------------------------------
  SUBROUTINE set_sym_bl ( )
  !-----------------------------------------------------------------------
	  !
	  ! Provides symmetry operations for all bravais lattices
	  ! Tests first the 24 proper rotations for the cubic lattice;
	  ! then the 8 rotations specific for the hexagonal axis (special axis c);
	  ! then inversion is added
	  !
	  USE matrix_inversion
	  
	  IMPLICIT NONE
	  !
	  ! sin3 = sin(pi/3), cos3 = cos(pi/3), msin3 = -sin(pi/3), mcos3 = -cos(pi/3)
	  !
	  real(RP), PARAMETER :: sin3 = 0.866025403784438597d0, cos3 = 0.5d0, &
	                        msin3 =-0.866025403784438597d0, mcos3 = -0.5d0, &
	                        eps1 = 1.0d-6
	  real(RP) :: s0(3, 3, 32), overlap (3, 3), rat (3), rot (3, 3), value
	  ! s0: the s matrices in cartesian axis
	  ! overlap: inverse overlap matrix between direct lattice
	  ! rat: the rotated of a direct vector ( cartesian )
	  ! rot: the rotated of a direct vector ( crystal axis )
	  ! value: component of the s matrix in axis basis
	  INTEGER :: jpol, kpol, mpol, irot, imat(24)
	  ! counters over the polarizations and the rotations
	
	  CHARACTER (len=45) :: s0name (64)
	  ! full name of the rotational part of each symmetry operation
	
	  data s0/ 1.d0,  0.d0,  0.d0,  0.d0,  1.d0,  0.d0,  0.d0,  0.d0,  1.d0, &
	          -1.d0,  0.d0,  0.d0,  0.d0, -1.d0,  0.d0,  0.d0,  0.d0,  1.d0, &
	          -1.d0,  0.d0,  0.d0,  0.d0,  1.d0,  0.d0,  0.d0,  0.d0, -1.d0, &
	           1.d0,  0.d0,  0.d0,  0.d0, -1.d0,  0.d0,  0.d0,  0.d0, -1.d0, &
	           0.d0,  1.d0,  0.d0,  1.d0,  0.d0,  0.d0,  0.d0,  0.d0, -1.d0, &
	           0.d0, -1.d0,  0.d0, -1.d0,  0.d0,  0.d0,  0.d0,  0.d0, -1.d0, &
	           0.d0, -1.d0,  0.d0,  1.d0,  0.d0,  0.d0,  0.d0,  0.d0,  1.d0, &
	           0.d0,  1.d0,  0.d0, -1.d0,  0.d0,  0.d0,  0.d0,  0.d0,  1.d0, &
	           0.d0,  0.d0,  1.d0,  0.d0, -1.d0,  0.d0,  1.d0,  0.d0,  0.d0, &
	           0.d0,  0.d0, -1.d0,  0.d0, -1.d0,  0.d0, -1.d0,  0.d0,  0.d0, &
	           0.d0,  0.d0, -1.d0,  0.d0,  1.d0,  0.d0,  1.d0,  0.d0,  0.d0, &
	           0.d0,  0.d0,  1.d0,  0.d0,  1.d0,  0.d0, -1.d0,  0.d0,  0.d0, &
	          -1.d0,  0.d0,  0.d0,  0.d0,  0.d0,  1.d0,  0.d0,  1.d0,  0.d0, &
	          -1.d0,  0.d0,  0.d0,  0.d0,  0.d0, -1.d0,  0.d0, -1.d0,  0.d0, &
	           1.d0,  0.d0,  0.d0,  0.d0,  0.d0, -1.d0,  0.d0,  1.d0,  0.d0, &
	           1.d0,  0.d0,  0.d0,  0.d0,  0.d0,  1.d0,  0.d0, -1.d0,  0.d0, &
	           0.d0,  0.d0,  1.d0,  1.d0,  0.d0,  0.d0,  0.d0,  1.d0,  0.d0, &
	           0.d0,  0.d0, -1.d0, -1.d0,  0.d0,  0.d0,  0.d0,  1.d0,  0.d0, &
	           0.d0,  0.d0, -1.d0,  1.d0,  0.d0,  0.d0,  0.d0, -1.d0,  0.d0, &
	           0.d0,  0.d0,  1.d0, -1.d0,  0.d0,  0.d0,  0.d0, -1.d0,  0.d0, &
	           0.d0,  1.d0,  0.d0,  0.d0,  0.d0,  1.d0,  1.d0,  0.d0,  0.d0, &
	           0.d0, -1.d0,  0.d0,  0.d0,  0.d0, -1.d0,  1.d0,  0.d0,  0.d0, &
	           0.d0, -1.d0,  0.d0,  0.d0,  0.d0,  1.d0, -1.d0,  0.d0,  0.d0, &
	           0.d0,  1.d0,  0.d0,  0.d0,  0.d0, -1.d0, -1.d0,  0.d0,  0.d0, &
	           cos3,  sin3, 0.d0, msin3,  cos3, 0.d0, 0.d0, 0.d0,  1.d0, &
	           cos3, msin3, 0.d0,  sin3,  cos3, 0.d0, 0.d0, 0.d0,  1.d0, &
	          mcos3,  sin3, 0.d0, msin3, mcos3, 0.d0, 0.d0, 0.d0,  1.d0, &
	          mcos3, msin3, 0.d0,  sin3, mcos3, 0.d0, 0.d0, 0.d0,  1.d0, &
	           cos3, msin3, 0.d0, msin3, mcos3, 0.d0, 0.d0, 0.d0, -1.d0, &
	           cos3,  sin3, 0.d0,  sin3, mcos3, 0.d0, 0.d0, 0.d0, -1.d0, &
	          mcos3, msin3, 0.d0, msin3,  cos3, 0.d0, 0.d0, 0.d0, -1.d0, &
	          mcos3,  sin3, 0.d0,  sin3,  cos3, 0.d0, 0.d0, 0.d0, -1.d0 /
	
	  data s0name/  'identity                                     ',&
	                '180 deg rotation - cart. axis [0,0,1]        ',&
	                '180 deg rotation - cart. axis [0,1,0]        ',&
	                '180 deg rotation - cart. axis [1,0,0]        ',&
	                '180 deg rotation - cart. axis [1,1,0]        ',&
	                '180 deg rotation - cart. axis [1,-1,0]       ',&
	                ' 90 deg rotation - cart. axis [0,0,-1]       ',&
	                ' 90 deg rotation - cart. axis [0,0,1]        ',&
	                '180 deg rotation - cart. axis [1,0,1]        ',&
	                '180 deg rotation - cart. axis [-1,0,1]       ',&
	                ' 90 deg rotation - cart. axis [0,1,0]        ',&
	                ' 90 deg rotation - cart. axis [0,-1,0]       ',&
	                '180 deg rotation - cart. axis [0,1,1]        ',&
	                '180 deg rotation - cart. axis [0,1,-1]       ',&
	                ' 90 deg rotation - cart. axis [-1,0,0]       ',&
	                ' 90 deg rotation - cart. axis [1,0,0]        ',&
	                '120 deg rotation - cart. axis [-1,-1,-1]     ',&
	                '120 deg rotation - cart. axis [-1,1,1]       ',&
	                '120 deg rotation - cart. axis [1,1,-1]       ',&
	                '120 deg rotation - cart. axis [1,-1,1]       ',&
	                '120 deg rotation - cart. axis [1,1,1]        ',&
	                '120 deg rotation - cart. axis [-1,1,-1]      ',&
	                '120 deg rotation - cart. axis [1,-1,-1]      ',&
	                '120 deg rotation - cart. axis [-1,-1,1]      ',&
	                ' 60 deg rotation - cryst. axis [0,0,1]       ',&
	                ' 60 deg rotation - cryst. axis [0,0,-1]      ',&
	                '120 deg rotation - cryst. axis [0,0,1]       ',&
	                '120 deg rotation - cryst. axis [0,0,-1]      ',&
	                '180 deg rotation - cryst. axis [1,-1,0]      ',&
	                '180 deg rotation - cryst. axis [2,1,0]       ',&
	                '180 deg rotation - cryst. axis [0,1,0]       ',&
	                '180 deg rotation - cryst. axis [1,1,0]       ',&
	                'inversion                                    ',&
	                'inv. 180 deg rotation - cart. axis [0,0,1]   ',&
	                'inv. 180 deg rotation - cart. axis [0,1,0]   ',&
	                'inv. 180 deg rotation - cart. axis [1,0,0]   ',&
	                'inv. 180 deg rotation - cart. axis [1,1,0]   ',&
	                'inv. 180 deg rotation - cart. axis [1,-1,0]  ',&
	                'inv.  90 deg rotation - cart. axis [0,0,-1]  ',&
	                'inv.  90 deg rotation - cart. axis [0,0,1]   ',&
	                'inv. 180 deg rotation - cart. axis [1,0,1]   ',&
	                'inv. 180 deg rotation - cart. axis [-1,0,1]  ',&
	                'inv.  90 deg rotation - cart. axis [0,1,0]   ',&
	                'inv.  90 deg rotation - cart. axis [0,-1,0]  ',&
	                'inv. 180 deg rotation - cart. axis [0,1,1]   ',&
	                'inv. 180 deg rotation - cart. axis [0,1,-1]  ',&
	                'inv.  90 deg rotation - cart. axis [-1,0,0]  ',&
	                'inv.  90 deg rotation - cart. axis [1,0,0]   ',&
	                'inv. 120 deg rotation - cart. axis [-1,-1,-1]',&
	                'inv. 120 deg rotation - cart. axis [-1,1,1]  ',&
	                'inv. 120 deg rotation - cart. axis [1,1,-1]  ',&
	                'inv. 120 deg rotation - cart. axis [1,-1,1]  ',&
	                'inv. 120 deg rotation - cart. axis [1,1,1]   ',&
	                'inv. 120 deg rotation - cart. axis [-1,1,-1] ',&
	                'inv. 120 deg rotation - cart. axis [1,-1,-1] ',&
	                'inv. 120 deg rotation - cart. axis [-1,-1,1] ',&
	                'inv.  60 deg rotation - cryst. axis [0,0,1]  ',&
	                'inv.  60 deg rotation - cryst. axis [0,0,-1] ',&
	                'inv. 120 deg rotation - cryst. axis [0,0,1]  ',&
	                'inv. 120 deg rotation - cryst. axis [0,0,-1] ',&
	                'inv. 180 deg rotation - cryst. axis [1,-1,0] ',&
	                'inv. 180 deg rotation - cryst. axis [2,1,0]  ',&
	                'inv. 180 deg rotation - cryst. axis [0,1,0]  ',&
	                'inv. 180 deg rotation - cryst. axis [1,1,0]  ' /
	
	  !    compute the overlap matrix for crystal axis
	
	  DO jpol = 1,3
	     DO kpol = 1,3
	        rot(kpol,jpol) = at(1,kpol)*at(1,jpol) +&
	                         at(2,kpol)*at(2,jpol) +&
	                         at(3,kpol)*at(3,jpol)
	     ENDDO
	  ENDDO
	  !
	  !    then its inverse (rot is used as work space)
	  !
	  CALL invmat (3, rot, overlap, value)
	
	  nrot = 1
	  DO irot = 1,32
	     !
	     !   for each possible symmetry
	     !
	     DO jpol = 1,3
	        DO mpol = 1,3
	           !
	           !   compute, in cartesian coordinates the rotated vector
	           !
	           rat(mpol) = s0(mpol,1,irot)*at(1,jpol) +&
	                       s0(mpol,2,irot)*at(2,jpol) +&
	                       s0(mpol,3,irot)*at(3,jpol)
	        ENDDO
	
	        DO kpol = 1,3
	           !
	           !   the rotated vector is projected on the direct lattice
	           !
	           rot(kpol,jpol) = at(1,kpol)*rat(1) +&
	                            at(2,kpol)*rat(2) +&
	                            at(3,kpol)*rat(3)
	        ENDDO
	     ENDDO
	     !
	     !  and the inverse of the overlap matrix is applied
	     !
	     DO jpol = 1,3
	        DO kpol = 1,3
	           value = overlap(jpol,1)*rot(1,kpol) +&
	           &       overlap(jpol,2)*rot(2,kpol) +&
	           &       overlap(jpol,3)*rot(3,kpol)
	           IF ( abs(dble(nint(value))-value) > eps1 ) THEN
	              !
	              ! if a noninteger is obtained, this implies that this operation
	              ! is not a symmetry operation for the given lattice
	              !
	              GOTO 10
	           ENDIF
	           s(kpol,jpol,nrot) = nint(value)
	        ENDDO
	     ENDDO
	     sname(nrot)=s0name(irot)
	     imat(nrot)=irot
	     nrot = nrot+1
	     IF (nrot > 25) THEN
			WRITE(stdout, '(a)') 'ERROR : some error in symmetry'
			STOP
		 ENDIF
 	10   CONTINUE
	  ENDDO
	  nrot = nrot-1
	  IF ( nrot /= 1 .AND. nrot /= 2 .AND. nrot /= 4 .AND. nrot /= 6 .AND. &
	       nrot /= 8 .AND. nrot /=12 .AND. nrot /=24 ) THEN
	       WRITE(stdout, '(a)') 'ERROR : some error in symmetry'
		   STOP
	  ENDIF
	!
	  !     set the inversion symmetry ( Bravais lattices have always inversion
	  !     symmetry )
	  !
	  DO irot = 1, nrot
	     sname(irot+nrot) = s0name(imat(irot)+32)
	     DO kpol = 1,3
	        DO jpol = 1,3
	           s(kpol,jpol,irot+nrot) = -s(kpol,jpol,irot)
	        ENDDO
	     ENDDO
	  ENDDO
	  nrot = 2*nrot
	  !
	  !    reset fractional translations to zero before checking the group
	  ! 
	  ft(:,:) = 0.0_RP
	  !
	  RETURN
	  !
  END SUBROUTINE set_sym_bl	
  
  !-----------------------------------------------------------------------
  SUBROUTINE kpoint_grid ( nrot, time_reversal, skip_equivalence, s, t_rev, &
                         bg, npk, k1,k2,k3, nk1,nk2,nk3, nks, xk, wk)
  !-----------------------------------------------------------------------
	!
	!  Automatic generation of a uniform grid of k-points
	!
	  USE kinds, ONLY: RP
	  IMPLICIT NONE
	  !
	  INTEGER, INTENT(in):: nrot, npk, k1, k2, k3, nk1, nk2, nk3, &
	                        t_rev(48), s(3,3,48)
	  LOGICAL, INTENT(in):: time_reversal, skip_equivalence
	  real(RP), INTENT(in):: bg(3,3)
	  !
	  INTEGER, INTENT(out) :: nks
	  real(RP), INTENT(out):: xk(3,npk)
	  real(RP), INTENT(out):: wk(npk)
	  ! LOCAL:
	  real(RP), PARAMETER :: eps=1.0d-5
	  real(RP) :: xkr(3), fact, xx, yy, zz
	  real(RP), ALLOCATABLE:: xkg(:,:), wkk(:)
	  INTEGER :: nkr, i,j,k, ns, n, nk
	  INTEGER, ALLOCATABLE :: equiv(:)
	  LOGICAL :: in_the_list
	  !
	  nkr=nk1*nk2*nk3
	  ALLOCATE (xkg( 3,nkr),wkk(nkr))
	  ALLOCATE (equiv( nkr))
	  !
	  DO i=1,nk1
	     DO j=1,nk2
	        DO k=1,nk3
	           !  this is nothing but consecutive ordering
	           n = (k-1) + (j-1)*nk3 + (i-1)*nk2*nk3 + 1
	           !  xkg are the components of the complete grid in crystal axis
	           xkg(1,n) = dble(i-1)/nk1 + dble(k1)/2/nk1
	           xkg(2,n) = dble(j-1)/nk2 + dble(k2)/2/nk2
	           xkg(3,n) = dble(k-1)/nk3 + dble(k3)/2/nk3
	        ENDDO
	     ENDDO
	  ENDDO
	
	  !  equiv(nk) =nk : k-point nk is not equivalent to any previous k-point
	  !  equiv(nk)!=nk : k-point nk is equivalent to k-point equiv(nk)
	
	  DO nk=1,nkr
	     equiv(nk)=nk
	  ENDDO
	
	  IF ( skip_equivalence ) THEN
	    WRITE(stdout, '(a)') 'WARNING: skipping k-point equivalence'
	    wkk = 1.d0
	  ELSE
	    DO nk=1,nkr
	    !  check if this k-point has already been found equivalent to another
	      IF (equiv(nk) == nk) THEN
	        wkk(nk)   = 1.0d0
	        !  check if there are equivalent k-point to this in the list
	        !  (excepted those previously found to be equivalent to another)
	        !  check both k and -k
	        DO ns=1,nrot
	           DO i=1,3
	              xkr(i) = s(i,1,ns) * xkg(1,nk) &
	                     + s(i,2,ns) * xkg(2,nk) &
	                     + s(i,3,ns) * xkg(3,nk)
	              xkr(i) = xkr(i) - nint( xkr(i) )
	           ENDDO
	           IF(t_rev(ns)==1) xkr = -xkr
	           xx = xkr(1)*nk1 - 0.5d0*k1
	           yy = xkr(2)*nk2 - 0.5d0*k2
	           zz = xkr(3)*nk3 - 0.5d0*k3
	           in_the_list = abs(xx-nint(xx))<=eps .and. &
	                         abs(yy-nint(yy))<=eps .and. &
	                         abs(zz-nint(zz))<=eps
	           IF (in_the_list) THEN
	              i = mod ( nint ( xkr(1)*nk1 - 0.5d0*k1 + 2*nk1), nk1 ) + 1
	              j = mod ( nint ( xkr(2)*nk2 - 0.5d0*k2 + 2*nk2), nk2 ) + 1
	              k = mod ( nint ( xkr(3)*nk3 - 0.5d0*k3 + 2*nk3), nk3 ) + 1
	              n = (k-1) + (j-1)*nk3 + (i-1)*nk2*nk3 + 1
	              IF (n>nk .and. equiv(n)==n) THEN
	                 equiv(n) = nk
	                 wkk(nk)=wkk(nk)+1.0d0
	              ELSE
	                 IF (equiv(n)/=nk .or. n<nk ) THEN
						WRITE(stdout, '(a)') 'ERROR : with kpoint checking algorithm'
						STOP
					 ENDIF
	              ENDIF
	           ENDIF
	           IF ( time_reversal ) THEN
	              xx =-xkr(1)*nk1 - 0.5d0*k1
	              yy =-xkr(2)*nk2 - 0.5d0*k2
	              zz =-xkr(3)*nk3 - 0.5d0*k3
	              in_the_list=abs(xx-nint(xx))<=eps.and.abs(yy-nint(yy))<=eps &
	                                                 .and. abs(zz-nint(zz))<=eps
	              IF (in_the_list) THEN
	                 i = mod ( nint (-xkr(1)*nk1 - 0.5d0 * k1 + 2*nk1), nk1 ) + 1
	                 j = mod ( nint (-xkr(2)*nk2 - 0.5d0 * k2 + 2*nk2), nk2 ) + 1
	                 k = mod ( nint (-xkr(3)*nk3 - 0.5d0 * k3 + 2*nk3), nk3 ) + 1
	                 n = (k-1) + (j-1)*nk3 + (i-1)*nk2*nk3 + 1
	                 IF (n>nk .and. equiv(n)==n) THEN
	                    equiv(n) = nk
	                    wkk(nk)=wkk(nk)+1.0d0
	                 ELSE
	                    IF (equiv(n)/=nk.or.n<nk) THEN
							WRITE(stdout, '(a)') 'ERROR : with kpoint checking algorithm'
							STOP
						ENDIF
	                 ENDIF
	              ENDIF
	           ENDIF
	        ENDDO
	      ENDIF
	    ENDDO
	  ENDIF
	
	  !  count irreducible points and order them
	
	  nks=0
	  fact=0.0d0
	  DO nk=1,nkr
	     IF (equiv(nk)==nk) THEN
	        nks=nks+1
	        IF (nks>npk) THEN
				WRITE(stdout, '(a)') 'ERROR: too many k-points'
				STOP
			ENDIF
	        wk(nks) = wkk(nk)
	        fact    = fact+wk(nks)
	        !  bring back into to the first BZ
	        DO i=1,3
	           xk(i,nks) = xkg(i,nk)-nint(xkg(i,nk))
	        ENDDO
	     ENDIF
	  ENDDO
	  !  go to cartesian axis (in units 2pi/a0)
	  CALL cryst_to_cart(nks,xk,bg,1)
	  !  normalize weights to one
	  DO nk=1,nks
	     wk(nk) = wk(nk)/fact
	  ENDDO
	
	  DEALLOCATE(equiv)
	  DEALLOCATE(xkg,wkk)
	
	  RETURN
  END SUBROUTINE kpoint_grid
  
	!-----------------------------------------------------------------------
	subroutine cryst_to_cart (nvec, vec, trmat, iflag)
	!-----------------------------------------------------------------------
	  !
	  !     This routine transforms the atomic positions or the k-point
	  !     components from crystallographic to cartesian coordinates 
	  !     ( iflag=1 ) and viceversa ( iflag=-1 ).
	  !     Output cartesian coordinates are stored in the input ('vec') array
	  !
	  !
	  USE kinds
	  implicit none
	  !
	  integer, intent(in) :: nvec, iflag
	  ! nvec:  number of vectors (atomic positions or k-points)
	  !        to be transformed from crystal to cartesian and vice versa
	  ! iflag: gives the direction of the transformation
	  real(RP), intent(in) :: trmat (3, 3)
	  ! trmat: transformation matrix
	  ! if iflag=1:
	  !    trmat = at ,  basis of the real-space lattice,       for atoms   or
	  !          = bg ,  basis of the reciprocal-space lattice, for k-points
	  ! if iflag=-1: the opposite
	  real(RP), intent(inout) :: vec (3, nvec)
	  ! coordinates of the vector (atomic positions or k-points) to be
	  ! transformed - overwritten on output
	  !
	  !    local variables
	  !
	  integer :: nv, kpol
	  ! counter on vectors
	  ! counter on polarizations
	  real(RP) :: vau (3)
	  ! workspace
	  !
	  !     Compute the cartesian coordinates of each vectors
	  !     (atomic positions or k-points components)
	  !
	  do nv = 1, nvec
	     if (iflag.eq.1) then
	        do kpol = 1, 3
	           vau (kpol) = trmat (kpol, 1) * vec (1, nv) + trmat (kpol, 2) &
	                * vec (2, nv) + trmat (kpol, 3) * vec (3, nv)
	        enddo
	     else
	        do kpol = 1, 3
	           vau (kpol) = trmat (1, kpol) * vec (1, nv) + trmat (2, kpol) &
	                * vec (2, nv) + trmat (3, kpol) * vec (3, nv)
	        enddo
	     endif
	     do kpol = 1, 3
	        vec (kpol, nv) = vau (kpol)
	     enddo
	  enddo
	  !
	  return
	end subroutine cryst_to_cart
  
END MODULE monkhorstpack
