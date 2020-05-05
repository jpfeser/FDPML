MODULE matrix_inversion

   use, intrinsic :: iso_fortran_env, only : stdin=>input_unit, &
											 stdout=>output_unit, &
											 stderr=>error_unit

   IMPLICIT NONE
   PRIVATE
   PUBLIC :: invmat

   INTERFACE invmat
      MODULE PROCEDURE invmat_r, invmat_c
   END INTERFACE

   CONTAINS

  SUBROUTINE invmat_r (n, a, a_inv, da)
  !-----------------------------------------------------------------------
  ! computes the inverse of a n*n real matrix "a" using LAPACK routines 
  ! if "a_inv" is not present, "a" contains the inverse on output
  ! if "a_inv" is present, it contains the inverse on output, "a" is unchanged
  ! if "da" is specified and if the matrix is dimensioned 3x3, 
  ! it also returns the determinant in "da"
  !
  USE kinds, ONLY : RP
  IMPLICIT NONE
  INTEGER, INTENT(in) :: n
  REAL(RP), DIMENSION (n,n), INTENT(inout)  :: a
  REAL(RP), DIMENSION (n,n), INTENT(out), OPTIONAL :: a_inv
  REAL(RP), OPTIONAL, INTENT(out) :: da
  !
  INTEGER :: info, lda, lwork
  ! info=0: inversion was successful
  ! lda   : leading dimension (the same as n)
  INTEGER, ALLOCATABLE :: ipiv (:)
  ! ipiv  : work space for pivoting
  REAL(RP), ALLOCATABLE :: work (:)
  ! more work space
  INTEGER, SAVE :: lworkfact = 64
  !
  IF ( PRESENT(da) ) THEN
     IF ( n == 3 ) THEN
        da = a(1,1)*(a(2,2)*a(3,3)-a(2,3)*a(3,2)) + &
             a(1,2)*(a(2,3)*a(3,1)-a(2,1)*a(3,3)) + &
             a(1,3)*(a(2,1)*a(3,2)-a(3,1)*a(2,2))
        IF (abs(da) < 1.d-10) THEN
			WRITE (stdout, '(a)') 'Matrix is singular'
			STOP
		ENDIF
     ELSE
        da = 0.0_RP
     ENDIF
  ENDIF
  !
  lda = n
  lwork=64*n
  ALLOCATE(ipiv(n), work(lwork) )
  !
  IF ( PRESENT(a_inv) ) THEN
     a_inv(:,:) = a(:,:)
     CALL dgetrf (n, n, a_inv, lda, ipiv, info)
  ELSE
     CALL dgetrf (n, n, a, lda, ipiv, info)
  END IF

  IF ( PRESENT(a_inv) ) THEN
     CALL dgetri (n, a_inv, lda, ipiv, work, lwork, info)
  ELSE
     CALL dgetri (n, a, lda, ipiv, work, lwork, info)
  END IF 

  !
  lworkfact = INT (work(1)/n)
  DEALLOCATE ( work, ipiv )

  END SUBROUTINE invmat_r

  SUBROUTINE invmat_c (n, a, a_inv, da)
  !-----------------------------------------------------------------------
  ! as invmat_r, for a complex matrix
  !
  USE kinds, ONLY : CP
  IMPLICIT NONE
  INTEGER, INTENT(IN) :: n
  COMPLEX (CP), DIMENSION (n,n), INTENT(INOUT)  :: a
  COMPLEX (CP), OPTIONAL, DIMENSION (n,n), INTENT(OUT) :: a_inv
  COMPLEX (CP), OPTIONAL, INTENT(OUT) :: da
  !
  INTEGER :: info, lda, lwork
  ! info=0: inversion was successful
  ! lda   : leading dimension (the same as n)
  INTEGER, ALLOCATABLE :: ipiv (:)
  ! ipiv  : work space for pivoting
  COMPLEX(CP), ALLOCATABLE :: work (:)
  ! more work space
  INTEGER, SAVE :: lworkfact = 64
  !
  IF ( PRESENT(da) ) THEN
     IF (n == 3) THEN
        da = a(1,1)*(a(2,2)*a(3,3)-a(2,3)*a(3,2)) + &
             a(1,2)*(a(2,3)*a(3,1)-a(2,1)*a(3,3)) + &
             a(1,3)*(a(2,1)*a(3,2)-a(3,1)*a(2,2))
        IF (abs(da) < 1.d-10) THEN
			WRITE (stdout, '(a)') 'Matrix is singular'
			STOP
		ENDIF
     ELSE
        da = (0.d0,0.d0)
     ENDIF
  ENDIF
  !
  lda = n
  lwork=64*n
  ALLOCATE(ipiv(n), work(lwork) )
  !
  IF ( PRESENT(a_inv) ) THEN
     a_inv(:,:) = a(:,:)
     CALL zgetrf (n, n, a_inv, lda, ipiv, info)
  ELSE
     CALL zgetrf (n, n, a, lda, ipiv, info)
  END IF

  IF ( PRESENT(a_inv) ) THEN
     CALL zgetri (n, a_inv, lda, ipiv, work, lwork, info)
  ELSE
     CALL zgetri (n, a, lda, ipiv, work, lwork, info)
  END IF

  !
  lworkfact = INT (work(1)/n)
  DEALLOCATE ( work, ipiv )
  !
  END SUBROUTINE invmat_c

END MODULE matrix_inversion
