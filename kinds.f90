MODULE kinds

use, intrinsic :: iso_fortran_env, ONLY : int32, int64

IMPLICIT NONE

INTEGER, PARAMETER :: RP = int64, IP = int64, CP = int64
REAL(KIND = RP), PARAMETER :: eps = EPSILON(0.0)

END MODULE kinds
