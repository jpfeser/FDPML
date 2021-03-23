MODULE kinds

use, intrinsic :: iso_fortran_env, ONLY : int32, int64, real32, real64

IMPLICIT NONE

INTEGER, PARAMETER :: RP = int64, IP = real64, CP = real64
REAL(KIND = RP), PARAMETER :: eps = EPSILON(0.0D0)

END MODULE kinds
