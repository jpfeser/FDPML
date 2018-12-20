MODULE dispersion
	
	USE io_module
	USE mp_module
	
	CONTAINS
	
	!----------------------------------------------------------------------------------------------------------------
	SUBROUTINE disp(frc, f_of_q, tau, zeu, m_loc, nr1, nr2, nr3, epsil, nat, ibrav, alat, at, ntyp, ityp, amass, & 
					omega, has_zstar, na_ifc, fd, asr, q, w2, z)
	!----------------------------------------------------------------------------------------------------------------
		USE rigid
		USE lapack95
		USE sum_rule
		USE essentials
		USE ws
		USE kinds
		
		IMPLICIT NONE
		
		CHARACTER(LEN=256) :: asr
		INTEGER, PARAMETER :: nrwsx = 200
		INTEGER :: nrws
		INTEGER :: nr1, nr2, nr3, nat, ibrav, ntyp
		REAL(KIND = RP), DIMENSION(nat) :: amass
		REAL(KIND = RP) :: omega, qq
		LOGICAL :: has_zstar, na_ifc, fd
		REAL(KIND=RP) :: at(3,3), atws(3,3), bg(3,3)
		REAL(KIND=RP) :: q(3), qhat(3)
		COMPLEX(KIND=CP) :: dyn(3,3,nat,nat)
		REAL(KIND = RP) :: f_of_q(3,3,nat,nat)
		INTEGER  :: ityp(nat), itau(nat)
		REAL(KIND=RP) :: alat, qh
		REAL(KIND=RP) :: rws(0:3,nrwsx)
		INTEGER :: it  !counters
		REAL(KIND = RP) :: epsil(3,3)
		REAL(KIND = RP) :: w2(3*nat)
		COMPLEX(KIND = RP) :: z(3*nat, 3*nat)
		REAL(KIND = RP) :: frc(nr1,nr2,nr3,3,3,nat,nat), zeu(3,3,nat), tau(3,nat)
		REAL(KIND = RP), ALLOCATABLE :: m_loc(:,:)
		
		! ityp_blk : atomic types for each atom of the original cell
		
		IF ( ntyp.lt.0 ) THEN 
			WRITE(*,*) 'Error in dispersion relationship'
			WRITE(*,*) 'ntyp < 0'
			WRITE(*,*) 'Check the IFC file.'
		ENDIF
		
		DO it=1,ntyp
			IF (amass(it) < 0.d0) THEN
				WRITE(*,*) 'Error in dispersion relationship'
				WRITE(*,*) 'Incorrect mass in IFC file'
				WRITE(*,*) 'Check the IFC file.'
			END IF
		END DO
		
!		CALL cell_volume(at, alat, omega)
		
		CALL recips(at, bg)
		!
		! build the WS cell corresponding to the force constant grid
		!
		atws(:,1) = at(:,1)*DBLE(nr1)
		atws(:,2) = at(:,2)*DBLE(nr2)
		atws(:,3) = at(:,3)*DBLE(nr3)
		! initialize WS r-vectors
		CALL wsinit(rws,nrwsx,nrws,atws)
		
		IF (asr /= 'no') THEN
			CALL set_asr (asr, nr1, nr2, nr3, frc, zeu, &
				nat, ibrav, tau)
		END IF
		
		dyn(:,:,:,:) = (0.D0, 0.D0)
		
!		lo_to_split = .FALSE.
		
!		f_of_q(:,:,:,:) = CMPLX(0.D0, 0.D0, KIND = RP)
!		IF(na_ifc) THEN

!          qq=sqrt(q(1)**2+q(2)**2+q(3)**3)
!          if(qq == 0.0) qq=1.0
!          qhat(1)=q(1)/qq
!          qhat(2)=q(2)/qq
!          qhat(3)=q(3)/qq

!          CALL nonanal_ifc (nat,nat,ityp,epsil,qhat,zeu,omega, &
!                          nr1, nr2, nr3,f_of_q)
!       END IF
        
        dyn(:,:,:,:) = CMPLX(0.d0,0.d0, KIND =RP)
                
		CALL frc_blk (dyn,q,tau,nat,nr1,nr2,nr3,frc,at,bg,rws,nrws,f_of_q,fd)
		
		!IF (has_zstar .and. .not.na_ifc) &
		!	CALL rgd_blk(nr1,nr2,nr3,nat,dyn,q,tau,epsil,zeu,bg,omega,+1.d0)
			
		!if (io_node) CALL write_dyn_on_file (q, dyn, nat)
						 
		! qhat(1) = q(1)*at(1,1)+q(2)*at(2,1)+q(3)*at(3,1)
        ! qhat(2) = q(1)*at(1,2)+q(2)*at(2,2)+q(3)*at(3,2)
        ! qhat(3) = q(1)*at(1,3)+q(2)*at(2,3)+q(3)*at(3,3)
        ! IF ( ABS( qhat(1) - NINT (qhat(1) ) ) <= eps .AND. &
             ! ABS( qhat(2) - NINT (qhat(2) ) ) <= eps .AND. &
             ! ABS( qhat(3) - NINT (qhat(3) ) ) <= eps ) THEN
           ! !
           ! ! q = 0 : we need the direction q => 0 for the non-analytic part
           ! !
! !           IF ( n == 1 ) THEN
! !              ! if q is the first point in the list
! !              IF ( nq > 1 ) THEN
! !                 ! one more point
! !                 qhat(:) = q(:,n) - q(:,n+1)
! !              ELSE
! !                 ! no more points
! !                 qhat(:) = 0.d0
! !              END IF
! !           ELSE IF ( n > 1 ) THEN
! !              ! if q is not the first point in the list
! !              IF ( q(1,n-1)==0.d0 .AND. &
! !                   q(2,n-1)==0.d0 .AND. &
! !                   q(3,n-1)==0.d0 .AND. n < nq ) THEN
! !                 ! if the preceding q is also 0 :
! !                 qhat(:) = q(:,n) - q(:,n+1)
! !              ELSE
! !                 ! if the preceding q is npt 0 :
! !                 qhat(:) = q(:,n) - q(:,n-1)
! !              END IF
! !           END IF
		   ! qhat(:) = 0.d0 
           ! qh = SQRT(qhat(1)**2+qhat(2)**2+qhat(3)**2)
           ! ! write(*,*) ' qh,  has_zstar ',qh,  has_zstar
           ! IF (qh /= 0.d0) qhat(:) = qhat(:) / qh
           ! IF (qh /= 0.d0 .AND. .NOT. has_zstar) THEN
                ! write(*,*) ' WARNING : Z* not found lo-to splitting will be absent'
           ! ELSE
              ! lo_to_split=.TRUE.
           ! ENDIF
           ! !
           ! CALL nonanal (nat, nat, itau, epsil, qhat, zeu, omega, dyn)
           ! !
        ! END IF
        
        ! if(io_node) call write_dyn_on_file(q,dyn,nat)
        
        CALL dyndiag(nat,ntyp,amass,ityp,dyn,w2,z)
		
	END SUBROUTINE disp
		
	!----------------------------------------------------------------------------
	SUBROUTINE frc_blk(dyn,q,tau,nat,nr1,nr2,nr3,frc,at,bg,rws,nrws,f_of_q,fd)
	!----------------------------------------------------------------------------
	! calculates the dynamical matrix at q from the (short-range part of the)
	! force constants
	!
	USE kinds,      ONLY : RP
	USE constants,  ONLY : tpi
	USE ws
	!
	IMPLICIT NONE
	INTEGER  :: nr1, nr2, nr3, nat, n1, n2, n3, &
			ipol, jpol, na, nb, m1, m2, m3, nint, i,j, nrws, nax
	COMPLEX(KIND = RP) :: dyn(3,3,nat,nat)
	REAL(KIND = RP) :: f_of_q(3,3,nat,nat)
	REAL(KIND = RP) :: frc(nr1,nr2,nr3,3,3,nat,nat), tau(3,nat), q(3), arg, &
				at(3,3), bg(3,3), r(3), weight, r_ws(3),  &
				total_weight, rws(0:3,nrws), alat
	REAL(KIND = RP),SAVE,ALLOCATABLE :: wscache(:,:,:,:,:)
	LOGICAL,SAVE :: first=.true.
	LOGICAL :: fd

	FIRST_TIME : IF (first) THEN
		first=.false.
		ALLOCATE( wscache(-2*nr3:2*nr3, -2*nr2:2*nr2, -2*nr1:2*nr1, nat,nat) )
		DO na=1, nat
			DO nb=1, nat
			total_weight=0.0d0
			!
				DO n1=-2*nr1,2*nr1
					DO n2=-2*nr2,2*nr2
						DO n3=-2*nr3,2*nr3
							DO i=1, 3
								r(i) = n1*at(i,1)+n2*at(i,2)+n3*at(i,3)
								r_ws(i) = r(i) + tau(i,na)-tau(i,nb)
								if (fd) r_ws(i) = r(i) + tau(i,nb)-tau(i,na)
							END DO
								wscache(n3,n2,n1,nb,na) = wsweight(r_ws,rws,nrws)
						ENDDO
					ENDDO
				ENDDO
			ENDDO
		ENDDO
	ENDIF FIRST_TIME
	
	!
	
	DO na=1, nat
		DO nb=1, nat
			total_weight=0.0d0
			DO n1=-2*nr1,2*nr1
				DO n2=-2*nr2,2*nr2
					DO n3=-2*nr3,2*nr3
					!
					! SUM OVER R VECTORS IN THE SUPERCELL - VERY VERY SAFE RANGE!
					!
					DO i=1, 3
						r(i) = n1*at(i,1)+n2*at(i,2)+n3*at(i,3)
					END DO

					weight = wscache(n3,n2,n1,nb,na) 
					IF (weight .GT. 0.0d0) THEN
                    !
                    ! FIND THE VECTOR CORRESPONDING TO R IN THE ORIGINAL CELL
                    !
                    m1 = MOD(n1+1,nr1)
                    IF(m1.LE.0) m1=m1+nr1
                    m2 = MOD(n2+1,nr2)
                    IF(m2.LE.0) m2=m2+nr2
                    m3 = MOD(n3+1,nr3)
                    IF(m3.LE.0) m3=m3+nr3
                
                    arg = tpi*(q(1)*r(1) + q(2)*r(2) + q(3)*r(3))
                    DO ipol=1, 3
                       DO jpol=1, 3
                          dyn(ipol,jpol,na,nb) =                 &
                               dyn(ipol,jpol,na,nb) +            &
                               (frc(m1,m2,m3,ipol,jpol,na,nb)+f_of_q(ipol,jpol,na,nb)	)     &
                               *CMPLX(COS(arg),-SIN(arg),kind=RP)*weight
                       END DO
                    END DO
					END IF ! Weight
					total_weight=total_weight + weight
					END DO
				END DO
			END DO
			IF (ABS(total_weight-DBLE(nr1*nr2*nr3)).GT.1.0d-8) THEN
				write(*,*) ' ERROR'
				write(*,*) ' Error: total weight is incorrect'
				STOP
			END IF
		END DO
	END DO
	
	RETURN
	
	END SUBROUTINE frc_blk

END MODULE
