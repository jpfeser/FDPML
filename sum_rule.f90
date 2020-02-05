MODULE sum_rule

CONTAINS

!----------------------------------------------------------------------
SUBROUTINE set_asr (asr, nr1, nr2, nr3, frc, zeu, nat, ibrav, tau)
  !-----------------------------------------------------------------------
  !
  USE kinds,      ONLY : RP
  !
  IMPLICIT NONE
  CHARACTER (LEN=10), intent(in) :: asr
  INTEGER, intent(in) :: nr1, nr2, nr3, nat, ibrav
  REAL(RP), intent(in) :: tau(3,nat)
  REAL(RP), intent(inout) :: frc(nr1,nr2,nr3,3,3,nat,nat), zeu(3,3,nat)
  !
  INTEGER :: axis, n, i, j, na, nb, n1,n2,n3, m,p,k,l,q,r, i1,j1,na1
  REAL(RP) :: zeu_new(3,3,nat)
  REAL(RP), ALLOCATABLE :: frc_new(:,:,:,:,:,:,:)
  type vector
     real(RP),pointer :: vec(:,:,:,:,:,:,:)
  end type vector
  !
  type (vector) u(6*3*nat)
  ! These are the "vectors" associated with the sum rules on force-constants
  !
  integer :: u_less(6*3*nat),n_less,i_less
  ! indices of the vectors u that are not independent to the preceding ones,
  ! n_less = number of such vectors, i_less = temporary parameter
  !
  integer, allocatable :: ind_v(:,:,:)
  real(RP), allocatable :: v(:,:)
  ! These are the "vectors" associated with symmetry conditions, coded by
  ! indicating the positions (i.e. the seven indices) of the non-zero elements (there
  ! should be only 2 of them) and the value of that element. We do so in order
  ! to limit the amount of memory used.
  !
  real(RP), allocatable :: w(:,:,:,:,:,:,:), x(:,:,:,:,:,:,:)
  ! temporary vectors and parameters
  real(RP) :: scal,norm2, sum
  !
  real(RP) :: zeu_u(6*3,3,3,nat)
  ! These are the "vectors" associated with the sum rules on effective charges
  !
  integer :: zeu_less(6*3),nzeu_less,izeu_less
  ! indices of the vectors zeu_u that are not independent to the preceding ones,
  ! nzeu_less = number of such vectors, izeu_less = temporary parameter
  !
  real(RP) :: zeu_w(3,3,nat), zeu_x(3,3,nat)
  ! temporary vectors

  ! Initialization. n is the number of sum rules to be considered (if asr.ne.'simple')
  ! and 'axis' is the rotation axis in the case of a 1D system
  ! (i.e. the rotation axis is (Ox) if axis='1', (Oy) if axis='2' and (Oz) if axis='3')
  !
  if((asr.ne.'simple').and.(asr.ne.'crystal').and.(asr.ne.'one-dim') &
                      .and.(asr.ne.'zero-dim')) then
     write(*,*) 'ERROR '
     write(*,*) 'Error in Acoustic Sum Rule'
     write(*,*) 'Incorrect acoustic sum rule'
     STOP
  endif
  !
  if(asr.eq.'simple') then
     !
     ! Simple Acoustic Sum Rule on effective charges
     !
     do i=1,3
        do j=1,3
           sum=0.0d0
           do na=1,nat
              sum = sum + zeu(i,j,na)
           end do
           do na=1,nat
              zeu(i,j,na) = zeu(i,j,na) - sum/nat
           end do
        end do
     end do
     !
     ! Simple Acoustic Sum Rule on force constants in real space
     !
     do i=1,3
        do j=1,3
           do na=1,nat
              sum=0.0d0
               do nb=1,nat
                  do n1=1,nr1
                     do n2=1,nr2
                        do n3=1,nr3
                           sum=sum+frc(n1,n2,n3,i,j,na,nb)
                        end do
                     end do
                  end do
               end do
               frc(1,1,1,i,j,na,na) = frc(1,1,1,i,j,na,na) - sum
               !               write(6,*) ' na, i, j, sum = ',na,i,j,sum
            end do
         end do
      end do
      !
      return
      !
   end if

  if(asr.eq.'crystal') n=3
  if(asr.eq.'one-dim') then
     ! the direction of periodicity is the rotation axis
     ! It will work only if the crystal axis considered is one of
     ! the cartesian axis (typically, ibrav=1, 6 or 8, or 4 along the
     ! z-direction)
     if (nr1*nr2*nr3.eq.1) axis=3
     if ((nr1.ne.1).and.(nr2*nr3.eq.1)) axis=1
     if ((nr2.ne.1).and.(nr1*nr3.eq.1)) axis=2
     if ((nr3.ne.1).and.(nr1*nr2.eq.1)) axis=3
     if (((nr1.ne.1).and.(nr2.ne.1)).or.((nr2.ne.1).and. &
          (nr3.ne.1)).or.((nr1.ne.1).and.(nr3.ne.1))) then
			write(*,*) 'ERROR '
			write(*,*) 'Error in Acoustic Sum Rule'
			write(*,*) 'Too many directions of periodicity in 1D system'
			STOP
     endif
     if ((ibrav.ne.1).and.(ibrav.ne.6).and.(ibrav.ne.8).and. &
          ((ibrav.ne.4).or.(axis.ne.3)) ) then
        write(*,*) 'asr: rotational axis may be wrong'
     endif
     write(*,'("asr rotation axis in 1D system= ",I4)') axis
     n=4
  endif
  if(asr.eq.'zero-dim') n=6
  !
  ! Acoustic Sum Rule on effective charges
  !
  ! generating the vectors of the orthogonal of the subspace to project
  ! the effective charges matrix on
  !
  zeu_u(:,:,:,:)=0.0d0
  do i=1,3
     do j=1,3
        do na=1,nat
           zeu_new(i,j,na)=zeu(i,j,na)
        enddo
     enddo
  enddo
  !
  p=0
  do i=1,3
     do j=1,3
        ! These are the 3*3 vectors associated with the
        ! translational acoustic sum rules
        p=p+1
        zeu_u(p,i,j,:)=1.0d0
        !
     enddo
  enddo
  !
  if (n.eq.4) then
     do i=1,3
        ! These are the 3 vectors associated with the
        ! single rotational sum rule (1D system)
        p=p+1
        do na=1,nat
           zeu_u(p,i,MOD(axis,3)+1,na)=-tau(MOD(axis+1,3)+1,na)
           zeu_u(p,i,MOD(axis+1,3)+1,na)=tau(MOD(axis,3)+1,na)
        enddo
        !
     enddo
  endif
  !
  if (n.eq.6) then
     do i=1,3
        do j=1,3
           ! These are the 3*3 vectors associated with the
           ! three rotational sum rules (0D system - typ. molecule)
           p=p+1
           do na=1,nat
              zeu_u(p,i,MOD(j,3)+1,na)=-tau(MOD(j+1,3)+1,na)
              zeu_u(p,i,MOD(j+1,3)+1,na)=tau(MOD(j,3)+1,na)
           enddo
           !
        enddo
     enddo
  endif
  !
  ! Gram-Schmidt orthonormalization of the set of vectors created.
  !
  nzeu_less=0
  do k=1,p
     zeu_w(:,:,:)=zeu_u(k,:,:,:)
     zeu_x(:,:,:)=zeu_u(k,:,:,:)
     do q=1,k-1
        r=1
        do izeu_less=1,nzeu_less
           if (zeu_less(izeu_less).eq.q) r=0
        enddo
        if (r.ne.0) then
           call sp_zeu(zeu_x,zeu_u(q,:,:,:),nat,scal)
           zeu_w(:,:,:) = zeu_w(:,:,:) - scal* zeu_u(q,:,:,:)
        endif
     enddo
     call sp_zeu(zeu_w,zeu_w,nat,norm2)
     if (norm2.gt.1.0d-16) then
        zeu_u(k,:,:,:) = zeu_w(:,:,:) / DSQRT(norm2)
     else
        nzeu_less=nzeu_less+1
        zeu_less(nzeu_less)=k
     endif
  enddo
  !
  ! Projection of the effective charge "vector" on the orthogonal of the
  ! subspace of the vectors verifying the sum rules
  !
  zeu_w(:,:,:)=0.0d0
  do k=1,p
     r=1
     do izeu_less=1,nzeu_less
        if (zeu_less(izeu_less).eq.k) r=0
     enddo
     if (r.ne.0) then
        zeu_x(:,:,:)=zeu_u(k,:,:,:)
        call sp_zeu(zeu_x,zeu_new,nat,scal)
        zeu_w(:,:,:) = zeu_w(:,:,:) + scal*zeu_u(k,:,:,:)
     endif
  enddo
  !
  ! Final substraction of the former projection to the initial zeu, to get
  ! the new "projected" zeu
  !
  zeu_new(:,:,:)=zeu_new(:,:,:) - zeu_w(:,:,:)
  call sp_zeu(zeu_w,zeu_w,nat,norm2)
  write(*,'("Norm of the difference between old and new effective ", &
       & "charges: ",F25.20)') SQRT(norm2)
  !
  ! Check projection
  !
  !write(6,'("Check projection of zeu")')
  !do k=1,p
  !  zeu_x(:,:,:)=zeu_u(k,:,:,:)
  !  call sp_zeu(zeu_x,zeu_new,nat,scal)
  !  if (DABS(scal).gt.1d-10) write(6,'("k= ",I8," zeu_new|zeu_u(k)= ",F15.10)') k,scal
  !enddo
  !
  do i=1,3
     do j=1,3
        do na=1,nat
           zeu(i,j,na)=zeu_new(i,j,na)
        enddo
     enddo
  enddo
  !
  ! Acoustic Sum Rule on force constants
  !
  !
  ! generating the vectors of the orthogonal of the subspace to project
  ! the force-constants matrix on
  !
  do k=1,18*nat
     allocate(u(k) % vec(nr1,nr2,nr3,3,3,nat,nat))
     u(k) % vec (:,:,:,:,:,:,:)=0.0d0
  enddo
  ALLOCATE (frc_new(nr1,nr2,nr3,3,3,nat,nat))
  do i=1,3
     do j=1,3
        do na=1,nat
           do nb=1,nat
              do n1=1,nr1
                 do n2=1,nr2
                    do n3=1,nr3
                       frc_new(n1,n2,n3,i,j,na,nb)=frc(n1,n2,n3,i,j,na,nb)
                    enddo
                 enddo
              enddo
           enddo
        enddo
     enddo
  enddo
  !
  p=0
  do i=1,3
     do j=1,3
        do na=1,nat
           ! These are the 3*3*nat vectors associated with the
           ! translational acoustic sum rules
           p=p+1
           u(p) % vec (:,:,:,i,j,na,:)=1.0d0
           !
        enddo
     enddo
  enddo
  !
  if (n.eq.4) then
     do i=1,3
        do na=1,nat
           ! These are the 3*nat vectors associated with the
           ! single rotational sum rule (1D system)
           p=p+1
           do nb=1,nat
              u(p) % vec (:,:,:,i,MOD(axis,3)+1,na,nb)=-tau(MOD(axis+1,3)+1,nb)
              u(p) % vec (:,:,:,i,MOD(axis+1,3)+1,na,nb)=tau(MOD(axis,3)+1,nb)
           enddo
           !
        enddo
     enddo
  endif
  !
  if (n.eq.6) then
     do i=1,3
        do j=1,3
           do na=1,nat
              ! These are the 3*3*nat vectors associated with the
              ! three rotational sum rules (0D system - typ. molecule)
              p=p+1
              do nb=1,nat
                 u(p) % vec (:,:,:,i,MOD(j,3)+1,na,nb)=-tau(MOD(j+1,3)+1,nb)
                 u(p) % vec (:,:,:,i,MOD(j+1,3)+1,na,nb)=tau(MOD(j,3)+1,nb)
              enddo
              !
           enddo
        enddo
     enddo
  endif
  !
  allocate (ind_v(9*nat*nat*nr1*nr2*nr3,2,7), v(9*nat*nat*nr1*nr2*nr3,2) )
  m=0
  do i=1,3
     do j=1,3
        do na=1,nat
           do nb=1,nat
              do n1=1,nr1
                 do n2=1,nr2
                    do n3=1,nr3
                       ! These are the vectors associated with the symmetry constraints
                       q=1
                       l=1
                       do while((l.le.m).and.(q.ne.0))
                          if ((ind_v(l,1,1).eq.n1).and.(ind_v(l,1,2).eq.n2).and. &
                               (ind_v(l,1,3).eq.n3).and.(ind_v(l,1,4).eq.i).and. &
                               (ind_v(l,1,5).eq.j).and.(ind_v(l,1,6).eq.na).and. &
                               (ind_v(l,1,7).eq.nb)) q=0
                          if ((ind_v(l,2,1).eq.n1).and.(ind_v(l,2,2).eq.n2).and. &
                               (ind_v(l,2,3).eq.n3).and.(ind_v(l,2,4).eq.i).and. &
                               (ind_v(l,2,5).eq.j).and.(ind_v(l,2,6).eq.na).and. &
                               (ind_v(l,2,7).eq.nb)) q=0
                          l=l+1
                       enddo
                       if ((n1.eq.MOD(nr1+1-n1,nr1)+1).and.(n2.eq.MOD(nr2+1-n2,nr2)+1) &
                            .and.(n3.eq.MOD(nr3+1-n3,nr3)+1).and.(i.eq.j).and.(na.eq.nb)) q=0
                       if (q.ne.0) then
                          m=m+1
                          ind_v(m,1,1)=n1
                          ind_v(m,1,2)=n2
                          ind_v(m,1,3)=n3
                          ind_v(m,1,4)=i
                          ind_v(m,1,5)=j
                          ind_v(m,1,6)=na
                          ind_v(m,1,7)=nb
                          v(m,1)=1.0d0/DSQRT(2.0d0)
                          ind_v(m,2,1)=MOD(nr1+1-n1,nr1)+1
                          ind_v(m,2,2)=MOD(nr2+1-n2,nr2)+1
                          ind_v(m,2,3)=MOD(nr3+1-n3,nr3)+1
                          ind_v(m,2,4)=j
                          ind_v(m,2,5)=i
                          ind_v(m,2,6)=nb
                          ind_v(m,2,7)=na
                          v(m,2)=-1.0d0/DSQRT(2.0d0)
                       endif
                    enddo
                 enddo
              enddo
           enddo
        enddo
     enddo
  enddo
  !
  ! Gram-Schmidt orthonormalization of the set of vectors created.
  ! Note that the vectors corresponding to symmetry constraints are already
  ! orthonormalized by construction.
  !
  n_less=0
  allocate (w(nr1,nr2,nr3,3,3,nat,nat), x(nr1,nr2,nr3,3,3,nat,nat))
  do k=1,p
     w(:,:,:,:,:,:,:)=u(k) % vec (:,:,:,:,:,:,:)
     x(:,:,:,:,:,:,:)=u(k) % vec (:,:,:,:,:,:,:)
     do l=1,m
        !
        call sp2(x,v(l,:),ind_v(l,:,:),nr1,nr2,nr3,nat,scal)
        do r=1,2
           n1=ind_v(l,r,1)
           n2=ind_v(l,r,2)
           n3=ind_v(l,r,3)
           i=ind_v(l,r,4)
           j=ind_v(l,r,5)
           na=ind_v(l,r,6)
           nb=ind_v(l,r,7)
           w(n1,n2,n3,i,j,na,nb)=w(n1,n2,n3,i,j,na,nb)-scal*v(l,r)
        enddo
     enddo
     if (k.le.(9*nat)) then
        na1=MOD(k,nat)
        if (na1.eq.0) na1=nat
        j1=MOD((k-na1)/nat,3)+1
        i1=MOD((((k-na1)/nat)-j1+1)/3,3)+1
     else
        q=k-9*nat
        if (n.eq.4) then
           na1=MOD(q,nat)
           if (na1.eq.0) na1=nat
           i1=MOD((q-na1)/nat,3)+1
        else
           na1=MOD(q,nat)
           if (na1.eq.0) na1=nat
           j1=MOD((q-na1)/nat,3)+1
           i1=MOD((((q-na1)/nat)-j1+1)/3,3)+1
        endif
     endif
     do q=1,k-1
        r=1
        do i_less=1,n_less
           if (u_less(i_less).eq.q) r=0
        enddo
        if (r.ne.0) then
           call sp3(x,u(q) % vec (:,:,:,:,:,:,:), i1,na1,nr1,nr2,nr3,nat,scal)
           w(:,:,:,:,:,:,:) = w(:,:,:,:,:,:,:) - scal* u(q) % vec (:,:,:,:,:,:,:)
        endif
     enddo
     call sp1(w,w,nr1,nr2,nr3,nat,norm2)
     if (norm2.gt.1.0d-16) then
        u(k) % vec (:,:,:,:,:,:,:) = w(:,:,:,:,:,:,:) / DSQRT(norm2)
     else
        n_less=n_less+1
        u_less(n_less)=k
     endif
  enddo
  !
  ! Projection of the force-constants "vector" on the orthogonal of the
  ! subspace of the vectors verifying the sum rules and symmetry contraints
  !
  w(:,:,:,:,:,:,:)=0.0d0
  do l=1,m
     call sp2(frc_new,v(l,:),ind_v(l,:,:),nr1,nr2,nr3,nat,scal)
     do r=1,2
        n1=ind_v(l,r,1)
        n2=ind_v(l,r,2)
        n3=ind_v(l,r,3)
        i=ind_v(l,r,4)
        j=ind_v(l,r,5)
        na=ind_v(l,r,6)
        nb=ind_v(l,r,7)
        w(n1,n2,n3,i,j,na,nb)=w(n1,n2,n3,i,j,na,nb)+scal*v(l,r)
     enddo
  enddo
  do k=1,p
     r=1
     do i_less=1,n_less
        if (u_less(i_less).eq.k) r=0
     enddo
     if (r.ne.0) then
        x(:,:,:,:,:,:,:)=u(k) % vec (:,:,:,:,:,:,:)
        call sp1(x,frc_new,nr1,nr2,nr3,nat,scal)
        w(:,:,:,:,:,:,:) = w(:,:,:,:,:,:,:) + scal*u(k)%vec(:,:,:,:,:,:,:)
     endif
     deallocate(u(k) % vec)
  enddo
  !
  ! Final substraction of the former projection to the initial frc, to get
  ! the new "projected" frc
  !
  frc_new(:,:,:,:,:,:,:)=frc_new(:,:,:,:,:,:,:) - w(:,:,:,:,:,:,:)
  call sp1(w,w,nr1,nr2,nr3,nat,norm2)
  write(*,'("Norm of the difference between old and new force-constants:",&
       &     F25.20)') SQRT(norm2)
  !
  ! Check projection
  !
  !write(6,'("Check projection IFC")')
  !do l=1,m
  !  call sp2(frc_new,v(l,:),ind_v(l,:,:),nr1,nr2,nr3,nat,scal)
  !  if (DABS(scal).gt.1d-10) write(6,'("l= ",I8," frc_new|v(l)= ",F15.10)') l,scal
  !enddo
  !do k=1,p
  !  x(:,:,:,:,:,:,:)=u(k) % vec (:,:,:,:,:,:,:)
  !  call sp1(x,frc_new,nr1,nr2,nr3,nat,scal)
  !  if (DABS(scal).gt.1d-10) write(6,'("k= ",I8," frc_new|u(k)= ",F15.10)') k,scal
  !  deallocate(u(k) % vec)
  !enddo
  !
  do i=1,3
     do j=1,3
        do na=1,nat
           do nb=1,nat
              do n1=1,nr1
                 do n2=1,nr2
                    do n3=1,nr3
                       frc(n1,n2,n3,i,j,na,nb)=frc_new(n1,n2,n3,i,j,na,nb)
                    enddo
                 enddo
              enddo
           enddo
        enddo
     enddo
  enddo
  deallocate (x, w)
  deallocate (v, ind_v)
  deallocate (frc_new)
  !
  return
end subroutine set_asr
!
!----------------------------------------------------------------------
subroutine sp_zeu(zeu_u,zeu_v,nat,scal)
  !-----------------------------------------------------------------------
  !
  ! does the scalar product of two effective charges matrices zeu_u and zeu_v
  ! (considered as vectors in the R^(3*3*nat) space, and coded in the usual way)
  !
  USE kinds, ONLY: RP
  implicit none
  integer i,j,na,nat
  real(RP) zeu_u(3,3,nat)
  real(RP) zeu_v(3,3,nat)
  real(RP) scal
  !
  !
  scal=0.0d0
  do i=1,3
    do j=1,3
      do na=1,nat
        scal=scal+zeu_u(i,j,na)*zeu_v(i,j,na)
      enddo
    enddo
  enddo
  !
  return
  !
end subroutine sp_zeu
!
!
!----------------------------------------------------------------------
subroutine sp1(u,v,nr1,nr2,nr3,nat,scal)
  !-----------------------------------------------------------------------
  !
  ! does the scalar product of two force-constants matrices u and v (considered as
  ! vectors in the R^(3*3*nat*nat*nr1*nr2*nr3) space, and coded in the usual way)
  !
  USE kinds, ONLY: RP
  implicit none
  integer nr1,nr2,nr3,i,j,na,nb,n1,n2,n3,nat
  real(RP) u(nr1,nr2,nr3,3,3,nat,nat)
  real(RP) v(nr1,nr2,nr3,3,3,nat,nat)
  real(RP) scal
  !
  !
  scal=0.0d0
  do i=1,3
    do j=1,3
      do na=1,nat
        do nb=1,nat
          do n1=1,nr1
            do n2=1,nr2
              do n3=1,nr3
                scal=scal+u(n1,n2,n3,i,j,na,nb)*v(n1,n2,n3,i,j,na,nb)
              enddo
            enddo
          enddo
        enddo
      enddo
    enddo
  enddo
  !
  return
  !
end subroutine sp1
!
!----------------------------------------------------------------------
subroutine sp2(u,v,ind_v,nr1,nr2,nr3,nat,scal)
  !-----------------------------------------------------------------------
  !
  ! does the scalar product of two force-constants matrices u and v (considered as
  ! vectors in the R^(3*3*nat*nat*nr1*nr2*nr3) space). u is coded in the usual way
  ! but v is coded as explained when defining the vectors corresponding to the
  ! symmetry constraints
  !
  USE kinds, ONLY: RP
  implicit none
  integer nr1,nr2,nr3,i,nat
  real(RP) u(nr1,nr2,nr3,3,3,nat,nat)
  integer ind_v(2,7)
  real(RP) v(2)
  real(RP) scal
  !
  !
  scal=0.0d0
  do i=1,2
    scal=scal+u(ind_v(i,1),ind_v(i,2),ind_v(i,3),ind_v(i,4),ind_v(i,5),ind_v(i,6), &
         ind_v(i,7))*v(i)
  enddo
  !
  return
  !
end subroutine sp2
!
!----------------------------------------------------------------------
subroutine sp3(u,v,i,na,nr1,nr2,nr3,nat,scal)
  !-----------------------------------------------------------------------
  !
  ! like sp1, but in the particular case when u is one of the u(k)%vec
  ! defined in set_asr (before orthonormalization). In this case most of the
  ! terms are zero (the ones that are not are characterized by i and na), so
  ! that a lot of computer time can be saved (during Gram-Schmidt).
  !
  USE kinds, ONLY: RP
  implicit none
  integer nr1,nr2,nr3,i,j,na,nb,n1,n2,n3,nat
  real(RP) u(nr1,nr2,nr3,3,3,nat,nat)
  real(RP) v(nr1,nr2,nr3,3,3,nat,nat)
  real(RP) scal
  !
  !
  scal=0.0d0
  do j=1,3
    do nb=1,nat
      do n1=1,nr1
        do n2=1,nr2
          do n3=1,nr3
            scal=scal+u(n1,n2,n3,i,j,na,nb)*v(n1,n2,n3,i,j,na,nb)
          enddo
        enddo
      enddo
    enddo
  enddo
  !
  return
  !
end subroutine sp3

END MODULE sum_rule
