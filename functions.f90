module functions
use parameters
implicit none

contains
! Returns the inverse of a complex GENERAL matrix calculated by finding the LU
! decomposition.  Depends on LAPACK. (adapted by AR from http://fortranwiki.org/fortran/show/inv)
function inv(A) result(Ainv)
  implicit none
  complex(kind=dp), dimension(:,:), intent(in) :: A
  complex(kind=dp), dimension(size(A,1),size(A,2)) :: Ainv

  complex(kind=dp), dimension(size(A,1)) :: work  ! work array for LAPACK
  integer, dimension(size(A,1)) :: ipiv   ! pivot indices
  integer :: n, info

  ! External procedures defined in LAPACK
  external ZGETRF
  external ZGETRI

  ! Store A in Ainv to prevent it from being overwritten by LAPACK
  Ainv = A
  n = size(A,1)

  ! ZGETRF computes an LU factorization of a general M-by-N matrix A
  ! using partial pivoting with row interchanges.
  call ZGETRF(n, n, Ainv, n, ipiv, info)

  if (info /= 0) then
     stop 'Matrix is numerically singular (ZGETRF)!'
  end if

  ! ZGETRI computes the inverse of a matrix using the LU factorization
  ! computed by ZGETRF.
  call ZGETRI(n, Ainv, n, ipiv, work, n, info)

  if (info /= 0) then
     stop 'Matrix inversion failed (ZGETRI)!'
  end if
end function inv
!
!
! Returns two arrays of eigenvalues and eigenvectors of a complex GENERAL matrix M
!
subroutine eigeng(M,eval,evec)
  implicit none
  complex(kind=dp), dimension(:,:), intent(in) :: M   ! our matrix
  complex(kind=dp), dimension(size(M,1),size(M,1)) :: A   ! our matrix (working array)
  complex(kind=dp), dimension(size(M,1)), intent(out) :: eval  ! eigenvalues
  complex(kind=dp), dimension(size(M,1),size(M,1)), intent(out) :: evec ! eigenvectors (second index numerates different eigenvalues)
  real(kind=dp), dimension(2*size(M,1)) :: rwork  ! work array for LAPACK
  complex(kind=dp), dimension(:), allocatable :: work  ! work array for LAPACK (optimal dimension is determined by ZGEEV)
  integer :: n, info, m1,m2
  !
  !
  if (size(M,1) /= size(M,2)) then
     stop 'Problem finding eigenvalues/eigenvectors (non-square matrix)!'
  endif
  !
  ! Store M in A to prevent it from being overwritten
  A = M
  !
  allocate(work(1))
  !
  call ZGEEV('N','V',size(M,1),A,size(M,1),eval,'N',1,evec,size(M,1),work,-1,rwork,info) ! calculation of the optimal dimension of the work array
  n = nint(real(work(1)))
  !
  deallocate(work)
  !
  if (info /= 0) then
     stop 'Problem finding optimal dimension of the work array (ZGEEV)!'
  end if
  !
  allocate(work(n))
  !
  call ZGEEV('N','V',size(M,1),A,size(M,1),eval,'N',1,evec,size(M,1),work,n,rwork,info) ! normal action of the routine
  !
  deallocate(work)
  !
  if (info /= 0) then
     stop 'Problem finding eigenvalues/eigenvectors (ZGEEV)!'
  end if
  !
  !
end subroutine eigeng

! Returns two arrays of eigenvalues and eigenvectors of a complex HERMITIAN matrix M
!
subroutine eigen(M,eval,A)
  implicit none
  complex(kind=dp), dimension(:,:), intent(in) :: M   ! our matrix
  real(kind=dp), dimension(size(M,1)), intent(out) :: eval  ! eigenvalues
  complex(kind=dp), dimension(size(M,1),size(M,1)), intent(out) :: A ! in: matrix, out: eigenvectors (second index numerates different eigenvalues)
  real(kind=dp), dimension(3*size(M,1)-2) :: rwork  ! work array for LAPACK
  complex(kind=dp), dimension(:), allocatable :: work  ! work array for LAPACK (optimal dimension is determined by ZHEEV)
  integer :: lwork, info
  !
  !
  if (size(M,1) /= size(M,2)) then
     stop 'Problem finding eigenvalues/eigenvectors (non-square matrix)!'
  endif
  !
  !
  ! Store M in A to prevent it from being overwritten
  A = M
  !
  !
  allocate(work(1))
  !
  call ZHEEV('V','U',size(M,1),A,size(M,1),eval,work,-1,rwork,info) ! calculation of the optimal dimension of the work array
  lwork = nint(real(work(1)))
  !
  deallocate(work)
  !
  if (info /= 0) then
     stop 'Problem finding optimal dimension of the work array (ZHEEV)!'
  end if
  !
  !
  allocate(work(lwork))
  !
  call ZHEEV('V','U',size(M,1),A,size(M,1),eval,work,lwork,rwork,info) ! normal action of the routine
  !
  deallocate(work)
  !
  if (info /= 0) then
     stop 'Problem finding eigenvalues/eigenvectors (ZHEEV)!'
  end if
  !
  !
end subroutine eigen

! GF by matrix inversion
subroutine green_inv(eta,epsil,Hk,Gk)
  !
  use parameters, only : ci,cmplx_0
  implicit none
  !
  complex(kind=dp), dimension(:,:), intent(in) :: Hk
  complex(kind=dp), dimension(size(Hk,1),size(Hk,2)), intent(out) :: Gk
  real(kind=dp), intent(in) :: eta, epsil
  integer :: m1, m2, numwann
  complex(kind=dp), dimension(size(Hk,1),size(Hk,2)) ::  Gk_tmp
  ! LAPACK variables
  integer :: info
  complex(kind=dp), dimension(size(Hk,1)) :: work  ! work array for LAPACK
  integer, dimension(size(Hk,1)) :: ipiv   ! pivot indices
  !
  numwann = size(Hk,1)
  !
  Gk_tmp(:,:) = cmplx_0
  Gk(:,:) = cmplx_0
  !
  do m1=1,numwann
  do m2=1,numwann
     !
     if (m1==m2) then   ! diagonal elements
        Gk_tmp(m1,m2) = ci*eta + epsil - Hk(m1,m2)  ! Denominator of the Green's function
     else               ! nondiagonal elements
        Gk_tmp(m1,m2) = -Hk(m1,m2)
     endif
     !
  enddo
  enddo
  !
  call ZGETRF(numwann, numwann, Gk_tmp, numwann, ipiv, info)
  if (info /= 0) stop 'ERROR in LU matrix factorization (ZGETRF)!'

  call ZGETRI(numwann, Gk_tmp, numwann, ipiv, work, numwann, info)
  if (info /= 0) stop 'ERROR in matrix inversion (ZGETRI)!'
  !
  Gk(:,:) = Gk_tmp(:,:)
  !
end subroutine green_inv

! Calculates the mn-component of a (retarded) Green's function from a given eigenvalues and eigenvectors of a Hamiltonian.
! (spectral representation of GF)
subroutine green_spec(eta,eps,G,kpt,up)
  !
  use parameters, only : ci,cmplx_0
  use parameters, only : nwann
  use parameters, only : eval_up, eval_dn, evec_up, evec_dn   ! previously calculated!
  implicit none
  !
!  complex(kind=dp), dimension(:,:) :: eigvec 
!  real(kind=dp), dimension(:) :: eigval  
  integer, intent(in) :: kpt
  real(kind=dp), intent(in) :: eta, eps
  complex(kind=dp), dimension(nwann,nwann), intent(out) :: G
  logical, intent(in) :: up
  !
  integer :: l, m1, m2
  complex(kind=dp), dimension(nwann) :: vec1, vec2
  complex(kind=dp) :: cdotc, denom
  !
  !
!  complex(dp), dimension(:,:), allocatable      ::  eigval_k ! eigenvalues and eigenvectors of the Hamiltonian in k-space
!  complex(dp), dimension(:,:,:), allocatable    ::  eigvec_k ! eigenvectors of the Hamiltonian in k-space
!  call eigen(ham_k(:,:,loop_kpt),eigval_k(:,loop_kpt),eigvec_k(:,:,loop_kpt))  ! Store eigval and eigvec of Hij(k) for all k-points
!  write(*,*) eigval_k(:,1) ! Eigenvalues at the Gamma-point                    ! (in order to avoid their recalculation in future)
  !
  G(:,:) = cmplx_0
  vec1(:) = cmplx_0
  vec2(:) = cmplx_0
  !
  do m1=1,nwann
  do m2=1,nwann
     if (up .EQV. .true.) then
        vec1(:) = evec_up(m2,:,kpt)
        vec2(:) = evec_up(m1,:,kpt)/(ci*eta + eps - eval_up(:,kpt))
        G(m1,m2) = dot_product(vec1(:),vec2(:))                    ! vec1 is conjugated automatically by dot_product since vec1&vec2 are complex!
!
!        do l=1,nwann
!            cdotc = evec_up(m1,l,kpt)*conjg(evec_up(m2,l,kpt))  ! numerator
!            denom = ci*eta + eps - eval_up(l,kpt) ! denominator
!            G(m1,m2) = G(m1,m2) + cdotc/denom
!         enddo
     else
        vec1(:) = evec_dn(m2,:,kpt)
        vec2(:) = evec_dn(m1,:,kpt)/(ci*eta + eps - eval_dn(:,kpt))
        G(m1,m2) = dot_product(vec1(:),vec2(:))            
!         do l=1,nwann
!            cdotc = evec_dn(m1,l,kpt)*conjg(evec_dn(m2,l,kpt))  ! numerator
!            denom = ci*eta + eps - eval_dn(l,kpt) ! denominator
!            G(m1,m2) = G(m1,m2) + cdotc/denom
!         enddo
     endif
  enddo
  enddo
  !
  !
end subroutine green_spec

! This is a simple function to search for an available unit.
! LUN_MIN and LUN_MAX define the range of possible LUNs to check.
! The UNIT value is returned by the function, and also by the optional
! argument. This allows the function to be used directly in an OPEN
! statement, and optionally save the result in a local variable.
! If no units are available, -1 is returned.
! This function is taken from fortranwiki.org (AR)
integer function newunit(unit)
  integer, intent(out), optional :: unit
! local
  integer, parameter :: LUN_MIN=10, LUN_MAX=10000
  logical :: opened
  integer :: lun
! begin
  newunit=-1
  do lun=LUN_MIN,LUN_MAX
    inquire(unit=lun,opened=opened)
    if (.not. opened) then
      newunit=lun
      exit
    end if
  end do
  if (present(unit)) unit=newunit
end function newunit

end module functions
