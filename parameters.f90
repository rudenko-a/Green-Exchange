module parameters

      implicit none
      include "mpif.h"
      integer :: iproc, nprocs, ierr
      !
      integer, parameter :: dp = selected_real_kind(15, 307)
!      integer, parameter :: dp = kind(1.0d0)
      real(dp), parameter :: pi=3.141592653589793238462643383279_dp
      real(dp), parameter :: twopi=2*pi
      real(dp), parameter :: fourpi=4*pi
      complex(dp), parameter :: ci=(0.0_dp,1.0_dp)
      complex(dp), parameter :: c1=(1.0_dp,0.0_dp)
      complex(dp), parameter :: cmplx_0=(0.0_dp,0.0_dp)
      real(dp), parameter :: kB=8.6173324*1e-5_dp
      integer :: nwann,nrpts_ws,nkpts,neps
      integer :: nrx,nry,nrz, nkx,nky,nkz
      real(dp) :: eta, depsil, eps_min, eps_max
      logical :: spin_polarization, dos_only, ham_restore
      character(len=50) :: ctrl_fname="control.inp"
      character(len=50) :: Hrup_fname="wannier90.up_hr.dat"
      character(len=50) :: Hrdn_fname="wannier90.dn_hr.dat"
      character(len=50) :: Hkup_fname="wannier90.up_hk.dat"
      character(len=50) :: Hkdn_fname="wannier90.dn_hk.dat"
      character(len=50) :: Exch_fname="exchange.dat"
      character(len=50) :: Dos_fname="dos"
      character(len=50) :: Delta_fname="delta.dat"
      integer, dimension(:), allocatable :: dos_file
      integer :: ctrl_file, Hrup_file, Hrdn_file, Hkup_file, Hkdn_file, Delta_file, Exch_file
      real(dp), dimension(:,:), allocatable          ::  eval_up, eval_dn ! eigenvalues of the Hamiltonian in k-space
      complex(dp), dimension(:,:,:), allocatable     ::  evec_up, evec_dn ! eigenvectors of the Hamiltonian in k-space
      real(dp) :: jcut   ! output cutoff for the exchange
      integer :: nsets                            ! the number of orbital sets
      integer, dimension(:), allocatable :: s_norb, s_first, s_last ! the number of orbitals in i-set
      integer, dimension(:,:), allocatable :: s_ind ! orbital index of j-orbital in i-set
      !

end module parameters
