!
!#### This program calculates exchange integrals (spin-polarized case) or magnetic susceptibility (non-polarized case) 
!#### in the Wannier function basis by using the Green's function formalism
!
program WannEx

use parameters
use functions
use fftw3
implicit none

  integer                                        ::  i,j,k,loop_rpt,loop_eps,loop_kpt,m1,m2,dummy_loop,m3,m4,ei,ej,m1_new,m2_new
  integer :: irx,iry,irz, ikx,iky,ikz
  integer, dimension(:), allocatable             ::  ndegen, ndegen_tmp  ! degeneracy of the real-space grid points
  integer, dimension(:,:),allocatable            ::  irvec, irvec_tmp  ! WS grid-points in lattice vectors
  real(dp), dimension(:,:),allocatable           ::  kpt_latt  ! kpoints in lattice vectors
  real(dp)                                       ::  epsil, rdotk  ! energy variable, R*k scalar product, magn. mom
  real(dp)                                       ::  time_begin, time_end, loop_begin, loop_end ! timings
  complex(dp)                                    ::  fac,cdotc,G
  character(len=10) :: ic ! dos file extension
  !
  complex(dp), dimension(:,:), allocatable       ::  Del ! \Delta=V_{\up}-V_{\down}, V_{\sigma}=sum_{k}H_{\sigma}(k)
  real(dp), dimension(:,:,:,:,:), allocatable    ::  J_r, J_buffer, Jshift ! Exchange integral (or susceptibility)
  complex(dp), dimension(:,:,:,:,:), allocatable ::  J_q
  real(dp)                                       ::  J_0
  !
  real(dp),dimension(:,:,:),allocatable          ::  Jav,Jav_grid
  !
  complex(dp), dimension(:,:), allocatable       ::  Theta ! Theta matrix for Tc
  complex(dp), dimension(:,:), allocatable       ::  Theta_evec 
  real(dp), dimension(:), allocatable            ::  Theta_eval
  !
  real(dp), dimension(:,:), allocatable          ::  pdos_up, pdos_dn ! DOS projected onto WFs + total DOS
  real(dp), dimension(:), allocatable            ::  dos_up, dos_dn, moment, buffer_up, buffer_dn
  complex(dp), dimension(:,:,:), allocatable     ::  Kern_J ! temporary array for calculating J_r, J_q
  !
  complex(dp), dimension(:,:,:), allocatable     ::  Hrup, Hrdn, Hkup, Hkdn, H_tmp1, H_tmp2 ! real-space and reciprocal-space Hamiltonians
  complex(dp), dimension(:,:,:,:,:), allocatable ::  Gk_up, Gk_dn      ! Green's functions in k-space
  complex(dp), dimension(:,:,:,:,:), allocatable ::  Gr_up, Gr_dn      ! Green's functions in R-space
  complex(dp), dimension(:,:), allocatable       ::  DelG1, DelG2      ! DelG_{ij} = \sum_k Delta_{ik}*G_{kj} - intermediate matrix products for J_{ij} calc.
  !
  real(dp) :: pdos_set1, pdos_set2, nelec
  !
  logical :: green_spectral=.true.
!
!
!#########################################
! 		INPUT PART
!#########################################
!
call cpu_time(time_begin)
!
call MPI_INIT(ierr)
call MPI_COMM_RANK(MPI_COMM_WORLD, iproc, ierr)
call MPI_COMM_SIZE(MPI_COMM_WORLD, nprocs, ierr)
!
if (iproc==0) write(*,'(a,i5,a)') 'Program is running on ',nprocs,' nodes'
!
!
!
call read_input()
!
if (iproc==0) then
    write(*,*) 'Version 10.10.2013'
    print *, '========================='
    print *, 'List of input parameters:'
    print *, 'nwann: ', nwann
    write(*,'(a,f10.5)')  ' eta: ', eta
    write(*,'(a,f10.5)')  ' depsil: ', depsil
    write(*,'(a,f10.5)')  ' eps_min: ', eps_min
    write(*,'(a,f10.5)')  ' eps_max: ', eps_max
    print *, 'dos_only: ', dos_only
    print *, 'spin_polarization: ', spin_polarization
    print *, 'ham_restore: ', ham_restore
    print *, 'nkx: ', nkx
    print *, 'nky: ', nky
    print *, 'nkz: ', nkz
    write(*,'(a,f10.5)') ' jcut: ', jcut
    print *, 'nsets: ', nsets
    do i=1,nsets
       write(*,*) (s_ind(i,j),j=1,6), s_norb(i)
    enddo
    print *, '========================='
endif
!
!
open(newunit(Hrup_file),file=trim(Hrup_fname),status='old',err=102) ! Open file with a real-space Hamiltonian in the WF basis
read(Hrup_file,*)
read(Hrup_file,*) i                                 
if (i /= nwann) stop 'ERROR: Number of orbitals in the spin-up Hamiltonian differs from nwann, aborting...'
read(Hrup_file,*) nrpts_ws                               ! Number of points in the real-space grid (Wigner-Seitz supercell)
allocate(ndegen(nrpts_ws))
read(Hrup_file,'(15I5)') (ndegen(i),i=1,nrpts_ws)        ! Degeneracy of a point in the grid
allocate(Hrup(nwann,nwann,nrpts_ws))     ;     Hrup(:,:,:)=cmplx_0
allocate(irvec(3,nrpts_ws))
if (iproc==0) write(*,'(a)') "Reading spin-up real-space Hamiltonian... "
do loop_rpt=1,nrpts_ws
   do m1=1,nwann
      do m2=1,nwann
!         read(Hrup_file,*) irvec(:,loop_rpt), j, i, Hrup(j,i,loop_rpt)         ! Read Hamiltonian
         read(Hrup_file,'(5I5,2F12.6)') irvec(:,loop_rpt), j, i, Hrup(j,i,loop_rpt)  ! Read Hamiltonian
      enddo
enddo
enddo
close(Hrup_file)
if (iproc==0) write(*,'(a)') "Done!"
!
!
    if (trim(Hrup_fname)==trim(Hrdn_fname)) then
       stop 'ERROR: Spin-up and spin-down Hamiltonian files must have different names, aborting...'
    endif
    open(newunit(Hrdn_file),file=trim(Hrdn_fname),status='old',err=103) ! Read spin-down Hamiltonian
    read(Hrdn_file,*)
    read(Hrdn_file,*) i
    if (i /= nwann) stop 'ERROR: Number of orbitals in the spin-down Hamiltonian differs from nwann, aborting...'
    read(Hrdn_file,*) j
    if (j /= nrpts_ws) stop 'ERROR: nrpts for spin-up Hamiltonian is not the same as for spin-down Hamiltonian, aborting...'
    allocate(ndegen_tmp(nrpts_ws))
    read(Hrdn_file,'(15I5)') (ndegen_tmp(i),i=1,nrpts_ws)
    do i=1,nrpts_ws
       if (ndegen_tmp(i) /= ndegen(i)) stop 'ERROR: ndegen for spin-up Hamiltonian is not the same as for spin-down Hamiltonian, aborting...'
    enddo
    deallocate(ndegen_tmp)
    allocate(Hrdn(nwann,nwann,nrpts_ws)) ;    Hrdn(:,:,:)=cmplx_0
    allocate(irvec_tmp(3,nrpts_ws))
 if (iproc==0) write(*,'(a)') "Reading spin-down real-space Hamiltonian... "
    do loop_rpt=1,nrpts_ws
      do m1=1,nwann
          do m2=1,nwann
!             read(Hrdn_file,*) irvec_tmp(:,loop_rpt), j, i, Hrdn(j,i,loop_rpt)  ! Read Hamiltonian
             read(Hrdn_file,'(5I5,2F12.6)') irvec_tmp(:,loop_rpt), j, i, Hrdn(j,i,loop_rpt)  ! Read Hamiltonian
          enddo
       enddo
       if ((irvec_tmp(1,loop_rpt) /= irvec(1,loop_rpt)) .or. (irvec_tmp(2,loop_rpt) /= irvec(2,loop_rpt)) .or. (irvec_tmp(3,loop_rpt) /= irvec(3,loop_rpt))) &
           stop 'Error: irvec for spin-up Hamiltonian is not the same as for spin-down Hamiltonian, aborting...'
    enddo
    deallocate(irvec_tmp)
    close(Hrdn_file)
 if (iproc==0) write(*,'(a)') "Done!"
!
!
allocate(kpt_latt(3,nkpts))
do ikx=1,nkx
   do iky=1,nky
      do ikz=1,nkz
         !
         loop_kpt = ikz + nkz*(iky-1) + nkz*nky*(ikx-1)
         kpt_latt(1,loop_kpt) = real(ikx-1,dp)/real(nkx)
         kpt_latt(2,loop_kpt) = real(iky-1,dp)/real(nky)
         kpt_latt(3,loop_kpt) = real(ikz-1,dp)/real(nkz)
         !
      enddo
   enddo
enddo
!
! Sorting orbitals in the Hamiltonians accorting to the defined sets
allocate(H_tmp1(nwann,nwann,nrpts_ws),H_tmp2(nwann,nwann,nrpts_ws))
H_tmp1=cmplx_0 ; H_tmp2=cmplx_0
m1_new=0 ; m2_new=0
do i=1,nsets
   do ei=1,s_norb(i)
      m1_new = m1_new + 1
      do m2=1,nwann
         H_tmp1(m1_new,m2,:) = Hrup(s_ind(i,ei),m2,:)
         H_tmp2(m1_new,m2,:) = Hrdn(s_ind(i,ei),m2,:)
      enddo
   enddo
enddo
!
do j=1,nsets
   do ej=1,s_norb(j)
      m2_new = m2_new + 1
      do m1=1,nwann
         Hrup(m1,m2_new,:) = H_tmp1(m1,s_ind(j,ej),:)
         Hrdn(m1,m2_new,:) = H_tmp2(m1,s_ind(j,ej),:)
      enddo
   enddo
enddo
!
deallocate(H_tmp1,H_tmp2)
!
!						 Write sorted Hamiltonians in stdout
!do loop_rpt=1,nrpts_ws
!       do m1=1,nwann
!          do m2=1,nwann
!             write(*,'(5i5,2f12.6)') irvec(:,loop_rpt), m2, m1, Hrup(m2,m1,loop_rpt)
!          enddo
!       enddo
!enddo
!do loop_rpt=1,nrpts_ws
!       do m1=1,nwann
!          do m2=1,nwann
!             write(*,'(5i5,2f12.6)') irvec(:,loop_rpt), m2, m1, Hrdn(m2,m1,loop_rpt)
!          enddo
!       enddo
!enddo
!
!
!#########################################
! 		MAIN PART
!#########################################
!
!
!#### Fourier transform Hij(R) -> Hij(k)
allocate(Hkup(nwann,nwann,nkpts),Hkdn(nwann,nwann,nkpts))
Hkup(:,:,:)=cmplx_0 ; Hkdn(:,:,:)=cmplx_0
!
!
if (ham_restore .EQV. .true.) then
    !
    if (iproc==0) write(*,*) "Restoring Hij(k) from file..."
    !
    open(newunit(Hkup_file),file=trim(Hkup_fname),status='old',form='unformatted',err=122)
    read(Hkup_file) dummy_loop
    read(Hkup_file) i
    read(Hkup_file) j
    read(Hkup_file) k
          if ((dummy_loop /= loop_kpt).or.(i /= nkx).or.(j /= nky).or.(k /= nkz)) then
              stop 'ERROR: input k-points grid does not match with the grid found in the spin-up Hamiltonian, aborting...'
          endif
       open(newunit(Hkdn_file),file=trim(Hkdn_fname),status='old',form='unformatted',err=133)
       read(Hkdn_file) dummy_loop
       read(Hkdn_file) i
       read(Hkdn_file) j
       read(Hkdn_file) k
          if ((dummy_loop /= loop_kpt).or.(i /= nkx).or.(j /= nky).or.(k /= nkz)) then
              stop 'ERROR: input k-points grid does not match with the grid found in the spin-down Hamiltonian, aborting...'
          endif
    !
    read(Hkup_file) Hkup(:,:,:)
    read(Hkdn_file) Hkdn(:,:,:)
    !       
else
    if (iproc==0) then
       !
       write(*,'(a)') "Slow fourier transform Hij(R) -> Hij(k) for both spins, Hij(k) is being stored... "
       !
       open(newunit(Hkup_file),file=trim(Hkup_fname),status='replace',form='unformatted',err=122)
       write(Hkup_file) loop_kpt
       write(Hkup_file) nkx
       write(Hkup_file) nky
       write(Hkup_file) nkz
          open(newunit(Hkdn_file),file=trim(Hkdn_fname),status='replace',form='unformatted',err=133)
          write(Hkdn_file) loop_kpt
          write(Hkdn_file) nkx
          write(Hkdn_file) nky
          write(Hkdn_file) nkz
       !
    endif
       !
    do loop_kpt=1,nkpts
       !
       do loop_rpt=1,nrpts_ws                        ! slow fourier transform on ONE processor !!! (to be parallelized)
          rdotk = twopi*dot_product(kpt_latt(:,loop_kpt),real(irvec(:,loop_rpt),dp))
          fac = exp(ci*rdotk)
          Hkup(:,:,loop_kpt) = Hkup(:,:,loop_kpt) + fac*Hrup(:,:,loop_rpt)/real(ndegen(loop_rpt),dp)
             Hkdn(:,:,loop_kpt) = Hkdn(:,:,loop_kpt) + fac*Hrdn(:,:,loop_rpt)/real(ndegen(loop_rpt),dp)
       enddo
       !
!   write(*,*) Hkup(:,:,loop_kpt)
    enddo
    !
    if (iproc==0) write(Hkup_file) Hkup(:,:,:)    ! store Hij(k) in a file
    if (iproc==0) write(Hkdn_file) Hkdn(:,:,:)   
    if (iproc==0) close(Hkup_file)
    if (iproc==0) close(Hkdn_file)
    if (iproc==0) write(*,'(a)') "Done!"
endif
!
deallocate(Hrup,Hrdn)
!
!
!#### Delta calculation: \Delta=V_{\up}-V_{\down}, where V_{\sigma}=sum_{k}H_{\sigma}(k)
allocate(Del(nwann,nwann))
Del(:,:)=cmplx_0
!
do loop_kpt=1,nkpts
   Del(:,:) = Del(:,:) + Hkup(:,:,loop_kpt) - Hkdn(:,:,loop_kpt)
enddo   
 !
Del(:,:) = Del(:,:)/real(nkpts,dp)
!
if (iproc==0) then
   !
   open(newunit(Delta_file),file=trim(Delta_fname),status='replace',err=108)
   do m1=1,nwann
      do m2=1,nwann
         write(Delta_file,'(2I5,2F12.6)') m2, m1, Del(m2,m1)
      enddo
   enddo
   !
   close(Delta_file)
endif
!
!
!#### Calculation of eigenvalues and eigenvectors of H_k
if (iproc==0) write(*,*) "Hamiltonian diagonalization (both spins)..."
!
if (green_spectral .EQV. .true.) then   ! calculate eigvec and eigval if the GF is calculated in spectral representation
   allocate(eval_up(nwann,nkpts),eval_dn(nwann,nkpts))
   allocate(evec_up(nwann,nwann,nkpts),evec_dn(nwann,nwann,nkpts))
   do loop_kpt=1,nkpts
      call eigen(Hkup(:,:,loop_kpt),eval_up(:,loop_kpt),evec_up(:,:,loop_kpt))  ! Store eigval and eigvec of Hij(k) for all k-points
      call eigen(Hkdn(:,:,loop_kpt),eval_dn(:,loop_kpt),evec_dn(:,:,loop_kpt))
   enddo
   deallocate(Hkup,Hkdn)
endif
!
!
epsil = eps_min - depsil
neps = nint((eps_max - eps_min)/depsil) + 1       ! number of points along the energy axis
!
!
allocate(Gk_up(nwann,nwann,nkx,nky,nkz),Gk_dn(nwann,nwann,nkx,nky,nkz))
allocate(Gr_up(nwann,nwann,nrx,nry,nrz),Gr_dn(nwann,nwann,nrx,nry,nrz))
allocate(pdos_up(nwann,neps),pdos_dn(nwann,neps))
allocate(dos_up(neps),dos_dn(neps))
allocate(DelG1(nwann,nwann),DelG2(nwann,nwann))
allocate(Kern_J(nwann,nwann,nrx*nry*nrz))
allocate(J_r(nsets,nsets,nrx,nry,nrz))
!
!
pdos_up(:,:)=cmplx_0  ;  pdos_dn(:,:)=cmplx_0
dos_up(:)=cmplx_0  ;  dos_dn(:)=cmplx_0
J_r(:,:,:,:,:)=0.0_dp
!
!
call MPI_Barrier(MPI_COMM_WORLD,ierr)
!
!
call cpu_time(loop_begin)
if (iproc==0) write(*,'(a,f10.3,a)') 'Time for initialization: ', loop_begin-time_begin, ' seconds'
!
! #### Start of a big loop over epsilon
if (iproc==0) write(*,*) "Starting of a big loop over epsilon..."
!
do loop_eps=floor((iproc*neps)/nprocs+1.d0),floor(((iproc+1.d0)*neps/nprocs))
!do loop_eps = 1,neps               
   !
   epsil = eps_min + depsil*(loop_eps-1)
   Gr_up(:,:,:,:,:) = cmplx_0  ;  Gr_dn(:,:,:,:,:) = cmplx_0
   Gk_up(:,:,:,:,:) = cmplx_0  ;  Gk_dn(:,:,:,:,:) = cmplx_0
   Kern_J(:,:,:) = cmplx_0
   !
   !#### Green's function Gmn(k) calculation
   write(*,'(a,f8.4,a,i5)') "Calculation of the Green's function for eps=",epsil,"..., rank=",iproc
   do loop_kpt = 1,nkpts                   ! loop over k
      !
      ikx = nint(nkx*kpt_latt(1,loop_kpt)) + 1
      iky = nint(nky*kpt_latt(2,loop_kpt)) + 1
      ikz = nint(nkz*kpt_latt(3,loop_kpt)) + 1
      !
      !
      if (green_spectral .EQV. .true.) then
         call green_spec(eta,epsil,Gk_up(:,:,ikx,iky,ikz),loop_kpt,.true.)   ! .true. = spin up
         call green_spec(eta,epsil,Gk_dn(:,:,ikx,iky,ikz),loop_kpt,.false.)  ! .false. = spin down
      else
         call green_inv(eta,epsil,Hkup(:,:,loop_kpt),Gk_up(:,:,ikx,iky,ikz))
         call green_inv(eta,epsil,Hkdn(:,:,loop_kpt),Gk_dn(:,:,ikx,iky,ikz))
      endif
      !
      !
      !#### DOS calculation for each epsilon - DOS(E) = -Im(\sum_k[Tr{Gmn(k,E)}]/(Nk*pi), similar for pdos.
      do m1=1,nwann   
         !
         pdos_up(m1,loop_eps) = pdos_up(m1,loop_eps) + aImag(Gk_up(m1,m1,ikx,iky,ikz))
         pdos_dn(m1,loop_eps) = pdos_dn(m1,loop_eps) + aImag(Gk_dn(m1,m1,ikx,iky,ikz))
         dos_up(loop_eps) = dos_up(loop_eps) + aImag(Gk_up(m1,m1,ikx,iky,ikz))
         dos_dn(loop_eps) = dos_dn(loop_eps) + aImag(Gk_dn(m1,m1,ikx,iky,ikz))
         !
      enddo
   !
   enddo                                   ! end of the loop over k-points
   !
   write(*,'(a,f8.4,a,i5)') "FINISHED Calculation of the Green's function for eps=",epsil,"..., rank=",iproc
   !
   !#### DOS normalization
   dos_up(loop_eps)=-dos_up(loop_eps)/(pi*real(nkpts,dp)) ; dos_dn(loop_eps)=-dos_dn(loop_eps)/(pi*real(nkpts,dp))               !Total DOS (normalized)
   do m1=1,nwann
      pdos_up(m1,loop_eps)=-pdos_up(m1,loop_eps)/(pi*real(nkpts,dp)) ; pdos_dn(m1,loop_eps)=-pdos_dn(m1,loop_eps)/(pi*real(nkpts,dp))  !DOS projected on each WF
   enddo
   !
   if (dos_only .EQV. .true.) cycle
   !
   !
   !#### Inverse Fourier transform Gij(k) -> Gij(R)
   if (iproc==0) write(*,'(a,f8.4,a)') "Inverse Fourier transform Gij(k) -> Gij(R) (both spins) for eps=",epsil,"..."
   do m1=1,nwann
   do m2=1,nwann
      call fft3d(Gk_up(m1,m2,:,:,:),Gr_up(m1,m2,:,:,:),-1)  ! -1 ------- K->R; +1 ------- R->K
      Gr_up(m1,m2,:,:,:) = Gr_up(m1,m2,:,:,:)/real(nkpts,dp)
      call fft3d(Gk_dn(m1,m2,:,:,:),Gr_dn(m1,m2,:,:,:),-1)
      Gr_dn(m1,m2,:,:,:) = Gr_dn(m1,m2,:,:,:)/real(nkpts,dp)
   enddo
   enddo
!
!
!#### J_r calculation
   write(*,'(a,f8.4,a)') "Calculation of the exchange contribution for eps=",epsil,"..."
   do irx=1,nrx
      do iry=1,nry
         do irz=1,nrz
         !
         DelG1(:,:)=cmplx_0
         DelG2(:,:)=cmplx_0
         !
         loop_rpt = irz + nrz*(iry-1) + nrz*nry*(irx-1)
         !
         !
         do i=1,nsets   ! Multiplication of block diagonal matrix Del_{ij} and general matrix G_{ij}. Each block corresponds to a predefined set of orbitals.
         do j=1,nsets
            !
            DelG1(s_first(i):s_last(i),s_first(j):s_last(j)) = matmul(Del(s_first(i):s_last(i),s_first(i):s_last(i)),Gr_up(s_first(i):s_last(i),s_first(j):s_last(j),irx,iry,irz))
            DelG2(s_first(i):s_last(i),s_first(j):s_last(j)) = matmul(Del(s_first(i):s_last(i),s_first(i):s_last(i)),Gr_dn(s_first(i):s_last(i),s_first(j):s_last(j),irx,iry,irz))
            !
            !                                                         
         enddo
         enddo
         !
         Kern_J(:,:,loop_rpt) = DelG1(:,:)*transpose(DelG2(:,:))
         !
         do i=1,nsets
         do j=1,nsets
            do ei=0,s_norb(i)-1
            do ej=0,s_norb(j)-1
               J_r(i,j,irx,iry,irz) = J_r(i,j,irx,iry,irz) + aImag(Kern_J(s_first(i)+ei,s_first(j)+ej,loop_rpt))
            enddo
            enddo
         enddo
         enddo

               !
         !
         enddo
      enddo
   enddo
   !
!   do m1=1,nwann
!      do m2=1,nwann
!         J_r(m1,m2,:) = J_r(m1,m2,:) + Imag(Kern_J(m1,m2,:))  ! sum over epsilon
!      enddo
!   enddo
   !
   !
enddo                                      ! end of the loop over epsilon
!
call MPI_Barrier(MPI_COMM_WORLD,ierr)
call cpu_time(loop_end)
if (iproc==0) write(*,'(a,f10.3,a)') 'Time for the loop over epsilon: ', loop_end-loop_begin, ' seconds'
!
!
if (iproc==0) write(*,*) 'Printing out DOS...'
!
!
!#### Opening DOS files
if (iproc==0) then
   allocate(dos_file(nsets+1))
   open(newunit(dos_file(nsets+1)),file=trim(Dos_fname)//'.dat',status='replace',err=107) ! total DOS
   write(dos_file(nsets+1),'(a)') '#  energy  dos_up      dos_dn        dos_up+dos_dn'
   do i=1,nsets
      write(ic,'(i4)') i
      ic=adjustl(ic)
      open(newunit(dos_file(i)),file=trim(Dos_fname)//trim(ic)//'.dat',status='replace',err=107) ! projected DOSs
      write(dos_file(i),'(a)') '#  energy  dos_up      dos_dn      dos_up-dos_dn'
   enddo
endif
!
!
!#### MPI_Reduce for the calculated DOS
allocate(buffer_up(neps),buffer_dn(neps))
buffer_up(:)=0.0_dp  ;  buffer_dn(:)=0.0_dp
Call MPI_Reduce(dos_up(:), buffer_up(:), neps, MPI_DOUBLE_PRECISION, MPI_SUM, 0, MPI_COMM_WORLD, ierr)
Call MPI_Reduce(dos_dn(:), buffer_dn(:), neps, MPI_DOUBLE_PRECISION, MPI_SUM, 0, MPI_COMM_WORLD, ierr)
dos_up(:) = buffer_up(:) ; dos_dn(:) = buffer_dn(:)
nelec=0.0_dp
if (iproc==0) then
   do loop_eps=1,neps
      epsil = eps_min + depsil*(loop_eps-1)
      if ((loop_eps/=1).or.(loop_eps/=neps)) then
         nelec = nelec + (dos_up(loop_eps) + dos_dn(loop_eps))*(eps_max-eps_min)/real(neps-1,dp)
      else
         nelec = nelec + (dos_up(loop_eps) + dos_dn(loop_eps))*(eps_max-eps_min)/(real(neps-1,dp)*2.0_dp)  ! trapezoidal integration
      endif
      write(dos_file(nsets+1),'(f9.4,2f12.8,1f14.8)') epsil, dos_up(loop_eps), dos_dn(loop_eps), dos_up(loop_eps)+dos_dn(loop_eps)
   enddo
   write(*,*)
   write(*,*) 'Number of electrons: ', nelec
endif
!
!
do m1=1,nwann
   buffer_up(:)=0.0_dp ; buffer_dn(:)=0.0_dp
   Call MPI_Reduce(pdos_up(m1,:), buffer_up(:), neps, MPI_DOUBLE_PRECISION, MPI_SUM, 0, MPI_COMM_WORLD, ierr)
   Call MPI_Reduce(pdos_dn(m1,:), buffer_dn(:), neps, MPI_DOUBLE_PRECISION, MPI_SUM, 0, MPI_COMM_WORLD, ierr)
   pdos_up(m1,:) = buffer_up(:) ; pdos_dn(m1,:) = buffer_dn(:)
enddo
!
allocate(moment(nsets))
moment(:)=0.0_dp
if (iproc==0) then
   do loop_eps=1,neps
     j = 0 
     epsil = eps_min + depsil*(loop_eps-1)
      do i=1,nsets
         pdos_set1=0.0_dp ; pdos_set2=0.0_dp
         do ei=1,s_norb(i)
            j = j + 1
            if ((loop_eps/=1).or.(loop_eps/=neps)) then
               moment(i) = moment(i) + (pdos_up(j,loop_eps) - pdos_dn(j,loop_eps))*(eps_max-eps_min)/real(neps-1,dp)
            else
               moment(i) = moment(i) + (pdos_up(j,loop_eps) - pdos_dn(j,loop_eps))*(eps_max-eps_min)/(real(neps-1,dp)*2.0_dp)
            endif        
            pdos_set1 = pdos_set1 + pdos_up(j,loop_eps)
            pdos_set2 = pdos_set2 + pdos_dn(j,loop_eps)
         enddo
         write(dos_file(i),'(f9.4,3f12.8)') epsil, pdos_set1, pdos_set2, pdos_set1-pdos_set2  ! pdos for the ith-set
      enddo
   enddo
   !
   write(*,*)
   write(*,*) 'Magnetic moment for each set of orbitals: '
   do i=1,nsets
      write(*,'(i5,f12.8)') i, moment(i)
   enddo
endif
!
!
call MPI_Barrier(MPI_COMM_WORLD,ierr)
!
deallocate(dos_up,dos_dn)
deallocate(buffer_up,buffer_dn)
deallocate(pdos_up,pdos_dn)
!
if (iproc==0) then
   close(dos_file(nsets+1))
   do i=1,nsets
      close(dos_file(i))
   enddo
   deallocate(dos_file)
endif
!
!
!#### Normalization of the calculated exchange (susceptibility)
J_r(:,:,:,:,:) = (eps_max - eps_min)*(J_r(:,:,:,:,:))/real(neps,dp) 
J_r(:,:,:,:,:) = J_r(:,:,:,:,:)/fourpi !twopi                  ! 1/(2*pi) because we do not sum over orbital indices (in such a case it would be 1/(4*pi))
!
! diagonal elements in the original unit cell [R=(0,0,0)] are meaningless => void them
do m1=1,nsets
   J_r(m1,m1,1,1,1)=0.0_dp
enddo
!
!
!#### MPI_Reduce for the calculated exchange
if (dos_only .NEQV. .true.) then
   if (iproc==0) write(*,*) 'Collecting parts of the calculated exchange from different processors...'
   allocate(J_buffer(nsets,nsets,nrx,nry,nrz))
    Call MPI_Reduce(J_r(:,:,:,:,:), J_buffer, nsets*nsets*nrx*nry*nrz, MPI_DOUBLE_PRECISION, MPI_SUM, 0, MPI_COMM_WORLD, ierr)
   if (iproc==0) J_r(:,:,:,:,:)=J_buffer(:,:,:,:,:)
   deallocate(J_buffer)
endif
!
!
!#### Sign correction to the exchange
if (iproc==0) then
   do m1=1,nsets
   do m2=1,nsets
      if (moment(m1)*moment(m2) < 0.0_dp) then   
         J_r(m1,m2,:,:,:) = -1.0_dp*J_r(m1,m2,:,:,:)
      endif
   enddo
   enddo
endif
!
!
if (green_spectral .EQV. .true.) then
   deallocate(eval_up,eval_dn)
   deallocate(evec_up,evec_dn)
else
   deallocate(Hkup)
   deallocate(Hkdn)
endif
!
deallocate(Kern_J)
deallocate(Gk_up,Gk_dn)
deallocate(Gr_up,Gr_dn)
deallocate(DelG1,DelG2)
deallocate(Del)
!
if (dos_only .EQV. .true.) then
    deallocate(J_r)
    deallocate(ndegen)
    deallocate(irvec)
    deallocate(kpt_latt)
    deallocate(s_norb,s_first,s_last)
    deallocate(s_ind)
    deallocate(moment)
    if (iproc==0) write(*,*) 'dos_only calculation finished!'
    if (iproc==0) call cpu_time(time_end)
    if (iproc==0) write(*,'(a,f10.3,a)') 'Time elapsed: ', time_end-time_begin, ' seconds'
    stop
endif
!
!#### J_q calculation (Fourier transform J(R) -> J(q))
allocate(J_q(nsets,nsets,nkx,nky,nkz))
J_q(:,:,:,:,:)=cmplx_0
do m1=1,nsets
do m2=1,nsets
   call fft3d(c1*J_r(m1,m2,:,:,:),J_q(m1,m2,:,:,:),1)
enddo
enddo
!
!#########################################
! 		OUTPUT PART
!#########################################
!
!
if (iproc==0) then
!
!#### Write J_r and J_q in the output file
open(newunit(Exch_file),file=trim(Exch_fname), status='replace',err=109)
write(Exch_file,*) 'Exchange in real space'
write(Exch_file,*) 'Number of interaction centers per cell: ', nsets
write(Exch_file,*) 'Total number of cells: ', nrx*nry*nrz
!
allocate(Jshift(nsets,nsets,nrx,nry,nrz))
!
do m1=1,nsets
do m2=1,nsets
   call fft_shift(J_r(m1,m2,:,:,:),Jshift(m1,m2,:,:,:))
enddo
enddo
!
do irx=1,nrx
do iry=1,nry
do irz=1,nrz
   do m1=1,nsets
   do m2=1,nsets
!      if (abs(Jshift(m2,m1,irx,iry,irz)) >= jcut) then
         write(Exch_file,'(5I5,1F12.8)') irx-nint(real(nrx,dp)/2.0_dp)-1,iry-nint(real(nry,dp)/2.0_dp)-1,irz-nint(real(nrz,dp)/2.0_dp)-1, m2, m1, Jshift(m2,m1,irx,iry,irz)
!      endif
   enddo
   enddo
enddo
enddo
enddo
!
deallocate(Jshift)
close(Exch_file)
!
!
!#### Write J_q(0) in the standard output
allocate(Theta(nsets,nsets),Theta_evec(nsets,nsets))
allocate(Theta_eval(nsets))
write(*,*) "Exchange at the Gamma point"
do ikx=1,nkx
do iky=1,nky
do ikz=1,nkz
   loop_kpt = ikz + nkz*(iky-1) + nkz*nky*(ikx-1)
   do m1=1,nsets
   do m2=1,nsets
!      if (abs(kpt_latt(1,loop_kpt)) < 1e-10_dp .and. abs(kpt_latt(2,loop_kpt)) < 1e-10_dp .and. abs(kpt_latt(3,loop_kpt)) < 1e-10_dp) then
         write(*,'(3F10.6,2I5,2F12.8)') kpt_latt(:,loop_kpt), m2, m1, J_q(m2,m1,ikx,iky,ikz)
         Theta(m2,m1) = c1 * J_q(m2,m1,ikx,iky,ikz) ! *(2/3)/kB
!      endif
   enddo
   enddo
enddo
enddo
enddo
!
!
!#### Write eigenvalues of Theta the standard output (that is mean-field Tc)
write(*,*) ""
call eigen(Theta(:,:),Theta_eval(:),Theta_evec(:,:))
write(*,*) "Positive Theta eigenvalues (in eV and in K)"
   do m1=1,nsets
      if (Theta_eval(m1) >= 0) then
         write(*,'(2F12.6)') (2.0_dp/3.0_dp)*Theta_eval(m1),(2.0_dp/3.0_dp)*Theta_eval(m1)/kB
      endif
   enddo
write(*,*) ""
deallocate(Theta,Theta_evec)
deallocate(Theta_eval)
!
!
!#### Write exchange averaged over unit cell in the standard output
write(*,*) "Exchange averaged over unit cell"  ! J_{R}=\sum_{nm}J_{R}^{nm}, n/=m for R=0
allocate(Jav(nrx,nry,nrz),Jav_grid(nrx,nry,nrz))
Jav(:,:,:)=0.0_dp   ;  Jav_grid(:,:,:)=0.0_dp   ;   J_0=0.0_dp
!
do irx=1,nrx
do iry=1,nry
do irz=1,nrz
   loop_rpt = irz + nrz*(iry-1) + nrz*nry*(irx-1)
   if ((irx == 1).and.(iry == 1).and.(irz == 1)) then
      Jav(1,1,1)=0.0_dp
   else
      do m1=1,nsets
      do m2=1,nsets
         Jav(irx,iry,irz) = Jav(irx,iry,irz) + J_r(m2,m1,irx,iry,irz)
      enddo
      enddo
   endif
J_0 = J_0 + Jav(irx,iry,irz)
enddo
enddo
enddo
write(*,*) ""
write(*,*) "Theta_zero (Tc from the averaged exchange - J_0)"
write(*,'(2F12.6)') (2.0_dp/3.0_dp)*J_0, (2.0_dp/3.0_dp)*J_0/kB 
write(*,*) ""
!
call fft_shift(Jav(:,:,:),Jav_grid(:,:,:))
!
do irx=1,nrx
do iry=1,nry
do irz=1,nrz
   if (abs(Jav_grid(irx,iry,irz)) >= jcut) then
      write(*,'(3I5,1F12.8)') irx-nint(real(nrx,dp)/2.0_dp)-1,iry-nint(real(nry,dp)/2.0_dp)-1,irz-nint(real(nrz,dp)/2.0_dp)-1, Jav_grid(irx,iry,irz) 
   endif
enddo
enddo
enddo
!
deallocate(Jav,Jav_grid)
endif !(iproc==0)
!
deallocate(J_r)
deallocate(J_q)
deallocate(ndegen)
deallocate(irvec)
deallocate(kpt_latt)
deallocate(s_norb,s_first,s_last)
deallocate(s_ind)
deallocate(moment)
!
!
if (iproc==0) call cpu_time(time_end)
!
if (iproc==0) write(*,'(a,f10.3,a)') 'Time elapsed: ', time_end-time_begin, ' seconds'
!
call MPI_FINALIZE(ierr)
!
stop
102 print *,'Problem opening R-Hamiltonian file for spin-up',trim(Hrup_fname),', aborting...'
stop
103 print *,'Problem opening R-Hamiltonian file for spin-down',trim(Hrdn_fname),', aborting...'
stop
122 print *,'Problem opening k-Hamiltonian file for spin-up',trim(Hkup_fname),', aborting...'
stop
133 print *,'Problem opening k-Hamiltonian file for spin-down',trim(Hkdn_fname),', aborting...'
stop
107 print *,'Problem opening DOS files, aborting...'
stop
108 print *,'Problem opening Delta file ',trim(Delta_fname),', aborting...'
stop
109 print *,'Cannot create output file for the exchange ',trim(Exch_fname),', aborting...'
end program WannEx

subroutine read_input()
use parameters, only : dp, ctrl_file
use parameters, only : nkpts, nkx, nky, nkz, nrx, nry, nrz, nwann
use parameters, only : eta, depsil, eps_min, eps_max
use parameters, only : dos_only, spin_polarization, ham_restore
use parameters, only : nsets, s_ind, s_norb, s_first, s_last, jcut
use functions, only : newunit, ctrl_fname
implicit none
character(len=255) :: buffer, label
character(len=1024) :: string
integer :: i,ei,j,s_io,lenstr,s_pos,s_cnt
integer :: pos
integer :: ios = 0
integer :: line = 0
logical :: found
!
open(newunit(ctrl_file),file=trim(ctrl_fname),status='old',err=113)
do while (ios == 0)
   !
   read(ctrl_file,'(A)', iostat=ios) buffer
   !
   if (ios == 0) then
      line = line + 1
      pos = scan(buffer, "=")
      label = trim(adjustl(buffer(1:pos-1)))
      buffer = buffer(pos+1:)
      !
      select case (label)
      case ('nwann')
           read(buffer, *, iostat=ios) nwann
      case ('eta')
           read(buffer, *, iostat=ios) eta
      case ('depsil')
           read(buffer, *, iostat=ios) depsil
      case ('eps_min')
           read(buffer, *, iostat=ios) eps_min
      case ('eps_max')
           read(buffer, *, iostat=ios) eps_max
      case ('dos_only')
           read(buffer, *, iostat=ios) dos_only
      case ('spin_polarization')
           read(buffer, *, iostat=ios) spin_polarization
      case ('ham_restore')
           read(buffer, *, iostat=ios) ham_restore
      case ('nkx')
           read(buffer, *, iostat=ios) nkx
      case ('nky')
           read(buffer, *, iostat=ios) nky
      case ('nkz')
           read(buffer, *, iostat=ios) nkz
      case ('jcut')
           read(buffer, *, iostat=ios) jcut
           jcut=abs(jcut)
      case ('nsets')
           read(buffer, *, iostat=ios) nsets
           !
           if (nsets > 0) then
              allocate(s_ind(nsets,255))         ! 255 - maximum number of orbitals for each set
              allocate(s_norb(nsets),s_first(nsets),s_last(nsets))
              s_ind(:,:)=0
              s_norb(:)=0
              do i=1,nsets                       ! read a sequence of orbital indices for each set (undefined length)
                 read(ctrl_file,'(a)',iostat=s_io) string
                 if (s_io>0) stop "Error in reading input file, aborting..."
                 if (s_io<0) stop "Error in reading input file: nsets is wrong, aborting..."
                 j=0
                 s_pos=-1
                 string=adjustl(string)
                 lenstr=len_trim(adjustl(string))
                 if (lenstr==0) then
                    write(*,*) 'ERROR in reading input file: parameter nsets is wrong, aborting...'
                    stop
                 endif
                 do while (s_pos/=0)
                    j=j+1
                    s_pos=scan(trim(adjustl(string)), " ")
                    lenstr=len_trim(adjustl(string))
                    if (s_pos==0) then
                       s_norb(i)=j    ! the number of orbitals in ith-set
                       read(string(1:lenstr),*,iostat=s_io) s_ind(i,j)
                       if (s_io>0) stop "ERROR in reading input file, aborting..."
                       exit
                    endif
                    read(string(1:s_pos-1),*,iostat=s_io) s_ind(i,j)
                    if (s_io>0) stop "ERROR in reading input file, aborting..."
                    string=trim(string(s_pos+1:))
                 enddo
              enddo
              !
              ! check whether the number of orbitals in all sets equal to the total number of orbitals
              s_cnt=0
              do i=1,nsets
                 s_cnt = s_cnt + s_norb(i)
              enddo
              if (s_cnt/=nwann) stop "ERROR: the number of orbitals in sets do not match the total number of orbitals, aborting..."
              !
              ! check whether every single orbital in present in sets
              do j=1,nwann
                 found=.false.
                 do i=1,nsets
                 do ei=1,s_norb(i)
                    if (j==s_ind(i,ei)) then
                       found=.true.
                       exit
                    endif
                 enddo
                 if (found) exit
                 enddo
                 if (.not.found) stop "ERROR: not all orbitals are found in sets, aborting..."
              enddo
              !
              ! generate two auxiliary arrays pointing to the first and last orbital of the ith-set
              j=0
              do i=1,nsets
                 s_first(i) = j + 1
                 do ei=1,s_norb(i)
                    j= j + 1
                 enddo
                 s_last(i) = s_first(i) + s_norb(i) - 1
              enddo
              !
           endif
           !
      end select
      !
      !
   endif
   !
enddo
close(ctrl_file)
!
if ((nkx <= 0).or.(nky <= 0).or.(nkz <= 0)) stop 'ERROR: FFT grid is not defined, aborting...'
!
nrx=nkx ; nry=nky ; nrz=nkz
nkpts = nkx*nky*nkz
!
if (nwann<=0) stop 'ERROR: nwann not defined, aborting...'
if (eta<=0) stop 'ERROR: eta not defined, aborting...'
if (depsil<=0) stop 'ERROR: depsil not defined, aborting...'
if (eps_min>=eps_max) stop 'ERROR: eps_min >= eps_max, aborting...'
if (jcut==0) stop 'ERROR: jcut not defined or equals zero, aborting...'
!
goto 999
113 print *,'ERROR in opening control file ',trim(ctrl_fname),', aborting...'
stop
999 continue

end subroutine read_input
