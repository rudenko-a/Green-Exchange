module fftw3

contains

subroutine fft3d(arr_in,arr_out,direction)
! Complex 3D Fast Fourier Transform by using the FTTW3 library
! --== FFTW computes an unnormalized transform! ==--
    use, intrinsic :: iso_c_binding
    use parameters, only : dp
implicit none
    include 'fftw3.f03'
complex(dp), dimension(:,:,:), intent(in) :: arr_in
complex(dp), dimension(size(arr_in,1),size(arr_in,2),size(arr_in,3)), intent(out) :: arr_out
integer, intent(in) :: direction ! -1 = forward, +1 = backward
integer :: nx,ny,nz
logical, save :: FirstCall_fw=.true.
logical, save :: FirstCall_bw=.true.
! intrinsic FFTW C variables
type(C_PTR), save :: plan_fw
type(C_PTR), save :: plan_bw
!integer*8 :: plan
!complex(dp), dimension(size(arr_in,1),size(arr_in,2),size(arr_in,3)) :: fft_arr
complex(C_DOUBLE_COMPLEX), dimension(size(arr_in,1),size(arr_in,2),size(arr_in,3)) :: fft_arr
!complex(C_DOUBLE_COMPLEX), dimension(size(arr_in,1),size(arr_in,2),size(arr_in,3)) :: fft_out
!
nx=size(arr_in,1) 
ny=size(arr_in,2) 
nz=size(arr_in,3) 

fft_arr(:,:,:) = arr_in(:,:,:)

if (direction == -1) then
   if (FirstCall_fw) then
      call dfftw_plan_dft_3d(plan_fw,nx,ny,nz,fft_arr,fft_arr,FFTW_FORWARD,FFTW_PATIENT)
      FirstCall_fw=.false.
      fft_arr(:,:,:) = arr_in(:,:,:)
   endif
   call dfftw_execute_dft(plan_fw,fft_arr,fft_arr)
!   call dfftw_destroy_plan(plan)
elseif (direction == 1) then
   if (FirstCall_bw) then
      call dfftw_plan_dft_3d(plan_bw,nx,ny,nz,fft_arr,fft_arr,FFTW_BACKWARD,FFTW_PATIENT)
      FirstCall_bw=.false.
      fft_arr(:,:,:) = arr_in(:,:,:)
   endif
   call dfftw_execute_dft(plan_bw,fft_arr,fft_arr)
!   call dfftw_destroy_plan(plan)
else
   stop 'FFT direction is not defined. Program stops!'
endif

arr_out(:,:,:) = fft_arr(:,:,:)

end subroutine fft3d
!
!
!
subroutine fft_shift(arr_in,arr_out)
! Shifts the origin of the FFT grid to the center
!
use parameters, only : dp, nrx, nry, nrz
implicit none
real(dp), dimension(:,:,:), intent(in) :: arr_in
real(dp), dimension(nrx,nry,nrz), intent(out) :: arr_out
integer :: irx, iry, irz !, loop_rpt
   !
   !
   do irx=nint(real(nrx,dp)/2.0_dp)+1,nrx
      do iry=nint(real(nry,dp)/2.0_dp)+1,nry
         do irz=nint(real(nrz,dp)/2.0_dp)+1,nrz
            !
!            loop_rpt = irz + nrz*(iry-1) + nrz*nry*(irx-1)
            arr_out(irx-nint(real(nrx,dp)/2.0_dp),iry-nint(real(nry,dp)/2.0_dp),irz-nint(real(nrz,dp)/2.0_dp)) = arr_in(irx,iry,irz)
            !
         enddo
         do irz=1,nint(real(nrz,dp)/2.0_dp)
            ! 
!            loop_rpt = irz + nrz*(iry-1) + nrz*nry*(irx-1)
            arr_out(irx-nint(real(nrx,dp)/2.0_dp),iry-nint(real(nry,dp)/2.0_dp),irz+nint(real(nrz,dp)/2.0_dp)) = arr_in(irx,iry,irz)
            !
         enddo
      enddo
      do iry=1,nint(real(nry,dp)/2.0_dp)
         do irz=nint(real(nrz,dp)/2.0_dp)+1,nrz
            !
!            loop_rpt = irz + nrz*(iry-1) + nrz*nry*(irx-1)
            arr_out(irx-nint(real(nrx,dp)/2.0_dp),iry+nint(real(nry,dp)/2.0_dp),irz-nint(real(nrz,dp)/2.0_dp)) = arr_in(irx,iry,irz)
            !
         enddo
         do irz=1,nint(real(nrz)/2.0_dp)
            !
!            loop_rpt = irz + nrz*(iry-1) + nrz*nry*(irx-1)
            arr_out(irx-nint(real(nrx,dp)/2.0_dp),iry+nint(real(nry,dp)/2.0_dp),irz+nint(real(nrz,dp)/2.0_dp)) = arr_in(irx,iry,irz)
            !
         enddo
      enddo
   enddo
   !
   do irx=1,nint(real(nrx,dp)/2.0_dp)
      do iry=nint(real(nry,dp)/2.0_dp)+1,nry
         do irz=nint(real(nrz,dp)/2.0_dp)+1,nrz
            !
!            loop_rpt = irz + nrz*(iry-1) + nrz*nry*(irx-1)
            arr_out(irx+nint(real(nrx,dp)/2.0_dp),iry-nint(real(nry,dp)/2.0_dp),irz-nint(real(nrz,dp)/2.0_dp)) = arr_in(irx,iry,irz)
            !
         enddo
         do irz=1,nint(real(nrz,dp)/2.0_dp)
            !
!            loop_rpt = irz + nrz*(iry-1) + nrz*nry*(irx-1)
            arr_out(irx+nint(real(nrx,dp)/2.0_dp),iry-nint(real(nry,dp)/2.0_dp),irz+nint(real(nrz,dp)/2.0_dp)) = arr_in(irx,iry,irz)
            !
         enddo
      enddo
      do iry=1,nint(real(nry,dp)/2.0_dp)
         do irz=nint(real(nrz,dp)/2.0_dp)+1,nrz
            !
!            loop_rpt = irz + nrz*(iry-1) + nrz*nry*(irx-1)
            arr_out(irx+nint(real(nrx,dp)/2.0_dp),iry+nint(real(nry,dp)/2.0_dp),irz-nint(real(nrz,dp)/2.0_dp)) = arr_in(irx,iry,irz)
            !
         enddo
         do irz=1,nint(real(nrz,dp)/2.0_dp)
            !
!            loop_rpt = irz + nrz*(iry-1) + nrz*nry*(irx-1)
            arr_out(irx+nint(real(nrx,dp)/2.0_dp),iry+nint(real(nry,dp)/2.0_dp),irz+nint(real(nrz,dp)/2.0_dp)) = arr_in(irx,iry,irz)
            !
         enddo
      enddo
   enddo
   !
   !
end subroutine fft_shift

end module fftw3
