!=============================================================
      module mainvar
!=============================================================
      implicit none

!     Constants:
      complex*8, parameter :: jcpx=(0.d0,1.d0) !imaginary unit
      real*8, parameter    :: pi = 4.d0*atan(1.d0)

!     Variables:
      integer :: nstep,ngrid,nwrite,nstate,krange_len
      real*8  :: xmin,xmax
      real*8  :: dt,dx
      real*8  :: E_tot, E_kin, E_pot
      real*8  :: redmass,k_ini,sigma_ini,xcenter_ini
      real*8, allocatable :: krange(:)
      real*8, allocatable :: xgrid(:)
      complex*8, allocatable :: psi(:,:),Op_Kin(:),psi_ad(:,:)
      complex*8, allocatable :: Op_Pot(:,:,:)

      real*8, allocatable :: eigenvals(:,:) !eigenvalues of the Hamiltonian
      real*8, allocatable :: u(:,:) !eigenvectors matrix for diagonalization
      real*8, allocatable :: u_matrices(:,:,:)
      
      contains
     !--------------------------------------------
      subroutine allocate_main_arrays(nalloc)
      implicit none
      integer :: nalloc

      if(nalloc == 1) then
         allocate(xgrid(ngrid))
         allocate(psi_ad(nstate, ngrid))
         allocate(psi(nstate, ngrid))
         allocate(Op_pot(nstate, nstate, ngrid))
         allocate(Op_kin(ngrid)) !always diagonal in diabatic representation
         allocate(u(nstate, nstate))
         allocate(u_matrices(nstate, nstate,ngrid))
         allocate(eigenvals(nstate, ngrid))
      else
         deallocate(xgrid)
         deallocate(psi)
         deallocate(psi_ad)
         deallocate(Op_pot)
         deallocate(Op_kin)
         deallocate(u)
         deallocate(u_matrices)
         deallocate(eigenvals)
       endif 

      end subroutine
     !--------------------------------------------

      endmodule
!=============================================================

