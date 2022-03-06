!=============================================================
      module mainmethods
!=============================================================
      implicit none

      contains
!======================================================================================
      subroutine Evaluate_Potential_Energy
!======================================================================================
!     Evaluate the kinetic energy of the wavepacket  
!              
!     The evaluation is performed in the Fourier space, as described in 
!     http://www.sam.math.ethz.ch/~hiptmair/StudentProjects/Bourquin/thesis.pdf
!======================================================================================
      use mainvar
      implicit none
      integer :: i,j,k
      real*8, allocatable :: VV(:,:)

      allocate(VV(nstate,nstate))

      E_pot = 0.d0
        
      do i =1,ngrid
        call Vpot_TullyA(xgrid(i), VV)
        write(*,*) i
        ! do j=1,2
        !     write(*,*) (VV(j,k),k=1,2)
        ! enddo
        E_pot = E_pot + dot_product(conjg(psi(:,i)), matmul(VV, psi(:,i)))*dx
      enddo

      deallocate(VV)

      end subroutine
!======================================================================================

!======================================================================================
      subroutine Evaluate_Kinetic_Energy
!======================================================================================
!     Evaluate the kinetic energy of the wavepacket  
!              
!     The evaluation is performed in the Fourier space, as described in 
!     http://www.sam.math.ethz.ch/~hiptmair/StudentProjects/Bourquin/thesis.pdf
!======================================================================================
      use mainvar
      implicit none
      integer :: i,j
      real*8 :: kx,kx_old,final_norm

      final_norm = 0.d0
      E_kin = 0.d0
      do i = 1,nstate
        call fourier(0,ngrid,psi(i,:)) !initialize
        call fourier(1,ngrid,psi(i,:)) !do forwad FFT
        call fourier(2,ngrid,psi(i,:)) !deallocate
      enddo
     ! Multiply the transformed the wavefunction by the 
     ! reciprocal vector \omega
      do i = 1,nstate
        j = 1
        kx_old=2.d0*pi*float(j)/(float(ngrid)*dx)
        
        do j=2,ngrid
           if(j < ngrid/2) then
             kx=2.d0*pi*float(j)/(float(ngrid)*dx)
           else
             kx=2.d0*pi*(float(j)-float(ngrid))/(float(ngrid)*dx)
           endif
           final_norm = final_norm+abs(psi(i,j))*abs(psi(i,j))*(kx-kx_old)
!           psi(i,:) = psi(i,:)/(60.d0/2048.d0) 
           E_kin = E_kin + 0.5d0*abs(psi(i,j))*kx*kx*abs(psi(i,j))*(kx-kx_old)/redmass
!           psi(i,:) = psi(i,:)*(60.d0/2048.d0) 
!           write(*,*) (kx_old-kx) 
           kx_old = kx   
        enddo
      enddo
     write(12,*) final_norm
      E_kin = E_kin/final_norm
     !FFT backward: Psi(k) --> Psi(x)
      ! useless to allocate and deallocate
      do i = 1,nstate
        call fourier(0,ngrid,psi(i,:)) !initialize
        call fourier(-1,ngrid,psi(i,:)) !do backward FFT
        call fourier(2,ngrid,psi(i,:)) !deallocate
      enddo
      end subroutine
!======================================================================================
      endmodule
!=============================================================
