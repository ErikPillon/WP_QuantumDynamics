!=============================================================
      program wp_propag
      use mainvar
      use mainmethods

      implicit none
      real*8 :: rescaling_factor
      integer :: i, ii, n
      logical :: ex
      character(10) :: ndir

     !Read input
      read(*,*) ngrid,xmin,xmax
      read(*,*) k_ini,xcenter_ini,sigma_ini,redmass
      read(*,*) nstep, dt, nwrite
      read(*,*) nstate
      read(*,*) krange_len
      allocate(krange(krange_len))
      read(*,*) (krange(i), i=1,krange_len)

     !grid resolution
      dx=(xmax-xmin)/dble(ngrid-1)

     !Allocate variables
      call allocate_main_arrays(1)

      do ii=1,krange_len
          k_ini = krange(ii)
          ! The following line needs to be used in the case of TullyB model
          ! k_ini = SQRT(redmass*2*exp(krange(ii)))
!          call create_folder(k_ini,ex)
!          if (ex .neqv. .TRUE.) then
            write (ndir, '(F10.1)') k_ini 
!            call chdir('./' // adjustl(trim(ndir)))
           !Initialize the Gaussian WP
            call Init_Psi
            if (ii == 1) then
                call calculate_eigenvalues
            endif 
           !Build the Pot and Kin energy part of the propagator
            call Set_Operators
            
           ! ------------------------------------------------------------
           ! Adjust the time step in such a way that the total dynamics covers
           ! approximately 20 units in nuclear coordinates
           ! ------------------------------------------------------------
            ! nstep = INT(20*redmass/(k_ini*dt))
            open(unit=10, status="new",file=trim(adjustl(ndir))//".data")
            write(10,"(2i10,3a5)") ngrid,nstep/nwrite,"2nd","1st Adiabatic","2nd Adiabatic"
           !Print the initial WP
            call write_psi(0)
            !---------------------------------------------------------
            !FFT time-propagation with symmetric Trotter-splitting
            !---------------------------------------------------------
            ! nstep_ext = min(2000,max(2000,int(10*nstep/redmass)))
             do i=1,nstep
                call Propagate_One_Step
                call Evaluate_Kinetic_Energy
                ! call Evaluate_Potential_Energy

                if(mod(i,nwrite) == 0) then
                    write(10,"(f12.5,a12)") i*dt,"time"
                    do n=1,ngrid
                       write(10,"(f12.5,4e15.5)")                              &
                           xgrid(n),abs(psi(1,n))*abs(psi(1,n))*dx, &!,real(psi(n)), aimag(psi(n))
                                    abs(psi(2,n))*abs(psi(2,n))*dx, &!,real(psi(n)), aimag(psi(n)) 
                                    abs(psi_ad(1,n))*abs(psi_ad(1,n))*dx, &!,real(psi(n)), aimag(psi(n))
                                    abs(psi_ad(2,n))*abs(psi_ad(2,n))*dx!,real(psi(n)), aimag(psi(n)) 
                    enddo
                endif    
                if(mod(i,nwrite) == 0) then
                    E_pot = 0.0
                    do n = 1,ngrid
                        E_pot = E_pot + eigenvals(1,n)*abs(psi_ad(1,n))*abs(psi_ad(1,n))*dx
                        E_pot = E_pot + eigenvals(2,n)*abs(psi_ad(2,n))*abs(psi_ad(2,n))*dx
                    enddo    
                    write(11,*) E_kin, E_pot
                endif    
!               !Print wavefunction
!                if(mod(i,nwrite) == 0) then
!                    call write_psi(i)
!                    write(*,*) "time step:",i
!                endif
             
             enddo
            !-----------------------------------------------
            !End Time-propagation
            !-----------------------------------------------
               close(10)
!              call chdir('../')
!          endif  
     enddo


     !Deallocate the variables
      call allocate_main_arrays(0)
      deallocate(krange)

      end program wp_propag


!======================================================================================
    subroutine create_folder(n,ex)
!======================================================================================
    implicit none
    real(kind=8), intent(in) :: n
    logical :: ex
    character(len=6) :: dirname
    
    write (dirname, '(F4.1)') n
    write(*,*) dirname
    INQUIRE(FILE=adjustl(trim(dirname)),EXIST=ex)
    IF (ex) THEN
        write(*,'(3X, "folder " F4.1 "already exist; I will skip this step")') n
    ELSE 
        call execute_command_line('mkdir -p ./' // adjustl(trim(dirname)))
    END IF    
    end subroutine create_folder
!======================================================================================

!======================================================================================
      subroutine Vpot_TullyA(x,V)
!======================================================================================
!      
!     Fewest surface hopping algorithm as introduced in the article
!     Tully, J.C. J. Chem. Phys. (1990) 93 1061. 
!======================================================================================
      use mainvar
      implicit none
      real*8, intent(in) :: x
      real*8 :: A=0.01, B=1.6, C=0.005, D=1.0
      real*8,dimension(nstate, nstate), intent(out) :: V 

      V = 0.0
      IF (x>0) THEN
          V(1,1) = A*(1.0-exp(-B*x))
      ELSE
          V(1,1) = -A*(1.0-exp(B*x))
      END IF
      V(2,2) = -V(1,1)
      V(1,2) = C*exp(-D*x*x)
      V(2,1) = V(1,2)
      end subroutine
!======================================================================================

!======================================================================================
      subroutine Vpot_TullyB(x,V)
!======================================================================================
!      
!     Fewest surface hopping algorithm as introduced in the article
!     Tully, J.C. J. Chem. Phys. (1990) 93 1061. 
! 
!     This is the so called "Double Avoided Crossing" model
!======================================================================================
      use mainvar
      implicit none
      real*8, intent(in) :: x
      real*8 :: A=0.1, B=0.28, C=0.015, D=0.06
      real*8 :: E0=0.05
      real*8,dimension(nstate, nstate), intent(out) :: V 

      V = 0.0
      
      V(1,1) = 0.0
      V(2,2) = -A*exp(-B*x*x)+E0
      V(1,2) = C*exp(-D*x*x)
      V(2,1) = V(1,2)
      end subroutine
!======================================================================================


!======================================================================================
      subroutine calculate_eigenvalues
!======================================================================================
!      
!     Fewest surface hopping algorithm as introduced in the article
!     Tully, J.C. J. Chem. Phys. (1990) 93 1061. 
! 
!     This is the so called "Double Avoided Crossing" model
!======================================================================================
      use mainvar
      implicit none
      real*8 :: A=0.01, B=1.6 , C=0.005, D=1.0 
      real*8 :: first_part, second_part, lambda12
      real*8 :: x
      integer :: i

      do i=1,ngrid
         ! Vv=Vpot_TullyA(xgrid(i))
         ! call Vpot_TullyA(xgrid(i),Vv)
         x = xgrid(i)
         first_part = A*(1-exp(-B*abs(x)))
         second_part = C*exp(-D*x*x)
         lambda12 = sqrt(first_part*first_part+second_part*second_part)
         eigenvals(1,i) = -lambda12 
         eigenvals(2,i) = lambda12 
      enddo

      end subroutine
!======================================================================================


!======================================================================================
      subroutine Propagate_One_Step
!======================================================================================
!     FFT time-propagation with symmetric Trotter-splitting              
!     Only one time step is done in this routine
!======================================================================================
      use mainvar
      implicit none
      integer :: i
     !Half-potential step: exp(-i*V/2)*Psi(x)
      do i = 1,ngrid
        psi(:,i)=matmul(Op_Pot(:,:,i),psi(:,i))
      enddo
     !FFT forward: Psi(x) --> Psi(k)



      do i = 1,nstate
        call fourier(0,ngrid,psi(i,:)) !initialize
        call fourier(1,ngrid,psi(i,:)) !do forwad FFT
        call fourier(2,ngrid,psi(i,:)) !deallocate
      enddo
     !Kinetic energy step: exp(-i*T)*Psi(k)
      do i = 1,nstate
        psi(i,:)=Op_Kin*psi(i,:)
      enddo
     !FFT backward: Psi(k) --> Psi(x)
      ! useless to allocate and deallocate
      do i = 1,nstate
        call fourier(0,ngrid,psi(i,:)) !initialize
        call fourier(-1,ngrid,psi(i,:)) !do backward FFT
        call fourier(2,ngrid,psi(i,:)) !deallocate
      enddo

     !Half-potential step: exp(-i*V/2)*Psi(x)
      do i = 1,ngrid
        psi(:,i)=matmul(Op_Pot(:,:,i),psi(:,i))
        psi_ad(:,i)=matmul(transpose(u_matrices(:,:,i)),psi(:,i))
      enddo
    
     !Normalize the wavefunction
      call normalize_psi

      end subroutine
!======================================================================================



!======================================================================================
      subroutine Init_Psi
!======================================================================================
!     Initialize the wave packet (Generating the Gaussian WP)             
!======================================================================================
      use mainvar
      implicit none
      integer :: i
      real*8  :: gaussian
      real*8  :: sigma_new
      ! real*8 :: stddev,xmean,k_0
      ! xcenter, sigma, k_0 already imported from main module

     !setting up grid points in x-space
      do i=1,ngrid
         xgrid(i) = xmin + dx*(i-1)
      enddo

     !Generating gaussian wavepacket
      sigma_new = 20.d0/k_ini
      if (sigma_new .ge. 3.0) then
          sigma_new = 3.0
      endif    
      psi=0.d0
      do i=1,ngrid
         gaussian=exp((-(xgrid(i)-xcenter_ini)**2)/(2.0*sigma_new**2))

        !in some case the direction is important 
        !(the sign can be changed depends on the studied problem)
         psi(1,i)=gaussian*exp(jcpx*k_ini*xgrid(i))
      enddo

      call normalize_psi

      end subroutine
!======================================================================================


!======================================================================================
      subroutine Set_Operators
!======================================================================================
!     Building up the Kin and Pot operator part of the Propagator              
!              
!     Since symmetric Trotter splitting is introduced in the 
!     short time propagator, the Op_Pot is only half time step is used
!======================================================================================
      use mainvar
      implicit none
      integer :: i
      real*8 :: Tt,kx
      complex*16 :: t
      integer :: m
      real*8, dimension(nstate,nstate) :: um
      real*8, dimension(nstate) :: eigvals 
      complex*16, allocatable :: e(:,:), cVv(:,:)
      real*8, allocatable :: Vv(:,:)

      allocate ( Vv ( nstate,nstate) )
      allocate ( e ( nstate,nstate) )
      allocate ( cVv ( nstate,nstate) )
      
      m = nstate
      t = -jcpx*dt/2.d0
     !exp[-i*V(x)*dt/2]
      do i=1,ngrid
         !Vv=Vpot_TullyA(xgrid(i))
         ! In the following we need to comment/uncomment the right model
         ! call Vpot_TullyA(xgrid(i),Vv)
         call Vpot_TullyA(xgrid(i),Vv)
         call jacobi(Vv,nstate,nstate,eigvals,um,1d-10,1)
         u_matrices(:,:,i) = um
         cVv = t*Vv
         call c8mat_expm1(m, cVv, e)
         ! write(56,*) Vv, cVv,e
         Op_Pot(:,:,i) = e
      enddo

     !exp[-iT*dt]
     !
     ! T = 1/2 kx*kx/m
      do i=1,ngrid
         if(i < ngrid/2) then
           kx=2.d0*pi*float(i)/(float(ngrid)*dx)
         else
           kx=2.d0*pi*(float(i)-float(ngrid))/(float(ngrid)*dx)
         endif

         Tt=0.5d0*kx*kx/redmass

         Op_kin(i)=exp(-jcpx*Tt*dt)
         ! write(88,*) i, u_matrices(:,:,i)
      enddo
      deallocate(cVv)
      deallocate(e)
      deallocate(Vv)
      write(*,*)"imag:",jcpx,"exp(-jcpx)",exp(-jcpx*pi)

      end subroutine
!======================================================================================
 

!======================================================================================
      subroutine normalize_psi
!======================================================================================
      use mainvar
      implicit none
      real*8         :: norm
      integer        :: i,j

      norm=0.0
      do j=1,nstate
        do i=1,ngrid
            norm = norm + dx*abs(psi(j,i))*abs(psi(j,i))
          enddo
      enddo

      psi=psi/sqrt(norm)

      end subroutine
!======================================================================================

!======================================================================================
      subroutine write_psi(istep)
!======================================================================================
      use mainvar
      implicit none
      integer :: n, istep
      integer, save :: call_first=1


      if(call_first == 1) then
         write(66,"(2i10,3a5)") ngrid,nstep/nwrite,"2nd","3rd","4th"
         call_first=99
      else

         write(66,"(f12.5,a12)") istep*dt,"time"
         do n=1,ngrid
            write(66,"(f12.5,4e15.5)")                              &
                xgrid(n),abs(psi(1,n))*abs(psi(1,n))*dx, &!,real(psi(n)), aimag(psi(n))
                         abs(psi(2,n))*abs(psi(2,n))*dx, &!,real(psi(n)), aimag(psi(n)) 
                         abs(psi_ad(1,n))*abs(psi_ad(1,n))*dx, &!,real(psi(n)), aimag(psi(n))
                         abs(psi_ad(2,n))*abs(psi_ad(2,n))*dx!,real(psi(n)), aimag(psi(n)) 
         enddo

      endif

      flush(66)

      end subroutine
!======================================================================================


!======================================================================
      subroutine fourier(dir,n,c)
!======================================================================
!     Interface to the FFTPACK5 library              
!     The employed FFTPACK5 is included in another file:
!     fftpack5.f90
!======================================================================
!     c --> input and output array (complex) to be FFT transformed
!            
!     n   --> the length of the complex sequence "psi".  The method is
!             more efficient when "n" is the product of small primes.
!             a complex array of length N which contains the sequence
! 
!     dir --> switch between initialization, forward and backward FFT:
!
!     first it has to be initialized:
!     call  fourier(0,n,c)
!
!     deallocation of FFT routine:
!     call  fourier(2,n,c)
!
!     forward FFT:
!     call  fourier(1,n,c) 
!
!     backward FFT:
!     call  fourier(-1,n,c)
!
!     normalization is contained in forward FFT calling
!======================================================================
      implicit none
      integer(kind=4) :: n,dir
      integer(kind=4) :: ier
      integer(kind=4), save ::lensav,lenwrk
      real(kind=8), allocatable, dimension(:),save :: work
      real(kind=8), allocatable, dimension(:),save :: wsave
      complex(kind=8) :: c(n)

      if(dir == 1) then !Forward FFT
         call cfft1f(n,1,c,n,wsave,lensav,work,lenwrk,ier)
      elseif(dir == -1) then !Backward FFT
         call cfft1b(n,1,c,n,wsave,lensav,work,lenwrk,ier)
      elseif(dir == 0) then !Allocate the work arrays
         lenwrk=2*n
         lensav=2*n + int(log(real(n,kind=8))/log(2.0D+00))+4
         allocate(work(1:lenwrk))  
         allocate(wsave(1:lensav)) 
         call cfft1i(n,wsave,lensav,ier)
      elseif(dir == 2) then !Deallocate work arrays
         deallocate(work)
         deallocate(wsave)
      endif
      
      end subroutine fourier
!======================================================================

