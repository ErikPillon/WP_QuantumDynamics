! ========================================================
!     ________    _______
!    /  ______|  /   __  \
!    | |         |  |  \  |
!    | |_____    |  |__/  |
!    |  _____|   |  _____/
!    | |         | |
!    | |______   | |
!    \________|  \_/
!
!    erik.pillon@gmail.com
!    http://ErikPillon.github.io
!
!    February 3rd, 2021
!
!    Documentation is like sex: when 
!    it is good, it is very, very good;
!    and when it is bad, it is better than nothing.
! ========================================================

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
