      subroutine calcEnergy(dpole, qpole, opole, hpole, d1v, d2v, d3v,
     $     d4v, nM, uTot)

      implicit none
      include '../commonblks/parameters.cmn'
      
c     Work multipoles. They start unpolarized and with the induction
c     loop we induce dipoles and quadrupoles.
      real*8 dpole(3,maxCoo/3), qpole(3,3,maxCoo/3)
      real*8 opole(3,3,3,maxCoo/3), hpole(3,3,3,3,maxCoo/3)

c     High order derivatives of the potential
      real*8 d1v(3,maxCoo/3), d2v(3,3,maxCoo/3), d3v(3,3,3,maxCoo/3)
      real*8 d4v(3,3,3,3,maxCoo/3), uTot, ud, uq, uo, uh, du

      integer nM

      integer n, i, j, k, l, s
      real*8 u

      uTot = 0.d0
      ud = 0.d0
      uq = 0.d0
      uo = 0.d0
      uh = 0.d0

      do n = 1, nM

         do i = 1, 3
c     Energy of the dipole
            du = d1v(i,n) * dpole(i,n)
            uTot = uTot + du
            ud = ud + du
            do j = 1, 3
c     Energy of the quadrupole
               du = d2v(j,i,n)*qpole(j,i,n) / 3.d0
               uTot = uTot + du
               uq = uq + du
               do k = 1, 3
c     Energy of the octopole
                  du = d3v(k,j,i,n)*opole(k,j,i,n) / 15.d0
                  uTot = uTot + du
                  uo = uo + du

                  do l = 1, 3
c     Energy of the hexadecapole
                     du = d4v(l,k,j,i,n)*hpole(l,k,j,i,n) / 105.d0 
                     uTot = uTot + du
                     uh = uh + du
                  end do
               end do
            end do

         end do
      end do
      uTot = uTot / 2.d0
      
      ud = ud / 2.d0
      uq = uq / 2.d0
      uo = uo / 2.d0
      uh = uh / 2.d0

      return
      end
c-----------------------------------------------------------------------
      subroutine calcEnergyI(dpole, dpole0, qpole, qpole0, opole, hpole,
     $     d1v, d2v, d3v, d4v, nM, uTot)  

      implicit none
      include '../commonblks/parameters.cmn'
      
c     Work multipoles. They start unpolarized and with the induction
c     loop we induce dipoles and quadrupoles.
      real*8 dpole(3,maxCoo/3), qpole(3,3,maxCoo/3)
      real*8 dpole0(3,maxCoo/3), qpole0(3,3,maxCoo/3)
      real*8 opole(3,3,3,maxCoo/3), hpole(3,3,3,3,maxCoo/3)

c     High order derivatives of the potential
      real*8 d1v(3,maxCoo/3), d2v(3,3,maxCoo/3), d3v(3,3,3,maxCoo/3)
      real*8 d4v(3,3,3,3,maxCoo/3), uTot

      integer nM

      integer n, i, j, k, l, s
      real*8 u

      uTot = 0.d0
      do n = 1, nM

         do i = 1, 3
c     Energy of the dipole
            uTot = uTot + d1v(i,n) * (dpole(i,n) - dpole0(i,n))
            do j = 1, 3
c     Energy of the quadrupole
               uTot = uTot + d2v(j,i,n)*(qpole(j,i,n)-qpole0(j,i,n)) / 
     $              3.d0 
            end do
            
         end do
      end do
      uTot = uTot / 2.d0

      return
      end
c-----------------------------------------------------------------------
      subroutine polarizationEnergy1(hp, dpole, dpole0, qpole, qpole0,
     $     d1v, d2v, nM, uPol)

      implicit none
      include '../commonblks/parameters.cmn'
      
c     Work multipoles. They start unpolarized and with the induction
c     loop we induce dipoles and quadrupoles.
      real*8 dpole(3,maxCoo/3), qpole(3,3,maxCoo/3)
      real*8 dpole0(3,maxCoo/3), qpole0(3,3,maxCoo/3)

c     High order derivatives of the potential
      real*8 d1v(3,maxCoo/3), d2v(3,3,maxCoo/3), uPol
      real*8 hp(3,3,3,maxCoo/3) 

      integer i, j, k, nM

      integer n
      real*8 u

      uPol = 0.d0
      do n = 1, nM
c     Dipole polarization
         u = -0.5d0 * ( (dpole(1,n)-dpole0(1,n)) * d1v(1,n) + 
     $        (dpole(2,n)-dpole0(2,n)) * d1v(2,n) +
     $        (dpole(3,n)-dpole0(3,n)) * d1v(3,n) )
         uPol = uPol + u

      end do

      return
      end
c-----------------------------------------------------------------------
      subroutine polarizationEnergy(dd, dq, qq, hp, d1v, d2v, nM, uPol)

      implicit none
      include '../commonblks/parameters.cmn'
      
c     Work multipoles. They start unpolarized and with the induction
c     loop we induce dipoles and quadrupoles.

c     High order derivatives of the potential
      real*8 d1v(3,maxCoo/3), d2v(3,3,maxCoo/3), uPol
      real*8 dd(3,3,maxCoo/3), dq(3,3,3,maxCoo/3), hp(3,3,3,maxCoo/3)
      real*8 qq(3,3,3,3,maxCoo/3)

      integer i, j, k, l, nM

      integer n
      real*8 u

      uPol = 0.d0
      do n = 1, nM

c     Dipole-dipole polarization
         do i = 1, 3
            do j = 1, 3
               uPol = uPol + 0.5d0 * dd(i,j,n) *d1v(i,n)*d1v(j,n)

c     Dipole-quadrupole polarization
               do k = 1, 3
                  uPol = uPol + dq(i,j,k,n) *d1v(i,n)*d2v(j,k,n) / 3.d0

c     Quadrupole-quadrupole polarization
                  do l = 1, 3
                     uPol = uPol + qq(i,j,k,l,n) *d2v(i,j,n)*d2v(k,l,n)
     $                    / 6.d0
                  end do

c     first hyperpolarization
c                  uPol = uPol - hp(i,j,k,n) * d1v(i,n)*d1v(j,n)*d1v(k,n)
c     $                 / 6.d0 
               end do
            end do
         end do
      end do

      return
      end
