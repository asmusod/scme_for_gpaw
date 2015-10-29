      subroutine calcRho(rho, ra, nM, a, a2, rMax2)

      implicit none
      include '../commonblks/parameters.cmn'
      include '../commonblks/compotent.cmn'

      real*8 ra(maxCoo), a(3), a2(3)
      integer nM, i, j, k, iOi, iOj
      real*8 rho(maxCoo/3), dr(3), dr1, dr2, f, rMax2
      
      do i = 1, nM
         rho(i) = 0.d0
      end do

      do i = 1, nM-1
         iOi = 3 * (2*nM + i - 1)
         do j = i+1, nM
            iOj = 3 * (2*nM + j - 1)

            do k = 1, 3
               dr(k) = ra(iOj+k) - ra(iOi+k)
               if (dr(k) .gt. a2(k)) then
                  dr(k) = dr(k) - a(k)
               else if (dr(k) .lt. -a2(k)) then
                  dr(k) = dr(k) + a(k)
               end if
               dr(k) = dr(k)
            end do

            dr2 = dr(1)**2 + dr(2)**2 + dr(3)**2
c            if (dr2 .gt. rMax2) cycle
            dr1 = sqrt(dr2)

c            if (dr1 .lt. potpar(12,1)) then
               call dens(dr1, f)
               
               rho(i) = rho(i) + f
               rho(j) = rho(j) + f
c            end if
c            if (i.eq.1) print '(3f20.5)', dr1, 1/dr1**3, rho(i)
         end do
      end do

      return
      end
c-----------------------------------------------------------------------
      subroutine calcAmp(rho, Amp, dAmp, nM)

      implicit none
      include '../commonblks/parameters.cmn'

      integer nM, i, j, k
      real*8 rho(maxCoo/3), dr(3), dr1, dr2, f, df, h, dh, rMax2
      real*8 dAmp(maxCoo/3), Amp(maxCoo/3), rhoAvg


      do i = 1, nM

c         rhoAvg = rhoAvg + rho(i)
         call hh(rho(i), Amp(i), dAmp(i))

c         print *, rho(i), Amp(i), dAmp(i)
c         stop

      end do

c      print *, rhoAvg / nM
c      stop

      return
      end
c-----------------------------------------------------------------------
c     Density
c
      subroutine dens(r, f)

      implicit none
      real*8 r, f, A, alpha, beta
      parameter (A=2.5d5, alpha=1.5d0, beta=3.d0)

      f = A * exp(-r/alpha) / r**beta

c      f = 53000.

      return
      end
c-----------------------------------------------------------------------
c     Derivative of the density
c
      subroutine dDens(r, f)

      implicit none
      real*8 r, f, A, alpha, beta
      parameter (A=2.5d5, alpha=1.5d0, beta=3.d0)

      f =- A * exp(-r/alpha) / r**(beta+1) * (beta + r/alpha)

c      f = 0.
      return
      end
c-----------------------------------------------------------------------
      subroutine hh(r, Amp, dAmp)

      implicit none
      include '../commonblks/parameters.cmn'
      include '../commonblks/compotent.cmn'

      integer n, i
      real*8 r, rp, Amp, dAmp, a(7), p1, p2
      character*40 envdata

c interpolate hexamer value
c$$$      data a, n / 0.004537832126562685, -3.96176669501045e-6,
c$$$     $     1.211268696438834e-9, -1.616683298632205e-13,
c$$$     $     1.083795945774495e-17, -2.748153055606425e-22,
c$$$     $     0., 5/

C Commented by Fer
C     data a, n / 0.005660034555240494, -9.78520357225129e-6,
C    -   6.075330570664922e-9, -1.681333815254198e-12,
C    -   2.197384518526545e-16, -1.047989953241298e-20,
C    $     0., 5/

C New values found for reparametrized potential
C     data a, n
C    &          /  1.5122785647D-02,
C    &            -1.9608006596D-05,
C    &             6.7252150745D-09,
C    &            -1.9673147232D-13,
C    &            -2.4522945269D-17,
C    &             0.0000000000D+00,
C    &             0.0000000000D+00,
C    &             4                 /
      data a, n
     &          /  0.0000000000D+00,
     &             0.0000000000D+00,
     &             0.0000000000D+00,
     &             0.0000000000D+00,
     &             0.0000000000D+00,
     &             0.0000000000D+00,
     &             0.0000000000D+00,
     &             0                 /


      INCLUDE '../commonblks/potparam.cmn' 
    

CAA C Debug
CAA C Read the parameters from the environment
CAA      call getenv('P_1',envdata)
CAA      read(envdata,*) a(1)
CAA      call getenv('P_2',envdata)
CAA      read(envdata,*) a(2)
CAA      call getenv('P_3',envdata)
CAA      read(envdata,*) a(3)
CAA      call getenv('P_4',envdata)
CAA      read(envdata,*) a(4)
CAA      call getenv('P_5',envdata)
CAA      read(envdata,*) a(5)
CAA      call getenv('P_6',envdata)
CAA      read(envdata,*) a(6)
CAA      call getenv('P_7',envdata)
CAA      read(envdata,*) a(7)
CAA      call getenv('N',envdata)
CAA      read(envdata,*) n

       a(1)=fsti(7)
       a(2)=fsti(8)
       a(3)=fsti(9)
       a(4)=fsti(10)
       a(5)=fsti(11)
       a(6)=fsti(12)
       a(7)=fsti(13)
       n=INT(fsti(14))

CAA       write(*,*) 'P_1 ',a(1) 
CAA       write(*,*) 'P_2 ',a(2)
CAA       write(*,*) 'P_3 ',a(3)
CAA       write(*,*) 'P_4 ',a(4)
CAA       write(*,*) 'P_5 ',a(5)
CAA       write(*,*) 'P_6 ',a(6)
CAA       write(*,*) 'P_7 ',a(7)
CAA       write(*,*) 'N   ',n
CAA       pause


c      if (r .lt. 1465.35283) then
      if (r .lt. 1600.d0) then
         Amp = 0.d0
         dAmp = 0.d0
c      else if (r .gt. 8275.18577201142) then
      else if (r .gt. 8000.) then
C     else if (r .gt. 10000.) then
C Commented by Fer
C        Amp = 0.024d0/2d0
C New values found for reparametrized potential
         Amp = 0.1750D0/2.0D0
C        Amp = 0.0610D0/2.0D0
         dAmp = 0.d0
      else
         Amp = a(1)
         dAmp = 0.d0
         rp = 1.d0

         do i = 1, n
            dAmp = dAmp + dble(i) * a(i+1) * rp
            rp = rp * r
            Amp = Amp + a(i+1) * rp
         end do
      end if

c      print *, r, Amp, dAmp
c      stop

      return
      end
