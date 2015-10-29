c-----------------------------------------------------------------------
      subroutine Tang_ToenniesN(r, f, n)

      implicit none
      real*8 r, f, b
c      parameter (b=1.92d0 / 0.529177249d0)
C     parameter (b=4.4d0)
      integer k, n
      real*8 x, t
      character*40 envdata

      INCLUDE '../commonblks/potparam.cmn'

CAA      call getenv('DL_1',envdata)
CAA      read(envdata,*) b
      b=fsti(1)

      x = b * r

      f = 1.d0
      
      t = exp(-x)
      do k = 0, n
         f = f - t
         t = t * x / (k+1)
      end do

      f = sqrt(f)

c      f = 1.d0
c      df = 0.d0
      return
      end
c-----------------------------------------------------------------------
      subroutine Tang_ToenniesNdF(r, f, df, n)

      implicit none
      real*8 r, f, df, b
c      parameter (b=1.92d0 / 0.529177249d0)
C     parameter (b=4.4d0)
      integer k, n
      real*8 x, t
      character*40 envdata

      INCLUDE '../commonblks/potparam.cmn'

CAA      call getenv('DL_1',envdata)
CAA      read(envdata,*) b
      b=fsti(1) 

      x = b * r

      f = 1.d0
      df = 0.d0
      
      t = exp(-x)
      do k = 0, n
         f = f - t
         df = df - t * b * (-1.d0 + k / x)
         t = t * x / (k+1)
      end do

      f = sqrt(f)
      df = df / (2.d0 * f)

c      f = 1.d0
c      df = 0.d0
      return
      end
