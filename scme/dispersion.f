      subroutine dispersion(ra, fa, uDisp, nM, a, a2)

      implicit none
      include '../commonblks/parameters.cmn'

      integer nM
      real*8 t1, t2, df, r, a(3), a2(3), uDisp
      real*8 r2, r6, r7, r8, r9, r10, r11
      real*8 f6, df6, f8, df8, f10, df10
      real*8 ra(maxCoo), fa(maxCoo), C6, C8, C10

      integer n, m, iOn, iOm, i
      real*8 dr(3), sc
      parameter (C6 = 46.4430d0 * 0.597527378d0, 
     $     C8  = 1141.7000d0 * 0.167324732d0, 
     $     C10 = 33441.0000d0 * 0.046855703d0)

c     Dispersion coefficients from Watts-Coker
c      Set I
c      parameter (C6 = 37.2484d0*1.15, C8  = 224.48d0*1.15, 
c     $     C10 = 1560.92d0*1.5 * 1.15)

c      parameter (C6 = 37.2484d0, C8  = 224.48d0, 
c     $     C10 = 1560.92d0*1.5)

      uDisp = 0.d0
      do n = 1, nM-1
c      do n = 1, 4
         iOn = 3 * (2*nM + n - 1)
         do m = n+1, nM
c         do m = 5, 8
            iOm = 3 * (2*nM + m - 1)

            do i = 1, 3
               dr(i) = ra(iOm+i) - ra(iOn+i)
               if (dr(i) .gt. a2(i)) then
                  dr(i) = dr(i) - a(i)
               else if (dr(i) .lt. -a2(i)) then
                  dr(i) = dr(i) + a(i)
               end if
            end do

            r2 = dr(1)**2 + dr(2)**2 + dr(3)**2
            r = sqrt(r2)
            call Tang_Toennies(r, f6, df6, f8, df8, f10, df10)

            r6 = r**6
            r7 = r6 * r
            r8 = r7 * r
            r9 = r8 * r
            r10 = r9 * r
            r11 = r10 * r

            uDisp = uDisp - C6 / r6 * f6 - C8 / r8 * f8 - C10 / r10
     $           * f10 

            do i = 1, 3
               df = -C6 * (6.d0 * f6 / r7 - df6 / r6)
               df = df - C8 * (8.d0 * f8 / r9 - df8 / r8)
               df = df - C10 * (10.d0 * f10 / r11 - df10 / r10)
               df = df * dr(i) / r
               fa(iOn+i) = fa(iOn+i) - df
               fa(iOm+i) = fa(iOm+i) + df
            end do
         end do
      end do

      return
      end
c-----------------------------------------------------------------------
      subroutine Tang_Toennies(r, f6, df6, f8, df8, f10, df10)

      implicit none
      real*8 b, r, f6, df6, f8, df8, f10, df10
      real*8 ff6, ff8, ff10, x2, x3, x4, x5, x6, x7, x8, x9, x10, x11 
      real*8 x12, x13, x14

      integer k, fact
      real*8 x, t
c      parameter (b = 2.48d0)
      parameter (b = 4.4d0)

      x = b * r

      f6 = 1.d0
      df6 = 0.d0

      t = exp(-x)
      do k = 0, 6
         f6 = f6 - t
         df6 = df6 - t * b * (-1.d0 + k / x)
         t = t * x / (k+1)
      end do

      f8 = f6
      df8 = df6
      do k = 7, 8
         f8 = f8 - t
         df8 = df8 - t * b * (-1.d0 + k / x)
         t = t * x / (k+1)
      end do

      f10 = f8
      df10 = df8
      do k = 9, 10
         f10 = f10 - t
         df10 = df10 - t * b * (-1.d0 + k / x)
         t = t * x / (k+1)
      end do

      return
      end


