      subroutine coreInt(ra, fa, uCore, nM, a, a2)

      implicit none
      include '../commonblks/parameters.cmn'
      include '../commonblks/compotent.cmn'

      integer nM, ai(10), kk, jj
      real*8 t1, t2, df, uCore, r, a(3), a2(3), avg
      real*8 ra(maxCoo), fa(maxCoo), c1, c2, c3, c4, c5,t11, t12, t13
      real*8 Amp(maxCoo/3), dAmp(maxCoo/3), rho(maxCoo/3), f, rMax2
      real*8 Ri1(3), Ri(3), Ri2(3), Rj1(3), Rj(3), Rj2(3)
      real*8 vRi1j1(3), vRi2j2(3), vRij(3)
      real*8 Rij
      real*8 S
      real*8 uR_An
      real*8 an_1, an_2, an_3, an_4
      real*8 df1, df2
      real*8 c5_r
      integer*4 Pt_i, Pt_i1, Pt_i2
      integer*4 Pt_j, Pt_j1, Pt_j2
      character*40 envdata

c     Using switching at short distances.
C     parameter(c1 = -2.14409d0,
C    $          c2 = 21773.8477335795*1.02d0,
C    $          c3 = -3.24106d0, 
C    $          c4 = c3/1.325d0)

      integer i, j, k, n, m, iOn, iOm, iOi, iOj, p
      real*8 dr(3)

      INCLUDE '../commonblks/potparam.cmn'
CAA
CAA C Read the core repulsion parameters from the environment
CAA      call getenv('CF_1',envdata)
CAA	 read(envdata,*) c1
CAA      call getenv('CF_2',envdata)
CAA      read(envdata,*) c2
CAA      call getenv('CF_3',envdata)
CAA      read(envdata,*) c3
CAA      call getenv('CF_4',envdata)
CAA      read(envdata,*) c4
CAA      call getenv('CF_5',envdata)
CAA      read(envdata,*) c5_r

      c1=fsti(2)
      c2=fsti(3)
      c3=fsti(4)
      c4=fsti(5)
      c5_r=fsti(6)

      call calcRho(rho, ra, nM, a, a2, rMax2)
      call calcAmp(rho, Amp, dAmp, nM)

C Debug
C     write(10,FMT='(A,I3,A,F16.10)')
C    &      'Rho(',int(nM/2),') = ', rho(int(nM/2))
C     do i=1,nM
C       write(10,FMT='(A,I3,A,F16.10)')
C    &        'Rho(',i,') = ', rho(i)
C     end do
C     write(10,FMT='(A,I3,A,F16.10)')
C    &      'Amp(',int(nM/2),') = ', Amp(int(nM/2))
C     if ( nM .le. 10 ) then
C       do i=1,nM
C         write(10,FMT='(I3,3X,2F16.10)')  i, rho(i), Amp(i)
C       end do
C     else
C       write(10,FMT='(I3,3X,2F16.10)') int(nM/2), rho(int(nM/2)),
C    &                                          Amp(int(nM/2))
C     end if

      uCore = 0.0D0

      do n = 1, nM-1

C Get the index of the first O atom
        iOn = 3 * (2*nM + n - 1)

        do m = n+1, nM

C Get the index of the second O atom
          iOm = 3 * (2*nM + m - 1)

C Adjust the O-O distance for the PBC's
          do i = 1, 3
            dr(i) = ra(iOm+i) - ra(iOn+i)
            if     (dr(i) .gt. a2(i)) then
              dr(i) = dr(i) - a(i)
            elseif (dr(i) .lt. -a2(i)) then
              dr(i) = dr(i) + a(i)
            end if
          end do

          r = dsqrt(dr(1)**2 + dr(2)**2 + dr(3)**2)

          t11 = dexp(c3*r)
          t12 = dexp(c3*r/c4)
          t13 = r**c1

C Debug
          if ( c5_r .ge. 0.0D0 ) then
            c5 = c5_r
          else
            c5 = (Amp(n)+Amp(m))
          end if

          t1 = c2*t13*(t11 + c5*t12)
C         t1 = c2*t13*t11

          uCore = uCore + t1

          t2 = c2*t13*t12

          df = (c2*(c1/r+c3)*t13*t11 + c2*c5*(c1/r+c3/c4)*t13*t12) / r
C         df = (c2*(c1/r+c3)*t13*t11) / r

          do i = 1, 3
            fa(iOn+i) = fa(iOn+i) + df * dr(i)
            fa(iOm+i) = fa(iOm+i) - df * dr(i)
          end do

c     Derivative of the embedding part (Amplitude)
c
C Debug
            if (c5_r .lt. 0.0D0 ) then
            do j = 1, nM
               if (j .eq. n) goto 11
               iOj = 3 * (2*nM + j - 1)

               do k = 1, 3
                  dr(k) = ra(iOj+k) - ra(iOn+k)
                  if (dr(k) .gt. a2(k)) then
                     dr(k) = dr(k) - a(k)
                  else if (dr(k) .lt. -a2(k)) then
                     dr(k) = dr(k) + a(k)
                  end if
               end do

               r = sqrt(dr(1)**2 + dr(2)**2 + dr(3)**2)
               
               call dDens(r, f)
                  
               df = dAmp(n) * f * t2 / r
               do k = 1, 3
                  fa(iOn+k) = fa(iOn+k) + df * dr(k)
                  fa(iOj+k) = fa(iOj+k) - df * dr(k)
               end do
 11         end do

            do j = 1, nM
               if (j .eq. m) goto 12
               iOj = 3 * (2*nM + j - 1)
               
               do k = 1, 3
                  dr(k) = ra(iOj+k) - ra(iOm+k)
                  if (dr(k) .gt. a2(k)) then
                     dr(k) = dr(k) - a(k)
                  else if (dr(k) .lt. -a2(k)) then
                     dr(k) = dr(k) + a(k)
                  end if
               end do

               r = sqrt(dr(1)**2 + dr(2)**2 + dr(3)**2)
               call dDens(r, f)

               df = dAmp(m) * f * t2 / r
               do k = 1, 3
                  fa(iOm+k) = fa(iOm+k) + df * dr(k)
                  fa(iOj+k) = fa(iOj+k) - df * dr(k)
               end do
 12         end do
            end if

         end do
      end do

C Read the core anisotropy parameters from the environment
C     call getenv('AN_1',envdata)
C     read(envdata,*) an_1
C     call getenv('AN_2',envdata)
C     read(envdata,*) an_2
C     call getenv('AN_3',envdata)
C     read(envdata,*) an_3
C     call getenv('AN_4',envdata)
C     read(envdata,*) an_4

C Adding the missing anisotropic exchange repulsion term

C     uR_An = 0.0D0

C     do i=1,nM-1

C Get the coordinates for molecule i
C       Pt_i1 = 6*(i-1)
C       Pt_i  = 3*(i-1+2*nM)
C       Pt_i2 = 3+6*(i-1)

C       do p=1,3
C         Ri1(p) = ra(Pt_i1 + p)
C         Ri(p)  = ra(Pt_i  + p)
C         Ri2(p) = ra(Pt_i2 + p)
C       end do

C       do j=i+1,nM

C Get the coordinates for molecule j
C         Pt_j1 = 6*(j-1)
C         Pt_j  = 3*(j-1+2*nM)
C         Pt_j2 = 3+6*(j-1)

C         do p=1,3
C           Rj1(p) = ra(Pt_j1 + p)
C           Rj(p)  = ra(Pt_j  + p)
C           Rj2(p) = ra(Pt_j2 + p)
C         end do

C Calculate the vectors
C         Rij   = 0.0D0

C         do p=1,3
C           vRi1j1(p) = Ri1(p) - Rj1(p)
C           vRi2j2(p) = Ri2(p) - Rj2(p)
C           vRij(p)   = Ri(p)  - Rj(p)
C           Rij   = Rij   + vRij(p)**2
C         end do

C         Rij = dsqrt(Rij)

C Calculate the energy contributions

C         S = 0.0D0

C         do p=1,3
C           S = S + 
C    &          ( vRi1j1(p) + vRi2j2(p) - 2.0D0*vRij(p) ) * vRij(p)
C         end do

C         uR_An = uR_An + an_1*S/Rij

C Calculate the derivatives

C With respect to the H positions:

C         do p=1,3
C           df1 = an_1*vRij(p)/Rij
C           fa(Pt_i1 + p) = fa(Pt_i1 + p) - df1
C           fa(Pt_i2 + p) = fa(Pt_i2 + p) - df1
C           fa(Pt_j1 + p) = fa(Pt_j1 + p) + df1
C           fa(Pt_j2 + p) = fa(Pt_j2 + p) + df1
C         end do

C With respect to the O positions:
          
C         do p=1,3
C           df2 = an_1*( -vRij(p)*S/Rij**2 + 
C    &                    vRi1j1(p) + vRi2j2(p) - 4.0D0*vRij(p) )/Rij
C           fa(Pt_i + p) = fa(Pt_i + p) - df2
C           fa(Pt_j + p) = fa(Pt_j + p) + df2
C         end do

C Debug
C     write(10,FMT='(F16.10)') uR_An
C     write(10,FMT='(F16.10)') 
C    &           (Ri1(1)-Rj1(1))*df11 + (Ri1(1)-Rj2(1))*df12
C       end do
C     end do

C     uCore = uCore + uR_An

      end

C***************************************************************************

