      subroutine shaker2()

      implicit real*8 (a-h,o-z)

      include '../commonblks/parameters.cmn'
      include '../commonblks/comtime.cmn'
      include '../commonblks/comconf.cmn'

      integer nMolecules, nOxygens, nHydrogens

      real*8 r1(3), r2(3), r3(3), v1(3), v2(3), v3(3), dt
      real*8 m1, m2, m3, invMass
      
C Debug
C     write(10,FMT='(1X,A)') ' Entering shaker2'

      if(nFPI .eq. 0) nimFPI(1)=1
cERB      nImages = 10
      nImages = nimFPI(1)
      if (nImages .gt. 2) then
         imageIni = 1
         imageFinal = nImages
      else
         imageIni = 1
         imageFinal = 1
      end if
      nOxygens = nAtms(2)
      nHydrogens = nAtms(1)
      nMolecules = nOxygens/nImages
      dt = stpsz

      do i = 1, nMolecules
c         print *, 'Shaking molecule # ', i, 'of ', nMolecules
cERB         do image = 1, nImages
         do image = imageIni, imageFinal
c$$$            if ((tag .ne. 1) .and. ((tag .lt. 500) .or. (tag .gt. 539)))
c$$$     $           then
            do j = 1, 3
               index1 = 6*nImages*(i-1) + 3*(image-1) + j
               index2 = 6*nImages*(i-1) + 3*nImages + 3*(image-1) + j
               indexO = 3*nImages*(i-1) + 3*nHydrogens + 
     $              3*(image-1) + j
               r1(j) = RA(index1)
               r2(j) = RA(index2)
               r3(j) = RA(indexO)
               v1(j) = VA(index1)
               v2(j) = VA(index2)
               v3(j) = VA(indexO)  
            end do

            call reconstruct(r1, r2, r3)

C Debug
C     write(10,FMT='(1X,A)') ' r1, r2, r3 '
C     write(10,FMT='(1X,3F20.14)') (r1(p),p=1,3)
C     write(10,FMT='(1X,3F20.14)') (r2(p),p=1,3)
C     write(10,FMT='(1X,3F20.14)') (r3(p),p=1,3)
C     write(10,FMT='(1X,A)') ' '

            m1 = sscrh(1)
            m2 = sscrh(1)
            itagOx = itag(nHydrogens+nImages*(i-1)+image)
            m3 = invMass(itagOx, sscrh(2))
              
C Debug
C     write(10,FMT='(1X,A)') ' m1, m2, m3 '
C     write(10,FMT='(1X,3F20.14)') m1, m2, m3
C     write(10,FMT='(1X,A)') ' '

c-----------------------------------------------------
c$$$            print '(A)', 'Before shaker2'
c$$$            print '(A,3g17.8)', 'r1 : ',(r1(j), j=1,3)
c$$$            print '(A,3g17.8)', 'r2 : ',(r2(j), j=1,3)
c$$$            print '(A,3g17.8)', 'r3 : ',(r3(j), j=1,3)
c$$$            
c$$$            print '(A,3g17.8)', 'v1 : ',(v1(j), j=1,3)
c$$$            print '(A,3g17.8)', 'v2 : ',(v2(j), j=1,3)
c$$$            print '(A,3g17.8)', 'v3 : ',(v3(j), j=1,3)
c$$$            print *
c$$$
c$$$         pause 'Shaker2 : Stopped at SHAKER2'
c-----------------------------------------------------
            
C Debug
C     write(10,FMT='(1X,A)') ' v1, v2, v3 '
C     write(10,FMT='(1X,3F20.14)') (v1(p),p=1,3)
C     write(10,FMT='(1X,3F20.14)') (v2(p),p=1,3)
C     write(10,FMT='(1X,3F20.14)') (v3(p),p=1,3)
C     write(10,FMT='(1X,A)') ' '

            call shake2(r1, r2, r3, v1, v2, v3, m1, m2, m3, dt, i)
c                ^^^^^ sscrh contains dt/(2 m) for each type of atom

C Debug
C     write(10,FMT='(1X,A)') ' v1, v2, v3 '
C     write(10,FMT='(1X,3F20.14)') (v1(p),p=1,3)
C     write(10,FMT='(1X,3F20.14)') (v2(p),p=1,3)
C     write(10,FMT='(1X,3F20.14)') (v3(p),p=1,3)
C     write(10,FMT='(1X,A)') ' '

            do j = 1, 3
               index1 = 6*nImages*(i-1) + 3*(image-1) + j
               index2 = index1 + 3*nImages
               indexO = 3*nImages*(i-1) + 3*nHydrogens + 
     $              3*(image-1) + j
               VA(index1) = VA(index1) + v1(j)
               VA(index2) = VA(index2) + v2(j)
               VA(indexO) = VA(indexO) + v3(j)

               FA(index1) = FA(index1) + v1(j) / sscrh(1)
               FA(index2) = FA(index2) + v2(j) / sscrh(1)
               FA(indexO) = FA(indexO) + v3(j) / sscrh(2)
            end do

c-----------------------------------------------------
c$$$            print '(A)', 'After shaker2'
c$$$            print '(A,3g17.8)', 'r1 : ',(r1(j), j=1,3)
c$$$            print '(A,3g17.8)', 'r2 : ',(r2(j), j=1,3)
c$$$            print '(A,3g17.8)', 'r3 : ',(r3(j), j=1,3)
c$$$         
c$$$            print '(A,3g17.8)', 'v1 : ',(v1(j), j=1,3)
c$$$            print '(A,3g17.8)', 'v2 : ',(v2(j), j=1,3)
c$$$            print '(A,3g17.8)', 'v3 : ',(v3(j), j=1,3)
c$$$            
c$$$            pause 'Shaker2 : Stopped at SHAKER2'
c-----------------------------------------------------

c$$$         end if
         end do
      end do

C Debug
C     write(10,FMT='(1X,A)') ' Leaving shaker2'

      end

c------------------------------------------------------------------------
      subroutine shake2(r1, r2, r3, v1, v2, v3, m1, m2, m3, dt, 
     $     moleNumber) 

      real*8 m(3,3), r1(3), r2(3), r3(3), r12(3), r23(3), r31(3) 
      real*8 v12(3), v23(3), v31(3), v1(3), v2(3), v3(3)
      real*8 b(3), dt, dummy, dot, bb(3), mm(3,3), c(3)
      real*8 m1, m2, m3
      integer indx(3)

C Debug
C     write(10,FMT='(1X,A)') ' Entering shake2'

C Debug
C     write(10,FMT='(1X,A)') ' r1, r2, r3 '
C     write(10,FMT='(1X,3F20.14)') (r1(p),p=1,3)
C     write(10,FMT='(1X,3F20.14)') (r2(p),p=1,3)
C     write(10,FMT='(1X,3F20.14)') (r3(p),p=1,3)
C     write(10,FMT='(1X,A)') ' '

C Debug
C     write(10,FMT='(1X,A)') ' v1, v2, v3 '
C     write(10,FMT='(1X,3F20.14)') (v1(p),p=1,3)
C     write(10,FMT='(1X,3F20.14)') (v2(p),p=1,3)
C     write(10,FMT='(1X,3F20.14)') (v3(p),p=1,3)
C     write(10,FMT='(1X,A)') ' '

C Debug
C     write(10,FMT='(1X,A)') ' m1, m2, m3 '
C     write(10,FMT='(1X,4F20.14,I5)') m1, m2, m3
C     write(10,FMT='(1X,A)') ' '

C Debug
C     write(10,FMT='(1X,A)') ' dt, moleNumber '
C     write(10,FMT='(1X,4F20.14,I5)') dt, moleNumber
C     write(10,FMT='(1X,A)') ' '

c$$$      print '(A,3g17.8)', 'r1 : ',(r1(j), j=1,3)
c$$$      print '(A,3g17.8)', 'r2 : ',(r2(j), j=1,3)
c$$$      print '(A,3g17.8)', 'r3 : ',(r3(j), j=1,3)
c$$$      
c$$$      print *
c$$$      print '(A,3g17.8)', 'v1 : ',(v1(j), j=1,3)
c$$$      print '(A,3g17.8)', 'v2 : ',(v2(j), j=1,3)
c$$$      print '(A,3g17.8)', 'v3 : ',(v3(j), j=1,3)
c$$$
c$$$      pause 'Shaker2 : Stopped at SHAKE 2'

      call calcDiffe(r1, r2, r3, r12, r23, r31)
      call calcDiffe(v1, v2, v3, v12, v23, v31)

C Debug
C     write(10,FMT='(1X,A)') ' r12, r23, r31 '
C     write(10,FMT='(1X,3F20.14)') (r12(p),p=1,3)
C     write(10,FMT='(1X,3F20.14)') (r23(p),p=1,3)
C     write(10,FMT='(1X,3F20.14)') (r31(p),p=1,3)
C     write(10,FMT='(1X,A)') ' '

C Debug
C     write(10,FMT='(1X,A)') ' v12, v23, v31 '
C     write(10,FMT='(1X,3F20.14)') (v12(p),p=1,3)
C     write(10,FMT='(1X,3F20.14)') (v23(p),p=1,3)
C     write(10,FMT='(1X,3F20.14)') (v31(p),p=1,3)
C     write(10,FMT='(1X,A)') ' '

c     *** Note that m1, m2 and m3 contain dt/(2 mi) that is the inverse of
c     twice the mass times dt ***
      m(1,1) =  dot(r12, r12) * (m1+m2)
      m(2,2) =  dot(r23, r23) * (m2+m3)
      m(3,3) =  dot(r31, r31) * (m3+m1)
      m(1,2) = -dot(r12, r23) * m2
      m(2,1) =  m(1,2)
      m(1,3) = -dot(r12, r31) * m1
      m(3,1) =  m(1,3)
      m(2,3) = -dot(r23, r31) * m3
      m(3,2) =  m(2,3)

      b(1) = -dot(r12, v12)
      b(2) = -dot(r23, v23)
      b(3) = -dot(r31, v31)

      do i = 1, 3
         bb(i) = b(i)
         do j = 1, 3
            mm(i,j) = m(i,j)
         end do
      end do

 
c$$$      print '(A)', 'M matrix'
c$$$      print *, 'Molecule number ', moleNumber
c$$$      do i = 1, 3
c$$$         print '(3g17.8)', (m(i,j), j=1,3)
c$$$      end do
c$$$      print *
c$$$
c$$$      print '(A)', 'Right side'
c$$$      print '(3g17.8)', (b(i), j=1,3)

C Debug
C     write(10,FMT='(1X,A)') ' m '
C     write(10,FMT='(1X,3F20.14)') (m(1,p),p=1,3)
C     write(10,FMT='(1X,3F20.14)') (m(2,p),p=1,3)
C     write(10,FMT='(1X,3F20.14)') (m(3,p),p=1,3)
C     write(10,FMT='(1X,A)') ' '

C Debug
C     write(10,FMT='(1X,A)') ' indx '
C     write(10,FMT='(1X,3I5)') indx
C     write(10,FMT='(1X,A)') ' '

C Debug
C     write(10,FMT='(1X,A)') ' dummy '
C     write(10,FMT='(1X,F20.14)') dummy
C     write(10,FMT='(1X,A)') ' '

c$$$      pause 'Shaker2 :  Stopped in SHAKE2 before ludcmp'
      call ludcmp(m, 3, 3, indx, dummy)

C Debug
C     write(10,FMT='(1X,A)') ' m '
C     write(10,FMT='(1X,3F20.14)') (m(1,p),p=1,3)
C     write(10,FMT='(1X,3F20.14)') (m(2,p),p=1,3)
C     write(10,FMT='(1X,3F20.14)') (m(3,p),p=1,3)
C     write(10,FMT='(1X,A)') ' '

C Debug
C     write(10,FMT='(1X,A)') ' indx '
C     write(10,FMT='(1X,3I5)') indx
C     write(10,FMT='(1X,A)') ' '

C Debug
C     write(10,FMT='(1X,A)') ' dummy '
C     write(10,FMT='(1X,F20.14)') dummy
C     write(10,FMT='(1X,A)') ' '

C Debug
C     write(10,FMT='(1X,A)') ' b '
C     write(10,FMT='(1X,3F20.14)') b
C     write(10,FMT='(1X,A)') ' '

c$$$      pause 'Shaker2 :  Stopped in SHAKE2 after ludcmp'
      call lubksb(m, 3, 3, indx, b)
    
C Debug
C     write(10,FMT='(1X,A)') ' b '
C     write(10,FMT='(1X,3F20.14)') b
C     write(10,FMT='(1X,A)') ' '

c$$$      print *
c$$$      print '(A)', 'In Shake2'
c$$$      do i = 1, 3
c$$$         c(i) = 0.d0
c$$$         do j = 1, 3
c$$$            c(i) = c(i) + mm(i,j) * b(j)
c$$$         end do
c$$$         print '(3g20.13)', c(i), bb(i), c(i)-bb(i)
c$$$      end do

c      print *
c      print '(A)', 'Corrections to the velocities in shaker2'

      do i = 1, 3
c$$$         v1(i) = v1(i) + m1 * ( b(1) * r12(i) - b(3) * r31(i))
c$$$         v2(i) = v2(i) + m2 * (-b(1) * r12(i) + b(2) * r23(i))
c$$$         v3(i) = v3(i) + m3 * (-b(2) * r23(i) + b(3) * r31(i))

         v1(i) = m1 * ( b(1) * r12(i) - b(3) * r31(i))
         v2(i) = m2 * (-b(1) * r12(i) + b(2) * r23(i))
         v3(i) = m3 * (-b(2) * r23(i) + b(3) * r31(i))
c$$$         print '(3g17.8)',
c$$$     $        m1 * ( b(1) * r12(i) - b(3) * r31(i)),
c$$$     $        m2 * (-b(1) * r12(i) + b(2) * r23(i)),
c$$$     $        m3 * (-b(2) * r23(i) + b(3) * r31(i))
      end do
c$$$      print *
c$$$
c$$$      call calcDiffe(v1, v2, v3, v12, v23, v31)
c$$$      print '(A,3g17.8)', 'dot :', dot(r12,v12), dot(r23,v23),
c$$$     $     dot(r31,v31)
c$$$      pause 'stopped at SHAKER2 to check orthogonality'

C Debug
C     write(10,FMT='(1X,A)') ' Leaving shake2'

      return
      end







