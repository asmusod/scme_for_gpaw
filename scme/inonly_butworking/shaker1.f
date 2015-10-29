      subroutine shaker1(RAold)

      implicit real*8 (a-h,o-z)

      include '../commonblks/parameters.cmn'
      include '../commonblks/comtime.cmn'
      include '../commonblks/comconf.cmn'
      include '../commonblks/comenergy.cmn'    


      real*8 RAold(maxCoo)
      integer nMolecules, nOxygens, nHydrogens

      real*8 r1N(3), r2N(3), r3N(3), r1O(3), r2O(3), r3O(3)
      real*8 dt, d(3), r12(3), r23(3), r31(3), dot
      real*8 r1(3), r2(3), r3(3), cf1(3), cf2(3), cf3(3)
      real*8 m1, m2, m3, invMass
      integer steps, index1, index2, indexO, nImages, tag

c     Number of oxygens and hydrogens
      if(nFPI .eq. 0) nimFPI(1)=1
cERB      nImages = 10
      nImages = nimFPI(1)
      if (nImages .gt. 2) then
c         imageIni = 1
         imageIni = 2
c         imageFinal = nImages
         imageFinal = nImages-1
      else
         imageIni = 1
         imageFinal = 1
      end if
      nOxygens = nAtms(2)
      nHydrogens = nAtms(1)
      nMolecules = nOxygens/nImages
      dt = stpsz

      do i = 1, nMolecules
         do image = imageIni, imageFinal

c     (tag #1 => permafrost) => skip calculation
            tag = itag(nHydrogens+i)
c$$$            if ((tag .ne. 1) .and. ((tag .lt. 500) .or. (tag .gt. 539)))
c$$$     $           then

c     Isolate new and old coordinates of molecule i.
               do j = 1, 3
                  index1 = 6*nImages*(i-1) + 3*(image-1) + j
                  index2 = index1 + 3*nImages
                  indexO = 3*nImages*(i-1) + 3*nHydrogens + 
     $                 3*(image-1) + j
                  r1O(j) = RAold(index1)
                  r1N(j) = RA(index1)
                  r2O(j) = RAold(index2)
                  r2N(j) = RA(index2)
                  r3O(j) = RAold(indexO)
                  r3N(j) = RA(indexO)
                  call reconstruct(r1O, r2O, r3O)
                  call reconstruct(r1N, r2N, r3N)
               end do

               call restoreMolecule(r1O, r2O, r3O, r1N, r2N, r3N)

C Debug
C       write(10,FMT='(1X,A)') ' '
C       write(10,FMT='(1X,3E23.15)') r1O
C       write(10,FMT='(1X,3E23.15)') r2O
C       write(10,FMT='(1X,3E23.15)') r3O
C       write(10,FMT='(1X,3E23.15)') r1N
C       write(10,FMT='(1X,3E23.15)') r2N
C       write(10,FMT='(1X,3E23.15)') r3N
C       write(10,FMT='(1X,A)') ' '

               m1 = sscrh(1)
               m2 = sscrh(1)
               itagOx = itag(nHydrogens+nImages*(i-1)+image)
               m3 = invMass(itagOx, sscrh(2))

C Debug
C       write(10,FMT='(1X,A)') ' '
C       write(10,FMT='(1X,3F20.10)') m1, m2, m3
C       write(10,FMT='(1X,A)') ' '

c-----------------------------------------------------------------------
c$$$            print '(A, 2i5)', 'Before shaker1', i, image
c$$$            print '(3i8)', nMolecules, nHydrogens, nOxygens
c$$$            print '(3i8)', index1, index2, indexO
c$$$            print '(A,3g17.8)', 'r1 : ',(r1N(j), j=1,3)
c$$$            print '(A,3g17.8)', 'r2 : ',(r2N(j), j=1,3)
c$$$            print '(A,3g17.8)', 'r3 : ',(r3N(j), j=1,3)
c$$$            print *
c$$$            print '(A,3g17.8)', 'r1 : ',(r1O(j), j=1,3)
c$$$            print '(A,3g17.8)', 'r2 : ',(r2O(j), j=1,3)
c$$$            print '(A,3g17.8)', 'r3 : ',(r3O(j), j=1,3)
c$$$            pause 'Shaker1 : Stopped at SHAKER 1 before SHAKE'
c-----------------------------------------------------------------------
                       
c           Make corrections in molecule i.
               call shake(r1N, r2N, r3N, r1O, r2O, r3O, cf1, cf2, cf3, 
     $              m1, m2, m3, dt, i)
c                   ^^^^^ sscrh contains dt/(2 m) for each type of atom

C Debug
C       write(10,FMT='(1X,A)') ' cf '
C       write(10,FMT='(1X,3E23.15)') cf1, cf2, cf3
C       write(10,FMT='(1X,A)') ' '

c$$$            print '(A, 2i5)', 'After shaker1 (mole, image)', i, image
c$$$            print '(3i8)', nMolecules, nHydrogens, nOxygens
c$$$            print '(3i8)', index1, index2, indexO
c$$$            print '(A,3g17.8)', 'r1 : ',(r1N(j), j=1,3)
c$$$            print '(A,3g17.8)', 'r2 : ',(r2N(j), j=1,3)
c$$$            print '(A,3g17.8)', 'r3 : ',(r3N(j), j=1,3)
c$$$            print *
c$$$            print '(A,3g17.8)', 'r1 : ',(r1O(j), j=1,3)
c$$$            print '(A,3g17.8)', 'r2 : ',(r2O(j), j=1,3)
c$$$            print '(A,3g17.8)', 'r3 : ',(r3O(j), j=1,3)

c           Re-store the new coordinates in array RA and VA
c           cf are the constraint forces times dt/2m
               do j = 1, 3
                  index1 = 6*nImages*(i-1) + 3*(image-1) + j
                  index2 = index1 + 3*nImages
                  indexO = 3*nHydrogens + 3*nImages*(i-1) +
     $                 3*(image-1) + j
                  RA(index1) = RA(index1) + cf1(j) * dt
                  RA(index2) = RA(index2) + cf2(j) * dt
                  RA(indexO) = RA(indexO) + cf3(j) * dt

                  VA(index1) = VA(index1) + cf1(j)
                  VA(index2) = VA(index2) + cf2(j)
                  VA(indexO) = VA(indexO) + cf3(j)

                  FA(index1) = FA(index1) + cf1(j) / sscrh(1)
                  FA(index2) = FA(index2) + cf2(j) / sscrh(1)
                  FA(indexO) = FA(indexO) + cf3(j) / sscrh(2)
               end do

c$$$            call calcDiffe(r1O, r2O, r3O, r12, r23, r31)
c$$$            print '(3f15.8)', sqrt(dot(r12,r12)), 
c$$$     $           sqrt(dot(r23,r23)), sqrt(dot(r31,r31))
c$$$            pause 'Stopped in Shaker1'
c$$$            end if
         end do
      end do

      return
      end
c------------------------------------------------------------------------
      function invMass(i, m)
      real*8 invMass, m

      if (i .eq. 1) then
         invMass = 0.0d0
      else
         invMass = m
      end if
      return
      end

c------------------------------------------------------------------------
      subroutine shake(r1N, r2N, r3N, r1O, r2O, r3O, cf1, cf2, cf3, m1, 
     $     m2, m3, dt, moleNumber)

      implicit real*8 (a-h,o-z)

      include '../commonblks/compotent.cmn'
      include '../commonblks/comgeom.cmn'    

      real*8 m(3,3), r1N(3), r2N(3), r3N(3), r12N(3), r23N(3), r31N(3) 
      real*8 r1O(3), r2O(3), r3O(3), r12O(3), r23O(3), r31O(3), d(3) 
      real*8 c(3,3), d1(3), d2(3), d3(3), q12(3), q23(3), q31(3)
      real*8 qSq, drSq, b(3), lambda(3), dt, dt2, dots(3,3), pi, g
      real*8 dummy, dot, epsilon, cf1(3), cf2(3), cf3(3), m1, m2, m3
      integer indx(3), N

C Debug
C     write(10,FMT='(1X,A)') ' Entering shake'

C Debug
C     write(10,FMT='(1X,A)') ' r1, r2, r3 '
C     write(10,FMT='(1X,3E23.15)') (r1N(p),p=1,3)
C     write(10,FMT='(1X,3E23.15)') (r2N(p),p=1,3)
C     write(10,FMT='(1X,3E23.15)') (r3N(p),p=1,3)
C     write(10,FMT='(1X,3E23.15)') (r1O(p),p=1,3)
C     write(10,FMT='(1X,3E23.15)') (r2O(p),p=1,3)
C     write(10,FMT='(1X,3E23.15)') (r3O(p),p=1,3)
C     write(10,FMT='(1X,A)') ' '

C Debug
C     write(10,FMT='(1X,A)') ' cf1, cf2, cf3 '
C     write(10,FMT='(1X,3E23.15)') (cf1(p),p=1,3)
C     write(10,FMT='(1X,3E23.15)') (cf2(p),p=1,3)
C     write(10,FMT='(1X,3E23.15)') (cf3(p),p=1,3)
C     write(10,FMT='(1X,A)') ' '

C Debug
C     write(10,FMT='(1X,A)') ' m1, m2, m3, '
C     write(10,FMT='(1X,3E23.15)') m1, m2, m3
C     write(10,FMT='(1X,A)') ' '

C Debug
C     write(10,FMT='(1X,A)') ' dt, moleNumber '
C     write(10,FMT='(1X,E23.15,I5)') dt, moleNumber
C     write(10,FMT='(1X,A)') ' '

      call calcDiffe(r1N, r2N, r3N, r12N, r23N, r31N)
      call calcDiffe(r1O, r2O, r3O, r12O, r23O, r31O)

C Debug
C     write(10,FMT='(1X,A)') ' r12N, r23N, r31N '
C     write(10,FMT='(1X,3E23.15)') (r12N(p),p=1,3)
C     write(10,FMT='(1X,3E23.15)') (r23N(p),p=1,3)
C     write(10,FMT='(1X,3E23.15)') (r31N(p),p=1,3)
C     write(10,FMT='(1X,A)') ' r12O, r23O, r31O '
C     write(10,FMT='(1X,3E23.15)') (r12O(p),p=1,3)
C     write(10,FMT='(1X,3E23.15)') (r23O(p),p=1,3)
C     write(10,FMT='(1X,3E23.15)') (r31O(p),p=1,3)
C     write(10,FMT='(1X,A)') ' '

c      print '(A, i5)', 'Shaking1 molecule ', moleNumber
      pi = 4.d0 * atan(1.d0)
      Ro = potpar(14, 1)
      theta0 = potpar(15, 1)*pi/180.d0

c      d(2) = Ro*Ro
      d(2) = Ro*Ro/(ax*ax)
      d(3) = d(2)
      d(1) = 4.d0 * d(2) * sin(theta0/2.d0)**2

C Debug
C     write(10,FMT='(1X,A)') ' d '
C     write(10,FMT='(1X,3E23.15)') d
C     write(10,FMT='(1X,A)') ' '

c     *** Note that m1, m2 and m3 contain dt/(2 mi) that is the inverse of
c     twice the mass times dt ***
      c(1,1) = (m1+m2) * dt 
      c(1,2) = -m2 * dt
      c(1,3) = -m1 * dt
      c(2,1) = c(1,2)
      c(2,2) = (m2+m3) * dt
      c(2,3) = -m3 * dt
      c(3,1) = c(1,3)
      c(3,2) = c(2,3)
      c(3,3) = (m3+m1) * dt
      
C Debug
C     write(10,FMT='(1X,A)') ' c '
C     write(10,FMT='(1X,3E23.15)') (c(1,p),p=1,3)
C     write(10,FMT='(1X,3E23.15)') (c(2,p),p=1,3)
C     write(10,FMT='(1X,3E23.15)') (c(3,p),p=1,3)
C     write(10,FMT='(1X,A)') ' '

      call calcDots(dots, r12O, r23O, r31O, r12N, r23N, r31N)

      call createM(c, dots, m)

c      print '(A)', 'Distances'
c      print '(3g17.8)', (d(i), i=1,3)
c      print *

c      print '(A)', 'dt/2mi'
c      print '(3g17.8)', m1, m2, m3
c      print *

c$$$      print '(A)', 'Mass matrix'
c$$$      do i = 1, 3
c$$$         print '(3g17.8)', (c(i,j), j=1,3)
c$$$      end do
c$$$      print *
c$$$
c$$$      print '(A)', 'Dots matrix'
c$$$      do i = 1, 3
c$$$         print '(3g17.8)', (dots(i,j), j=1,3)
c$$$      end do
c$$$      print *
c$$$
c$$$      print '(A)', 'M matrix'
c$$$      print *, 'Molecule number ', moleNumber
c$$$      do i = 1, 3
c$$$         print '(3g17.8)', (m(i,j), j=1,3)
c$$$      end do
c$$$      print *
c$$$         
c$$$      print *, 'Shaker1 :  Stopped in SHAKE before ludcmp'
c$$$      pause 'Shaker1 :  Stopped in SHAKE before ludcmp'

C Debug
C     write(10,FMT='(1X,A)') ' m before ludcmp '
C     write(10,FMT='(1X,3E23.15)') (m(1,p),p=1,3)
C     write(10,FMT='(1X,3E23.15)') (m(2,p),p=1,3)
C     write(10,FMT='(1X,3E23.15)') (m(3,p),p=1,3)
C     write(10,FMT='(1X,A)') ' '

      call ludcmp(m, 3, 3, indx, dummy)

C Debug
C     write(10,FMT='(1X,A)') ' m after ludcmp '
C     write(10,FMT='(1X,3E23.15)') (m(1,p),p=1,3)
C     write(10,FMT='(1X,3E23.15)') (m(2,p),p=1,3)
C     write(10,FMT='(1X,3E23.15)') (m(3,p),p=1,3)
C     write(10,FMT='(1X,A)') ' '

C Debug
C     write(10,FMT='(1X,A)') ' indx '
C     write(10,FMT='(1X,3I5)') indx
C     write(10,FMT='(1X,A)') ' '

c$$$      print '(A)', 'M matrix'
c$$$      print *, 'Molecule number ', moleNumber
c$$$      do i = 1, 3
c$$$         print '(3g17.8)', (m(i,j), j=1,3)
c$$$      end do
c$$$      print *
c$$$
c$$$      pause 'Shaker1 :  Stopped in SHAKE right after ludcmp'

      lambda(1) = 0.d0
      lambda(2) = 0.d0
      lambda(3) = 0.d0

      epsilon = 1.d0
      iteration = 0
      do while(epsilon .gt. 1e-5)
         iteration = iteration + 1

c$$$         print '(A,i5,A,f15.10)', 'Iteration: ', iteration, 
c$$$     $        'epsilon: ', epsilon
         if (iteration .gt. 10) then
            print '(A, i5, A)', 'Shake did not converge in molecule ', 
     $           moleNumber,'. Hit any key to do 10'
            print '(A)', 'more iterations or Ctrl+C to stop'

            print '(3f15.10)', (r12N(ii), ii=1,3)
            print '(3f15.10)', (r31N(ii), ii=1,3)
            print '(3f15.10)', (r23N(ii), ii=1,3)

            stop
c$$$            pause 'shaking ...'
            iteration = 0
         end if

         do j = 1, 3
            d1(j) = c(1,j) * lambda(j)
            d2(j) = c(2,j) * lambda(j)
            d3(j) = c(3,j) * lambda(j)
         end do
         
C Debug
C     write(10,FMT='(1X,A)') ' iteration d '
C     write(10,FMT='(1X,I5,3E23.15)') iteration, d1
C     write(10,FMT='(1X,I5,3E23.15)') iteration, d2
C     write(10,FMT='(1X,I5,3E23.15)') iteration, d3
C     write(10,FMT='(1X,A)') ' '

         do i = 1, 3
            q12(i) = r12N(i) + (d1(1) * r12O(i) + d1(2) * r23O(i) + 
     $           d1(3) * r31O(i))            
            q23(i) = r23N(i) + (d2(1) * r12O(i) + d2(2) * r23O(i) + 
     $           d2(3) * r31O(i))
            q31(i) = r31N(i) + (d3(1) * r12O(i) + d3(2) * r23O(i) + 
     $           d3(3) * r31O(i))
         end do

C Debug
C     write(10,FMT='(1X,A)') ' iteration q '
C     write(10,FMT='(1X,I5,3E23.15)') iteration, q12
C     write(10,FMT='(1X,I5,3E23.15)') iteration, q23
C     write(10,FMT='(1X,I5,3E23.15)') iteration, q31
C     write(10,FMT='(1X,A)') ' '

         epsilon = 0.d0
         qSq = dot(q12, q12)
         drSq = 2.d0 * (d1(1) * dots(1,1) + d1(2) * dots(1,2) + 
     $        d1(3) * dots(1,3))
         b(1) = d(1) - qSq 
         epsilon = epsilon + abs(b(1))/d(1)
         b(1) = b(1) + drSq
c         print '(A,g15.7, $)', 'q^2 : ', qSq
         
         qSq = dot(q23, q23)
         drSq = 2.d0 * (d2(1) * dots(2,1) + d2(2) * dots(2,2) + 
     $        d2(3) * dots(2,3))  
         b(2) = d(2) - qSq 
         epsilon = epsilon + abs(b(2))/d(2)
         b(2) = b(2) + drSq
c         print '(g15.7, $)', qSq

         qSq = dot(q31, q31)
         drSq = 2.d0 * (d3(1) * dots(3,1) + d3(2) * dots(3,2) + 
     $        d3(3) * dots(3,3))  
         b(3) = d(3) - qSq 
         epsilon = epsilon + abs(b(3))/d(3)
         b(3) = b(3) + drSq
c         print '(g15.7,$)', qSq

C Debug
C     write(10,FMT='(1X,A)') ' m before lubksb '
C     write(10,FMT='(1X,3E23.15)') (m(1,p),p=1,3)
C     write(10,FMT='(1X,3E23.15)') (m(2,p),p=1,3)
C     write(10,FMT='(1X,3E23.15)') (m(3,p),p=1,3)
C     write(10,FMT='(1X,A)') ' '

C Debug
C     write(10,FMT='(1X,A)') ' indx before lubksb '
C     write(10,FMT='(1X,3I5)') indx
C     write(10,FMT='(1X,A)') ' '

C Debug
C     write(10,FMT='(1X,A)') ' iteration b before lubksb'
C     write(10,FMT='(1X,I5,3E23.15)') iteration, b
C     write(10,FMT='(1X,A)') ' '

         call lubksb(m, 3, 3, indx, b)

C Debug
C     write(10,FMT='(1X,A)') ' m after lubksb '
C     write(10,FMT='(1X,3E23.15)') (m(1,p),p=1,3)
C     write(10,FMT='(1X,3E23.15)') (m(2,p),p=1,3)
C     write(10,FMT='(1X,3E23.15)') (m(3,p),p=1,3)
C     write(10,FMT='(1X,A)') ' '

C Debug
C     write(10,FMT='(1X,A)') ' indx after lubksb '
C     write(10,FMT='(1X,3I5)') indx
C     write(10,FMT='(1X,A)') ' '

C Debug
C     write(10,FMT='(1X,A)') ' iteration b after lubksb'
C     write(10,FMT='(1X,I5,3E23.15)') iteration, b
C     write(10,FMT='(1X,A)') ' '

         do i = 1, 3
            lambda(i) = b(i)
            b(i) = 0.d0
         end do
         epsilon = epsilon/3.d0
         
c         iter = iter+1
c         print '(i4,g13.8)', iter, epsilon

      end do

c$$$      print '(A,i3,A,A,2x,g13.8)',' after', iter, ' iterations, ',
c$$$     $     'epsilon =', epsilon
c      iter = 0

c      print '(A)', 'Constraint forces in shaker1:'

C Debug
C     write(10,FMT='(1X,A)') ' lambda '
C     write(10,FMT='(1X,3F20.14)') lambda
C     write(10,FMT='(1X,A)') ' '

      do i = 1, 3
         g = m1 * (lambda(1) * r12O(i) - lambda(3) * r31O(i))
         cf1(i) = g
c         print '(i4,g17.8,$)', i, g 

         g = m2 * (-lambda(1) * r12O(i) + lambda(2) * r23O(i))
         cf2(i) = g
c         print '(g17.8,$)', g

         g = m3 * (-lambda(2) * r23O(i) + lambda(3) * r31O(i))
         cf3(i) = g
c         print '(g17.8,f8.4)', g 
      end do
c      pause ' Check the constraint forces Sir!'

C Debug
C     write(10,FMT='(1X,A)') ' Leaving shake'

      end

c------------------------------------
      subroutine createM(c, dots, m)
      real*8 c(3,3), dots(3,3), m(3,3)

      do i = 1, 3
         do j = 1, 3
            m(i,j) = 2.d0 * c(i,j) * dots(i,j)
         end do
      end do

      return
      end
c------------------------------------
      subroutine calcDots(dots, r12O, r23O, r31O, r12N, r23N, r31N)
      real*8 dots(3,3), r12O(3), r23O(3), r31O(3), dot
      real*8 r12N(3), r23N(3), r31N(3)

      dots(1,1) = dot(r12O, r12N)
      dots(1,2) = dot(r23O, r12N)
      dots(1,3) = dot(r31O, r12N)

      dots(2,1) = dot(r12O, r23N)
      dots(2,2) = dot(r23O, r23N)
      dots(2,3) = dot(r31O, r23N)

      dots(3,1) = dot(r12O, r31N)
      dots(3,2) = dot(r23O, r31N)
      dots(3,3) = dot(r31O, r31N)

      return
      end
c------------------------------------------------------------------------
      subroutine calcDiffe(r1, r2, r3, r12, r23, r31)
      real*8 r1(3), r2(3), r3(3)
      real*8 r12(3), r23(3), r31(3)

      call add(r1, r2, -1, r12)
      call add(r2, r3, -1, r23)
      call add(r3, r1, -1, r31)

      return
      end
c-----------------------------------
      subroutine add(v1, v2, a, v3)
      real*8 v1(3), v2(3), v3(3)
      integer a

      do i = 1, 3
         v3(i) = v1(i) + a * v2(i)
      end do

      return
      end
c------------------------------------------------------------------------

      subroutine restoreMolecule(r1O, r2O, r3O, r1N, r2N, r3N)

      implicit real*8 (a-h,o-z)
      include '../commonblks/parameters.cmn'

      real*8 r1O(3), r2O(3), r3O(3), r1N(3), r2N(3), r3N(3)

      call resetHydrogen(r1O, r3O)
      call resetHydrogen(r2O, r3O)

      call resetHydrogen(r1N, r3N)
      call resetHydrogen(r2N, r3N)

      return
      end


      subroutine resetHydrogen(r, rO)

      implicit real*8 (a-h,o-z)
      include '../commonblks/comgeom.cmn'    

      real*8 r(3), rO(3), dist

      dist = r(1)-rO(1)
      if (dist .gt. 1.d0/2.0d0) then
         r(1) = r(1) - 1.d0
      elseif(dist .lt. -1.d0/2.0d0) then
         r(1) = r(1) + 1.d0
      end if
      
      dist = r(2)-rO(2)
      if (dist .gt. ay/2.0d0/ax) then
         r(2) = r(2) - ay/ax
      elseif(dist .lt. -ay/2.0d0/ax) then
         r(2) = r(2) + ay/ax
      end if
      
      dist = r(3)-rO(3)
      if (dist .gt. az/2.0d0/ax) then
         r(3) = r(3) - az/ax
      elseif(dist .lt. -az/2.0d0) then
         r(3) = r(3) + az/ax
      end if

      return
      end
