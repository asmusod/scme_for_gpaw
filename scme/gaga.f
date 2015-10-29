c-----------------------------------------------------------------------
c                  New potential based on multipole moments.
c     Each molecule is represented as a multipole expansion up to
c     hexadecapole moment.
c-----------------------------------------------------------------------
c     An arrays with positions enters as argument ra(). 
c     It is assumed that the array of positions ra() has a length equal
c     to 9 times the number of molecules and that its structure is the
c     following: 
c           ra(l+6*(i-1)) stores the l-th coordinate of the first
c                         hydrogen in the i-th molecule
c           ra(l+3+6*(i-1)) stores the l-th coordinate of the second
c                         hydrogen in the i-th molecule
c           ra(l+3*(i-1+nHydrogens)) stores the l-th coordinate of the
c                         oxygen in the i-th molecule. (nHydrogens is
c                         the total number of hydrogen atoms)
c     The routine will return an array fa() with the forces on each atom.
c     The array fa() has the same structure as ra().
c     The total potential energy of the configuration is also calculated
c     and returned in 'uTot'.
c
c     This routine was written by Enrique Batista.
c-----------------------------------------------------------------------

c23456789012345678901234567890123456789012345678901234567890123456789012
c        10        20        30        40        50        60        70

      subroutine gagafe(nAtms, raOri, itagl, fa, uTot, 
     $                virial,ETOUT,eQM,natm)

c      implicit real*8 (a-h,o-z)


      implicit none
      integer natm
      include '../commonblks/parameters.cmn'
      include '../commonblks/compotent.cmn'
      include '../commonblks/comtime.cmn'
      include '../commonblks/comgeom.cmn'    

ctimming
      integer*8 ti, tf, irtc
      real*8 t1, t2, t3, t4, t5, t6, t7, rMax, rMax2
ctimming
      real*8 d1, d2

      real*8 pi, raOri(maxCoo), a(3), a2(3), fa(maxCoo)
      real*8 fCM(3,maxCoo/3), fsf(3,maxCoo/3), tau(3,maxCoo/3)
      real*8 uD, uQ, uH, uPol, uDisp, uRho, uCore, uES

c     atomic position, centers of mass, principal axes.
      real*8 ra(maxCoo), rCM(3,maxCoo/3), x(3,3,maxCoo/3)

c     Electric fields
      real*8 eD(3,maxCoo/3), eQ(3,maxCoo/3), eH(3,maxCoo/3) 
      real*8 eT(3,maxCoo/3)
      real*8 eTQM(3,maxCoo/3)
c     QM field from gpaw edens in from ASE
      real*8 eQM
      intent(in) eQM
                
c     Derivatives of E
      real*8 dEddr(3,3,maxCoo/3), dEqdr(3,3,maxCoo/3)
      real*8 dEhdr(3,3,maxCoo/3), dEtdr(3,3,maxCoo/3)

c     High order derivatives of the potential
      real*8 d1v(3,maxCoo/3), d2v(3,3,maxCoo/3), d3v(3,3,3,maxCoo/3)
      real*8 d4v(3,3,3,3,maxCoo/3), d5v(3,3,3,3,3,maxCoo/3) 

      real*8 uTot, faB(maxCoo), u, virial, convFactor
      real*8 uCov, uAng, Amp, rMin, lambda, gam, rOHmax
      integer nM, nO, nH, nAtms(maxComp), i, ii, j, k, l, NC, nst, Np

c     Unpolarized multipoles
      real*8 d0(3), q0(3,3), o0(3,3,3), h0(3,3,3,3)

c     Work multipoles. They start unpolarized and with the induction
c     loop we induce dipoles and quadrupoles.
      real*8 dpole(3,maxCoo/3), qpole(3,3,maxCoo/3)
      real*8 opole(3,3,3,maxCoo/3), hpole(3,3,3,3,maxCoo/3)

      real*8 dpole0(3,maxCoo/3), qpole0(3,3,maxCoo/3)

c     Polarizabilities
      real*8 dd0(3,3), dq0(3,3,3), hp0(3,3,3), qq0(3,3,3,3)
      real*8 dd(3,3,maxCoo/3), dq(3,3,3,maxCoo/3), hp(3,3,3,maxCoo/3)
      real*8 qq(3,3,3,3,maxCoo/3)

c     Move along tetrahedral axis
      integer oi
      real*8 m(3,3), mI(3,3), ft(3), dft(3), df(3), ftddr, dot, ddd

c     Study bulk modulus
      real*8 raOld(maxCoo), scale, dr(3)
      integer iO, iH1, iH2
c     Check derivatives
      integer at, coord, iMol, iComp
      real*8 r1(3), r2(3), ro(3)
      real*8 dx, di(3), fexa, fnum, uOld, r0(3)

      integer itagl(maxatoms)
      logical*1 firstVisit, converged, debug, iSlab
      logical*1 addCore

C     Out to ASE
      real*8 ETOUT(3,nAtms(2))
      integer ct
      intent(out) ETOUT


C     Debug
      integer     p, q, r, s
      real*8      kk1, kk2
      character*3 Axis
      real*8 weirdcheck(nAtms(1)+nAtms(2),3) 

      data firstVisit / .true. /

      save

      debug = .false.
      if (firstVisit) then
         pi = 3.14159265358979312d0


c     Size of the simulation cell
         a(1) = ax
         a(2) = ay
         a(3) = az
         a2(1) = ax/2
         a2(2) = ay/2
         a2(3) = az/2

         nO = nAtms(2)          ! Number of OXYGENS 
         nH = nAtms(1)          ! Number of HYDROGENS  
         nM = nO                ! Number of molecules
         call readPoles(d0, q0, o0, h0)
         call readPolariz(dd0, dq0, hp0, qq0)

         NC = potpar(2,1)       ! Number of cell in each direction

         rMax = 11.d0
         rMax2 = rMax*rMax
c     Convert from Debye^2/angstrom^3 to eV
         convFactor = 14.39975841d0 / 4.803206799d0**2
         
         iSlab = .FALSE.
c         iSlab = .TRUE.

         addCore = .TRUE.
c         addCore = .FALSE.

      end if
      
C Debug
C     do i=1,3*nM
C       write(10,FMT='(1X,3F20.10)') (RaOri(3*(i-1)+p),p=1,3)
C     end do

      uTot = 0.0d0

c$$$      ti = irtc()

c     Recover broken molecules due to PBC
      call recoverMolecules(raOri, ra, nH, nO, a, a2)

c     Aaand then set them back since nothanks
c      ra = raOri        
c      do i=1,(nAtms(1)+nAtms(2))
c         do j=1,3
c            weirdcheck(i,j) ra(i*3)-raOri()
c         end do
c      end do

c***********************
c     Shift one molecule to check F = -grad U

c$$$      iMol = 1
c$$$      iComp = 3
c$$$      iH1 = 2*(iMol -1)
c$$$      iH2 = iH1 + 1
c$$$      iO = nAtms(1) + iMol - 1
c$$$      dx = 0.00025
c$$$      do i = 1, 3
c$$$         r1(i) = ra(3*iH1 + i)
c$$$         r2(i) = ra(3*iH2 + i)
c$$$         rO(i) = ra(3*iO + i)
c$$$      end do
c$$$
c$$$      Np = 160
c$$$
c$$$      open(61, file='forces', status='unknown')
c$$$      do ii = -Np/2, Np/2
c$$$         uTot = 0.d0
c$$$         do i = 1, nAtms(1)+nAtms(2)
c$$$            do j = 1, 3
c$$$               fa(3*(i-1)+j) = 0.d0
c$$$            end do
c$$$         end do
c$$$         ra(3*iH1 + iComp) = r1(iComp) + ii * dx
c$$$         ra(3*iH2 + iComp) = r2(iComp) + ii * dx
c$$$         ra(3*iO + iComp) = rO(iComp) + ii * dx
c************************         


c$$$      goto 113

c$$$      print *, 'In gagafe'

      call calcCentersOfMass(ra, nM, rCM)

      call findPpalAxes(ra, nM, x)

      call rotatePoles(d0, q0, o0, h0, dpole0, qpole0, opole, hpole, nM
     $     ,x) 
      call setUnpolPoles(dpole, qpole, dpole0, qpole0, nM) 
      call rotatePolariz(dd0, dq0, qq0, hp0, dd, dq, qq, hp, nM, x)

C Added by Fer.
      if ( firstVisit ) then

        Axis = 'xyz'
        kk1 = 2.5417709D0
        kk2 = 1.88972666351031921149D0

C Print out the multipole moments
C       write(10,FMT='(1X,A)') '--------------------------------'

C       write(10,FMT='(1X,A)') 'Dipole Moment'
C       do p=1,3
C         write(10,FMT='(1X,3X,A,2F11.5)')
C    &          Axis(p:p),
C    &          d0(p), d0(p)/kk1
C       end do

C       write(10,FMT='(1X,A)') 'Quadrupole Moment'
C       do p=1,3
C         do q=p,3
C           write(10,FMT='(1X,2X,2A,2F11.5)')
C    &          Axis(p:p), Axis(q:q),
C    &          q0(p,q), q0(p,q)/kk1*kk2
C         end do
C       end do

C       write(10,FMT='(1X,A)') 'Octupole Moment'
C       do p=1,3
C         do q=p,3
C           do r=q,3
C             write(10,FMT='(1X,1X,3A,2F11.5)')
C    &              Axis(p:p), Axis(q:q), Axis(r:r),
C    &              o0(p,q,r), o0(p,q,r)/kk1*(kk2)**2
C           end do
C         end do
C       end do

C       write(10,FMT='(1X,A)') 'Hexadecapole Moment'
C       do p=1,3
C         do q=p,3
C           do r=q,3
C             do s=r,3
C               write(10,FMT='(1X,4A,2F11.5)')
C    &                Axis(p:p), Axis(q:q), Axis(r:r), Axis(s:s),
C    &                h0(p,q,r,s), h0(p,q,r,s)/kk1*(kk2)**3
C             end do
C           end do
C         end do
C       end do

C Print out the polarizabilities
C       write(10,FMT='(1X,A)') '--------------------------------'

C       write(10,FMT='(1X,A)') 'DD Polarizability'
C       do p=1,3
C         do q=p,3
C           write(10,FMT='(1X,3X,2A,2F11.5)')
C    &          Axis(p:p), Axis(q:q),
C    &          dd0(p,q), dd0(p,q)*(kk2)**3
C         end do
C       end do

C       write(10,FMT='(1X,A)') 'DQ Polarizability'
C       do p=1,3
C         do q=1,3
C           do r=q,3
C             write(10,FMT='(1X,1X,A,1X,2A,2F11.5)')
C    &              Axis(p:p), Axis(q:q), Axis(r:r),
C    &              dq0(p,q,r), dq0(p,q,r)*(kk2)**4
C           end do
C         end do
C       end do

C       write(10,FMT='(1X,A)') 'QQ Polarizability'
C       do p=1,3
C         do q=p,3
C           do r=p,3
C             do s=r,3
C               if ( q .le. s ) then
C                 write(10,FMT='(1X,2A,1X,2A,2F11.5)')
C    &                  Axis(p:p), Axis(q:q), Axis(r:r), Axis(s:s),
C    &                  qq0(p,q,r,s), qq0(p,q,r,s)*(kk2)**5
C               end if
C               if ( (p.eq.1) .and. (q.eq.3) .and.
C    &               (r.eq.2) .and. (s.eq.2)       ) then
C                 write(10,FMT='(1X,2A,1X,2A,2F11.5)')
C    &                  Axis(p:p), Axis(q:q), Axis(r:r), Axis(s:s),
C    &                  qq0(p,q,r,s), qq0(p,q,r,s)*(kk2)**5
C               end if
C             end do
C           end do
C         end do
C       end do

C       write(10,FMT='(1X,A)') '--------------------------------'

      end if

      call calcEhigh(rCM, opole, hpole, nM, NC, a, a2, uH, eH, dEhdr,
     $     rMax2, iSlab)

c     Here's where the induction loop begins
      converged = .false.

c      converged = .true.
c      ii = 0

      do while (.not. converged)
c         ii = ii + 1
         call calcEdip_quad(rCM, dpole, qpole, nM, NC, a, a2, uD, uQ, eD
     $        , dEddr, rMax2, iSlab)

         call addFields(eH, eD, eT, nM)
         call addFieldsQM(eH,eD,eTQM,eQM,nM) 
c        FOR QM/MM embedding via dpole

         call addDfields(dEhdr, dEddr, dEtdr, nM)

c     Induce dipoles and quadrupoles
         converged = .true.
         call induceDipole(dPole, dpole0, eTQM, dEtdr, dd, dq, hp, nM
     $        , converged)
         call induceQpole(qPole, qpole0, eT, dEtdr, dq, qq, nM,
     $        converged) 

      end do                    !End of the induction loop

c    get out eT for ase - 2015 hack get out DIPOLE
      do ct = 1,nAtms(2)
        do i = 1,3
          ETOUT(i,ct) = dPole(i,ct)
        end do
      end do

c      print *,ETOUT

c      d1 = sqrt(dpole(1,1)**2+dpole(2,1)**2+dpole(3,1)**2)
c      print *, d1
c      stop

c$$$      d1 = 0.d0
c$$$      do i = 1, nAtms(2)
c$$$         d1 = d1 + sqrt(dpole(1,i)**2+dpole(2,i)**2+dpole(3,i)**2)
c$$$      end do
c$$$      print *, 'Avg. dipole: ', d1 / nAtms(2)
c$$$      stop


c     With the polarized multipoles, calculate the derivarives of the
c     electrostatic potential, up to 5th order.

      call calcDv(rCM, dpole, qpole, opole, hpole, nM, NC, a, a2, 
     $     d1v, d2v, d3v, d4v, d5v, rMax2, fsf, iSlab)  

c     Compute the force on the center of mass
      call forceCM(dpole, qpole, opole, hpole, d2v, d3v, d4v, d5v, nM,
     $     fsf, fCM)  

c     Compute the torque on the molecule
      call torqueCM(dpole, qpole, opole, hpole, d1v, d2v, d3v, d4v, nM,
     $     tau)  

c     Find 3 forces, one at oxygen and one at each hydrogen such that
c     the total force and total torque agree with the ones calculated
c     for the multipoles.

      call atomicForces(fCM, tau, ra, rCM, nM, fa)

c     Calculate the energy of interaction between the multipole moments
c     and the electric field of the other molecules.
      call calcEnergy(dpole0, qpole0, opole, hpole, d1v, d2v, d3v, d4v, 
     $     nM, uTot)

c      call calcEnergy(dpole, qpole, opole, hpole, d1v, d2v, d3v, d4v, 
c     $     nM, uTot)  
c      uPol = 0.d0
c      call calcEnergyI(dpole, dpole0, qpole, qpole0, opole,
c     $     hpole, d1v, d2v, d3v, d4v, nM, uPol)
c      print *, uTot*convFactor, uPol*convFactor, (uPol+uTot)*Convfactor
c      stop

c      print '(2f15.7,$)', uTot*convFactor, uPol*convFactor

c      uTot = uTot - uPol
    
      virial = 0.0d0
      do i = 1, 3*(nH+nO)
         fa(i) = convFactor * fa(i)
c         virial = virial - ra(i)*fa(i)
      end do
      uTot = uTot * convFactor

      uES = uTot

c 112  call dispersion2(ra, fa, uDisp, nM, a, a2, NC, rMax2, iSlab)
      call dispersion(ra, fa, uDisp, nM, a, a2)
      uTot = uTot + uDisp
      
      if (addCore) then
 111     call coreInt(ra, fa, uCore, nM, a, a2)
         uTot = uTot + uCore
      end if

C Added by Fer. -------------------------------------------------------------
C     write(10,FMT='(1X,A)')
C    &      '--------------------------------------------'
C     write(10,FMT='(1X,A,2F20.10)')
C    &      ' Electrostatic Energy: ', uES,   627.5100402*uES/27.211396
C     write(10,FMT='(1X,A,2F20.10)')
C    &      ' Dispersion Energy:    ', uDisp, 627.5100402*uDisp/27.211396
C     write(10,FMT='(1X,A,2F20.10)')
C    &      ' Repulsion Energy:     ', uCore, 627.5100402*uCore/27.211396
C     write(10,FMT='(1X,A,2F20.10)')
C    &      ' Total Energy:         ', uTot,  627.5100402*uTot/27.211396
C     write(10,FMT='(1X,A)')
C    &      '--------------------------------------------'
C ---------------------------------------------------------------------------

c      print *, uES, uDisp, uCore
c      stop

c****************************
c     Check F = - grad U

c      fexa = fa(3*iH1+iComp) + fa(3*iH2+iComp) + fa(3*iO+iComp)
c      fnum = -(uTot-uOld) / (dx)
c      if (ii .gt. -15) print '(i6,3f15.6)', ii, ii*dx, fexa, fnum
c      uOld = uTot

c      print '(i6,3f15.6)', ii, ii*dx, uTot, ra(3*iO + iComp)

c$$$      if (mod(ii,2) .eq. 0) then
c$$$         fexa = fa(3*iH1+iComp) + fa(3*iH2+iComp) + fa(3*iO+iComp)
c$$$      else
c$$$         if (ii .gt. -Np/2+1) then
c$$$            fnum = -(uTot - uOld)/(2.d0*dx)
c$$$            print '(i6,3f20.12)', ii, ii*dx, fexa, fnum
c$$$            write(61, '(i6,3f20.12)') ii, ii*dx, fexa, fnum
c$$$         end if
c$$$         uOld = uTot
c$$$      end if
c$$$      end do
c$$$      stop
c*******************************


c$$$ 112  open(44, file='allPoles', status='unknown')
c$$$c     Dipole
c$$$      write(44,'(A)') '#'
c$$$      write(44,'(A)') '# **** DIPOLES ****'
c$$$      do i = 1, nAtms(2)
c$$$         d1 = sqrt(dpole(1,i)**2+dpole(2,i)**2+dpole(3,i)**2)
c$$$         write(44,'(i5, f5.2,4f15.10)') i, 8.d0, dpole(1,i), dpole(2,i),
c$$$     $        dpole(3,i), d1
c$$$      end do
c$$$      
c$$$c     quadrupole
c$$$      write(44,'(A)') '#'
c$$$      write(44,'(A)') '# **** QUADRUPOLES ****'
c$$$      do i = 1, nAtms(2)
c$$$         write(44, '(A,i5)') '# Molecule ', i
c$$$         do j = 1, 3
c$$$            write(44,'(3f20.10)') qpole(1,j,i), qpole(2,j,i), qpole(3,j
c$$$     $           ,i)
c$$$         end do
c$$$         write(44, '(A)') '#'
c$$$      end do
c$$$
c$$$c     octopole
c$$$      write(44,'(A)') '# **** OCTOPOLES ****'
c$$$      do i = 1, nAtms(2)
c$$$         write(44, '(A,i5)') '# Molecule ', i
c$$$         do j = 1, 3
c$$$            do k = 1, 3
c$$$               write(44,'(3f20.10)') opole(1,j,k,i)*6.d0, opole(2,j,k,i)
c$$$     $              *6.d0, opole(3,j,k,i)*6.d0
c$$$            end do
c$$$            write(44, '(A)') '#'
c$$$         end do
c$$$      end do
c$$$
c$$$c     hexadecapole
c$$$      write(44,'(A)') '# **** HEXADECAPOLES ****'
c$$$      do i = 1, nAtms(2)
c$$$         write(44, '(A,i5)') '# Molecule ', i
c$$$         do j = 1, 3
c$$$            do k = 1, 3
c$$$               do l = 1, 3
c$$$                  write(44,'(3f20.10)') hpole(1,j,k,l,i)*24.d0, hpole(2
c$$$     $                 ,j,k,l,i)*24.d0, hpole(3,j,k,l,i)*24.d0
c$$$               end do
c$$$               write(44, '(A)') '#'
c$$$            end do
c$$$         end do
c$$$      end do
c$$$      close(44)
c$$$      stop
c$$$
c$$$      if (firstVisit) firstVisit = .false.
c$$$
c$$$      if (.false.) then
c$$$         write(97, '(A)') '=========================================='
c$$$         do i = 1, 288
c$$$            write(97, '(3f16.9,2i6)') fa(3*(i-1)+1), fa(3*(i-1)+2), fa(3
c$$$     $           *(i-1)+3)
c$$$         end do
c$$$      end if

      print '(4f16.10)', uTot, uES, uDisp, uCore
c$$$      stop

C Debug
C     do i=1,nM
C       write(10,FMT='(A,I2,3F16.10)') ' Force: ', i,
C    &  (fa(6*(i-1)+p)+fa(3*(i-1+2*nM)+p)+fa(3+6*(i-1)+p),p=1,3)
C     end do
C     write(10,FMT='(A)') ' '

      return
      end

