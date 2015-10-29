c-----------------------------------------------------------------------
      subroutine readPoles(dpole, qpole, opole, hpole)

C Modified by Fer
      implicit none

      integer   IBUFFER
      parameter (IBUFFER = 256)
      character(LEN=IBUFFER) line
      integer*4 iUnit
      integer*4 ios

      integer*4 i
      integer*4 p, q, r, s
      integer*4 nC
      real*8    kk1, kk2, Val
      real*8 dpole(3), qpole(3,3), opole(3,3,3), hpole(3,3,3,3)

C Initialize the multipole moments
      do p=1,3
        dpole(p) = 0.0D0
        do q=1,3
          qpole(p,q) = 0.0D0
          do r=1,3
            opole(p,q,r) = 0.0D0
            do s=1,3
              hpole(p,q,r,s) = 0.0D0
            end do
          end do
        end do
      end do

C Open the multipoles file
      iUnit = 55
      open(iUnit, file="multipoles", status="old")
      ios = 0

C au to Debye constant
      kk1 = 2.5417709D0

C Ang to au constant
      kk2 = 1.88972666351031921149

C Read the dipole moment components
      call getLine(iUnit, line, ios)
      read(line, *) nC

      do i=1,nC
        call getLine(iUnit, line, ios)
        read(line, *) p, Val
        dpole(p) = Val*kk1
      end do

C Read the quadrupole moment components
      call getLine(iUnit, line, ios)
      read(line, *) nC

      do i=1,nC
        call getLine(iUnit, line, ios)
        read(line, *) p, q, Val
        qpole(p,q) = Val*kk1/kk2
      end do

C Read the octupole moment components
      call getLine(iUnit, line, ios)
      read(line, *) nC

      do i=1,nC
        call getLine(iUnit, line, ios)
        read(line, *) p, q, r, Val
        opole(p,q,r) = Val*kk1/(kk2)**2
      end do

C Read the hexadecapole moment components
      call getLine(iUnit, line, ios)
      read(line, *) nC

      do i=1,nC
        call getLine(iUnit, line, ios)
        read(line, *) p, q, r, s, Val
        hpole(p,q,r,s) = Val*kk1/(kk2)**3
      end do

C Enforce the symmetry constraints
C Note: Symmetry constraints are only enforced for non-zero
C       components in water

      opole(1,3,1) = opole(1,1,3)
      opole(3,1,1) = opole(1,1,3)
      opole(2,3,2) = opole(2,2,3)
      opole(3,2,2) = opole(2,2,3)

      hpole(1,2,1,2) = hpole(1,1,2,2)
      hpole(1,2,2,1) = hpole(1,1,2,2)
      hpole(2,1,1,2) = hpole(1,1,2,2)
      hpole(2,1,2,1) = hpole(1,1,2,2)
      hpole(2,2,1,1) = hpole(1,1,2,2)
      hpole(1,3,1,3) = hpole(1,1,3,3)
      hpole(1,3,3,1) = hpole(1,1,3,3)
      hpole(3,1,1,3) = hpole(1,1,3,3)
      hpole(3,1,3,1) = hpole(1,1,3,3)
      hpole(3,3,1,1) = hpole(1,1,3,3)
      hpole(2,3,2,3) = hpole(2,2,3,3)
      hpole(2,3,3,2) = hpole(2,2,3,3)
      hpole(3,2,2,3) = hpole(2,2,3,3)
      hpole(3,2,3,2) = hpole(2,2,3,3)
      hpole(3,3,2,2) = hpole(2,2,3,3)
        
C Close the multipoles file
      close(iUnit)

      end

c-----------------------------------------------------------------------
c     This routine extracts the next uncommented line from a file. The
c     comments are specified by starting the line with a '#' sign

      subroutine getLine(unit, line, ios)

      parameter (IBUFFER = 256)
      integer unit, ios
c      character(LEN=IBUFFER) line, fmt
      character*256 line, fmt
      character*1 first

c      write(fmt, 100) IBUFFER
c 100  format('(A',i3.3,')')

      write(fmt, '("(A",i3.3,")")') IBUFFER

      read(unit, fmt, IOSTAT=ios, end=200, err=200) line
      read(line, '(A1)') first
      do while (first .eq. '#')
         read(unit, fmt, IOSTAT=ios, end=200, err=200) line
         read(line, '(A1)') first
      end do

 200  return
      end
