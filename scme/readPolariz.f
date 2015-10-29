c-----------------------------------------------------------------------
      subroutine readPolariz(dd, dq, hp, qq)

C Modified by Fer
      implicit none

      integer   IBUFFER
      parameter (IBUFFER = 256)
      character(LEN=IBUFFER) line
      integer*4 iUnit
      integer*4 ios

      integer i
      integer*4 p, q, r, s
      integer*4 nC
      real*8    kk1, kk2, Val
      real*8 dd(3,3), dq(3,3,3), hp(3,3,3), qq(3,3,3,3), scale

C Initialize the polarizabilities
      do p=1,3
        do q=1,3
          dd(p,q) = 0.0D0
          do r=1,3
            hp(p,q,r) = 0.0D0
            dq(p,q,r) = 0.0D0
            do s=1,3
              qq(p,q,r,s) = 0.0D0
            end do
          end do
        end do
      end do

C Open the polarizabilities file
      iUnit = 55
      open(iUnit, file="polar", status="old")
      ios = 0

C au to Debye constant
      kk1 = 2.5417709D0

C Ang to au constant
      kk2 = 1.88972666351031921149

C Read the dipole-dipole polarizability components
      call getLine(iUnit, line, ios)
      read(line, *) nC

      do i=1,nC
        call getLine(iUnit, line, ios)
        read(line, *) p, q, Val
        dd(p,q) = Val/(kk2)**3
      end do

C Read the dipole-quadrupole polarizability components
      call getLine(iUnit, line, ios)
      read(line, *) nC

      do i=1,nC
        call getLine(iUnit, line, ios)
        read(line, *) p, q, r, Val
        dq(p,q,r) = Val/(kk2)**4
      end do

C Read the quadrupole-quadrupole polarizability components
      call getLine(iUnit, line, ios)
      read(line, *) nC

      do i=1,nC
        call getLine(iUnit, line, ios)
        read(line, *) p, q, r, s, Val
        qq(p,q,r,s) = Val/(kk2)**5
      end do

C Enforce the symmetry constraints
C Note: Symmetry constraints are only enforced for non-zero
C       components in water

      dq(1,3,1) = dq(1,1,3)
      dq(2,3,2) = dq(2,2,3)

      qq(2,2,1,1) = qq(1,1,2,2)
      qq(3,3,1,1) = qq(1,1,3,3)
      qq(3,3,2,2) = qq(2,2,3,3)
      qq(1,2,2,1) = qq(1,2,1,2)
      qq(2,1,1,2) = qq(1,2,1,2)
      qq(2,1,2,1) = qq(1,2,1,2)
      qq(1,3,3,1) = qq(1,3,1,3)
      qq(3,1,1,3) = qq(1,3,1,3)
      qq(3,1,3,1) = qq(1,3,1,3)
      qq(2,3,3,2) = qq(2,3,2,3)
      qq(3,2,2,3) = qq(2,3,2,3)
      qq(3,2,3,2) = qq(2,3,2,3)

C Close the polarizabilities file
      close(iUnit)

      end

