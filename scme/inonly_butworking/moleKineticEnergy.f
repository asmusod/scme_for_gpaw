      subroutine moleKineticEnergy(v1x, v1y, v1z, v2x, v2y, v2z, vOx, 
     $     vOy, vOz, Ktot, Krot, Ktra)

c     Variables used:
c     v1x  : x component of the velocity of hydrogen number 1
c     v1sq : v1 square = v1.v1
c     v12  : v1.v2
c     VcmX : x component of the velocity of the center of mass
c     Mt   : Total mass of the molecule = 1+1+16 = 18

      real*8 v1x, v1y, v1z, v2x, v2y, v2z, vOx, vOy, vOz 
      real*8 Ktot, Krot, Ktra, Ktra2
      real*8 v1sq, v2sq, vOsq, v12, vO1, vO2
      real*8 VcmX, VcmY, VcmZ, vCmSq

c     Total kinetic energy of the molecule: Sum[1/2 Mi Vi.Vi]
      v1sq = v1x*v1x + v1y*v1y + v1z*v1z
      v2sq = v2x*v2x + v2y*v2y + v2z*v2z
      vOsq = vOx*vOx + vOy*vOy + vOz*vOz

      Ktot = (v1sq + v2sq + 16.0d0 * vOsq) / 2.0d0

c     Translational kinetic energy: 1/2 Mt Vcm.Vcm
      VcmX = (v1x + v2x + 16.0d0 * vOx) / 18.0d0
      VcmY = (v1y + v2y + 16.0d0 * vOy) / 18.0d0
      VcmZ = (v1z + v2z + 16.0d0 * vOz) / 18.0d0
      vCmSq = VcmX*VcmX + VcmY*VcmY + VcmZ*VcmZ

      Ktra = 18.0d0 * vCmSq / 2.0d0

c     Rotational kinetic energy
      Krot = Ktot - Ktra

cERB      print '(4f15.10)', Ktra, Krot, Ktot, Ktra2

      return
      end
