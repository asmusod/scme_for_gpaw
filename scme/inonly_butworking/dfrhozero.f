c  Just a bluff, for Tersoff and pair potentials, such as LJ.
c  Dblexp:  June 92
c  Generalized to arbitrary order   Feb 92
c  This routine evaluates teh derivative of the embedding function with respect
c  to its argument (the density).

	SUBROUTINE DFRHODRHO(rho1,natm1,ipot,fprho)
	
	implicit real*8 (a-h,o-z)
	
        include '../commonblks/parameters.cmn'
        include '../commonblks/comconf.cmn'
        include '../commonblks/compotent.cmn'

        dimension fprho(MAXATOMS,3),natm1(3),rho1(MAXATOMS,3)

c     ---------------------------------------------------------------

	do i = 1,natm1(ipot)
	  fprho(i,ipot) = 0.0
	end do
		 	
	RETURN
	END

