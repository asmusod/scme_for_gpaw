c                                                                                         
C  Evaluate the error function      (routine from Bill Swope)                             
c                                                                                         
      FUNCTION ERF(X)

      implicit real*8 (a-h,o-z)

      DATA P,A1,A2,A3,A4,A5/0.3275911,0.254829592,-0.284496736,
     +                      1.421413741,-1.453152027,1.061405429/

      T=1./(1.+P*X)
      ERF=1.-(A1+T*(A2+T*(A3+T*(A4+T*A5))))*T*EXP(-X**2)

      RETURN
      END
