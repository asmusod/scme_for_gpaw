c                                                                                         
      FUNCTION ERFINV(R)
c                                                                                         
      implicit real*8 (a-h,o-z)
c                                                                                         
      DIMENSION D(5)
c                                                                                         
      DATA RT1,RT2/0.707107,1.224745/
      DATA RTPI2,D(1)/0.8862269,1.0/
c                                                                                         
      ERFINV=SIGN(1.0D+00,R)
      R=ABS(R)
      IF(R.GT..5) GO TO 10
        XG=R*RTPI2
        GO TO 20
10    XG=3.0189*R**2-2.5241*R+.9833
c                                                                                         
20    RG=ERF(XG)
      DELR=R-RG
      IF(ABS(DELR).LT..0001) GO TO 40
      IF(XG-0.0001) 33,31,31
31    IF(ABS(XG-RT1)-0.0001) 34,32,32
32    IF(ABS(XG-RT2)-0.0001) 34,35,35
33    XG=R*RTPI2
      GOTO 20
34    X=DELR
      GO TO 39
35    D(2)=-2.*XG
      D(3)=4.*XG**2-2.
      D(4)=12.*XG-8.*XG**3
      D(5)=16.*XG**4-48.*XG**2+12.
      ADD=1.
      X=0.
      DO 30 I=1,5
        ADD=ADD*DELR/FLOAT(I)
30      X=X+ADD/D(I)
39    XG=XG+RTPI2*EXP(XG**2)*X
      GO TO 20
40    ERFINV=ERFINV*XG
c                                                                                         
      RETURN
      END
                                                                                   

