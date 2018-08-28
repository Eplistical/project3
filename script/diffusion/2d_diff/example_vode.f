      PROGRAM TEST
      EXTERNAL CAL_DUDT, CAL_JAC
      DOUBLE PRECISION ATOL, RPAR, RTOL, RWORK, T, TOUT, U
      DIMENSION U(3), RWORK(67), IWORK(33)
      NEQ = 3
      U(1) = 1.0D0
      U(2) = 0.0D0
      U(3) = 0.0D0
      T = 0.0D0
      TOUT = 0.4D0
      ITOL = 1
      RTOL = 1.D-4
      ATOL = 1.D-14
      ITASK = 1
      ISTATE = 1
      IOPT = 0
      LRW = 67
      LIW = 33
      MF = 21

      DO 40 IOUT = 1,11
        CALL DVODE(CAL_DUDT,NEQ,U,T,TOUT,ITOL,RTOL,ATOL,ITASK,ISTATE,
     1            IOPT,RWORK,LRW,IWORK,LIW,CAL_JAC,MF,RPAR,IPAR)
        WRITE(6,20)T,U(1),U(2),U(3)
  20    FORMAT(' At t =',D12.4,'   y =',3D14.6)
        IF (ISTATE .LT. 0) GO TO 80
  40    TOUT = TOUT*10.

      STOP
  80  WRITE(6,90)ISTATE
  90  FORMAT(///' Error halt: ISTATE =',I3)
      END


      SUBROUTINE CAL_DUDT (NEQ, T, U, DUDT, RPAR, IPAR)

      DOUBLE PRECISION RPAR, T, U, DUDT
      DIMENSION U(NEQ), DUDT(NEQ)

      DUDT(1) = -.04D0*U(1) + 1.D4*U(2)*U(3)
      DUDT(3) = 3.D7*U(2)*U(2)
      DUDT(2) = -DUDT(1) - DUDT(3)

      END


      SUBROUTINE CAL_JAC (NEQ, T, U, ML, MU, PD, NRPD, RPAR, IPAR)
      DOUBLE PRECISION PD, RPAR, T, U

      DIMENSION U(NEQ), PD(NRPD,NEQ)

      PD(1,1) = -.04D0
      PD(1,2) = 1.D4*U(3)
      PD(1,3) = 1.D4*U(2)
      PD(2,1) = .04D0
      PD(2,3) = -PD(1,3)
      PD(3,2) = 6.D7*U(2)
      PD(2,2) = -PD(1,2) - PD(3,2)

      END
