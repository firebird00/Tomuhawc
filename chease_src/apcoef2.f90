!*DECK C2SP06
!*CALL PROCESS
SUBROUTINE APCOEF2(KN,PX,PAP,PT)
  !        ################################
  !
  !                                        AUTHORS:
  !                                        H. LUTJENS,  CRPP-EPFL
  !                                        A. BONDESON, CRPP-EPFL
  !**********************************************************************
  !                                                                     *
  !  C2SP06  EVALUATE P' IF GIVEN AS POLYNOMIALS IN SEVERAL SECTIONS    *
  !                                                                     *
  !**********************************************************************
  !
  USE globals
  IMPLICIT NONE
  REAL(RKIND)      ::     PT
  REAL(RKIND)      ::     ZS6
  REAL(RKIND)      ::     ZS5
  REAL(RKIND)      ::     ZS4
  REAL(RKIND)      ::     ZS3
  REAL(RKIND)      ::     ZS2
  REAL(RKIND)      ::     ZS1
  REAL(RKIND)      ::     ZSG6
  REAL(RKIND)      ::     ZSG5
  REAL(RKIND)      ::     ZSG4
  REAL(RKIND)      ::     ZSG3
  REAL(RKIND)      ::     ZSG2
  REAL(RKIND)      ::     PX
  REAL(RKIND)      ::     ZSG1
  INTEGER          ::     J1
  REAL(RKIND)      ::     ZF2
  REAL(RKIND)      ::     ZF1
  REAL(RKIND)      ::     ZF0
  REAL(RKIND)      ::     ZE1
  REAL(RKIND)      ::     ZE0
  REAL(RKIND)      ::     ZD3
  REAL(RKIND)      ::     ZD2
  REAL(RKIND)      ::     ZD1
  REAL(RKIND)      ::     ZD0
  REAL(RKIND)      ::     ZC2
  REAL(RKIND)      ::     ZC1
  REAL(RKIND)      ::     ZC0
  REAL(RKIND)      ::     ZB3
  REAL(RKIND)      ::     ZB2
  REAL(RKIND)      ::     ZB1
  REAL(RKIND)      ::     ZB0
  REAL(RKIND)      ::     ZA1
  REAL(RKIND)      ::     ZA0
  REAL(RKIND)      ::     PAP
  INTEGER          ::     KN
  DIMENSION &
       &   PAP(*),  PT(KN),   PX(KN)
  !
  !----*----*----*---*----*----*----*----*----*----*----*----*----*----*-
  !
  CALL COPYAPP(PAP,ZA0,ZA1,ZB0,ZB1,ZB2,ZB3,ZC0,ZC1,ZC2, &
       &                    ZD0,ZD1,ZD2,ZD3,ZE0,ZE1,ZF0,ZF1,ZF2)
  !
  DO J1=1,KN
     !
     ZSG1 = SIGN(RC1P,PAP(1) - PX(J1))
     ZSG2 = SIGN(RC1P,(PX(J1) - PAP(1)) * (PAP(2) - PX(J1)))
     ZSG3 = SIGN(RC1P,(PX(J1) - PAP(2)) * (PAP(3) - PX(J1)))
     ZSG4 = SIGN(RC1P,(PX(J1) - PAP(3)) * (PAP(4) - PX(J1)))
     ZSG5 = SIGN(RC1P,(PX(J1) - PAP(4)) * (PAP(5) - PX(J1)))
     ZSG6 = SIGN(RC1P,PX(J1) - PAP(5))
     !
     ZS1 = MAX(RC0P,ZSG1)
     ZS2 = MAX(RC0P,ZSG2)
     ZS3 = MAX(RC0P,ZSG3)
     ZS4 = MAX(RC0P,ZSG4)
     ZS5 = MAX(RC0P,ZSG5)
     ZS6 = MAX(RC0P,ZSG6)
     !
     PT(J1) = PT(J1) - ZS1 * (ZA0 + ZA1 * PX(J1)) - &
          &                     ZS2 * (ZB0 + PX(J1) * (ZB1 + PX(J1) * &
          &                                     (ZB2 + PX(J1) * ZB3))) - &
          &                     ZS3 * (ZC0 + PX(J1) * (ZC1 + PX(J1) * ZC2)) - &
          &                     ZS4 * (ZD0 + PX(J1) * (ZD1 + PX(J1) * &
          &                                     (ZD2 + PX(J1) * ZD3))) - &
          &                     ZS5 * (ZE0 + ZE1 * PX(J1)) - &
          &                     ZS6 * (ZF0 + PX(J1) * (ZF1 + PX(J1) * ZF2))
     !
  END DO
  !
  RETURN
END SUBROUTINE APCOEF2
