!*DECK C2SA06
!*CALL PROCESS
SUBROUTINE PACKMEP(KN,KPOID,PMESH,PPLACE,PWIDTH,PSOLPD)
  !        =======================================================
  !
  !                                        AUTHORS:
  !                                        H. LUTJENS,  CRPP-EPFL
  !                                        A. BONDESON, CRPP-EPFL
  !
  USE globals
  IMPLICIT NONE
  INTEGER          ::     J5
  INTEGER          ::     I
  REAL(RKIND)      ::     ZF
  REAL(RKIND)      ::     ZDP
  REAL(RKIND)      ::     PMESH
  INTEGER          ::     J3
  REAL(RKIND)      ::     PSOLPD
  REAL(RKIND)      ::     ZC
  REAL(RKIND)      ::     PPLACE
  REAL(RKIND)      ::     PWIDTH
  REAL(RKIND)      ::     ZZ
  INTEGER          ::     J1
  REAL(RKIND)      ::     ZW
  REAL(RKIND)      ::     ZS
  INTEGER          ::     J2
  REAL(RKIND)      ::     ZM
  INTEGER          ::     KPOID
  INTEGER          ::     KN
  INTEGER          ::     IM
  PARAMETER (IM = 401)
  !
  DIMENSION &
       &   PMESH(KN),   PPLACE(KPOID),   PWIDTH(KPOID),   ZW(IM)
  !
  !---*----*----*----*----*----*----*----*----*----*----*----*----*----*
  !
  !**********************************************************************
  !                                                                     *
  ! 1. STEP FOR EQUIDISTANT THETA'-MESH                                 *
  !                                                                     *
  !**********************************************************************
  !
  ZM = 2._RKIND * CPI / REAL(IM - 1,RKIND)
  !
  !**********************************************************************
  !                                                                     *
  ! 2. FILL IN DENSITY FUNCTION                                         *
  !                                                                     *
  !**********************************************************************
  !     
  DO J2=1,IM
     !
     ZS     = (J2 - 1) * ZM
     ZW(J2) = 0._RKIND
     !
     DO J1=1,KPOID
        !
        ZZ = SQRT(PWIDTH(J1)**2+1)
        !
        ZW(J2) = ZW(J2) + 2._RKIND/(ZZ*PWIDTH(J1))*( &
             &            ATAN2(PWIDTH(J1)*TAN(.25_RKIND*(ZS-PPLACE(J1))),ZZ+1)+ &
             &            ATAN2(PWIDTH(J1)*TAN(.25_RKIND*(ZS-PPLACE(J1))),ZZ-1)+ &
             &            ATAN2(PWIDTH(J1)*TAN(.25_RKIND*(   PPLACE(J1))),ZZ+1)+ &
             &            ATAN2(PWIDTH(J1)*TAN(.25_RKIND*(   PPLACE(J1))),ZZ-1))
        !
     END DO
  END DO
  !     
  !**********************************************************************
  !                                                                     *
  ! 3. NORMALIZE IT TO ONE                                              *
  !                                                                     *
  !**********************************************************************
  !
  ZC = 2._RKIND* CPI * (1 - PSOLPD) / ZW(IM)
  !
  DO J3=1,IM
     !
     ZS     = (J3 - 1) * ZM
     ZW(J3) = ZS * PSOLPD + ZC * ZW(J3)
     !
  END DO
  !
  !**********************************************************************
  !                                                                     *
  ! 4. FIND MESH POSITIONS                                              *
  !                                                                     *
  !**********************************************************************
  !
  PMESH( 1) = 0._RKIND
  PMESH(KN) = 2._RKIND*CPI
  !
  ZDP = 2._RKIND*CPI / REAL(KN - 1,RKIND)
  ZF  = ZDP
  I   = 1
  !
  DO J5=2,IM
     !
4    CONTINUE 
     !
     IF (ZW(J5) .LE. ZF) GOTO 5
     !
     I        = I + 1
     ZS       = (J5 - 2) * ZM
     PMESH(I) = ZS + (ZF - ZW(J5-1)) * ZM / (ZW(J5) - ZW(J5-1))
     ZF       = ZF + ZDP
     !
     GOTO 4
     !
5    CONTINUE 
  END DO
  !
  RETURN
END SUBROUTINE PACKMEP
