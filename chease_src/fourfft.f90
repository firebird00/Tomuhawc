!*DECK C2SM22
!*CALL PROCESS
SUBROUTINE FOURFFT(KPSI,PTETA,KMMAX)
  !        ####################################
  !
  !                                        AUTHOR O. SAUTER, CRPP-EPFL
  !**********************************************************************
  !                                                                     *
  ! C2SM22 COMPUTE FAST FOURIER TRANSFORM OF THE DIFFERENT EQL, EQI     *
  !        AND EQ3 TERMS                                                *
  !                                                                     *
  !**********************************************************************
  !
  USE globals
  IMPLICIT NONE
  !
  REAL(RKIND)      ::     PTETA
  INTEGER          ::     IM2P
  INTEGER          ::     IM2
  INTEGER          ::     IMP1
  INTEGER          ::     IP
  REAL(RKIND)      ::     ZN
  REAL(RKIND)      ::     ZFFTEQ3
  REAL(RKIND)      ::     ZFFTEQI
  REAL(RKIND)      ::     ZFFTEQL
  REAL(RKIND)      ::     ZWORK
  REAL(RKIND)      ::     ZD2FUN
  REAL(RKIND)      ::     ZD
  REAL(RKIND)      ::     ZC
  REAL(RKIND)      ::     ZB
  REAL(RKIND)      ::     ZA
  REAL(RKIND)      ::     ZH
  INTEGER          ::     ICHIISO
  INTEGER          ::     ISRCHFGE
  INTEGER          ::     ICHISO
  INTEGER          ::     IGCHISO
  REAL(RKIND)      ::     ZCHIFFT
  INTEGER          ::     I
  REAL(RKIND)      ::     ZDCHI
  INTEGER          ::     KMMAX
  INTEGER          ::     KPSI
  INTEGER          ::     INCHI
  REAL(RKIND)      ::     ZLOG2
  DIMENSION &
       &     PTETA(*), ICHIISO(2*NPCHI1), &
       &     ZCHIFFT(2*NPCHI1), ZD2FUN(NPMGS*NTP1), ZWORK(NPMGS*NTP1,3), &
       &     ZA(2*NPCHI1), ZB(2*NPCHI1), ZC(2*NPCHI1), ZD(2*NPCHI1), &
       &     ZFFTEQL(2*NPCHI1,19), ZFFTEQI(2*NPCHI1,32), &
       &     ZFFTEQ3(2*NPCHI1,17)
  !
  !----*----*----*---*----*----*----*----*----*----*----*----*----*----*-
  !
  !     1. COMPUTE EQUIDISTANT CHI MESH, 2**L INTERVALS
  !
  ZLOG2 = LOG(1._RKIND*NCHI)/LOG(2._RKIND)
  INCHI = 2**NINT(ZLOG2)
  IF (INCHI .GE. 2*NPCHI1) INCHI = INCHI/2
  !
  IF ((KPSI .EQ. 1) .AND. (NVERBOSE .GE. 1)) THEN
    PRINT *,' '
    PRINT *,' NUMBER OF CHI POINTS USED FOR FOURIER= ',INCHI
    PRINT *,' '
  ENDIF

  IF (KMMAX .GT. INCHI/2) THEN
    IF (NVERBOSE .GE. 0) PRINT *,' NOT ENOUGH POINTS FOR FFT: MMAX= ',KMMAX, &
      &       ' INCHI/2= ',INCHI/2
    STOP 'fourfft'
  ENDIF

  ZDCHI = 2._RKIND * CPI / REAL(INCHI,RKIND)
  DO I=1,INCHI+1
     ZCHIFFT(I) = REAL(I-1,RKIND)*ZDCHI
  ENDDO
  !
  !     2. PREPARE COEFFICIENTS FOR THE CUBIC SPLINE FIT DEPENDING
  !     ONLY ON RELATIVE POSITION OF ZCHIFFT(I) WITH RESPECT TO CHIISO
  !
  IGCHISO = NMGAUS*NT1
  CALL GCHI(KPSI)
  DO I=1,INCHI
     ICHISO = ISRCHFGE(IGCHISO,CHIISO,1,ZCHIFFT(I)) - 1
     !
     IF (ICHISO .LT. 1) THEN
        ICHISO = 1
     ELSE IF (ICHISO .GT. IGCHISO) THEN
        ICHISO = IGCHISO
     ENDIF
     ICHIISO(I) = ICHISO
     !
     ZH = CHIISO(ICHISO+1) - CHIISO(ICHISO)
     ZA(I) = (CHIISO(ICHISO+1) - ZCHIFFT(I)) / ZH
     ZB(I) = (ZCHIFFT(I) - CHIISO(ICHISO)) / ZH
     ZC(I) = (ZA(I) + 1) * (ZA(I) - 1) * ZH * &
          &       (CHIISO(ICHISO+1) - ZCHIFFT(I)) / 6._RKIND
     ZD(I) = (ZB(I) + 1) * (ZB(I) - 1) * ZH * &
          &       (ZCHIFFT(I) - CHIISO(ICHISO)) / 6._RKIND
     !
  ENDDO
  !
  !     3. FOR EACH ARRAY: COMPUTE VALUES ON ZCHIFFT USING A PERIODIC
  !     CUBIC SPLINE FIT AND COMPUTE FULL FOURIER TRANSFORM
  !
  !     EQL
  !
  CALL SPLIFFT(CHIISO,EQL(1,1),IGCHISO,RC2PI,ZD2FUN, &
       &     ZWORK,ZFFTEQL(1,1),INCHI,ICHIISO,ZA,ZB,ZC,ZD)
  CALL SPLIFFT(CHIISO,EQL(1,2),IGCHISO,RC2PI,ZD2FUN, &
       &     ZWORK,ZFFTEQL(1,2),INCHI,ICHIISO,ZA,ZB,ZC,ZD)
  CALL SPLIFFT(CHIISO,EQL(1,3),IGCHISO,RC2PI,ZD2FUN, &
       &     ZWORK,ZFFTEQL(1,3),INCHI,ICHIISO,ZA,ZB,ZC,ZD)
  CALL SPLIFFT(CHIISO,EQL(1,4),IGCHISO,RC2PI,ZD2FUN, &
       &     ZWORK,ZFFTEQL(1,4),INCHI,ICHIISO,ZA,ZB,ZC,ZD)
  CALL SPLIFFT(CHIISO,EQL(1,5),IGCHISO,RC2PI,ZD2FUN, &
       &     ZWORK,ZFFTEQL(1,5),INCHI,ICHIISO,ZA,ZB,ZC,ZD)
  CALL SPLIFFT(CHIISO,EQL(1,6),IGCHISO,RC2PI,ZD2FUN, &
       &     ZWORK,ZFFTEQL(1,6),INCHI,ICHIISO,ZA,ZB,ZC,ZD)
  CALL SPLIFFT(CHIISO,EQL(1,7),IGCHISO,RC2PI,ZD2FUN, &
       &     ZWORK,ZFFTEQL(1,7),INCHI,ICHIISO,ZA,ZB,ZC,ZD)
  CALL SPLIFFT(CHIISO,EQL(1,8),IGCHISO,RC2PI,ZD2FUN, &
       &     ZWORK,ZFFTEQL(1,8),INCHI,ICHIISO,ZA,ZB,ZC,ZD)
  CALL SPLIFFT(CHIISO,EQL(1,9),IGCHISO,RC2PI,ZD2FUN, &
       &     ZWORK,ZFFTEQL(1,9),INCHI,ICHIISO,ZA,ZB,ZC,ZD)
  CALL SPLIFFT(CHIISO,EQL(1,10),IGCHISO,RC2PI,ZD2FUN, &
       &     ZWORK,ZFFTEQL(1,10),INCHI,ICHIISO,ZA,ZB,ZC,ZD)
  CALL SPLIFFT(CHIISO,EQL(1,11),IGCHISO,RC2PI,ZD2FUN, &
       &     ZWORK,ZFFTEQL(1,11),INCHI,ICHIISO,ZA,ZB,ZC,ZD)
  CALL SPLIFFT(CHIISO,EQL(1,12),IGCHISO,RC2PI,ZD2FUN, &
       &     ZWORK,ZFFTEQL(1,12),INCHI,ICHIISO,ZA,ZB,ZC,ZD)
  CALL SPLIFFT(CHIISO,EQL(1,13),IGCHISO,RC2PI,ZD2FUN, &
       &     ZWORK,ZFFTEQL(1,13),INCHI,ICHIISO,ZA,ZB,ZC,ZD)
  CALL SPLIFFT(CHIISO,EQL(1,14),IGCHISO,RC2PI,ZD2FUN, &
       &     ZWORK,ZFFTEQL(1,14),INCHI,ICHIISO,ZA,ZB,ZC,ZD)
  CALL SPLIFFT(CHIISO,EQL(1,15),IGCHISO,RC2PI,ZD2FUN, &
       &     ZWORK,ZFFTEQL(1,15),INCHI,ICHIISO,ZA,ZB,ZC,ZD)
  CALL SPLIFFT(CHIISO,EQL(1,16),IGCHISO,RC2PI,ZD2FUN, &
       &     ZWORK,ZFFTEQL(1,16),INCHI,ICHIISO,ZA,ZB,ZC,ZD)
  CALL SPLIFFT(CHIISO,EQL(1,17),IGCHISO,RC2PI,ZD2FUN, &
       &     ZWORK,ZFFTEQL(1,17),INCHI,ICHIISO,ZA,ZB,ZC,ZD)
  CALL SPLIFFT(CHIISO,EQL(1,18),IGCHISO,RC2PI,ZD2FUN, &
       &     ZWORK,ZFFTEQL(1,18),INCHI,ICHIISO,ZA,ZB,ZC,ZD)
  CALL SPLIFFT(CHIISO,EQL(1,19),IGCHISO,RC2PI,ZD2FUN, &
       &     ZWORK,ZFFTEQL(1,19),INCHI,ICHIISO,ZA,ZB,ZC,ZD)
  !
  !     EQI
  !
  CALL SPLIFFT(CHIISO,EQI(1,1),IGCHISO,RC2PI,ZD2FUN, &
       &     ZWORK,ZFFTEQI(1,1),INCHI,ICHIISO,ZA,ZB,ZC,ZD)
  CALL SPLIFFT(CHIISO,EQI(1,2),IGCHISO,RC2PI,ZD2FUN, &
       &     ZWORK,ZFFTEQI(1,2),INCHI,ICHIISO,ZA,ZB,ZC,ZD)
  CALL SPLIFFT(CHIISO,EQI(1,3),IGCHISO,RC2PI,ZD2FUN, &
       &     ZWORK,ZFFTEQI(1,3),INCHI,ICHIISO,ZA,ZB,ZC,ZD)
  CALL SPLIFFT(CHIISO,EQI(1,4),IGCHISO,RC2PI,ZD2FUN, &
       &     ZWORK,ZFFTEQI(1,4),INCHI,ICHIISO,ZA,ZB,ZC,ZD)
  CALL SPLIFFT(CHIISO,EQI(1,5),IGCHISO,RC2PI,ZD2FUN, &
       &     ZWORK,ZFFTEQI(1,5),INCHI,ICHIISO,ZA,ZB,ZC,ZD)
  CALL SPLIFFT(CHIISO,EQI(1,6),IGCHISO,RC2PI,ZD2FUN, &
       &     ZWORK,ZFFTEQI(1,6),INCHI,ICHIISO,ZA,ZB,ZC,ZD)
  CALL SPLIFFT(CHIISO,EQI(1,7),IGCHISO,RC2PI,ZD2FUN, &
       &     ZWORK,ZFFTEQI(1,7),INCHI,ICHIISO,ZA,ZB,ZC,ZD)
  CALL SPLIFFT(CHIISO,EQI(1,8),IGCHISO,RC2PI,ZD2FUN, &
       &     ZWORK,ZFFTEQI(1,8),INCHI,ICHIISO,ZA,ZB,ZC,ZD)
  CALL SPLIFFT(CHIISO,EQI(1,9),IGCHISO,RC2PI,ZD2FUN, &
       &     ZWORK,ZFFTEQI(1,9),INCHI,ICHIISO,ZA,ZB,ZC,ZD)
  CALL SPLIFFT(CHIISO,EQI(1,10),IGCHISO,RC2PI,ZD2FUN, &
       &     ZWORK,ZFFTEQI(1,10),INCHI,ICHIISO,ZA,ZB,ZC,ZD)
  CALL SPLIFFT(CHIISO,EQI(1,11),IGCHISO,RC2PI,ZD2FUN, &
       &     ZWORK,ZFFTEQI(1,11),INCHI,ICHIISO,ZA,ZB,ZC,ZD)
  CALL SPLIFFT(CHIISO,EQI(1,12),IGCHISO,RC2PI,ZD2FUN, &
       &     ZWORK,ZFFTEQI(1,12),INCHI,ICHIISO,ZA,ZB,ZC,ZD)
  CALL SPLIFFT(CHIISO,EQI(1,13),IGCHISO,RC2PI,ZD2FUN, &
       &     ZWORK,ZFFTEQI(1,13),INCHI,ICHIISO,ZA,ZB,ZC,ZD)
  CALL SPLIFFT(CHIISO,EQI(1,14),IGCHISO,RC2PI,ZD2FUN, &
       &     ZWORK,ZFFTEQI(1,14),INCHI,ICHIISO,ZA,ZB,ZC,ZD)
  CALL SPLIFFT(CHIISO,EQI(1,15),IGCHISO,RC2PI,ZD2FUN, &
       &     ZWORK,ZFFTEQI(1,15),INCHI,ICHIISO,ZA,ZB,ZC,ZD)
  CALL SPLIFFT(CHIISO,EQI(1,16),IGCHISO,RC2PI,ZD2FUN, &
       &     ZWORK,ZFFTEQI(1,16),INCHI,ICHIISO,ZA,ZB,ZC,ZD)
  CALL SPLIFFT(CHIISO,EQI(1,17),IGCHISO,RC2PI,ZD2FUN, &
       &     ZWORK,ZFFTEQI(1,17),INCHI,ICHIISO,ZA,ZB,ZC,ZD)
  CALL SPLIFFT(CHIISO,EQI(1,18),IGCHISO,RC2PI,ZD2FUN, &
       &     ZWORK,ZFFTEQI(1,18),INCHI,ICHIISO,ZA,ZB,ZC,ZD)
  CALL SPLIFFT(CHIISO,EQI(1,19),IGCHISO,RC2PI,ZD2FUN, &
       &     ZWORK,ZFFTEQI(1,19),INCHI,ICHIISO,ZA,ZB,ZC,ZD)
  CALL SPLIFFT(CHIISO,EQI(1,20),IGCHISO,RC2PI,ZD2FUN, &
       &     ZWORK,ZFFTEQI(1,20),INCHI,ICHIISO,ZA,ZB,ZC,ZD)
  CALL SPLIFFT(CHIISO,EQI(1,21),IGCHISO,RC2PI,ZD2FUN, &
       &     ZWORK,ZFFTEQI(1,21),INCHI,ICHIISO,ZA,ZB,ZC,ZD)
  CALL SPLIFFT(CHIISO,EQI(1,22),IGCHISO,RC2PI,ZD2FUN, &
       &     ZWORK,ZFFTEQI(1,22),INCHI,ICHIISO,ZA,ZB,ZC,ZD)
  CALL SPLIFFT(CHIISO,EQI(1,23),IGCHISO,RC2PI,ZD2FUN, &
       &     ZWORK,ZFFTEQI(1,23),INCHI,ICHIISO,ZA,ZB,ZC,ZD)
  CALL SPLIFFT(CHIISO,EQI(1,24),IGCHISO,RC2PI,ZD2FUN, &
       &     ZWORK,ZFFTEQI(1,24),INCHI,ICHIISO,ZA,ZB,ZC,ZD)
  CALL SPLIFFT(CHIISO,EQI(1,25),IGCHISO,RC2PI,ZD2FUN, &
       &     ZWORK,ZFFTEQI(1,25),INCHI,ICHIISO,ZA,ZB,ZC,ZD)
  CALL SPLIFFT(CHIISO,EQI(1,26),IGCHISO,RC2PI,ZD2FUN, &
       &     ZWORK,ZFFTEQI(1,26),INCHI,ICHIISO,ZA,ZB,ZC,ZD)
  CALL SPLIFFT(CHIISO,EQI(1,27),IGCHISO,RC2PI,ZD2FUN, &
       &     ZWORK,ZFFTEQI(1,27),INCHI,ICHIISO,ZA,ZB,ZC,ZD)
  CALL SPLIFFT(CHIISO,EQI(1,28),IGCHISO,RC2PI,ZD2FUN, &
       &     ZWORK,ZFFTEQI(1,28),INCHI,ICHIISO,ZA,ZB,ZC,ZD)
  CALL SPLIFFT(CHIISO,EQI(1,29),IGCHISO,RC2PI,ZD2FUN, &
       &     ZWORK,ZFFTEQI(1,29),INCHI,ICHIISO,ZA,ZB,ZC,ZD)
  CALL SPLIFFT(CHIISO,EQI(1,30),IGCHISO,RC2PI,ZD2FUN, &
       &     ZWORK,ZFFTEQI(1,30),INCHI,ICHIISO,ZA,ZB,ZC,ZD)
  CALL SPLIFFT(CHIISO,EQI(1,31),IGCHISO,RC2PI,ZD2FUN, &
       &     ZWORK,ZFFTEQI(1,31),INCHI,ICHIISO,ZA,ZB,ZC,ZD)
  CALL SPLIFFT(CHIISO,EQI(1,32),IGCHISO,RC2PI,ZD2FUN, &
       &     ZWORK,ZFFTEQI(1,32),INCHI,ICHIISO,ZA,ZB,ZC,ZD)
  !
  !     EQ3
  !
  CALL SPLIFFT(CHIISO,EQ3(1,1),IGCHISO,RC2PI,ZD2FUN, &
       &     ZWORK,ZFFTEQ3(1,1),INCHI,ICHIISO,ZA,ZB,ZC,ZD)
  CALL SPLIFFT(CHIISO,EQ3(1,2),IGCHISO,RC2PI,ZD2FUN, &
       &     ZWORK,ZFFTEQ3(1,2),INCHI,ICHIISO,ZA,ZB,ZC,ZD)
  CALL SPLIFFT(CHIISO,EQ3(1,3),IGCHISO,RC2PI,ZD2FUN, &
       &     ZWORK,ZFFTEQ3(1,3),INCHI,ICHIISO,ZA,ZB,ZC,ZD)
  CALL SPLIFFT(CHIISO,EQ3(1,4),IGCHISO,RC2PI,ZD2FUN, &
       &     ZWORK,ZFFTEQ3(1,4),INCHI,ICHIISO,ZA,ZB,ZC,ZD)
  CALL SPLIFFT(CHIISO,EQ3(1,5),IGCHISO,RC2PI,ZD2FUN, &
       &     ZWORK,ZFFTEQ3(1,5),INCHI,ICHIISO,ZA,ZB,ZC,ZD)
  CALL SPLIFFT(CHIISO,EQ3(1,6),IGCHISO,RC2PI,ZD2FUN, &
       &     ZWORK,ZFFTEQ3(1,6),INCHI,ICHIISO,ZA,ZB,ZC,ZD)
  CALL SPLIFFT(CHIISO,EQ3(1,7),IGCHISO,RC2PI,ZD2FUN, &
       &     ZWORK,ZFFTEQ3(1,7),INCHI,ICHIISO,ZA,ZB,ZC,ZD)
  CALL SPLIFFT(CHIISO,EQ3(1,8),IGCHISO,RC2PI,ZD2FUN, &
       &     ZWORK,ZFFTEQ3(1,8),INCHI,ICHIISO,ZA,ZB,ZC,ZD)
  CALL SPLIFFT(CHIISO,EQ3(1,9),IGCHISO,RC2PI,ZD2FUN, &
       &     ZWORK,ZFFTEQ3(1,9),INCHI,ICHIISO,ZA,ZB,ZC,ZD)
  CALL SPLIFFT(CHIISO,EQ3(1,10),IGCHISO,RC2PI,ZD2FUN, &
       &     ZWORK,ZFFTEQ3(1,10),INCHI,ICHIISO,ZA,ZB,ZC,ZD)
  CALL SPLIFFT(CHIISO,EQ3(1,11),IGCHISO,RC2PI,ZD2FUN, &
       &     ZWORK,ZFFTEQ3(1,11),INCHI,ICHIISO,ZA,ZB,ZC,ZD)
  CALL SPLIFFT(CHIISO,EQ3(1,12),IGCHISO,RC2PI,ZD2FUN, &
       &     ZWORK,ZFFTEQ3(1,12),INCHI,ICHIISO,ZA,ZB,ZC,ZD)
  CALL SPLIFFT(CHIISO,EQ3(1,13),IGCHISO,RC2PI,ZD2FUN, &
       &     ZWORK,ZFFTEQ3(1,13),INCHI,ICHIISO,ZA,ZB,ZC,ZD)
  CALL SPLIFFT(CHIISO,EQ3(1,14),IGCHISO,RC2PI,ZD2FUN, &
       &     ZWORK,ZFFTEQ3(1,14),INCHI,ICHIISO,ZA,ZB,ZC,ZD)
  CALL SPLIFFT(CHIISO,EQ3(1,15),IGCHISO,RC2PI,ZD2FUN, &
       &     ZWORK,ZFFTEQ3(1,15),INCHI,ICHIISO,ZA,ZB,ZC,ZD)
  CALL SPLIFFT(CHIISO,EQ3(1,16),IGCHISO,RC2PI,ZD2FUN, &
       &     ZWORK,ZFFTEQ3(1,16),INCHI,ICHIISO,ZA,ZB,ZC,ZD)
  CALL SPLIFFT(CHIISO,EQ3(1,17),IGCHISO,RC2PI,ZD2FUN, &
       &     ZWORK,ZFFTEQ3(1,17),INCHI,ICHIISO,ZA,ZB,ZC,ZD)
  !-----------------------------------------------------------------------
  !
  !     4. COPY SPECTRUM INTO WANTED ARRAYS
  !     F_M = CMPLX(FFT(2*M+1), -FFT(2*M+2)), THUS MMAX<INCHI
  !     AS USE EXP(-I*K*X) INSTEAD OF EXP(I*K*X) IN REALFT, MINUS
  !     SIGN FOR IMAGINARY PART
  !
  ZN = 1._RKIND / REAL(INCHI,RKIND)
  !
  IF (MOD(KPSI,2) .EQ. 0) THEN
     IP = KPSI / 2 + 1
     DO IMP1=1,KMMAX+1
        IM2 = 2*IMP1 - 1
        IM2P = IM2 + 1
        !     EQL
        DG11L(IP,IMP1) =CMPLX(ZFFTEQL(IM2, 1),-ZFFTEQL(IM2P, 1))*ZN
        DG22L(IP,IMP1) =CMPLX(ZFFTEQL(IM2, 2),-ZFFTEQL(IM2P, 2))*ZN
        DG33L(IP,IMP1) =CMPLX(ZFFTEQL(IM2, 3),-ZFFTEQL(IM2P, 3))*ZN
        DG12L(IP,IMP1) =CMPLX(ZFFTEQL(IM2, 4),-ZFFTEQL(IM2P, 4))*ZN
        JG11L(IP,IMP1) =CMPLX(ZFFTEQL(IM2, 5),-ZFFTEQL(IM2P, 5))*ZN
        JG22L(IP,IMP1) =CMPLX(ZFFTEQL(IM2, 6),-ZFFTEQL(IM2P, 6))*ZN
        JG33L(IP,IMP1) =CMPLX(ZFFTEQL(IM2, 7),-ZFFTEQL(IM2P, 7))*ZN
        JG12L(IP,IMP1) =CMPLX(ZFFTEQL(IM2, 8),-ZFFTEQL(IM2P, 8))*ZN
        JACOBI(IP,IMP1)=CMPLX(ZFFTEQL(IM2, 9),-ZFFTEQL(IM2P, 9))*ZN
        J2U(IP,IMP1,1) =CMPLX(ZFFTEQL(IM2,10),-ZFFTEQL(IM2P,10))*ZN
        J3U(IP,IMP1,1) =CMPLX(ZFFTEQL(IM2,11),-ZFFTEQL(IM2P,11))*ZN
        B2E(IP,IMP1,1) =CMPLX(ZFFTEQL(IM2,12),-ZFFTEQL(IM2P,12))*ZN
        B3E(IP,IMP1,1) =CMPLX(ZFFTEQL(IM2,13),-ZFFTEQL(IM2P,13))*ZN
        PEQ(IP,IMP1,1) =CMPLX(ZFFTEQL(IM2,14),-ZFFTEQL(IM2P,14))*ZN
        DPEDS(IP,IMP1,1)=CMPLX(ZFFTEQL(IM2,15),-ZFFTEQL(IM2P,15)) &
             &         *ZN
        GCHDZ(IP,IMP1) =CMPLX(ZFFTEQL(IM2,16),-ZFFTEQL(IM2P,16))*ZN
        GSDZ (IP,IMP1) =CMPLX(ZFFTEQL(IM2,17),-ZFFTEQL(IM2P,17))*ZN
        GBZ(IP,IMP1)   =CMPLX(ZFFTEQL(IM2,18),-ZFFTEQL(IM2P,18))*ZN
        GBR(IP,IMP1)   =CMPLX(ZFFTEQL(IM2,19),-ZFFTEQL(IM2P,19))*ZN
        !     EQI
        IDIY2(IP,IMP1) =CMPLX(ZFFTEQI(IM2, 1),-ZFFTEQI(IM2P, 1))*ZN
        IDIY3(IP,IMP1) =CMPLX(ZFFTEQI(IM2, 2),-ZFFTEQI(IM2P, 2))*ZN
        IG122(IP,IMP1) =CMPLX(ZFFTEQI(IM2, 3),-ZFFTEQI(IM2P, 3))*ZN
        IG123(IP,IMP1) =CMPLX(ZFFTEQI(IM2, 4),-ZFFTEQI(IM2P, 4))*ZN
        INXX(IP,IMP1)  =CMPLX(ZFFTEQI(IM2, 5),-ZFFTEQI(IM2P, 5))*ZN
        INXY(IP,IMP1)  =CMPLX(ZFFTEQI(IM2, 6),-ZFFTEQI(IM2P, 6))*ZN
        INYY(IP,IMP1)  =CMPLX(ZFFTEQI(IM2, 7),-ZFFTEQI(IM2P, 7))*ZN
        INZZ(IP,IMP1)  =CMPLX(ZFFTEQI(IM2, 8),-ZFFTEQI(IM2P, 8))*ZN
        IJ0QX(IP,IMP1) =CMPLX(ZFFTEQI(IM2, 9),-ZFFTEQI(IM2P, 9))*ZN
        IJ0QY(IP,IMP1) =CMPLX(ZFFTEQI(IM2,10),-ZFFTEQI(IM2P,10))*ZN
        IGPX2(IP,IMP1) =CMPLX(ZFFTEQI(IM2,11),-ZFFTEQI(IM2P,11))*ZN
        IGPX3(IP,IMP1) =CMPLX(ZFFTEQI(IM2,12),-ZFFTEQI(IM2P,12))*ZN
        IGPY2(IP,IMP1) =CMPLX(ZFFTEQI(IM2,13),-ZFFTEQI(IM2P,13))*ZN
        IGPY3(IP,IMP1) =CMPLX(ZFFTEQI(IM2,14),-ZFFTEQI(IM2P,14))*ZN
        IDRXX(IP,IMP1) =CMPLX(ZFFTEQI(IM2,15),-ZFFTEQI(IM2P,15))*ZN
        IRXZ(IP,IMP1)  =CMPLX(ZFFTEQI(IM2,16),-ZFFTEQI(IM2P,16))*ZN
        IDRYX(IP,IMP1) =CMPLX(ZFFTEQI(IM2,17),-ZFFTEQI(IM2P,17))*ZN
        IRYX(IP,IMP1)  =CMPLX(ZFFTEQI(IM2,18),-ZFFTEQI(IM2P,18))*ZN
        IDRZX(IP,IMP1) =CMPLX(ZFFTEQI(IM2,19),-ZFFTEQI(IM2P,19))*ZN
        IRZY(IP,IMP1)  =CMPLX(ZFFTEQI(IM2,20),-ZFFTEQI(IM2P,20))*ZN
        VISXZ(IP,IMP1) =CMPLX(ZFFTEQI(IM2,21),-ZFFTEQI(IM2P,21))*ZN
        VISYZ(IP,IMP1) =CMPLX(ZFFTEQI(IM2,22),-ZFFTEQI(IM2P,22))*ZN
        IVS11(IP,IMP1) =CMPLX(ZFFTEQI(IM2,23),-ZFFTEQI(IM2P,23))*ZN
        IVS12(IP,IMP1) =CMPLX(ZFFTEQI(IM2,24),-ZFFTEQI(IM2P,24))*ZN
        IVS21(IP,IMP1) =CMPLX(ZFFTEQI(IM2,25),-ZFFTEQI(IM2P,25))*ZN
        IVS22(IP,IMP1) =CMPLX(ZFFTEQI(IM2,26),-ZFFTEQI(IM2P,26))*ZN
        GSFC(IP,IMP1)  =CMPLX(ZFFTEQI(IM2,27),-ZFFTEQI(IM2P,27))*ZN
        GSCC(IP,IMP1)  =CMPLX(ZFFTEQI(IM2,28),-ZFFTEQI(IM2P,28))*ZN
        GSFS(IP,IMP1)  =CMPLX(ZFFTEQI(IM2,29),-ZFFTEQI(IM2P,29))*ZN
        GSCS(IP,IMP1)  =CMPLX(ZFFTEQI(IM2,30),-ZFFTEQI(IM2P,30))*ZN
        GCFC(IP,IMP1)  =CMPLX(ZFFTEQI(IM2,31),-ZFFTEQI(IM2P,31))*ZN
        GCFS(IP,IMP1)  =CMPLX(ZFFTEQI(IM2,32),-ZFFTEQI(IM2P,32))*ZN
        !     EQ3
        EQRHO(IP,IMP1) =CMPLX(ZFFTEQ3(IM2, 1),-ZFFTEQ3(IM2P, 1))*ZN
        DRHOS(IP,IMP1) =CMPLX(ZFFTEQ3(IM2, 2),-ZFFTEQ3(IM2P, 2))*ZN
        EQROT(IP,IMP1) =CMPLX(ZFFTEQ3(IM2, 3),-ZFFTEQ3(IM2P, 3))*ZN
        DROT(IP,IMP1)  =CMPLX(ZFFTEQ3(IM2, 4),-ZFFTEQ3(IM2P, 4))*ZN
        FEQ(IP,IMP1)   =CMPLX(ZFFTEQ3(IM2, 5),-ZFFTEQ3(IM2P, 5))*ZN
        IWSQ1(IP,IMP1) =CMPLX(ZFFTEQ3(IM2, 6),-ZFFTEQ3(IM2P, 6))*ZN
        IWSQ2(IP,IMP1) =CMPLX(ZFFTEQ3(IM2, 7),-ZFFTEQ3(IM2P, 7))*ZN
        IWSQ3(IP,IMP1) =CMPLX(ZFFTEQ3(IM2, 8),-ZFFTEQ3(IM2P, 8))*ZN
        IJ0QZ(IP,IMP1) =CMPLX(ZFFTEQ3(IM2, 9),-ZFFTEQ3(IM2P, 9))*ZN
        JACOF(IP,IMP1) =CMPLX(ZFFTEQ3(IM2,10),-ZFFTEQ3(IM2P,10))*ZN
        B2F(IP,IMP1)   =CMPLX(ZFFTEQ3(IM2,11),-ZFFTEQ3(IM2P,11))*ZN
        B3F(IP,IMP1)   =CMPLX(ZFFTEQ3(IM2,12),-ZFFTEQ3(IM2P,12))*ZN
        JACOS(IP,IMP1) =CMPLX(ZFFTEQ3(IM2,13),-ZFFTEQ3(IM2P,13))*ZN
        IGF22(IP,IMP1) =CMPLX(ZFFTEQ3(IM2,14),-ZFFTEQ3(IM2P,14))*ZN
        B3FC(IP,IMP1)  =CMPLX(ZFFTEQ3(IM2,15),-ZFFTEQ3(IM2P,15))*ZN
        B2FC(IP,IMP1)  =CMPLX(ZFFTEQ3(IM2,16),-ZFFTEQ3(IM2P,16))*ZN
        DJCOF(IP,IMP1) =CMPLX(ZFFTEQ3(IM2,17),-ZFFTEQ3(IM2P,17))*ZN
     ENDDO
  ELSE IF (MOD(KPSI,2) .EQ. 1) THEN
     IP = (KPSI + 1) / 2
     DO IMP1=1,KMMAX+1
        IM2 = 2*IMP1 - 1
        IM2P = IM2 + 1
        !     EQL
        DG11LM(IP,IMP1)=CMPLX(ZFFTEQL(IM2, 1),-ZFFTEQL(IM2P, 1))*ZN
        DG22LM(IP,IMP1)=CMPLX(ZFFTEQL(IM2, 2),-ZFFTEQL(IM2P, 2))*ZN
        DG33LM(IP,IMP1)=CMPLX(ZFFTEQL(IM2, 3),-ZFFTEQL(IM2P, 3))*ZN
        DG12LM(IP,IMP1)=CMPLX(ZFFTEQL(IM2, 4),-ZFFTEQL(IM2P, 4))*ZN
        JG11LM(IP,IMP1)=CMPLX(ZFFTEQL(IM2, 5),-ZFFTEQL(IM2P, 5))*ZN
        JG22LM(IP,IMP1)=CMPLX(ZFFTEQL(IM2, 6),-ZFFTEQL(IM2P, 6))*ZN
        JG33LM(IP,IMP1)=CMPLX(ZFFTEQL(IM2, 7),-ZFFTEQL(IM2P, 7))*ZN
        JG12LM(IP,IMP1)=CMPLX(ZFFTEQL(IM2, 8),-ZFFTEQL(IM2P, 8))*ZN
        JACOBM(IP,IMP1)=CMPLX(ZFFTEQL(IM2, 9),-ZFFTEQL(IM2P, 9))*ZN
        J2E(IP,IMP1,1) =CMPLX(ZFFTEQL(IM2,10),-ZFFTEQL(IM2P,10))*ZN
        J3E(IP,IMP1,1) =CMPLX(ZFFTEQL(IM2,11),-ZFFTEQL(IM2P,11))*ZN
        B2U(IP,IMP1,1) =CMPLX(ZFFTEQL(IM2,12),-ZFFTEQL(IM2P,12))*ZN
        B3U(IP,IMP1,1) =CMPLX(ZFFTEQL(IM2,13),-ZFFTEQL(IM2P,13))*ZN
        PRE(IP,IMP1,1) =CMPLX(ZFFTEQL(IM2,14),-ZFFTEQL(IM2P,14))*ZN
        DPEDSM(IP,IMP1,1)=CMPLX(ZFFTEQL(IM2,15),-ZFFTEQL(IM2P,15)) &
             &         *ZN
        GCHDZM(IP,IMP1)=CMPLX(ZFFTEQL(IM2,16),-ZFFTEQL(IM2P,16))*ZN
        GSDZM (IP,IMP1)=CMPLX(ZFFTEQL(IM2,17),-ZFFTEQL(IM2P,17))*ZN
        GBZM(IP,IMP1)  =CMPLX(ZFFTEQL(IM2,18),-ZFFTEQL(IM2P,18))*ZN
        GBRM(IP,IMP1)  =CMPLX(ZFFTEQL(IM2,19),-ZFFTEQL(IM2P,19))*ZN
        !     EQI
        IDIY2M(IP,IMP1)=CMPLX(ZFFTEQI(IM2, 1),-ZFFTEQI(IM2P, 1))*ZN
        IDIY3M(IP,IMP1)=CMPLX(ZFFTEQI(IM2, 2),-ZFFTEQI(IM2P, 2))*ZN
        IG122M(IP,IMP1)=CMPLX(ZFFTEQI(IM2, 3),-ZFFTEQI(IM2P, 3))*ZN
        IG123M(IP,IMP1)=CMPLX(ZFFTEQI(IM2, 4),-ZFFTEQI(IM2P, 4))*ZN
        INXXM(IP,IMP1) =CMPLX(ZFFTEQI(IM2, 5),-ZFFTEQI(IM2P, 5))*ZN
        INXYM(IP,IMP1) =CMPLX(ZFFTEQI(IM2, 6),-ZFFTEQI(IM2P, 6))*ZN
        INYYM(IP,IMP1) =CMPLX(ZFFTEQI(IM2, 7),-ZFFTEQI(IM2P, 7))*ZN
        INZZM(IP,IMP1) =CMPLX(ZFFTEQI(IM2, 8),-ZFFTEQI(IM2P, 8))*ZN
        IJ0QXM(IP,IMP1)=CMPLX(ZFFTEQI(IM2, 9),-ZFFTEQI(IM2P, 9))*ZN
        IJ0QYM(IP,IMP1)=CMPLX(ZFFTEQI(IM2,10),-ZFFTEQI(IM2P,10))*ZN
        IGPX2M(IP,IMP1)=CMPLX(ZFFTEQI(IM2,11),-ZFFTEQI(IM2P,11))*ZN
        IGPX3M(IP,IMP1)=CMPLX(ZFFTEQI(IM2,12),-ZFFTEQI(IM2P,12))*ZN
        IGPY2M(IP,IMP1)=CMPLX(ZFFTEQI(IM2,13),-ZFFTEQI(IM2P,13))*ZN
        IGPY3M(IP,IMP1)=CMPLX(ZFFTEQI(IM2,14),-ZFFTEQI(IM2P,14))*ZN
        IDRXXM(IP,IMP1)=CMPLX(ZFFTEQI(IM2,15),-ZFFTEQI(IM2P,15))*ZN
        IRXZM(IP,IMP1) =CMPLX(ZFFTEQI(IM2,16),-ZFFTEQI(IM2P,16))*ZN
        IDRYXM(IP,IMP1)=CMPLX(ZFFTEQI(IM2,17),-ZFFTEQI(IM2P,17))*ZN
        IRYXM(IP,IMP1) =CMPLX(ZFFTEQI(IM2,18),-ZFFTEQI(IM2P,18))*ZN
        IDRZXM(IP,IMP1)=CMPLX(ZFFTEQI(IM2,19),-ZFFTEQI(IM2P,19))*ZN
        IRZYM(IP,IMP1) =CMPLX(ZFFTEQI(IM2,20),-ZFFTEQI(IM2P,20))*ZN
        VISXZM(IP,IMP1)=CMPLX(ZFFTEQI(IM2,21),-ZFFTEQI(IM2P,21))*ZN
        VISYZM(IP,IMP1)=CMPLX(ZFFTEQI(IM2,22),-ZFFTEQI(IM2P,22))*ZN
        IVS11M(IP,IMP1)=CMPLX(ZFFTEQI(IM2,23),-ZFFTEQI(IM2P,23))*ZN
        IVS12M(IP,IMP1)=CMPLX(ZFFTEQI(IM2,24),-ZFFTEQI(IM2P,24))*ZN
        IVS21M(IP,IMP1)=CMPLX(ZFFTEQI(IM2,25),-ZFFTEQI(IM2P,25))*ZN
        IVS22M(IP,IMP1)=CMPLX(ZFFTEQI(IM2,26),-ZFFTEQI(IM2P,26))*ZN
        GSFCM(IP,IMP1) =CMPLX(ZFFTEQI(IM2,27),-ZFFTEQI(IM2P,27))*ZN
        GSCCM(IP,IMP1) =CMPLX(ZFFTEQI(IM2,28),-ZFFTEQI(IM2P,28))*ZN
        GSFSM(IP,IMP1) =CMPLX(ZFFTEQI(IM2,29),-ZFFTEQI(IM2P,29))*ZN
        GSCSM(IP,IMP1) =CMPLX(ZFFTEQI(IM2,30),-ZFFTEQI(IM2P,30))*ZN
        GCFCM(IP,IMP1) =CMPLX(ZFFTEQI(IM2,31),-ZFFTEQI(IM2P,31))*ZN
        GCFSM(IP,IMP1) =CMPLX(ZFFTEQI(IM2,32),-ZFFTEQI(IM2P,32))*ZN
        !     EQ3
        EQRHOM(IP,IMP1)=CMPLX(ZFFTEQ3(IM2, 1),-ZFFTEQ3(IM2P, 1))*ZN
        DRHOSM(IP,IMP1)=CMPLX(ZFFTEQ3(IM2, 2),-ZFFTEQ3(IM2P, 2))*ZN
        EQROTM(IP,IMP1)=CMPLX(ZFFTEQ3(IM2, 3),-ZFFTEQ3(IM2P, 3))*ZN
        DROTM(IP,IMP1) =CMPLX(ZFFTEQ3(IM2, 4),-ZFFTEQ3(IM2P, 4))*ZN
        FEQM(IP,IMP1)  =CMPLX(ZFFTEQ3(IM2, 5),-ZFFTEQ3(IM2P, 5))*ZN
        IWSQ1M(IP,IMP1)=CMPLX(ZFFTEQ3(IM2, 6),-ZFFTEQ3(IM2P, 6))*ZN
        IWSQ2M(IP,IMP1)=CMPLX(ZFFTEQ3(IM2, 7),-ZFFTEQ3(IM2P, 7))*ZN
        IWSQ3M(IP,IMP1)=CMPLX(ZFFTEQ3(IM2, 8),-ZFFTEQ3(IM2P, 8))*ZN
        IJ0QZM(IP,IMP1)=CMPLX(ZFFTEQ3(IM2, 9),-ZFFTEQ3(IM2P, 9))*ZN
        JACOFM(IP,IMP1)=CMPLX(ZFFTEQ3(IM2,10),-ZFFTEQ3(IM2P,10))*ZN
        B2FM(IP,IMP1)  =CMPLX(ZFFTEQ3(IM2,11),-ZFFTEQ3(IM2P,11))*ZN
        B3FM(IP,IMP1)  =CMPLX(ZFFTEQ3(IM2,12),-ZFFTEQ3(IM2P,12))*ZN
        JACOSM(IP,IMP1)=CMPLX(ZFFTEQ3(IM2,13),-ZFFTEQ3(IM2P,13))*ZN
        IGF22M(IP,IMP1)=CMPLX(ZFFTEQ3(IM2,14),-ZFFTEQ3(IM2P,14))*ZN
        B3FCM(IP,IMP1) =CMPLX(ZFFTEQ3(IM2,15),-ZFFTEQ3(IM2P,15))*ZN
        B2FCM(IP,IMP1) =CMPLX(ZFFTEQ3(IM2,16),-ZFFTEQ3(IM2P,16))*ZN
        DJCOFM(IP,IMP1)=CMPLX(ZFFTEQ3(IM2,17),-ZFFTEQ3(IM2P,17))*ZN
     ENDDO
  ENDIF
  !
  RETURN
END SUBROUTINE FOURFFT
