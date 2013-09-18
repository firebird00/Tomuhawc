!*DECK C2SM09
!*CALL PROCESS
SUBROUTINE XTINIT(NPMAX,PSMISO)
  !        ##############################################
  !
  !                                        AUTHORS:
  !                                        H. LUTJENS,  CRPP-EPFL
  !                                        A. BONDESON, CRPP-EPFL
  !**********************************************************************
  !                                                                     *
  ! C2SM09 FOURIER TRANSFORM THE EQ'S COMPUTED IN GIJLIN ACCORDING TO   *
  !        [1], EQ. (22)                                                *
  !                                                                     *
  !**********************************************************************
  !
  USE globals
  IMPLICIT NONE
  !
  REAL(RKIND)      ::     PSMISO
  REAL(RKIND)      ::     ZBNORM
  REAL(RKIND)      ::     ZF
  REAL(RKIND)      ::     ZJAC
  REAL(RKIND)      ::     ZBND
  REAL(RKIND)      ::     ZRHO
  REAL(RKIND)      ::     ZGP
  REAL(RKIND)      ::     ZPSI
  REAL(RKIND)      ::     ZR
  REAL(RKIND)      ::     ZZ
  REAL(RKIND)      ::     ZCP1
  REAL(RKIND)      ::     ZD
  REAL(RKIND)      ::     ZC
  REAL(RKIND)      ::     ZB
  REAL(RKIND)      ::     ZA
  REAL(RKIND)      ::     ZH
  REAL(RKIND)      ::     ZC1
  REAL(RKIND)      ::     ZB1
  REAL(RKIND)      ::     ZA1
  REAL(RKIND)      ::     equilibriumdensity
  REAL(RKIND)      ::     ZDPDS
  REAL(RKIND)      ::     ZDPDT
  REAL(RKIND)      ::     ZDRSDT
  REAL(RKIND)      ::     ZFP
  REAL(RKIND)      ::     ZEPS
  INTEGER          ::     IC
  INTEGER          ::     ICHIM
  INTEGER          ::     IS0
  INTEGER          ::     IT0
  INTEGER          ::     J
  INTEGER          ::     JG
  INTEGER          ::     JPOP
  INTEGER          ::     JS
  INTEGER          ::     JT
  INTEGER          ::     KGAUS
  INTEGER          ::     KPSI
  INTEGER          ::     NPMAX
  REAL(RKIND)      ::     ZDTHETAHAT
  REAL(RKIND)      ::     THETAHAT
  REAL(RKIND)      ::     ZD2RHO
  REAL(RKIND)      ::     ZD2SIG
  REAL(RKIND)      ::     ZD2TET
  REAL(RKIND)      ::     ZTET
  REAL(RKIND)      ::     RPARXT
  REAL(RKIND)      ::     ZPARXT
  REAL(RKIND)      ::     ZBNDPARXT
  REAL(RKIND)      ::     RHOMAP
  REAL(RKIND)      ::     ZTETA
  DIMENSION &
       &   IC((NPMGS+1)*NPCHI),        IS0((NPMGS+1)*NPCHI),  IT0((NPMGS+1)*NPCHI), &
       &   PSMISO(NPMAX),  &
       &   ZDTHETAHAT((NPMGS+1)*NTP2), THETAHAT(NTP2), &
       &   RPARXT(NPCHI),   ZPARXT(NPCHI),  & 
       &   ZA1(NTP2),    ZB1(NTP2),    ZC1(NTP2),    ZD2RHO(NTP2), &
       &   ZD2SIG(NTP2), ZD2TET(NTP2), ZTET(NTP2),   ZBNDPARXT(NPCHI), & 
       &   RHOMAP(NTP2), ZTETA(NTP2,5), ZBND(NTP2,5), &
       &   ZJAC(NTP2),   ZBNORM(NTP2)
  REAL(RKIND)  :: &
       &   ZDBDS(NTP1,16),  ZDBDT(NTP1,16),  &
       &   ZPCEL(NTP1,16),  ZS(NTP1),        ZS1(NTP1), &
       &   ZS2(NTP1),       ZT(NTP1),        ZT1(NTP1), &
       &   ZT2(NTP1)
  !
  !----*----*----*---*----*----*----*----*----*----*----*----*----*----*-
  !
  DO KPSI=1,NPMAX
     !
     DO J=1,NMGAUS*NT1
        !
        ZDTHETAHAT(J) = RRISO(J,KPSI)*RHOISO(J,KPSI)*BNDISO(J,KPSI)/DPSISO(J,KPSI)
        !
     END DO
     !
     THETAHAT(1) = 0._RKIND
     !
     DO J=1,NT1
        !
        JG = (J - 1) * NMGAUS
        !
        THETAHAT(J+1) = THETAHAT(J)
        !
        DO KGAUS=1,NMGAUS
           THETAHAT(J+1) = THETAHAT(J+1) + WGTPSI(JG+KGAUS,KPSI)*ZDTHETAHAT(JG+KGAUS)
        END DO
        !
     END DO
     !
     ZCP1 = THETAHAT(NT2) / TWOPI
     !
     CALL DSCAL(NT2,RC1P/ZCP1,THETAHAT,1)
     !
     CALL DCOPY(NT2,TETMAP(1,KPSI),1,ZTET,1)
     !
     DO J=2,NT2
        !
        IF (ZTET(J) .LT. ZTET(J-1)) THEN
           !
           ZTET(J) = ZTET(J) + 2._RKIND * CPI * (1._RKIND + &
                &                 INT(.5_RKIND * ABS(ZTET(J) - ZTET(J-1)) / CPI))       
           !
        ENDIF
        !
     END DO
     !
     ZEPS   = 1.E-3_RKIND
     !
     DO J=1,NT1
        ZTETA(J,1) = ZTET(J)
        ZTETA(J,2) = ZTET(J) - 2._RKIND * ZEPS
        ZTETA(J,3) = ZTET(J) -      ZEPS
        ZTETA(J,4) = ZTET(J) +      ZEPS
        ZTETA(J,5) = ZTET(J) + 2._RKIND * ZEPS
     END DO
     !
     CALL BOUND(NT1,ZTETA(1,1),ZBND(1,1))
     CALL BOUND(NT1,ZTETA(1,2),ZBND(1,2))
     CALL BOUND(NT1,ZTETA(1,3),ZBND(1,3))
     CALL BOUND(NT1,ZTETA(1,4),ZBND(1,4))
     CALL BOUND(NT1,ZTETA(1,5),ZBND(1,5))
     !
     CALL RESETI(IC,NT1,1)
     DO JT = 1,NT1
        DO JG=1,NT1
           IF (IC(JG).EQ.1) THEN
              IT0(JG) = JT-1
              IF (ZTET(JG).LE.CT(JT)) IC(JG)  = 0
           ENDIF
        ENDDO
     ENDDO
     CALL RESETI(IC,NT1,1)
     DO JS = 1,NS1
        DO JG=1,NT1
           IF (IC(JG).EQ.1) THEN
              IS0(JG) = JS-1
              IF (SIGMAP(JG,KPSI).LE.CSIG(JS)) IC(JG)  = 0
           ENDIF
        ENDDO
     ENDDO
     !
     DO J=1,NT1
        IF (IS0(J) .GT. NS) IS0(J) = NS
        IF (IS0(J) .LT. 1)  IS0(J) = 1
        IF (IT0(J) .GT. NT) IT0(J) = NT
        IF (IT0(J) .LT. 1)  IT0(J) = 1
        !
        ZT(J)  = ZTET(J)
        ZS(J)  = SIGMAP(J,KPSI)
        ZS1(J) = CSIG(IS0(J))
        ZS2(J) = CSIG(IS0(J)+1)
        ZT1(J) = CT(IT0(J))
        ZT2(J) = CT(IT0(J)+1)
     END DO
     !
     CALL PSICEL(IS0,IT0,NT1,NTP1,ZPCEL,CPSICL)
     CALL BASIS2(NT1,NTP1,ZS1,ZS2,ZT1,ZT2,ZS,ZT,ZDBDS,ZDBDT)
     !
     DO J=1,NT1
        !
        ZDRSDT = (ZBND(J,2) + 8._RKIND*(ZBND(J,4) - ZBND(J,3)) - &
             & ZBND(J,5)) / (12._RKIND * ZEPS)
        !
        ZDPDS = ZDBDS(J, 1) * ZPCEL(J, 1) + &
             &           ZDBDS(J, 2) * ZPCEL(J, 2) + &
             &           ZDBDS(J, 3) * ZPCEL(J, 3) + &
             &           ZDBDS(J, 4) * ZPCEL(J, 4) + &
             &           ZDBDS(J, 5) * ZPCEL(J, 5) + &
             &           ZDBDS(J, 6) * ZPCEL(J, 6) + &
             &           ZDBDS(J, 7) * ZPCEL(J, 7) + &
             &           ZDBDS(J, 8) * ZPCEL(J, 8) + &
             &           ZDBDS(J, 9) * ZPCEL(J, 9) + &
             &           ZDBDS(J,10) * ZPCEL(J,10) + &
             &           ZDBDS(J,11) * ZPCEL(J,11) + &
             &           ZDBDS(J,12) * ZPCEL(J,12) + &
             &           ZDBDS(J,13) * ZPCEL(J,13) + &
             &           ZDBDS(J,14) * ZPCEL(J,14) + &
             &           ZDBDS(J,15) * ZPCEL(J,15) + &
             &           ZDBDS(J,16) * ZPCEL(J,16)
        !
        ZDPDT = ZDBDT(J, 1) * ZPCEL(J, 1) + &
             &           ZDBDT(J, 2) * ZPCEL(J, 2) + &
             &           ZDBDT(J, 3) * ZPCEL(J, 3) + &
             &           ZDBDT(J, 4) * ZPCEL(J, 4) + &
             &           ZDBDT(J, 5) * ZPCEL(J, 5) + &
             &           ZDBDT(J, 6) * ZPCEL(J, 6) + &
             &           ZDBDT(J, 7) * ZPCEL(J, 7) + &
             &           ZDBDT(J, 8) * ZPCEL(J, 8) + &
             &           ZDBDT(J, 9) * ZPCEL(J, 9) + &
             &           ZDBDT(J,10) * ZPCEL(J,10) + &
             &           ZDBDT(J,11) * ZPCEL(J,11) + &
             &           ZDBDT(J,12) * ZPCEL(J,12) + &
             &           ZDBDT(J,13) * ZPCEL(J,13) + &
             &           ZDBDT(J,14) * ZPCEL(J,14) + &
             &           ZDBDT(J,15) * ZPCEL(J,15) + &
             &           ZDBDT(J,16) * ZPCEL(J,16)
        !
        ZRHO   = SIGMAP(J,KPSI) * ZBND(J,1)
        ZR = ZRHO * COS(ZTETA(J,1)) + R0
        !
        ZFP    = (ZDPDS**2 + (ZDPDT / SIGMAP(J,KPSI) - ZDPDS * ZDRSDT / &
             &            ZBND(J,1))**2) / ZBND(J,1)**2
        ZBNORM(J)= SQRT(TMF(KPSI)**2 + ZFP) / ZR
        !
        ! JACOBIAN FOR GAMMA INTEGRATION IN S**2 (GAMMA IN LINEAR IN S**2 CLOS TO AXIS)
        !
        ZJAC(J)   = ZCP1
     END DO
     !
     DO JPOP=1,NPOPULATIONS
        !
        DO J=1,NT1
           !
           RHOMAP(J) = equilibriumdensity(PSIISO(KPSI),ZBNORM(J),JPOP)*ZJAC(J)
           !
        END DO
        RHOMAP(NT2) = RHOMAP(1)
        !
        CALL SPLCY(THETAHAT,SIGMAP(1,KPSI),NT1,RC2PI, &
             &              ZD2SIG,ZA1,ZB1,ZC1)
        CALL SPLCYP(THETAHAT,ZTET,NT1,RC2PI,RC2PI, &
             &               ZD2TET,ZA1,ZB1,ZC1)
        CALL SPLCY(THETAHAT,RHOMAP,NT1,RC2PI, &
             &              ZD2RHO,ZA1,ZB1,ZC1)
        !
        ZD2SIG(NT2) = ZD2SIG(1) 
        ZD2TET(NT2) = ZD2TET(1) 
        ZD2RHO(NT2) = ZD2RHO(1) 
        !
        CALL RESETI(IC,NCHI,1)
        DO JG=1,NCHI
           DO JT = 1,NT2
              IF (IC(JG).EQ.1) THEN
                 IT0(JG) = JT-1
                 IF (THETAHAT(JT).GE.CHIM(JG)) IC(JG)  = 0
              ENDIF
           ENDDO
        ENDDO
        !
        DO J=1,NCHI
           !
           ICHIM = IT0(J)
           !
           IF (ICHIM .LT. 1)   ICHIM = 1
           IF (ICHIM .GT. NT1) ICHIM = NT1
           !
           ZH = THETAHAT(ICHIM+1) - THETAHAT(ICHIM)
           ZA = (THETAHAT(ICHIM+1) - CHIM(J)) / ZH
           ZB = (CHIM(J) - THETAHAT(ICHIM)) / ZH
           ZC = (ZA + 1) * (ZA - 1) * ZH * &
                &        (THETAHAT(ICHIM+1) - CHIM(J)) / 6._RKIND
           ZD = (ZB + 1) * (ZB - 1) * ZH * &
                &        (CHIM(J) - THETAHAT(ICHIM)) / 6._RKIND
           ! 
           RHOPARXT(J,KPSI,JPOP) = ZA*RHOMAP(ICHIM) + ZB*RHOMAP(ICHIM+1) + &
                &                  ZC*ZD2RHO(ICHIM) + ZD*ZD2RHO(ICHIM+1)
           TETPARXT(J,KPSI,JPOP) = ZA*ZTET(ICHIM)   + ZB*ZTET(ICHIM+1) + &
                &                  ZC*ZD2TET(ICHIM) + ZD*ZD2TET(ICHIM+1)
           !
           IF (TETPARXT(J,KPSI,JPOP) .LT. CT(1)) &
                &                   TETPARXT(J,KPSI,JPOP) = TETPARXT(J,KPSI,JPOP) + 2._RKIND*CPI
           IF (TETPARXT(J,KPSI,JPOP) .GT. CT(NT1)) &
                &                   TETPARXT(J,KPSI,JPOP) = TETPARXT(J,KPSI,JPOP) - 2._RKIND*CPI
           !
           IF (KPSI .EQ. NPMAX) THEN
              !
              SIGPARXT(J,KPSI,JPOP) = 1._RKIND
              !
           ELSE
              !
              SIGPARXT(J,KPSI,JPOP) = ZA*SIGMAP(ICHIM,KPSI) + ZB*SIGMAP(ICHIM+1,KPSI) + &
                   &                  ZC*ZD2SIG(ICHIM)      + ZD*ZD2SIG(ICHIM+1)
              !
           ENDIF
        END DO
     ENDDO
  ENDDO
  !
  GAMMAPARXT=0._RKIND
  DO JPOP=1,NPOPULATIONS
     DO KPSI=1,NPMAX
        !
        CALL BOUND(NCHI,TETPARXT(1,KPSI,JPOP),ZBNDPARXT)
        !
        DO J=1,NCHI
           ZRHO   = SIGPARXT(J,KPSI,JPOP)*ZBNDPARXT(J)
           RPARXT(J) = ZRHO * COS(TETPARXT(J,KPSI,JPOP)) + R0
           ZPARXT(J) = ZRHO * SIN(TETPARXT(J,KPSI,JPOP)) + RZ0
        ENDDO
        DO J=1,NCHI
           IF (KPSI==1) THEN
              GAMMAPARXT(J,1,JPOP)=RHOPARXT(J,KPSI,JPOP)*PSIISO(KPSI)
           ELSE
              GAMMAPARXT(J,KPSI,JPOP)=GAMMAPARXT(J,KPSI-1,JPOP)+ &
                                      RHOPARXT(J,KPSI,JPOP)*(PSIISO(KPSI) - PSIISO(KPSI-1))
           ENDIF
        ENDDO
        !
        WRITE(NXTOR) PSMISO(KPSI),PSIISO(KPSI)
        WRITE(NXTOR) (RPARXT(J),J=1,NCHI)
        WRITE(NXTOR) (ZPARXT(J),J=1,NCHI)
        WRITE(NXTOR) (CHIM(J),J=1,NCHI)
        WRITE(NXTOR) (GAMMAPARXT(J,KPSI,JPOP),J=1,NCHI)
        !
     ENDDO
  ENDDO
  !
  RETURN
END SUBROUTINE XTINIT
function equilibriumdensity(psi,bnorm,kpop)
  !-----------------------------------------------------------------------------
  !
  !     Particle equilibrium distribution function (SI)
  !
  ! Currently, density & temperature are centered-Maxwellian.
  ! Energy distribution is Maxwellian or Slowing-Down.
  ! Anisotropy is Maxwellian around a given value lambda0.
  ! Convention : if sd_vmax <= 0.0, Energy is maxwellian. Else : Slowing-down.
  !              if dlambda >= 100.0, energy isotropy is assumed.
  !
  use globals
  !
  real(rkind), parameter :: accuracy = 1e-5 ! arbitrary for the moment.
  real(rkind) :: psi, bnorm, equilibriumdensity
  real(rkind) :: isotropicequilibriumdensity,anisotropyfactor
  integer     :: kpop
  !
  equilibriumdensity = isotropicequilibriumdensity(psi,kpop)   &
       *anisotropyfactor(1.0_rkind,accuracy,bnorm,kpop)
  !       
  return
end function equilibriumdensity
!-----------------------------------------------------------------------------
function isotropicequilibriumdensity(psi,kpop)
  !-----------------------------------------------------------------------------
  use globals
  !
  implicit none
  real(rkind) :: psi, isotropicequilibriumdensity
  integer     :: kpop

  integer, parameter :: distribtype = 2
  integer :: i,ISRCHFGT
  real(rkind) :: t
  
  real(rkind) :: psimax
  REAL(RKIND) :: PP,ZS1
  
  if (distribtype.eq.1) then
     !
     ! 1st type of distribution
     !
     psimax = 4.5_rkind*dpsi_ni(kpop)*cpsrf
     if (psi > psimax) then
        isotropicequilibriumdensity = 0
     else
        isotropicequilibriumdensity = (1-(psi/psimax)**2)**2
     end if
     
  elseif(distribtype.eq.2) then
     !
     ! 2nd type of distribution
     !
     psimax = 4.5_rkind*dpsi_ni(kpop)*cpsrf
     if (psi>psimax) then
        isotropicequilibriumdensity = 0.0_rkind
     else
        isotropicequilibriumdensity = exp(-psi/dpsi_ni(kpop)/cpsrf)
     end if
     
  elseif(distribtype.eq.3) then
     !
     ! 3rd type of distribution
     !
     if ((nppfun.ne.7).and.(npp.eq.1)) then
        print *, 'In isotropicequilibriumdensity, nppfun must be 7.'
        print *, '                                npp    must be 1.'
        stop 
     end if

     isotropicequilibriumdensity = 1.0_rkind
     PP = (psi - cpsrf)*(-SPSIM/cpsrf)
     ZS1 = 1._RKIND - PP / SPSIM
     if (ZS1>=1) then 
        isotropicequilibriumdensity = 0.0_rkind
     else  
        isotropicequilibriumdensity = (1-ZS1**AP(3))**AP(2)
     end if

  else
     !
     ! ATTENTION : check cpr(0)
     !

     psimax = 0.95*cpsrf

     if (psi>psimax) then
        isotropicequilibriumdensity = 0._rkind
     else        
        i = ISRCHFGT(NISO1EFF,psiiso,1,psi) - 1
        t = (sqrt(psi) - sqrt(psiiso(i)))/(sqrt(psiiso(i+1)) - sqrt(psiiso(i)))
        isotropicequilibriumdensity = cpr(i)*(1._rkind-t) + cpr(i+1)*t
     end if
  end if
  
end function isotropicequilibriumdensity
!---------------------------------------------------------------------------------
function anisotropyfactor(x,accuracy,bnorm,kpop)
  !---------------------------------------------------------------------------------
  ! If dlambda -> infty : - we recover the isotropic case
  !                       - the function returns 1.
  use globals         
  implicit none
  real(rkind) :: x, anisotropyfactor
  real(rkind) :: bnorm
  integer     :: kpop
  !
  real(rkind) :: accuracy
  real(rkind) :: intervalsize
  integer     :: nintervals
  integer     :: i
  !
  if (dlambda(kpop) > 100._rkind) then       ! isotropic
     anisotropyfactor = 1.0_rkind
  else 
     intervalsize = min(accuracy, 2.0_rkind*bnorm*accuracy*dlambda(kpop)/lambda0(kpop))
     nintervals   = int(x/intervalsize) + 1
     intervalsize = x/real(nintervals)
     !
     anisotropyfactor = 0.0_rkind
     do i=0,nintervals-1
        anisotropyfactor = anisotropyfactor                       &
             + exp(-((1._rkind-(intervalsize*real(i))**2-lambda0(kpop)) &
             /(bnorm*dlambda(kpop)))**2 )
        ! bnorm and dlambda not equal to zero
     end do
     anisotropyfactor = anisotropyfactor*intervalsize
  end if
  !
  return
end function anisotropyfactor
!---------------------------------------------------------------------------------
