!*DECK C2SE02
!*CALL PROCESS
SUBROUTINE DIRECT
  !        #################
  !
  !                                        AUTHORS:
  !                                        H. LUTJENS,  CRPP-EPFL
  !                                        A. BONDESON, CRPP-EPFL
  !
  !**********************************************************************
  !                                                                     *
  ! C2SE02 GAUSS SEIDEL ELIMINATION                                     *
  !                                                                     *
  !**********************************************************************
  !
  USE globals
  IMPLICIT NONE
  !
  !**********************************************************************
  !                                                                     *
  ! 1. SOLVE L * Y = B                                                  *
  !                                                                     *
  !**********************************************************************
  !
  CALL LYV(A,B,N4NSNT,NP4NST,NBAND,NPBAND)
  !
  !**********************************************************************
  !                                                                     *
  ! 2. SOLVE D * W = Y                                                  *
  !                                                                     *
  !**********************************************************************
  !
  CALL DWY(A,B,N4NSNT,NP4NST,NBAND,NPBAND)
  !
  !**********************************************************************
  !                                                                     *
  ! 3. SOLVE LT * X = W                                                 *
  !                                                                     *
  !**********************************************************************
  !
  CALL LTXW(A,B,N4NSNT,NP4NST,NBAND,NPBAND)
  !
  !**********************************************************************
  !                                                                     *
  ! 4. RESULT IN B                                                      *
  !                                                                     *
  !**********************************************************************
  !
  RETURN
END SUBROUTINE DIRECT