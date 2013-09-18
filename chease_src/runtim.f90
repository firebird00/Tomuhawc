SUBROUTINE runtim
  !
  !   Print elapsed time
  !
  USE prec_const
  !
  REAL, SAVE :: stime
  INTEGER, SAVE :: icall = 0
  REAL :: eltime
  INTEGER :: c,t,ic_max
  !----------------------------------------------------------------------
  CALL system_clock(c,t,ic_max)
  IF( t .EQ. 0 ) RETURN   ! No clock on this system
  !
  eltime = REAL(c) / REAL(t)
  IF( icall .EQ. 0 ) THEN  ! Initial call
     stime = eltime
     print *,'clock step= ',1._rkind/real(t,rkind),' max cpu before restart clock= ', &
          & real(ic_max,rkind)/real(t,rkind)/3600._rkind,' hours'
     CALL system_CLOCK(c,t)
     eltime = REAL(c) / REAL(t)
     icall = 1
     return
  END IF
  !    if has gone once around clock adapt stime
  IF (eltime .LE. stime) stime = stime - REAL(ic_max)/REAL(t)
  WRITE(*,'(/a,1pe14.6,a)') 'CPU TIME USED SO FAR =',eltime-stime,' SECS'
  !
END SUBROUTINE runtim
