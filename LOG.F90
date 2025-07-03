
	  !********************************************************
	  !*                                                      *
	  !********************************************************
	  !************************* MEGA *************************
	  !********************************************************
	  !*                                                      *
	  !********* Meshfree Explicit Galerkin Analysis **********
	  !*                                                      *
	  !*                                                      *
	  !*                                                      *
	  !*     Copyright 2016 Michael C Hillman                 *
	  !*                                                      *
	  !********************************************************
	  
     SUBROUTINE LOG_HEADER
     OPEN(50,FILE='MEGA.log',STATUS='REPLACE')
     WRITE(50,*)    '********************************************************'
     WRITE(50,*)    '*                                                      *'
     WRITE(50,*)    '********************************************************'
     WRITE(50,*)    '************************* MEGA *************************'
     WRITE(50,*)    '********************************************************'
     WRITE(50,*)    '*                                                      *'
     WRITE(50,*)    '********* Meshfree Explicit Galerkin Analysis **********'
     WRITE(50,*)    '*                                                      *'
     WRITE(50,*)    '*                                                      *'
     WRITE(50,*)    '*                                                      *'
     WRITE(50,*)    '*     Copyright 2016 Michael C Hillman                 *'
     WRITE(50,*)    '*                                                      *'
     WRITE(50,*)    '********************************************************'
     CLOSE(50)
     OPEN(50,FILE='WARN.log',STATUS='REPLACE')
     WRITE(50,*)    '********************************************************'
     WRITE(50,*)    '********************************************************'
     WRITE(50,*)    '************* MEGA WARNING LOG *************************'
     WRITE(50,*)    '********************************************************'
     WRITE(50,*)    '********************************************************'
     CLOSE(50)
     RETURN
	 END SUBROUTINE
     
     
     
     
      SUBROUTINE WARN(STRING_LINE)
<<<<<<< HEAD


=======
>>>>>>> a5f824d (問題 2：子程式需要 ACC ROUTINE 指令)
	  !
	  ! FUNCTION OF THIS SUBROUTINE:
	  !
	  ! APPEND THE LOG FILE WITH A STRING
	  !
	  IMPLICIT NONE
	  CHARACTER(*):: STRING_LINE
      OPEN(50,FILE='MEGA.log',STATUS='OLD',POSITION='APPEND')
      WRITE(50,*) '********************************************'
      WRITE(50,*) '************    WARNING   ******************'
      WRITE(50,*) '********************************************'
      WRITE(50,*) 
      WRITE(50,*) STRING_LINE
      WRITE(50,*) 
      WRITE(50,*) '*******************************************'
      CLOSE(50)
      
      OPEN(50,FILE='WARN.log',STATUS='OLD',POSITION='APPEND')
      WRITE(50,*) '********************************************'
      WRITE(50,*) '************    WARNING   ******************'
      WRITE(50,*) '********************************************'
      WRITE(50,*) 
      WRITE(50,*) STRING_LINE
      WRITE(50,*) 
      WRITE(50,*) '*******************************************'
      CLOSE(50)
	  
      WRITE(*,*) '********************************************'
      WRITE(*,*) '************    WARNING   ******************'
      WRITE(*,*) '********************************************'
      WRITE(*,*) 
      WRITE(*,*) STRING_LINE
      WRITE(*,*) 
      WRITE(*,*) '*******************************************'
      
      RETURN
	  END SUBROUTINE
      
      
      
     SUBROUTINE STRONG_WARN(STRING_LINE)
	  !
	  ! FUNCTION OF THIS SUBROUTINE:
	  !
	  ! APPEND THE LOG FILE WITH A STRING
	  !
	  IMPLICIT NONE
	  CHARACTER(*):: STRING_LINE
      OPEN(50,FILE='MEGA.log',STATUS='OLD',POSITION='APPEND')
      WRITE(50,*) 
      WRITE(50,*) '*********************************************************'
      WRITE(50,*) '************    !!!STRONG WARNING!!!   ******************'
      WRITE(50,*) '*********************************************************'
      WRITE(50,*) 
      WRITE(50,*) STRING_LINE
      WRITE(50,*) 
      WRITE(50,*) '********************************************************'
      WRITE(50,*) '* YOU ARE STRONGLY ENCOURAGED TO CHECK YOUR INPUT FILE *'
      WRITE(50,*) '********************************************************'
      WRITE(50,*) 
      CLOSE(50)
      
      WRITE(*,*) 
      WRITE(*,*) '*********************************************************'
      WRITE(*,*) '************    !!!STRONG WARNING!!!   ******************'
      WRITE(*,*) '*********************************************************'
      WRITE(*,*) 
      WRITE(*,*) STRING_LINE
      WRITE(*,*) 
      WRITE(*,*) '********************************************************'
      WRITE(*,*) '* YOU ARE STRONGLY ENCOURAGED TO CHECK YOUR INPUT FILE *'
      WRITE(*,*) '********************************************************'
      WRITE(*,*) 
      
      RETURN
	  END SUBROUTINE
      
     SUBROUTINE LOG_APPEND_SPACE(STRING_LINE)
	  !
	  ! FUNCTION OF THIS SUBROUTINE:
	  !
	  ! APPEND THE LOG FILE WITH A STRING WITH SPACES ABOVE AND BELOW
	  !
	  IMPLICIT NONE
	  CHARACTER(*):: STRING_LINE
      OPEN(50,FILE='MEGA.log',STATUS='OLD',POSITION='APPEND')
      WRITE(50,*) 
      WRITE(50,*) STRING_LINE
      WRITE(50,*) 
      CLOSE(50)
      RETURN
	  END SUBROUTINE
      
      
     SUBROUTINE WRITE_OUT(STRING_LINE)
	  !
	  ! FUNCTION OF THIS SUBROUTINE:
	  !
	  ! APPEND THE LOG FILE WITH A STRING WITH SPACES ABOVE AND BELOW,
      ! AND ALSO WRITE THE SAME TO SCREEN
	  !
	  IMPLICIT NONE
	  CHARACTER(*):: STRING_LINE
      OPEN(50,FILE='MEGA.log',STATUS='OLD',POSITION='APPEND')
      WRITE(50,*) 
      WRITE(50,*) STRING_LINE
      WRITE(50,*) 
      WRITE(*,*) 
      WRITE(*,*) STRING_LINE
      WRITE(*,*) 
      CLOSE(50)
      RETURN
	  END SUBROUTINE
      
      
      
     SUBROUTINE LOG_APPEND(STRING_LINE)
	  !
	  ! FUNCTION OF THIS SUBROUTINE:
	  !
	  ! APPEND THE LOG FILE WITH A STRING
	  !
	  IMPLICIT NONE
	  CHARACTER(*):: STRING_LINE
      OPEN(50,FILE='MEGA.log',STATUS='OLD',POSITION='APPEND')
      WRITE(50,*) STRING_LINE
      CLOSE(50)
      RETURN
	  END SUBROUTINE