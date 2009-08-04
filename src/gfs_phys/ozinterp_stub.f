!**********************************************************************
!
      SUBROUTINE ozinterpol(me,latd,nlats,IDATE,FHOUR &
      ,jindx1,jindx2,ozplin,ozplout,ddy)
!
      implicit none
      integer             latd
!
 
      integer  JINDX1(LATD), JINDX2(LATD)
      integer  me,idate(4),nlats
      integer  IDAT(8),JDAT(8)
!
      real(kind=8) ozplin(1,1,1,1),fhour
      real(kind=8) DDY(1)
      real(kind=8) ozplout(1,1,1)
!
!
      RETURN
      END
