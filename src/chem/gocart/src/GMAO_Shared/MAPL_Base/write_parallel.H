! $Id: write_parallel.H,v 1.2 2007/12/19 14:17:25 trayanov Exp $

#ifdef NAME_
#undef NAME_
#endif

#define NAME_ WRITE_PARALLEL_

#include "overload.macro"


subroutine SUB_ ( data, UNIT, format, RC)

#if (VARTYPE_ > 0)
  TYPE_(kind=EKIND_ ), intent(in   )            :: data DIMENSIONS_
#else
  character(LEN=*),    intent(in   )            :: data DIMENSIONS_
#endif
  integer,             intent(in   ),  optional :: UNIT
  character(len=*),    intent(in   ),  optional :: format
  integer         ,    intent(  out),  optional :: RC

  character(len=ESMF_MAXSTR) :: FORMATTED
  character(len=ESMF_MAXSTR) :: IAM='WRITE_PARALLEL'
  integer :: status

  if (present(unit)) then
     if (unit == UNDEF) then
        RETURN_(ESMF_SUCCESS)
     endif
     if(unit < 0) then

#if (RANK_ > 2 || VARTYPE_ <= 0)

        ASSERT_(.FALSE.)

#else

        ASSERT_(-UNIT<=LAST_UNIT)
        munit => MEM_units(-unit)
        munit%prevrec = munit%prevrec + 1

        if(.not.associated(munit%Records)) then
           allocate(Rec(16),stat=status)
           munit%Records => REC
           VERIFY_(STATUS)
        elseif(size(munit%Records)< munit%prevrec) then
           allocate(REC(munit%prevrec*2),stat=status)
           VERIFY_(STATUS)
           REC(:munit%prevrec-1) = munit%Records
           deallocate(munit%Records)
           munit%Records => REC
        endif

#if (RANK_ == 0)
        call dealloc_(munit%Records(munit%prevrec),rc=status)	
        VERIFY_(STATUS)
#endif

#if (RANK_ == 1)
        call alloc_(munit%Records(munit%prevrec),TKR_,size(data,1),rc=status)	
        VERIFY_(STATUS)
#endif

#if (RANK_ == 2)
        call alloc_(munit%Records(munit%prevrec),TKR_,size(data,1),size(data,2),rc=status)	
        VERIFY_(STATUS)
#endif

        munit%Records(munit%prevrec)%TKR_  = data

#endif
     else ! unit is > 0
        if (MAPL_AM_I_ROOT()) then
           inquire(unit=UNIT, formatted=FORMATTED)
           if   (FORMATTED == "YES") then
              if (present(format)) then;    write(UNIT, format) data
              else;    write(UNIT, *     ) data
              end if
           elseif(FORMATTED == "NO") then; write(UNIT        ) data
           end if
        endif
     end if
  else

     if (MAPL_AM_I_ROOT()) then
        if (present(format)) then; write(*, format) data
        else; write(*,      *) data
        end if
     end if
  end if

  RETURN_(ESMF_SUCCESS)

end subroutine SUB_

!---------------------------
#undef NAME_

#undef DIMENSIONS_
#undef RANK_
#undef VARTYPE_
