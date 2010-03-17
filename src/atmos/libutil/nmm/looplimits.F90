!===============================================================
!  From Jim Abeles.
!  Block partition the loop bounds (lb...ub) -> (i1...i2).
!  The number of tasks is ntasks;  taskid = 0, 1, ..., ntasks-1.
!  The first nt1 tasks get a chunk one bigger than the rest.
!  The counts and displacements arrays range from 1 to ntasks.
!===============================================================
!
      subroutine looplimits(taskid, ntasks, lb, ub, i1, i2)
      implicit none
      integer taskid, ntasks, lb, ub, i1, i2
      integer chunk, nwork, nt1, nt2
      integer itask, netdisp
      integer counts(ntasks), displacements(ntasks)

      nwork = ub - lb + 1
      chunk = nwork/ntasks
      nt1 = nwork - ntasks*chunk
      nt2 = ntasks - nt1

      netdisp = lb
      do itask = 1, nt1
         counts(itask) = chunk + 1
         displacements(itask) = netdisp  
         netdisp = min(ub,netdisp+chunk+1)
      end do
      do itask = nt1 + 1 , ntasks
         counts(itask) = chunk
         displacements(itask) = netdisp  
         netdisp = min(ub,netdisp+chunk)
      end do

      i1 = displacements(taskid+1)
      i2 = min(ub,i1+counts(taskid+1)-1)

      return
      end
