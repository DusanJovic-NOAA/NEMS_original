!$$$  Module Documentation Block
!
! Module: mersenne_twister   Modern random number generator
!   Prgmmr: Iredell          Org: W/NX23     date: 2005-06-14
!
! Abstract: This module calculates random numbers using the Mersenne twister.
!   (It has been adapted to a Fortran 90 module from open source software.
!   The comments from the original software are given below in the remarks.)
!   The Mersenne twister (aka MT19937) is a state-of-the-art random number
!   generator based on Mersenne primes and originally developed in 1997 by
!   Matsumoto and Nishimura. It has a period before repeating of 2^19937-1,
!   which certainly should be good enough for geophysical purposes. :-)
!   Considering the algorithm's robustness, it runs fairly speedily.
!   (Some timing statistics are given below in the remarks.)
!   This adaptation uses the standard Fortran 90 random number interface,
!   which can generate an arbitrary number of random numbers at one time.
!   The random numbers generated are uniformly distributed between 0 and 1.
!   The module also can generate random numbers from a Gaussian distribution
!   with mean 0 and standard deviation 1, using a Numerical Recipes algorithm.
!   The module also can generate uniformly random integer indices.
!   There are also thread-safe versions of the generators in this adaptation,
!   necessitating the passing of generator states which must be kept private.
!
! Program History Log:
!   2005-06-14  Mark Iredell
!
! Usage: 
!   The module can be compiled with 4-byte reals or with 8-byte reals, but
!   4-byte integers are required. The module should be endian-independent.
!   The Fortran 90 interfaces random_seed and random_number are overloaded
!   and can be used as in the standard by adding the appropriate use statement
!     use mersenne_twister
!   In the below use cases, harvest is a real array of arbitrary size,
!   and iharvest is an integer array of arbitrary size.
!   To generate uniformly distributed random numbers between 0 and 1,
!     call random_number(harvest)
!   To generate Gaussian distributed random numbers with 0 mean and 1 sigma,
!     call random_gauss(harvest)
!   To generate uniformly distributed random integer indices between 0 and n,
!     call random_index(n,iharvest)
!   In standard "saved" mode, the random number generator can be used without
!   setting a seed. But to set a seed, only 1 non-zero integer is required, e.g.
!     call random_setseed(4357) ! set default seed
!   The full generator state can be set via the standard interface random_seed,
!   but it is recommended to use this method only to restore saved states, e.g.
!     call random_seed(size=lsave)  ! get size of generator state seed array
!     allocate isave(lsave)         ! allocate seed array
!     call random_seed(get=isave)   ! fill seed array (then maybe save to disk)
!     call random_seed(put=isave)   ! restore state (after read from disk maybe)
!   Locally kept generator states can also be saved in a seed array, e.g.
!     type(random_stat):: stat
!     call random_seed(get=isave,stat=stat)  ! fill seed array
!     call random_seed(put=isave,stat=stat)  ! restore state
!   To generate random numbers in a threaded region, the "thread-safe" mode
!   must be used where generator states of type random_state are passed, e.g.
!     type(random_stat):: stat(8)
!     do i=1,8                               ! threadable loop
!       call random_setseed(7171*i,stat(i))  ! thread-safe call
!     enddo
!     do i=1,8                               ! threadable loop
!       call random_number(harvest,stat(i))  ! thread-safe call
!     enddo
!     do i=1,8                               ! threadable loop
!       call random_gauss(harvest,stat(i))   ! thread-safe call
!     enddo
!     do i=1,8                               ! threadable loop
!       call random_index(n,iharvest,stat(i))! thread-safe call
!     enddo
!   There is also a relatively inefficient "interactive" mode available, where
!   setting seeds and generating random numbers are done in the same call.
!   There is also a functional mode available, returning one value at a time.
!   
! Public Defined Types:
!   random_stat       Generator state (private contents)
!   
! Public Subprograms:
!   random_seed         determine size or put or get state
!     size              optional integer output size of seed array
!     put               optional integer(:) input seed array
!     get               optional integer(:) output seed array
!     stat              optional type(random_stat) (thread-safe mode)
!   random_setseed      set seed (thread-safe mode)
!     inseed            integer seed input
!     stat              type(random_stat) output
!   random_setseed      set seed (saved mode)
!     inseed            integer seed input
!   random_number       get mersenne twister random numbers (thread-safe mode)
!     harvest           real(:) numbers output
!     stat              type(random_stat) input
!   random_number       get mersenne twister random numbers (saved mode)
!     harvest           real(:) numbers output
!   random_number       get mersenne twister random numbers (interactive mode)
!     harvest           real(:) numbers output
!     inseed            integer seed input
!   random_number_f     get mersenne twister random number (functional mode)
!     harvest           real number output
!   random_gauss        get gaussian random numbers (thread-safe mode)
!     harvest           real(:) numbers output
!     stat              type(random_stat) input
!   random_gauss        get gaussian random numbers (saved mode)
!     harvest           real(:) numbers output
!   random_gauss        get gaussian random numbers (interactive mode)
!     harvest           real(:) numbers output
!     inseed            integer seed input
!   random_gauss_f      get gaussian random number (functional mode)
!     harvest           real number output
!   random_index        get random indices (thread-safe mode)
!     imax              integer maximum index input
!     iharvest          integer(:) numbers output
!     stat              type(random_stat) input
!   random_index        get random indices (saved mode)
!     imax              integer maximum index input
!     iharvest          integer(:) numbers output
!   random_index        get random indices (interactive mode)
!     imax              integer maximum index input
!     iharvest          integer(:) numbers output
!     inseed            integer seed input
!   random_index_f      get random index (functional mode)
!     imax              integer maximum index input
!     iharvest          integer number output
!
! Remarks:
!   (1) Here are the comments in the original open source code:
!     A C-program for MT19937: Real number version
!     genrand() generates one pseudorandom real number (double)
!     which is uniformly distributed on [0,1]-interval, for each
!     call. sgenrand(seed) set initial values to the working area
!     of 624 words. Before genrand(), sgenrand(seed) must be
!     called once. (seed is any 32-bit integer except for 0).
!     Integer generator is obtained by modifying two lines.
!     Coded by Takuji Nishimura, considering the suggestions by
!     Topher Cooper and Marc Rieffel in July-Aug. 1997.
!     This library is free software; you can redistribute it and/or
!     modify it under the terms of the GNU Library General Public
!     License as published by the Free Software Foundation; either
!     version 2 of the License, or (at your option) any later
!     version.
!     This library is distributed in the hope that it will be useful,
!     but WITHOUT ANY WARRANTY; without even the implied warranty of
!     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
!     See the GNU Library General Public License for more details.
!     You should have received a copy of the GNU Library General
!     Public License along with this library; if not, write to the
!     Free Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA
!     02111-1307  USA
!     Copyright (C) 1997 Makoto Matsumoto and Takuji Nishimura.
!     When you use this, send an email to: matumoto@math.keio.ac.jp
!     with an appropriate reference to your work.
!     Fortran translation by Hiroshi Takano.  Jan. 13, 1999.
!
!   (2) On a single IBM Power4 processor on the NCEP operational cluster (2005)
!     each Mersenne twister random number takes less than 30 ns, about 3 times
!     slower than the default random number generator, and each random number
!     from a Gaussian distribution takes less than 150 ns.
!     
! Attributes:
!   Language: Fortran 90
!
!$$$
      module mersenne_twister
!        private
        public random_stat
	public random_number
        public random_setseed
	integer,parameter:: ttt=624
        type random_stat
          private
          integer:: mti=ttt+1
          integer:: mt(0:ttt-1)
          integer:: iset
          real:: gset
        end type
        interface random_setseed
          module procedure random_setseed_s
          module procedure random_setseed_t
        end interface
        interface random_number
          module procedure random_number_i
          module procedure random_number_s
          module procedure random_number_t
        end interface
        contains
        subroutine random_number_i(harvest,inseed)
          implicit none
          real,intent(out):: harvest(:)
          integer,intent(in):: inseed
        end subroutine
!  Subprogram random_number_s
!  Generates random numbers in saved mode; overloads Fortran 90 standard.
        subroutine random_number_s(harvest)
          implicit none
          real,intent(out):: harvest(:)
        end subroutine
!  Subprogram random_number_t
!  Generates random numbers in thread-safe mode.
        subroutine random_number_t(harvest,stat)
          implicit none
          real,intent(out):: harvest(:)
          type(random_stat),intent(inout):: stat
        end subroutine
        subroutine random_setseed_s(inseed)
          implicit none
          integer,intent(in):: inseed
        end subroutine
!  Subprogram random_setseed_t
!  Sets seed in thread-safe mode.
        subroutine random_setseed_t(inseed,stat)
          implicit none
          integer,intent(in):: inseed
          type(random_stat),intent(out):: stat
        end subroutine
!  All the subprograms
!  Subprogram random_seed
!  Sets and gets state; overloads Fortran 90 standard.
      end module mersenne_twister
