! wrapper is supposed to stop test.f getting too messy, but will be
! a slight overhead.

module rng_wrapper_mod

  use precision_mod, only : ik4,rk8,ik8
  use, intrinsic :: iso_c_binding

  implicit none

  ! num_vals is number of calls to rng routines - so vals produced is (eg threefry4x32) four times this - yes i should rename it
  ! so if using threefry4x32 - each call produces four values. so if vector length is 8 - actually get 32 values
  ! essentially vector length - or more - can do as many as you want. multiple of vector length is best.
  ! important: integer vector lengths -  sse/avx=4 avx2=8 avx512=16
  ! it is not to hard to set up to get this value from random123 - and for many systems from preproc defines - eg. __avx512f__

  !set this to your system
  integer(kind=ik4),parameter   ::  num_vals  = 16 ! vector length (or multiple of)
  

  integer(kind=ik4),parameter   ::  buffersize_single = 4_ik4    !eg. for 4x32 this should be 4
  !integer(kind=ik4),parameter   ::  buffersize_single = 2_ik4

  !pain setting here and in random123
  !integer(kind=ik4),parameter   ::  num_vals  = 4 ! 4/8/16 should be num. calls as vals is four times this.
  integer(kind=ik4),parameter   ::  buffersize_multi  = num_vals*buffersize_single ! 12 = 3 calls to 4x64... 16 = 4 calls

  integer(kind=ik4),parameter   ::  rng_seed= 12345_ik4
  integer(kind=ik4),parameter   ::  rng_key = 54321_ik4

  character(*),parameter        ::  threefry_name = "threefry4x32"
  !character(*),parameter        ::  philox_name   = "philox4x32"

  !character(*),parameter        ::  threefry_name = "threefry4x64"
  !character(*),parameter        ::  philox_name   = "philox4x64"

  !character(*),parameter        ::  threefry_name = "threefry2x64"
  !character(*),parameter        ::  philox_name   = "philox2x64"

  !externals - quick way - before doing interfaces
  integer(kind=ik4)             ::  threefry4x64f
  integer(kind=ik4)             ::  threefry2x64f
  integer(kind=ik4)             ::  philox4x64f
  integer(kind=ik4)             ::  philox4x32f

  external                          threefry4x64f
  external                          threefry2x64f
  external                          philox4x64f
  external                          philox4x32f

    
  ! interface block
  interface
    
    ! Utility routine - currently just prints rather than returning.
    integer(c_int) function v123_get_vec_sizes() bind(c,name='v123_get_vec_sizes')
      import :: c_int
    end function v123_get_vec_sizes   
              
    ! scalar variants

    !original threefry - should be built -no-vec as faster
    integer(c_int) function threefry4x32f(ctr4, key4, buff) bind(c,name='threefry4x32f')
      import :: c_int, c_double, num_vals
      integer(kind=c_int),dimension(4) ,intent(in)             :: ctr4
      integer(kind=c_int),dimension(4) ,intent(in)             :: key4        
      real(kind=c_double),dimension(4) ,intent(inout)   :: buff
      
      !note - alt. can just use an asterisk for size if dont know here - only passing pointer - works for me
      !real(kind=c_double),dimension(*) ,intent(inout)          :: buff
    end function threefry4x32f    

    !scalar - built with vec
    integer(c_int) function threefry4x32f_scalar(ctr4, key4, buff) bind(c,name='threefry4x32f_scalar')
      import :: c_int, c_double, num_vals
      integer(kind=c_int),dimension(4) ,intent(in)             :: ctr4
      integer(kind=c_int),dimension(4) ,intent(in)             :: key4        
      real(kind=c_double),dimension(4) ,intent(inout)   :: buff
    end function threefry4x32f_scalar 

    !scalar - built with novec
    integer(c_int) function threefry4x32f_novec(ctr4, key4, buff) bind(c,name='threefry4x32f_novec')
      import :: c_int, c_double, num_vals
      integer(kind=c_int),dimension(4) ,intent(in)             :: ctr4
      integer(kind=c_int),dimension(4) ,intent(in)             :: key4        
      real(kind=c_double),dimension(4) ,intent(inout)   :: buff
    end function threefry4x32f_novec 


    !vector vairants
    
    !sh - working notes on alternative routines and performance.
    !sh - need testing to see if performance can be gained from 1 and 2 due to simpler interface/calculation. 
    !sh - note that one consideration is the number of arguments - may go into stack passing - one possibility
    !     is to return a single buffer of random numbers rather than four.
    !     register pressure in general is also a consideration.

    ! alternative routines
    !  1. threefry4x32f_multi_ss                7  arguments 
    !  2. threefry4x32f_multi_ctrkey1           8  arguments 
    !  3. threefry4x32f_multi_ctrkey_all       12 arguments
    
    ! all arguments are pointers (pass by reference) *update added increment and in_num_vals to some - pass by value
    
    ! 1. single stream version/s: 
    ! provide one ctr-set and one key-set like scalar version. ctr1 is incremented across vector dimesion with supplied increment.
    ! an increment of 1 is as good as any for quality of random numbers
    ! values 2/3/4 of each ctr/key set are replicated across vector dimension.
    ! Four vectors of random number are returned.

    !fixed vector length version (#define set in c file)
    integer(c_int) function threefry4x32f_multi_ss_fix(ctr4, key4, increment, buff) bind(c,name='threefry4x32f_multi_ss_fix')
      import :: c_int, c_double, num_vals
      integer(kind=c_int),dimension(4) ,intent(in)             :: ctr4
      integer(kind=c_int),dimension(4) ,intent(in)             :: key4        
      integer(kind=c_int),value                                :: increment        
      real(kind=c_double),dimension(num_vals) ,intent(inout)   :: buff
      
      !note - can just use an asterisk for size - only passing pointer
      !real(kind=c_double),dimension(*) ,intent(inout)   :: buff
      
      !integer(kind=c_int),value                                :: in_num_vals
    end function threefry4x32f_multi_ss_fix    
    
    !provided length version (in_num_vals)
    integer(c_int) function threefry4x32f_multi_ss(ctr4, key4, increment, buff, in_num_vals) bind(c,name='threefry4x32f_multi_ss')
      import :: c_int, c_double, num_vals
      integer(kind=c_int),dimension(4) ,intent(in)             :: ctr4
      integer(kind=c_int),dimension(4) ,intent(in)             :: key4        
      integer(kind=c_int),value                                :: increment        
      real(kind=c_double),dimension(num_vals) ,intent(inout)   :: buff
      integer(kind=c_int),value                                :: in_num_vals
    end function threefry4x32f_multi_ss    


    ! 2. single counters/key vector version
    !    supply vector of values for ctr1 and key1 - a default ctr/key is used for ctr/key 2/3/4
    ! four vectors of random number are returned.      .
    !integer function threefry4x32f_multi_ctrkey1(ctr1, key1, def_ctr, def_key, ran1, ran2, ran3, ran4) bind(c,name='threefry4x32f_multi_ctrkey1')
    !  import :: c_int, c_double
    !  integer(kind=c_int),dimension(*) ,intent(in)    :: ctr1
    !  integer(kind=c_int),dimension(*) ,intent(in)    :: key1
    !  integer(kind=c_int)              ,intent(in)    :: def_ctr, def_key
    !  real(kind=c_double),dimension(*) ,intent(inout) :: ran1, ran2, ran3, ran4
    !end function threefry4x32f_multi_ctrkey1    
    
    ! 3. any version - all counters/keys supplied
    ! any configuration of counters and keys is possible
    ! either four vectors of random number are returned. 
    !or    
    ! buffer of random number are returned in order as if n scalar calls (for vector length n)      

    integer(c_int) function threefry4x32f_multi_ctrkey_all(ctr1, ctr2, ctr3, ctr4, key1, key2, key3, key4, buff, in_num_vals) bind(c,name='threefry4x32f_multi_ctrkey_all')
      import :: c_int, c_double, num_vals
      integer(kind=c_int),dimension(*) ,intent(in)             :: ctr1, ctr2, ctr3, ctr4
      integer(kind=c_int),dimension(*) ,intent(in)             :: key1, key2, key3, key4
      real(kind=c_double),dimension(num_vals) ,intent(inout)   :: buff
      integer(kind=c_int),value                                :: in_num_vals        
    end function threefry4x32f_multi_ctrkey_all    


  end interface


 contains

  subroutine test_vec_sizes()
    implicit none    
    integer(kind=ik4) :: ierror ! a return value 
       
    ierror = v123_get_vec_sizes()
    
  end subroutine test_vec_sizes

!=========================================================================================================
!single call versions
!=========================================================================================================


  subroutine rng_threefry(number,seed,key)
    !use random123

    implicit none

    real(kind=rk8),dimension(buffersize_single),intent(out)    ::  number

    !integer(kind=ik4),intent(inout)            ::  seed
    integer(kind=ik4),intent(in)               ::  seed

    integer(kind=ik4),intent(in)               ::  key
    !------- for random123

    integer(kind=ik4),dimension(4) :: key4
    integer(kind=ik4),dimension(4) :: sed4
    !integer(kind=ik8),dimension(4) :: key4
    !integer(kind=ik8),dimension(4) :: sed4


    real(kind=rk8)   ,dimension(buffersize_single) :: ran4

    integer(kind=ik4) :: ierror ! a return value
    !------- for random123

    sed4(1) = seed
    sed4(2) = rng_seed
    sed4(3) = rng_seed
    sed4(4) = rng_seed

    key4(1)  = key
    key4(2)  = rng_key
    key4(3)  = rng_key
    key4(4)  = rng_key

    !uses original random123 source - should be compiled no-vec - to prevent disastrous vectorisation
    ierror = threefry4x32f(sed4,key4,ran4)

    !my c version - file should be compiled -no-vec - should be same performance as original
    !ierror = threefry4x32f_novec(sed4,key4,ran4)

    !my c version - identical to novec but file is compiled with vec for comparison (usually much slower)  
    !ierror = threefry4x32f_scalar(sed4,key4,ran4)

    number = ran4

  end subroutine



!=========================================================================================================
!multi call versions
!=========================================================================================================


  subroutine rng_threefry_multi_ss_fix(number,seed,key)
    !use random123

    implicit none

    real(kind=rk8),dimension(buffersize_multi),intent(inout)    ::  number

    !integer(kind=ik4),intent(inout)            ::  seed
    integer(kind=ik4),intent(in)               ::  seed

    integer(kind=ik4),intent(in)               ::  key
    !------- for random123

    integer(kind=ik4),dimension(4) :: key4
    integer(kind=ik4),dimension(4) :: sed4
    !integer(kind=ik8),dimension(4) :: key4
    !integer(kind=ik8),dimension(4) :: sed4


    real(kind=rk8)   ,dimension(buffersize_multi) :: ran4

    integer(kind=ik4) :: ierror ! a return value
    !------- for random123

    sed4(1) = seed
    sed4(2) = rng_seed
    sed4(3) = rng_seed
    sed4(4) = rng_seed

    key4(1)  = key
    key4(2)  = rng_key
    key4(3)  = rng_key
    key4(4)  = rng_key

    !increment is four
    ierror = threefry4x32f_multi_ss_fix(sed4,key4,4,ran4)

    number = ran4

  end subroutine 
  
  !-----------------------------------------------------------------------------

  subroutine rng_threefry_multi_ss(number,seed,key)
    !use random123

    implicit none

    real(kind=rk8),dimension(buffersize_multi),intent(inout)    ::  number

    !integer(kind=ik4),intent(inout)            ::  seed
    integer(kind=ik4),intent(in)               ::  seed

    integer(kind=ik4),intent(in)               ::  key
    !------- for random123

    integer(kind=ik4),dimension(4) :: key4
    integer(kind=ik4),dimension(4) :: sed4
    !integer(kind=ik8),dimension(4) :: key4
    !integer(kind=ik8),dimension(4) :: sed4


    real(kind=rk8)   ,dimension(buffersize_multi) :: ran4

    integer(kind=ik4) :: ierror ! a return value
    !------- for random123

    sed4(1) = seed
    sed4(2) = rng_seed
    sed4(3) = rng_seed
    sed4(4) = rng_seed

    key4(1)  = key
    key4(2)  = rng_key
    key4(3)  = rng_key
    key4(4)  = rng_key

    !increment is four
    ierror = threefry4x32f_multi_ss(sed4,key4,4,ran4,num_vals)

    number = ran4

  end subroutine 

  !-----------------------------------------------------------------------------

  subroutine rng_threefry_multi_all(number,seed,key)
    !use random123

    implicit none

    real(kind=rk8),dimension(buffersize_multi),intent(inout)    ::  number

    !integer(kind=ik4),intent(inout)            ::  seed
    integer(kind=ik4),intent(in)               ::  seed

    integer(kind=ik4),intent(in)               ::  key
    !------- for random123

    integer(kind=ik4),dimension(num_vals)      :: key1, key2, key3, key4 
    integer(kind=ik4),dimension(num_vals)      :: sed1, sed2, sed3, sed4 

    real(kind=rk8)   ,dimension(buffersize_multi) :: ran4

    integer(kind=ik4) :: ierror ! a return value
    integer(kind=ik4) :: i
    !------- for random123

    do i=1, num_vals
      sed1(i)=seed+(i-1)*4
    enddo

    !array notation
    sed2 = rng_seed
    sed3 = rng_seed
    sed4 = rng_seed

    key1  = key
    key2  = rng_key
    key3  = rng_key
    key4  = rng_key

    ierror = threefry4x32f_multi_ctrkey_all(sed1, sed2, sed3, sed4, key1, key2, key3, key4, ran4 ,num_vals)

    number = ran4

  end subroutine 

  !-----------------------------------------------------------------------------


end module rng_wrapper_mod
