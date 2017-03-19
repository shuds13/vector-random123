program test_r123
  use precision_mod
  use rng_wrapper_mod
  use mpi
  implicit none


  integer(kind=ik4)                       :: buffersize !Number per r123 call x no of r123 call per rng call
  
  integer(kind=ik4)                       :: n = 10000000
  integer(kind=ik4)                       :: num_rounds = 2
  
  integer(kind=ik4)                       :: i,position,round
  
  real(kind=rk8)                          :: t1, t2, time1, time2
  real(kind=rk8)                          :: time_tfry_single,time_tfry_multi_ss_fix,time_tfry_multi_ss,time_tfry_multi_all
  !real(kind=rk8)                          :: time_phil_single,time_phil_multi

  real(kind=rk8),allocatable              :: r(:)
  !real(kind=rk8),dimension(buffersize)    :: r
  
  real(kind=rk8)                          :: rval
  real(kind=rk8)                          :: rsum1,rsum2,rsum3,rsum4,rsum5
  
  integer(kind=ik4)                       :: fixed_key = 10
  integer(kind=ik4)                       :: mpierr,ierr


!   INTEGER(KIND=ik4)                       ::  test_vec_sizes
!   EXTERNAL                                    test_vec_sizes


  time1 = 0.0_rk8 
  
  
  !For mpi_wtime
  call mpi_init(mpierr)

  write(*,*) ''   
  write(*,111) THREEFRY_NAME
  !write(*,111) THREEFRY_NAME,PHILOX_NAME
  
  !Check
  write(*,*) ''   
  !ierr = test_vec_sizes()
  call test_vec_sizes()
  write(*,3) NUM_VALS
  


  !Note to avoid confusion
  !Using i as counter when generate random123 - but as each call generates X values - i'm skipping values
  !Could simply set separate generator value and use all, but this is way i've done for now.

!=========================================================================================================
!SINGLE CALL VERSIONS
!=========================================================================================================
  
  write(*,*) '' 
  write(*,20)
  write(*,201) 
  write(*,*) '' 
  
  buffersize=BUFFERSIZE_SINGLE
  allocate(r(buffersize))
  
  rsum1 = 0.0_rk8 
  do round=1,num_rounds
  
    write(*,123) round
   
    t1     = mpi_wtime()     
    do i=1,n
      !Note - same thing of - skipping values - to keep N as seed (as each seed generates 4 values)
      position=mod(i-1,buffersize)
      !print*,'i,position', i,position !test

      if (position==0) then
        !In rirst instance doing one r123call with orig routine so 4 values 
        !call rng_threefry(r,int(i,KIND=8),INT(fixed_key,KIND=8))
        call rng_threefry(r,i,fixed_key)
        rval=r(1)
      else 
        rval=r(position+1)      
      endif

!       !Debug  - print some values
!       if (i.lt.10) then
!         print*,'rval',rval
!       endif

      
      rsum1 = rsum1 + rval !accum a sum

    enddo  
    t2     = mpi_wtime()
    time1 = t2 - t1
    write(*,11) round,time1
    
    time_tfry_single=time1 !Added later for analysis 

    write(*,*) ''

  enddo !rounds  
  
  write(*,124) rsum1
  write(*,*) ''
  deallocate(r)

!=========================================================================================================
!MULTI CALL VERSIONS - SS1 VARIANTS
!=========================================================================================================
  
  write(*,*) '' 
  write(*,20)
  write(*,202) NUM_VALS
  write(*,*) '' 
  
  buffersize=BUFFERSIZE_MULTI
  allocate(r(buffersize))
  
  rsum2 = 0.0_rk8 
  do round=1,num_rounds
  
    write(*,123) round
   
    t1     = mpi_wtime()     
    do i=1,n
      !Note - same thing of - skipping values - to keep N as seed (as each seed generates 4 values)
      position=mod(i-1,buffersize)
      !print*,'i,position', i,position !test

      if (position==0) then
        !In rirst instance doing one r123call with orig routine so 4 values 
        !call rng_threefry(r,int(i,KIND=8),INT(fixed_key,KIND=8))
        call rng_threefry_multi_ss_fix(r,i,fixed_key)
        rval=r(1)
      else 
        rval=r(position+1)      
      endif

!       !Debug  - print some values
!       if (i.lt.10) then
!         print*,'rval',rval
!       endif
      
      rsum2 = rsum2 + rval !accum a sum

    enddo  
    t2     = mpi_wtime()
    time1 = t2 - t1 
    write(*,12) round,time1
    time_tfry_multi_ss_fix=time1 !Added later for analysis 


    write(*,*) ''

  enddo !rounds  
  
  write(*,124) rsum2
  write(*,*) ''
  deallocate(r)
  
  !----------------------------------------------------------------------------------------
  
  write(*,*) '' 
  write(*,20)
  write(*,203) NUM_VALS
  write(*,*) '' 
  
  buffersize=BUFFERSIZE_MULTI
  allocate(r(buffersize))
  
  rsum3 = 0.0_rk8 
  do round=1,num_rounds
  
    write(*,123) round
   
    t1     = mpi_wtime()     
    do i=1,n
      !Note - same thing of - skipping values - to keep N as seed (as each seed generates 4 values)
      position=mod(i-1,buffersize)
      !print*,'i,position', i,position !test

      if (position==0) then
        !In rirst instance doing one r123call with orig routine so 4 values 
        !call rng_threefry(r,int(i,KIND=8),INT(fixed_key,KIND=8))
        call rng_threefry_multi_ss(r,i,fixed_key)
        rval=r(1)
      else 
        rval=r(position+1)      
      endif

!       !Debug  - print some values
!       if (i.lt.10) then
!         print*,'rval',rval
!       endif
      
      rsum3 = rsum3 + rval !accum a sum

    enddo  
    t2     = mpi_wtime()
    time1 = t2 - t1 
    write(*,13) round,time1
    time_tfry_multi_ss=time1 !Added later for analysis 


    write(*,*) ''

  enddo !rounds  
  
  write(*,124) rsum3
  write(*,*) ''
  deallocate(r)

!=========================================================================================================
!MULTI CALL VERSIONS - ALL VARIANT
!=========================================================================================================
  
  write(*,*) '' 
  write(*,20)
  write(*,205) NUM_VALS
  write(*,*) '' 
  
  buffersize=BUFFERSIZE_MULTI
  allocate(r(buffersize))
  
  rsum5 = 0.0_rk8 
  do round=1,num_rounds
  
    write(*,123) round
   
    t1     = mpi_wtime()     
    do i=1,n
      !Note - same thing of - skipping values - to keep N as seed (as each seed generates 4 values)
      position=mod(i-1,buffersize)
      !print*,'i,position', i,position !test

      if (position==0) then
        !In rirst instance doing one r123call with orig routine so 4 values 
        !call rng_threefry(r,int(i,KIND=8),INT(fixed_key,KIND=8))
        call rng_threefry_multi_all(r,i,fixed_key)
        rval=r(1)
      else 
        rval=r(position+1)      
      endif

!       !Debug  - print some values
!       if (i.lt.10) then
!         print*,'rval',rval
!       endif
      
      rsum5 = rsum5 + rval !accum a sum

    enddo  
    t2     = mpi_wtime()
    time1 = t2 - t1 
    write(*,15) round,time1
    time_tfry_multi_all=time1 !Added later for analysis 

    write(*,*) ''

  enddo !rounds  
  
  write(*,124) rsum5
  write(*,*) ''
  deallocate(r)
!=========================================================================================================
  
  ! Analysis
  write(*,223)  
  if (rsum1 .eq. rsum2) then
     write(*,*) 'Rsum single/multi fixlen (ss1) matches'
  else  
     write(*,*) 'Rsum single/multi fixlen (ss1) =', rsum1/rsum2
  endif

  if (rsum1 .eq. rsum3) then
     write(*,*) 'Rsum single/multi arglen (ss1) matches'
  else  
     write(*,*) 'Rsum single/multi arglen (ss1) =', rsum1/rsum3
  endif

  if (rsum1 .eq. rsum5) then
     write(*,*) 'Rsum single/multi arglen (all) matches'
  else  
     write(*,*) 'Rsum single/multi arglen (all) =', rsum1/rsum5
  endif

  write(*,*) ''  
  write(*,*) 'Round 2 times: single v multi(ss1)'
  write(*,224) time_tfry_single/time_tfry_multi_ss_fix
  write(*,225) time_tfry_single/time_tfry_multi_ss
  write(*,228) time_tfry_single/time_tfry_multi_all
  write(*,229) time_tfry_multi_ss/time_tfry_multi_all

  !write(*,205) time_phil_single/time_phil_multi
  write(*,20)  
  
  write(*,*) ''
  
  
  call mpi_finalize(mpierr)
  
!1    FORMAT('Threefry  Time:',F10.6)

!For grepping have diff prints
11    FORMAT('Threefry  Time single           round ',i3,':',F10.6)
12    FORMAT('Threefry  Time multi ss1 fixlen round ',i3,':',F10.6)
13    FORMAT('Threefry  Time multi ss1 arglen round ',i3,':',F10.6)
15    FORMAT('Threefry  Time multi all arglen round ',i3,':',F10.6)

20   FORMAT('----------------------------------------------')

!2    FORMAT('Philox    Time:',F10.6)
3    FORMAT('Fortran vector length: NUM_VALS', i4 , '  - should match 32 or 64 value above')

123  FORMAT('Round:  ',i3)
124  FORMAT('Rsum:   ',F17.6)

201  FORMAT('Single Call version')
202  FORMAT('Multiple Call version Single stream fixlen (SSF):',i3,' calls at a time')
203  FORMAT('Multiple Call version Single stream arglen (SS1):',i3,' calls at a time')
205  FORMAT('Multiple Call version Any config    arglen (ALL):',i3,' calls at a time')

223  FORMAT('Analysis -------------------------------------')
224  FORMAT(' Threefry multi(ss1) fixlen speedup: ',F5.2)
225  FORMAT(' Threefry multi(ss1) arglen speedup: ',F5.2)
228  FORMAT(' Threefry multi(all) arglen speedup: ',F5.2)
229  FORMAT(' Threefry multi(all) slowdown (<1) relative to multi(ss1): ',F5.2)

!205  FORMAT(' Philox   multi speedup:  ',F5.2)
!206  FORMAT('----------------------------------------')

111 FORMAT('Using ',a12)
!111 FORMAT('Using ',a12, ' and ',a10)

! 1    FORMAT('Threefry ',i3, '   Time:',F10.6)
! 2    FORMAT('Philox   ',i3, '   Time:',F10.6)

end program test_r123








