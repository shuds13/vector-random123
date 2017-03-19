## Directory Contents

Currently, these are standalone files where the appropriate file can be merged into users source
code, as it is likely users will want to only use one variation of routines. Example makefile
options can be seen in the fortran driver directory. Thus the files in here include replication
of definitions.

An alternative library format which will employ shared header files may be implemented in the future.
Also routines may be renamed.

### File variants

The b and r variants are based on the format in which the random numbers are supplied.

eg 
v123_threefry4x32b_dble.c  Routines that return a single array (buffer) of RNs.
v123_threefry4x32r_dble.c: Routines that return 4 arrays corresponding to the 4 input counters.

By default the b version should be easy to compare values with the original scalar version called
repeatedly writing values into a buffer. 

For example, if the routine threefry4x32f_multi_ss (in v123_threefry4x32b_dble.c) is called with a given
counter and key set, an increment of one, and an array/vector size of 8, then the buffer returned
should contain the same 32 values, as if the original threefry4x32 is called 8 times (incrementing
the counter by one each time), writing  the 4 values of each call into one buffer consecutively.

The r variant may lend itself more naturally to, for example, an array/vector of particles are being
processed in parallel.

The suffix _dble deontes that a double is returned. For floats modify the buff argument to a float and
adjust the conversion according to the Random123 u01.h or u01fixedpt.h header files.

More details of each routine are given in the files.


### List of current variants


v123_threefry4x32b_dble.c 
v123_threefry4x32r_dble.c 

v123_threefry4x64b_dble.c 

v123_threefry2x64b_dble.c




### Proposed renaming

The routines may be renamed in the following format.

v123_threefry4x32b_dble.c:  threefry4x32f_multi_ss_fix         to  v123_threefry4x32b_ss_fix
v123_threefry4x32b_dble.c:  threefry4x32f_multi_ss             to  v123_threefry4x32b_ss
v123_threefry4x32b_dble.c:  threefry4x32f_multi_ctrkey_all     to  v123_threefry4x32b_ctrkey_all

v123_threefry4x32r_dble.c:  threefry4x32f_multi_ss_fix         to  v123_threefry4x32r_ss_fix
v123_threefry4x32r_dble.c:  threefry4x32f_multi_ss             to  v123_threefry4x32r_ss
v123_threefry4x32r_dble.c:  threefry4x32f_multi_ctrkey_all     to  v123_threefry4x32r_ctrkey_all

etc...

----------------------------------------------------------------------------------------------------
