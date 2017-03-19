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


### Setting Parameters ###

Modify the following paramaters (macros) according to the vector length of system.

VECTOR_LENGTH_BYTES  
NUM_VALS_32  
NUM_VALS_64  

VECTOR_LENGTH_BYTES is used for alignment of data and should be at least the vector length in bytes for your system.

NUM_VALS_32 and NUM_VALS_64 are only used for the function variants that have a "_fix" suffix. In others the vector size is passed as an argument.

NUM_VALS_32 is the number of sets produced (loop iteration count) for *x32 routines (eg threefry4x32f_multi_ss_fix). This should be the integer vector size of system in terms of data type or a multiple thereof.
 
eg. For threefry4x32f_multi_ss_fix - if NUM_VALS_32 is 8 - that is equivalent to doing 8 calls of threefry4x32 in random123. As each threefry4x32 produces a set of 4 values - in total 32 values will be returned.

Examples:

eg. AVX-512  
 #define VECTOR_LENGTH_BYTES 64  
 #define NUM_VALS_32 16  
 #define NUM_VALS_64 8  

eg. AVX2  
 #define VECTOR_LENGTH_BYTES 32  
 #define NUM_VALS_32 8  
 #define NUM_VALS_64 4  

eg. AVX/SSE  
 #define VECTOR_LENGTH_BYTES 32  
 #define NUM_VALS_32 4  
 #define NUM_VALS_64 2  

Note integer vector operations in AVX are limited to AVX-128 - so the same as SSE - however the VECTOR_LENGTH_BYTES should be aligned to full vector size.

As many values as you want can be produced - with multiples of vector size preferable.
Eg. You could have:

eg. AVX2  
 #define VECTOR_LENGTH_BYTES 32  
 #define NUM_VALS_32 8000  
 #define NUM_VALS_64 4000  
 
It is not difficult to use these values for initialisation of data structures in the calling code - and for many systems the processors capability can be detected automatically - eg. `#if defined(__AVX512F__)` To use automatic detection set the macro AUTO_SELECT to 1 in the source files. The values can be checked with a call to v123_get_vec_sizes().

Alternatively you can use the functions that pass in the vector size (or number of sets). Depending on the system, using a compile time constant for the vector length may be either the same or faster than passing in as an argument. It is worth always checking that the loop has succesfully vectorised.

### OpenMP ###

It is preferable to have OpenMP enabled, as the OpenMP SIMD directive is used on the main loop. This makes it more likely for the compiler to successfully vectorise the loop, as well as indicating alignment of work arrays. This requires that OpenMP 4.0 or higher is supported. The loop should vectorise without these pragmas enabled, but it is always worth checking. It is expected that threading will normally be implemented above this level, but threaded versions of the routines can be found under the loop-only directory. The functions could also be structured as openmp SIMD subroutines using scalar arguments. 

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
