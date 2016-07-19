Vector-random123 (a.k.a Vector123) is a vectorizable port of the counter based random 
number generators in Random123. The reformatted funtions are written in plain C and
optimised for performance and clarity.

Random123 can be found *here*


--------------------------------------------------------------------------------

Note: There is no guarrantee of correctness - user should test their code 

Answers for vectorised functions should be identical to originals.


--------------------------------------------------------------------------------
**Note - current status 
--------------------------------------------------------------------------------

-Only threefry4x32 fully implemented

-Though answers tested have been identical to the crush resistant Random123, a thorough QA
 has not been carried out and no RNG test batteries have yet been run.

-currently set to creates values between 0 < x <= 1 by default
 **this will be made generic in future.
 
-documentation - including routines available will be added soon.
 Currently see vector123.c files and in fortran harness see interface section 
 in rng_wrapper.F for some info.
 
-Results for different platforms - also to follow
 
-Current version contains vectorisable loop inside the function. A vectorisable version
that can be called inside a loop is to follow.

-Currently the vector123 functions are set up to be standalone. However, see the comments at top
of vector123.c files to see how to pick up the Random123/ directory. This provides access
to some of the features files for different systems.

 
--------------------------------------------------------------------------------
Instructions for use
--------------------------------------------------------------------------------

Note: The fortran-test-driver/ directory contains an example of using the vector versions of
threefry4x32. It also tests performance of different implementations. Running this may be the
quickest way to go.

-------------------------

In  vector123.c set the following for vector length of system.

VECTOR_LENGTH_BYTES
NUM_VALS_32
NUM_VALS_64

Default values will be set

VECTOR_LENGTH_BYTES is used for alignment of data and should be at least the vector length
in bytes for your system.

NUM_VALS_32 and NUM_VALS_64 are only used for "fix" routines - where vector size is not passed
as an argument.

NUM_VALS_32 is the number of sets produced (loop iteration count) for *x32 routines
 (eg threefry4x32f_multi_ss_fix). Should be integer vector size of system in terms of data type
 or a multiple thereof.
 
eg. For threefry4x32f_multi_ss_fix - if NUM_VALS_32 is 8 - that is equivalent to doing 8 calls
    of threefry4x32 in random123.
    
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

Note integer vector operations in AVX are limited to AVX-128 - so same as SSE - however
the VECTOR_LENGTH_BYTES should be aligned to full vector size.

As many values as you want can be produced - with multiples of vector size preferable.
eg could have.

eg. AVX2
 #define VECTOR_LENGTH_BYTES 32
 #define NUM_VALS_32 8000
 #define NUM_VALS_64 4000


It is not to hard to set up to get these values from vector123 - and for many systems
from preproc defines - eg. __avx512f__. See commented out lines at top of vector123.c





