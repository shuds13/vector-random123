## vector-random123

vector-random123 (abbr. below as vector123) is a vectorizable port of the counter based random number generators in Random123. The reformatted functions are written in plain C and optimized for performance and clarity.

Random123 resources can be found [here](http://www.thesalmons.org/john/random123)

Note that while the random numbers generated for a given set of counters/keys in the vectorised functions should be identical to the originals, there is no guarantee of correctness. This project is still in development. Initial focus is on the Threefry generator.


## Current status 


* Only threefry4x32 fully implemented

* Though answers tested have been identical to the crush resistant Random123, a thorough QA has not been carried out and no RNG test batteries have yet been run.

* __Priority__ The functions currently output doubles between 0 <= x < 1 (closed/open) by default. Different options will be provided soon.
  
* Documentation - including routines available will be added soon. Currently see vector123.c files and in fortran harness see interface section in rng_wrapper.F for some info.
 
* Results for different platforms - also to follow
 
* Current version contains vectorisable loop inside the function. A vectorisable version that can be called inside a loop is to follow.

* Currently the vector123 functions are set up to be standalone. However, see the comments at top of vector123.c files to see how to pick up the Random123/ directory. This provides access to some of the features files for different systems.

 
## Instructions for use

Note: The fortran-test-driver/ directory contains an example of using the vector versions of threefry4x32. It also tests performance of different implementations. Running this may be the quickest way to go.


### Setting Parameters ###

In vector123.c modify the following paramaters (macros) according to the vector length of system.

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
 
It is not difficult to use these values for initialisation of data structures in the calling code - and for many systems the processors capability can be detected automatically - eg. `#if defined(__AVX512F__)` See commented out lines at top of vector123.c Alternatively you can use the functions that pass in the vector size (or number of sets). Depending on the system, using a compile time constant for the vector length may be either the same or faster than passing in as an argument. It is worth always checking that the loop has succesfully vectorised.

### OpenMP ###

It is preferable to have OpenMP enabled, as the OpenMP SIMD directive is used on the main loop. This makes it more likely for the compiler to successfully vectorise the loop, as well as indicating alignment of work arrays. This requires that OpenMP 4.0 or higher is supported. The loop should vectorise without these pragmas enabled, but it is always worth checking.
