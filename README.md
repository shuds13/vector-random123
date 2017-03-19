## vector-random123

vector-random123 (abbr. below as vector123) is a vectorizable port of the counter based random number generators in Random123. The reformatted functions are written in plain C and optimized for performance, clarity and portability.

Random123 resources can be found [here](http://www.thesalmons.org/john/random123)

Note that while the random numbers generated for a given set of counters/keys in the vectorised functions should be identical to the originals, there is no guarantee of correctness. This project is still in development. Initial focus is on the Threefry generator.


## Current status 


* Only threefry4x32 fully implemented and tested

* threefry4x32, threefry4x64, threefry2x64 have been implemented.


* Though answers tested have been identical to the crush resistant Random123, a thorough QA has not been carried out and no RNG test batteries have yet been run.

* __Priority__ The functions currently output doubles between 0 <= x < 1 (closed/open) by default. Different options will be provided soon.
  
* Documentation - including routines available will be added soon. Some information is provided in the source files (see source directory) files and in fortran harness see interface section in rng_wrapper.F for some info.
 
* Results for different platforms - also to follow
 
* Current version contains vectorisable loop inside the function. A vectorisable version that can be called inside a loop is to follow.

* Currently the vector123 functions are set up to be standalone. However, see the comments at top of source files to see how to pick up the Random123/ directory. This provides access to some of the features files for different systems.

 
## Instructions for use

Note: The fortran-test-driver/ directory contains an example of using the vector versions of threefry4x32. It also tests performance of different implementations. Running this may be the quickest way to go.


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
 
It is not difficult to use these values for initialisation of data structures in the calling code - and for many systems the processors capability can be detected automatically - eg. `#if defined(__AVX512F__)` See commented out lines at top of the source files Alternatively you can use the functions that pass in the vector size (or number of sets). Depending on the system, using a compile time constant for the vector length may be either the same or faster than passing in as an argument. It is worth always checking that the loop has succesfully vectorised.

### OpenMP ###

It is preferable to have OpenMP enabled, as the OpenMP SIMD directive is used on the main loop. This makes it more likely for the compiler to successfully vectorise the loop, as well as indicating alignment of work arrays. This requires that OpenMP 4.0 or higher is supported. The loop should vectorise without these pragmas enabled, but it is always worth checking.


<br />

## Benchmarking ##

Results from benchmarking on Ivy Bridge, Haswell and Xeon Phi(KNL)

See loop-only directory for more details - including source code used and full results.

All results were based on threefry4x32


### Single threaded performance ###

These tests show that the threefry generator on its own gains very impressive vectorisation
performance - near perfect when taking into account the relative number of ALUs and vector units
on each system. In particular KNL, when using the fast rol intrinsic has 16.3x over non-vectorised
version - 100% vector efficiency. This should be supported from source code in future compiler versions.
The results below are all based on the threefry4x32 generator with the default 20 rounds.

Comparison producing 64 million random numbers.

Vector v Scalar for pure threefry4x32 generator:
* Ivy Bridge ~ 2.9x
* Haswell    ~ 5.2x
* KNL        ~ 12.8 / 16.4x(AVX512 intrinsics)

![image](https://cloud.githubusercontent.com/assets/16457059/24082629/ade23b6a-0cc0-11e7-95cd-03a791b16d1a.png)  

The conversion to floats or double, however, takes longer than the generator and reduces speed up significantly.
In particular, converting 32-bit ints to doubles will require a split of the vector. The resuts below do not
include the intrinsics version.

Vector v Scalar with conversion to double:
* Ivy Bridge ~ 1.8x
* Haswell    ~ 2.3x
* KNL        ~ 6.8x

![image](https://cloud.githubusercontent.com/assets/16457059/24082662/41add746-0cc1-11e7-9367-4481cad143fd.png)  

Vector v scalar with conversion to float:
* Ivy Bridge ~ 2.1x
* Haswell    ~ 3.1x
* KNL        ~ 9.4x

![image](https://cloud.githubusercontent.com/assets/16457059/24082647/12ec9280-0cc1-11e7-96cb-bb1256246e74.png)  

#### Haswell v KNL ####

The results below are based on the source code versions.

While KNL has the best speed up from vectorisation - Haswell is 1.6x faster than KNL for single thread performance
on pure threefry (0.078sec to 0.123). Note: The AVX512 intrinsics version on KNL achieves 0.096sec

KNL seems to do better with the conversions, however.  
With conversion to floats KNL has parity with Haswell at single thread.  0.179 v 0.183(knl)

Likewise with conversion to doubles - 0.274 v 0.264(knl)

Given the lower clockrate KNLs cycles per byte (cpB) is signficantly better than Haswell, taking advantage of better vector performance. The graph below shows a comparison of cpB for single thread runs on each of the platforms. This includes the KNL intrinsics version which achieves a cpB of 0.48.

![image](https://cloud.githubusercontent.com/assets/16457059/24082740/6860d78e-0cc2-11e7-9d2c-c48ce1441949.png)  

### Fully populated processors - using OpenMP threads ###

Comparison producing 640 million random integers (no float/double conversions).

E5-2680v3 Haswell 12-core 2.5GHz  
v  
KNL-7210 (64 cores - 4HW threads per core) 1.3GHz

Best results on one 12-core Haswell:
* 24 threads (2 threads/core): = 0.077792 seconds (3.1 x non-vec code)

Best results on one KNL-7210:
* 128 threads (2HW threads/core): = 0.015215 seconds (14.0 x non-vec code)

KNL = 5.1x over Haswell

![image](https://cloud.githubusercontent.com/assets/16457059/24082837/c326f9f4-0cc3-11e7-984f-f6b94194716a.png)  

These results suggest that the vectorised threefry generator, on KNL in particular, is producing close to peak integer operations per second. Quantification of this to follow.

