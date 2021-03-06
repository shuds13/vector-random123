# vector-random123

## Introduction

vector-random123 (abbr. below as vector123) is a vectorizable port of the counter based random number generators in Random123. The reformatted functions are written in plain C and optimized for performance, clarity and portability.

Random123 resources can be found [here](http://www.thesalmons.org/john/random123)

Note that while the random numbers generated for a given set of counters/keys in the vectorised functions should be identical to the originals, there is no guarantee of correctness. This project is still in development. 

Initial focus is on the Threefry generator. This generator is composed mostly of simple operations which are single cycle latency on many platforms. This makes it more performance portable than Philox which requires implementation of mulhilo. It also tends to perform well for both on-the-fly and batch generation approaches. Philox may be included in the future, however.

Current performance results can be viewed under [Benchmarking](#benchmarking) below.


## Current status 


* Only threefry4x32 fully implemented and tested

* threefry4x32, threefry4x64, threefry2x64 have also been implemented.

* Though answers tested have been identical to the crush resistant Random123, a thorough QA has not been carried out and no RNG test batteries have yet been run.

* The functions currently output doubles between 0 <= x < 1 (closed/open) by default. Different options can be found in the Random123 file u01.h or in more recent versions u01fixedpt.h.
  
* Documentation is limied but some information is provided in the source directory and source files themselves. 

* Currently the vector123 functions are set up to be standalone. However, see the comments at top of source files to see how to pick up the Random123/ directory. This provides access to some of the features files for different systems.

 
## Instructions for use

Instructions for incorporating the source files are in the readme in the source directory. The source files themselves also contain explanations of the functions.

The loop-only programs can be run to test scalar and vector performance on a given system. These programs simply generate a large array of random numbers (either ints, floats or doubles) along with timing information. The existing build scripts can be easily modified to create either scalar or vector builds for intel systems. Of course these will need to be modified for other platforms. A version written using AVX512 vector intrinsics is included which, in current benchmarking, is about 30% faster than the source code version on the KNL (Intel Knight's Landing Xeon Phi). Threaded versions of these routines are also available.

The fortran-test-driver/ directory contains an example of calling the vector versions of threefry4x32 within a Fortran progam. It also tests performance of different implementations and compares with the original scalar routines, including the calling overheads that maybe incurred in a typical application.

### OpenMP ###

It is preferable to have OpenMP enabled, as the OpenMP SIMD directive is used on the main loop. This makes it more likely for the compiler to successfully vectorise the loop, as well as indicating alignment of work arrays. This requires that OpenMP 4.0 or higher is supported. The loop should vectorise without these pragmas enabled, but it is always worth checking. It is expected that threading will normally be implemented above this level, but threaded versions of the code can be found under the loop-only directory. The functions could also be structured as OpenMP SIMD functions using scalar arguments.
<br />

## Benchmarking

Results from benchmarking on Ivy Bridge, Haswell and Xeon Phi (KNL).

See loop-only directory for more details - including source code used and full results.

The results below are all based on the threefry4x32 generator with the default 20 rounds.


### Single threaded performance ###

These tests show that the threefry generator on its own gains very impressive vectorisation
performance - near perfect when taking into account the relative number of ALUs and vector units
on each system. In particular KNL, when using the fast rol intrinsic has 16.3x over non-vectorised
version - 100% vector efficiency. This should be supported from source code in future compiler versions.


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

Vector v Scalar with conversion to float:
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

A popular performance metric for RNGs is cycles per byte (cpB). KNL has a signficantly better cpB than Haswell,  taking advantage of better vector performance, while this metric is not affected by KNLs lower clockrate. The graph below shows a comparison of cpB for single thread runs on each of the platforms. This includes the KNL intrinsics version which achieves a cpB of 0.48.

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

