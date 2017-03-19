//--------------------------------------------------------------------------------------------------
// S. Hudson: Vector_Random123: threefry4x32
// Standard C vectorizable implementation of Random123 threefry4x32 functions.
// Based on Random123 functions - see copyright notice.
//--------------------------------------------------------------------------------------------------

/*
Copyright 2010-2011, D. E. Shaw Research.
All rights reserved.

Redistribution and use in source and binary forms, with or without
modification, are permitted provided that the following conditions are
met:

* Redistributions of source code must retain the above copyright
  notice, this list of conditions, and the following disclaimer.

* Redistributions in binary form must reproduce the above copyright
  notice, this list of conditions, and the following disclaimer in the
  documentation and/or other materials provided with the distribution.

* Neither the name of D. E. Shaw Research nor the names of its
  contributors may be used to endorse or promote products derived from
  this software without specific prior written permission.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
"AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
(INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
*/

//---------------------------------------------------------------------------------------------------
// S. Hudson: Standalone version - This file should work standalone by incorporating into code as is.
//                                 See "Notes" for picking up random123 configuration settings.
//---------------------------------------------------------------------------------------------------

/*

Overview:

  Below is a list of available routines which implement by default 20 rounds of threefry4x32. Instead of a single
  set of 4 values, there routines produce multiple sets as described below. The variants allow for simplified
  interfaces for single streams (ss).

  threefry4x32f_multi_ss_fix:     Single stream version using vector/array length specified in this file (NUM_VALS_32)
  threefry4x32f_multi_ss:         Single stream version using provided vector/array length
  threefry4x32f_multi_ctrkey_all: Multi-stream version - user supplies explicit counters/keys for all RNs (Random Numbers).

  Utility Routines
  v123_get_vec_sizes:             Returns parameter/macro values for vector and array sizes
                                  VECTOR_LENGTH_BYTES / NUM_VALS_32 / NUM_VALS_64
  
  --------------------------------------------------------------------------------------------------------------------------------
  
  Original scalar code: 
  
  For comparison, the original threefry4x32 code passes one set of counters (as array) and one set of
  keys (as array) and produce one set of RNs.
  (ctr1, ctr2, ctr3, ctr4) and (key1, key2, key3, key4) - produces (ran1, ran2, ran3, ran4)

  --------------------------------------------------------------------------------------------------------------------------------
  
  Vector code:

  Single stream variants:
  threefry4x32f_multi_ss_fix and threefry4x32f_multi_ss:
         
  Input: Pass one set of counters (as array) and one set of keys (as array) as with scalar code.
  CTR4 = (ctr1, ctr2, ctr3, ctr4) and KEY4 = (key1, key2, key3, key4)
  
  Further sets are generated as below - given a supplied increment of 1

  Array length = 8 ---------------> 
                  
   0       1      2       3       4       5       6       7     
  ctr1   ctr1+1  ctr1+2  ctr1+3  ctr1+4  ctr1+5  ctr1+6  ctr1+7  
  ctr2   ctr2    ctr2    ctr2    ctr2    ctr2    ctr2    ctr2  
  ctr3   ctr3    ctr3    ctr3    ctr3    ctr3    ctr3    ctr3  
  ctr4   ctr4    ctr4    ctr4    ctr4    ctr4    ctr4    ctr4  
  
  Note: Only one ctr value need be changed per set to produce four different RNs as output.

  In this version, the code returns four vectors of random numbers corresponding to the array length above.
  E.g. Array length = 8 (32 RNs will be produced)
  Returned arrays:
  RAN1 = (ran1(0), ran1(1), ran1(2), ran1(3), ran1(4), ran1(5), ran1(6), ran1(7))
  RAN2 = (ran2(0), ran2(1), ran2(2), ran2(3), ran2(4), ran2(5), ran2(6), ran2(7))
  RAN3 = (ran3(0), ran3(1), ran3(2), ran3(3), ran3(4), ran3(5), ran3(6), ran3(7))
  RAN4 = (ran4(0), ran4(1), ran4(2), ran4(3), ran4(4), ran4(5), ran4(6), ran4(7))
  
  The array size can be anything, but the machine vector length or a multiple thereof is preferable for performance.
  
  For threefry4x32f_multi_ss_fix, the array size is set in this file (NUM_VALS_32) which is either a specified
  default or, if a vector instruction set is recognised at compilation by the listed predefined macros, it is set to
  the machine vector length. This must, of course, be consistent with calling code. See NUM_VALS_32 below. Note that
  this constant is not used by the other routines which take the array length as an argument.


  threefry4x32f_multi_ctrkey_all:
  The routine threefry4x32f_multi_ctrkey_all enables the user to supply any permutation of input counters and keys.
  The array arguments supplied for counter and keys are of array/vector length with a seperate array argument for each
  counter/key.
  
  CTR1 = (ctr1_v0, ctr1_v1, ctr1_v2, ctr1_v3, ctr1_v4, ctr1_v5, ctr1_v6, ctr1_v7 ...)
  CTR2 = (ctr2_v0, ctr2_v1, ctr2_v2, ctr2_v3, ctr2_v4, ctr2_v5, ctr2_v6, ctr2_v7 ...)
  etc...
  
  --------------------------------------------------------------------------------------------------------------------------------


Notes:

  This file should work standalone by incorporating into code as is.
  Alternatively, Random123 include files (threefry.h) can be used
  and will pick up the compiler feature files etc.
  which sets some macros  (eg.directives for inlining)
  eg #define R123_FORCE_INLINE(decl) decl __attribute__((__always_inline__))
  Also if use threefry.h  - comment out section indicated below to avoid duplicate definitions
  -> so far no difference has been found on platforms tested
  Not tried this on CUDA platforms..

  By default these routines produce RNs which are doubles uniformly distributed between 0 and 1 inclusive/exlusive.
  Alternative conversions can be found in Random123 header files u01.h (up to 1.08) or u01fixedpt.h

  *Alternative routines exist which produces a single buffer of RNs as output. 
  
  The current code is laid out for clarity and ease of debugging/profiling - rather than minimal code size or coverage
  of all possible variations. The code can easily be modified for variants as desired.
  
  OpenMP pragmas: It is not essential to enable OpenMP and there is no threading inside this routine. However,
  the OpenMP pragmas help to ensure loops are vectorized.

*/


// ------------------------------------------------------------------------------------------------/
// ------------- Automatic determination of vector length for threefry4x32f_multi_ss_fix ----------/
// ------------------------------------------------------------------------------------------------/

// The macro VECTOR_LENGTH_BYTES is used for alignment of arrays in all routines (on supported platforms).
// This should be at least the machine vector length in bytes.

// NUM_VALS_32 is used only in routines with the _fix suffix (which refers to a fixed array size). This is actually
// the number of sets of values (see above).

// For supported platforms the conditional preprocessor code below will set these macros to architecture vector lengths.
// These can be overidden here by the user. 

// Generally I dont think should matter if too big - as long as fits with calling routines and multiple of correct size.
// VECTOR_LENGTH_BYTES // For alignment if using - should not matter if too big (avx256=32 avx512=64) 
// NUM_VALS_32         // Number of R123 sets of 32-bit values. Should be  vector length (or multiple of) - LOOP SIZE
// NUM_VALS_64         // Number of R123 sets of 64-bit values. Should be  vector length (or multiple of) - LOOP SIZE


// Enable/disable auto-selection
#define AUTO_SELECT 0
//#define AUTO_SELECT 1

// ----------------------------------------------------------
// Set values for supported platforms (*may move to sep. header file - shared by diff generators)
// Initially this is focussed on intel platforms - for others it is recommended to set manually
#if AUTO_SELECT

// AVX macros etc.. are defined in many non-intel compilers also (including gcc) - so do not specify __ICC here.
//#if defined(__ICC)

#if   defined(__AVX512F__)
#define VECTOR_LENGTH_BYTES 64
#define NUM_VALS_32 16
#define NUM_VALS_64 8
#elif   defined(__AVX2__)
#define VECTOR_LENGTH_BYTES 32
#define NUM_VALS_32 8
#define NUM_VALS_64 4
#elif   defined(__AVX__)
#define VECTOR_LENGTH_BYTES 32
#define NUM_VALS_32 4 //AVX1 - integers limited to AVX-128
#define NUM_VALS_64 2
#elif   defined(__SSE__)
#define VECTOR_LENGTH_BYTES 32
#define NUM_VALS_32 4
#define NUM_VALS_64 2
#else
//ICC DEFAULTS
#define VECTOR_LENGTH_BYTES 64
#define NUM_VALS_32 8
#define NUM_VALS_64 4
#endif

//#endif // (__ICC)

#else // If not using AUTO_SELECT
//DEFAULTS
#define VECTOR_LENGTH_BYTES 64
#define NUM_VALS_32 16
#define NUM_VALS_64 8

#endif //AUTO_SELECT
// ----------------------------------------------------------


// ------------------------------------------------------------------------------------------------/
// ------------------------------ Include Header Files --------------------------------------------/
// ------------------------------------------------------------------------------------------------/

#include <stdio.h>

// Allows use of short forms like uint32_t
#include <stdint.h>


//Random123 - To pick up macros from R123 - uncomment here and comment out section below
//#include <threefry.h>

//#include <u01.h>
//OR
//#include <u01fixedpt.h>


// ------------------------------------------------------------------------------------------------/
// ------------- For standalome version - comment out if using threefry.h etc ---------------------/
// ------------------------------------------------------------------------------------------------/

//*

// Set macros for use in routines below (*may move to sep. file - shared by diff generators)

#define SKEIN_MK_64(hi32,lo32)  ((lo32) + (((uint64_t) (hi32)) << 32))
#define SKEIN_KS_PARITY64         SKEIN_MK_64(0x1BD11BDA,0xA9FC1A22)
#define SKEIN_KS_PARITY32         0x1BD11BDA

//For u01 - divide a single result by range of unsigned integer to get value 0 to 1
//Variants for open/closed - 32/64 etc... - only subset here - see u01.h
#define R123_0x1p_32f               (1.f/4294967296.f)                        // for u01_closed_closed_32_24 int4 to s.p
#define R123_0x1p_32                (1./4294967296.)                          // for 32_53 CO,OC,OO int4 to d.p
#define R123_0x1p_64                (1./(4294967296.*4294967296.))            // for u01_closed_closed_64_53 int8 to d.p


// u01_closed_closed_32_24 is i*R123_0x1p_32f
// u01_closed_open_32_53   is i*R123_0x1p_32

// u01_closed_open_32_53   is i*R123_0x1p_32
// u01_open_closed_32_53   is (1.+i)*R123_0x1p_32
// u01_open_open_32_53     is (0.5+i)*R123_0x1p_32

// u01_closed_closed_64_53 is i*R123_0x1p_64



// Original contains some macro definitions

// If use with features files - which define the declaration macros eg. R123_CUDA_DEVICE
//------------------
// R123_CUDA_DEVICE R123_STATIC_INLINE R123_FORCE_INLINE(uint32_t RotL_32(uint32_t x, unsigned int N));
// R123_CUDA_DEVICE R123_STATIC_INLINE uint32_t RotL_32(uint32_t x, unsigned int N)
// {
//     return (x << (N & 31)) | (x >> ((32-N) & 31));
// }
// R123_CUDA_DEVICE R123_STATIC_INLINE R123_FORCE_INLINE(uint32_t RotL_32(uint32_t x, unsigned int N));
// R123_CUDA_DEVICE R123_STATIC_INLINE uint32_t RotL_32(uint32_t x, unsigned int N)
// {
//     return (x << (N & 31)) | (x >> ((32-N) & 31));
// }
//------------------

static inline uint32_t RotL_32(uint32_t x, unsigned int N)
{
    return (x << (N & 31)) | (x >> ((32-N) & 31));
}

// Rotation constants:
enum r123_enum_threefry32x4 {
    // Output from skein_rot_search: (srs-B128-X5000.out)
    // Random seed = 1. BlockSize = 64 bits. sampleCnt =  1024. rounds =  8, minHW_or=28
    // Start: Mon Aug 24 22:41:36 2009
    // ...
    // rMin = 0.472. #0A4B[*33] [CRC=DD1ECE0F. hw_OR=31. cnt=16384. blkSize= 128].format 
    R_32x4_0_0=10, R_32x4_0_1=26,
    R_32x4_1_0=11, R_32x4_1_1=21,
    R_32x4_2_0=13, R_32x4_2_1=27,
    R_32x4_3_0=23, R_32x4_3_1= 5,
    R_32x4_4_0= 6, R_32x4_4_1=20,
    R_32x4_5_0=17, R_32x4_5_1=11,
    R_32x4_6_0=25, R_32x4_6_1=10,
    R_32x4_7_0=18, R_32x4_7_1=20

};

enum r123_enum_threefry_wcnt {
    WCNT2=2,
    WCNT4=4
};


//*/

// ------------------------------------------------------------------------------------------------/
// ------------- End standalome version - comment out if using threefry.h etc ---------------------/
// -------------------------------------------------------------------------------------------------


//=========================================================================================================

//Utility call
//Returns parameter/macro values for vector and array sizes
//Plan to move to header file along with macro definitions
//int test_vec_sizes_() {
int v123_get_vec_sizes() {
  printf("VECTOR_LENGTH_BYTES %d\n",VECTOR_LENGTH_BYTES);
  printf("NUM_VALS_32 %d\n",NUM_VALS_32);
  printf("NUM_VALS_64 %d\n",NUM_VALS_64);  
  return 0;
}


//=========================================================================================================
// Threefry functions
//=========================================================================================================

// Note: These may be renamed in future - eg. v123_threefry4x32r_ss_fix

// Single stream version (with constant number of sets - NUM_VALS_32): 
// Provide one ctr-set and one key-set like scalar version. ctr1 is incremented
// across vector dimesion with supplied increment.
// An increment of 1 is as good as any for quality of random numbers
// Values 2/3/4 of each ctr/key set are replicated across vector dimension.
// Four vectors of random number are returned.
int threefry4x32f_multi_ss_fix(unsigned int* CTR4,unsigned int* KEY4, int INCREMENT, 
                               double* RAN1, double* RAN2, double* RAN3, double* RAN4) {
  
    int ivec;    
    int ictr;  // For iteration over ctk/keys 1,2,3,4
    
    __attribute__((aligned(VECTOR_LENGTH_BYTES))) uint32_t ks[4+1]; //Assuming one set of keys
              
    //Use 1D arrays for aligned accesses when supply vector size
    __attribute__((aligned(VECTOR_LENGTH_BYTES))) uint32_t X0[NUM_VALS_32];
    __attribute__((aligned(VECTOR_LENGTH_BYTES))) uint32_t X1[NUM_VALS_32];
    __attribute__((aligned(VECTOR_LENGTH_BYTES))) uint32_t X2[NUM_VALS_32];
    __attribute__((aligned(VECTOR_LENGTH_BYTES))) uint32_t X3[NUM_VALS_32];
    
    ks[4] =  SKEIN_KS_PARITY32;

    //Assumes only first counter varies

    for (ivec=0;ivec <NUM_VALS_32 ; ivec++) {
      X0[ivec]  = CTR4[0] + INCREMENT*ivec; // increment counter
      X1[ivec]  = CTR4[1];
      X2[ivec]  = CTR4[2];
      X3[ivec]  = CTR4[3];               
    }
    

    //loop over 4 vals (not vector dimension)
    for (ictr=0;ictr < 4; ictr++) {
      ks[ictr] = KEY4[ictr];
      ks[4] ^= KEY4[ictr]; 
    }

    //loop over vector length
    #pragma omp simd aligned(X0,X1,X2,X3,ks)
    for (ivec=0;ivec < NUM_VALS_32; ivec++) {

      X0[ivec] += ks[0]; X1[ivec] += ks[1]; X2[ivec] += ks[2]; X3[ivec] += ks[3];


      X0[ivec] += X1[ivec]; X1[ivec] = RotL_32(X1[ivec],R_32x4_0_0); X1[ivec] ^= X0[ivec]; 
      X2[ivec] += X3[ivec]; X3[ivec] = RotL_32(X3[ivec],R_32x4_0_1); X3[ivec] ^= X2[ivec]; 

      X0[ivec] += X3[ivec]; X3[ivec] = RotL_32(X3[ivec],R_32x4_1_0); X3[ivec] ^= X0[ivec]; 
      X2[ivec] += X1[ivec]; X1[ivec] = RotL_32(X1[ivec],R_32x4_1_1); X1[ivec] ^= X2[ivec]; 

      X0[ivec] += X1[ivec]; X1[ivec] = RotL_32(X1[ivec],R_32x4_2_0); X1[ivec] ^= X0[ivec]; 
      X2[ivec] += X3[ivec]; X3[ivec] = RotL_32(X3[ivec],R_32x4_2_1); X3[ivec] ^= X2[ivec]; 

      X0[ivec] += X3[ivec]; X3[ivec] = RotL_32(X3[ivec],R_32x4_3_0); X3[ivec] ^= X0[ivec]; 
      X2[ivec] += X1[ivec]; X1[ivec] = RotL_32(X1[ivec],R_32x4_3_1); X1[ivec] ^= X2[ivec]; 

      /* InjectKey(r=1) */                                            
      X0[ivec] += ks[1]; X1[ivec] += ks[2]; X2[ivec] += ks[3]; X3[ivec] += ks[4]; 
      X3[ivec] += 1;     /* X[WCNT4-1] += r  */                 


      X0[ivec] += X1[ivec]; X1[ivec] = RotL_32(X1[ivec],R_32x4_4_0); X1[ivec] ^= X0[ivec]; 
      X2[ivec] += X3[ivec]; X3[ivec] = RotL_32(X3[ivec],R_32x4_4_1); X3[ivec] ^= X2[ivec]; 

      X0[ivec] += X3[ivec]; X3[ivec] = RotL_32(X3[ivec],R_32x4_5_0); X3[ivec] ^= X0[ivec]; 
      X2[ivec] += X1[ivec]; X1[ivec] = RotL_32(X1[ivec],R_32x4_5_1); X1[ivec] ^= X2[ivec]; 

      X0[ivec] += X1[ivec]; X1[ivec] = RotL_32(X1[ivec],R_32x4_6_0); X1[ivec] ^= X0[ivec]; 
      X2[ivec] += X3[ivec]; X3[ivec] = RotL_32(X3[ivec],R_32x4_6_1); X3[ivec] ^= X2[ivec]; 

      X0[ivec] += X3[ivec]; X3[ivec] = RotL_32(X3[ivec],R_32x4_7_0); X3[ivec] ^= X0[ivec]; 
      X2[ivec] += X1[ivec]; X1[ivec] = RotL_32(X1[ivec],R_32x4_7_1); X1[ivec] ^= X2[ivec]; 

      /* InjectKey(r=2) */                                            
      X0[ivec] += ks[2]; X1[ivec] += ks[3]; X2[ivec] += ks[4]; X3[ivec] += ks[0]; 
      X3[ivec] += 2;     /* X[WCNT4-1] += r  */                 


      X0[ivec] += X1[ivec]; X1[ivec] = RotL_32(X1[ivec],R_32x4_0_0); X1[ivec] ^= X0[ivec]; 
      X2[ivec] += X3[ivec]; X3[ivec] = RotL_32(X3[ivec],R_32x4_0_1); X3[ivec] ^= X2[ivec]; 

      X0[ivec] += X3[ivec]; X3[ivec] = RotL_32(X3[ivec],R_32x4_1_0); X3[ivec] ^= X0[ivec]; 
      X2[ivec] += X1[ivec]; X1[ivec] = RotL_32(X1[ivec],R_32x4_1_1); X1[ivec] ^= X2[ivec]; 

      X0[ivec] += X1[ivec]; X1[ivec] = RotL_32(X1[ivec],R_32x4_2_0); X1[ivec] ^= X0[ivec]; 
      X2[ivec] += X3[ivec]; X3[ivec] = RotL_32(X3[ivec],R_32x4_2_1); X3[ivec] ^= X2[ivec]; 

      X0[ivec] += X3[ivec]; X3[ivec] = RotL_32(X3[ivec],R_32x4_3_0); X3[ivec] ^= X0[ivec]; 
      X2[ivec] += X1[ivec]; X1[ivec] = RotL_32(X1[ivec],R_32x4_3_1); X1[ivec] ^= X2[ivec]; 

      /* InjectKey(r=3) */                                            
      X0[ivec] += ks[3]; X1[ivec] += ks[4]; X2[ivec] += ks[0]; X3[ivec] += ks[1]; 
      X3[ivec] += 3;     /* X[WCNT4-1] += r  */                 


      X0[ivec] += X1[ivec]; X1[ivec] = RotL_32(X1[ivec],R_32x4_4_0); X1[ivec] ^= X0[ivec]; 
      X2[ivec] += X3[ivec]; X3[ivec] = RotL_32(X3[ivec],R_32x4_4_1); X3[ivec] ^= X2[ivec]; 

      X0[ivec] += X3[ivec]; X3[ivec] = RotL_32(X3[ivec],R_32x4_5_0); X3[ivec] ^= X0[ivec]; 
      X2[ivec] += X1[ivec]; X1[ivec] = RotL_32(X1[ivec],R_32x4_5_1); X1[ivec] ^= X2[ivec]; 

      X0[ivec] += X1[ivec]; X1[ivec] = RotL_32(X1[ivec],R_32x4_6_0); X1[ivec] ^= X0[ivec]; 
      X2[ivec] += X3[ivec]; X3[ivec] = RotL_32(X3[ivec],R_32x4_6_1); X3[ivec] ^= X2[ivec]; 

      X0[ivec] += X3[ivec]; X3[ivec] = RotL_32(X3[ivec],R_32x4_7_0); X3[ivec] ^= X0[ivec]; 
      X2[ivec] += X1[ivec]; X1[ivec] = RotL_32(X1[ivec],R_32x4_7_1); X1[ivec] ^= X2[ivec]; 

      /* InjectKey(r=1) */                                            
      X0[ivec] += ks[4]; X1[ivec] += ks[0]; X2[ivec] += ks[1]; X3[ivec] += ks[2]; 
      X3[ivec] += 4;     /* X[WCNT4-1] += r  */                 


      X0[ivec] += X1[ivec]; X1[ivec] = RotL_32(X1[ivec],R_32x4_0_0); X1[ivec] ^= X0[ivec]; 
      X2[ivec] += X3[ivec]; X3[ivec] = RotL_32(X3[ivec],R_32x4_0_1); X3[ivec] ^= X2[ivec]; 

      X0[ivec] += X3[ivec]; X3[ivec] = RotL_32(X3[ivec],R_32x4_1_0); X3[ivec] ^= X0[ivec]; 
      X2[ivec] += X1[ivec]; X1[ivec] = RotL_32(X1[ivec],R_32x4_1_1); X1[ivec] ^= X2[ivec]; 

      X0[ivec] += X1[ivec]; X1[ivec] = RotL_32(X1[ivec],R_32x4_2_0); X1[ivec] ^= X0[ivec]; 
      X2[ivec] += X3[ivec]; X3[ivec] = RotL_32(X3[ivec],R_32x4_2_1); X3[ivec] ^= X2[ivec]; 

      X0[ivec] += X3[ivec]; X3[ivec] = RotL_32(X3[ivec],R_32x4_3_0); X3[ivec] ^= X0[ivec]; 
      X2[ivec] += X1[ivec]; X1[ivec] = RotL_32(X1[ivec],R_32x4_3_1); X1[ivec] ^= X2[ivec]; 

      /* InjectKey(r=1) */                                            
      X0[ivec] += ks[0]; X1[ivec] += ks[1]; X2[ivec] += ks[2]; X3[ivec] += ks[3]; 
      X3[ivec] += 5;     /* X[WCNT4-1] += r  */                 


      //Convert integers to doubles uniformly distributed between 0 and 1 inclusive/exlusive
      //u01_closed_open_32_53
      RAN1[ivec] = X0[ivec]*R123_0x1p_32;
      RAN2[ivec] = X1[ivec]*R123_0x1p_32;
      RAN3[ivec] = X2[ivec]*R123_0x1p_32;
      RAN4[ivec] = X3[ivec]*R123_0x1p_32;

    }
    

    return 0; 
}

//=========================================================================================================


// Single stream version (with number of sets provided - IN_NUM_VALS): 
// Provide one ctr-set and one key-set like scalar version. ctr1 is incremented
// across vector dimesion with supplied increment.
// An increment of 1 is as good as any for quality of random numbers
// Values 2/3/4 of each ctr/key set are replicated across vector dimension.
// Four vectors of random number are returned.
int threefry4x32f_multi_ss(unsigned int* CTR4,unsigned int* KEY4, int INCREMENT, 
                           double* RAN1, double* RAN2, double* RAN3, double* RAN4, const int IN_NUM_VALS) {
  
    int ivec;    
    int ictr;  // For iteration over ctk/keys 1,2,3,4
    
    __attribute__((aligned(VECTOR_LENGTH_BYTES))) uint32_t ks[4+1]; //Assuming one set of keys
              
    //Use 1D arrays for aligned accesses when supply vector size
    __attribute__((aligned(VECTOR_LENGTH_BYTES))) uint32_t X0[IN_NUM_VALS];
    __attribute__((aligned(VECTOR_LENGTH_BYTES))) uint32_t X1[IN_NUM_VALS];
    __attribute__((aligned(VECTOR_LENGTH_BYTES))) uint32_t X2[IN_NUM_VALS];
    __attribute__((aligned(VECTOR_LENGTH_BYTES))) uint32_t X3[IN_NUM_VALS];
    
    ks[4] =  SKEIN_KS_PARITY32;

    //Assumes only first counter varies

    for (ivec=0;ivec <IN_NUM_VALS ; ivec++) {
      X0[ivec]  = CTR4[0] + INCREMENT*ivec; // increment counter
      X1[ivec]  = CTR4[1];
      X2[ivec]  = CTR4[2];
      X3[ivec]  = CTR4[3];               
    }
    

    //loop over 4 vals (not vector dimension)
    for (ictr=0;ictr < 4; ictr++) {
      ks[ictr] = KEY4[ictr];
      ks[4] ^= KEY4[ictr]; 
    }

    //loop over vector length
    #pragma omp simd aligned(X0,X1,X2,X3,ks)
    for (ivec=0;ivec < IN_NUM_VALS; ivec++) {

      X0[ivec] += ks[0]; X1[ivec] += ks[1]; X2[ivec] += ks[2]; X3[ivec] += ks[3];


      X0[ivec] += X1[ivec]; X1[ivec] = RotL_32(X1[ivec],R_32x4_0_0); X1[ivec] ^= X0[ivec]; 
      X2[ivec] += X3[ivec]; X3[ivec] = RotL_32(X3[ivec],R_32x4_0_1); X3[ivec] ^= X2[ivec]; 

      X0[ivec] += X3[ivec]; X3[ivec] = RotL_32(X3[ivec],R_32x4_1_0); X3[ivec] ^= X0[ivec]; 
      X2[ivec] += X1[ivec]; X1[ivec] = RotL_32(X1[ivec],R_32x4_1_1); X1[ivec] ^= X2[ivec]; 

      X0[ivec] += X1[ivec]; X1[ivec] = RotL_32(X1[ivec],R_32x4_2_0); X1[ivec] ^= X0[ivec]; 
      X2[ivec] += X3[ivec]; X3[ivec] = RotL_32(X3[ivec],R_32x4_2_1); X3[ivec] ^= X2[ivec]; 

      X0[ivec] += X3[ivec]; X3[ivec] = RotL_32(X3[ivec],R_32x4_3_0); X3[ivec] ^= X0[ivec]; 
      X2[ivec] += X1[ivec]; X1[ivec] = RotL_32(X1[ivec],R_32x4_3_1); X1[ivec] ^= X2[ivec]; 

      /* InjectKey(r=1) */                                            
      X0[ivec] += ks[1]; X1[ivec] += ks[2]; X2[ivec] += ks[3]; X3[ivec] += ks[4]; 
      X3[ivec] += 1;     /* X[WCNT4-1] += r  */                 


      X0[ivec] += X1[ivec]; X1[ivec] = RotL_32(X1[ivec],R_32x4_4_0); X1[ivec] ^= X0[ivec]; 
      X2[ivec] += X3[ivec]; X3[ivec] = RotL_32(X3[ivec],R_32x4_4_1); X3[ivec] ^= X2[ivec]; 

      X0[ivec] += X3[ivec]; X3[ivec] = RotL_32(X3[ivec],R_32x4_5_0); X3[ivec] ^= X0[ivec]; 
      X2[ivec] += X1[ivec]; X1[ivec] = RotL_32(X1[ivec],R_32x4_5_1); X1[ivec] ^= X2[ivec]; 

      X0[ivec] += X1[ivec]; X1[ivec] = RotL_32(X1[ivec],R_32x4_6_0); X1[ivec] ^= X0[ivec]; 
      X2[ivec] += X3[ivec]; X3[ivec] = RotL_32(X3[ivec],R_32x4_6_1); X3[ivec] ^= X2[ivec]; 

      X0[ivec] += X3[ivec]; X3[ivec] = RotL_32(X3[ivec],R_32x4_7_0); X3[ivec] ^= X0[ivec]; 
      X2[ivec] += X1[ivec]; X1[ivec] = RotL_32(X1[ivec],R_32x4_7_1); X1[ivec] ^= X2[ivec]; 

      /* InjectKey(r=2) */                                            
      X0[ivec] += ks[2]; X1[ivec] += ks[3]; X2[ivec] += ks[4]; X3[ivec] += ks[0]; 
      X3[ivec] += 2;     /* X[WCNT4-1] += r  */                 


      X0[ivec] += X1[ivec]; X1[ivec] = RotL_32(X1[ivec],R_32x4_0_0); X1[ivec] ^= X0[ivec]; 
      X2[ivec] += X3[ivec]; X3[ivec] = RotL_32(X3[ivec],R_32x4_0_1); X3[ivec] ^= X2[ivec]; 

      X0[ivec] += X3[ivec]; X3[ivec] = RotL_32(X3[ivec],R_32x4_1_0); X3[ivec] ^= X0[ivec]; 
      X2[ivec] += X1[ivec]; X1[ivec] = RotL_32(X1[ivec],R_32x4_1_1); X1[ivec] ^= X2[ivec]; 

      X0[ivec] += X1[ivec]; X1[ivec] = RotL_32(X1[ivec],R_32x4_2_0); X1[ivec] ^= X0[ivec]; 
      X2[ivec] += X3[ivec]; X3[ivec] = RotL_32(X3[ivec],R_32x4_2_1); X3[ivec] ^= X2[ivec]; 

      X0[ivec] += X3[ivec]; X3[ivec] = RotL_32(X3[ivec],R_32x4_3_0); X3[ivec] ^= X0[ivec]; 
      X2[ivec] += X1[ivec]; X1[ivec] = RotL_32(X1[ivec],R_32x4_3_1); X1[ivec] ^= X2[ivec]; 

      /* InjectKey(r=3) */                                            
      X0[ivec] += ks[3]; X1[ivec] += ks[4]; X2[ivec] += ks[0]; X3[ivec] += ks[1]; 
      X3[ivec] += 3;     /* X[WCNT4-1] += r  */                 


      X0[ivec] += X1[ivec]; X1[ivec] = RotL_32(X1[ivec],R_32x4_4_0); X1[ivec] ^= X0[ivec]; 
      X2[ivec] += X3[ivec]; X3[ivec] = RotL_32(X3[ivec],R_32x4_4_1); X3[ivec] ^= X2[ivec]; 

      X0[ivec] += X3[ivec]; X3[ivec] = RotL_32(X3[ivec],R_32x4_5_0); X3[ivec] ^= X0[ivec]; 
      X2[ivec] += X1[ivec]; X1[ivec] = RotL_32(X1[ivec],R_32x4_5_1); X1[ivec] ^= X2[ivec]; 

      X0[ivec] += X1[ivec]; X1[ivec] = RotL_32(X1[ivec],R_32x4_6_0); X1[ivec] ^= X0[ivec]; 
      X2[ivec] += X3[ivec]; X3[ivec] = RotL_32(X3[ivec],R_32x4_6_1); X3[ivec] ^= X2[ivec]; 

      X0[ivec] += X3[ivec]; X3[ivec] = RotL_32(X3[ivec],R_32x4_7_0); X3[ivec] ^= X0[ivec]; 
      X2[ivec] += X1[ivec]; X1[ivec] = RotL_32(X1[ivec],R_32x4_7_1); X1[ivec] ^= X2[ivec]; 

      /* InjectKey(r=1) */                                            
      X0[ivec] += ks[4]; X1[ivec] += ks[0]; X2[ivec] += ks[1]; X3[ivec] += ks[2]; 
      X3[ivec] += 4;     /* X[WCNT4-1] += r  */                 


      X0[ivec] += X1[ivec]; X1[ivec] = RotL_32(X1[ivec],R_32x4_0_0); X1[ivec] ^= X0[ivec]; 
      X2[ivec] += X3[ivec]; X3[ivec] = RotL_32(X3[ivec],R_32x4_0_1); X3[ivec] ^= X2[ivec]; 

      X0[ivec] += X3[ivec]; X3[ivec] = RotL_32(X3[ivec],R_32x4_1_0); X3[ivec] ^= X0[ivec]; 
      X2[ivec] += X1[ivec]; X1[ivec] = RotL_32(X1[ivec],R_32x4_1_1); X1[ivec] ^= X2[ivec]; 

      X0[ivec] += X1[ivec]; X1[ivec] = RotL_32(X1[ivec],R_32x4_2_0); X1[ivec] ^= X0[ivec]; 
      X2[ivec] += X3[ivec]; X3[ivec] = RotL_32(X3[ivec],R_32x4_2_1); X3[ivec] ^= X2[ivec]; 

      X0[ivec] += X3[ivec]; X3[ivec] = RotL_32(X3[ivec],R_32x4_3_0); X3[ivec] ^= X0[ivec]; 
      X2[ivec] += X1[ivec]; X1[ivec] = RotL_32(X1[ivec],R_32x4_3_1); X1[ivec] ^= X2[ivec]; 

      /* InjectKey(r=1) */                                            
      X0[ivec] += ks[0]; X1[ivec] += ks[1]; X2[ivec] += ks[2]; X3[ivec] += ks[3]; 
      X3[ivec] += 5;     /* X[WCNT4-1] += r  */                 


      //Convert integers to doubles uniformly distributed between 0 and 1 inclusive/exlusive
      //u01_closed_open_32_53
      RAN1[ivec] = X0[ivec]*R123_0x1p_32;
      RAN2[ivec] = X1[ivec]*R123_0x1p_32;
      RAN3[ivec] = X2[ivec]*R123_0x1p_32;
      RAN4[ivec] = X3[ivec]*R123_0x1p_32;

    }
    

    return 0; 
}

//=========================================================================================================

// Interface 
// All = Supply all counters/keys
// Four vectors of random number are returned.      
int threefry4x32f_multi_ctrkey_all(unsigned int* CTR1,unsigned int* CTR2,unsigned int* CTR3,unsigned int* CTR4,  
                                   unsigned int* KEY1,unsigned int* KEY2,unsigned int* KEY3,unsigned int* KEY4, 
                                   double* RAN1, double* RAN2, double* RAN3, double* RAN4, const int IN_NUM_VALS) {
  
    int ivec;  // For iteration over vector dimesion
    //int ictr;  // For iteration over ctk/keys 1,2,3,4
          
    //Use 1D arrays for aligned accesses when supply vector size
    __attribute__((aligned(VECTOR_LENGTH_BYTES))) uint32_t ks0[IN_NUM_VALS];    
    __attribute__((aligned(VECTOR_LENGTH_BYTES))) uint32_t ks1[IN_NUM_VALS];
    __attribute__((aligned(VECTOR_LENGTH_BYTES))) uint32_t ks2[IN_NUM_VALS];
    __attribute__((aligned(VECTOR_LENGTH_BYTES))) uint32_t ks3[IN_NUM_VALS];
    __attribute__((aligned(VECTOR_LENGTH_BYTES))) uint32_t ks4[IN_NUM_VALS];

              
    //Use 1D arrays for aligned accesses when supply vector size
    __attribute__((aligned(VECTOR_LENGTH_BYTES))) uint32_t X0[IN_NUM_VALS];
    __attribute__((aligned(VECTOR_LENGTH_BYTES))) uint32_t X1[IN_NUM_VALS];
    __attribute__((aligned(VECTOR_LENGTH_BYTES))) uint32_t X2[IN_NUM_VALS];
    __attribute__((aligned(VECTOR_LENGTH_BYTES))) uint32_t X3[IN_NUM_VALS];


    //__attribute__((aligned(VECTOR_LENGTH_BYTES))) double buff[IN_NUM_VALS*4];
    

    for (ivec=0;ivec <IN_NUM_VALS ; ivec++) {
      X0[ivec]  = CTR1[ivec];
      X1[ivec]  = CTR2[ivec];
      X2[ivec]  = CTR3[ivec];
      X3[ivec]  = CTR4[ivec];       
      ks4[ivec] =  SKEIN_KS_PARITY32;                          
    }


    for (ivec=0;ivec <IN_NUM_VALS ; ivec++) {
        ks0[ivec] = KEY1[ivec];
        ks1[ivec] = KEY2[ivec];
        ks2[ivec] = KEY3[ivec];
        ks3[ivec] = KEY4[ivec];
        
        //could interleave or here
        ks4[ivec] ^= KEY1[ivec]; 
        ks4[ivec] ^= KEY2[ivec]; 
        ks4[ivec] ^= KEY3[ivec]; 
        ks4[ivec] ^= KEY4[ivec]; 
    }



    //Loop over vector length (or multiple thereof)
    #pragma omp simd aligned(X0,X1,X2,X3,ks0,ks1,ks2,ks3,ks4)
    for (ivec=0;ivec < IN_NUM_VALS; ivec++) {

      X0[ivec] += ks0[ivec]; X1[ivec] += ks1[ivec]; X2[ivec] += ks2[ivec]; X3[ivec] += ks3[ivec];


      X0[ivec] += X1[ivec]; X1[ivec] = RotL_32(X1[ivec],R_32x4_0_0); X1[ivec] ^= X0[ivec]; 
      X2[ivec] += X3[ivec]; X3[ivec] = RotL_32(X3[ivec],R_32x4_0_1); X3[ivec] ^= X2[ivec]; 

      X0[ivec] += X3[ivec]; X3[ivec] = RotL_32(X3[ivec],R_32x4_1_0); X3[ivec] ^= X0[ivec]; 
      X2[ivec] += X1[ivec]; X1[ivec] = RotL_32(X1[ivec],R_32x4_1_1); X1[ivec] ^= X2[ivec]; 

      X0[ivec] += X1[ivec]; X1[ivec] = RotL_32(X1[ivec],R_32x4_2_0); X1[ivec] ^= X0[ivec]; 
      X2[ivec] += X3[ivec]; X3[ivec] = RotL_32(X3[ivec],R_32x4_2_1); X3[ivec] ^= X2[ivec]; 

      X0[ivec] += X3[ivec]; X3[ivec] = RotL_32(X3[ivec],R_32x4_3_0); X3[ivec] ^= X0[ivec]; 
      X2[ivec] += X1[ivec]; X1[ivec] = RotL_32(X1[ivec],R_32x4_3_1); X1[ivec] ^= X2[ivec]; 

      /* InjectKey(r=1) */                                            
      X0[ivec] += ks1[ivec]; X1[ivec] += ks2[ivec]; X2[ivec] += ks3[ivec]; X3[ivec] += ks4[ivec]; 
      X3[ivec] += 1;     /* X[WCNT4-1] += r  */                 


      X0[ivec] += X1[ivec]; X1[ivec] = RotL_32(X1[ivec],R_32x4_4_0); X1[ivec] ^= X0[ivec]; 
      X2[ivec] += X3[ivec]; X3[ivec] = RotL_32(X3[ivec],R_32x4_4_1); X3[ivec] ^= X2[ivec]; 

      X0[ivec] += X3[ivec]; X3[ivec] = RotL_32(X3[ivec],R_32x4_5_0); X3[ivec] ^= X0[ivec]; 
      X2[ivec] += X1[ivec]; X1[ivec] = RotL_32(X1[ivec],R_32x4_5_1); X1[ivec] ^= X2[ivec]; 

      X0[ivec] += X1[ivec]; X1[ivec] = RotL_32(X1[ivec],R_32x4_6_0); X1[ivec] ^= X0[ivec]; 
      X2[ivec] += X3[ivec]; X3[ivec] = RotL_32(X3[ivec],R_32x4_6_1); X3[ivec] ^= X2[ivec]; 

      X0[ivec] += X3[ivec]; X3[ivec] = RotL_32(X3[ivec],R_32x4_7_0); X3[ivec] ^= X0[ivec]; 
      X2[ivec] += X1[ivec]; X1[ivec] = RotL_32(X1[ivec],R_32x4_7_1); X1[ivec] ^= X2[ivec]; 

      /* InjectKey(r=2) */                                            
      X0[ivec] += ks2[ivec]; X1[ivec] += ks3[ivec]; X2[ivec] += ks4[ivec]; X3[ivec] += ks0[ivec]; 
      X3[ivec] += 2;     /* X[WCNT4-1] += r  */                 


      X0[ivec] += X1[ivec]; X1[ivec] = RotL_32(X1[ivec],R_32x4_0_0); X1[ivec] ^= X0[ivec]; 
      X2[ivec] += X3[ivec]; X3[ivec] = RotL_32(X3[ivec],R_32x4_0_1); X3[ivec] ^= X2[ivec]; 

      X0[ivec] += X3[ivec]; X3[ivec] = RotL_32(X3[ivec],R_32x4_1_0); X3[ivec] ^= X0[ivec]; 
      X2[ivec] += X1[ivec]; X1[ivec] = RotL_32(X1[ivec],R_32x4_1_1); X1[ivec] ^= X2[ivec]; 

      X0[ivec] += X1[ivec]; X1[ivec] = RotL_32(X1[ivec],R_32x4_2_0); X1[ivec] ^= X0[ivec]; 
      X2[ivec] += X3[ivec]; X3[ivec] = RotL_32(X3[ivec],R_32x4_2_1); X3[ivec] ^= X2[ivec]; 

      X0[ivec] += X3[ivec]; X3[ivec] = RotL_32(X3[ivec],R_32x4_3_0); X3[ivec] ^= X0[ivec]; 
      X2[ivec] += X1[ivec]; X1[ivec] = RotL_32(X1[ivec],R_32x4_3_1); X1[ivec] ^= X2[ivec]; 

      /* InjectKey(r=3) */                                            
      X0[ivec] += ks3[ivec]; X1[ivec] += ks4[ivec]; X2[ivec] += ks0[ivec]; X3[ivec] += ks1[ivec]; 
      X3[ivec] += 3;     /* X[WCNT4-1] += r  */                 


      X0[ivec] += X1[ivec]; X1[ivec] = RotL_32(X1[ivec],R_32x4_4_0); X1[ivec] ^= X0[ivec]; 
      X2[ivec] += X3[ivec]; X3[ivec] = RotL_32(X3[ivec],R_32x4_4_1); X3[ivec] ^= X2[ivec]; 

      X0[ivec] += X3[ivec]; X3[ivec] = RotL_32(X3[ivec],R_32x4_5_0); X3[ivec] ^= X0[ivec]; 
      X2[ivec] += X1[ivec]; X1[ivec] = RotL_32(X1[ivec],R_32x4_5_1); X1[ivec] ^= X2[ivec]; 

      X0[ivec] += X1[ivec]; X1[ivec] = RotL_32(X1[ivec],R_32x4_6_0); X1[ivec] ^= X0[ivec]; 
      X2[ivec] += X3[ivec]; X3[ivec] = RotL_32(X3[ivec],R_32x4_6_1); X3[ivec] ^= X2[ivec]; 

      X0[ivec] += X3[ivec]; X3[ivec] = RotL_32(X3[ivec],R_32x4_7_0); X3[ivec] ^= X0[ivec]; 
      X2[ivec] += X1[ivec]; X1[ivec] = RotL_32(X1[ivec],R_32x4_7_1); X1[ivec] ^= X2[ivec]; 

      /* InjectKey(r=1) */                                            
      X0[ivec] += ks4[ivec]; X1[ivec] += ks0[ivec]; X2[ivec] += ks1[ivec]; X3[ivec] += ks2[ivec]; 
      X3[ivec] += 4;     /* X[WCNT4-1] += r  */                 


      X0[ivec] += X1[ivec]; X1[ivec] = RotL_32(X1[ivec],R_32x4_0_0); X1[ivec] ^= X0[ivec]; 
      X2[ivec] += X3[ivec]; X3[ivec] = RotL_32(X3[ivec],R_32x4_0_1); X3[ivec] ^= X2[ivec]; 

      X0[ivec] += X3[ivec]; X3[ivec] = RotL_32(X3[ivec],R_32x4_1_0); X3[ivec] ^= X0[ivec]; 
      X2[ivec] += X1[ivec]; X1[ivec] = RotL_32(X1[ivec],R_32x4_1_1); X1[ivec] ^= X2[ivec]; 

      X0[ivec] += X1[ivec]; X1[ivec] = RotL_32(X1[ivec],R_32x4_2_0); X1[ivec] ^= X0[ivec]; 
      X2[ivec] += X3[ivec]; X3[ivec] = RotL_32(X3[ivec],R_32x4_2_1); X3[ivec] ^= X2[ivec]; 

      X0[ivec] += X3[ivec]; X3[ivec] = RotL_32(X3[ivec],R_32x4_3_0); X3[ivec] ^= X0[ivec]; 
      X2[ivec] += X1[ivec]; X1[ivec] = RotL_32(X1[ivec],R_32x4_3_1); X1[ivec] ^= X2[ivec]; 

      /* InjectKey(r=1) */                                            
      X0[ivec] += ks0[ivec]; X1[ivec] += ks1[ivec]; X2[ivec] += ks2[ivec]; X3[ivec] += ks3[ivec]; 
      X3[ivec] += 5;     /* X[WCNT4-1] += r  */                 


      //Convert integers to doubles uniformly distributed between 0 and 1 inclusive/exlusive
      //u01_closed_open_32_53
      RAN1[ivec] = X0[ivec]*R123_0x1p_32;
      RAN2[ivec] = X1[ivec]*R123_0x1p_32;
      RAN3[ivec] = X2[ivec]*R123_0x1p_32;
      RAN4[ivec] = X3[ivec]*R123_0x1p_32;
      

    }
    
    return 0; 
}

