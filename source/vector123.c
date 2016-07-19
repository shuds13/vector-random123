//--------------------------------------------------------------------------------------------------
// S. Hudson: Vector123
// Standard C vectorizable implementation of Random123 functions.
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

//--------------------------------------------------------------------------------------------------
// S. Hudson: Standalone version - threefry vectorised.
//--------------------------------------------------------------------------------------------------
// At time of writing - restricted interface/options - this is not a full random123 port.
// Alternative is to use this with Random123 library and will pick up the compiler feature files etc.
// which sets some macros  (eg.directives for inlining)
// eg #define R123_FORCE_INLINE(decl) decl __attribute__((__always_inline__))
// Also if use threefry.h  - comment out section indicated below as will pick up from there.
// -> so far I've found makes no difference on platforms tested for threefry.
// Not tried this on CUDA platforms..
//--------------------------------------------------------------------------------------------------


#include <stdio.h>

// Allows use of short forms like uint32_t
#include <stdint.h>


//Random123 - To pick up macros from R123 - uncomment
//#include <threefry.h>
//#include <philox.h>

//#include <u01.h>
//OR
//#include <u01fixedpt.h>


//Try and automatic where can for vector length
//Note a sep. scalar/autovec version may handle automatically -though careful that its approp for integers
// - ie. avx1 only avx-128 for integers!! That said I guess if too long with break into two or whatever.

// Generally I dont think should matter if too big - as long as fits with calling routines and multiple of correct size.
// VECTOR_LENGTH_BYTES // For alignment if using - should not matter if too big (avx256=32 avx512=64) 
// NUM_VALS_32         // For multi-val version - should be  vector length in words (or multiple of) - LOOP SIZE
// NUM_VALS_64         // For multi-val version - should be  vector length in words (or multiple of) - LOOP SIZE

//This works but not too helpful as has to match calling routine - unless do that same way - make same.
// #if defined(__ICC)
// #if   defined(__AVX512F__)
// #define VECTOR_LENGTH_BYTES 64
// #define NUM_VALS_32 16
// #define NUM_VALS_64 8
// #elif   defined(__AVX2__)
// #define VECTOR_LENGTH_BYTES 32
// #define NUM_VALS_32 8
// #define NUM_VALS_64 4
// #elif   defined(__AVX__)
// #define VECTOR_LENGTH_BYTES 32
// #define NUM_VALS_32 4 //AVX1 - integers limited to AVX-128
// #define NUM_VALS_64 2
// #elif   defined(__SSE__)
// #define VECTOR_LENGTH_BYTES 32
// #define NUM_VALS_32 4 //AVX1 - integers limited to AVX-128
// #define NUM_VALS_64 2
// #else
// //ICC DEFAULTS
// #define VECTOR_LENGTH_BYTES 32
// #define NUM_VALS_32 8
// #define NUM_VALS_64 4
// #endif
// #else
// //DEFAULTS
// #define VECTOR_LENGTH_BYTES 32
// #define NUM_VALS_32 8
// #define NUM_VALS_64 4
// #endif

//eg. AVX-512
//#define VECTOR_LENGTH_BYTES 64
//#define NUM_VALS_32 16
//#define NUM_VALS_64 8

//eg. AVX2
// #define VECTOR_LENGTH_BYTES 32
// #define NUM_VALS_32 8
// #define NUM_VALS_64 4

//eg. AVX/SSE - may keep length 32 for AVX1
 #define VECTOR_LENGTH_BYTES 32
// //#define VECTOR_LENGTH_BYTES 16
 #define NUM_VALS_32 4
 #define NUM_VALS_64 2

// ------------------------------------------------------------------------------------------------/
// ------------- For standalome version - comment out if using threefry.h etc ---------------------/
// ------------------------------------------------------------------------------------------------/

//If want to do a formal vectorised random123 including all options
// the prob. more sensible to "use" than write all options out here.
//Just I like the idea of not having to have a sub-dir - so keep this to options 
//I'm likely to use.

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

static inline uint64_t RotL_64(uint64_t x, unsigned int N)
{
    return (x << (N & 63)) | (x >> ((64-N) & 63));
}

/* Rotation constants: */
enum r123_enum_threefry32x4 {
    /* Output from skein_rot_search: (srs-B128-X5000.out)
    // Random seed = 1. BlockSize = 64 bits. sampleCnt =  1024. rounds =  8, minHW_or=28
    // Start: Mon Aug 24 22:41:36 2009
    // ...
    // rMin = 0.472. #0A4B[*33] [CRC=DD1ECE0F. hw_OR=31. cnt=16384. blkSize= 128].format    */
    R_32x4_0_0=10, R_32x4_0_1=26,
    R_32x4_1_0=11, R_32x4_1_1=21,
    R_32x4_2_0=13, R_32x4_2_1=27,
    R_32x4_3_0=23, R_32x4_3_1= 5,
    R_32x4_4_0= 6, R_32x4_4_1=20,
    R_32x4_5_0=17, R_32x4_5_1=11,
    R_32x4_6_0=25, R_32x4_6_1=10,
    R_32x4_7_0=18, R_32x4_7_1=20

};

enum r123_enum_threefry64x4 {
    /* These are the R_256 constants from the Threefish reference sources
       with names changed to R_64x4... */
    R_64x4_0_0=14, R_64x4_0_1=16,
    R_64x4_1_0=52, R_64x4_1_1=57,
    R_64x4_2_0=23, R_64x4_2_1=40,
    R_64x4_3_0= 5, R_64x4_3_1=37,
    R_64x4_4_0=25, R_64x4_4_1=33,
    R_64x4_5_0=46, R_64x4_5_1=12,
    R_64x4_6_0=58, R_64x4_6_1=22,
    R_64x4_7_0=32, R_64x4_7_1=32
};

enum r123_enum_threefry64x2 {
    /*
    // Output from skein_rot_search: (srs64_B64-X1000)
    // Random seed = 1. BlockSize = 128 bits. sampleCnt =  1024. rounds =  8, minHW_or=57
    // Start: Tue Mar  1 10:07:48 2011
    // rMin = 0.136. #0325[*15] [CRC=455A682F. hw_OR=64. cnt=16384. blkSize= 128].format   
    */
    R_64x2_0_0=16,
    R_64x2_1_0=42,
    R_64x2_2_0=12,
    R_64x2_3_0=31,
    R_64x2_4_0=16,
    R_64x2_5_0=32,
    R_64x2_6_0=24,
    R_64x2_7_0=21
};

enum r123_enum_threefry32x2 {
    /* Output from skein_rot_search (srs32x2-X5000.out)
    // Random seed = 1. BlockSize = 64 bits. sampleCnt =  1024. rounds =  8, minHW_or=28
    // Start: Tue Jul 12 11:11:33 2011
    // rMin = 0.334. #0206[*07] [CRC=1D9765C0. hw_OR=32. cnt=16384. blkSize=  64].format   */
    R_32x2_0_0=13,
    R_32x2_1_0=15,
    R_32x2_2_0=26,
    R_32x2_3_0= 6,
    R_32x2_4_0=17,
    R_32x2_5_0=29,
    R_32x2_6_0=16,
    R_32x2_7_0=24

    /* 4 rounds: minHW =  4  [  4  4  4  4 ]
    // 5 rounds: minHW =  6  [  6  8  6  8 ]
    // 6 rounds: minHW =  9  [  9 12  9 12 ]
    // 7 rounds: minHW = 16  [ 16 24 16 24 ]
    // 8 rounds: minHW = 32  [ 32 32 32 32 ]
    // 9 rounds: minHW = 32  [ 32 32 32 32 ]
    //10 rounds: minHW = 32  [ 32 32 32 32 ]
    //11 rounds: minHW = 32  [ 32 32 32 32 ] */
    };

enum r123_enum_threefry_wcnt {
    WCNT2=2,
    WCNT4=4
};

// ------------------------------------------------------------------------------------------------/
// ------------- End standalome version - comment out if using threefry.h etc ---------------------/
// -------------------------------------------------------------------------------------------------

//=========================================================================================================

//Utility call
int test_vec_sizes_() {
  //Default values vector sizes NUM_VALS_32 and NUM_VALS_64
  //In some functions these can be overwritten with provided value
  printf("--------------------------------------------------------------------- \n");
  printf("Default vector lengths from vector123.c: Note some routines overwrite \n");
  printf("VECTOR_LENGTH_BYTES %d\n",VECTOR_LENGTH_BYTES);
  printf("NUM_VALS_32 %d\n",NUM_VALS_32);
  printf("NUM_VALS_64 %d\n",NUM_VALS_64);  
  return 0;
}


//=========================================================================================================
//MULTI CALL VERSIONS
//=========================================================================================================

// ------------------------------------ threefry4x32f_multi_ss -----------------------------------------

//extern "C"
//1D array version - with fixed length
int threefry4x32f_multi_ss_fix(unsigned int* CTR4,unsigned int* KEY4, int INCREMENT, double* buff) {
  
    int ivec;    
    int ictr;  // For iteration over ctk/keys 1,2,3,4
    
    //uint32_t ks[4+1]; //Assuming one set of keys
    __attribute__((aligned(VECTOR_LENGTH_BYTES))) uint32_t ks[4+1]; //Assuming one set of keys
              
    //uint32_t X[4][NUM_VALS_32];
    __attribute__((aligned(VECTOR_LENGTH_BYTES))) uint32_t X[4][NUM_VALS_32];

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
      buff[ivec*4+0] = X0[ivec]*R123_0x1p_32;
      buff[ivec*4+1] = X1[ivec]*R123_0x1p_32;
      buff[ivec*4+2] = X2[ivec]*R123_0x1p_32;
      buff[ivec*4+3] = X3[ivec]*R123_0x1p_32;

    }
    

    return 0; 
}

//1D array version - with provided length
int threefry4x32f_multi_ss(unsigned int* CTR4,unsigned int* KEY4, int INCREMENT, double* buff, const int IN_NUM_VALS) {
  
    int ivec;    
    int ictr;  // For iteration over ctk/keys 1,2,3,4
    
    //uint32_t ks[4+1]; //Assuming one set of keys
    __attribute__((aligned(VECTOR_LENGTH_BYTES))) uint32_t ks[4+1]; //Assuming one set of keys
              
    //uint32_t X[4][IN_NUM_VALS];
    __attribute__((aligned(VECTOR_LENGTH_BYTES))) uint32_t X[4][IN_NUM_VALS];

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
      buff[ivec*4+0] = X0[ivec]*R123_0x1p_32;
      buff[ivec*4+1] = X1[ivec]*R123_0x1p_32;
      buff[ivec*4+2] = X2[ivec]*R123_0x1p_32;
      buff[ivec*4+3] = X3[ivec]*R123_0x1p_32;

    }
    

    return 0; 
}



// ------------------------------------ Other lengths - old needs updating ----------------------------------

//old
int threefry4x64f_multi(unsigned int* CTR4,unsigned int* KEY4, double* buff) {
  
    int i;    
    
    //uint64_t ks[4+1]; //Assuming one set of keys
    __attribute__((aligned(VECTOR_LENGTH_BYTES))) uint64_t ks[4+1]; //Assuming one set of keys
              
    //uint64_t X[4][NUM_VALS_64];
    __attribute__((aligned(VECTOR_LENGTH_BYTES))) uint64_t X[4][NUM_VALS_64];
    
    ks[4] =  SKEIN_KS_PARITY64;

    //Assumes only first counter varies......... big assumption 
    
    for (i=0;i <NUM_VALS_64 ; i++) {
      X[0][i]  = *CTR4 + 4*i; // increment counter
      X[1][i]  = CTR4[1];
      X[2][i]  = CTR4[2];
      X[3][i]  = CTR4[3];               
    }
    

    //loop over 4 vals (not vector dimension)
    for (i=0;i < 4; i++) {
      ks[i] = KEY4[i];
      ks[4] ^= KEY4[i]; 
    }

    // ks[4] ^= KEY4[0];
    // ks[4] ^= KEY4[1];
    // ks[4] ^= KEY4[2];
    // ks[4] ^= KEY4[3];

    //loop over vector length
    //#pragma simd
    
    for (i=0;i < NUM_VALS_64; i++) {

      X[0][i] += ks[0]; X[1][i] += ks[1]; X[2][i] += ks[2]; X[3][i] += ks[3];


      X[0][i] += X[1][i]; X[1][i] = RotL_64(X[1][i],R_64x4_0_0); X[1][i] ^= X[0][i]; 
      X[2][i] += X[3][i]; X[3][i] = RotL_64(X[3][i],R_64x4_0_1); X[3][i] ^= X[2][i]; 

      X[0][i] += X[3][i]; X[3][i] = RotL_64(X[3][i],R_64x4_1_0); X[3][i] ^= X[0][i]; 
      X[2][i] += X[1][i]; X[1][i] = RotL_64(X[1][i],R_64x4_1_1); X[1][i] ^= X[2][i]; 

      X[0][i] += X[1][i]; X[1][i] = RotL_64(X[1][i],R_64x4_2_0); X[1][i] ^= X[0][i]; 
      X[2][i] += X[3][i]; X[3][i] = RotL_64(X[3][i],R_64x4_2_1); X[3][i] ^= X[2][i]; 

      X[0][i] += X[3][i]; X[3][i] = RotL_64(X[3][i],R_64x4_3_0); X[3][i] ^= X[0][i]; 
      X[2][i] += X[1][i]; X[1][i] = RotL_64(X[1][i],R_64x4_3_1); X[1][i] ^= X[2][i]; 

      /* InjectKey(r=1) */                                            
      X[0][i] += ks[1]; X[1][i] += ks[2]; X[2][i] += ks[3]; X[3][i] += ks[4]; 
      X[4-1][i] += 1;     /* X[WCNT4-1] += r  */                 


      X[0][i] += X[1][i]; X[1][i] = RotL_64(X[1][i],R_64x4_4_0); X[1][i] ^= X[0][i]; 
      X[2][i] += X[3][i]; X[3][i] = RotL_64(X[3][i],R_64x4_4_1); X[3][i] ^= X[2][i]; 

      X[0][i] += X[3][i]; X[3][i] = RotL_64(X[3][i],R_64x4_5_0); X[3][i] ^= X[0][i]; 
      X[2][i] += X[1][i]; X[1][i] = RotL_64(X[1][i],R_64x4_5_1); X[1][i] ^= X[2][i]; 

      X[0][i] += X[1][i]; X[1][i] = RotL_64(X[1][i],R_64x4_6_0); X[1][i] ^= X[0][i]; 
      X[2][i] += X[3][i]; X[3][i] = RotL_64(X[3][i],R_64x4_6_1); X[3][i] ^= X[2][i]; 

      X[0][i] += X[3][i]; X[3][i] = RotL_64(X[3][i],R_64x4_7_0); X[3][i] ^= X[0][i]; 
      X[2][i] += X[1][i]; X[1][i] = RotL_64(X[1][i],R_64x4_7_1); X[1][i] ^= X[2][i]; 

      /* InjectKey(r=2) */                                            
      X[0][i] += ks[2]; X[1][i] += ks[3]; X[2][i] += ks[4]; X[3][i] += ks[0]; 
      X[4-1][i] += 2;     /* X[WCNT4-1] += r  */                 


      X[0][i] += X[1][i]; X[1][i] = RotL_64(X[1][i],R_64x4_0_0); X[1][i] ^= X[0][i]; 
      X[2][i] += X[3][i]; X[3][i] = RotL_64(X[3][i],R_64x4_0_1); X[3][i] ^= X[2][i]; 

      X[0][i] += X[3][i]; X[3][i] = RotL_64(X[3][i],R_64x4_1_0); X[3][i] ^= X[0][i]; 
      X[2][i] += X[1][i]; X[1][i] = RotL_64(X[1][i],R_64x4_1_1); X[1][i] ^= X[2][i]; 

      X[0][i] += X[1][i]; X[1][i] = RotL_64(X[1][i],R_64x4_2_0); X[1][i] ^= X[0][i]; 
      X[2][i] += X[3][i]; X[3][i] = RotL_64(X[3][i],R_64x4_2_1); X[3][i] ^= X[2][i]; 

      X[0][i] += X[3][i]; X[3][i] = RotL_64(X[3][i],R_64x4_3_0); X[3][i] ^= X[0][i]; 
      X[2][i] += X[1][i]; X[1][i] = RotL_64(X[1][i],R_64x4_3_1); X[1][i] ^= X[2][i]; 

      /* InjectKey(r=3) */                                            
      X[0][i] += ks[3]; X[1][i] += ks[4]; X[2][i] += ks[0]; X[3][i] += ks[1]; 
      X[4-1][i] += 3;     /* X[WCNT4-1] += r  */                 


      X[0][i] += X[1][i]; X[1][i] = RotL_64(X[1][i],R_64x4_4_0); X[1][i] ^= X[0][i]; 
      X[2][i] += X[3][i]; X[3][i] = RotL_64(X[3][i],R_64x4_4_1); X[3][i] ^= X[2][i]; 

      X[0][i] += X[3][i]; X[3][i] = RotL_64(X[3][i],R_64x4_5_0); X[3][i] ^= X[0][i]; 
      X[2][i] += X[1][i]; X[1][i] = RotL_64(X[1][i],R_64x4_5_1); X[1][i] ^= X[2][i]; 

      X[0][i] += X[1][i]; X[1][i] = RotL_64(X[1][i],R_64x4_6_0); X[1][i] ^= X[0][i]; 
      X[2][i] += X[3][i]; X[3][i] = RotL_64(X[3][i],R_64x4_6_1); X[3][i] ^= X[2][i]; 

      X[0][i] += X[3][i]; X[3][i] = RotL_64(X[3][i],R_64x4_7_0); X[3][i] ^= X[0][i]; 
      X[2][i] += X[1][i]; X[1][i] = RotL_64(X[1][i],R_64x4_7_1); X[1][i] ^= X[2][i]; 

      /* InjectKey(r=1) */                                            
      X[0][i] += ks[4]; X[1][i] += ks[0]; X[2][i] += ks[1]; X[3][i] += ks[2]; 
      X[4-1][i] += 4;     /* X[WCNT4-1] += r  */                 


      X[0][i] += X[1][i]; X[1][i] = RotL_64(X[1][i],R_64x4_0_0); X[1][i] ^= X[0][i]; 
      X[2][i] += X[3][i]; X[3][i] = RotL_64(X[3][i],R_64x4_0_1); X[3][i] ^= X[2][i]; 

      X[0][i] += X[3][i]; X[3][i] = RotL_64(X[3][i],R_64x4_1_0); X[3][i] ^= X[0][i]; 
      X[2][i] += X[1][i]; X[1][i] = RotL_64(X[1][i],R_64x4_1_1); X[1][i] ^= X[2][i]; 

      X[0][i] += X[1][i]; X[1][i] = RotL_64(X[1][i],R_64x4_2_0); X[1][i] ^= X[0][i]; 
      X[2][i] += X[3][i]; X[3][i] = RotL_64(X[3][i],R_64x4_2_1); X[3][i] ^= X[2][i]; 

      X[0][i] += X[3][i]; X[3][i] = RotL_64(X[3][i],R_64x4_3_0); X[3][i] ^= X[0][i]; 
      X[2][i] += X[1][i]; X[1][i] = RotL_64(X[1][i],R_64x4_3_1); X[1][i] ^= X[2][i]; 

      /* InjectKey(r=1) */                                            
      X[0][i] += ks[0]; X[1][i] += ks[1]; X[2][i] += ks[2]; X[3][i] += ks[3]; 
      X[4-1][i] += 5;     /* X[WCNT4-1] += r  */                 


      //Convert integers to doubles uniformly distributed between 0 and 1 inclusive/exlusive
      //u01_closed_open_64_53
      buff[i*4+0] = X[0][i]*R123_0x1p_64;
      buff[i*4+1] = X[1][i]*R123_0x1p_64;
      buff[i*4+2] = X[2][i]*R123_0x1p_64;
      buff[i*4+3] = X[3][i]*R123_0x1p_64;

    }

    return 0; 
}



//old - may not work
int threefry2x64f_multi(unsigned int* CTR4,unsigned int* KEY4, double* buff) {
  
    int i;    
    
    //uint64_t ks[2+1]; //Assuming one set of keys
    __attribute__((aligned(VECTOR_LENGTH_BYTES))) uint64_t ks[2+1]; //Assuming one set of keys
              
    //uint64_t X[2][NUM_VALS_64];
    __attribute__((aligned(VECTOR_LENGTH_BYTES))) uint64_t X[2][NUM_VALS_64];
    
    ks[2] =  SKEIN_KS_PARITY64;

    //Assumes only first counter varies......... big assumption 
    
    for (i=0;i <NUM_VALS_64 ; i++) {
      X[0][i]  = *CTR4 + 2*i; // increment counter
      X[1][i]  = CTR4[1];
    }
    

    //loop over 2 vals (not vector dimension)
    for (i=0;i < 2; i++) {
      ks[i] = KEY4[i];
      ks[2] ^= KEY4[i]; 
    }

    //loop over vector length
    //#pragma simd
    
    for (i=0;i < NUM_VALS_64; i++) {

      X[0][i] += ks[0]; X[1][i] += ks[1];


      X[0][i] += X[1][i]; X[1][i] = RotL_64(X[1][i],R_64x2_0_0); X[1][i] ^= X[0][i]; 
      X[0][i] += X[1][i]; X[1][i] = RotL_64(X[1][i],R_64x2_1_0); X[1][i] ^= X[0][i]; 
      X[0][i] += X[1][i]; X[1][i] = RotL_64(X[1][i],R_64x2_2_0); X[1][i] ^= X[0][i]; 
      X[0][i] += X[1][i]; X[1][i] = RotL_64(X[1][i],R_64x2_3_0); X[1][i] ^= X[0][i]; 


      /* InjectKey(r=1) */                                            
      X[0][i] += ks[1]; X[1][i] += ks[2]; 
      X[1][i] += 1;

      X[0][i] += X[1][i]; X[1][i] = RotL_64(X[1][i],R_64x2_4_0); X[1][i] ^= X[0][i]; 
      X[0][i] += X[1][i]; X[1][i] = RotL_64(X[1][i],R_64x2_5_0); X[1][i] ^= X[0][i]; 
      X[0][i] += X[1][i]; X[1][i] = RotL_64(X[1][i],R_64x2_6_0); X[1][i] ^= X[0][i]; 
      X[0][i] += X[1][i]; X[1][i] = RotL_64(X[1][i],R_64x2_7_0); X[1][i] ^= X[0][i]; 

      /* InjectKey(r=2) */                                            
      X[0][i] += ks[2]; X[1][i] += ks[0]; 
      X[1][i] += 2;


      X[0][i] += X[1][i]; X[1][i] = RotL_64(X[1][i],R_64x2_0_0); X[1][i] ^= X[0][i]; 
      X[0][i] += X[1][i]; X[1][i] = RotL_64(X[1][i],R_64x2_1_0); X[1][i] ^= X[0][i]; 
      X[0][i] += X[1][i]; X[1][i] = RotL_64(X[1][i],R_64x2_2_0); X[1][i] ^= X[0][i]; 
      X[0][i] += X[1][i]; X[1][i] = RotL_64(X[1][i],R_64x2_3_0); X[1][i] ^= X[0][i]; 

      /* InjectKey(r=3) */                                            
      X[0][i] += ks[0]; X[1][i] += ks[1]; 
      X[1][i] += 3;

      X[0][i] += X[1][i]; X[1][i] = RotL_64(X[1][i],R_64x2_4_0); X[1][i] ^= X[0][i]; 
      X[0][i] += X[1][i]; X[1][i] = RotL_64(X[1][i],R_64x2_5_0); X[1][i] ^= X[0][i]; 
      X[0][i] += X[1][i]; X[1][i] = RotL_64(X[1][i],R_64x2_6_0); X[1][i] ^= X[0][i]; 
      X[0][i] += X[1][i]; X[1][i] = RotL_64(X[1][i],R_64x2_7_0); X[1][i] ^= X[0][i]; 


      /* InjectKey(r=4) */                                            
      X[0][i] += ks[1]; X[1][i] += ks[2];  
      X[1][i] += 4;


      X[0][i] += X[1][i]; X[1][i] = RotL_64(X[1][i],R_64x2_0_0); X[1][i] ^= X[0][i]; 
      X[0][i] += X[1][i]; X[1][i] = RotL_64(X[1][i],R_64x2_1_0); X[1][i] ^= X[0][i]; 
      X[0][i] += X[1][i]; X[1][i] = RotL_64(X[1][i],R_64x2_2_0); X[1][i] ^= X[0][i]; 
      X[0][i] += X[1][i]; X[1][i] = RotL_64(X[1][i],R_64x2_3_0); X[1][i] ^= X[0][i]; 

      /* InjectKey(r=4) */                                            
      X[0][i] += ks[2]; X[1][i] += ks[0]; 
      X[1][i] += 5;


      //Convert integers to doubles uniformly distributed between 0 and 1 inclusive/exlusive
      //u01_closed_open_64_53
      buff[i*4+0] = X[0][i]*R123_0x1p_64;
      buff[i*4+1] = X[1][i]*R123_0x1p_64;

    }

    return 0; 
}



// ------------------------------------ threefry4x32f_multi_ctrkey_all -----------------------------------------

// Version with output format double* RAN1, double* RAN2, double* RAN3, double* RAN4
// Version with provided IN_NUM_VALS - prob slower on some systems... to fixed value

// //extern "C"
// // Interface 
// // any = Supply all seeds/keys
// // --- Return 4 vectors of values as output 
// int threefry4x32f_multi_ctrkey_all(unsigned int* CTR1,unsigned int* CTR2,unsigned int* CTR3,unsigned int* CTR4,  
//                                    unsigned int* KEY1,unsigned int* KEY2,unsigned int* KEY3,unsigned int* KEY4, 
//                                    double* RAN1, double* RAN2, double* RAN3, double* RAN4, const int IN_NUM_VALS) {
//   
//     int ivec;  // For iteration over vector dimesion
//     //int ictr;  // For iteration over ctk/keys 1,2,3,4
//           
//     //Use 1D arrays for aligned accesses when supply vector size
//     __attribute__((aligned(VECTOR_LENGTH_BYTES))) uint32_t ks0[IN_NUM_VALS];    
//     __attribute__((aligned(VECTOR_LENGTH_BYTES))) uint32_t ks1[IN_NUM_VALS];
//     __attribute__((aligned(VECTOR_LENGTH_BYTES))) uint32_t ks2[IN_NUM_VALS];
//     __attribute__((aligned(VECTOR_LENGTH_BYTES))) uint32_t ks3[IN_NUM_VALS];
//     __attribute__((aligned(VECTOR_LENGTH_BYTES))) uint32_t ks4[IN_NUM_VALS];
// 
//               
//     //Use 1D arrays for aligned accesses when supply vector size
//     __attribute__((aligned(VECTOR_LENGTH_BYTES))) uint32_t X0[IN_NUM_VALS];
//     __attribute__((aligned(VECTOR_LENGTH_BYTES))) uint32_t X1[IN_NUM_VALS];
//     __attribute__((aligned(VECTOR_LENGTH_BYTES))) uint32_t X2[IN_NUM_VALS];
//     __attribute__((aligned(VECTOR_LENGTH_BYTES))) uint32_t X3[IN_NUM_VALS];
// 
// 
//     //__attribute__((aligned(VECTOR_LENGTH_BYTES))) double buff[IN_NUM_VALS*4];
//     
// 
//     for (ivec=0;ivec <IN_NUM_VALS ; ivec++) {
//       X0[ivec]  = CTR1[ivec];
//       X1[ivec]  = CTR2[ivec];
//       X2[ivec]  = CTR3[ivec];
//       X3[ivec]  = CTR4[ivec];       
//       ks4[ivec] =  SKEIN_KS_PARITY32;                          
//     }
// 
// 
//     for (ivec=0;ivec <IN_NUM_VALS ; ivec++) {
//         ks0[ivec] = KEY1[ivec];
//         ks1[ivec] = KEY2[ivec];
//         ks2[ivec] = KEY3[ivec];
//         ks3[ivec] = KEY4[ivec];
//         
//         //could interleave or here
//         ks4[ivec] ^= KEY1[ivec]; 
//         ks4[ivec] ^= KEY2[ivec]; 
//         ks4[ivec] ^= KEY3[ivec]; 
//         ks4[ivec] ^= KEY4[ivec]; 
//     }
// 
// 
// 
//     //Loop over vector length (or multiple thereof)
//     #pragma omp simd aligned(X0,X1,X2,X3,ks0,ks1,ks2,ks3,ks4)
//     for (ivec=0;ivec < IN_NUM_VALS; ivec++) {
// 
//       X0[ivec] += ks0[ivec]; X1[ivec] += ks1[ivec]; X2[ivec] += ks2[ivec]; X3[ivec] += ks3[ivec];
// 
// 
//       X0[ivec] += X1[ivec]; X1[ivec] = RotL_32(X1[ivec],R_32x4_0_0); X1[ivec] ^= X0[ivec]; 
//       X2[ivec] += X3[ivec]; X3[ivec] = RotL_32(X3[ivec],R_32x4_0_1); X3[ivec] ^= X2[ivec]; 
// 
//       X0[ivec] += X3[ivec]; X3[ivec] = RotL_32(X3[ivec],R_32x4_1_0); X3[ivec] ^= X0[ivec]; 
//       X2[ivec] += X1[ivec]; X1[ivec] = RotL_32(X1[ivec],R_32x4_1_1); X1[ivec] ^= X2[ivec]; 
// 
//       X0[ivec] += X1[ivec]; X1[ivec] = RotL_32(X1[ivec],R_32x4_2_0); X1[ivec] ^= X0[ivec]; 
//       X2[ivec] += X3[ivec]; X3[ivec] = RotL_32(X3[ivec],R_32x4_2_1); X3[ivec] ^= X2[ivec]; 
// 
//       X0[ivec] += X3[ivec]; X3[ivec] = RotL_32(X3[ivec],R_32x4_3_0); X3[ivec] ^= X0[ivec]; 
//       X2[ivec] += X1[ivec]; X1[ivec] = RotL_32(X1[ivec],R_32x4_3_1); X1[ivec] ^= X2[ivec]; 
// 
//       /* InjectKey(r=1) */                                            
//       X0[ivec] += ks1[ivec]; X1[ivec] += ks2[ivec]; X2[ivec] += ks3[ivec]; X3[ivec] += ks4[ivec]; 
//       X3[ivec] += 1;     /* X[WCNT4-1] += r  */                 
// 
// 
//       X0[ivec] += X1[ivec]; X1[ivec] = RotL_32(X1[ivec],R_32x4_4_0); X1[ivec] ^= X0[ivec]; 
//       X2[ivec] += X3[ivec]; X3[ivec] = RotL_32(X3[ivec],R_32x4_4_1); X3[ivec] ^= X2[ivec]; 
// 
//       X0[ivec] += X3[ivec]; X3[ivec] = RotL_32(X3[ivec],R_32x4_5_0); X3[ivec] ^= X0[ivec]; 
//       X2[ivec] += X1[ivec]; X1[ivec] = RotL_32(X1[ivec],R_32x4_5_1); X1[ivec] ^= X2[ivec]; 
// 
//       X0[ivec] += X1[ivec]; X1[ivec] = RotL_32(X1[ivec],R_32x4_6_0); X1[ivec] ^= X0[ivec]; 
//       X2[ivec] += X3[ivec]; X3[ivec] = RotL_32(X3[ivec],R_32x4_6_1); X3[ivec] ^= X2[ivec]; 
// 
//       X0[ivec] += X3[ivec]; X3[ivec] = RotL_32(X3[ivec],R_32x4_7_0); X3[ivec] ^= X0[ivec]; 
//       X2[ivec] += X1[ivec]; X1[ivec] = RotL_32(X1[ivec],R_32x4_7_1); X1[ivec] ^= X2[ivec]; 
// 
//       /* InjectKey(r=2) */                                            
//       X0[ivec] += ks2[ivec]; X1[ivec] += ks3[ivec]; X2[ivec] += ks4[ivec]; X3[ivec] += ks0[ivec]; 
//       X3[ivec] += 2;     /* X[WCNT4-1] += r  */                 
// 
// 
//       X0[ivec] += X1[ivec]; X1[ivec] = RotL_32(X1[ivec],R_32x4_0_0); X1[ivec] ^= X0[ivec]; 
//       X2[ivec] += X3[ivec]; X3[ivec] = RotL_32(X3[ivec],R_32x4_0_1); X3[ivec] ^= X2[ivec]; 
// 
//       X0[ivec] += X3[ivec]; X3[ivec] = RotL_32(X3[ivec],R_32x4_1_0); X3[ivec] ^= X0[ivec]; 
//       X2[ivec] += X1[ivec]; X1[ivec] = RotL_32(X1[ivec],R_32x4_1_1); X1[ivec] ^= X2[ivec]; 
// 
//       X0[ivec] += X1[ivec]; X1[ivec] = RotL_32(X1[ivec],R_32x4_2_0); X1[ivec] ^= X0[ivec]; 
//       X2[ivec] += X3[ivec]; X3[ivec] = RotL_32(X3[ivec],R_32x4_2_1); X3[ivec] ^= X2[ivec]; 
// 
//       X0[ivec] += X3[ivec]; X3[ivec] = RotL_32(X3[ivec],R_32x4_3_0); X3[ivec] ^= X0[ivec]; 
//       X2[ivec] += X1[ivec]; X1[ivec] = RotL_32(X1[ivec],R_32x4_3_1); X1[ivec] ^= X2[ivec]; 
// 
//       /* InjectKey(r=3) */                                            
//       X0[ivec] += ks3[ivec]; X1[ivec] += ks4[ivec]; X2[ivec] += ks0[ivec]; X3[ivec] += ks1[ivec]; 
//       X3[ivec] += 3;     /* X[WCNT4-1] += r  */                 
// 
// 
//       X0[ivec] += X1[ivec]; X1[ivec] = RotL_32(X1[ivec],R_32x4_4_0); X1[ivec] ^= X0[ivec]; 
//       X2[ivec] += X3[ivec]; X3[ivec] = RotL_32(X3[ivec],R_32x4_4_1); X3[ivec] ^= X2[ivec]; 
// 
//       X0[ivec] += X3[ivec]; X3[ivec] = RotL_32(X3[ivec],R_32x4_5_0); X3[ivec] ^= X0[ivec]; 
//       X2[ivec] += X1[ivec]; X1[ivec] = RotL_32(X1[ivec],R_32x4_5_1); X1[ivec] ^= X2[ivec]; 
// 
//       X0[ivec] += X1[ivec]; X1[ivec] = RotL_32(X1[ivec],R_32x4_6_0); X1[ivec] ^= X0[ivec]; 
//       X2[ivec] += X3[ivec]; X3[ivec] = RotL_32(X3[ivec],R_32x4_6_1); X3[ivec] ^= X2[ivec]; 
// 
//       X0[ivec] += X3[ivec]; X3[ivec] = RotL_32(X3[ivec],R_32x4_7_0); X3[ivec] ^= X0[ivec]; 
//       X2[ivec] += X1[ivec]; X1[ivec] = RotL_32(X1[ivec],R_32x4_7_1); X1[ivec] ^= X2[ivec]; 
// 
//       /* InjectKey(r=1) */                                            
//       X0[ivec] += ks4[ivec]; X1[ivec] += ks0[ivec]; X2[ivec] += ks1[ivec]; X3[ivec] += ks2[ivec]; 
//       X3[ivec] += 4;     /* X[WCNT4-1] += r  */                 
// 
// 
//       X0[ivec] += X1[ivec]; X1[ivec] = RotL_32(X1[ivec],R_32x4_0_0); X1[ivec] ^= X0[ivec]; 
//       X2[ivec] += X3[ivec]; X3[ivec] = RotL_32(X3[ivec],R_32x4_0_1); X3[ivec] ^= X2[ivec]; 
// 
//       X0[ivec] += X3[ivec]; X3[ivec] = RotL_32(X3[ivec],R_32x4_1_0); X3[ivec] ^= X0[ivec]; 
//       X2[ivec] += X1[ivec]; X1[ivec] = RotL_32(X1[ivec],R_32x4_1_1); X1[ivec] ^= X2[ivec]; 
// 
//       X0[ivec] += X1[ivec]; X1[ivec] = RotL_32(X1[ivec],R_32x4_2_0); X1[ivec] ^= X0[ivec]; 
//       X2[ivec] += X3[ivec]; X3[ivec] = RotL_32(X3[ivec],R_32x4_2_1); X3[ivec] ^= X2[ivec]; 
// 
//       X0[ivec] += X3[ivec]; X3[ivec] = RotL_32(X3[ivec],R_32x4_3_0); X3[ivec] ^= X0[ivec]; 
//       X2[ivec] += X1[ivec]; X1[ivec] = RotL_32(X1[ivec],R_32x4_3_1); X1[ivec] ^= X2[ivec]; 
// 
//       /* InjectKey(r=1) */                                            
//       X0[ivec] += ks0[ivec]; X1[ivec] += ks1[ivec]; X2[ivec] += ks2[ivec]; X3[ivec] += ks3[ivec]; 
//       X3[ivec] += 5;     /* X[WCNT4-1] += r  */                 
// 
// 
//       //Convert integers to doubles uniformly distributed between 0 and 1 inclusive/exlusive
//       //u01_closed_open_32_53
//       RAN1[ivec] = X0[ivec]*R123_0x1p_32;
//       RAN2[ivec] = X1[ivec]*R123_0x1p_32;
//       RAN3[ivec] = X2[ivec]*R123_0x1p_32;
//       RAN4[ivec] = X3[ivec]*R123_0x1p_32;
//       
// 
//     }
//     
//     return 0; 
// }





//-----------------------------------------------------------------------------------------------
// Version with output format double* buff - same order as original - means scattered writes...
// Version with provided IN_NUM_VALS - prob slower on some systems... to fixed value

//extern "C"
// Interface 
// any = Supply all seeds/keys
// --- Return 4 vectors of values as output 
int threefry4x32f_multi_ctrkey_all(unsigned int* CTR1,unsigned int* CTR2,unsigned int* CTR3,unsigned int* CTR4,  
                                   unsigned int* KEY1,unsigned int* KEY2,unsigned int* KEY3,unsigned int* KEY4, 
                                   double* buff, const int IN_NUM_VALS) {
  
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
//       RAN1[ivec] = X0[ivec]*R123_0x1p_32;
//       RAN2[ivec] = X1[ivec]*R123_0x1p_32;
//       RAN3[ivec] = X2[ivec]*R123_0x1p_32;
//       RAN4[ivec] = X3[ivec]*R123_0x1p_32;
      
      buff[ivec*4+0] = X0[ivec]*R123_0x1p_32;
      buff[ivec*4+1] = X1[ivec]*R123_0x1p_32;
      buff[ivec*4+2] = X2[ivec]*R123_0x1p_32;
      buff[ivec*4+3] = X3[ivec]*R123_0x1p_32;

    }
    
    return 0; 
}

