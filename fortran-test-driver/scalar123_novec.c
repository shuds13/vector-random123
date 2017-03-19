//--------------------------------------------------------------------------------------------------
// S. Hudson: Standard C implementation of Random123 functions.
// Based on Random123 functions - see copyright notice.
// scalar123_novec is essentially a duplicate file to scalar123.c but
// with renamed routines - which are then compiled with novec for comparison
// with vector compiled code in one program
// This is just for test harness
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
// S. Hudson: Standalone version - threefry
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

//=========================================================================================================
//SINGLE CALL VERSIONS
//=========================================================================================================

//extern "C"
int threefry4x32f_novec(unsigned int* CTR4,unsigned int* KEY4, double* buff) {
  
    
    //use ks[4] as fits in rotation of injection keys
    //uint32_t ks4
    //ks4 =  SKEIN_KS_PARITY32;
    int i;
    uint32_t ks[4+1]; 
    uint32_t X[4]; 
       
    ks[4] =  SKEIN_KS_PARITY32; 
       
    for (i=0;i < 4; i++) {
      ks[i] = KEY4[i];
      X[i]  = CTR4[i];
      ks[4] ^= KEY4[i]; 
    }
    
    //printf("Hello from vector 123 - scalar version\n");
    
    //See how instructions are interleaved with pairs in each
    //So already some independence - atleast 2 way
    //Maybe why multi calls does not make that much difference without vect (ie from ALUs).
                            
   //InjectKey bit can be vectorised when contingous but not when wraps around
   //Might want to force alignment also.
    X[0] += ks[0]; X[1] += ks[1]; X[2] += ks[2]; X[3] += ks[3];
    
    
    X[0] += X[1]; X[1] = RotL_32(X[1],R_32x4_0_0); X[1] ^= X[0]; 
    X[2] += X[3]; X[3] = RotL_32(X[3],R_32x4_0_1); X[3] ^= X[2]; 

    X[0] += X[3]; X[3] = RotL_32(X[3],R_32x4_1_0); X[3] ^= X[0]; 
    X[2] += X[1]; X[1] = RotL_32(X[1],R_32x4_1_1); X[1] ^= X[2]; 

    X[0] += X[1]; X[1] = RotL_32(X[1],R_32x4_2_0); X[1] ^= X[0]; 
    X[2] += X[3]; X[3] = RotL_32(X[3],R_32x4_2_1); X[3] ^= X[2]; 

    X[0] += X[3]; X[3] = RotL_32(X[3],R_32x4_3_0); X[3] ^= X[0]; 
    X[2] += X[1]; X[1] = RotL_32(X[1],R_32x4_3_1); X[1] ^= X[2]; 

    /* InjectKey(r=1) */                                            
    X[0] += ks[1]; X[1] += ks[2]; X[2] += ks[3]; X[3] += ks[4]; 
    X[4-1] += 1;     /* X[WCNT4-1] += r  */                 

                                                                    
    X[0] += X[1]; X[1] = RotL_32(X[1],R_32x4_4_0); X[1] ^= X[0]; 
    X[2] += X[3]; X[3] = RotL_32(X[3],R_32x4_4_1); X[3] ^= X[2]; 

    X[0] += X[3]; X[3] = RotL_32(X[3],R_32x4_5_0); X[3] ^= X[0]; 
    X[2] += X[1]; X[1] = RotL_32(X[1],R_32x4_5_1); X[1] ^= X[2]; 

    X[0] += X[1]; X[1] = RotL_32(X[1],R_32x4_6_0); X[1] ^= X[0]; 
    X[2] += X[3]; X[3] = RotL_32(X[3],R_32x4_6_1); X[3] ^= X[2]; 

    X[0] += X[3]; X[3] = RotL_32(X[3],R_32x4_7_0); X[3] ^= X[0]; 
    X[2] += X[1]; X[1] = RotL_32(X[1],R_32x4_7_1); X[1] ^= X[2]; 

    /* InjectKey(r=2) */                                            
    X[0] += ks[2]; X[1] += ks[3]; X[2] += ks[4]; X[3] += ks[0]; 
    X[4-1] += 2;     /* X[WCNT4-1] += r  */                 


    X[0] += X[1]; X[1] = RotL_32(X[1],R_32x4_0_0); X[1] ^= X[0]; 
    X[2] += X[3]; X[3] = RotL_32(X[3],R_32x4_0_1); X[3] ^= X[2]; 

    X[0] += X[3]; X[3] = RotL_32(X[3],R_32x4_1_0); X[3] ^= X[0]; 
    X[2] += X[1]; X[1] = RotL_32(X[1],R_32x4_1_1); X[1] ^= X[2]; 

    X[0] += X[1]; X[1] = RotL_32(X[1],R_32x4_2_0); X[1] ^= X[0]; 
    X[2] += X[3]; X[3] = RotL_32(X[3],R_32x4_2_1); X[3] ^= X[2]; 

    X[0] += X[3]; X[3] = RotL_32(X[3],R_32x4_3_0); X[3] ^= X[0]; 
    X[2] += X[1]; X[1] = RotL_32(X[1],R_32x4_3_1); X[1] ^= X[2]; 

    /* InjectKey(r=3) */                                            
    X[0] += ks[3]; X[1] += ks[4]; X[2] += ks[0]; X[3] += ks[1]; 
    X[4-1] += 3;     /* X[WCNT4-1] += r  */                 


    X[0] += X[1]; X[1] = RotL_32(X[1],R_32x4_4_0); X[1] ^= X[0]; 
    X[2] += X[3]; X[3] = RotL_32(X[3],R_32x4_4_1); X[3] ^= X[2]; 

    X[0] += X[3]; X[3] = RotL_32(X[3],R_32x4_5_0); X[3] ^= X[0]; 
    X[2] += X[1]; X[1] = RotL_32(X[1],R_32x4_5_1); X[1] ^= X[2]; 

    X[0] += X[1]; X[1] = RotL_32(X[1],R_32x4_6_0); X[1] ^= X[0]; 
    X[2] += X[3]; X[3] = RotL_32(X[3],R_32x4_6_1); X[3] ^= X[2]; 

    X[0] += X[3]; X[3] = RotL_32(X[3],R_32x4_7_0); X[3] ^= X[0]; 
    X[2] += X[1]; X[1] = RotL_32(X[1],R_32x4_7_1); X[1] ^= X[2]; 

    /* InjectKey(r=1) */                                            
    X[0] += ks[4]; X[1] += ks[0]; X[2] += ks[1]; X[3] += ks[2]; 
    X[4-1] += 4;     /* X[WCNT4-1] += r  */                 


    X[0] += X[1]; X[1] = RotL_32(X[1],R_32x4_0_0); X[1] ^= X[0]; 
    X[2] += X[3]; X[3] = RotL_32(X[3],R_32x4_0_1); X[3] ^= X[2]; 

    X[0] += X[3]; X[3] = RotL_32(X[3],R_32x4_1_0); X[3] ^= X[0]; 
    X[2] += X[1]; X[1] = RotL_32(X[1],R_32x4_1_1); X[1] ^= X[2]; 

    X[0] += X[1]; X[1] = RotL_32(X[1],R_32x4_2_0); X[1] ^= X[0]; 
    X[2] += X[3]; X[3] = RotL_32(X[3],R_32x4_2_1); X[3] ^= X[2]; 

    X[0] += X[3]; X[3] = RotL_32(X[3],R_32x4_3_0); X[3] ^= X[0]; 
    X[2] += X[1]; X[1] = RotL_32(X[1],R_32x4_3_1); X[1] ^= X[2]; 

    /* InjectKey(r=1) */                                            
    X[0] += ks[0]; X[1] += ks[1]; X[2] += ks[2]; X[3] += ks[3]; 
    X[4-1] += 5;     /* X[WCNT4-1] += r  */                 
    
    //Convert integers to doubles uniformly distributed between 0 and 1 inclusive/exlusive 
    //u01_closed_open_32_53
    buff[0] = X[0]*R123_0x1p_32;
    buff[1] = X[1]*R123_0x1p_32;
    buff[2] = X[2]*R123_0x1p_32;
    buff[3] = X[3]*R123_0x1p_32;
    
    
  // Again Fortran will expect a pointer to the value
  //return randval; 
  return 0; 
}





//extern "C"
int threefry4x64f_novec(unsigned int* CTR4,unsigned int* KEY4, double* buff) {
  
    
    int i;
    uint64_t ks[4+1]; 
    uint64_t X[4]; 
       
    ks[4] =  SKEIN_KS_PARITY64; 
       
    for (i=0;i < 4; i++) {
      ks[i] = KEY4[i];
      X[i]  = CTR4[i];
      ks[4] ^= KEY4[i]; 
    }
    
    //See how instructions are interleaved with pairs in each
    //So already some independence - atleast 2 way
    //Maybe why multi calls does not make that much difference without vect (ie from ALUs).
                            
   //InjectKey bit can be vectorised when contingous but not when wraps around
   //Might want to force alignment also.
    X[0] += ks[0]; X[1] += ks[1]; X[2] += ks[2]; X[3] += ks[3];
    
    
    X[0] += X[1]; X[1] = RotL_64(X[1],R_64x4_0_0); X[1] ^= X[0]; 
    X[2] += X[3]; X[3] = RotL_64(X[3],R_64x4_0_1); X[3] ^= X[2]; 

    X[0] += X[3]; X[3] = RotL_64(X[3],R_64x4_1_0); X[3] ^= X[0]; 
    X[2] += X[1]; X[1] = RotL_64(X[1],R_64x4_1_1); X[1] ^= X[2]; 

    X[0] += X[1]; X[1] = RotL_64(X[1],R_64x4_2_0); X[1] ^= X[0]; 
    X[2] += X[3]; X[3] = RotL_64(X[3],R_64x4_2_1); X[3] ^= X[2]; 

    X[0] += X[3]; X[3] = RotL_64(X[3],R_64x4_3_0); X[3] ^= X[0]; 
    X[2] += X[1]; X[1] = RotL_64(X[1],R_64x4_3_1); X[1] ^= X[2]; 

    /* InjectKey(r=1) */                                            
    X[0] += ks[1]; X[1] += ks[2]; X[2] += ks[3]; X[3] += ks[4]; 
    X[4-1] += 1;     /* X[WCNT4-1] += r  */                 

                                                                    
    X[0] += X[1]; X[1] = RotL_64(X[1],R_64x4_4_0); X[1] ^= X[0]; 
    X[2] += X[3]; X[3] = RotL_64(X[3],R_64x4_4_1); X[3] ^= X[2]; 

    X[0] += X[3]; X[3] = RotL_64(X[3],R_64x4_5_0); X[3] ^= X[0]; 
    X[2] += X[1]; X[1] = RotL_64(X[1],R_64x4_5_1); X[1] ^= X[2]; 

    X[0] += X[1]; X[1] = RotL_64(X[1],R_64x4_6_0); X[1] ^= X[0]; 
    X[2] += X[3]; X[3] = RotL_64(X[3],R_64x4_6_1); X[3] ^= X[2]; 

    X[0] += X[3]; X[3] = RotL_64(X[3],R_64x4_7_0); X[3] ^= X[0]; 
    X[2] += X[1]; X[1] = RotL_64(X[1],R_64x4_7_1); X[1] ^= X[2]; 

    /* InjectKey(r=2) */                                            
    X[0] += ks[2]; X[1] += ks[3]; X[2] += ks[4]; X[3] += ks[0]; 
    X[4-1] += 2;     /* X[WCNT4-1] += r  */                 


    X[0] += X[1]; X[1] = RotL_64(X[1],R_64x4_0_0); X[1] ^= X[0]; 
    X[2] += X[3]; X[3] = RotL_64(X[3],R_64x4_0_1); X[3] ^= X[2]; 

    X[0] += X[3]; X[3] = RotL_64(X[3],R_64x4_1_0); X[3] ^= X[0]; 
    X[2] += X[1]; X[1] = RotL_64(X[1],R_64x4_1_1); X[1] ^= X[2]; 

    X[0] += X[1]; X[1] = RotL_64(X[1],R_64x4_2_0); X[1] ^= X[0]; 
    X[2] += X[3]; X[3] = RotL_64(X[3],R_64x4_2_1); X[3] ^= X[2]; 

    X[0] += X[3]; X[3] = RotL_64(X[3],R_64x4_3_0); X[3] ^= X[0]; 
    X[2] += X[1]; X[1] = RotL_64(X[1],R_64x4_3_1); X[1] ^= X[2]; 

    /* InjectKey(r=3) */                                            
    X[0] += ks[3]; X[1] += ks[4]; X[2] += ks[0]; X[3] += ks[1]; 
    X[4-1] += 3;     /* X[WCNT4-1] += r  */                 


    X[0] += X[1]; X[1] = RotL_64(X[1],R_64x4_4_0); X[1] ^= X[0]; 
    X[2] += X[3]; X[3] = RotL_64(X[3],R_64x4_4_1); X[3] ^= X[2]; 

    X[0] += X[3]; X[3] = RotL_64(X[3],R_64x4_5_0); X[3] ^= X[0]; 
    X[2] += X[1]; X[1] = RotL_64(X[1],R_64x4_5_1); X[1] ^= X[2]; 

    X[0] += X[1]; X[1] = RotL_64(X[1],R_64x4_6_0); X[1] ^= X[0]; 
    X[2] += X[3]; X[3] = RotL_64(X[3],R_64x4_6_1); X[3] ^= X[2]; 

    X[0] += X[3]; X[3] = RotL_64(X[3],R_64x4_7_0); X[3] ^= X[0]; 
    X[2] += X[1]; X[1] = RotL_64(X[1],R_64x4_7_1); X[1] ^= X[2]; 

    /* InjectKey(r=1) */                                            
    X[0] += ks[4]; X[1] += ks[0]; X[2] += ks[1]; X[3] += ks[2]; 
    X[4-1] += 4;     /* X[WCNT4-1] += r  */                 


    X[0] += X[1]; X[1] = RotL_64(X[1],R_64x4_0_0); X[1] ^= X[0]; 
    X[2] += X[3]; X[3] = RotL_64(X[3],R_64x4_0_1); X[3] ^= X[2]; 

    X[0] += X[3]; X[3] = RotL_64(X[3],R_64x4_1_0); X[3] ^= X[0]; 
    X[2] += X[1]; X[1] = RotL_64(X[1],R_64x4_1_1); X[1] ^= X[2]; 

    X[0] += X[1]; X[1] = RotL_64(X[1],R_64x4_2_0); X[1] ^= X[0]; 
    X[2] += X[3]; X[3] = RotL_64(X[3],R_64x4_2_1); X[3] ^= X[2]; 

    X[0] += X[3]; X[3] = RotL_64(X[3],R_64x4_3_0); X[3] ^= X[0]; 
    X[2] += X[1]; X[1] = RotL_64(X[1],R_64x4_3_1); X[1] ^= X[2]; 

    /* InjectKey(r=1) */                                            
    X[0] += ks[0]; X[1] += ks[1]; X[2] += ks[2]; X[3] += ks[3]; 
    X[4-1] += 5;     /* X[WCNT4-1] += r  */                 
    
    //Convert integers to doubles uniformly distributed between 0 and 1 inclusive/exlusive 
    //u01_closed_open_64_53
    buff[0] = X[0]*R123_0x1p_64;
    buff[1] = X[1]*R123_0x1p_64;
    buff[2] = X[2]*R123_0x1p_64;
    buff[3] = X[3]*R123_0x1p_64;
    
    
  // Again Fortran will expect a pointer to the value
  //return randval; 
  return 0; 
}



//extern "C"
int threefry2x64f_novec(unsigned int* CTR4,unsigned int* KEY4, double* buff) {
  
    
    int i;
    uint64_t ks[2+1]; 
    uint64_t X[2]; 
       
    ks[2] =  SKEIN_KS_PARITY64; 
       
    for (i=0;i < 2; i++) {
      ks[i] = KEY4[i];
      X[i]  = CTR4[i];
      ks[2] ^= KEY4[i]; 
    }
    
    //See how instructions are interleaved with pairs in each
    //So already some independence - atleast 2 way
    //Maybe why multi calls does not make that much difference without vect (ie from ALUs).
                            
   //InjectKey bit can be vectorised when contingous but not when wraps around
   //Might want to force alignment also.
    X[0] += ks[0]; X[1] += ks[1];
    
    
    X[0] += X[1]; X[1] = RotL_64(X[1],R_64x2_0_0); X[1] ^= X[0]; 
    X[0] += X[1]; X[1] = RotL_64(X[1],R_64x2_1_0); X[1] ^= X[0]; 
    X[0] += X[1]; X[1] = RotL_64(X[1],R_64x2_2_0); X[1] ^= X[0]; 
    X[0] += X[1]; X[1] = RotL_64(X[1],R_64x2_3_0); X[1] ^= X[0]; 

    /* InjectKey(r=1) */                                            
    X[0] += ks[1]; X[1] += ks[2]; 
    X[1] += 1;

                                                                    
    X[0] += X[1]; X[1] = RotL_64(X[1],R_64x2_4_0); X[1] ^= X[0]; 
    X[0] += X[1]; X[1] = RotL_64(X[1],R_64x2_5_0); X[1] ^= X[0]; 
    X[0] += X[1]; X[1] = RotL_64(X[1],R_64x2_6_0); X[1] ^= X[0]; 
    X[0] += X[1]; X[1] = RotL_64(X[1],R_64x2_7_0); X[1] ^= X[0]; 


    /* InjectKey(r=2) */                                            
    X[0] += ks[2]; X[1] += ks[0];
    X[1] += 2;

    X[0] += X[1]; X[1] = RotL_64(X[1],R_64x2_0_0); X[1] ^= X[0]; 
    X[0] += X[1]; X[1] = RotL_64(X[1],R_64x2_1_0); X[1] ^= X[0]; 
    X[0] += X[1]; X[1] = RotL_64(X[1],R_64x2_2_0); X[1] ^= X[0]; 
    X[0] += X[1]; X[1] = RotL_64(X[1],R_64x2_3_0); X[1] ^= X[0]; 

    /* InjectKey(r=3) */  
    X[0] += ks[0]; X[1] += ks[1];
    X[1] += 3;

    X[0] += X[1]; X[1] = RotL_64(X[1],R_64x2_4_0); X[1] ^= X[0]; 
    X[0] += X[1]; X[1] = RotL_64(X[1],R_64x2_5_0); X[1] ^= X[0]; 
    X[0] += X[1]; X[1] = RotL_64(X[1],R_64x2_6_0); X[1] ^= X[0]; 
    X[0] += X[1]; X[1] = RotL_64(X[1],R_64x2_7_0); X[1] ^= X[0]; 

    /* InjectKey(r=4) */  
    X[0] += ks[1]; X[1] += ks[2];
    X[1] += 4;

    X[0] += X[1]; X[1] = RotL_64(X[1],R_64x2_0_0); X[1] ^= X[0]; 
    X[0] += X[1]; X[1] = RotL_64(X[1],R_64x2_1_0); X[1] ^= X[0]; 
    X[0] += X[1]; X[1] = RotL_64(X[1],R_64x2_2_0); X[1] ^= X[0]; 
    X[0] += X[1]; X[1] = RotL_64(X[1],R_64x2_3_0); X[1] ^= X[0]; 

    /* InjectKey(r=4) */  
    X[0] += ks[2]; X[1] += ks[0];
    X[1] += 5;


    //Convert integers to doubles uniformly distributed between 0 and 1 inclusive/exlusive 
    //u01_closed_open_64_53
    buff[0] = X[0]*R123_0x1p_64;
    buff[1] = X[1]*R123_0x1p_64;
    
    
  // Again Fortran will expect a pointer to the value
  //return randval; 
  return 0; 
}
