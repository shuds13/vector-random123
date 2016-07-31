//--------------------------------------------------------------------------------------------------
// S. Hudson: Vector123 loop only test - Intrinsics plus portable
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

#include <time.h>
#include <sys/time.h>
#include <stdio.h>
#include <stdlib.h>
#include "immintrin.h"

// Allows use of short forms like uint32_t
#include <stdint.h>


//Random123 - To pick up macros from R123 - uncomment
//#include <threefry.h>
//#include <philox.h>

//#include <u01.h>
//OR
//#include <u01fixedpt.h>

// -----------------------------------------------------------------------------
// Is supported vector instruction set detected
// -----------------------------------------------------------------------------

//Comment out this line to disable intrinsic versions
#define ENABLE_INSTRINSIC

#ifdef ENABLE_INSTRINSIC

#if    defined(__AVX512F__)
#define INTRINSICS
#define VLEN 512
#define VEC_ROL
#define VECTOR_LENGTH_BYTES 64
#define VSTRING "AVX512"

#elif  defined(__AVX2__)
#define INTRINSICS
#define VLEN 256
#define VECTOR_LENGTH_BYTES 32
#define VSTRING "AVX2"

#elif  defined(__AVX__)
#define INTRINSICS
#define VLEN 128   #Integer vlen only avx-128
#define VECTOR_LENGTH_BYTES 32
#define VSTRING "AVX"

#else
#define VECTOR_LENGTH_BYTES 64 //Default high
#define VSTRING "None"
#endif

#else // ENABLE_INSTRINSIC not set
#define VECTOR_LENGTH_BYTES 64 //Default high
#define VSTRING "None"
#define VLEN Unknown
#endif

// -----------------------------------------------------------------------------

#define NUM_VALS_32 16000000
#define NUM_VALS_64 8000000

#define SKEIN_MK_64(hi32,lo32)  ((lo32) + (((uint64_t) (hi32)) << 32))
#define SKEIN_KS_PARITY64         SKEIN_MK_64(0x1BD11BDA,0xA9FC1A22)
#define SKEIN_KS_PARITY32         0x1BD11BDA

#define R123_0x1p_32                (1./4294967296.)   // for 32_53 CO,OC,OO int4 to d.p


// Define halfround
#ifdef VEC_ROL

#define halfround(X,Y,ROL_VAL) \
      X = _mmVLEN_add_epi32(X, Y); \
      Y = _mmVLEN_rol_epi32(Y, ROL_VAL); \
      Y = _mmVLEN_xor_epi32(Y, X);

#else 

#define halfround(X,Y,ROL_VAL) \
      X = _mmVLEN_add_epi32(X, Y); \
      Y = _mmVLEN_myrol_epi32(Y, ROL_VAL);  \
      Y = _mmVLEN_xor_epi32(Y, X);

//Could put in func directly
#define _mmVLEN_myrol_epi32(Y,ROL_VAL)  \
	Y = _mmVLEN_slli_epi32(Y, ROL_VAL);  \
	__mVLENi epi32L = _mmVLEN_srli_epi32(Y, 32 - ROL_VAL);  \
  Y = _mmVLEN_or_epi32(Y, epi32L)
#endif


static inline uint32_t RotL_32(uint32_t x, unsigned int N)
{
    return (x << (N & 31)) | (x >> ((32-N) & 31));
}

static inline uint64_t RotL_64(uint64_t x, unsigned int N)
{
    return (x << (N & 63)) | (x >> ((64-N) & 63));
}

//*check - not sure why epi32A needs to be a constant
//template <int nCount>
//__mVLENi _mmVLEN_myrol_epi32(__mVLENi const & epi32A)

// static inline __mVLENi _mmVLEN_myrol_epi32(__mVLENi epi32A,unsigned int nCount)
// {
// 	__mVLENi epi32H = _mmVLEN_slli_epi32(epi32A, nCount);
// 	__mVLENi epi32L = _mmVLEN_srli_epi32(epi32A, 32 - nCount);
// 	return _mmVLEN_or_epi32(epi32H, epi32L);
// 	//return _mmVLEN_or_siVLEN(epi32AH, epi32L);	//what would be the difference?
// }


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


//Utility call
void test_vec_sizes() {
  //Default values vector sizes NUM_VALS_32 and NUM_VALS_64
  //In some functions these can be overwritten with provided value
  printf("--------------------------------------------------------------------- \n");
  printf("Default vector lengths from vector123.c: Note some routines overwrite \n");
  printf("VECTOR_LENGTH_BYTES %d\n",VECTOR_LENGTH_BYTES);
  printf("NUM_VALS_32 %d\n",NUM_VALS_32);
  printf("NUM_VALS_64 %d\n",NUM_VALS_64); 
#ifdef INTRINSICS  
  printf("Intrinsics supported vector instruction set detected %s\n",VSTRING); 
  printf("VECTOR_LENGTH_BITS %d\n",VLEN);   
#else
  printf("Intrinsics supported vector instruction set NOT detected OR intrinsics disabled\n");
  printf("Running portable code\n");    
#endif
}


#ifdef INTRINSICS
int main (void)
{

    int ivec;    
    int ictr;  // For iteration over ctk/keys 1,2,3,4
    
    //clock - cpu timing
    clock_t start, diff;
    int msec;

    //gettimeofday - wall clock timing
    struct timeval  tv1, tv2;
    
    //*check using omp to tell compiler can assume aligned on the loop
    // so really want to use omp to force alignement on declaration
        
    //uint32_t ks[4+1]; //Assuming one set of keys
     //   __attribute__((aligned(VECTOR_LENGTH_BYTES))) uint32_t ks[4+1];
     __attribute__((aligned(VECTOR_LENGTH_BYTES))) uint32_t ks0;
     __attribute__((aligned(VECTOR_LENGTH_BYTES))) uint32_t ks1;
     __attribute__((aligned(VECTOR_LENGTH_BYTES))) uint32_t ks2;
     __attribute__((aligned(VECTOR_LENGTH_BYTES))) uint32_t ks3;
     __attribute__((aligned(VECTOR_LENGTH_BYTES))) uint32_t ks4;

    //Use 1D arrays for aligned accesses when supply vector size
    //__attribute__((aligned(VECTOR_LENGTH_BYTES))) uint32_t X0[NUM_VALS_32];
    //__attribute__((aligned(VECTOR_LENGTH_BYTES))) uint32_t X1[NUM_VALS_32];
    //__attribute__((aligned(VECTOR_LENGTH_BYTES))) uint32_t X2[NUM_VALS_32];
    //__attribute__((aligned(VECTOR_LENGTH_BYTES))) uint32_t X3[NUM_VALS_32];

    __attribute__((aligned(VECTOR_LENGTH_BYTES))) uint32_t *X0;
    __attribute__((aligned(VECTOR_LENGTH_BYTES))) uint32_t *X1;
    __attribute__((aligned(VECTOR_LENGTH_BYTES))) uint32_t *X2;
    __attribute__((aligned(VECTOR_LENGTH_BYTES))) uint32_t *X3;
    
    
    //Test
    test_vec_sizes();
    

    X0 = (uint32_t*)malloc(NUM_VALS_32 * sizeof(uint32_t));
    X1 = (uint32_t*)malloc(NUM_VALS_32 * sizeof(uint32_t));
    X2 = (uint32_t*)malloc(NUM_VALS_32 * sizeof(uint32_t));
    X3 = (uint32_t*)malloc(NUM_VALS_32 * sizeof(uint32_t));
        
    //start = clock();
    //gettimeofday(&tv1, NULL);
    
    ks4 =  SKEIN_KS_PARITY32;

    //Assumes only first counter varies

    for (ivec=0;ivec <NUM_VALS_32 ; ivec++) {
      X0[ivec]  = ivec; // increment counter
      X1[ivec]  = 12345;
      X2[ivec]  = 12345;
      X3[ivec]  = 12345;               
    }
    
    //loop over 4 vals (not vector dimension)
    //for (ictr=0;ictr < 4; ictr++) {
    //  ks[ictr] = KEY4[ictr];
    //  ks[4] ^= KEY4[ictr]; 
    //}
    
    ks0 = 10;
    ks1 = 54321;
    ks2 = 54321;
    ks3 = 54321;
    
    //for (ictr=0;ictr < 4; ictr++) {
    //  ks[4] ^= ks[ictr]; 
   // }
   
   ks4 ^= ks0;
   ks4 ^= ks1;
   ks4 ^= ks2;
   ks4 ^= ks3;


    start = clock();
    gettimeofday(&tv1, NULL);

    //Ok so if use intrinsics need sep vector length from num vals
   
    //load X array vectors - can load all but X0 just once for threefry4x32f_multi_ss_fix 
//     __mVLENi tempX0 = _mmVLEN_load_siVLEN(&X0[ivec]);
    __mVLENi tempX1 = _mmVLEN_load_siVLEN(&X1[0]);
    __mVLENi tempX2 = _mmVLEN_load_siVLEN(&X2[0]);
    __mVLENi tempX3 = _mmVLEN_load_siVLEN(&X3[0]);
    
    //load ks vector
    __mVLENi tempks0 = _mmVLEN_set1_epi32(ks0);
    __mVLENi tempks1 = _mmVLEN_set1_epi32(ks1);
    __mVLENi tempks2 = _mmVLEN_set1_epi32(ks2);
    __mVLENi tempks3 = _mmVLEN_set1_epi32(ks3);
    __mVLENi tempks4 = _mmVLEN_set1_epi32(ks4);
    
    
    //Loop over vector length
    //May be benefit to unrolling (post-vec) to gain more independent ops (trade-off:reg pressure)
    #pragma omp simd aligned(X0,X1,X2,X3)
    for (ivec=0;ivec < NUM_VALS_32; ivec+=16) {
             
/*      
      X0[ivec] += ks[0];
      X1[ivec] += ks[1];
      X2[ivec] += ks[2];
      X3[ivec] += ks[3];
*/
      // load X array vectors - load all here for threefry4x32f_multi_ctrkey_all 
      __mVLENi tempX0 = _mmVLEN_load_siVLEN(&X0[ivec]);
//       __mVLENi tempX1 = _mmVLEN_load_siVLEN(&X1[ivec]);
//       __mVLENi tempX2 = _mmVLEN_load_siVLEN(&X2[ivec]);
//       __mVLENi tempX3 = _mmVLEN_load_siVLEN(&X3[ivec]);
//       
//       //load ks vector
//       __mVLENi tempks0 = _mmVLEN_set1_epi32(ks0);
//       __mVLENi tempks1 = _mmVLEN_set1_epi32(ks1);
//       __mVLENi tempks2 = _mmVLEN_set1_epi32(ks2);
//       __mVLENi tempks3 = _mmVLEN_set1_epi32(ks3);
//       __mVLENi tempks4 = _mmVLEN_set1_epi32(ks4);
      
      //-----------------------------------------------------
      //Injection 0. (from ks0)
      tempX0 = _mmVLEN_add_epi32(tempX0, tempks0);
      tempX1 = _mmVLEN_add_epi32(tempX1, tempks1);
      tempX2 = _mmVLEN_add_epi32(tempX2, tempks2);
      tempX3 = _mmVLEN_add_epi32(tempX3, tempks3);
      
      //sh - so far used 8 vector registers.
      //-----------------------------------------------------
      
      //Rounds 1-4
      halfround(tempX0,tempX1,R_32x4_0_0)
      halfround(tempX2,tempX3,R_32x4_0_1)
      
      halfround(tempX0,tempX3,R_32x4_1_0)
      halfround(tempX2,tempX1,R_32x4_1_1)
      
      halfround(tempX0,tempX1,R_32x4_2_0)
      halfround(tempX2,tempX3,R_32x4_2_1)
      
      halfround(tempX0,tempX3,R_32x4_3_0)
      halfround(tempX2,tempX1,R_32x4_3_1)
     
      //-----------------------------------------------------
      //Injection 1. (from ks1)
      tempX0 = _mmVLEN_add_epi32(tempX0, tempks1);
      tempX1 = _mmVLEN_add_epi32(tempX1, tempks2);
      tempX2 = _mmVLEN_add_epi32(tempX2, tempks3);
      tempX3 = _mmVLEN_add_epi32(tempX3, tempks4);
      
      //-----------------------------------------------------
      
      //Rounds 5-8
      halfround(tempX0,tempX1,R_32x4_4_0)
      halfround(tempX2,tempX3,R_32x4_4_1)
      
      halfround(tempX0,tempX3,R_32x4_5_0)
      halfround(tempX2,tempX1,R_32x4_5_1)
      
      halfround(tempX0,tempX1,R_32x4_6_0)
      halfround(tempX2,tempX3,R_32x4_6_1)
      
      halfround(tempX0,tempX3,R_32x4_7_0)
      halfround(tempX2,tempX1,R_32x4_7_1)

      //-----------------------------------------------------
      //Injection 2. (from ks2)
      tempX0 = _mmVLEN_add_epi32(tempX0, tempks2);
      tempX1 = _mmVLEN_add_epi32(tempX1, tempks3);
      tempX2 = _mmVLEN_add_epi32(tempX2, tempks4);
      tempX3 = _mmVLEN_add_epi32(tempX3, tempks0);
      
      //-----------------------------------------------------
      
      //Rounds 9-12
      halfround(tempX0,tempX1,R_32x4_0_0)
      halfround(tempX2,tempX3,R_32x4_0_1)
      
      halfround(tempX0,tempX3,R_32x4_1_0)
      halfround(tempX2,tempX1,R_32x4_1_1)
      
      halfround(tempX0,tempX1,R_32x4_2_0)
      halfround(tempX2,tempX3,R_32x4_2_1)
      
      halfround(tempX0,tempX3,R_32x4_3_0)
      halfround(tempX2,tempX1,R_32x4_3_1)
      
      //-----------------------------------------------------
      //Injection 3. (from ks3)
      tempX0 = _mmVLEN_add_epi32(tempX0, tempks3);
      tempX1 = _mmVLEN_add_epi32(tempX1, tempks4);
      tempX2 = _mmVLEN_add_epi32(tempX2, tempks0);
      tempX3 = _mmVLEN_add_epi32(tempX3, tempks1);
      
      //-----------------------------------------------------
      
      //Rounds 13-16
      halfround(tempX0,tempX1,R_32x4_4_0)
      halfround(tempX2,tempX3,R_32x4_4_1)
      
      halfround(tempX0,tempX3,R_32x4_5_0)
      halfround(tempX2,tempX1,R_32x4_5_1)
      
      halfround(tempX0,tempX1,R_32x4_6_0)
      halfround(tempX2,tempX3,R_32x4_6_1)
      
      halfround(tempX0,tempX3,R_32x4_7_0)
      halfround(tempX2,tempX1,R_32x4_7_1)
      
      //-----------------------------------------------------
      //Injection 4. (from ks4)
      tempX0 = _mmVLEN_add_epi32(tempX0, tempks4);
      tempX1 = _mmVLEN_add_epi32(tempX1, tempks0);
      tempX2 = _mmVLEN_add_epi32(tempX2, tempks1);
      tempX3 = _mmVLEN_add_epi32(tempX3, tempks2);
      
      //-----------------------------------------------------
      
      //Rounds 17-20
      halfround(tempX0,tempX1,R_32x4_0_0)
      halfround(tempX2,tempX3,R_32x4_0_1)
      
      halfround(tempX0,tempX3,R_32x4_1_0)
      halfround(tempX2,tempX1,R_32x4_1_1)
      
      halfround(tempX0,tempX1,R_32x4_2_0)
      halfround(tempX2,tempX3,R_32x4_2_1)
      
      halfround(tempX0,tempX3,R_32x4_3_0)
      halfround(tempX2,tempX1,R_32x4_3_1)
      
      //-----------------------------------------------------
    
      //Injection 5. (from ks0)
      tempX0 = _mmVLEN_add_epi32(tempX0, tempks0);
      tempX1 = _mmVLEN_add_epi32(tempX1, tempks1);
      tempX2 = _mmVLEN_add_epi32(tempX2, tempks2);
      tempX3 = _mmVLEN_add_epi32(tempX3, tempks3);
      
      //-----------------------------------------------------

//Store X values
      _mmVLEN_store_epi32(&X0[ivec],tempX0);
      _mmVLEN_store_epi32(&X1[ivec],tempX1);
      _mmVLEN_store_epi32(&X2[ivec],tempX2);
      _mmVLEN_store_epi32(&X3[ivec],tempX3);

    }

    diff = clock() - start;
    gettimeofday(&tv2, NULL);

    msec = diff * 1000 / CLOCKS_PER_SEC;
//    printf("Time taken %d seconds %d milliseconds\n", msec/1000, msec%1000);
    printf("Time taken (CPU)  = %d.%d seconds\n", msec/1000, msec%1000);
    
    printf ("Time taken (Wall) = %f seconds\n",
         (double) (tv2.tv_usec - tv1.tv_usec) / 1000000 +
         (double) (tv2.tv_sec - tv1.tv_sec));
    
    printf("X0[100]  = %u\n",X0[100]);
    printf("X3[1000] = %u\n",X3[1000]);
    printf("\nX2[NUM_VALS_32-1] = %u\n\n",X2[NUM_VALS_32-1]);
        
    free(X0);
    free(X1);
    free(X2);
    free(X3);
    

    return 0; 
}


//--------------------------------------------------------------------------------------------------
#else

int main (void)
{

    int ivec;    
    int ictr;  // For iteration over ctk/keys 1,2,3,4
    
    //clock - cpu timing
    clock_t start, diff;
    int msec;

    //gettimeofday - wall clock timing
    struct timeval  tv1, tv2;
    
    //uint32_t ks[4+1]; //Assuming one set of keys
    __attribute__((aligned(VECTOR_LENGTH_BYTES))) uint32_t ks[4+1]; //Assuming one set of keys
               
    //Use 1D arrays for aligned accesses when supply vector size
    //__attribute__((aligned(VECTOR_LENGTH_BYTES))) uint32_t X0[NUM_VALS_32];
    //__attribute__((aligned(VECTOR_LENGTH_BYTES))) uint32_t X1[NUM_VALS_32];
    //__attribute__((aligned(VECTOR_LENGTH_BYTES))) uint32_t X2[NUM_VALS_32];
    //__attribute__((aligned(VECTOR_LENGTH_BYTES))) uint32_t X3[NUM_VALS_32];

    __attribute__((aligned(VECTOR_LENGTH_BYTES))) uint32_t *X0;
    __attribute__((aligned(VECTOR_LENGTH_BYTES))) uint32_t *X1;
    __attribute__((aligned(VECTOR_LENGTH_BYTES))) uint32_t *X2;
    __attribute__((aligned(VECTOR_LENGTH_BYTES))) uint32_t *X3;
    
    X0 = (uint32_t*)malloc(NUM_VALS_32 * sizeof(uint32_t));
    X1 = (uint32_t*)malloc(NUM_VALS_32 * sizeof(uint32_t));
    X2 = (uint32_t*)malloc(NUM_VALS_32 * sizeof(uint32_t));
    X3 = (uint32_t*)malloc(NUM_VALS_32 * sizeof(uint32_t));
    
    
    //start = clock();
    //gettimeofday(&tv1, NULL);
    
    ks[4] =  SKEIN_KS_PARITY32;

    //Assumes only first counter varies

    for (ivec=0;ivec <NUM_VALS_32 ; ivec++) {
      X0[ivec]  = ivec; // increment counter
      X1[ivec]  = 12345;
      X2[ivec]  = 12345;
      X3[ivec]  = 12345;               
    }
    
    //loop over 4 vals (not vector dimension)
    //for (ictr=0;ictr < 4; ictr++) {
    //  ks[ictr] = KEY4[ictr];
    //  ks[4] ^= KEY4[ictr]; 
    //}
    
    ks[0] = 10;
    ks[1] = 54321;
    ks[2] = 54321;
    ks[3] = 54321;
    
    for (ictr=0;ictr < 4; ictr++) {
      ks[4] ^= ks[ictr]; 
    }
   
    start = clock();
    gettimeofday(&tv1, NULL);

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
      //buff[ivec*4+0] = X0[ivec]*R123_0x1p_32;
      //buff[ivec*4+1] = X1[ivec]*R123_0x1p_32;
      //buff[ivec*4+2] = X2[ivec]*R123_0x1p_32;
      //buff[ivec*4+3] = X3[ivec]*R123_0x1p_32;

    }

    diff = clock() - start;
    gettimeofday(&tv2, NULL);

    msec = diff * 1000 / CLOCKS_PER_SEC;
//    printf("Time taken %d seconds %d milliseconds\n", msec/1000, msec%1000);
    printf("Time taken (CPU)  = %d.%d seconds\n", msec/1000, msec%1000);
    
    printf ("Time taken (Wall) = %f seconds\n",
         (double) (tv2.tv_usec - tv1.tv_usec) / 1000000 +
         (double) (tv2.tv_sec - tv1.tv_sec));
    
    printf("X0[100]  = %u\n",X0[100]);
    printf("X3[1000] = %u\n",X3[1000]);
    printf("\nX2[NUM_VALS_32-1] = %u\n\n",X2[NUM_VALS_32-1]);
        
    free(X0);
    free(X1);
    free(X2);
    free(X3);
    

    return 0; 
}


#endif







