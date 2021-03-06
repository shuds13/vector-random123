//--------------------------------------------------------------------------------------------------
// S. Hudson: Vector123 loop-only test - using AVX512 instrinsics
//            Apply vectorizable threefry4x32 to generate large array of random floats between 0 and 1
//            Aims to test performance of threefry and conversion to float - scalar v vectorized
//            Number of values produced is 4*NUM_VALS_32
//            AVX512 vector intrinsics version.
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


#define VECTOR_LENGTH_BYTES 64
#define NUM_VALS_32 16000000
#define NUM_VALS_64 8000000

#define SKEIN_MK_64(hi32,lo32)  ((lo32) + (((uint64_t) (hi32)) << 32))
#define SKEIN_KS_PARITY64       SKEIN_MK_64(0x1BD11BDA,0xA9FC1A22)
#define SKEIN_KS_PARITY32       0x1BD11BDA


//#define R123_0x1p_32                (1./4294967296.)   // for 32_53 CO,OC,OO int4 to d.p
#define R123_0x1p_32f               (1.f/4294967296.f) //float version for CC
#define R123_0x1p_24f               (1.f/16777216.f)   //float version for CO



//*
//Generate vector rol instruction - fast
#define halfround(X,Y,ROL_VAL) \
      X = _mm512_add_epi32(X, Y); \
      Y = _mm512_rol_epi32(Y, ROL_VAL); \
      Y = _mm512_xor_epi32(Y, X);
//*/

/*
//Use a composed rol - slow
#define halfround(X,Y,ROL_VAL) \
      X = _mm512_add_epi32(X, Y); \
      Y = _mm512_myrol_epi32(Y, ROL_VAL); \
      Y = _mm512_xor_epi32(Y, X);
//*/


static inline uint64_t rdtsc(){
  unsigned int lo,hi;
  __asm__ __volatile__ ("rdtsc" : "=a" (lo) "=d" (hi));
  return ((uint64_t)hi << 32) | lo;
}


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
//__m512i _mm512_myrol_epi32(__m512i const & epi32A)
static inline __m512i _mm512_myrol_epi32(__m512i epi32A,unsigned int nCount)
{
	__m512i epi32H = _mm512_slli_epi32(epi32A, nCount);
	__m512i epi32L = _mm512_srli_epi32(epi32A, 32 - nCount);
	return _mm512_or_epi32(epi32H, epi32L);
	//return _mm512_or_si512(epi32AH, epi32L);	//what would be the difference?
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

int main (void)
{

    int ivec;    
    int ictr;  // For iteration over ctk/keys 1,2,3,4
    
    //clock - cpu timing
    clock_t start, diff;
    int msec;
    uint64_t rdtsc_count1, rdtsc_count2;

    //gettimeofday - wall clock timing
    struct timeval  tv1, tv2;
    
    //*check using omp to tell compiler can assume aligned on the loop
    // so really want to use omp to force alignment on declaration
        
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

    //__attribute__((aligned(VECTOR_LENGTH_BYTES))) float *buff0;
    //__attribute__((aligned(VECTOR_LENGTH_BYTES))) float *buff1;
    //__attribute__((aligned(VECTOR_LENGTH_BYTES))) float *buff2;
    //__attribute__((aligned(VECTOR_LENGTH_BYTES))) float *buff3;

    __attribute__((aligned(VECTOR_LENGTH_BYTES))) float *buff;

    X0 = (uint32_t*)malloc(NUM_VALS_32 * sizeof(uint32_t));
    X1 = (uint32_t*)malloc(NUM_VALS_32 * sizeof(uint32_t));
    X2 = (uint32_t*)malloc(NUM_VALS_32 * sizeof(uint32_t));
    X3 = (uint32_t*)malloc(NUM_VALS_32 * sizeof(uint32_t));
    
    //buff0 = (float*)malloc(NUM_VALS_32 * sizeof(float));
    //buff1 = (float*)malloc(NUM_VALS_32 * sizeof(float));
    //buff2 = (float*)malloc(NUM_VALS_32 * sizeof(float));
    //buff3 = (float*)malloc(NUM_VALS_32 * sizeof(float));
    
    buff = (float*)malloc(NUM_VALS_32 * 4 * sizeof(float));
	
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
    rdtsc_count1 = rdtsc();

    //Ok so if use intrinsics need sep vector length from num vals
    
    //loop over vector length        
    for (ivec=0;ivec < NUM_VALS_32; ivec+=16) {
             
/*      
      X0[ivec] += ks[0];
      X1[ivec] += ks[1];
      X2[ivec] += ks[2];
      X3[ivec] += ks[3];
*/
      // load X array vectors   
      __m512i tempX0 = _mm512_load_si512(&X0[ivec]);
      __m512i tempX1 = _mm512_load_si512(&X1[ivec]);
      __m512i tempX2 = _mm512_load_si512(&X2[ivec]);
      __m512i tempX3 = _mm512_load_si512(&X3[ivec]);
      
      //load ks vector
      __m512i tempks0 = _mm512_set1_epi32(ks0);
      __m512i tempks1 = _mm512_set1_epi32(ks1);
      __m512i tempks2 = _mm512_set1_epi32(ks2);
      __m512i tempks3 = _mm512_set1_epi32(ks3);
      __m512i tempks4 = _mm512_set1_epi32(ks4);
      
      //__m512i tempks1 = _mm512_broadcastd_epi32(&ks[1]);
      //__m512i tempks2 = _mm512_broadcastd_epi32(&ks[2]);
      //__m512i tempks3 = _mm512_broadcastd_epi32(&ks[3]);

      //-----------------------------------------------------
      //Injection 0. (from ks0)
      tempX0 = _mm512_add_epi32(tempX0, tempks0);
      tempX1 = _mm512_add_epi32(tempX1, tempks1);
      tempX2 = _mm512_add_epi32(tempX2, tempks2);
      tempX3 = _mm512_add_epi32(tempX3, tempks3);
      
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
      tempX0 = _mm512_add_epi32(tempX0, tempks1);
      tempX1 = _mm512_add_epi32(tempX1, tempks2);
      tempX2 = _mm512_add_epi32(tempX2, tempks3);
      tempX3 = _mm512_add_epi32(tempX3, tempks4);

      __m512i incrm = _mm512_set1_epi32(1);
      tempX3 = _mm512_add_epi32(tempX3, incrm);
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
      tempX0 = _mm512_add_epi32(tempX0, tempks2);
      tempX1 = _mm512_add_epi32(tempX1, tempks3);
      tempX2 = _mm512_add_epi32(tempX2, tempks4);
      tempX3 = _mm512_add_epi32(tempX3, tempks0);
	    
      __m512i incrm = _mm512_set1_epi32(2);
      tempX3 = _mm512_add_epi32(tempX3, incrm);	    
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
      tempX0 = _mm512_add_epi32(tempX0, tempks3);
      tempX1 = _mm512_add_epi32(tempX1, tempks4);
      tempX2 = _mm512_add_epi32(tempX2, tempks0);
      tempX3 = _mm512_add_epi32(tempX3, tempks1);
	    
      __m512i incrm = _mm512_set1_epi32(3);
      tempX3 = _mm512_add_epi32(tempX3, incrm);	    
      
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
      tempX0 = _mm512_add_epi32(tempX0, tempks4);
      tempX1 = _mm512_add_epi32(tempX1, tempks0);
      tempX2 = _mm512_add_epi32(tempX2, tempks1);
      tempX3 = _mm512_add_epi32(tempX3, tempks2);
	    
      __m512i incrm = _mm512_set1_epi32(4);
      tempX3 = _mm512_add_epi32(tempX3, incrm);	    
      
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
      tempX0 = _mm512_add_epi32(tempX0, tempks0);
      tempX1 = _mm512_add_epi32(tempX1, tempks1);
      tempX2 = _mm512_add_epi32(tempX2, tempks2);
      tempX3 = _mm512_add_epi32(tempX3, tempks3);
	    
      __m512i incrm = _mm512_set1_epi32(5);
      tempX3 = _mm512_add_epi32(tempX3, incrm);	    
      
      //-----------------------------------------------------
      //x0[ivec]>>8
      tempX0 = _mm512_srli_epi32(tempX0, 8);
      tempX1 = _mm512_srli_epi32(tempX1, 8);
      tempX2 = _mm512_srli_epi32(tempX2, 8);
      tempX3 = _mm512_srli_epi32(tempX3, 8);
      
      //x0[ivec]>>8 convert to float
      __mm512 temp_buff0 = _mm512_cvtepi32_ps(tempX0);
      __mm512 temp_buff1 = _mm512_cvtepi32_ps(tempX1);
      __mm512 temp_buff2 = _mm512_cvtepi32_ps(tempX2);
      __mm512 temp_buff3 = _mm512_cvtepi32_ps(tempX3);
	    
      __mm512 temp_vec_0x1p_24f = _mm512_set1_ps(R123_0x1p_24f);
      
      //(x0[ivec]>>8)*R123_0x1p_24f
      temp_buff0 = _mm512_mul_ps(tempX0,temp_vec_0x1p_24f);
      temp_buff1 = _mm512_mul_ps(tempX1,temp_vec_0x1p_24f);     
      temp_buff2 = _mm512_mul_ps(tempX2,temp_vec_0x1p_24f);
      temp_buff3 = _mm512_mul_ps(tempX3,temp_vec_0x1p_24f);

      _mm512_storeu_ps(&buff[ivec])
	      
      //Store floats in buff     
      _mm512_storeu_ps(&buff[NUM_VALS_32*0 + ivec],temp_buff0);
      _mm512_storeu_ps(&buff[NUM_VALS_32*1 + ivec],temp_buff1);
      _mm512_storeu_ps(&buff[NUM_VALS_32*2 + ivec],temp_buff2);
      _mm512_storeu_ps(&buff[NUM_VALS_32*3 + ivec],temp_buff3);

    }

    rdtsc_count2 = rdtsc();

    diff = clock() - start;
    gettimeofday(&tv2, NULL);
    
    printf("Num sets: %d   Tot. num randoms: %d\n",NUM_VALS_32,NUM_VALS_32*4);
    printf("cycles: %llu\n",rdtsc_count2 - rdtsc_count1);

    msec = diff * 1000 / CLOCKS_PER_SEC;
    printf("Time taken (CPU)  = %d.%d seconds\n", msec/1000, msec%1000);
    
    printf ("Time taken (Wall) = %f seconds\n",
         (double) (tv2.tv_usec - tv1.tv_usec) / 1000000 +
         (double) (tv2.tv_sec - tv1.tv_sec));
      
    printf("buff[100]  = %f\n",buff[100]);
    printf("buff[1000] = %f\n",buff[1000]);
    printf("\nbuff[3*NUM_VALS_32-1] = %f\n\n",buff[3*NUM_VALS_32-1]);
        
    free(X0);
    free(X1);
    free(X2);
    free(X3);

    free(buff);
	
    //free(buff0);    
    //free(buff1);    
    //free(buff2);    
    //free(buff3);    
    

    return 0; 
}
