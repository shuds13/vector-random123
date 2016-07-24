//--------------------------------------------------------------------------------------------------
// S. Hudson: Vector123 loop only test
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
#define NUM_VALS_32 1600000
#define NUM_VALS_64 800000

#define SKEIN_MK_64(hi32,lo32)  ((lo32) + (((uint64_t) (hi32)) << 32))
#define SKEIN_KS_PARITY64         SKEIN_MK_64(0x1BD11BDA,0xA9FC1A22)
#define SKEIN_KS_PARITY32         0x1BD11BDA

#define R123_0x1p_32                (1./4294967296.)   // for 32_53 CO,OC,OO int4 to d.p


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
static inline __m512i _mm512_myrol_epi32(__m512i const & epi32A,unsigned int nCount)
{
	__m512i const epi32H = _mm512_slli_epi32(epi32A, nCount);
	__m512i const epi32L = _mm512_srli_epi32(epi32A, 32 - nCount);
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
    
    clock_t start, diff;
    int msec;
    
    //*check using omp to tell compiler can assume aligned on the loop
    // so really want to use omp to force alignement on declaration
    
    
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

//*check - example code - remove
/*
	__m256i zero;
	zero = _mm256_xor_si256(zero,zero);
	for(int i=0;i<count;++i) {
	  __m256i temp = _mm256_load_si256(&array[i]);
	  __m256i NegCarry = _mm256_cmpgt_epi32(zero,temp);
	  _mm256_store_si256(&array[i], _mm256_sub_epi32(_mm256_add_epi32(temp, temp), NegCarry));
	}
*/

    //Ok so if use intrinsics need sep vector length from num vals
    
    //loop over vector length
    #pragma omp simd aligned(X0,X1,X2,X3,ks)         
    for (ivec=0;ivec < NUM_VALS_32; ivec+=16) {
    
    // Replace 256 with 512 equiv where nec. (maybe can do with 128 instructions also?)
    // Where poss could use a macro to define length - 256/512 etc...
    // under #ifdef - then #ifdef directly for the rest
    // Also could cert do this with auto-detection - ifdef __AVX512F__ ETC..
    // maybe inc knc - looks like that doesnt have the rol though - though does landing
    // really have it?
    // Also could even combine with portable code - if detects AVX2/512 ETC
    // it runs the "special" version!!!
    
    //*check
    // Also i've written ev out in general code version - but this would get really
    // long here - so could just say every line is basically a function/tempate (inline)
    // inputs are: 2 X values and one rotation constant
    // then just key injections - i could try doing that with the original code
    // see if performs same - and whether same clarity for profiling/debugging etc..
    
/*      
      X0[ivec] += ks[0];
      X1[ivec] += ks[1];
      X2[ivec] += ks[2];
      X3[ivec] += ks[3];
*/

      // load X array vectors   
      __m512i tempX0 = _mm512_load_si512(&X0[i]);
      __m512i tempX1 = _mm512_load_si512(&X1[i]);
      __m512i tempX2 = _mm512_load_si512(&X2[i]);
      __m512i tempX3 = _mm512_load_si512(&X3[i]);
      
      //load ks vector
      //*check - do i have to load and broadcast???
      
      //__m512i tempks = _mm512_load_si512(&ks[i]);
      //__m512i tempks _mm512_broadcastd_epi32(__m128i a)
      __m512i tempks0 = _mm512_broadcastd_epi32(&ks[0]);
      __m512i tempks1 = _mm512_broadcastd_epi32(&ks[1]);
      __m512i tempks2 = _mm512_broadcastd_epi32(&ks[2]);
      __m512i tempks3 = _mm512_broadcastd_epi32(&ks[3]);
      
      tempX0 = _mm512_add_epi32(tempX0, tempks0);
      tempX1 = _mm512_add_epi32(tempX1, tempks1);
      tempX2 = _mm512_add_epi32(tempX2, tempks2);
      tempX3 = _mm512_add_epi32(tempX3, tempks3);
      
      //sh - so far used 8 vector registers.
      
      //X0[ivec] += X1[ivec];
      tempX0 = _mm512_add_epi32(tempX0, tempX1);
            
      //X1[ivec] = RotL_32(X1[ivec],R_32x4_0_0);
#if defined(__AVX512F__)
     // AVX512 introduces vectorised rotate left intrinsic
     //__m512i _mm512_rol_epi32 (__m512i a, const int imm8);
     
     //check* works - if not try just const - rather than enum type...
     __m512i _mm512_rol_epi32 (tempX1, R_32x4_0_0); 
#else
     //Use shift left/shift right and OR
     __m512i _mm512_myrol_epi32(tempX1, R_32x4_0_0);      
#endif      
      
      
      //X1[ivec] ^= X0[ivec];
      tempX1 = _mm512_xor_epi32(tempX1, tempX0)
       


    }

    diff = clock() - start;

    msec = diff * 1000 / CLOCKS_PER_SEC;
//    printf("Time taken %d seconds %d milliseconds\n", msec/1000, msec%1000);
    printf("Time taken %d.%d seconds\n", msec/1000, msec%1000);
    printf("X0[100]  = %d\n",X0[100]);
    printf("X3[1000] = %d\n",X3[1000]);
    
    free(X0);
    free(X1);
    free(X2);
    free(X3);
    

    return 0; 
}
