//--------------------------------------------------------------------------------------------------
// S. Hudson: Vector123 loop-only test
//            Apply vectorizable threefry4x32 to generate large array of random floats between 0 and 1
//            Aims to test performance of pure conversion loop - scalar v vectorized
//            Number of values produced is 4*NUM_VALS_32
//            Standard source code version.
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
#include <sys/time.h>

// Allows use of short forms like uint32_t
#include <stdint.h>


//Random123 - To pick up macros from R123 - uncomment
//#include <threefry.h>
//#include <philox.h>

//#include <u01.h>
//OR
//#include <u01fixedpt.h>


#define VECTOR_LENGTH_BYTES 64
#define NUM_VALS_32 160000000
#define NUM_VALS_64 80000000

#define SKEIN_MK_64(hi32,lo32)  ((lo32) + (((uint64_t) (hi32)) << 32))
#define SKEIN_KS_PARITY64         SKEIN_MK_64(0x1BD11BDA,0xA9FC1A22)
#define SKEIN_KS_PARITY32         0x1BD11BDA

//#define R123_0x1p_32                (1./4294967296.)   // for 32_53 CO,OC,OO int4 to d.p
#define R123_0x1p_32f               (1.f/4294967296.f) //float version for CC
#define R123_0x1p_24f               (1.f/16777216.f)   //float version for CO


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

    __attribute__((aligned(VECTOR_LENGTH_BYTES))) float *buff;
    

    X0 = (uint32_t*)malloc(NUM_VALS_32 * sizeof(uint32_t));
    X1 = (uint32_t*)malloc(NUM_VALS_32 * sizeof(uint32_t));
    X2 = (uint32_t*)malloc(NUM_VALS_32 * sizeof(uint32_t));
    X3 = (uint32_t*)malloc(NUM_VALS_32 * sizeof(uint32_t));
    
    buff = (float*)malloc(NUM_VALS_32 * 4 * sizeof(float));
    
    //start = clock();
    //gettimeofday(&tv1, NULL);
    
    ks[4] =  SKEIN_KS_PARITY32;

    //Assumes only first counter varies

    #pragma omp parallel for default(none), private(ivec), shared(X0,X1,X2,X3)
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
    rdtsc_count1 = rdtsc();

    //loop over vector length
    //#pragma omp simd aligned(X0,X1,X2,X3,ks)         
    #pragma omp parallel for default(none), private(ivec), shared(X0,X1,X2,X3,ks,buff)
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


      //Convert integers to floats uniformly distributed between 0 and 1 inclusive/exlusive
      //u01_closed_open_32_53
      //buff[ivec*4+0] = X0[ivec]*R123_0x1p_32;
      //buff[ivec*4+1] = X1[ivec]*R123_0x1p_32;
      //buff[ivec*4+2] = X2[ivec]*R123_0x1p_32;
      //buff[ivec*4+3] = X3[ivec]*R123_0x1p_32;
      
      //Vector order stores - u01_closed_closed_32_24
      //buff[NUM_VALS_32*0 + ivec] = X0[ivec]*R123_0x1p_32f;
      //buff[NUM_VALS_32*1 + ivec] = X1[ivec]*R123_0x1p_32f;
      //buff[NUM_VALS_32*2 + ivec] = X2[ivec]*R123_0x1p_32f;
      //buff[NUM_VALS_32*3 + ivec] = X3[ivec]*R123_0x1p_32f;
      
      //Vector order stores - u01_closed_open_32_24
      buff[NUM_VALS_32*0 + ivec] = (X0[ivec]>>8)*R123_0x1p_24f;
      buff[NUM_VALS_32*1 + ivec] = (X1[ivec]>>8)*R123_0x1p_24f;
      buff[NUM_VALS_32*2 + ivec] = (X2[ivec]>>8)*R123_0x1p_24f;
      buff[NUM_VALS_32*3 + ivec] = (X3[ivec]>>8)*R123_0x1p_24f;

    }

    rdtsc_count2 = rdtsc();

    diff = clock() - start;
    gettimeofday(&tv2, NULL);

#pragma omp parallel
#pragma omp master
    printf("Num threads:  %d\n",omp_get_num_threads());

    printf("Num sets: %d   Tot. num randoms: %d\n",NUM_VALS_32,NUM_VALS_32*4);
    printf("cycles: %llu\n",rdtsc_count2 - rdtsc_count1);

    msec = diff * 1000 / CLOCKS_PER_SEC;
//    printf("Time taken %d seconds %d milliseconds\n", msec/1000, msec%1000);
    printf("Time taken (CPU)  = %d.%d seconds\n", msec/1000, msec%1000);
    
    printf ("Time taken (Wall) = %f seconds\n",
         (double) (tv2.tv_usec - tv1.tv_usec) / 1000000 +
         (double) (tv2.tv_sec - tv1.tv_sec));
    
    printf("buff[100]  = %f\n",buff[100]);
    printf("buff[1000] = %f\n",buff[1000]);
    printf("buff[3*NUM_VALS_32-1] = %f\n\n",buff[3*NUM_VALS_32-1]);
        
    free(X0);
    free(X1);
    free(X2);
    free(X3);
    free(buff);    
    

    return 0; 
}
