//#include <Random123/threefry.h>

#include <stdio.h>

//Random123
#include <threefry.h>
#include <philox.h>
#include <u01.h>

//=========================================================================================================
//SINGLE CALL VERSIONS
//=========================================================================================================

//extern "C"
//int threefry4x32f(unsigned long long int* SED4,unsigned long long int* KEY4, double* buff) {
int threefry4x32f(unsigned int* SED4,unsigned int* KEY4, double* buff) {
  
    threefry4x32_ctr_t rand_int;

    threefry4x32_ctr_t ctr; 
    threefry4x32_key_t key;
    
    ctr.v[0] = *SED4; 
    ctr.v[1] = *(SED4+1);
    ctr.v[2] = *(SED4+2); 
    ctr.v[3] = *(SED4+3);
    
    key.v[0] = *KEY4; 
    key.v[1] = *(KEY4+1);
    key.v[2] = *(KEY4+2); 
    key.v[3] = *(KEY4+3);
    
    
    //Produce 4 randomised integers
    rand_int = threefry4x32(ctr, key);
    
    //Convert integers to doubles uniformly distributed between 0 and 1 inclusive/exlusive
    //Note: Second value is size of mantissa: 24 for s.p or 53 for d.p)
    //Note: In version where use only first value - only first of these calls is strictly necessary
    *buff     = u01_closed_open_32_53(rand_int.v[0]);
    *(buff+1) = u01_closed_open_32_53(rand_int.v[1]);
    *(buff+2) = u01_closed_open_32_53(rand_int.v[2]);
    *(buff+3) = u01_closed_open_32_53(rand_int.v[3]);


  // Again Fortran will expect a pointer to the value
  //return randval; 
  return 0; 
}

//extern "C"
//int threefry4x64f(unsigned long long int* SED4,unsigned long long int* KEY4, double* buff) {
int threefry4x64f(unsigned int* SED4,unsigned int* KEY4, double* buff) {
  
    threefry4x64_ctr_t rand_int;

    threefry4x64_ctr_t ctr; 
    threefry4x64_key_t key;
    
    ctr.v[0] = *SED4; 
    ctr.v[1] = *(SED4+1);
    ctr.v[2] = *(SED4+2); 
    ctr.v[3] = *(SED4+3);
    
    key.v[0] = *KEY4; 
    key.v[1] = *(KEY4+1);
    key.v[2] = *(KEY4+2); 
    key.v[3] = *(KEY4+3);
    
    
    //Produce 4 randomised integers
    rand_int = threefry4x64(ctr, key);
    
    //Convert integers to doubles uniformly distributed between 0 and 1 inclusive/exlusive
    //Note: Second value is size of mantissa: 24 for s.p or 53 for d.p)
    //Note: In version where use only first value - only first of these calls is strictly necessary
    *buff     = u01_closed_open_64_53(rand_int.v[0]);
    *(buff+1) = u01_closed_open_64_53(rand_int.v[1]);
    *(buff+2) = u01_closed_open_64_53(rand_int.v[2]);
    *(buff+3) = u01_closed_open_64_53(rand_int.v[3]);


  // Again Fortran will expect a pointer to the value
  //return randval; 
  return 0; 
}

//extern "C"
//int philox4x64f(unsigned long long int* SED4,unsigned long long int* KEY4, double* buff) {
int philox4x64f(unsigned int* SED4,unsigned int* KEY4, double* buff) {
  
    philox4x64_ctr_t rand_int;

    philox4x64_ctr_t ctr; 
    philox4x64_key_t key;
    
    ctr.v[0] = *SED4; 
    ctr.v[1] = *(SED4+1);
    ctr.v[2] = *(SED4+2); 
    ctr.v[3] = *(SED4+3);
    
    key.v[0] = *KEY4; 
    key.v[1] = *(KEY4+1);
    
    // *note philox uses only 2 keys
    //key.v[2] = *(KEY4+2); 
    //key.v[3] = *(KEY4+3);
    
    
    //Produce 4 randomised integers
    rand_int = philox4x64(ctr, key);
    
    //Convert integers to doubles uniformly distributed between 0 and 1 inclusive/exlusive
    //Note: Second value is size of mantissa: 24 for s.p or 53 for d.p)
    //Note: In version where use only first value - only first of these calls is strictly necessary
    *buff     = u01_closed_open_64_53(rand_int.v[0]);
    *(buff+1) = u01_closed_open_64_53(rand_int.v[1]);
    *(buff+2) = u01_closed_open_64_53(rand_int.v[2]);
    *(buff+3) = u01_closed_open_64_53(rand_int.v[3]);


  // Again Fortran will expect a pointer to the value
  //return randval; 
  return 0; 
}



//=========================================================================================================
//MULTI CALL VERSIONS
//=========================================================================================================

//extern "C"
//int threefry4x64f_multi_(unsigned long long int* SED4,unsigned long long int* KEY4, double* buff) {
int threefry4x64f_multi_(unsigned int* SED4,unsigned int* KEY4, double* buff) {
  
    threefry4x64_ctr_t rand_int1;
    threefry4x64_ctr_t rand_int2;
    threefry4x64_ctr_t rand_int3;

    threefry4x64_ctr_t ctr1; 
    threefry4x64_ctr_t ctr2; 
    threefry4x64_ctr_t ctr3; 

    //For this I can use same key
    threefry4x64_key_t key;
    
    ctr1.v[0] = *SED4; 
    ctr1.v[1] = *(SED4+1);
    ctr1.v[2] = *(SED4+2); 
    ctr1.v[3] = *(SED4+3);

    ctr2.v[0] = *SED4 + 4; //Add 4 to counter (to match one call version - rem only using ev 4th counter)
    ctr2.v[1] = *(SED4+1);
    ctr2.v[2] = *(SED4+2); 
    ctr2.v[3] = *(SED4+3);

    ctr3.v[0] = *SED4 + 8;
    ctr3.v[1] = *(SED4+1);
    ctr3.v[2] = *(SED4+2); 
    ctr3.v[3] = *(SED4+3);

    
    key.v[0] = *KEY4; 
    key.v[1] = *(KEY4+1);
    key.v[2] = *(KEY4+2); 
    key.v[3] = *(KEY4+3);
    
    
    //Produce 4 randomised integers per call
    rand_int1 = threefry4x64(ctr1, key);
    rand_int2 = threefry4x64(ctr2, key);
    rand_int3 = threefry4x64(ctr3, key);
    
    //Convert integers to doubles uniformly distributed between 0 and 1 inclusive/exlusive
    //Note: Second value is size of mantissa: 24 for s.p or 53 for d.p)
    //Note: In version where use only first value - only first of these calls is strictly necessary
    *buff     = u01_closed_open_64_53(rand_int1.v[0]);
    *(buff+1) = u01_closed_open_64_53(rand_int1.v[1]);
    *(buff+2) = u01_closed_open_64_53(rand_int1.v[2]);
    *(buff+3) = u01_closed_open_64_53(rand_int1.v[3]);

    *(buff+4) = u01_closed_open_64_53(rand_int2.v[0]);
    *(buff+5) = u01_closed_open_64_53(rand_int2.v[1]);
    *(buff+6) = u01_closed_open_64_53(rand_int2.v[2]);
    *(buff+7) = u01_closed_open_64_53(rand_int2.v[3]);

    *(buff+8)  = u01_closed_open_64_53(rand_int3.v[0]);
    *(buff+9)  = u01_closed_open_64_53(rand_int3.v[1]);
    *(buff+10) = u01_closed_open_64_53(rand_int3.v[2]);
    *(buff+11) = u01_closed_open_64_53(rand_int3.v[3]);

  // Again Fortran will expect a pointer to the value
  //return randval; 
  return 0; 
}

//extern "C"
//int philox4x64f_multi_(unsigned long long int* SED4,unsigned long long int* KEY4, double* buff) {
int philox4x64f_multi_(unsigned int* SED4,unsigned int* KEY4, double* buff) {
  
    philox4x64_ctr_t rand_int1;
    philox4x64_ctr_t rand_int2;
    philox4x64_ctr_t rand_int3;

    philox4x64_ctr_t ctr1; 
    philox4x64_ctr_t ctr2; 
    philox4x64_ctr_t ctr3; 

    //For this I can use same key
    philox4x64_key_t key;
    
    ctr1.v[0] = *SED4; 
    ctr1.v[1] = *(SED4+1);
    ctr1.v[2] = *(SED4+2); 
    ctr1.v[3] = *(SED4+3);

    ctr2.v[0] = *SED4 + 4; //Add 4 to counter (to match one call version - rem only using ev 4th counter)
    ctr2.v[1] = *(SED4+1);
    ctr2.v[2] = *(SED4+2); 
    ctr2.v[3] = *(SED4+3);

    ctr3.v[0] = *SED4 + 8;
    ctr3.v[1] = *(SED4+1);
    ctr3.v[2] = *(SED4+2); 
    ctr3.v[3] = *(SED4+3);
    
    key.v[0] = *KEY4; 
    key.v[1] = *(KEY4+1);
    
    // *note philox uses only 2 keys
    //key.v[2] = *(KEY4+2); 
    //key.v[3] = *(KEY4+3);
    
    
    //Produce 4 randomised integers
    rand_int1 = philox4x64(ctr1, key);
    rand_int2 = philox4x64(ctr2, key);
    rand_int3 = philox4x64(ctr3, key);
    
    //Convert integers to doubles uniformly distributed between 0 and 1 inclusive/exlusive
    //Note: Second value is size of mantissa: 24 for s.p or 53 for d.p)
    //Note: In version where use only first value - only first of these calls is strictly necessary
    *buff     = u01_closed_open_64_53(rand_int1.v[0]);
    *(buff+1) = u01_closed_open_64_53(rand_int1.v[1]);
    *(buff+2) = u01_closed_open_64_53(rand_int1.v[2]);
    *(buff+3) = u01_closed_open_64_53(rand_int1.v[3]);

    *(buff+4) = u01_closed_open_64_53(rand_int2.v[0]);
    *(buff+5) = u01_closed_open_64_53(rand_int2.v[1]);
    *(buff+6) = u01_closed_open_64_53(rand_int2.v[2]);
    *(buff+7) = u01_closed_open_64_53(rand_int2.v[3]);

    *(buff+8)  = u01_closed_open_64_53(rand_int3.v[0]);
    *(buff+9)  = u01_closed_open_64_53(rand_int3.v[1]);
    *(buff+10) = u01_closed_open_64_53(rand_int3.v[2]);
    *(buff+11) = u01_closed_open_64_53(rand_int3.v[3]);


  // Again Fortran will expect a pointer to the value
  //return randval; 
  return 0; 
}





//On some platforms (power8) philox4x64 is not working.
//Try 4x32

//extern "C"
//int philox4x32f(unsigned long long int* SED4,unsigned long long int* KEY4, double* buff) {
int philox4x32f(unsigned int* SED4,unsigned int* KEY4, double* buff) {
  
    philox4x32_ctr_t rand_int;

    philox4x32_ctr_t ctr; 
    philox4x32_key_t key;
    
    ctr.v[0] = *SED4; 
    ctr.v[1] = *(SED4+1);
    ctr.v[2] = *(SED4+2); 
    ctr.v[3] = *(SED4+3);
    
    key.v[0] = *KEY4; 
    key.v[1] = *(KEY4+1);
    
    // *note philox uses only 2 keys
    //key.v[2] = *(KEY4+2); 
    //key.v[3] = *(KEY4+3);
    
    
    //Produce 4 randomised integers
    rand_int = philox4x32(ctr, key);
    
    //Convert integers to doubles uniformly distributed between 0 and 1 inclusive/exlusive
    //Note: Second value is size of mantissa: 24 for s.p or 53 for d.p)
    //Note: In version where use only first value - only first of these calls is strictly necessary
    *buff     = u01_closed_open_64_53(rand_int.v[0]);
    *(buff+1) = u01_closed_open_64_53(rand_int.v[1]);
    *(buff+2) = u01_closed_open_64_53(rand_int.v[2]);
    *(buff+3) = u01_closed_open_64_53(rand_int.v[3]);


  // Again Fortran will expect a pointer to the value
  //return randval; 
  return 0; 
}



//extern "C"
//int philox4x32f_multi_(unsigned long long int* SED4,unsigned long long int* KEY4, double* buff) {
int philox4x32f_multi_(unsigned int* SED4,unsigned int* KEY4, double* buff) {
  
    philox4x32_ctr_t rand_int1;
    philox4x32_ctr_t rand_int2;
    philox4x32_ctr_t rand_int3;

    philox4x32_ctr_t ctr1; 
    philox4x32_ctr_t ctr2; 
    philox4x32_ctr_t ctr3; 

    //For this I can use same key
    philox4x32_key_t key;
    
    ctr1.v[0] = *SED4; 
    ctr1.v[1] = *(SED4+1);
    ctr1.v[2] = *(SED4+2); 
    ctr1.v[3] = *(SED4+3);

    ctr2.v[0] = *SED4 + 4; //Add 4 to counter (to match one call version - rem only using ev 4th counter)
    ctr2.v[1] = *(SED4+1);
    ctr2.v[2] = *(SED4+2); 
    ctr2.v[3] = *(SED4+3);

    ctr3.v[0] = *SED4 + 8;
    ctr3.v[1] = *(SED4+1);
    ctr3.v[2] = *(SED4+2); 
    ctr3.v[3] = *(SED4+3);
    
    key.v[0] = *KEY4; 
    key.v[1] = *(KEY4+1);
    
    // *note philox uses only 2 keys
    //key.v[2] = *(KEY4+2); 
    //key.v[3] = *(KEY4+3);
    
    
    //Produce 4 randomised integers
//     rand_int1 = philox4x32(ctr1, key);
//     rand_int2 = philox4x32(ctr2, key);
//     rand_int3 = philox4x32(ctr3, key);

//4z64 not working on power8 - try 4x32
    rand_int1 = philox4x32(ctr1, key);
    rand_int2 = philox4x32(ctr2, key);
    rand_int3 = philox4x32(ctr3, key);
    
    //Convert integers to doubles uniformly distributed between 0 and 1 inclusive/exlusive
    //Note: Second value is size of mantissa: 24 for s.p or 53 for d.p)
    //Note: In version where use only first value - only first of these calls is strictly necessary
    *buff     = u01_closed_open_64_53(rand_int1.v[0]);
    *(buff+1) = u01_closed_open_64_53(rand_int1.v[1]);
    *(buff+2) = u01_closed_open_64_53(rand_int1.v[2]);
    *(buff+3) = u01_closed_open_64_53(rand_int1.v[3]);

    *(buff+4) = u01_closed_open_64_53(rand_int2.v[0]);
    *(buff+5) = u01_closed_open_64_53(rand_int2.v[1]);
    *(buff+6) = u01_closed_open_64_53(rand_int2.v[2]);
    *(buff+7) = u01_closed_open_64_53(rand_int2.v[3]);

    *(buff+8)  = u01_closed_open_64_53(rand_int3.v[0]);
    *(buff+9)  = u01_closed_open_64_53(rand_int3.v[1]);
    *(buff+10) = u01_closed_open_64_53(rand_int3.v[2]);
    *(buff+11) = u01_closed_open_64_53(rand_int3.v[3]);


  // Again Fortran will expect a pointer to the value
  //return randval; 
  return 0; 
}



