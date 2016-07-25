#TRY 02
#iCC -O2 -xhost -o test vector123.c

#TRY 03
#iCC -O3 -xhost -o test vector123.c

#KNL 02 - PROB NO DIFF THAN XHOST
icc -O2 -qopenmp -xMIC-AVX512 -o test-vec.x vector123.c

#KNL 03 - PROB NO DIFF THAN XHOST
#iCC -O3 -xMIC-AVX512 -o test vector123.c
