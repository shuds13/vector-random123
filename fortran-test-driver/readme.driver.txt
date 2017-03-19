The driver runs several versions for timing comparison - each produces 10 million random
numbers. The first is the original serial versions. These are followed by

single stream fixed - a single stream of numbers varying only one counter with const vector length.
single stream arg - provides vector length as an argument
ctrkey_all - Provides counter/keys for all Random number sets. Also provides vector length as argument.

The driver is written in Fortran and uses a wrapper to the C routines. 

----------------------------------------------------------------
In  vector123.c set the following for vector length of system.

VECTOR_LENGTH_BYTES
NUM_VALS_32
NUM_VALS_64

Default values will be set

*see readme.txt in source directory for more details.

Also for test harness must set num_vals at top of rng_wrapper.F
This should be same as NUM_VALS_32 or NUM_VALS_64 and is the number of sets
to produce in each all the the vectorised routines.

Note:
It is not to hard to set up to get these values from vector123 - and for many systems
from preproc defines - eg. __avx512f__. See commented out lines at top of vector123.c

The test harness runs several versions for timing comparison - each produces 10 million random
numbers. The first is the original serial versions. These are followed by

single stream fixed - a single stream of numbers varying only one counter with const vector length.
single stream arg - provides vector length as an argument
ctrkey_all - Provides counter/keys for all Random number sets. Also provides vector length as argument.


----------------------------------------------------------------
Building
----------------------------------------------------------------
Select desired makefile from list eg.

>cp makefile.intel makefile

then
> make test


Note: The makefiles use mpi versions of the compilers. Currently MPI is only used to pick up
mpi_wtime. This could easily be removed/changed to different timer if MPI is not available.

