#Note: No threaded version of the intrinsics.
#Note: Default - 640 million RNs produced - 10x non-threaded version.

#Pure integer transformations
REG_SRC_FILE=v123_standalone_loop_nocnvrt_source_thrd.c
#AVX_SRC_FILE=v123_standalone_loop_nocnvrt_avx512_thrd.c
BINARY_BASE_NAME=v123_loop_nocnvrt_thrd

#Producing floats between 0 and 1 (inclusive/exlusive)
#REG_SRC_FILE=v123_standalone_loop_cnvrt-float_source_thrd.c
##AVX_SRC_FILE=v123_standalone_loop_cnvrt-float_avx512_thrd.c
#BINARY_BASE_NAME=v123_loop_floats_thrd

#Producing doubles between 0 and 1 (inclusive/exlusive)
#REG_SRC_FILE=v123_standalone_loop_cnvrt-dble_source_thrd.c
##AVX_SRC_FILE=v123_standalone_loop_cnvrt-dble_avx512_thrd.c
#BINARY_BASE_NAME=v123_loop_doubles_thrd

#----------------------------------------------------
#Compiler Flags

#02
FLAGS_INTEL_NOVEC="-O2 -xhost -qopenmp -no-vec"
FLAGS_INTEL="-O2 -xhost -qopenmp"
FLAGS_INTEL_KNL_NOVEC="-O2 -qopenmp -xMIC-AVX512 -no-vec"
FLAGS_INTEL_KNL="-O2 -qopenmp -xMIC-AVX512"

#03
# FLAGS_INTEL_NOVEC="-O3 -xhost -qopenmp -no-vec"
# FLAGS_INTEL="-O3 -xhost -qopenmp"
# FLAGS_INTEL_KNL_NOVEC="-O3 -qopenmp -xMIC-AVX512 -no-vec"
# FLAGS_INTEL_KNL="-O3 -qopenmp -xMIC-AVX512"

EXTRAS=""
#EXTRAS="-qopt-report=2"
#EXTRAS="-qopt-report=2 -qopt-report-file="

#------------------------------------------------------------------------------------------
#Build Intel Xeon
icc ${FLAGS_INTEL_NOVEC} ${EXTRAS} -o ${BINARY_BASE_NAME}_novec.x  ${REG_SRC_FILE}
icc ${FLAGS_INTEL}       ${EXTRAS} -o ${BINARY_BASE_NAME}_vec.x    ${REG_SRC_FILE}

#------------------------------------------------------------------------------------------
#Build Intel KNL
icc ${FLAGS_INTEL_KNL_NOVEC}  ${EXTRAS} -o ${BINARY_BASE_NAME}_knl_novec.x       ${REG_SRC_FILE}
icc ${FLAGS_INTEL_KNL}        ${EXTRAS} -o ${BINARY_BASE_NAME}_knl_vec.x         ${REG_SRC_FILE}
#icc ${FLAGS_INTEL_KNL}        ${EXTRAS} -o ${BINARY_BASE_NAME}_knl_vec_avx512.x  ${AVX_SRC_FILE}

#------------------------------------------------------------------------------------------
