Single threaded performance
===========================

These tests show that the threefry generator on its own gains very impressive vectorisation
performance - near perfect when taking into account the relative number of ALUs and vector units
on each system. In particular KNL, when using the fast rol intrinsic has 16.3x over non-vectorised
version - 100% vector efficiency. This should be supported from source code in future compiler verssions.

Comparison producing 64 million random numbers.

Vector v scalar for pure threefry4x32 generator:
Ivy Bridge ~ 2.9x
Haswell    ~ 5.2x
KNL        ~ 12.8/16.4x(AVX512 intrinsics)

The conversion to floats or double, however, takes longer than the generator and reduces speed up significantly.
In particular, if converting 32-bit ints to doubles will require a split of the vector.

Vector v scalar with conversion to double:
Ivy Bridge ~ 1.8x
Haswell    ~ 2.3x
KNL        ~ 6.8x


Vector v scalar with conversion to float:
Ivy Bridge ~ 2.1x
Haswell    ~ 3.1x
KNL        ~ 9.4x

------------------------------

Haswell v KNL:

The results below are based on the source code versions.

While KNL has the best speed up from vectorisation - Haswell is 1.6x faster than KNL for single thread performance
on pure threefry (0.078sec to 0.123). Note: The AVX512 intrinsics version on KNL achieves 0.096sec

KNL seems to do better with the conversions, however.
With conversion to floats KNL has parity with Haswell at single thread.  0.179 v 0.183(knl)

Likewise with conversion to doubles - 0.274 v 0.264(knl)

Thus KNLs cycles per byte is signficantly better than Haswell, taking advantage of better vector performance.
