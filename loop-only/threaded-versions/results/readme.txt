Fully populated processors - using OpenMP threads
=================================================

Comparison producing 640 million random integers (no float/double conversions).

E5-2680v3 Haswell 12-core 2.5GHz 
v
KNL-7210 (64 cores - 4HW threads per core) 1.3GHz

================================================================================

Best results on one 12-core Haswell:
24 threads (2 threads/core): = 0.077792 seconds (3.1 x non-vec code)

Best results on one KNL-7210:
128 threads (2HW threads/core): = 0.015215 seconds (14.0 x non-vec code)

KNL = 5.1x Haswell

Note: 
Threaded KNL results used the source code version
(does not include fast rol instrinsic as of intel 16.0):
