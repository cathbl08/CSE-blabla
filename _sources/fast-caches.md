# Caches

Large matrices cannot fit entirely into the CPU's fastest memory, the cache. Without optimization, the CPU spends more time waiting for data from the RAM than performing calculations.
CSE-blabla solves this using a blocked approach to ensure maximum data reuse.

We implement this through a two-level blocking system found in `matrix.hpp` and `matexpr.hpp`.
At the first level, the library partitions large matrices into macro-blocks that fit within the Level 2 cache.
At the second level, these blocks are further divided into tiny tiles processed by our micro-kernels, which are small enough to be stored directly in the Level 1 cache and CPU registers.
This ensures that every piece of data loaded from the main memory is used for as many calculations as possible before being replaced.

Users do not need to manage these memory levels manually, as the library is designed to handle cache-blocking automatically through our Expression Template system. By simply using standard matrix operators, the library selects the most efficient evaluation path based on the specific dimensions and storage order of your data.