# Performance

## Overview

This section details the key HPC concepts implemented, including vectorization (SIMD), cache optimization, and pipelining, to ensure efficient linear algebra operations.

## Vectorization (SIMD)

To unleash the full power of modern CPUs, which support Single Instruction, Multiple Data (SIMD) operations, we implemented a custom wrapper for vector data types and intrinsic functions.

Our solution is the `SIMD<T, S>` template class. It's the engine that handles this parallel magic to get Fast-CSErious!

We built the SIMD class to be flexible and fast, using a split-and-conquer strategy.

If you ask for a vector size S, our class automatically breaks it down into smaller pieces until it reaches the native size your CPU supports (for example a 4-wide AVX register). By using masking, we avoid slow single-number loops for the few remaining elements at the end of an array.
This means our class can handle any vector length you throw at it.

This is achieved by splitting the vector into a low part (m_lo) and a high part (m_hi) recursively until the size is 1 or a native SIMD size.

You don't have to worry about Intel AVX or ARM NEON. We define the basic arithmetic operators (+, *, -) once in `simd.hpp`.
Our system automatically passes the job down to the low-level intrinsic functions available for your machine via files like `simd_avx.hpp` or `simd_arm64.hpp`.
   
The specific, highly-optimized intrinsic functions are only contained in the included files:

```cpp
#ifdef __AVX__
#include "simd_avx.hpp" // Contains SIMD<double, 4> implementation for Intel/AMD
#endif

#if defined(__aarch64__) || defined(_M_ARM64)
#include "simd_arm64.hpp" // Contains SIMD implementation for ARM processors
#endif
```

Now, a simple C++ vector addition is automatically translated into a single parallel instruction on your CPU.

```cpp
SIMD<double,4> va(a);
SIMD<double,4> vb(b);

// this simple addition translates directly into a fast, parallel CPU instruction
SIMD<double,4> vc = va + vb;
```

## Vectorizing mathematical functions

Why implement our own sin or cos? Because the standard library's versions are designed for accuracy across a massive range, which makes them slow. For high-performance simulations, we need speed.

We replaced slow standard calls with custom SIMD-friendly functions implemented in `simd.hpp` that compute exp, sin, and cos in an efficient, parallel manner.

We implemented sin and cos using Argument Reduction to shrink the input angle to a small, central interval, and then used Chebyshev Polynomials for fast approximation.

The vectorized exponential was made fast using Argument Reduction based on ln(2) and a quick calculation of 2^q by directly manipulating the exponent bits of the double-precision number.

To use these vectorized functions, simply pass a SIMD vector to them:

```cpp
// Create a SIMD vector with four exponents: -1, 0, 1, and 2
SIMD<double,4> exponents(-1.0, 0.0, 1.0, 2.0);

// Compute e^x for all four values simultaneously
auto result = exponential(exponents);
```
