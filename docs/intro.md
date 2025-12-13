# Welcome to ASC-bla's documentation! 
# Hello from Group CSE-blabla!

ASC-bla is a C++ library for basic linear algebra operations.
The library provides template classes **Vector** and **Matrix**.

The development of the library is performed by the group CSE-blabla.

Group members:
- Andreas Kickenweitz, 11771225
- Catherine Bley, 12002266
- Marcus Konrath, 01300688
- Natalia Tylek, 12332258

## Overview

ASC-bla is a C++ dense linear algebra library providing efficient, template-based classes for basic operations on vectors and matrices.
It combines expression templates, Lapack interface, and pybind11 binding for use from both C++ and Python.

## Installation

install it via git-clone:

    git clone https://github.com/cathbl08/CSE-blabla.git

upstream repo:
    git clone https://github.com/TUWien-ASC/ASC-bla.git


To configure and build some tests do

    cd CSE-blabla
    mkdir build
    cd build
    cmake ..
    make

on Windows systems replace only **make** with

    cmake --build .
    

## Using ASC-bla

To use ASC-bla in your code, set the compiler include path properly, and include the header files

    #include <vector.hpp>
    #include <matrix.hpp>

All objects are implemented in the namespace ASC_bla. To use them with less typing, you can set

    namespace bla = ASC_bla;

or
    
    using namespace ASC_bla;

    

You can create vectors and compute with vectors like:

                 
```cpp
Vector<double> x(5), y(5), z(5);
for (int i = 0; i < x.Size(); i++)
   x(i) = i;
y = 5.0
z = x+3*y;
cout << "z = " << z << endl;
```

For matrices you can choose between row-major (`RowMajor`) or column-major (`ColMajor`) storage,
default is row-major.

```cpp
Matrix<double,RowMajor> m1(5,3), m2(3,3);
for (int i = 0; i < m1.Height(); i++)
  for (int j = 0; j < m1.Width(); j++)
    m1(i,j) = i+j;
m2 = 3.
Matrix product = m1 * m2;
```

You can extract a row or a column from a matrix:

```cpp
auto r2 = M.row(2);   
auto c1 = M.col(1);
```

Lapack interface

Vector BLAS-1 routine example:

```cpp
Vector<double> x(5), y(5);
for (size_t i = 0; i < x.size(); i++) { x(i) = i; y(i) = 2; }
addVectorLapack(2.0, x, y);
``` 
Matrix-Matrix multiplication BLAS-3 routine example:

either by direct call:

```cpp
#include <lapack_interface.hpp>

Matrix<double, RowMajor> A(2,3), B(3,2);
Matrix<double, RowMajor> C(2,2);
multMatMatLapack(A, B, C);   
``` 
or by syntactic sugar using | Lapack with expression templates:

```cpp
#include <matexpr.hpp>

Matrix<double, RowMajor> A(2,3), B(3,2), C(2,2);
MatrixView<double, RowMajor> Av = A, Bv = B, Cv = C;
Cv = (Av * Bv) | Lapack;
``` 



## Use the libary from Python

You can also use ASC-bla from Python. After building the library with cmake, go to the `python` subdirectory and install the python package with

    pip install .

Another, more user friendly way with no need to build from source is to install directly from the git repository:

    pip install git+https://github.com/cathbl08/CSE-blabla.git@main

Then you can use the library from Python:

```python
from ASCsoft.bla import Vector
from ASCsoft.bla import Matrix
from ASCsoft.bla import matmul_lapack
```

with the use of expression templates you can:

```python
x = Vector(5)
y = Vector(5)
for i in range(5):
    x[i] = i
y = 5.0
z = x + 3 * y
print("z =", z)
``` 
both the addition and the scalar multiplication are done efficiently without temporary objects and there is no need for explicit loops in your code

similarly with matrices:

```python
m1 = Matrix(5,3)
m2 = Matrix(3,3)
for i in range(5):
    for j in range(3):
        m1[i,j] = i + j
m2 = 3.7
matmul = m1 * m2
```

you only need to write the loops for initialization, the matrix multiplication is perform with efficient memory access


Another way to compute matrix-matrix products is to use the LAPACK interface:

```python
n = 4
A = Matrix(n, n)
B = Matrix(n, n)
for i in range(n):
    for j in range(n):
        A[i, j] = 0.5*i + 0.25*j
        B[i, j] = 0.2*i - 0.1*j

# Expression-template 
C_et = A * B

# LAPACK 
C_la = matmul_lapack(A, B)
```

Performance benchmark can be found in the `py_demos` folder.


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

## Branch-free sorting??

