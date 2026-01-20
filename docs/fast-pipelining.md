# Pipelining

In high-performance computing, a CPU can stall if it must wait for a calculation to finish before starting the next.
To prevent this, CSE-blabla utilizes tile-based accumulation to hide instruction latency and maximize throughput.

Our implementation in `matexpr.hpp`, specifically MultMatExpr::GetTile, optimizes the Arithmetic Intensity of matrix operations.
Instead of a single dependency chain, we process a tile of height H to maintain multiple independent accumulators. While the CPU waits for the result of one row, it simultaneously executes the next, ensuring the arithmetic units are never idle.

To push performance toward the hardware's theoretical limit, we use three key strategies.
We load a single scalar A_ik and broadcast it across all lanes of a SIMD vector.
This scalar is multiplied against an entire SIMD block of Matrix B, performing multiple operations per load.
Organizing data into tiles reduces the ratio of slow memory transfers to fast floating-point operations.

Benchmarks in `Fast-CSErious/demos/simd_timings.cpp` confirm that this pipelined approach outperforms standard loops by saturating the CPU's pipeline. As a user, these optimizations are triggered automatically via expression templates. 

Simply writing a standard matrix product utilizes the optimized kernels:

```cpp
#include <matrix.hpp>
using namespace ASC_bla;

// The Expression Template system automatically breaks this 
// operation into pipelined SIMD tiles:
Matrix<double> C = A * B;
```



