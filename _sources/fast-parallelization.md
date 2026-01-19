# Parallelization

For lightweight thread synchronization, we utilize the Compare-and-Swap (CAS) operation.
Unlike standard mutexes that put threads to sleep, this implementation uses `std::atomic<T>::compare_exchange_strong` to create a spinlock. 
This allows a thread to poll a memory location and update it only if the state is exactly as expected. This mechanism is the engine behind our `SimpleLockFreeQueue` (found in `Fast-CSErious/concurrentqueue/benchmarks/simplelockfree.h`), allowing workers to claim tasks with near-zero latency.

We also parallelized the matrix-matrix multiplication by using the RunParallel task manager to distribute work across the available CPU cores. 
In `demo_tasks.cpp`, we demonstrated a row-based partitioning strategy where the result matrix is divided into independent blocks. Each thread calculates a specific range of rows, ensuring thread safety by design since no two workers write to the same memory location.

```cpp
StartWorkers(7); // Main thread + 7 workers
const size_t num_tasks = 8;

RunParallel(num_tasks, [N, &A, &B, &C](int task_id, int size) {
    size_t first = (N * task_id) / size;
    size_t next = (N * (task_id + 1)) / size;
    
    // Each thread computes its unique slice of the result matrix
    for (size_t i = first; i < next; i++)
        for (size_t j = 0; j < N; j++)
            for (size_t k = 0; k < N; k++)
                C(i,j) += A(i,k) * B(k,j);
});
StopWorkers();
```

Performance can be visualized using the Vite Trace Explorer to ensure all CPU cores are equally utilized during parallel execution.

By running the parallel matrix multiplication for a 500 x 500 system, the Vite tracer visualizes how the workload is distributed across 4 total threads.

<img src="images/vite_tracer_4.png" alt="Parallel Vite Trace (4 threads)"  width="800px">

Here, we used a 2000 x 2000 system and 8 threads to run the matrix multiplication.

<img src="images/vite_tracer_8.png" alt="Parallel Vite Trace (8 threads)"  width="800px">


