from ASCsoft.bla import Matrix, matmul_lapack
# import numpy as np 

from time import time

# measure our matrix multiplication time

n = 1

data_et = []
data_lapack = []
while n <= 512: #keep 1024 out for a second 
    n = 2*n

    A = Matrix(n,n)
    B = Matrix(n,n)
    runs =  1+int(min( 1e8 / n**3, 1000))

    ts = time()
    for i in range(runs):
        C = A*B
    te = time()
    t_et = (te-ts)/runs

    ts = time()
    for i in range(runs):
        C = matmul_lapack(A,B)
    te = time()
    t_lapack = (te-ts)/runs

    # print(f"n = {n:4d}   ET={t_et:.6f}s   LAPACK={t_lapack:.6f}s")
    # data_et.append((n, t_et))
    # data_lapack.append((n, t_lapack))

    flops = 2 * (n ** 3)  
    gflops_et = flops / (t_et * 1e9)
    gflops_lapack = flops / (t_lapack * 1e9)

    print(f"n = {n:4d}   "
          f"ET = {t_et:.6f}s  ({gflops_et:.3f} GFLOP/s)   "
          f"LAPACK = {t_lapack:.6f}s  ({gflops_lapack:.3f} GFLOP/s)")

    data_et.append((n, t_et, gflops_et))
    data_lapack.append((n, t_lapack, gflops_lapack))


import matplotlib.pyplot as plt

# sizes_et, times_et = zip(*data_et)
# sizes_la, times_la = zip(*data_lapack)

sizes, times_et, gflops_et = zip(*data_et)
sizes, times_la, gflops_la = zip(*data_lapack)

plt.figure(figsize=(8,6))
plt.plot(sizes, times_et, marker='o', linestyle='-', color='b', label='Execution time (ExpTemplate)')
plt.plot(sizes, times_la, marker='o', linestyle='-', color='r', label='Execution time (LAPACK)')

plt.xscale('log') # , base=2)
plt.yscale('log')

plt.xlabel('Matrix size (N)')
plt.ylabel('Time (seconds)')
plt.title('Execution Time vs. Matrix Size; ExpTemplate vs. LAPACK')
plt.grid(True, which='both', linestyle='--', linewidth=0.5)
plt.legend()
plt.tight_layout()
plt.show()


plt.figure(figsize=(8,6))
plt.plot(sizes, gflops_et, marker='o', linestyle='-', color='b', label='ExpTemplate GFLOP/s')
plt.plot(sizes, gflops_la, marker='o', linestyle='-', color='r', label='LAPACK GFLOP/s')
plt.xscale('log')
plt.xlabel('Matrix size (N)')
plt.ylabel('GFLOP/s')
plt.title('Performance: ExpTemplate vs. LAPACK')
plt.grid(True, which='both', linestyle='--', linewidth=0.5)
plt.legend(); plt.tight_layout()
plt.show()
    
