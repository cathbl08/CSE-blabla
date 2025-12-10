#include <iostream>
#include <sstream>
// include matrix library
#include <matrix.hpp>


#include <taskmanager.hpp>
#include <timer.hpp>
#include <cmath>

using namespace ASC_HPC;
using std::cout, std::endl;

// check namespace for matrix
using namespace ASC_bla;

int main()
{
  timeline = std::make_unique<TimeLine>("demo.trace");

  StartWorkers(3);
  
  RunParallel(10, [] (int i, int size)
  {
    static Timer t("timer one");
    RegionTimer reg(t);
    cout << "I am task " << i << " out of " << size << endl;
  });

  
  RunParallel(6, [] (int i, int s)
  {
    RunParallel(6, [i] (int j, int s2)
    {
      static Timer t("timer nested", {1,0,0});
      RegionTimer reg(t);
      
      std::stringstream str;
      str << "nested, i,j = " << i << "," << j << "\n";
      cout << str.str();
    });
  });



  
  RunParallel(100,  [] (int i, int size)
  {
    static Timer t("one of 100", { 0, 0, 1});
    RegionTimer reg(t);
  });

  {
    static Timer t("100x10 parallel runs", { 0, 0, 1});
    RegionTimer reg(t);

    for (int k = 0; k < 100; k++)
      RunParallel(10,  [] (int i, int size)
      {
        ;
      });
  }

  
  {
    static Timer t("10x10x10 parallel runs", { 0, 0, 1});
    RegionTimer reg(t);

    for (int k = 0; k < 10; k++)
      RunParallel(10, [] (int i, int size)
      {
        RunParallel (10, [] (int j, int size)
        {
        });
      });
  }

  
  RunParallel(1000,  [] (int i, int size)
  {
    static Timer t("timer 1000 tasks", { 1, 0, 0});
    RegionTimer reg(t);
  });

  RunParallel(100, [] (int i, int s)
  {
    static Timer t("timer 4", { 1, 1, 0});
    RegionTimer reg(t);    
    RunParallel(100, [i](size_t j, size_t s2)
    {
      ;
    });
  });

  // use RunParallel for Matrix-Matrix multiplication
  // matrix size
  const size_t N = 500;
  // use 4 tasks (3 workers and 1 main)
  const size_t num_tasks = 4;

  // define matrices A, B and C
  Matrix<double, ColMajor> A(N,N);
  Matrix<double, ColMajor> B(N,N);
  Matrix<double, ColMajor> C(N,N);

  // initialize matrices A and B
  for (size_t i = 0; i < N; i++)
    for (size_t j = 0; j < N; j++)
    {
      A(i,j) = static_cast<double>(i+j);
      B(i,j) = static_cast<double>(i-j);
    }

  // timer for the entire parallel matrix multiplication
  static Timer t_matmul("parallel matrix multiplication", {10,0,0});
  RegionTimer reg_matmul(t_matmul);

  // use RunParallel to perform matrix multiplication C = A * B
  RunParallel(num_tasks, [N, &A, &B, &C](int task_id, int size)
  {
    // determine the range of rows for this task
    size_t first = (N * task_id) / size;
    size_t next = (N * (task_id + 1)) / size;

    // perform multiplication for assigned rows
    /*
    for (size_t i = row_start; i < row_end; i++)
      for (size_t j = 0; j < N; j++)
      {
        double sum = 0.0;
        for (size_t k = 0; k < N; k++)
          sum += A(i,k) * B(k,j);
        C(i,j) = sum;
      }
    */
   C.rows(first, next) = A.rows(first, next) * B;
  });


  
  StopWorkers();
}

