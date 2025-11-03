#include <iostream>

#include <vector.hpp>
#include <matrix.hpp>
#include <lapack_interface.hpp>
#include "matexpr.hpp"  

using namespace ASC_bla;
using namespace std;


int main()
{
  Vector<double> x(5);
  Vector<double> y(5);

  for (int i = 0; i < x.size(); i++)
    {
      x(i) = i;
      y(i) = 2;
    }

  cout << "x = " << x << endl;
  cout << "y = " << y << endl;
  
  addVectorLapack (2, x, y);  
  cout << "y+2*x = " << y << endl;

  Matrix<double, RowMajor> A(2, 3);
  Matrix<double, RowMajor> B(3, 2);
  Matrix<double, ColMajor> C(2, 2);

  A(0,0) = 1; A(0,1) = 2; A(0,2) = 3;
  A(1,0) = 4; A(1,1) = 5; A(1,2) = 6; 

  B(0,0) = 7; B(0,1) = 8;
  B(1,0) = 9; B(1,1) = 10;
  B(2,0) = 11; B(2,1) = 12;

  multMatMatLapack 
    (A,
     B,
     C);  

  cout << "C = A*B = " << endl;
  cout << C << endl;

  // Row-major result C
  Matrix<double, RowMajor> C_rm(2,2);
  multMatMatLapack(A, B, C_rm);

  cout << "C_rm (RowMajor) = A*B = " << endl;
  cout << C_rm << endl;

  double expected[2][2] = {{58, 64}, {139, 154}};
  bool ok = true;
  for (size_t i = 0; i < 2; ++i)
    for (size_t j = 0; j < 2; ++j)
      if (std::abs(C_rm(i,j) - expected[i][j]) > 1e-12) ok = false;

  cout << (ok ? "RowMajor matmul OK " : "RowMajor matmul mismatch ehh )))") << endl;


  // ------- Same matrices, now via syntactic sugar -------
  Matrix<double, RowMajor> C_sugar(2,2);
  MatrixView<double, RowMajor> Av = A, Bv = B, Cv = C_sugar;

  Cv = (Av * Bv) | Lapack;

  cout << "C_sugar (RowMajor) = (A*B)|Lapack =" << endl;
  cout << C_sugar << endl;

  bool ok_sugar = true;
  for (size_t i = 0; i < 2; ++i)
    for (size_t j = 0; j < 2; ++j)
      if (std::abs(C_sugar(i,j) - expected[i][j]) > 1e-12) ok_sugar = false;

  cout << (ok_sugar ? "Sugar matmul OK" : "Sugar matmul mismatch") << endl;

  return 0;
  
}

  
