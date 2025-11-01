#include <iostream>

#include <vector.hpp>
#include <matrix.hpp>
#include <lapack_interface.hpp>


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
}

  
