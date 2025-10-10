#include <iostream>

#include <vector.hpp>
#include <matrix.hpp>

namespace bla = ASC_bla;

int main()
{
  size_t n = 5, m = 4;
  bla::Vector<double> x(n), y(n);
  bla::Vector<std::complex<double>> compl(n);
  bla::Matrix<double> A(n,m), B(n,m), C(n,n), D(n,m);

  // init vectors
   for (size_t i = 0; i < x.Size(); i++){
     x(i) = i;
     y(i) = 10;
     compl(i).real(i);
     compl(i).imag(i*10);
  }
  
  // init A
  for (size_t i = 0; i < A.Rows(); i++){
    for (size_t j = 0; j < A.Cols(); j++){
      A(i,j) = i*10 + j;
    }
  }

  // init C
  for (size_t i = 0; i < C.Rows(); i++){
    for (size_t j = 0; j < C.Cols(); j++){
      C(i,j) = 1;
    }
  }

  bla::Vector<double> z1 = x+y;
  bla::Vector<double> z2 = x-y;
  bla::Vector<double> z3 = -x;
  double z4 = x*y;
  
  std::cout << "x = " << x << std::endl;
  std::cout << "y = " << y << std::endl;
  std::cout << "x+y = " << z1 << std::endl;
  std::cout << "x-y = " << z2 << std::endl;
  std::cout << "-x = " << z3 << std::endl;
  std::cout << "x*y = " << z4 << std::endl;

  std::cout << std::endl;

  std::cout << "compl = " << compl << std::endl;
  std::cout << "compl + x = " << compl+x << std::endl;
  std::cout << "y - compl = " << y-compl << std::endl;
  std::cout << "compl + compl = " << compl+compl << std::endl;

  std::cout << std::endl;

  std::cout << "---- Matrix A ----" << std::endl;
  std::cout << A << std::endl;
  std::cout << "A(1,2) = " << A(1,2) << std::endl;

  std::cout << std::endl;

  std::cout << "---- Matrix A transpose ----" << std::endl;
  std::cout << Transpose(A) << std::endl;

  std::cout << std::endl;
  
  B = A;

  std::cout << "---- Matrix B ----" << std::endl;
  std::cout << B << std::endl;
  std::cout << "B(4,3) = " << B(4,3) << std::endl;

  std::cout << std::endl;

  D = C*x;

  std::cout << "---- Matrix D=C*x ----" << std::endl;
  std::cout << D << std::endl;
  std::cout << std::endl;

  bla::Matrix<double> E(x);
  std::cout << "---- Matrix E ----" << std::endl;
  std::cout << E << std::endl;

  std::cout << "---- Matrix C*E ----" << std::endl;
  std::cout << C*E << std::endl;
  std::cout << std::endl;  

}
