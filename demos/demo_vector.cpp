#include <iostream>

#include <vector.hpp>
#include <matrix.hpp>

namespace bla = ASC_bla;


int main()
{
  size_t n = 5, m = 4;
  bla::Vector<double> x(n), y(n);
  bla::Matrix<double> A(n,m), B(n,m), C(m,n);

  for (size_t i = 0; i < x.Size(); i++)
    {
      x(i) = i;
      y(i) = 10;
    }
  
  for (size_t i = 1; i <= A.Rows(); i++){
    for (size_t j = 1; j <= A.Cols(); j++){
      A(i,j) = i*10 + j;
    }
  }

  bla::Vector<double> z1 = x+y;
  bla::Vector<double> z2 = x-y;
  bla::Vector<double> z3 = -x;
  double z4 = x*y;
  
  std::cout << "x+y = " << z1 << std::endl;
  std::cout << "x-y = " << z2 << std::endl;
  std::cout << "-x = " << z3 << std::endl;
  std::cout << "x*y = " << z4 << std::endl;

  std::cout << std:: endl;

  std::cout << "---- Matrix A ----" << std::endl;
  std::cout << A << std::endl;
  std::cout << "A(1,2) = " << A(1,2) << std::endl;

  std::cout << std::endl;
  
  B = A;

  std::cout << "---- Matrix B ----" << std::endl;
  std::cout << B << std::endl;
  std::cout << "B(4,3) = " << B(4,3) << std::endl;

  std::cout << std::endl;
}
