#include <iostream>
#include <random>

#include "vector.hpp"
#include "vecexpr.hpp"

namespace bla = ASC_bla;

int main()
{
  size_t n = 5, m = 4;
  bla::Vector<double> x(n), y(n);
  bla::Vector<std::complex<double>> comple(n);

  // init vectors
  for (size_t i = 0; i < x.Size(); i++)
  {
    x(i) = i;
    y(i) = 10;
    comple(i).real(i);
    comple(i).imag(i * 10);
  }

  bla::Vector<double> z1 = x + y;
  // bla::Vector<double> z2 = x-y;
  // bla::Vector<double> z3 = -x;
  // double z4 = x*y;

  std::cout << "x = " << x << std::endl;
  std::cout << "y = " << y << std::endl;
  std::cout << "x+y = " << z1 << std::endl;
  // std::cout << "x-y = " << z2 << std::endl;
  // std::cout << "-x = " << z3 << std::endl;
  // std::cout << "x*y = " << z4 << std::endl;

  std::cout << std::endl;

  std::cout << "comple = " << comple << std::endl;
  // std::cout << "comple + x = " << comple+x << std::endl;
  // std::cout << "y - comple = " << y-comple << std::endl;
  // std::cout << "comple + comple = " << comple+comple << std::endl;
}
