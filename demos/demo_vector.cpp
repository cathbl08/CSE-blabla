#include <iostream>

#include <vector.hpp>

namespace bla = ASC_bla;


int main()
{
  size_t n = 5;
  bla::Vector<double> x(n), y(n);

  for (size_t i = 0; i < x.Size(); i++)
    {
      x(i) = i;
      y(i) = 10;
    }

  bla::Vector<double> z1 = x+y;
  bla::Vector<double> z2 = x-y;
  bla::Vector<double> z3 = -x;
  double z4 = x*y;
  
  std::cout << "x+y = " << z1 << std::endl;
  std::cout << "x-y = " << z2 << std::endl;
  std::cout << "-x = " << z3 << std::endl;
  std::cout << "x*y = " << z4 << std::endl;
}
