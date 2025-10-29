#include <iostream>

#include <vector.hpp>
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
}

  
