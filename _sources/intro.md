# Welcome to ASC-bla's documentation!


ASC-bla is a C++ library for basic linear algebra operations.
The library provides template classes **Vector** and **Matrix**.

## Installation

install it via git-clone:

    git clone https://github.com/TUWien-ASC/ASC-bla.git


To configure and build some tests do

    cd ASC-bla
    mkdir build
    cd build
    cmake ..
    make
    

## Using ASC-bla

To use ASC-bla in your code, set the compiler include path properly, and include the header files

    #include <vector.h>
    #include <matrix.h>

All objects are implemented in the namespace ASC_bla. To use them with less typing, you can set

    namespace bla = ASC_bla;

or even

    
    using namespace ASC_bla;

    

You can create vectors and compute with vectors like:

                 
```cpp
Vector<double> x(5), y(5), z(5);
for (int i = 0; i < x.Size(); i++)
   x(i) = i;
y = 5.0
z = x+3*y;
cout << "z = " << z << endl;
```

For matrices you can choose between row-major (`RowMajor`) or column-major (`ColMajor`) storage,
default is row-major.

```cpp
Matrix<double,RowMajor> m1(5,3), m2(3,3);
for (int i = 0; i < m1.Height(); i++)
  for (int j = 0; j < m1.Width(); j++)
    m1(i,j) = i+j;
m2 = 3.7;
Matrix product = m1 * m2;
```

You can extract a row or a column from a matrix:

```cpp
Vector col1 = product.Col(1);
```

some changes ...  

   
