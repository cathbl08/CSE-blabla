#ifndef FILE_MATRIX
#define FILE_MATRIX

#include <iostream>
//#include <vector.hpp>

namespace ASC_bla
{
  
  enum ORDERING { ColMajor, RowMajor };
  template <typename T, ORDERING ORD=RowMajor>
  class Matrix {
    size_t rows;
    size_t cols;
    T* data;
    
  public:
    Matrix (size_t _row, size_t _col) 
      : rows(_row), cols(_col), data(new T[_row * _col]) { ; }

    Matrix (const Matrix & a)
      : Matrix(a.Size())
    {
      *this = a;
    }

    Matrix (Matrix && a)
      : size(0), data(nullptr)
    {
      std::swap(rows, a.rows);
      std::swap(cols, a.cols);
      std::swap(data, a.data);
    }

    ~Matrix () { delete [] data; }
    
    Vector & operator=(const Vector & v2)
    {
      for (size_t i = 0; i < size; i++)
        data[i] = v2(i);
      return *this;
    }

    Vector & operator= (Vector && v2)
    {
      std::swap(size, v2.size);
      std::swap(data, v2.data);
      return *this;
    }
    
    size_t Size() const { return size; }
    T & operator()(size_t i) { return data[i]; }
    const T & operator()(size_t i) const { return data[i]; }
  };


  template <typename T, ORDERING ORD>
  Matrix<T, ORD> operator+ (const Matrix<T, ORD> & a, const Matrix<T, ORD> & b)
  {
    Matrix<T, ORD> sum(a.Row(),a.Col());
    for (size_t i = 0; i < a.Row(); i++)
      for (size_t j = 0; j < a.Col(); j++)
        sum(i,j) = a(i,j)+b(i,j);
    return sum;
  }
  
  template <typename T, ORDERING ORD>
  Matrix<T, ORD> operator* (const Matrix<T, ORD> & a, const Matrix<T, ORD> & b)
  {
    if constexpr (ORD == ColMajor)
      {
       
      }
    else
      {
        
      }
  }

  template <typename T>
  std::ostream & operator<< (std::ostream & ost, const Vector<T> & v)
  {
    if (v.Size() > 0)
      ost << v(0);
    for (size_t i = 1; i < v.Size(); i++)
      ost << ", " << v(i);
    return ost;
  }
  
}

#endif
