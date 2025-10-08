// prevent multiple inclusions of headerfile
#ifndef FILE_MATRIX
#define FILE_MATRIX

#include <iostream>
//#include <vector.hpp>

namespace ASC_bla
{
  
  enum ORDERING { ColMajor, RowMajor };
  template <typename T, ORDERING ORD=RowMajor>
  class Matrix {
  private:
    size_t rows;
    size_t cols;
    T* data;
    
  public:
    // constructor
    Matrix (size_t _row, size_t _col) 
      : rows(_row), cols(_col), data(new T[_row * _col]) { ; }

    // copy constructor
    Matrix (const Matrix & A)
      : Matrix(A.Rows(), A.Cols())
    {
      *this = A;
    }

    // move constructor
    Matrix (Matrix && A)
      : size(0), data(nullptr)
    {
      std::swap(rows, A.Rows());
      std::swap(cols, A.Cols());
      std::swap(data, A.data);
    }

    // destructor
    ~Matrix () { delete [] data; }
    
    // assignment operator (copy)
    Matrix & operator=(const Matrix & A2)
    {
      if constexpr (ORD == RowMajor){
        for (size_t i = 1; i <= A2.Rows(); i++){
          for (size_t j = 1; j <= A2.Cols(); j++){
            data[(i-1)*cols + (j-1)] = A2(i,j);
          }
        }
      }
      else{
        for (size_t i = 1; i <= A2.Rows(); i++){
          for (size_t j = 1; j <= A2.Cols(); j++){
            data[(j-1)*rows + (i-1)] = A2(i,j);
          }
        }
      }

      return *this;
    }

    // assignment operator (move)
    Matrix & operator= (Matrix && A2)
    {
      std::swap(rows, A2.rows);
      std::swap(cols, A2.cols);
      std::swap(data, A2.data);
      return *this;
    }
    
    size_t Rows() const { return rows; }
    size_t Cols() const { return cols; }

    // access operator
    T & operator()(size_t i, size_t j) {
      if constexpr (ORD == RowMajor)
        return data[(i-1)*cols + (j-1)];
      return data[(j-1)*rows + (i-1)];
    }

    // access operator (for const objects)
    const T & operator()(size_t i, size_t j) const {
      if constexpr (ORD == RowMajor)
        return data[(i-1)*cols + (j-1)];
      return data[(j-1)*rows + (i-1)];
    }
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

  template <typename T, ORDERING ORD>
  std::ostream & operator<< (std::ostream & ost, const Matrix<T,ORD> & A)
  {
    for (size_t i = 1; i <= A.Rows(); i++){
      for (size_t j = 1; j <= A.Cols(); j++){
        ost << A(i, j);
        j == A.Cols() ? ost << std::endl : ost << ", ";
      }
    }
    return ost;
  }
  
}

#endif
