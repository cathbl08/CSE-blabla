// prevent multiple inclusions of headerfile
#ifndef FILE_MATRIX
#define FILE_MATRIX

#include <iostream>
#include <vector.hpp>

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

    // constructor for converting a Vector object to a Matrix object
    Matrix (const ASC_bla::Vector<T> & a)
      : Matrix(a.Size(), 1)
    {
      for (size_t i = 0; i < a.Size(); i++){
        data[i] = a(i);
      }
    }

    // copy constructor
    Matrix (const Matrix & A)
      : Matrix(A.Rows(), A.Cols())
    {
      *this = A;
    }

    // move constructor
    Matrix (Matrix && A)
      : rows(0), cols(0), data(nullptr)
    {
      std::swap(rows, A.rows);
      std::swap(cols, A.cols);
      std::swap(data, A.data);
    }

    // destructor
    ~Matrix () { delete [] data; }
    
    // assignment operator (copy)
    Matrix & operator=(const Matrix & A2)
    {
      if constexpr (ORD == RowMajor){
        for (size_t i = 0; i < A2.Rows(); i++){
          for (size_t j = 0; j < A2.Cols(); j++){
            data[i*cols + j] = A2(i,j);
          }
        }
      }
      else{
        for (size_t i = 0; i < A2.Rows(); i++){
          for (size_t j = 0; j < A2.Cols(); j++){
            data[j*rows + i] = A2(i,j);
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
        return data[i*cols + j];
      return data[j*rows + i];
    }

    // access operator (for const objects)
    const T & operator()(size_t i, size_t j) const {
      if constexpr (ORD == RowMajor)
        return data[i*cols + j];
      return data[j*rows + i];
    }
  };

  // matrix-matrix addition
  template <typename T, ORDERING ORD>
  Matrix<T, ORD> operator+ (const Matrix<T, ORD> & A, const Matrix<T, ORD> & B)
  {
    Matrix<T, ORD> sum(A.Rows(),A.Cols());
    for (size_t i = 1; i <= A.Rows(); i++)
      for (size_t j = 1; j <= A.Cols(); j++)
        sum(i,j) = A(i,j)+B(i,j);
    return sum;
  }
  
  // matrix-matrix multiplication
  template <typename T, ORDERING ORD>
  Matrix<T, ORD> operator* (const Matrix<T, ORD> & A, const Matrix<T, ORD> & B)
  { 
    if (A.Cols() != B.Rows())
      throw std::invalid_argument("Matrix multiplication is not defined.");
    else{
      Matrix<T, ORD> C(A.Rows(),B.Cols());
      for (size_t i = 0; i < A.Rows(); i++){
        for (size_t j = 0; j < B.Cols(); j++){
          T sum = T{};
          for (size_t k = 0; k < A.Cols(); k++)
            sum += A(i,k)*B(k,j);
          C(i,j) = sum;
        }
      }
      return C;
    }
  }

  // matrix-vector multiplication
  // this function converts the vector to a matrix first, in order to utilize matrix-matrix multiplication function
  template <typename T, ORDERING ORD>
  Matrix<T, ORD> operator* (const Matrix<T, ORD> & A, const ASC_bla::Vector<T> & b)
  {
    Matrix<T, ORD> B(b);
    return A*B;
  }

  template <typename T, ORDERING ORD>
  std::ostream & operator<< (std::ostream & ost, const Matrix<T,ORD> & A)
  {
    for (size_t i = 0; i < A.Rows(); i++){
      for (size_t j = 0; j < A.Cols(); j++){
        ost << A(i, j);
        j == A.Cols()-1 ? ost << std::endl : ost << ", ";
      }
    }
    return ost;
  }
  
}

#endif
