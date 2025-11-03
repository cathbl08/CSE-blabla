// prevent multiple inclusions of headerfile
#ifndef FILE_MATRIX
#define FILE_MATRIX

#include <string>
#include <algorithm>
#include <iostream>
#include <vector.hpp>
#include "matexpr.hpp"

// we need a matrix from matrix expression that call matrix view for that we need a constructor 
namespace ASC_bla
{
  enum ORDERING { ColMajor, RowMajor };

// starting with MatrixView class
  template <typename T, ORDERING ORD=RowMajor, typename TDIST=std::integral_constant<size_t,1>>
  class MatrixView : public MatExpr<MatrixView<T,ORD, TDIST>>
  {
  protected:
    size_t m_rows, m_cols, m_dist;
    T * m_data;  

    // implement ordering for MatrixView
    constexpr size_t Index(std::size_t i, std::size_t j) const noexcept {
  if constexpr (ORD == RowMajor) return i * m_dist + j;  
  else return i + j * m_dist;   
}

  public:

    size_t rows() const noexcept { return m_rows; }
    size_t cols() const noexcept { return m_cols; }
    size_t dist() const noexcept { return m_dist; }
    T* data() noexcept { return m_data; }
    const T* data() const noexcept { return m_data; }
    
    MatrixView() = default;
    MatrixView(const MatrixView &) = default;

    T& operator()(size_t i, size_t j) {
    return m_data[Index(i,j)];
    }
    const T& operator()(size_t i, size_t j) const {
      return m_data[Index(i,j)];
    }

    // T*       data()       { return m_data; }
    // const T* data() const { return m_data; }


    // matrixview from another matrixview
    template <typename TDIST2, ORDERING ORD2>
    MatrixView (const MatrixView<T, ORD2, TDIST2> & A)
      : m_data(A.data()), m_rows(A.rows()), m_cols(A.cols()), m_dist (ORD == RowMajor ? A.cols() : A.rows()) { }
    
    // view for contiguous data
    MatrixView (size_t rows, size_t cols, T * data)
      : m_data(data), m_rows(rows), m_cols(cols) { 
        m_dist = (ORD == RowMajor) ? cols : rows;
      }
    

    // view for strided data
    MatrixView (size_t rows, size_t cols,  size_t dist, T * data)
      : m_data(data), m_rows(rows), m_cols(cols), m_dist(dist) { 
      
      }

    // crtp 
    // size_t rows() const { return m_rows; }
    // size_t cols() const { return m_cols; }

    // obtain specific row as a VectorView
    auto row(size_t i) const
    {
      if constexpr (ORD == RowMajor)
        return VectorView<T>(m_cols, &m_data[Index(i,0)]);
      else
        // Column major: elements are strided by m_dist (which is m_rows)
        return VectorView<T, size_t>(m_cols, m_dist, &m_data[i]);
    }
    // obtain specific column as a VectorView
    auto col(size_t j) const
    {
      if constexpr (ORD == ColMajor)
        return VectorView<T>(m_rows, &m_data[Index(0,j)]);
      else
        // Row major: elements are strided by m_dist (which is m_cols)
        return VectorView<T, size_t>(m_rows, m_dist, &m_data[j]);
    }

    // obtain a specific range of rows as a MatrixView; row 'last' is not included in the result
    auto rows(size_t first, size_t last) const
    {
      size_t range = last - first;

      if constexpr (ORD == RowMajor)
        return MatrixView<T>(range, m_cols, &m_data[Index(first,0)]);
      else
        // Column major: elements are strided by m_dist (which is m_rows)
        return MatrixView<T, ORD, size_t>(range, m_cols, m_dist, &m_data[first]);
    }

    // obtain a specific range of columns as a MatrixView; column 'last' is not included in the result
    auto cols(size_t first, size_t last) const
    {
      size_t range = last - first;

      if constexpr (ORD == ColMajor)
        return MatrixView<T, ColMajor>(m_rows, range, &m_data[Index(0,first)]);
      else
        // Column major: elements are strided by m_dist (which is m_rows)
        return MatrixView<T, ORD, size_t>(m_rows, range, m_dist, &m_data[first]);
    }

    // assignment operator but from MatExpr
    template <typename TB>
    MatrixView & operator= (const MatExpr<TB> & A)
    {
      assert ((m_rows == A.rows()) && (m_cols == A.cols()));
      if constexpr (ORD == RowMajor){
        for (size_t i = 0; i < m_rows; i++)
          for (size_t j = 0; j < m_cols; j++)
            m_data[i * m_dist + j] = A(i,j);
        
      }
      else{
        for (size_t j = 0; j < m_cols; ++j)
          for (size_t i = 0; i < m_rows; ++i)
            m_data[i + j * m_dist] = A(i,j);
      }
      return *this;
    }

    // matrix transpose
    auto transpose()
    {
      if constexpr (ORD == RowMajor)
        return MatrixView<T, ColMajor>(m_cols, m_rows, m_dist, m_data);
      else
        return MatrixView<T, RowMajor>(m_cols, m_rows, m_dist, m_data);
    }

    template <typename TA, typename TB>
    MatrixView& operator=(const LapackMultExpr<TA,TB>& L)
    {
      static_assert(std::is_same_v<T,double>,
                    "Lapack path currently implemented only for double");
      multMatMatLapack(L.left(), L.right(), *this);
      return *this;
    }


  }; // end class MatrixView

  // matmatmult LAPACK row majow transposition
  template <class T, ORDERING ORD, class TDIST>
  inline auto trans(MatrixView<T,ORD,TDIST> v)
  {
      return v.transpose();
  }

  template <typename T, ORDERING ORD=RowMajor, typename TDIST=std::integral_constant<size_t,1>>
  class Matrix : public MatrixView<T,ORD, TDIST>
  {
    // typedef MatrixView<T,ORD> BASE;
    // using BASE:: 
    using BASE = MatrixView<T,ORD, TDIST>;
    using BASE::m_rows; using BASE::m_cols; using BASE::m_dist; using BASE::m_data;
    using BASE::operator=; // bring in MatrixView assignment from expressions

  public:
    Matrix (size_t rows, size_t cols) : BASE() {
      m_rows = rows;
      m_cols = cols;
      m_dist = (ORD == RowMajor) ? cols : rows;
      m_data = new T[rows*cols];
      for (size_t i = 0; i < rows*cols; i++)
        m_data[i] = T{};
    }

    template <typename TB>
    Matrix(const MatExpr<TB>& A)
      : Matrix(A.rows(), A.cols())  
    {
      BASE::operator=(A);          
    }

    // constructor for creating the identity matrix (to be adapted in case, e.g., the zero matrix is needed, etc.)
    Matrix (size_t _dim)
      : Matrix(_dim,_dim){
       for (size_t i = 0; i < _dim; ++i) (*this)(i,i) = T(1);
    }
 
    // constructor for converting a Vector object to a Matrix object
    Matrix (const ASC_bla::Vector<T> & a)
      : Matrix(a.size(), 1)
    {
      for (size_t i = 0; i < a.size(); i++){
        m_data[i] = a(i);
      }
    }

    // copy constructor
    Matrix (const Matrix & A)
      : Matrix(A.rows(), A.cols())
    {
      *this = A;
    }

    // move constructor changed
    Matrix (Matrix && A) : BASE() {
      m_rows = A.m_rows;
      m_cols = A.m_cols;
      m_dist = A.m_dist;
      m_data = A.m_data;
      A.m_rows = 0; A.m_cols = 0; A.m_dist = 0; A.m_data = nullptr;
    }

    // destructor
    ~Matrix () { delete [] m_data; }


    // row manipulation
    Matrix& swapRows(size_t i1, size_t i2) {
      std::swap(*(this->row(i1).data()), *(this->row(i2).data()));
      return *this;
    }

    // copy-assign
    Matrix& operator=(const Matrix& A2) {
      assert(m_rows == A2.rows() && m_cols == A2.cols());
      for (size_t i = 0; i < m_rows; ++i)
        for (size_t j = 0; j < m_cols; ++j)
          m_data[this->Index(i,j)] = A2(i,j);
      return *this;
    }

    // resize matrix, to be expanded! As of now, basic functionality to remove rows/cols
    void resize(size_t new_rows, size_t new_cols, const std::string& keep = "upper left")
    {
      if (new_rows == m_rows && new_cols == m_cols) {
        return;
      }

      T* new_data = new T[new_rows * new_cols]{};
      size_t new_dist = (ORD == RowMajor) ? new_cols : new_rows;

      size_t rows_to_copy = std::min(m_rows, new_rows);
      size_t cols_to_copy = std::min(m_cols, new_cols);

      if (keep == "upper left")
      {
        for (size_t i = 0; i < rows_to_copy; ++i) {
          for (size_t j = 0; j < cols_to_copy; ++j) {
            if constexpr (ORD == RowMajor)
              new_data[i * new_dist + j] = (*this)(i, j);
            else // ColMajor
              new_data[i + j * new_dist] = (*this)(i, j);
          }
        }
      }
      else if (keep == "lower right")
      {
        size_t old_row_start = m_rows - rows_to_copy;
        size_t old_col_start = m_cols - cols_to_copy;
        size_t new_row_start = new_rows - rows_to_copy;
        size_t new_col_start = new_cols - cols_to_copy;

        for (size_t i = 0; i < rows_to_copy; ++i) {
          for (size_t j = 0; j < cols_to_copy; ++j) {
            if constexpr (ORD == RowMajor)
              new_data[(new_row_start + i) * new_dist + (new_col_start + j)] = (*this)(old_row_start + i, old_col_start + j);
            else // ColMajor
              new_data[(new_row_start + i) + (new_col_start + j) * new_dist] = (*this)(old_row_start + i, old_col_start + j);
          }
        }
      }

      delete[] m_data;
      m_data = new_data;
      m_rows = new_rows;
      m_cols = new_cols;
      m_dist = new_dist;
    }

    // matrix inverse using Gauss-Jordan: [A,I] -> [I, A^-1]
    Matrix inv() const {
      if (m_rows != m_cols){
        throw std::invalid_argument("Matrix inversion not defined. Matrix must be quadratic.");
      }
      else
      {
        Matrix<T, ORD> inverse(m_cols), tab(m_rows, 2*m_cols); // initializes as identity matrix
        tab = *this<inverse;
        for (size_t i = 0; i < m_rows; i++){
          int pivot = i;
          while (pivot < m_rows && std::fabs(tab(pivot,i)) < 1e-9){
            pivot++;
          }
          if (pivot == m_rows){
            throw std::invalid_argument("Matrix is singular.");
          }
          
          tab.swapRows(i,pivot);
          T pivot_val = tab(i,i);
          
          // normalize pivot, row transformation
          for (int j = i; j < tab.cols(); j++) {
            tab(i,j) /= pivot_val;
          }
          
          // zeros below pivot
          for (size_t k = i+1; k < m_rows; k++){
            T factor = tab(k,i);
            for (size_t j = i; j < tab.cols(); j++){
              tab(k,j) -= factor * tab(i,j);
            }
          }
        }
        
        // zeros above pivot
        for (int i = m_rows-1; i >= 0; i--){
          for (int k = i - 1; k >= 0; k--){
            T factor = tab(k,i);
            for (size_t j = i; j < tab.cols(); ++j){
              tab(k,j) -= factor * tab(i,j);
            }
          }
          
        }
        // std::cout << tab << std::endl; // DOES NOTHING
        Matrix<T, ORD> res(m_cols,m_cols);
        for (size_t i = 0; i < m_rows; i++){
          for (size_t j = 0; j < m_rows; j++){
            res(i,j) = tab(i,j+m_cols);
            }
          }
        return res;
      }
    }

    // helper for QR decompositon
    static int signum(double x) {
        return (0.0 < x) - (x < 0.0);
    }

    std::tuple<Matrix<T, ORD>, Matrix<T, ORD>> qr_decomp(){
    // Matrix qr_decomp()
    // {
      
      int m = this->rows();
      int n = this->cols();
        
      Matrix<T, ORD> Q(m,m), Ak(m, n), R(m, n);
      Ak = *this;

      for (size_t k = 0; k < std::min(m,n); k++){

        // std::cout << "-----------" << std::endl;
        // std::cout << "iteration " << k << std::endl;
        
        // if (k == 0){
        Matrix<T, ORD> Hk(Ak.rows(), Ak.rows()), Ik(Ak.rows()), vk_tp(1, Ak.rows());
        Vector<T> tmp_vk = Ak.col(0);
        double norm_val = tmp_vk.norm();
        
        // std::cout << "xk = " << tmp_vk << std::endl;
        
        tmp_vk = tmp_vk + signum(tmp_vk(0))*norm_val*Ik.col(0);
        norm_val = pow(tmp_vk.norm(),2); // new value

        // std::cout << "vk = " << tmp_vk << std::endl;
        
        auto factor = -2.0/norm_val;

        Matrix<T, ORD> vk(tmp_vk); // Matrix(Vector)
        vk_tp = vk.transpose();

        // std::cout << "vk = " << vk << std::endl;

        Hk = Ik + (factor * (vk * vk_tp));
        // std::cout << "Hk = " << std::endl;
        // std::cout << Hk << std::endl;
        // std::cout << std::endl;
        // std::cout << "factor = " << std::endl;
        // std::cout << factor << std::endl;

        Ak = Hk*Ak;
        std::cout << "Ak = " << Ak << std::endl;

        for (size_t j = 0; j < Ak.cols(); j++){
          R(k, k+j) = Ak(0, j);
        }
        // For the next iteration, work on the lower-right submatrix
        Ak.resize(Ak.rows()-1, Ak.cols()-1, "lower right");
        // std::cout << "---------------" << std::endl;
        // std::cout << "Ak = " << Ak << std::endl;
      }
      Q = *this*R.inv();
      return {Q, R};
    }

    // assignment operator (move)
    // Matrix & operator= (Matrix && A2)
    // {
    //   std::swap(m_rows, A2.rows);
    //   std::swap(m_cols, A2.cols);
    //   std::swap(m_data, A2.data);
    //   return *this;
    // }

    Matrix & operator = (Matrix&& A2)
    {
      if (this == &A2) return *this; // self-assignment check
      delete [] m_data; // free existing resource
      m_rows = A2.m_rows;
      m_cols = A2.m_cols;
      m_dist = A2.m_dist;
      m_data = A2.m_data;
      A2.m_rows = 0; A2.m_cols = 0; A2.m_dist = 0; A2.m_data = nullptr;
      return *this;
    }
  };

  // matrix-matrix addition
  template <typename T, ORDERING ORD>
  Matrix<T, ORD> operator+ (const Matrix<T, ORD> & A, const Matrix<T, ORD> & B)
  {
    Matrix<T, ORD> sum(A.rows(),A.cols());
    for (size_t i = 0; i < A.rows(); i++)
      for (size_t j = 0; j < A.cols(); j++)
        sum(i,j) = A(i,j)+B(i,j);
    return sum;
  }
  
  // matrix-matrix multiplication
  template <typename T, ORDERING ORD>
  Matrix<T, ORD> operator* (const Matrix<T, ORD> & A, const Matrix<T, ORD> & B)
  { 
    if (A.cols() != B.rows())
      throw std::invalid_argument("Matrix multiplication is not defined.");
    else{
      Matrix<T, ORD> C(A.rows(),B.cols());
      for (size_t i = 0; i < A.rows(); i++){
        for (size_t j = 0; j < B.cols(); j++){
          T sum = T{};
          for (size_t k = 0; k < A.cols(); k++)
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

  template <typename T, ORDERING ORD, typename TDIST>
  std::ostream & operator<< (std::ostream & ost, const MatrixView<T,ORD,TDIST> & A)
  {
    for (size_t i = 0; i < A.rows(); i++){
      for (size_t j = 0; j < A.cols(); j++){
        ost << A(i, j);
        j == A.cols()-1 ? ost << std::endl : ost << ", ";
      }
    }
    return ost;
  }

  // sideway concatenation
  template <typename T, ORDERING ORD>
  Matrix<T, ORD> operator< (const Matrix<T, ORD> & A, const Matrix<T, ORD> & B)
  {
    Matrix<T, ORD> C(A.rows(), A.cols() + B.cols());
    for (size_t i = 0; i < C.rows(); i++){
      for (size_t j = 0; j < A.cols(); j++){
        C(i,j) = A(i,j);
      }
      for (size_t j = 0; j < B.cols(); j++){
        C(i,j + A.cols()) = B(i,j);    
      }
    }
    return C;
  }

}

#endif
