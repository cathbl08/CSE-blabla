#ifndef FILE_EXPRESSION_MATRIX
#define FILE_EXPRESSION_MATRIX

#include <cassert>
#include <iostream>

namespace ASC_bla
{
  template <typename T>
  class MatExpr
  {
  public:
    auto derived() const { return static_cast<const T&> (*this); }
    // size_t size() const { return derived().size(); }
    size_t rows() const { return derived().rows(); }
    size_t cols() const { return derived().cols(); }
    auto operator() (size_t i, size_t j) const { return derived()(i,j); }
  };
  
 // ***************** Sum of two matrices *****************

  template <typename TA, typename TB>
  class SumMatExpr : public MatExpr<SumMatExpr<TA,TB>>
  {
    TA A;
    TB B;
  public:
    SumVecExpr (TA _A, TB _B) : A(_A), B(_B) { }
    auto operator() (size_t i, size_t j) const { return A(i,j)+B(i,j); }
    size_t rows() const { return A.rows(); }
    size_t cols() const { return B.cols(); }      
  };
  
  template <typename TA, typename TB>
  auto operator+ (const MatExpr<TA> & A, const MatExpr<TB> & B)
  {
    assert (A.rows() == B.rows()) && (A.cols() == B.cols());
    return SumMatExpr(A.derived(), B.derived());
  }


// ***************** Scaling a matrix *****************
  
  template <typename TSCAL, typename TV>
  class ScaleMatExpr : public MatExpr<ScaleMatExpr<TSCAL,TV>>
  {
    TSCAL scal;
    TV mat;
  public:
    ScaleMatExpr (TSCAL _scal, TV _mat) : scal(_scal), mat(_mat) { }
    auto operator() (size_t i, size_t j) const { return scal*mat(i,j); }
    size_t rows() const { return mat.rows(); }
    size_t cols() const { return mat.cols(); }      
  };
  
  template <typename T>
  auto operator* (double scal, const MatExpr<T> & m)
  {
    return ScaleMatExpr(scal, m.derived());
  }

  template <typename TA, typename TB>
  class MultMatExpr : public MatExpr<MultMatExpr<TA,TB>>
  {
    TA A;
    TB B;

  public:
    MultMatExpr(TA _A, TB _B) : A(_A), B(_B){}
    auto operator() (size_t i, size_t j) const
    {
      using R = decltype (A(0,0) * B(0,0));
      R multsum{};
      for (size_t k = 0; k<A._cols(); k++)
        multsum += A_(i,k) * B_(k,j);
      return multsum;
    
    }

    size_t rows() const { return A_.rows(); }
    size_t cols() const { return B_.cols(); } 
  };

  // **************** get i'th row of matrix *****************
 
  // template <typename TA, typename TB>
  // auto dot (const VecExpr<TA> & a, const VecExpr<TB> & b)
  // {
  //   assert (a.size() == b.size());

  //   using elemtypeA = typename std::invoke_result<TA,size_t>::type;
  //   using elemtypeB = typename std::invoke_result<TB,size_t>::type;
  //   using TSUM = decltype(std::declval<elemtypeA>()*std::declval<elemtypeB>());

  //   TSUM sum = 0;
  //   for (size_t i = 0; i < a.size(); i++)
  //     sum += a(i)*b(i);
  //   return sum;
  // }

  // ***************** Output operator *****************

  // template <typename T>
  // std::ostream & operator<< (std::ostream & ost, const VecExpr<T> & v)
  // {
  //   if (v.size() > 0)
  //     ost << v(0);
  //   for (size_t i = 1; i < v.size(); i++)
  //     ost << ", " << v(i);
  //   return ost;
  // }

  template <typename T, ORDERING ORD>
  std::ostream & operator<< (std::ostream & ost, const Matrix<T,ORD> & A)
  {
    for (size_t i = 0; i < A.rows(); i++){
      for (size_t j = 0; j < A.cols(); j++){
        ost << A(i,j);
        j == A.cols() - 1 ? ost << std::endl : ost << ", ";
      }
    }
    return ost;
  }
  
}
 
#endif
