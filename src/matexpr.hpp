#ifndef FILE_EXPRESSION_MATRIX
#define FILE_EXPRESSION_MATRIX

#include <cassert>
#include <iostream>

#include <type_traits>
#include <utility>

#include "vecexpr.hpp"

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
    SumMatExpr (TA _A, TB _B) : A(_A), B(_B) { }
    auto operator() (size_t i, size_t j) const { return A(i,j)+B(i,j); }
    size_t rows() const { return A.rows(); }
    size_t cols() const { return A.cols(); }      
  };
  
  template <typename TA, typename TB>
  auto operator+ (const MatExpr<TA> & A, const MatExpr<TB> & B)
  {
    assert (A.rows() == B.rows()) && (A.cols() == B.cols());
    return SumMatExpr(A.derived(), B.derived());
  }


// ***************** Scaling a matrix *****************
  
  template <typename TSCAL, typename TM>
  class ScaleMatExpr : public MatExpr<ScaleMatExpr<TSCAL,TM>>
  {
    TSCAL scal;
    TM mat;
  public:
    ScaleMatExpr (TSCAL _scal, TM _mat) : scal(_scal), mat(_mat) { }
    auto operator() (size_t i, size_t j) const { return scal*mat(i,j); }
    size_t rows() const { return mat.rows(); }
    size_t cols() const { return mat.cols(); }      
  };
  
  // template <typename S, typename TM>
  // auto operator* (S scal, const MatExpr<TM> & m)
  // {
  //   return ScaleMatExpr<S, TM>(scal, m.derived());
  // }

  template <typename S, typename TM,
          std::enable_if_t<std::is_arithmetic_v<S>, int> = 0>
  auto operator*(S scal, const MatExpr<TM>& m)
  {
    return ScaleMatExpr<S, TM>(scal, m.derived());
  }

  template <typename TM, typename S,
            std::enable_if_t<std::is_arithmetic_v<S>, int> = 0>
  auto operator*(const MatExpr<TM>& m, S scal)
  {
    return ScaleMatExpr<S, TM>(scal, m.derived());
  }


  // ***************** Product of two matrices *****************

  template <typename TA, typename TB>
  class MultMatExpr : public MatExpr<MultMatExpr<TA,TB>>
  {
    
    TA A;
    TB B;

  public:

    const TA& left()  const { return A; }
    const TB& right() const { return B; }

    MultMatExpr(TA _A, TB _B) : A(_A), B(_B){}
    auto operator() (size_t i, size_t j) const
    {
      using elemtypeA = typename std::invoke_result<TA, size_t, size_t>::type;
      using elemtypeB = typename std::invoke_result<TB, size_t, size_t>::type;
      using TSUM = decltype(std::declval<elemtypeA>() * std::declval<elemtypeB>());

      TSUM multsum{};                         
      for (size_t k = 0; k < A.cols(); ++k)
        multsum += A(i,k) * B(k,j);
      return multsum;
    }

  
    size_t rows() const { return A.rows(); }
    size_t cols() const { return B.cols(); } 
  };

  template <typename TA, typename TB>
  auto operator* (const MatExpr<TA>& A, const MatExpr<TB>& B)
  {
    assert( A.cols() == B.rows() );
    return MultMatExpr(A.derived(), B.derived());
  }
    

  // **************** Product of matrix and vector  *****************

  template <typename TM, typename TV>
  class MultMatVecExpr : public VecExpr<MultMatVecExpr<TM,TV>>
  {
    TM A; 
    TV x;

  public:
    MultMatVecExpr(TM _A, TV _x) : A(_A), x(_x) {}
    size_t size() const { return A.rows(); }

    auto operator()(size_t i) const {
      using elemtypeA = typename std::invoke_result<TM, size_t, size_t>::type;
      using elemtypeB = typename std::invoke_result<TV, size_t>::type;
      using TSUM = decltype(std::declval<elemtypeA>() * std::declval<elemtypeB>());

      TSUM sum{};
      for (size_t k = 0; k < A.cols(); ++k)
        sum += A(i,k) * x(k);
      return sum;
    }
  };

  template <typename TA, typename TV>
  auto operator* (const MatExpr<TA>& A, const VecExpr<TV>& x)
  {
    assert( A.cols() == x.size() );
    return MultMatVecExpr(A.derived(), x.derived());
  }

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

  // Note: Matrix stream operator<< is defined in matrix.hpp to avoid include cycles.
  
  struct T_Lapack {};
  inline constexpr T_Lapack Lapack{};

  // The “decorated” expression
  template <class TA, class TB>
  class LapackMultExpr : public MatExpr<LapackMultExpr<TA,TB>> {
    TA A; TB B;
  public:
    LapackMultExpr(TA a, TB b) : A(std::move(a)), B(std::move(b)) {}
    size_t rows() const { return A.rows(); }
    size_t cols() const { return B.cols(); }
    auto operator()(size_t i, size_t j) const { return MultMatExpr<TA,TB>(A,B)(i,j); }
    const TA& left()  const { return A; }
    const TB& right() const { return B; }
  };

  // Generic operator| for ANY MultMatExpr
  template <class TA, class TB>
  inline auto operator|(const MultMatExpr<TA,TB>& mult, T_Lapack) {
    // We need accessors on MultMatExpr:
    //   const TA& left() const;  const TB& right() const;
    return LapackMultExpr<TA,TB>(mult.left(), mult.right());
  }

}
 
#endif
