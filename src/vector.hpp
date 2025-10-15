// prevent multiple inclusions of headerfile
#ifndef FILE_VECTOR
#define FILE_VECTOR

#include <iostream>
#include <complex>

template <typename TA, typename TB>
// capture the result type of the addition of two vectors (this allows for addition/subtraction of complex and real vectors)
using TRES = decltype(std::declval<TA>() + std::declval<TB>());
// typedef decltype(std::declval<TA>()+std::declval<TB>()) TRES;
// the above line was mentioned on the course's homepage, but did not work

#include "vecexpr.hpp"


namespace ASC_bla
{
 
  template <typename T, typename TDIST = std::integral_constant<size_t,1> >
  class VectorView : public VecExpr<VectorView<T,TDIST>>
  {
  protected:
    T * m_data;
    size_t m_size;
    TDIST m_dist;
  public:
    VectorView() = default;
    VectorView(const VectorView &) = default;
    
    template <typename TDIST2>
    VectorView (const VectorView<T,TDIST2> & v2)
      : m_data(v2.data()), m_size(v2.Size()), m_dist(v2.dist()) { }
    
    VectorView (size_t size, T * data)
      : m_data(data), m_size(size) { }
    
    VectorView (size_t size, TDIST dist, T * data)
      : m_data(data), m_size(size), m_dist(dist) { }
    
    template <typename TB>
    VectorView & operator= (const VecExpr<TB> & v2)
    {
      assert (m_size == v2.size());
      for (size_t i = 0; i < m_size; i++)
        m_data[m_dist*i] = v2(i);
      return *this;
    }

    VectorView & operator= (T scal)
    {
      for (size_t i = 0; i < m_size; i++)
        m_data[m_dist*i] = scal;
      return *this;
    }

    T * data() const { return m_data; }
    size_t size() const { return m_size; }
    auto dist() const { return m_dist; }
    
    T & operator()(size_t i) { return m_data[m_dist*i]; }
    const T & operator()(size_t i) const { return m_data[m_dist*i]; }
    
    auto range(size_t first, size_t next) const {
      assert(first <= next && next <= m_size);
      return VectorView(next-first, m_dist, m_data+first*m_dist);
    }

    auto slice(size_t first, size_t slice) const {
      return VectorView<T,size_t> (m_size/slice, m_dist*slice, m_data+first*m_dist);
    }
      
  };
  
  

  
  template <typename T>
  class Vector : public VectorView<T>
  {
    typedef VectorView<T> BASE;
    using BASE::m_size;
    using BASE::m_data;
  public:
    Vector (size_t size) 
      : VectorView<T> (size, new T[size]) { ; }
    
    // copy constructor
    Vector (const Vector & v)
      : Vector(v.size())
    {
      *this = v;
    }

    // move constructor
    Vector (Vector && v)
      : VectorView<T> (0, nullptr)
    {
      std::swap(m_size, v.m_size);
      std::swap(m_data, v.m_data);
    }

    template <typename TB>
    Vector (const VecExpr<TB> & v)
      : Vector(v.size())
    {
      *this = v;
    }
    
    ~Vector () { delete [] m_data; }

    using BASE::operator=;
    Vector & operator=(const Vector & v2)
    {
      assert (m_size == v2.m_size);    
      for (size_t i = 0; i < m_size; i++)
        m_data[i] = v2(i);
      return *this;
    }

    // assignment operator (move)
    Vector & operator= (Vector && v2)
    {
      std::swap(m_size, v2.m_size);
      std::swap(m_data, v2.m_data);
      return *this;
    }
    
    size_t Size() const { return size; }
    // access operator
    T & operator()(size_t i) { return data[i]; }
    // access operator (for const objects)
    const T & operator()(size_t i) const { return data[i]; }
  };

  /*
  // vector addition - also handles addition of complex-complex or complex-real vectors
  template <typename TA, typename TB>
  Vector<TRES<TA,TB>> operator+ (const Vector<TA> & a, const Vector<TB> & b)
  {
    Vector<TRES<TA,TB>> sum(a.Size());
      for (size_t i = 0; i < a.Size(); i++)
        sum(i) = a(i)+b(i);
      return sum;
  }

  // vector subtraction - also handles subtraction of complex-complex or complex-real vectors
  template <typename TA, typename TB>
  Vector<TRES<TA,TB>> operator- (const Vector<TA> & a, const Vector<TB> & b)
  {
    Vector<TRES<TA,TB>> diff(a.Size());
    for (size_t i = 0; i < a.Size(); i++)
      diff(i) = a(i)-b(i);
    return diff;
  }

  template <typename T>
  Vector<T> operator- (const Vector<T> & a)
  {
    Vector<T> diff(a.Size());
    for (size_t i = 0; i < a.Size(); i++){
      if (a(i) != 0)
        diff(i) = -a(i);
      else
        diff(i) = a(i);
    }
    return diff;
  }

  template <typename T>
  T operator* (const Vector<T> & a, const Vector<T> & b)
  {
    T mult=0;
    for (size_t i = 0; i < a.Size(); i++) {
      mult += a(i)*b(i);
    }
      

    return mult;
  }
*/  
  template <typename T>
  std::ostream & operator<< (std::ostream & ost, const Vector<T> & v)
  {
    if (v.size() > 0)
      ost << v(0);
    for (size_t i = 1; i < v.size(); i++)
      ost << ", " << v(i);
    return ost;
  }
  
}

#endif
