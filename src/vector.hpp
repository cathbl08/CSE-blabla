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

namespace ASC_bla
{
  template <typename T>
  class Vector
  {
    size_t size;
    T * data;
    
  public:
    // constructor
    Vector (size_t _size) 
      : size(_size), data(new T[size]) { ; }
    
    // copy constructor
    Vector (const Vector & v)
      : Vector(v.Size())
    {
      *this = v;
    }

    // move constructor
    Vector (Vector && v)
      : size(0), data(nullptr)
    {
      std::swap(size, v.size);
      std::swap(data, v.data);
    }

    // destructor
    ~Vector () { delete [] data; }

    // assignment operator (copy)
    Vector & operator=(const Vector & v2)
    {
      for (size_t i = 0; i < size; i++)
        data[i] = v2(i);
      return *this;
    }

    // assignment operator (move)
    Vector & operator= (Vector && v2)
    {
      std::swap(size, v2.size);
      std::swap(data, v2.data);
      return *this;
    }
    
    size_t Size() const { return size; }
    // access operator
    T & operator()(size_t i) { return data[i]; }
    // access operator (for const objects)
    const T & operator()(size_t i) const { return data[i]; }
  };

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
