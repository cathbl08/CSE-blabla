#ifndef FILE_VECTOR
#define FILE_VECTOR

#include <iostream>

namespace ASC_bla
{
  
  template <typename T>
  class Vector
  {
    size_t size;
    T * data;
    
  public:
    Vector (size_t _size) 
      : size(_size), data(new T[size]) { ; }
    
    Vector (const Vector & v)
      : Vector(v.Size())
    {
      *this = v;
    }

    Vector (Vector && v)
      : size(0), data(nullptr)
    {
      std::swap(size, v.size);
      std::swap(data, v.data);
    }

    ~Vector () { delete [] data; }
    
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


  template <typename T>
  Vector<T> operator+ (const Vector<T> & a, const Vector<T> & b)
  {
    Vector<T> sum(a.Size());
    for (size_t i = 0; i < a.Size(); i++)
      sum(i) = a(i)+b(i);
    return sum;
  }

  template <typename T>
  Vector<T> operator- (const Vector<T> & a, const Vector<T> & b)
  {
    Vector<T> diff(a.Size());
    for (size_t i = 0; i < a.Size(); i++)
      diff(i) = a(i)-b(i);
    return diff;
  }

  template <typename T>
  Vector<T> operator- (const Vector<T> & a)
  {
    Vector<T> diff(a.Size());
    for (size_t i = 0; i < a.Size(); i++)
      diff(i) = -a(i);
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
