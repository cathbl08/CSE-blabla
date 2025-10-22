#include <iostream>

#include <vector.hpp>
#include <matrix.hpp>

#include "vecexpr.hpp"
#include "matexpr.hpp"
#include <random>

namespace bla = ASC_bla;

int main()
{
  size_t n = 5, m = 4;
  bla::Vector<double> x(n), y(n);
  bla::Vector<std::complex<double>> comple(n);
  bla::Matrix<double> A(n,m), B(n,m), C(n,n), D(n,m);

  // init vectors
   for (size_t i = 0; i < x.Size(); i++){
     x(i) = i;
     y(i) = 10;
     comple(i).real(i);
     comple(i).imag(i*10);
  
  
  // init A
  for (size_t i = 0; i < A.Rows(); i++){
    for (size_t j = 0; j < A.Cols(); j++){
      A(i,j) = i*10 + j;
    }
  }

  // init C
  for (size_t i = 0; i < C.Rows(); i++){
    for (size_t j = 0; j < C.Cols(); j++){
      C(i,j) = 1;
    }
  }

  bla::Vector<double> z1 = x+y;
  bla::Vector<double> z2 = x-y;
  bla::Vector<double> z3 = -x;
  double z4 = x*y;
  
  std::cout << "x = " << x << std::endl;
  std::cout << "y = " << y << std::endl;
  std::cout << "x+y = " << z1 << std::endl;
  std::cout << "x-y = " << z2 << std::endl;
  std::cout << "-x = " << z3 << std::endl;
  std::cout << "x*y = " << z4 << std::endl;

  std::cout << std::endl;

  std::cout << "comple = " << comple << std::endl;
  std::cout << "comple + x = " << comple+x << std::endl;
  std::cout << "y - comple = " << y-comple << std::endl;
  std::cout << "comple + comple = " << comple+comple << std::endl;

  std::cout << std::endl;

  std::cout << "---- Matrix A ----" << std::endl;
  std::cout << A << std::endl;
  std::cout << "A(1,2) = " << A(1,2) << std::endl;

  std::cout << std::endl;

  // std::cout << "---- Matrix A transpose ----" << std::endl;
  // std::cout << Transpose(A) << std::endl;

  std::cout << std::endl;
  
  B = A;

  std::cout << "---- Matrix B ----" << std::endl;
  std::cout << B << std::endl;
  std::cout << "B(4,3) = " << B(4,3) << std::endl;

  std::cout << std::endl;

  D = C*x;

  std::cout << "---- Matrix D=C*x ----" << std::endl;
  std::cout << D << std::endl;
  std::cout << std::endl;

  bla::Matrix<double> E(x);
  std::cout << "---- Matrix E ----" << std::endl;
  std::cout << E << std::endl;

  std::cout << "---- Matrix C*E ----" << std::endl;
  std::cout << C*E << std::endl;
  std::cout << std::endl; 
  
  // std::cout << "---- Concetenated Matrix [A,B] ----" << std::endl;
  // std::cout << (A<B) << std::endl;
  // std::cout << std::endl;
  
  // int i1 = 0, i2 = 1;

  // std::cout << "---- 1st & 2nd rows of [A,B] swapped ----" << std::endl;
  // std::cout << (A<B).swapRows(i1,i2) << std::endl;
  // std::cout << std::endl;

  //------------------------------------------------------------------------//

  int dim = 3;
  std::random_device rd;
  std::mt19937 gen(rd());
  std::uniform_int_distribution<> distrib1(1, 4), distrib2(1, 4);
  bla::Matrix<double> mat(dim,dim), inv(dim,dim);
  
  for (size_t i = 0; i < dim; i++)
    for (size_t j = 0; j < dim; j++)
      mat(i,j) = distrib1(gen)+distrib2(gen);

  // std::cout << "---- inverse Matrix test ----" << std::endl;
  // std::cout << "mat = " << std::endl;
  // std::cout << mat << std::endl;
  // std::cout << std::endl;

  // std::cout << "inverse(mat) = " << std::endl;
  // std::cout << mat.inv() << std::endl;
  // std::cout << std::endl;

  // std::cout << "---- mat*inv ?= I ---- " << std::endl;
  // std::cout << mat*mat.inv() << std::endl;
  // std::cout << std::endl;


}
}
