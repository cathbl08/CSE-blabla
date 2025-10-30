#include <iostream>
#include <random>

#include "matrix.hpp"
#include "matexpr.hpp"

namespace bla = ASC_bla;

int main()
{
    size_t n = 5, m = 4;
    bla::Vector<double> x(n), y(n);
    bla::Matrix<double> A(n, m), B(n, m), C(n, n), D(n, m);

    // init vector
    for (size_t i = 0; i < x.size(); i++)
        x(i) = i;

    // init A
    for (size_t i = 0; i < A.rows(); i++)
    {
        for (size_t j = 0; j < A.cols(); j++)
            A(i, j) = i * 10 + j;
    }
    // init C
    for (size_t i = 0; i < C.rows(); i++)
    {
        for (size_t j = 0; j < C.cols(); j++)
            C(i, j) = 1;
    }

    std::cout << "---- Matrix A ----" << std::endl;
    std::cout << A << std::endl;
    std::cout << "A(1,2) = " << A(1, 2) << std::endl;
    std::cout << std::endl;

    std::cout << "---- A.row(2) ----" << std::endl;
    std::cout << A.row(2) << std::endl;
    std::cout << "---- A.col(1) ----" << std::endl;
    std::cout << A.col(1) << std::endl;
    std::cout << std::endl;

    std::cout << "---- A.rows(2,4) ----" << std::endl;
    std::cout << A.rows(2,4) << std::endl;
    std::cout << "---- A.cols(1,3) ----" << std::endl;
    std::cout << A.cols(1,3) << std::endl;
    std::cout << std::endl;

    bla::Matrix<double, ASC_bla::ColMajor> AA(A);
    std::cout << "---- Matrix AA (Column Major from A) ----" << std::endl;
    std::cout << AA << std::endl;

    std::cout << "---- AA.row(2) ----" << std::endl;
    std::cout << AA.row(2) << std::endl;
    std::cout << "---- AA.col(1) ----" << std::endl;
    std::cout << AA.col(1) << std::endl;
    std::cout << std::endl;

    std::cout << "---- AA.rows(2,4) ----" << std::endl;
    std::cout << AA.rows(2,4) << std::endl;
    std::cout << "---- AA.cols(1,3) ----" << std::endl;
    std::cout << AA.cols(1,3) << std::endl;
    std::cout << std::endl;


    std::cout << "---- Matrix A transpose ----" << std::endl;
    std::cout << A.transpose() << std::endl;
    std::cout << std::endl;

    std::cout << "---- Matrix AA transpose ----" << std::endl;
    std::cout << AA.transpose() << std::endl;
    std::cout << std::endl;

    B = A;

    std::cout << "---- Matrix B ----" << std::endl;
    std::cout << B << std::endl;
    std::cout << "B(4,3) = " << B(4, 3) << std::endl;
    std::cout << std::endl;

    D = C * x;

    std::cout << "---- Matrix D=C*x ----" << std::endl;
    std::cout << D << std::endl;
    std::cout << std::endl;

    bla::Matrix<double> E(x);
    std::cout << "---- Matrix E ----" << std::endl;
    std::cout << E << std::endl;

    std::cout << "---- Matrix C*E ----" << std::endl;
    std::cout << C * E << std::endl;
    std::cout << std::endl;

    std::cout << "---- Concetenated Matrix [A,B] ----" << std::endl;
    std::cout << (A < B) << std::endl;
    std::cout << std::endl;

    int i1 = 0, i2 = 1;

    std::cout << "---- 1st & 2nd rows of [A,B] swapped ----" << std::endl;
    std::cout << (A < B).swapRows(i1, i2) << std::endl;
    std::cout << std::endl;

    //------------------------------------------------------------------------//

    int dim = 3;
    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_int_distribution<> distrib1(1, 4), distrib2(1, 4);
    bla::Matrix<double> mat(dim, dim), inv(dim, dim);

    for (size_t i = 0; i < dim; i++)
        for (size_t j = 0; j < dim; j++)
            mat(i, j) = distrib1(gen) + distrib2(gen);

    std::cout << "---- inverse Matrix test ----" << std::endl;
    std::cout << "mat = " << std::endl;
    std::cout << mat << std::endl;
    std::cout << std::endl;

    std::cout << "inverse(mat) = " << std::endl;
    std::cout << mat.inv() << std::endl;
    std::cout << std::endl;

    std::cout << "---- mat*inv ?= I ---- " << std::endl;
    std::cout << mat * mat.inv() << std::endl;
    std::cout << std::endl;
}
