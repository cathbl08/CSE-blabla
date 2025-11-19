#ifndef FILE_CACHEMAT_HPP
#define FILE_CACHEMAT_HPP

#include "matrix.hpp"
#include <iostream>

namespace ASC_bla
{
    // TODO use ordering from matrix.hpp
    enum ORDERING { ColMajor, RowMajor };

    template <typename T, ORDERING ORD, typename TDIST>
    void addMatMat(MatrixView<> A, MatrixView<> B, MatrixView<> C)
    {
        constexpr size_t BH = 96;
        constexpr size_t BW = 96;
        alignas(64) double memBA[BH * BW];
        for (size_t i1 = 0; i1 < A.rows(); i1 += BH)
            for (size_t j1 = 0; j1 < A.cols(); j1 += BW)
            {
                size_t i2 = min(A.rows(), i1 + BH);
                size_t j2 = min(A.cols(), j1 + BW);

                MatrixView<> Ablock(i2 - i1, j2 - j1, BW, memBA);
                Ablock = A.rows(i1, i2).cols(j1, j2);
                addMatMat2(Ablock, B.rows(j1, j2), C.rows(i1, i2));
            }
    }
    template <size_t H, size_t W>
    void AddMatMatKernel<H, W>(size_t K,
                         const double* A, size_t distA,
                         const double* B, size_t distB,
                         double* C, size_t distC)
    {
        for (size_t k = 0; k < K; ++k)
        {
            for (size_t w = 0; w < 12; ++w)
            {
                C[0 * distC + w] += A[0 * distA + k] * B[k * distB + w];
                C[1 * distC + w] += A[1 * distA + k] * B[k * distB + w];
                C[2 * distC + w] += A[2 * distA + k] * B[k * distB + w];
                C[3 * distC + w] += A[3 * distA + k] * B[k * distB + w];
            }
        }
    }

    void addMatMat2(MatrixView<> A, MatrixView<> B, MatrixView<> C)
    {
        constexpr size_t H = 4;
        constexpr size_t W = 12;

        for (size_t j = 0; j + W <= C.cols(); j += W)
            for (size_t i = 0; i + H <= C.rows(); i += H)
                AddMatMatKernel<H, W>(A.cols(), &A(i, 0), A.dist(),
                                      &B(0, j), B.dist(), &C(i, j), C.dist());
        // leftover rows and cols
    }
}

#endif