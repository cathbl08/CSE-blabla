#include <cstdint>
#include <catch2/catch_test_macros.hpp>
#include <catch2/catch_approx.hpp> 


#include <matrix.hpp>
#include <vector.hpp>
#include <matexpr.hpp>

using namespace ASC_bla;
using Catch::Approx;     


TEST_CASE("Matrix basic size test", "[matrix]") {
    Matrix<double> A(3, 4);
    REQUIRE(A.rows() == 3);
    REQUIRE(A.cols() == 4);
}

TEST_CASE("Matrix access and assignment", "[matrix]") {
    Matrix<double> A(2, 2);
    A(0,0) = 1.0; A(0,1) = 2.0;
    A(1,0) = 3.0; A(1,1) = 4.0;

    REQUIRE(A(0,0) == 1.0);
    REQUIRE(A(1,1) == 4.0);
}

TEST_CASE("Matrix ET multiplication", "[matrix][ET]") {
    // Matrix multiplication using your expression template path

    Matrix<double> A(2,2);
    Matrix<double> B(2,2);

    A(0,0)=1; A(0,1)=2;
    A(1,0)=3; A(1,1)=4;

    B(0,0)=5; B(0,1)=6;
    B(1,0)=7; B(1,1)=8;

    // This should call MultMatExpr through the MatrixView operator=
    Matrix<double> C = A * B;

    REQUIRE(C.rows() == 2);
    REQUIRE(C.cols() == 2);

    REQUIRE(C(0,0) == Approx(19.0));  // 1*5 + 2*7
    REQUIRE(C(0,1) == Approx(22.0));  // 1*6 + 2*8
    REQUIRE(C(1,0) == Approx(43.0));  // 3*5 + 4*7
    REQUIRE(C(1,1) == Approx(50.0));  // 3*6 + 4*8
}

TEST_CASE("Matrix ET addition", "[matrix][add]") {
    Matrix<double> A(2,2);
    Matrix<double> B(2,2);

    A(0,0)=1; A(0,1)=1;
    A(1,0)=1; A(1,1)=1;

    B(0,0)=2; B(0,1)=2;
    B(1,0)=2; B(1,1)=2;

    Matrix<double> C = A + B;

    REQUIRE(C(0,0) == Approx(3.0));
    REQUIRE(C(1,1) == Approx(3.0));
}

TEST_CASE("Matrix scaling ET", "[matrix][scale]") {
    Matrix<double> A(2,2);

    A(0,0)=1; A(0,1)=2;
    A(1,0)=3; A(1,1)=4;

    Matrix<double> C = 2.0 * A;

    REQUIRE(C(0,0) == Approx(2.0));
    REQUIRE(C(0,1) == Approx(4.0));
    REQUIRE(C(1,0) == Approx(6.0));
    REQUIRE(C(1,1) == Approx(8.0));
}
