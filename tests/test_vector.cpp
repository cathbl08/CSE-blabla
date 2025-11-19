#include <cstdint>

#define CATCH_CONFIG_MAIN
#include <catch2/catch_test_macros.hpp>

#include <vector.hpp>
using namespace ASC_bla;


TEST_CASE( "vector-test", "[vector size test]" ) {
  Vector<double> x(5);
  Vector<double> y = x;
  REQUIRE(y.size() == 5);
}

TEST_CASE( "vector-test 2", "[vector add and scale test]" ) {
  Vector<double> x(5);
  x = 3;
  Vector<double> y = x + x + (-2) * x + x;
  REQUIRE(y(0) == 3); 
}

TEST_CASE( "vector-test 3", "[buggy vector add and scale test]" ) {
  Vector<double> x(5);
  x = 3;
  Vector<double> y = x + x + (-2) * x + x;
  REQUIRE(y(0) == 4);  // wrong test !!
}