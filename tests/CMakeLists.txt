include(FetchContent)

FetchContent_Declare(
  Catch2
  GIT_REPOSITORY https://github.com/catchorg/Catch2.git
  GIT_TAG        v3.4.0  # Use the latest stable version
)

FetchContent_MakeAvailable(Catch2)


add_executable (tests_vector tests_vector.cpp tests_matrix.cpp)
target_sources (tests_vector PUBLIC ../src/vector.hpp ../src/vecexpr.hpp ../src/matrix.hpp ../src/matexpr.hpp)

target_link_libraries(tests_vector Catch2::Catch2WithMain)



add_test(NAME catchtest COMMAND tests_vector)

