add_executable (demo_vector demo_vector.cpp)
add_executable (demo_matrix demo_matrix.cpp)
target_sources (demo_vector PUBLIC ../src/vector.hpp ../src/vecexpr.hpp)
target_sources (demo_matrix PUBLIC ../src/matrix.hpp ../src/matexpr.hpp)
