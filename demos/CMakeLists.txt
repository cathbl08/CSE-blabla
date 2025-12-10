add_executable (demo_vector demo_vector.cpp)
add_executable (demo_matrix demo_matrix.cpp)

# added
add_executable (demo_tasks demo_tasks.cpp 
                ../Fast-CSErious/src/taskmanager.cpp 
                ../Fast-CSErious/src/timer.cpp)

target_sources (demo_vector PUBLIC ../src/vector.hpp ../src/vecexpr.hpp)
target_sources (demo_matrix PUBLIC ../src/matrix.hpp ../src/matexpr.hpp)

# added
target_link_libraries(demo_tasks PRIVATE /opt/homebrew/opt/libomp/lib/libomp.dylib)
target_compile_options(demo_tasks PRIVATE -Xpreprocessor -fopenmp)
target_include_directories(demo_tasks PRIVATE /opt/homebrew/opt/libomp/include)