# ODE-to-CSE

## Introduction

ODE-to-CSE is a C++ library specifically designed for solving nonlinear ODEs.
The library functions by defining the equation through a right-hand side function object. It offers a variety of time-stepping algorithms that can be utilized to simulate dynamic systems, such as multi-body mass-spring systems.

You can find the underlying theory for these methods here: [Solving ODEs Theory](https://jschoeberl.github.io/IntroSC/ODEs/ODEs.html)

### How to Build and Run the Library

The library depends on pybind11 and lapack, which must be installed on your system. To get started, clone the repository and initialize the nanoblas linear algebra submodule:

```console
git clone https://github.com/my_fork_on_github/ASC-ODE.git
cd ASC-ODE
git submodule update --init
```

To build the library, navigate to the project root and execute the standard CMake workflow:

```console
mkdir build
cd build
cmake ..
make
```

After the build is complete, you can run the demonstration executable from the build directory:

```console
./test_ode
```

This generates simulation results in a text file (e.g., `output_test_ode.txt`).
To visualize these results and generate plots for the mass-spring system, navigate to the demos folder and use the Python script:

```console
python3 plotmassspring.py
```

### Changing the simulation method

You can modify the simulation behavior in `demos/test_ode.cpp`.
The library supports multiple methods including Explicit Euler, Improved Euler, and Implicit Euler.

To switch methods or adjust the resolution, modify the main() function:

```cpp
// From demos/test_ode.cpp
int steps = 200; // Change steps for higher resolution
...
// Switch the stepper method
// ExplicitEuler stepper(rhs);
ImprovedEuler stepper(rhs);
```

Rebuild the project using `make` after any changes to `test_ode.cpp` for the changes to take effect.
