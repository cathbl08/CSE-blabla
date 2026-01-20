# Automatic Differentiation

In scientific computing, particularly when using implicit time-stepping methods like Implicit Euler or Crank-Nicolson, we frequently need to solve non-linear systems using Newton's method. This requires computing the Jacobian matrix of the system.

Our library provides two automated approaches:

Numerical Differentiation: Approximating derivatives using finite differences.

Automatic Differentiation (AutoDiff): Computing exact derivatives algorithmically using the chain rule.

You can find the implementation for these features in `autodiff.hpp`.

## Numerical Differentiation

Numerical differentiation approximates the derivative without requiring changes to the function's implementation. Our library employs a central difference quotient for the i-th partial derivative:
\begin{equation}
    \frac{\partial f}{\partial x_i} (x) \approx \frac{f(x + \epsilon e_i) - f(x - \epsilon e_i)}{2\epsilon}
\end{equation}
where $e_i$ is the unit vector in the i-th direction and $\epsilon$ is a small perturbation value.

## Automatic Differentiation

Automatic Differentiation (AD) is a technique where the computer tracks both the value of a function and its derivative simultaneously during evaluation.
Unlike numerical differentiation, AD does not use a perturbation parameter $\epsilon$ and does not suffer from approximation errors. Instead, it applies the rules of calculus to every elementary operation.

The core of our AD implementation is the AutoDiff class template found in `autodiff.hpp`. It represents a variable that carries both a scalar value and a gradient vector.

```cpp
template <size_t N, typename T = double>
class AutoDiff
{
  private:
    T m_val;                 // The function value
    std::array<T, N> m_deriv; // The gradient vector
  // ...
};
```

Where N is the number of independent variables, and T is the data type (usually double).

We have overloaded standard arithmetic operators (+, -, *, /) and mathematical functions (sin, cos, exp, log) to work with AutoDiff types, ensuring the derivatives are propagated correctly according to calculus rules.

## Example Usage

You can explore a full working example in `demos/demo_autodiff.cpp`.

### Basic usage

To calculate gradients, you must first define your independent variables using the Variable helper class. Variable<I>(value) initializes the I-th independent variable, setting its partial derivative at index I to 1.0.

```cpp
#include <autodiff.hpp>
#include <iostream>

using namespace ASC_ode;

int main() {
  double x = 1.0, y = 2.0;

  // 1. Define independent variables x (index 0) and y (index 1)
  // We use N=2 because there are two independent variables.
  AutoDiff<2> adx = Variable<0>(x);
  AutoDiff<2> ady = Variable<1>(y);

  // 2. Define the function: f(x,y) = x * sin(y)
  AutoDiff<2> result = adx * sin(ady);

  // 3. Access results
  std::cout << "Value: " << result.value() << std::endl;       // Output: 0.909297
  std::cout << "df/dx: " << result.deriv()[0] << std::endl;    // Output: 0.909297 (sin(2))
  std::cout << "df/dy: " << result.deriv()[1] << std::endl;    // Output: -0.416147 (x*cos(2))
  
  return 0;
}
```

### Higher-order derivatives

The AutoDiff class is templated on the type T. By nesting AutoDiff types, you can compute higher-order derivatives (e.g., Hessians). This works because the value of the outer AutoDiff becomes an AutoDiff itself, tracking the derivative of the derivative.

```cpp
// Computing 2nd derivatives for f(x) = x^2
// func  = x^2
// func' = 2x
// func''= 2

// Define a variable with 1 dimension, where the underlying type is also AutoDiff<1>
AutoDiff<1, AutoDiff<1>> addx{Variable<0>(2.0)}; 

auto result = addx * addx;

std::cout << "Value: " << result.value().value() << std::endl; // 4
std::cout << "1st Deriv: " << result.value().deriv()[0] << std::endl; // 4
std::cout << "2nd Deriv: " << result.deriv()[0].deriv()[0] << std::endl; // 2
```

### Integration with ODE solvers

To use AutoDiff with our ODE solvers (e.g., for the Newton solver), you do not need to manually compute the Jacobian matrix. Instead, you can write your ODE system using a templated evaluate function.

Inherit from NonlinearFunctionAutoDiff (this adapter class automatically handles the Jacobian logic).

Implement a template function T_evaluate. This allows the function to be called with double (for standard evaluation) or AutoDiff (for Jacobian computation).

### Example: the pendulum

The pendulum equation $\alpha''(t) + \frac{g}{L} \sin(\alpha(t)) = 0$ can be rewritten as a first-order system. 
Here is how to implement it to support AutoDiff:

```cpp
// Inherit from NonlinearFunctionAutoDiff using CRTP
class PendulumAD : public NonlinearFunctionAutoDiff<PendulumAD> 
{
    double m_length;
    double m_gravity;
public:
    PendulumAD(double length, double gravity) 
      : m_length(length), m_gravity(gravity) {}

    size_t dimX() const override { return 2; }
    size_t dimF() const override { return 2; }

    // Templated evaluation function
    // T can be 'double' or 'AutoDiff<2>'
    template <typename T>
    void T_evaluate(VectorView<T> x, VectorView<T> f) const
    {
        f(0) = x(1);
        f(1) = - (m_gravity / m_length) * sin(x(0)); 
    }
};
```

By implementing T_evaluate, the library automatically generates the evaluateDeriv function. It creates AutoDiff variables for the input vector x, calls your function, and extracts the derivatives from the result to fill the Jacobian matrix.

### Legendre Polynomials

The library also includes a utility to compute Legendre polynomials, which serves as a test for the AutoDiff class capabilities. This is demonstrated in `demo_autodiff.cpp`, where polynomials and their derivatives are computed and exported to legendre.txt for plotting.

```cpp
int N = 100;
for (int i = 0; i <= N; ++i)
{
  double x = -1.0 + 2.0 * i / N;
  AutoDiff<1> adx = Variable<0>(x);
  std::vector<AutoDiff<1>> P;
  
  // Compute polynomials up to degree 5 at point x
  LegendrePolynomials(5, adx, P);
  
  // P[k].deriv()[0] contains the derivative of the k-th polynomial
}
```

