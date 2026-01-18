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

## Improved Euler Method

The Improved Euler method is an explicit second-order Runge-Kutta scheme. It uses a predictor-corrector approach to achieve better numerical accuracy than the standard Euler method by interpreting the difference quotient as an approximation to the derivative at the mid-point of the interval.

### Mathematical Overview
The method interprets the difference quotient as an approximation to the derivative at the mid-point of the interval.
Given an autonomous ODE $\dot{y} = f(y)$, the solution is updated at each time step using a two-step predictor-corrector process:

An intermediate state $\tilde{y}$ is estimated at the midpoint of the interval:
$$ \tilde{y} = y_n + \frac{\tau}{2} f(y_n) $$

The final update uses the derivative evaluated at this intermediate point:

$$ y_{n+1} = y_n + \tau f(\tilde{y}) $$


### Implementation
In our library, the ExplicitEuler and ImprovedEuler class are defined in `timestepper.hpp` and inherit from the base TimeStepper class.

The constructor requires a shared pointer to a NonlinearFunction which represents the right-hand side $f(y)$ of the ODE:
```cpp
ExplicitEuler(std::shared_ptr<NonlinearFunction> rhs);
ImprovedEuler(std::shared_ptr<NonlinearFunction> rhs);
```
### Parameters
$\tau$ : Time step size
$y$ : The current solution vector $y_n$, modified in-place to contain $y_{n+1}$ after the step.

The core logic calculates the intermediate $y\_tilde$ before performing the final update:

```cpp
void doStep(double tau, VectorView<double> y) override
{
  // Evaluate f at current y
  this->m_rhs->evaluate(y, m_vecf);
  
  // Compute intermediate predictor point y_tilde
  Vector<> y_tilde = y + tau / 2 * m_vecf;
  
  // Evaluate f at y_tilde and perform final update
  this->m_rhs->evaluate(y_tilde, m_vecf);
  y += tau * m_vecf;
}
```

### How can you use this implementation?

As demonstrated in `test_ode.cpp`, you can initialize the stepper with a right-hand side object (such as a MassSpring system) and iterate through the desired number of steps:

```cpp
auto rhs = std::make_shared<MassSpring>(1.0, 1.0);
ImprovedEuler stepper2(rhs);

for (int i = 0; i < steps; i++)
{
  stepper2.doStep(tau, y);
  outfile2 << (i+1) * tau << "  " << y(0) << " " << y(1) << std::endl;
}
```

## Implicit Euler Method

The Implicit Euler method, also known as the backward Euler method, is an implicit first-order time-stepping scheme. Unlike explicit methods, the new step $y_{n+1}$ appears on both sides of the equation.

### Mathematical Overview

For a constant mesh size $\tau$, the numerical update is defined by the following equation:
\begin{equation}
    y_{n+1} = y_n + \tau f(y_{n+1})
\end{equation}

To determine the value of $y_{n+1}$, we must solve a non-linear system of equations at every time step:
\begin{equation}
    G(y_{n+1}) = y_{n+1} - y_n - \tau f(y_{n+1}) = 0
\end{equation}

In this library, the root of $G(y_{n+1})$ is found using a Newton solver.

### Implementation

The ImplicitEuler class is defined in `timestepper.hpp` and inherits from TimeStepper.
The constructor initializes the internal equation $G(y)$ by combining several NonlinearFunction objects:

```cpp
ImplicitEuler(std::shared_ptr<NonlinearFunction> rhs) 
: TimeStepper(rhs), m_tau(std::make_shared<Parameter>(0.0)) 
{
  m_yold = std::make_shared<ConstantFunction>(rhs->dimX());
  auto ynew = std::make_shared<IdentityFunction>(rhs->dimX());
  // Equation: G(y) = y_new - y_old - tau * f(y_new)
  m_equ = ynew - m_yold - m_tau * m_rhs;
}
```

### Parameters

$\tau$ : Time step size
$y$ : The current solution vector $y_n$, modified in-place by the Newton solver to become $y_{n+1}$.

The `doStep` method sets the previous state and time step size before invoking the Newton solver to find $y_{n+1}$:

```cpp
void doStep(double tau, VectorView<double> y) override
{
  m_yold->set(y);      // Update y_n for the current step
  m_tau->set(tau);     // Update current time step size
  NewtonSolver(m_equ, y); // Solve G(y_{n+1}) = 0
}
```

### How can you use this implementation?

As shown in `test_ode.cpp`, you can initialize the stepper with a right-hand side object. Note that since this is an implicit method, the NonlinearFunction (such as MassSpring) should ideally provide an implementation of evaluateDeriv to assist the Newton solver.

```cpp
auto rhs = std::make_shared<MassSpring>(1.0, 1.0);
ImplicitEuler stepper3(rhs);

for (int i = 0; i < steps; i++)
{
  stepper3.doStep(tau, y);
  outfile3 << (i+1) * tau << "  " << y(0) << " " << y(1) << std::endl;
}
```

## Crank-Nicolson Method

The Crank-Nicolson method is an implicit second-order time-stepping scheme. It can be viewed as the average of the Explicit and Implicit Euler methods or derived by applying the trapezoidal rule to the equivalent integral equation. It provides a better rate of convergence compared to first-order implicit schemes.

### Mathematical Overview

For a constant mesh size $\tau$, the update is defined by:
\begin{equation}
    y_{n+1} = y_n + \frac{\tau}{2} (f(y_n) + f(y_{n+1}))
\end{equation}

Finding $y_{n+1}$ requires solving the following non-linear system at each step:
\begin{equation}
    G(y_{n+1}) = y_{n+1} - y_n - \frac{\tau}{2} (f(y_n) + f(y_{n+1})) = 0
\end{equation}

This is achieved in our library using a Newton solver.

### Implementation

The Crank-Nicolson class is defined in `timestepper.hpp` and inherits from TimeStepper. It constructs the residual equation $G(y)$ using internal storage for previous values.

The constructor sets up the identity and constant functions needed to define the implicit equation:

```cpp
CrankNicholson(std::shared_ptr<NonlinearFunction> rhs) 
: TimeStepper(rhs), m_tau(std::make_shared<Parameter>(0.0)) 
{
  m_yold = std::make_shared<ConstantFunction>(rhs->dimX());
  m_fold = std::make_shared<ConstantFunction>(rhs->dimF());
  auto ynew = std::make_shared<IdentityFunction>(rhs->dimX());
  
  // Equation: G(y) = y_new - y_old - 0.5 * tau * (f(y_old) + f(y_new))
  m_equ = ynew - m_yold - m_tau * (0.5 * (m_fold + m_rhs));
}
```

### Parameters

$\tau$ : Time step size
$y$ : The current solution vector $y_n$, modified in-place by the Newton solver to become $y_{n+1}$.

The doStep implementation first evaluates the function at the current point to store $f(y_n)$ before triggering the solver:

```cpp
void doStep(double tau, VectorView<double> y) override
{
  Vector<> fold_val(m_rhs->dimF());
  m_rhs->evaluate(y, fold_val); // Store f(y_n)
  m_fold->set(fold_val);

  m_yold->set(y);
  m_tau->set(tau);
  NewtonSolver(m_equ, y); // Solve for y_{n+1}
}
```

### How can you use this implementation?

As shown in test_ode.cpp, you can initialize the Crank-Nicolson stepper and run the simulation similarly to other methods:

```cpp
auto rhs = std::make_shared<MassSpring>(1.0, 1.0);
CrankNicholson stepper4(rhs);

for (int i = 0; i < steps; i++)
{
  stepper4.doStep(tau, y);
  outfile4 << (i+1) * tau << "  " << y(0) << " " << y(1) << std::endl;
}
```

## Newton Solver

The Newton solver is the core numerical engine used by the implicit time-stepping schemes (Implicit Euler and Crank-Nicolson) in the ODE-to-CSE library. It is a powerful iterative method for finding the roots of nonlinear systems of equations.

### Mathematical Overview

The Newton solver finds a root $x$ such that $F(x) = 0$ by iteratively applying the following formula:
\begin{equation}
    x_{k+1} = x_k - J_F(x_k)^{-1} F(x_k)
\end{equation}
where $J_F(x_k)$ is the Jacobian matrix evaluated at $x_k$. The iteration terminates when the residual norm satisfies:
\begin{equation}
    \|F(x_{k+1})\| < \text{tol}
\end{equation}
or after reaching a maximum number of steps.

### Implementation
The NewtonSolver is implemented as a standalone function that operates on objects inheriting from the NonlinearFunction base class. This base class ensures that the solver can access both the function evaluation and its Jacobian.

```cpp
void NewtonSolver(std::shared_ptr<NonlinearFunction> func, 
                  VectorView<double> x, 
                  double tol = 1e-10, 
                  int maxsteps = 10, 
                  std::function<void(int,double,VectorView<double>)> callback = nullptr)
```

### Parameters

func: A NonlinearFunction object that provides evaluate() and evaluateDeriv() methods.
x: The initial guess for the root, modified in-place to contain the solution.
tol: The tolerance for convergence based on the residual norm.
maxsteps: The maximum number of iterations to perform.
callback: An optional function called at each iteration.

### How does it work?

At each iteration, the solver evaluates the Jacobian at the current point $x_k$.
In our implementation, this typically involves solving a linear system or inverting the matrix to find the update direction. If the residual norm does not reach the required tol within maxsteps, the solver throws a std::domain_error to indicate that convergence was not achieved.