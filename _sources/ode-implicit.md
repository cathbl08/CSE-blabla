# Implicit Euler Method

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

