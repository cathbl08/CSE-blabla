# Crank-Nicolson Method

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

