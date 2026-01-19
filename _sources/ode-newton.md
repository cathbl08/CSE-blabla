# Newton Solver

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