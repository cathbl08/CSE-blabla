# Runge-Kutta Methods

Runge-Kutta (RK) methods are a family of one-step time integration schemes that achieve higher order accuracy by combining multiple evaluations (“stages”) of the right-hand side. This project covers an explicit and an implicit Runge-Kutta solver.

## Mathematical Overview

Recall the integral formulation of the ODE on time interval $[t_i, t_{i+1}]$:

$$
y(t_{i+1}) = y(t_i) + \int_{t_i}^{t_{i+1}} f(s, y(s))\, ds.
$$

Runge-Kutta methods obtain higher accuracy by approximating this integral with a quadrature rule:

$$
y_{i+1} = y_i + \tau \sum_{j=0}^{s-1} b_j\, f\big(t_i+c_j \tau,\, y_i^j\big),
$$

where

- $s$ is the number of stages,
- $c_j, b_j$ are nodes and weights of a quadrature rule on $[0,1]$,
- $y_i^j \approx y(t_i + c_j\tau)$ are the (unknown) stage values.

Stage values are defined by a second quadrature on $[0,c_j]$ (skipping the index $i$):

$$
y^j = y_i + \tau \sum_{l=0}^{s-1} a_{jl}\, f\big(t_i+c_l\tau,\, y^l\big), \qquad 0\le j < s.
$$

A particular RK method is defined by the coefficients $(A,b,c)$, typically given by the Butcher tableau

$$
\begin{array}{c|ccc}
c_0 & a_{0,0} & \dots & a_{0,s-1} \\
\vdots & \vdots & & \vdots \\
c_{s-1} & a_{s-1,0} & \dots & a_{s-1,s-1} \\
\hline
 & b_0 & \dots & b_{s-1}
\end{array}.
$$

Instead of solving for the stage values $y^j$, many implementations solve for stage slopes $k^j$:

$$
k^j = f\Big(t_i+c_j\tau,\, y_i + \tau\sum_{l=0}^{s-1} a_{jl} k^l\Big),
$$

and update

$$
y_{i+1} = y_i + \tau \sum_{l=0}^{s-1} b_l k^l.
$$

## Explicit RK methods

For explicit RK methods, $A$ is strictly lower triangular ($a_{jl}=0$ for $l\ge j$). This makes the stages sequential: $k^0$ can be computed first, then $k^1$, and so on.

## Implementation: `ExplicitRungeKutta`

The class `ExplicitRungeKutta` derives from `TimeStepper` and takes a right-hand side `NonlinearFunction` plus a Butcher tableau $(A,b,c)$.

Key design points in the implementation:

- **Dimension checks:** it requires `dimX == dimF` (ODE system $y' = f(y)$ with matching input/output dimension).
- **Explicitness check:** it validates that $A$ is strictly lower triangular (up to a small tolerance).
- **Stage storage:** all stage vectors $k^j\in\mathbb{R}^n$ are stored contiguously in a single vector `m_k` of length $s\,n$; stage $j$ is accessed via `m_k.range(j*n, (j+1)*n)`.

Algorithmically, each step computes

$$
y^{(i)} = y_i + \tau\sum_{j=0}^{i-1} a_{ij} k^j,\qquad k^i = f\big(y^{(i)}\big),
$$

and finally applies

$$
y_{i+1} = y_i + \tau\sum_{j=0}^{s-1} b_j k^j.
$$

The implementation mirrors this directly: for each stage $i$ it forms a temporary vector `m_y_tmp` and evaluates the RHS into the stage slice.

### Convenience tableau: classical RK4

The helper `ERK_RK4_Tableau()` returns the classical 4th order explicit RK coefficients $(A,b,c)$:

$$
\begin{array}{c|cccc}
0 \\
	\frac{1}{2} & \tfrac{1}{2} \\
	\frac{1}{2} & 0 & \tfrac{1}{2} \\
1 & 0 & 0 & 1 \\
\hline
 & \tfrac{1}{6} & \tfrac{1}{3} & \tfrac{1}{3} & \tfrac{1}{6}
\end{array}.
$$

Usage pattern:

```cpp
auto rhs = /* std::shared_ptr<NonlinearFunction> */;
auto [A, b, c] = ASC_ode::ERK_RK4_Tableau();
ASC_ode::ExplicitRungeKutta stepper(rhs, A, b, c);

stepper.doStep(tau, y);
```

## Implicit Runge-Kutta methods

Implicit RK methods use a general matrix $A$ and therefore require solving a coupled nonlinear system for the stage slopes $(k^0,\ldots,k^{s-1})$.

Gaussian collocation methods (Gauss-Legendre) are a common choice for high order accuracy; Radau IIA variants enforce a node at $c_{s-1}=1$ and are often preferred for stiff problems due to better damping (L-stability).

## Implementation: `ImplicitRungeKutta`

The class `ImplicitRungeKutta` also derives from `TimeStepper`, but it builds a nonlinear equation in $k\in\mathbb{R}^{sn}$ and solves it with the library Newton solver.

The coupled system is

$$
k - \widetilde f\big(\tilde y_i + \tau\, (A\otimes k)\big) = 0,
$$

where

- $k=(k^0,\ldots,k^{s-1})$ with each $k^j\in\mathbb{R}^n$,
- $\tilde y_i=(y_i,\ldots,y_i)\in\mathbb{R}^{sn}$ is the stacked current state,
- $\widetilde f(x_0,\ldots,x_{s-1})=(f(x_0),\ldots,f(x_{s-1}))$ applies the RHS stage-wise,
- $A\otimes$ is the block operator that mixes stage vectors using the RK matrix $A$.

In code, the constructor composes this equation from reusable `NonlinearFunction` building blocks:

- `MultipleFunc(rhs, s)` implements $\widetilde f$.
- `ConstantFunction(sn)` holds $\tilde y_i$ (updated each step).
- `Parameter` stores the current step size $\tau$.
- `MatVecFunc(A, n)` implements $A\otimes$ on a stacked stage vector.
- `IdentityFunction(sn)` represents the unknown $k$ (so the final equation has the form “identity minus something”).

During `doStep(tau, y)` the implementation:

1. Stacks the current state $y$ into $\tilde y_i$ and stores it in the constant function.
2. Sets the parameter $\tau$.
3. Solves for `m_k` with `NewtonSolver(m_equ, m_k)` (initial guess `m_k = 0`).
4. Updates $y \leftarrow y + \tau\sum_j b_j k^j$.

### Built-in tableaus in `implicitRK.hpp`

The header provides example nodes/weights:

- **Gauss-Legendre 2-stage**: `Gauss2a`, `Gauss2b`, `Gauss2c`.
- **Gauss-Legendre 3 nodes**: `Gauss3c` (nodes only; see the coefficient generation below).

If you include these headers in multiple translation units, prefer making such global objects `inline` to avoid ODR/linker issues.

## Computing tableaus from nodes (arbitrary order)

The helper `computeABfromC(c)` constructs $(A,b)$ from a chosen node vector $c$ by enforcing polynomial exactness up to degree $s-1$.

It builds a Vandermonde-like matrix

$$
M_{ij} = c_j^i\qquad (0\le i,j < s)
$$

and solves the moment conditions

$$
\sum_{l=0}^{s-1} b_l c_l^k = \int_0^1 x^k\,dx = \frac{1}{k+1},\qquad 0\le k < s,
$$

and for each stage $j$:

$$
\sum_{l=0}^{s-1} a_{jl} c_l^k = \int_0^{c_j} x^k\,dx = \frac{c_j^{k+1}}{k+1},\qquad 0\le k < s.
$$

In the implementation, the inverse of $M$ is computed (via `calcInverse`) and then used to form $b$ and each row of $A$.

## Quadrature helpers for Gauss / Radau nodes

The file `implicitRK.hpp` includes routines from [Numerical Recipes](https://numerical.recipes/book.html) to generate nodes and weights:

- `GaussLegendre(x, w)` for Gauss-Legendre nodes on $[0,1]$.
- `GaussJacobi(x, w, alf, bet)` for Gauss-Jacobi nodes (used to generate Radau nodes).
- `GaussRadau(x, w)` which computes Radau nodes on $[0,1]$ by calling `GaussJacobi` with parameters $(\alpha,\beta)=(1,0)$ and then appending the endpoint $1$.

These utilities let you generate $c$ (and optionally weights) and then use `computeABfromC(c)` to build the full tableau.

## How can you use this implementation?

As for the other time-steppers, you create a right-hand side object (derived from `NonlinearFunction`), construct the stepper, and call `doStep(tau, y)` in a loop.

### Explicit RK (example: RK4)

Explicit RK methods only need RHS evaluations `evaluate(y, f)`.

```cpp
#include <explicitRK.hpp>

using namespace ASC_ode;

auto rhs = std::make_shared<MassSpring>(1.0, 1.0);
auto [A, b, c] = ERK_RK4_Tableau();
ExplicitRungeKutta stepper(rhs, A, b, c);

for (int i = 0; i < steps; i++)
{
  stepper.doStep(tau, y);
  outfile << (i+1) * tau << "  " << y(0) << " " << y(1) << std::endl;
}
```

### Implicit RK (example: Gauss-Legendre 2)

Implicit RK methods solve a nonlinear system in each step using `NewtonSolver`. For good performance and robustness, your RHS should provide a Jacobian (either by implementing `evaluateDeriv` or by using the AutoDiff adapter discussed on the previous page).

```cpp
#include <implicitRK.hpp>

using namespace ASC_ode;

auto [Gauss3a,Gauss3b] = computeABfromC (Gauss3c); // Gauss3c is predefined in implicitRK.hpp
auto rhs = std::make_shared<MassSpring>(1.0, 1.0);
ImplicitRungeKutta stepper(rhs, Gauss3a, Gauss3b, Gauss3c);

for (int i = 0; i < steps; ++i)
{
    stepper.doStep(tau, y);
    outfile << (i+1) * tau << "  " << y(0) << " " << y(1) << std::endl;
}
```

If you want a different implicit method, generate nodes $c$ (e.g. using `GaussLegendre` or `GaussRadau`), compute `(A,b)` via `computeABfromC(c)`, and then construct `ImplicitRungeKutta(rhs, A, b, c)`.


## Plots of mass-spring system and RC circuit

### Explicit Runge-Kutta (RK4)
<img src="images/ExplicitRungeKutta.png" alt="Time Evolution and Phase Plot of Mass-Spring System solved with Explicit Runge-Kutta RK4"  width="800px">

<img src="images/RC_ExplicitRungeKutta.png" alt="Time Evolution of RC Circuit solved with Explicit Runge-Kutta RK4"  width="600px">

### Implicit Runge-Kutta (Gauss-Legendre 2)
<img src="images/ImplicitRungeKutta.png" alt="Time Evolution and Phase Plot of Mass-Spring System solved with Implicit Runge-Kutta Gauss-Legendre2"  width="800px">

<img src="images/RC_ImplicitRungeKutta.png" alt="Time Evolution of RC Circuit solved with Implicit Runge-Kutta Gauss-Legendre2"  width="600px">
