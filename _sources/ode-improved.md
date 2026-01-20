# Improved Euler Method

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

### Plots of mass-spring system and RC circuit

<img src="images/ImprovedEuler.png" alt="Time Evolution and Phase Plot of Mass-Spring System solved with Improved Euler Method"  width="800px">

<img src="images/RC_ImprovedEuler.png" alt="Time Evolution of RC Circuit solved with Improved Euler Method"  width="600px">


