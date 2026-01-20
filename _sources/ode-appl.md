# Applications

This page showcases a few use-cases of the ODE solvers from the ODE-to-CSE part of this project.
The focus is on *results*: plots and short simulation videos.

## Minimal background (why ODE solvers show up)

Many dynamical systems can be written as

$$
\dot y(t) = f\big(t, y(t)\big),
$$

which is exactly the interface used by the implemented time-steppers (Explicit/Improved/Implicit Euler, Crank–Nicolson, explicit and implicit Runge–Kutta).

For mechanical systems, a common modeling route is

- start from a potential energy $U(x)$,
- compute forces $F_i(x) = -\partial U/\partial x_i$,
- apply Newton’s law $m_i \ddot x_i = F_i(x)$,
- reduce to first order by introducing velocities $v_i = \dot x_i$.

In practice, the interesting part is what you can simulate with the resulting ODE (or DAE (differential-algebraic equation), if constraints are present).

## 1) Mass–spring system

The following clip demonstrates a mass–spring scene generated via the Python bindings.

<video controls width="900">
	<source src="_static/videos/mass_spring.mp4" type="video/mp4">
	Your browser does not support the video tag.
</video>

## 2) Mechanical systems with constraints (distance constraints)

Constraints turn many “geometry-driven” problems into differential-algebraic equations (DAEs).
A typical example is a pendulum: instead of describing the motion by an angle, you can enforce a fixed length constraint

$$
g(x) = \lVert x - x_0 \rVert^2 - l^2 = 0.
$$

In the implementation, such constraints can be enforced using Lagrange multipliers. With Python bindings and a scene graph (masses, springs, constraints), one can also assemble larger systems (e.g. crane-like structure).

### Double pendulum

<video controls width="900">
	<source src="_static/videos/double_pendulum.mp4" type="video/mp4">
	Your browser does not support the video tag.
</video>

### Chain

<video controls width="900">
	<source src="_static/videos/chain.mp4" type="video/mp4">
	Your browser does not support the video tag.
</video>

### Crane-like structure (vibration)

<video controls width="900">
	<source src="_static/videos/crane.mp4" type="video/mp4">
	Your browser does not support the video tag.
</video>

### Spinning top / Kreisel

<video controls width="900">
	<source src="_static/videos/kreisel.mp4" type="video/mp4">
	Your browser does not support the video tag.
</video>

### Spinning top / Kreisel (variant)

<video controls width="900">
	<source src="_static/videos/kreisel_notreally.mp4" type="video/mp4">
	Your browser does not support the video tag.
</video>

## How were the animations created?

The animations were produced from the Jupyter notebooks used to set up and simulate the mechanical systems (mass–spring, double pendulum, chain, crane, spinning top / Kreisel).

To reproduce them:

1. Navigate to the corresponding notebook files (the `.ipynb` notebooks found in the folder `mechsystem`).
2. Open a notebook in VS Code (or Jupyter) and execute the cells from top to bottom.
   - The first cells typically extend `sys.path` to load the compiled Python bindings from `../build/mechsystem`.
   - The subsequent cells build the scene (masses, springs, distance constraints) and run a simulation loop while updating a `pythreejs` visualization.
3. (Optional) Record the interactive 3D viewer (screen recording) to obtain the `.mp4` clips shown above.
