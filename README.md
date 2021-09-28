# dune-elastodynamics
!!! PROTOTYPE !!!
Created during my master's thesis for testing partitioned coupling schemes with
the library preCICE.

This DUNE module is used to solve the linear equations of elastodynamics.
Assemblers for the calculation of stiffness, mass and lumped mass operators are provided.
Different explicit and implicit time-stepping methods for the solution of the corresponding
system of second order ordinary differential equations are implemented. 

## Installation

It is assumed the dune core modules and all other necessary libraries are allready installed.

## Governing equations

```latex
\begin{cases}
\rho \frac{\partial^2 u}{\partial t^2} + \nabla \cdot \sigma(u) = f & \text{(Newton's second law)} \\[10pt]
\sigma(u) = C : \epsilon(u) & \text{(Hooke's law)} \\[10pt]
\epsilon(u) = \frac{1}{2}(\nabla u + (\nabla u)^T) & \text{(Strain-displacement relation)} \\[10pt]
\end{cases} 
```

## Time-stepping methods

For the time-stepping, direct methods are implemented for solving second order
ordinary differential equations.

Explicit Runge-Kutta-Nyström methods with fixed timestep size:

- `nyström4`: fourth order Nyström method with fixed timestep
- `nyström5`: fitfh order Nyström method with fixed timestep

Explicit Runge-Kutta-Nyström methods with adaptive timestepping:

- `bettisrkn45`:
- `dprkn64`:
- `dprkn86`:

A popular approach in structural dynamics is the family of Newmark methods:

- `stoermer`:
- `foxgoodwin`:
- `linearacceleration`:
- `constantacceleration`:

