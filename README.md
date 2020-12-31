# dune-elastodynamics
This DUNE module is used to solve the linear equations of elastodynamics.
Assemblers for the calculation of stiffness, mass and lumped mass operators are provided.
Different explicit and implicit time-stepping methods for the solution of the corresponding
system of second order ordinary differential equations are implemented.

## Governing equations

$
    \begin{cases}
      \rho \frac{\partial^2 u}{\partial t^2} + \nabla \cdot \sigma(u) = f & \text{(Newton's second law)} \\[10pt]
      \sigma(u) = C : \epsilon(u) & \text{(Hooke's law)} \\[10pt]
      \epsilon(u) = \frac{1}{2}(\nabla u + (\nabla u)^T) & \text{(Strain-displacement relation)} \\[10pt]
    \end{cases} 
$

## Assemblers





The discretized system can be written as:


## Time-stepping methods

Second-order ordinary differential equations of the form

can be solved directly.

Explicit Runge-Kutta-Nyström methods with fixed timestep size:

- `Nyström4`:
- `Nyström5`:

Explicit Runge-Kutta-Nyström methods with adaptive timestepping:

- `BettisRKN45`:
- `DPRKN64`:
- `DPRKN86`:

A popular approach in structural dynamics is the family of Newmark methods:

- `Stoermer`:
- `FoxGoodwin`:
- `LinearAcceleration`:
- `ConstantAcceleration`:

