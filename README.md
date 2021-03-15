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

For the discretization in space, finite elements are used. The assemblers work with
function space basis trees, e.g. a lagragian basis of order p and dimension k is
constructed as:

```
auto basis = makeBasis(gridView, power<dim>(lagrange<p>()));
```

Stiffness assembler:

- `stiffness`

Mass assembler:

- `consistentmass`:
- `hrzlumpedmass`:
- `lobattolumpedmass`:

## Time-stepping methods

For the time-stepping, direct methods are implemented for solving second order
ordinary differential equations.

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

