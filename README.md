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
```C
// create basis
int p = 2, k = 3;
auto basis = makeBasis(gridView, power<dim>(lagrange<p>()));
```
This basis will be used by an operator assembler, sorting local computations into
the global matrix (stiffness, mass, ...). For e.g. a sparse stiffness matrix, the operator
assembler is first created with an appropriate basis. Afterwards the datastructure, which
should hold the entries is initialized with the operator assembler. In a third step the
operator is assembled by using a local assembler, the datastructure and a boolean indicating
if the matrix is lumped or not. The local assembler in this case for the stiffness is set
up with a given Young's modulus and Poisson ratio.
```C
// assemble stiffness
operatorType stiffnessMatrix;
double E = 1098500.0, nu = 0.3;
  
OperatorAssembler operatorAssembler(basis);
StiffnessAssembler stiffnessAssembler(E, nu);
operatorAssembler.initialize(stiffnessMatrix);
operatorAssembler.assemble(stiffnessAssembler, stiffnessMatrix, false);
```

Global assembler:

- `operator`:

Stiffness assembler:

- `stiffness`:

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

