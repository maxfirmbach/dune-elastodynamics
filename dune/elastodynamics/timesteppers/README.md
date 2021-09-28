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
