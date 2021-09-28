## Time-stepping methods

For the time-stepping, direct methods are implemented for solving second order
ordinary differential equations without converting them into a system of first
order equations.



Explicit Runge-Kutta-Nyström methods with fixed timestep size (expect a lumped
mass matrix $`M`$):

- `nyström4`: fourth order method with fixed timestep
- `nyström5`: fitfh order method with fixed timestep

Explicit Runge-Kutta-Nyström methods with adaptive timestepping (expect a lumped
mass matrix $`M`$):

- `bettisrkn45`: fitfh order method with fourth order error estimation [[1]](#1)
- `dprkn64`: sixth order method with fourth order error estimation
- `dprkn86`: eigth order method with sixth order error estimation

A popular approach in structural dynamics is the family of Newmark methods:

- `stoermer`: explicit central difference method
- `foxgoodwin`: fourth order conditionally stable implicit method
- `linearacceleration`: second order conditionally stable implicit method
- `constantacceleration`: second order unconditionally stable implicit method

## References

<a id="1">[1]</a> 
Bettis, D. G. (1973). 
A Runge-Kutta Nyström algorithm.
Celestial mechanics, 8, 229-233.
