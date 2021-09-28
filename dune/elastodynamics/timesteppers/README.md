## Time-stepping methods

For the time-stepping, direct methods are implemented for solving second order
ordinary differential equations without converting them into a system of first
order equations.



Explicit Runge-Kutta-Nyström methods with fixed timestep size (expect a lumped
mass matrix $`M`$):

- `nyström4`: fourth order method with fixed timestep [[1]](#1)
- `nyström5`: fitfh order method with fixed timestep [[1]](#1)

Explicit Runge-Kutta-Nyström methods with adaptive timestepping (expect a lumped
mass matrix $`M`$):

- `bettisrkn45`: fitfh order method with fourth order error estimation [[2]](#2)
- `dprkn64`: sixth order method with fourth order error estimation [[3]](#3)
- `dprkn86`: eigth order method with sixth order error estimation [[3]](#3)

A popular approach in structural dynamics is the family of Newmark methods:

- `stoermer`: explicit central difference method
- `foxgoodwin`: fourth order conditionally stable implicit method
- `linearacceleration`: second order conditionally stable implicit method
- `constantacceleration`: second order unconditionally stable implicit method

## References

<a id="1">[1]</a> 
Hairer E., Norsett S. P., Wanner G. (1993). 
Solving Ordinary Differential Equations I - Nonstiff Problems.

<a id="2">[2]</a> 
Bettis D. G. (1973). 
A Runge-Kutta Nyström algorithm.
Celestial mechanics, 8, 229-233.

<a id="3">[3]</a> 
Dormand J. R., El-Mikkawy M. E., Prince P. J. (1987). 
High-Order Embedded Runge-Kutta-Nystrom Formulae.
IMA Journal of Numerical Analysis, 7(4), 423-430.

<a id="3">[3]</a> 
Dormand J. R., El-Mikkawy M. E., Prince P. J. (1987). 
High-Order Embedded Runge-Kutta-Nystrom Formulae.
IMA Journal of Numerical Analysis, 7(4), 423-430.

<a id="4">[4]</a> 
Newmark N. M. (1959). 
A method of computation for structural dynamics.
Journal of the Engineering Mechanics Division, 85(EM3), 67-94.
