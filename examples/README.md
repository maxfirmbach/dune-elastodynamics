# Examples
The dune-elastodynamics module contains several examples for static and dynamic linear
elasticity cases.

## Example01
In this example, a slender cantilever beam is put under a unit load at the free end.
The overall displacement of the beam is calculated. One can also check the analytical
solution of the Euler-Bernoulli or Timoschenko beam theory.

### Setup
The structured mesh is generated with gmsh and consists of .
<figure>
	<center>
		<img src='Example01/Setup.png'>
	</center>
</figure>

### Calculation
The main new features are:
- setting up the necessary data structures
- load the grid
- creating the function space basis tree
- assembling the problem and boundary conditions
- solving the linear system
- getting output for paraview

### Result
<figure>
	<center>
		<img src='Example01/Result.png'>
	</center>
</figure>

## Example02

## Example03
