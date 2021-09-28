## Assemblers

For the assembly of the different operators (stiffness, mass) a two layered approach is chosen.
The operator assembler itself collects the local contributions and add's them into the global
system matrix. The local operator assembler calculates the contribution of each element.

Local assembler:

- `stiffness`: computes the stiffness contribution in terms of linear elasticity
- `consistentmass`: computes the full/consistent mass matrix contributions
- `hrzlumpedmass`: computes a lumped mass contribution by scaling the diagonal terms [[1]](#1)
- `lobattolumpedmass`: computes a lumped mass contribution based on a special quadrature

## Example

Constructing the stiffness operator:

```cpp
operatorType stiffnessMatrix;
double E = 1000000, nu = 0.3;
  
Elastodynamics::OperatorAssembler<Basis> operatorAssembler(basis);
operatorAssembler.initialize(stiffnessMatrix);
Elastodynamics::StiffnessAssembler stiffnessAssembler(E, nu);
operatorAssembler.assemble(stiffnessAssembler, stiffnessMatrix, false);
```

## References

<a id="1">[1]</a> 
Hinton E., Rock T., Zienkiewicz O. C. (1976). 
A note on mass lumping and related processes in the finite element method.
Earthquake Engin. and Struct. Dyn., 4
