// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:

// This is an example file for a simulation script with truss
// elements to solve static loadcases ... 
// can easily be advanced for dynamic loadcases

#include <config.h>

#include <stdio.h>

#include <dune/grid/io/file/vtk.hh>
#include <dune/grid/io/file/gmshreader.hh>

#include <dune/foamgrid/foamgrid.hh>

#include <dune/istl/matrix.hh>
#include <dune/istl/bvector.hh>
#include <dune/istl/bcrsmatrix.hh>
#include <dune/istl/bdmatrix.hh>
#include <dune/istl/io.hh>

#include <dune/functions/functionspacebases/powerbasis.hh>
#include <dune/functions/functionspacebases/lagrangebasis.hh>
#include <dune/functions/gridfunctions/discreteglobalbasisfunction.hh>

#include <dune/elastodynamics/assemblers/operatorassembler.hh>
// here include your stiffness/mass element assembler

using namespace Dune;

int main (int argc, char *argv[]) {

  // set up MPI (is always needed)
  const MPIHelper& mpiHelper = MPIHelper::instance(argc, argv);

  const int p        = 1; // element order
  const int dim      = 1; // element dimension
  const int dimworld = 2; // world dimension

  using Grid = FoamGrid<dim, dimworld>;
  using GridView = Grid::LeafGridView;
  
  // load and generate grid / gridview
  GridFactory<Grid> factory;
  GmshReader<Grid>::read(factory, "framework_2.msh");
  std::shared_ptr<Grid> grid(factory.createGrid());
  GridView gridView = grid->leafGridView();
  
  // generate Basis
  using namespace Functions::BasisFactory;
  auto basis = makeBasis(gridView, power<dimworld>(lagrange<p>()));
  using Basis = decltype(basis);
  
  // define operators needed
  using operatorType = BCRSMatrix<FieldMatrix<double, dimworld, dimworld>>;
  using blockVector  = BlockVector<FieldVector<double, dimworld>>;
  
  // assemble problem
  //
  // This would be student part
  //
  // Here you should assemble the stiffness matrix of the problem
  // and maybe also the mass matrix for dynamic problems
  //
  // Example for 2D/3D solid elements for stiffness and mass matrix:
  //
  // Elastodynamics::OperatorAssembler<Basis> operatorAssembler(basis);
  //
  // double E = 4000000.0, nu = 0.3;
  // operatorType stiffnessMatrix;
  // operatorAssembler.initialize(stiffnessMatrix);
  // Elastodynamics::StiffnessAssembler stiffnessAssembler(E, nu);
  // operatorAssembler.assemble(stiffnessAssembler, stiffnessMatrix, false);
  //
  // double rho = 1.0;
  // operatorType massMatrix;
  // operatorAssembler.initialize(massMatrix);
  // Elastodynamics::ConsistentMassAssembler massAssembler(rho);
  // operatorAssembler.assemble(massAssembler, massMatrix, false);
  //
  
  blockVector loadVector(basis.size()), displacementVector(basis.size());
  loadVector = 0.0, displacementVector = 0.0;
  
  // set dirichlet boundary
  //
  // This would be student part
  //
  
  // set neumann boundary
  //
  // This would be student part
  //
  
  // solve the arising linear system Ku=f
  // with an appropriate method for static cases
  //
  // if dynamic set a time loop and an appropriate
  // integrator to evolve the ODE in time
  //
  // This would be student part
  //
  
  // generate global displacement function
  using displacementRange = FieldVector<double, dimworld>;
  auto displacementFunction = Functions::makeDiscreteGlobalBasisFunction<displacementRange> (basis, displacementVector);
  
  // generate output writer (refine output mesh depending on element order)
  SubsamplingVTKWriter<GridView> vtkWriter(gridView, refinementLevels(p));
  vtkWriter.addVertexData(displacementFunction, VTK::FieldInfo("displacement", VTK::FieldInfo::Type::vector, dimworld));
  vtkWriter.write("framework-result");

}
