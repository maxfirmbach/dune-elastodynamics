// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:

#include <config.h>

#include <dune/common/parallel/mpihelper.hh>

#include <dune/grid/uggrid.hh>
#include <dune/grid/io/file/gmshreader.hh>

#include <dune/istl/matrix.hh>
#include <dune/istl/bcrsmatrix.hh>
#include <dune/istl/bvector.hh>
#include <dune/istl/preconditioners.hh>
#include <dune/istl/solvers.hh>

#include <dune/functions/functionspacebases/basistags.hh>
#include <dune/functions/functionspacebases/powerbasis.hh>
#include <dune/functions/functionspacebases/lagrangebasis.hh>

#include <dune/elastodynamics/assemblers/operatorassembler.hh>
#include <dune/elastodynamics/assemblers/stiffnessassembler.hh>

#include <dune/elastodynamics/utilities/boundaryindexbcassembler.hh>

using namespace Dune;
const int dim = 2;
const int p = 2;

int main(int argc, char** argv) {

  const MPIHelper& mpiHelper = MPIHelper::instance(argc, argv);
  bool passed = true;
  
  // generate Grid
  using Grid = UGGrid<dim>;
  using GridView = Grid::LeafGridView;
  
  auto mesh = "beam.msh";
  std::vector<int> materialIndex, boundaryIndex;
  GridFactory<Grid> factory;
  GmshReader<Grid>::read(factory, mesh, boundaryIndex, materialIndex, true);
  std::shared_ptr<Grid> grid(factory.createGrid());    
  auto gridView = grid->leafGridView();
  
  // generate Basis
  using namespace Functions::BasisBuilder;
  auto basis = makeBasis(gridView, power<dim>(lagrange<p>()));
  using Basis = decltype(basis);

  // define operators needed
  using operatorType = BCRSMatrix<FieldMatrix<double, dim, dim>>;
  using blockVector  = BlockVector<FieldVector<double, dim>>;

  // assemble problem
  operatorType stiffnessMatrix;
  double E = 1000000, nu = 0.3;
  
  Elastodynamics::OperatorAssembler<Basis> operatorAssembler(basis);
  operatorAssembler.initialize(stiffnessMatrix);
  Elastodynamics::StiffnessAssembler stiffnessAssembler(E, nu);
  operatorAssembler.assemble(stiffnessAssembler, stiffnessMatrix, false);

  // set boundary conditions
  blockVector loadVector(basis.size());
  loadVector = 0.0;
  
  // two nodes at free end, force should sum up to 1N
  FieldVector<double, dim> force = {0.0, 0.5};
  
  Elastodynamics::BoundaryIndexBCAssembler<Basis> bcAssembler(basis, boundaryIndex);
  bcAssembler.assembleMatrix(stiffnessMatrix);
  bcAssembler.assembleVector(loadVector, force);

  // solve linear system
  blockVector x(basis.size());
  x = 0.0;

  MatrixAdapter<operatorType, blockVector, blockVector> stiffnessOperator(stiffnessMatrix);
  SeqILU<operatorType, blockVector, blockVector> preconditioner(stiffnessMatrix, 1.0);
  CGSolver<blockVector> solver(stiffnessOperator, preconditioner, 1e-16, 10000, 1);
  
  InverseOperatorResult statistics;
  solver.apply(x, loadVector, statistics);

  // exact displacement at beam tip x = F*L^3 / 3*E*I, I = b*h^3 / 12
  // thus get maxium displacement  
  double max = 0.0;
  for( int i=0; i<basis.size(); i++) {
    for( int k=0; k<dim; k++) {
      max = std::max(max, x[i][k]);
    }
  }

  // analytical solution 0.108
  passed = passed and std::abs(0.1080-max) < 1e-4;

  return passed ? 0 : 1;

}
