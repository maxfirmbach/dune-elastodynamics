// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:

#include <config.h>

#include <dune/common/parallel/mpihelper.hh>

#include <dune/grid/io/file/vtk.hh>
#include <dune/grid/io/file/gmshreader.hh>

#include <dune/foamgrid/foamgrid.hh>

#include <dune/istl/matrix.hh>
#include <dune/istl/bvector.hh>
#include <dune/istl/bcrsmatrix.hh>
#include <dune/istl/bdmatrix.hh>
#include <dune/istl/io.hh>
//#include <dune/istl/solvers.hh> //algebraic system solve
//#include <dune/istl/preconditioners.hh> //algebraic system solve
#include <dune/istl/umfpack.hh>

#include <dune/functions/functionspacebases/powerbasis.hh>
#include <dune/functions/functionspacebases/lagrangebasis.hh>
#include <dune/functions/functionspacebases/interpolate.hh>
#include <dune/functions/gridfunctions/discreteglobalbasisfunction.hh>

#include <dune/elastodynamics/assemblers/operatorassembler.hh>
#include <dune/elastodynamics/assemblers/stiffnessassembler_truss.hh>
// here include your stiffness/mass element assembler

using namespace Dune;

int main (int argc, char *argv[]) {

  // set up MPI (is always needed)
  const MPIHelper& mpiHelper = MPIHelper::instance(argc, argv);
  bool passed = true;

  // load and generate grid / gridview
  const int p        = 1; // element order
  const int dim      = 1; // element dimension
  const int dimworld = 2; // world dimension

  using Grid = FoamGrid<dim, dimworld>;
  using GridView = Grid::LeafGridView;
  
  GridFactory<Grid> factory;
  GmshReader<Grid>::read(factory, "framework_1.msh");
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
  Elastodynamics::OperatorAssembler<Basis> operatorAssembler(basis);
  double E = 2e8, A = 0.005;
  operatorType stiffnessMatrix;
  operatorAssembler.initialize(stiffnessMatrix);
  Elastodynamics::StiffnessAssemblerTruss stiffnessAssembler(E, A);
  operatorAssembler.assemble(stiffnessAssembler, stiffnessMatrix, false);
  
  blockVector loadVector(basis.size()), displacementVector(basis.size());
  loadVector = 0.0, displacementVector = 0.0;
  
  // set dirichlet boundary
  auto verschiebungPredicate = [](auto coordinate) { return coordinate[1] > -1; };
  blockVector dirichletNodes(basis.size());
  Functions::interpolate(basis, dirichletNodes, verschiebungPredicate);
   
  FieldMatrix<double, dimworld, dimworld> I = {{1, 0}, {0, 1}};  
  FieldMatrix<double, dimworld, dimworld> O = {{0, 0}, {0, 0}};

  for (int i=0; i < stiffnessMatrix.N(); i++) {
    if (dirichletNodes[i][0]) {
      auto cIt = stiffnessMatrix[i].begin();
      auto cEndIt = stiffnessMatrix[i].end();
      for (; cIt!=cEndIt; ++cIt)
      *cIt = (cIt.index()==i) ? I : O;
    }
  }  

  // set neumann boundary
  auto loadPredicate = [](auto coordinate) { return coordinate[1] < -1; };
  blockVector neumannNodes(basis.size());
  Functions::interpolate(basis, neumannNodes, loadPredicate);
  
  FieldVector<double, dimworld> L = {0, -1732};
  
  for (int i=0; i < neumannNodes.size(); i++)
  {
    if (neumannNodes[i][0])
    {
      loadVector[i] = L;
    }
  }
  
  UMFPack<operatorType> solver(stiffnessMatrix);
  
  InverseOperatorResult statistics;
  solver.apply(displacementVector, loadVector, statistics);
  
  // get maxium displacement  
  double max = 0.0;
  for( int i=0; i<basis.size(); i++) {
    for( int k=0; k<dimworld; k++) {
      max = std::max(max, std::abs(displacementVector[i][k]));
    }
  }

  // analytical solution 
  passed = passed and std::abs(max-0.01155) < 1e-5;

  return passed ? 0 : 1;

}
