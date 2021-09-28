// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:

#include <config.h>

#include <dune/common/parallel/mpihelper.hh>

#include <dune/grid/uggrid.hh>
#include <dune/grid/io/file/gmshreader.hh>

#include <dune/istl/matrix.hh>
#include <dune/istl/matrixmatrix.hh>
#include <dune/istl/bcrsmatrix.hh>
#include <dune/istl/bdmatrix.hh>
#include <dune/istl/bvector.hh>
#include <dune/istl/preconditioners.hh>
#include <dune/istl/solvers.hh>
#include <dune/istl/eigenvalue/poweriteration.hh>
#include <dune/istl/umfpack.hh>

#include <dune/functions/functionspacebases/basistags.hh>
#include <dune/functions/functionspacebases/powerbasis.hh>
#include <dune/functions/functionspacebases/lagrangebasis.hh>

#include <dune/elastodynamics/assemblers/operatorassembler.hh>
#include <dune/elastodynamics/assemblers/stiffnessassembler.hh>
#include <dune/elastodynamics/assemblers/hrzlumpedmassassembler.hh>

#include <dune/elastodynamics/utilities/boundaryindexbcassembler.hh>

using namespace Dune;
const int dim = 2;
const int p = 2;

int main(int argc, char** argv) {

  Dune::MPIHelper& mpiHelper = Dune::MPIHelper::instance(argc, argv);
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
  using diagonalType = BDMatrix<FieldMatrix<double, dim, dim>>;
  using blockVector  = BlockVector<FieldVector<double, dim>>;

  {
    // assemble problem
    Elastodynamics::OperatorAssembler<Basis> operatorAssembler(basis);
  
    double E = 1000000, nu = 0.3;
    operatorType stiffnessMatrix;
    operatorAssembler.initialize(stiffnessMatrix);
    Elastodynamics::StiffnessAssembler stiffnessAssembler(E, nu);
    operatorAssembler.assemble(stiffnessAssembler, stiffnessMatrix, false);
  
    double rho = 1.0;
    diagonalType massMatrix(basis.size());
    Elastodynamics::HRZLumpedMassAssembler massAssembler(rho);
    operatorAssembler.assemble(massAssembler, massMatrix, true);

    // set boundary conditions
    Elastodynamics::BoundaryIndexBCAssembler<Basis> bcAssembler(basis, boundaryIndex);
    bcAssembler.assembleMatrix(stiffnessMatrix);
    bcAssembler.assembleMatrix(massMatrix);

    // calculate eigenvalue problem of dynamic equation Kx=c^2Mx
    // analytical solution ~5.6387
    double lambda;
    blockVector x(basis.size());
    operatorType dynamicMatrix;
    massMatrix.invert();
    matMultMat(dynamicMatrix, massMatrix, stiffnessMatrix);
  
    PowerIteration_Algorithms<operatorType, blockVector> eigenvalueAlgorithm(dynamicMatrix, 10000, 1);  
    auto iterationMatrix = eigenvalueAlgorithm.getIterationMatrix();
    UMFPack<operatorType> solver(iterationMatrix);
  
    x = 1.0;
    eigenvalueAlgorithm.applyInverseIteration(20.0, 1e-6, solver, x, lambda);
    std::cout << "smallest eigenvalue: " << std::sqrt(lambda) << std::endl;
  
    passed = passed and std::abs(5.638-std::sqrt(lambda)) < 1e-2;
  }

  return passed ? 0 : 1;

}
