// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:

#include <config.h>

#include <dune/common/parallel/mpihelper.hh>

#include <dune/grid/uggrid.hh>
#include <dune/grid/io/file/gmshreader.hh>
#include <dune/grid/io/file/vtk.hh>

#include <dune/istl/matrixmarket.hh>
#include <dune/istl/matrix.hh>
#include <dune/istl/bvector.hh>
#include <dune/istl/bcrsmatrix.hh>
#include <dune/istl/preconditioners.hh>
#include <dune/istl/solvers.hh>

#include <dune/functions/functionspacebases/basistags.hh>
#include <dune/functions/functionspacebases/powerbasis.hh>
#include <dune/functions/functionspacebases/lagrangebasis.hh>

#include <dune/functions/functionspacebases/basistags.hh>
#include <dune/functions/functionspacebases/powerbasis.hh>
#include <dune/functions/functionspacebases/lagrangebasis.hh>
#include <dune/functions/gridfunctions/discreteglobalbasisfunction.hh>

#include <dune/elastodynamics/assemblers/operatorassembler.hh>
#include <dune/elastodynamics/assemblers/stiffnessassembler.hh>
#include <dune/elastodynamics/assemblers/consistentmassassembler.hh>

#include <dune/elastodynamics/utilities/boundaryindexbcassembler.hh>

#include <dune/elastodynamics/timesteppers/coefficients.hh>
#include <dune/elastodynamics/timesteppers/timestepcontroller.hh>
#include <dune/elastodynamics/timesteppers/newmark.hh>

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
  GmshReader<Grid>::read(factory, mesh, boundaryIndex, materialIndex);
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
  Elastodynamics::OperatorAssembler<Basis> operatorAssembler(basis);
  
  double E = 3.0e7, nu = 0.3;
  operatorType stiffnessMatrix;
  operatorAssembler.initialize(stiffnessMatrix);
  Elastodynamics::StiffnessAssembler stiffnessAssembler(E, nu);
  operatorAssembler.assemble(stiffnessAssembler, stiffnessMatrix, false);
  
  double rho = 0.3/386.0;
  operatorType massMatrix;
  operatorAssembler.initialize(massMatrix);
  Elastodynamics::ConsistentMassAssembler massAssembler(rho);
  operatorAssembler.assemble(massAssembler, massMatrix, false);

  blockVector loadVector(basis.size()), displacement(basis.size()), velocity(basis.size()), acceleration(basis.size());
  loadVector = 0.0, displacement = 0.0, velocity = 0.0, acceleration = 0.0;

  // set boundary and initial conditions
  Elastodynamics::BoundaryIndexBCAssembler<Basis> bcAssembler(basis, boundaryIndex);
  bcAssembler.assembleMatrix(stiffnessMatrix);
  bcAssembler.assembleMatrix(massMatrix);
  loadMatrixMarket(displacement, "Displacement.mm");
  
  // generate output writer
  using displacementRange = Dune::FieldVector<double, dim>;
  Functions::LagrangeBasis<GridView, p> Pbasis(gridView);
  auto displacementFunction = Functions::makeDiscreteGlobalBasisFunction<displacementRange> (Pbasis, displacement);
    
  auto vtkWriter = std::make_shared<SubsamplingVTKWriter<GridView>> (gridView, refinementLevels(2));
  VTKSequenceWriter<GridView> vtkSequenceWriter(vtkWriter, "solid");
  vtkWriter->addVertexData(displacementFunction, VTK::FieldInfo("displacement", VTK::FieldInfo::Type::vector, dim));
  vtkSequenceWriter.write(0.0);
    
  // setup looping
  double t = 0.0;
  double dt = 0.00025;
  int iteration_count = 0;
  
  // set up Newmark-beta method
  FixedStepController fixed(t, dt);
  NewmarkCoefficients coefficients = ConstantAcceleration();
  Newmark<GridView, operatorType, blockVector> newmark(gridView, massMatrix, stiffnessMatrix, coefficients, fixed);
  newmark.initialize(acceleration, loadVector);
    
  while(t < 0.025) {
    std::cout << "time " << t << std::endl;
    std::cout << "iter " << iteration_count << std::endl; 
    
    newmark.step(displacement, velocity, acceleration, loadVector);
    vtkSequenceWriter.write(t);

    t += dt;
    iteration_count += 1;

  }

  return 1;
} 
