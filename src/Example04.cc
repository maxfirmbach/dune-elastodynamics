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
#include <dune/istl/umfpack.hh>

#include <dune/functions/functionspacebases/powerbasis.hh>
#include <dune/functions/functionspacebases/lagrangebasis.hh>
#include <dune/functions/functionspacebases/interpolate.hh>
#include <dune/functions/gridfunctions/discreteglobalbasisfunction.hh>

#include <dune/elastodynamics/assemblers/operatorassembler.hh>
#include <dune/elastodynamics/assemblers/stiffnessassembler_truss.hh>

using namespace Dune;

int main (int argc, char *argv[]) {

  // set up MPI (is always needed)
  const MPIHelper& mpiHelper = MPIHelper::instance(argc, argv);

  // load and generate grid / gridview
  const int p        = 2; // element order
  const int dim      = 1; // element dimension
  const int dimworld = 2; // world dimension

  using Grid = FoamGrid<dim, dimworld>;
  using GridView = Grid::LeafGridView;
  
  GridFactory<Grid> factory;
  GmshReader<Grid>::read(factory, "framework_3.msh");
  std::shared_ptr<Grid> grid(factory.createGrid());
  GridView gridView = grid->leafGridView();
  
  // generate Basis
  std::cout << "Basis Order = " << p << std::endl;
  using namespace Functions::BasisFactory;
  auto basis = makeBasis(gridView, power<dimworld>(lagrange<p>()));
  using Basis = decltype(basis);
  
  // define operators needed
  using operatorType = BCRSMatrix<FieldMatrix<double, dimworld, dimworld>>;
  using blockVector  = BlockVector<FieldVector<double, dimworld>>;

  // assemble problem
  Elastodynamics::OperatorAssembler<Basis> operatorAssembler(basis);
  double E = 1e7, A = 0.1;

  operatorType stiffnessMatrix;
  operatorAssembler.initialize(stiffnessMatrix);
  Elastodynamics::StiffnessAssemblerTruss stiffnessAssembler(E, A);
  operatorAssembler.assemble(stiffnessAssembler, stiffnessMatrix, false);

  blockVector loadVector(basis.size()), displacementVector(basis.size());
  loadVector = 0.0, displacementVector = 0.0;
  
  // set dirichlet boundary
  auto displacementPredicate = [](auto coordinate) { return coordinate[1] > -0.001; };
  blockVector dirichletNodes(basis.size());
  Functions::interpolate(basis, dirichletNodes, displacementPredicate);

  FieldMatrix<double, dimworld, dimworld> I = {{1,0},{0,1}};
  FieldMatrix<double, dimworld, dimworld> O = {{0,0},{0,0}};

  for (int i=0; i < stiffnessMatrix.N(); i++) {
    if (dirichletNodes[i][0]) {
      auto cIt = stiffnessMatrix[i].begin();
      auto cEndIt = stiffnessMatrix[i].end();
      for (; cIt!=cEndIt; ++cIt)
        *cIt = (cIt.index()==i) ? I : O;
    }
  }
  
  // set neumann boundary
  auto loadPredicate = [](auto coordinate) { return coordinate[1] < -8.500; };
  blockVector neumannNodes(basis.size());
  Functions::interpolate(basis, neumannNodes, loadPredicate);
  
  FieldVector<double, dimworld> L = {0, -1732};
  
  for (int i=0; i < neumannNodes.size(); i++) {
    if (neumannNodes[i][0]) {
      loadVector[i] = L;
    }
  }

  // solve the arising linear system Ku=f 
  UMFPack<operatorType> umfpack(stiffnessMatrix);
  umfpack.setVerbosity(2);
  
  InverseOperatorResult statistics;
  umfpack.apply(displacementVector, loadVector, statistics);

  // generate global displacement function
  using displacementRange = FieldVector<double, dimworld>;
  auto displacementFunction = Functions::makeDiscreteGlobalBasisFunction<displacementRange> (basis, displacementVector);
  
  // generate output writer (refine output mesh depending on element order)
  SubsamplingVTKWriter<GridView> vtkWriter(gridView, refinementLevels(p));
  vtkWriter.addVertexData(displacementFunction, VTK::FieldInfo("displacement", VTK::FieldInfo::Type::vector, dimworld));
  vtkWriter.write("framework-result");
}
