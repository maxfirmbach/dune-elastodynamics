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
#include <dune/istl/solvers.hh>
#include <dune/istl/preconditioners.hh>
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
  const int p        = 1; // element order
  const int dim      = 1; // element dimension
  const int dimworld = 3; // world dimension

  using Grid = FoamGrid<dim, dimworld>;
  using GridView = Grid::LeafGridView;
  
  GridFactory<Grid> factory;
  GmshReader<Grid>::read(factory, "bridge.msh");
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
  double E = 210000000, A = 0.01;
  operatorType stiffnessMatrix;
  operatorAssembler.initialize(stiffnessMatrix);
  Elastodynamics::StiffnessAssemblerTrussQuadrature stiffnessAssembler(E, A);
  operatorAssembler.assemble(stiffnessAssembler, stiffnessMatrix, false);
  
  blockVector loadVector(basis.size()), displacementVector(basis.size());
  loadVector = 0.0, displacementVector = 0.0;

  // set dirichlet boundary
  auto verschiebungPredicate = [](auto coordinate)
  {
    return coordinate[0] < 3.16
        || coordinate[0] > 741.1
        || coordinate[1] < -48.6;
  };
  blockVector dirichletNodes(basis.size());
  Functions::interpolate(basis, dirichletNodes, verschiebungPredicate);
   
  FieldMatrix<double, dimworld, dimworld> I = {{1, 0, 0}, {0, 1, 0}, {0, 0, 1}};   
  FieldMatrix<double, dimworld, dimworld> O = {{0, 0, 0}, {0, 0, 0}, {0, 0, 0}};

  for (int i=0; i < stiffnessMatrix.N(); i++) {
    if (dirichletNodes[i][0]) {
      auto cIt = stiffnessMatrix[i].begin();
      auto cEndIt = stiffnessMatrix[i].end();
      for (; cIt!=cEndIt; ++cIt)
      *cIt = (cIt.index()==i) ? I : O;
    }
  }  

  // set neumann boundary
  auto loadonePredicate = [](auto coordinate) //Function erstellen
  {
    return coordinate[0] > 3.15 && coordinate[0] < 741.2 && coordinate[1] > 0 && coordinate[1] < 70 ; 
    
    //markierung von dem Kraft am Ende wirken soll coordinate[0] > 3.15 && coordinate[0] < 741.2 && coordinate[1] > 0 && coordinate[1] < 17 verkehrlast   
    //coordinate[0] > 3.15 && coordinate[0] < 741.2 && coordinate[1] > 0 && coordinate[1] < 70 windlast  
    //coordinate[0] > 100 && coordinate[0] < 150 && coordinate[1] > 0 && coordinate[1] < 17 last auf linken seite
  };
  blockVector neumannNodes(basis.size()); //bool vector erstellen
  Functions::interpolate(basis, neumannNodes, loadonePredicate);

  FieldVector<double, dimworld> Lone = {0, 0, -0.5}; //windlast
 //FieldVector<double, dimworld> Lone = {0, -10000, 0};//verkehrlast
 // FieldVector<double, dimworld> Ltwo = {-500, 0};
  
  for (int i=0; i < neumannNodes.size(); i++)
  {
    if (neumannNodes[i][0])
    {
      loadVector[i] = Lone;
    }
  }
 
  UMFPack<operatorType> solver(stiffnessMatrix);
  solver.setVerbosity(2);
  
  InverseOperatorResult statistics;
  solver.apply(displacementVector, loadVector, statistics);

  // generate global displacement function
  using displacementRange = FieldVector<double, dimworld>;
  auto displacementFunction = Functions::makeDiscreteGlobalBasisFunction<displacementRange> (basis, displacementVector);
  
  // generate output writer (refine output mesh depending on element order)
  SubsamplingVTKWriter<GridView> vtkWriter(gridView, refinementLevels(p));
  vtkWriter.addVertexData(displacementFunction, VTK::FieldInfo("displacement", VTK::FieldInfo::Type::vector, dimworld));
  vtkWriter.write("framework-result");

}
