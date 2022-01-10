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
  Elastodynamics::StiffnessAssemblerTrussQuadrature stiffnessAssembler(E, A);
  operatorAssembler.assemble(stiffnessAssembler, stiffnessMatrix, false);
  
  printSparseMatrix(std::cout, stiffnessMatrix, "", "", 3, 2);
  
  blockVector loadVector(basis.size()), displacementVector(basis.size()); //basis nennt die Groesse von loadVector und displacementVector
  loadVector = 0.0, displacementVector = 0.0;
  
  // set dirichlet boundary
  auto verschiebungPredicate = [](auto coordinate) //Function erstellen
  {
    return coordinate[1] > -0.00001; //Markierung von den Auflagern
  };
  blockVector dirichletNodes(basis.size()); //bool vector erstellen
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
  auto loadPredicate = [](auto coordinate) //Function erstellen
  {
    return coordinate[1] < -8.655555; //markierung von dem Kraft am Ende wirken soll
  };
  blockVector neumannNodes(basis.size()); //bool vector erstellen
  Functions::interpolate(basis, neumannNodes, loadPredicate);
  
  FieldVector<double, dimworld> L = {0, -1732};
  
  for (int i=0; i < neumannNodes.size(); i++)
  {
    if (neumannNodes[i][0])
    {
      loadVector[i] = L;
    }
  }
  
  MatrixAdapter<operatorType, blockVector, blockVector> linearOperator(stiffnessMatrix);
  SeqILU<operatorType, blockVector, blockVector> preconditioner(stiffnessMatrix, 1.0);
  CGSolver<blockVector> cg(linearOperator,
                      preconditioner,
                      1e-12,    // residual reduction factor
                      10000,    // Maximum number of iterations
                      2);       // Verbosity of the solver
  
  InverseOperatorResult statistics;
  cg.apply(displacementVector, loadVector, statistics);

  // get maxium displacement    
  double max = 0.0;
  for( int i=0; i<basis.size(); i++) {
    for( int k=0; k<dimworld; k++) {
      max = std::max(max, std::abs(displacementVector[i][k]));
    }
  }
  std::cout << max-0.01155 << std::endl;

  // generate global displacement function
  using displacementRange = FieldVector<double, dimworld>;
  auto displacementFunction = Functions::makeDiscreteGlobalBasisFunction<displacementRange> (basis, displacementVector);
  
  // generate output writer (refine output mesh depending on element order)
  VTKWriter<GridView> vtkWriter(gridView);
  vtkWriter.addVertexData(displacementFunction, VTK::FieldInfo("displacement", VTK::FieldInfo::Type::vector, dimworld));
  vtkWriter.write("framework-result");

} 
