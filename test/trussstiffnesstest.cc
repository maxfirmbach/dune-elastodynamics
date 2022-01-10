// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:

#include <config.h>
#include <stdio.h>

#include <dune/common/parallel/mpihelper.hh>

#include <dune/grid/uggrid.hh>
#include <dune/foamgrid/foamgrid.hh>
#include <dune/grid/io/file/gmshreader.hh>

#include <dune/istl/matrix.hh>
#include <dune/istl/bcrsmatrix.hh>
#include <dune/istl/bdmatrix.hh>
#include <dune/istl/bvector.hh>
#include <dune/istl/io.hh>

#include <dune/functions/functionspacebases/basistags.hh>
#include <dune/functions/functionspacebases/powerbasis.hh>
#include <dune/functions/functionspacebases/lagrangebasis.hh>

#include <dune/elastodynamics/assemblers/operatorassembler.hh>
#include <dune/elastodynamics/assemblers/stiffnessassembler_truss.hh>

using namespace Dune;

int main(int argc, char** argv) {

  Dune::MPIHelper& mpiHelper = Dune::MPIHelper::instance(argc, argv);
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

  // 2 node bar (simple approach)
  {
    Elastodynamics::OperatorAssembler<Basis> operatorAssembler(basis);
    double E = 5.0, A = 2.0;
    operatorType stiffnessMatrix;
    operatorAssembler.initialize(stiffnessMatrix);
    Elastodynamics::StiffnessAssemblerTrussSimple stiffnessAssembler(E, A);
    operatorAssembler.assemble(stiffnessAssembler, stiffnessMatrix, false);
  
    passed = passed and std::abs(stiffnessMatrix.frobenius_norm()-3.1234752377721211047) < 1e-10;
  }

  // 2 node bar (quadrature approach)
  {
    Elastodynamics::OperatorAssembler<Basis> operatorAssembler(basis);
    double E = 5.0, A = 2.0;
    operatorType stiffnessMatrix;
    operatorAssembler.initialize(stiffnessMatrix);
    Elastodynamics::StiffnessAssemblerTruss stiffnessAssembler(E, A);
    operatorAssembler.assemble(stiffnessAssembler, stiffnessMatrix, false);
  
    passed = passed and std::abs(stiffnessMatrix.frobenius_norm()-3.1234752377721211047) < 1e-10;
  }

  return passed ? 0 : 1;

}
