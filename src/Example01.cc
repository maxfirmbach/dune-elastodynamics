// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:

#ifdef HAVE_CONFIG_H
# include "config.h"
#endif

#include <vector>
#include <string>
#include <stdio.h>
#include <math.h>
#include <iostream>

#include <dune/common/parallel/mpihelper.hh>

#include <dune/common/function.hh>

#include <dune/grid/uggrid.hh>
#include <dune/grid/io/file/gmshreader.hh>
#include <dune/grid/io/file/vtk/vtksequencewriter.hh>
#include <dune/grid/io/file/vtk/subsamplingvtkwriter.hh>

#include <dune/istl/matrix.hh>
#include <dune/istl/bcrsmatrix.hh>
#include <dune/istl/bvector.hh>
#include <dune/istl/eigenvalue/poweriteration.hh>
#include <dune/istl/umfpack.hh>

#include <dune/functions/functionspacebases/basistags.hh>
#include <dune/functions/functionspacebases/powerbasis.hh>
#include <dune/functions/functionspacebases/lagrangebasis.hh>

#include <dune/elastodynamics/assemblers/operatorassembler.hh>
#include <dune/elastodynamics/assemblers/stiffnessassembler.hh>

using namespace Dune;

const int dim = 2;
const int p = 1;

int main(int argc, char** argv)
{
  bool passed = true;
  Dune::MPIHelper& mpiHelper = Dune::MPIHelper::instance(argc, argv);;
  
  using operatorType = BCRSMatrix<FieldMatrix<double, dim, dim>>;
  using blockVector  = BlockVector<FieldVector<double, dim>>;
  
  // get grid from .msh file
  using Grid = UGGrid<dim>;
  using GridView = Grid::LeafGridView;

  auto mesh = "triangle.msh";
  std::vector<int> materialEntity, boundaryEntity;
  GridFactory<Grid> factory;
  GmshReader<Grid>::read(factory, mesh, boundaryEntity, materialEntity, true, true);
  shared_ptr<Grid> grid(factory.createGrid());    
  auto gridView = grid->leafGridView();
  
  using namespace Functions::BasisBuilder;
  auto basis = makeBasis(gridView, power<dim>(lagrange<p>()));
  using Basis = decltype(basis);
  
  Elastodynamics::OperatorAssembler<Basis> operatorAssembler(basis);
  
  double E = 25.0, nu = 0.2;
  operatorType stiffnessMatrix;
  operatorAssembler.initialize(stiffnessMatrix);
  Elastodynamics::StiffnessAssembler stiffnessAssembler(E, nu);
  operatorAssembler.assemble(stiffnessAssembler, stiffnessMatrix, false);
  
  blockVector x(basis.size());
  double lambda;
  PowerIteration_Algorithms<operatorType, blockVector> eigenvalueAlgorithm(stiffnessMatrix, 1000, 1);
  auto iterationMatrix = eigenvalueAlgorithm.getIterationMatrix();
  UMFPack<operatorType> solver(iterationMatrix);
  
  x = 1.0;
  eigenvalueAlgorithm.applyPowerIteration(1e-10, x, lambda);
  passed = passed and std::abs(lambda-69.664793948382650) < 1e-10;
  
  x = 1.0;
  eigenvalueAlgorithm.applyInverseIteration(29.78, 1e-10, solver, x, lambda);
  passed = passed and std::abs(lambda-30.0) < 1e-10;
  
  
  double solution[6] = {69.664793948382650, 30.0, 10.335206051617348, 0.0, 0.0, 0.0};
  

  return passed ? 0 : 1;

}
