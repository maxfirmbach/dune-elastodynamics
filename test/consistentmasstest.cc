// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:

#include <config.h>
#include <stdio.h>

#include <dune/common/parallel/mpihelper.hh>

#include <dune/grid/uggrid.hh>
#include <dune/foamgrid/foamgrid.hh>

#include <dune/istl/matrix.hh>
#include <dune/istl/bcrsmatrix.hh>

#include <dune/functions/functionspacebases/basistags.hh>
#include <dune/functions/functionspacebases/powerbasis.hh>
#include <dune/functions/functionspacebases/lagrangebasis.hh>

#include <dune/elastodynamics/assemblers/operatorassembler.hh>
#include <dune/elastodynamics/assemblers/consistentmassassembler.hh>

// test for the correct implementation of the hrzmasslumping
// for different element types and basis function orders

using namespace Dune;

struct ConsistentMassTestSuit {

  template <class Basis, class Matrix>
  void assemble(const Basis& basis, Matrix& matrix) {

    double rho = 1.0;
    Elastodynamics::OperatorAssembler<Basis> operatorAssembler(basis);
    operatorAssembler.initialize(matrix);
    Elastodynamics::ConsistentMassAssembler consistentmassAssembler(rho);
    operatorAssembler.assemble(consistentmassAssembler, matrix, false);

  }
    
};

int main(int argc, char** argv) {

  bool passed = true;
  Dune::MPIHelper& mpiHelper = Dune::MPIHelper::instance(argc, argv);
  using namespace Functions::BasisBuilder;

  // test consistent mass for bar element
  {
    ConsistentMassTestSuit barTest;
    
    using Grid = FoamGrid<1, 2>;
    using GridView = Grid::LeafGridView;
  
    GridFactory<Grid> factory;
    factory.insertVertex({0, 0});
    factory.insertVertex({1, 0});
    factory.insertElement(GeometryTypes::simplex(1), {0, 1}); 
    std::shared_ptr<Grid> grid(factory.createGrid());
    auto gridView = grid->leafGridView();
    
    using operatorType = BCRSMatrix<FieldMatrix<double, 2, 2>>;
    using blockVector  = BlockVector<FieldVector<double, 2>>;
  
    {
      // 2 node bar
      auto basis = makeBasis(gridView, power<2>(lagrange<1>()));
      using Basis = decltype(basis);
    
      std::cout << "Test: Consistent mass for linear bar" << std::endl;
      operatorType matrix;
      barTest.assemble<Basis>(basis, matrix);
      passed = passed and std::abs(matrix.frobenius_norm()-0.74535599249993) < 1e-10;
    }
  }

  // test consistent mass for simplex element
  {
    ConsistentMassTestSuit simplexTest;
    
    using Grid = UGGrid<2>;
    using GridView = Grid::LeafGridView;
  
    GridFactory<Grid> factory;
    factory.insertVertex({0, 0});
    factory.insertVertex({1, 0});
    factory.insertVertex({0, 1});
    factory.insertElement(GeometryTypes::simplex(2), {0, 1, 2}); 
    std::shared_ptr<Grid> grid(factory.createGrid());
    auto gridView = grid->leafGridView();
    
    using operatorType = BCRSMatrix<FieldMatrix<double, 2, 2>>;
    using blockVector  = BlockVector<FieldVector<double, 2>>;
    
    {
      // 3 node simplex
      auto basis = makeBasis(gridView, power<2>(lagrange<1>()));
      using Basis = decltype(basis);
    
      std::cout << "Test: Consistent mass for linear simplex" << std::endl;
      operatorType matrix;
      simplexTest.assemble<Basis>(basis, matrix);
      passed = passed and std::abs(matrix.frobenius_norm()-0.250000000000000) < 1e-10;
    
    }
    
    {
      // 6 node simplex
      auto basis = makeBasis(gridView, power<2>(lagrange<2>()));
      using Basis = decltype(basis);
      
      std::cout << "Test: Consistent mass for quadratic simplex" << std::endl;
      operatorType matrix;
      simplexTest.assemble<Basis>(basis, matrix);
      passed = passed and std::abs(matrix.frobenius_norm()-0.272675359817956) < 1e-10;
    }
  }
  
  // test consistent mass for cube element
  // ...
   
  return passed ? 0 : 1;

}
