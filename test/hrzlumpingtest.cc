// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:

#include <config.h>
#include <stdio.h>

#include <dune/common/parallel/mpihelper.hh>

#include <dune/grid/uggrid.hh>
#include <dune/foamgrid/foamgrid.hh>

#include <dune/istl/matrix.hh>
#include <dune/istl/bdmatrix.hh>

#include <dune/functions/functionspacebases/basistags.hh>
#include <dune/functions/functionspacebases/powerbasis.hh>
#include <dune/functions/functionspacebases/lagrangebasis.hh>

#include <dune/elastodynamics/assemblers/operatorassembler.hh>
#include <dune/elastodynamics/assemblers/hrzlumpedmassassembler.hh>

// test for the correct implementation of the hrzmasslumping
// for different element types and basis function orders

using namespace Dune;

struct HRZLumpingTestSuit {

  template <class Basis, class Matrix>
  void assemble(const Basis& basis, Matrix& lumpedMatrix) {

    double rho = 1.0;
    lumpedMatrix = 0.0;
    Elastodynamics::OperatorAssembler<Basis> operatorAssembler(basis);
    Elastodynamics::HRZLumpedMassAssembler hrzlumpedAssembler(rho);
    operatorAssembler.assemble(hrzlumpedAssembler, lumpedMatrix, true);

  }
  
  template <class Vector, class Matrix>
  bool check(Vector& solution, Matrix& lumpedMatrix) {
  
    bool passed = true;

    for( int i=0; i<solution.N(); i++) {
      for( int k=0; k<solution[1].N(); k++) {
        passed = passed and std::abs(lumpedMatrix[i][i][k][k]-solution[i][k]) < 1e-10;
      }
    }
    
    return passed;
  }
  
};

int main(int argc, char** argv) {

  bool passed = true;
  Dune::MPIHelper& mpiHelper = Dune::MPIHelper::instance(argc, argv);
  using namespace Functions::BasisBuilder;

  // test hrzlumping for bar elements of different order
  {
    HRZLumpingTestSuit barTest;
    
    using Grid = FoamGrid<1, 2>;
    using GridView = Grid::LeafGridView;
  
    GridFactory<Grid> factory;
    factory.insertVertex({0, 0});
    factory.insertVertex({1, 0});
    factory.insertElement(GeometryTypes::simplex(1), {0, 1}); 
    std::shared_ptr<Grid> grid(factory.createGrid());
    auto gridView = grid->leafGridView();
    
    using operatorType = BDMatrix<FieldMatrix<double, 2, 2>>;
    using blockVector  = BlockVector<FieldVector<double, 2>>;
      
    {
      // 2 node bar in 2 dim world
      auto basis = makeBasis(gridView, power<2>(lagrange<1>()));
      using Basis = decltype(basis);
    
      std::cout << "Test: HRZ lumping for linear bar" << std::endl;
      operatorType lumpedMatrix(basis.size());
      barTest.assemble<Basis>(basis, lumpedMatrix);
      blockVector solution = {{1.0, 1.0}, {1.0, 1.0}};
      solution *= 1.0/2.0;
      passed = passed and barTest.check(solution, lumpedMatrix);
    }
  }

  // test hrzlumping for simplex elements of different order
  {
    HRZLumpingTestSuit simplexTest;
    
    using Grid = UGGrid<2>;
    using GridView = Grid::LeafGridView;
  
    GridFactory<Grid> factory;
    factory.insertVertex({0, 0});
    factory.insertVertex({1, 0});
    factory.insertVertex({0, 1});
    factory.insertElement(GeometryTypes::simplex(2), {0, 1, 2}); 
    std::shared_ptr<Grid> grid(factory.createGrid());
    auto gridView = grid->leafGridView();
    
    using operatorType = BDMatrix<FieldMatrix<double, 2, 2>>;
    using blockVector  = BlockVector<FieldVector<double, 2>>;
    
    {
      // 3 node simplex
      auto basis = makeBasis(gridView, power<2>(lagrange<1>()));
      using Basis = decltype(basis);
    
      std::cout << "Test: HRZ lumping for linear simplex" << std::endl;
      operatorType lumpedMatrix(basis.size());
      simplexTest.assemble<Basis>(basis, lumpedMatrix);
      blockVector solution = {{1.0, 1.0}, {1.0, 1.0}, {1.0, 1.0}};
      solution *= 1.0/6.0;
      passed = passed and simplexTest.check(solution, lumpedMatrix);
    }
  
    {
      // 6 node simplex
      auto basis = makeBasis(gridView, power<2>(lagrange<2>()));
      using Basis = decltype(basis);
      
      std::cout << "Test: HRZ lumping for quadratic simplex" << std::endl;
      operatorType lumpedMatrix(basis.size());
      simplexTest.assemble<Basis>(basis, lumpedMatrix);
      blockVector solution = {{3.0, 3.0}, {3.0, 3.0}, {3.0, 3.0}, 
                              {16.0, 16.0}, {16.0, 16.0}, {16.0, 16.0}};
      solution *= 1.0/114.0;
      passed = passed and simplexTest.check(solution, lumpedMatrix);
    }
  
    {
      // 10 node simplex
      auto basis = makeBasis(gridView, power<2>(lagrange<3>()));
      using Basis = decltype(basis);
      
      std::cout << "Test: HRZ lumping for cubic simplex" << std::endl;
      operatorType lumpedMatrix(basis.size());
      simplexTest.assemble<Basis>(basis, lumpedMatrix);
      blockVector solution =  {{19.0, 19.0}, {19.0, 19.0}, {19.0, 19.0},
                               {135.0, 135.0}, {135.0, 135.0}, {135.0, 135.0},
                               {135.0, 135.0}, {135.0, 135.0}, {135.0, 135.0},
                               {486.0, 486.0}};
      solution *= 1.0/2706.0;
      passed = passed and simplexTest.check(solution, lumpedMatrix);
    }
  }
  
  // test hrzlumping for cube elements of different order
  {
    HRZLumpingTestSuit cubeTest;
    
    using Grid = UGGrid<2>;
    using GridView = Grid::LeafGridView;
  
    GridFactory<Grid> factory;
    factory.insertVertex({0, 0});
    factory.insertVertex({1, 0});
    factory.insertVertex({0, 1});
    factory.insertVertex({1, 1});
    factory.insertElement(GeometryTypes::cube(2), {0, 1, 2, 3}); 
    std::shared_ptr<Grid> grid(factory.createGrid());
    auto gridView = grid->leafGridView();
    
    using operatorType = BDMatrix<FieldMatrix<double, 2, 2>>;
    using blockVector  = BlockVector<FieldVector<double, 2>>;
  
    {
      // 4 node cube
      auto basis = makeBasis(gridView, power<2>(lagrange<1>()));
      using Basis = decltype(basis);
      
      std::cout << "Test: HRZ lumping for linear cube" << std::endl;
      operatorType lumpedMatrix(basis.size());
      cubeTest.assemble<Basis>(basis, lumpedMatrix);
      blockVector solution = {{1.0, 1.0}, {1.0, 1.0}, {1.0, 1.0}, {1.0, 1.0}};
      solution *= 1.0/4.0;
      passed = passed and cubeTest.check(solution, lumpedMatrix);
    }
    
    {
      // 9 node cube
      auto basis = makeBasis(gridView, power<2>(lagrange<2>()));
      using Basis = decltype(basis);
      
      std::cout << "Test: HRZ lumping for quadratic cube" << std::endl;
      operatorType lumpedMatrix(basis.size());
      cubeTest.assemble<Basis>(basis, lumpedMatrix);
      blockVector solution = {{1.0, 1.0}, {1.0, 1.0}, {1.0, 1.0}, {1.0, 1.0},
                              {4.0, 4.0}, {4.0, 4.0}, {4.0, 4.0}, {4.0, 4.0},
                              {16.0, 16.0}};
      solution *= 1.0/36.0;
      passed = passed and cubeTest.check(solution, lumpedMatrix);
    }
  }
  
  return passed ? 0 : 1;

}
