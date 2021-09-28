// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:

/*
This file is taken from the Book: DUNE-The Distributed and Unified
Numerics Environment by Oliver Sander and is modified!
*/

#ifndef DISTRIBUTED_SCALAR_PRODCUT_HH
#define DISTRIBUTED_SCALAR_PRODCUT_HH

#include <dune/istl/scalarproducts.hh>
#include <dune/istl/solvercategory.hh>
#include <dune/elastodynamics/parallel/vectordatahandle.hh>

template<class GridView, class Vector>
class DistributedScalarProduct : public Dune::ScalarProduct<Vector> {
    
  using typename Dune::ScalarProduct<Vector>::field_type;
  using typename Dune::ScalarProduct<Vector>::real_type;
    
  private:
    
    const GridView& gridView_;
    
  public:
    
  DistributedScalarProduct(const GridView& gridView)
    : gridView_(gridView)
  {}
  
  virtual field_type dot(const Vector& x, const Vector& y) const override
  { return gridView_.comm().sum(x.dot(y)); }
      
  virtual real_type norm(const Vector& x) const override
  {
    auto xConsistent = x;
    VectorExchangeAdd<GridView, Vector> datahandle(gridView_, xConsistent);
    gridView_.communicate(datahandle, Dune::InteriorBorder_InteriorBorder_Interface, Dune::ForwardCommunication);
      
    auto localNorm2 = x.dot(xConsistent);
    return std::sqrt(gridView_.comm().sum(localNorm2));
  }
      
  virtual Dune::SolverCategory::Category category() const override
  { return Dune::SolverCategory::sequential; }

};

#endif
