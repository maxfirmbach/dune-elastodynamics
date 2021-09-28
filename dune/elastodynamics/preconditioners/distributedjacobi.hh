// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:

/*
This file is taken from the Book: DUNE-The Distributed and Unified
Numerics Environment by Oliver Sander and is modified!
*/

#ifndef DISTRIBUTED_JACOBI_HH
#define DISTRIBUTED_JACOBI_HH

#include <dune/istl/preconditioner.hh>
#include <dune/istl/solvercategory.hh>
#include <dune/elastodynamics/parallel/vectordatahandle.hh>

template<class GridView, class Matrix, class Vector>
class DistributedJacobi : public Dune::Preconditioner<Vector, Vector> {

  private:

    const GridView& gridView_;
    const Matrix& matrix_;
    Vector consistentDiagonal_;
    
  public:

    DistributedJacobi(const GridView& gridView, const Matrix& matrix)
      : gridView_(gridView),
        matrix_(matrix)
    {
      Vector diagonal(matrix_.N());
      for(int i=0; i<diagonal.size(); i++) {
        for(int j=0; j<GridView::dimension; j++) {
          diagonal[i][j] = matrix_[i][i][j][j];
        }
      }
      consistentDiagonal_ = diagonal;
      VectorExchangeAdd<GridView, Vector> datahandle(gridView_, consistentDiagonal_);
      gridView_.communicate(datahandle, Dune::All_All_Interface, Dune::ForwardCommunication);
    }

    virtual void pre (Vector& x, Vector& b) {}

    virtual void apply (Vector& v, const Vector& r)
    {
      auto rConsistent = r;
      VectorExchangeAdd<GridView, Vector> datahandle(gridView_, rConsistent);
      gridView_.communicate(datahandle, Dune::InteriorBorder_InteriorBorder_Interface, Dune::ForwardCommunication);
      
      for(int i=0; i<matrix_.N(); i++) {
        for(int j=0; j<GridView::dimension; j++) {
          v[i][j] = rConsistent[i][j]/consistentDiagonal_[i][j];
        }
      }
    }

    virtual void post (Vector& x) {}

    virtual Dune::SolverCategory::Category category() const
    { return Dune::SolverCategory::sequential; }

};

#endif
