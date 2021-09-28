// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:

#ifndef DIRICHLET_BOUNDARY_HH
#define DIRICHLET_BOUNDARY_HH

namespace Dune::Elastodynamics {

  template<class Function>
  class DirichletBoundaryAssembler {

    private:
    
      double value_;

    public:

      typedef typename Dune::BlockVector<Dune::FieldVector<double, 1>> LocalVector;
      
      DirichletBoundaryAssembler(double value)
        : value_(value)
      {}
      
      template<class BoundaryIterator, class LocalView>
      void assemble(const BoundaryIterator& it, LocalVector& localVector, LocalView& localView) {

        localVector.resize(localView.size());
        localVector = value_;
             
      }
  };
}     

#endif
