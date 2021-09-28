// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:

#ifndef BOUNDARY_ASSEMBLER_HH
#define BOUNDARY_ASSEMBLER_HH

#include <dune/geometry/referenceelements.hh>

namespace Dune::Elastodynamics {

  template<class Basis>
  class BoundaryAssembler {
  
    private:

      using GridView = typename Basis::GridView;
      
      const Basis& basis_;
      const BoundaryPatch<GridView>& boundaryPatch_;
      
      template <class LocalBoundaryAssemblerType, class GlobalVectorType>
	  void addEntries(LocalBoundaryAssemblerType& localAssembler, GlobalVectorType& b) {

        auto localView = basis_.localView();

        typedef typename LocalBoundaryAssemblerType::LocalVector LocalVector;
 
        for( auto it = boundaryPatch_.begin(); it != boundaryPatch_.end(); it++) {
 
          const auto& element = it->inside();
          localView.bind(element);
                   
          LocalVector localVector;
          localAssembler.assemble(it, localVector, localView);
 
          for( int i=0; i<localView.size(); i++) {
            auto row = localView.index(i);
            b[row[0]][row[1]] = localVector[i];
          }
        }
      }

    public:
    
      BoundaryAssembler(const Basis& basis, const BoundaryPatch<GridView>& boundaryPatch)
        : basis_(basis), boundaryPatch_(boundaryPatch)
      {}
      
      template<class LocalBoundaryAssemblerType, class GlobalVectorType>
      void assemble(LocalBoundaryAssemblerType& localAssembler, GlobalVectorType& b)
      {
        addEntries(localAssembler, b);
      }
      
  }; 
}     

#endif
