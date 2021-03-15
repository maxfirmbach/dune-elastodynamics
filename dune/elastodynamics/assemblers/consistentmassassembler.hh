// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:

#ifndef CONSISTENT_MASS_ASSEMBLER_HH
#define CONSISTENT_MASS_ASSEMBLER_HH

#include <dune/geometry/quadraturerules.hh>

namespace Dune::Elastodynamics {

  class ConsistentMassAssembler {

    private:
    
      double rho_;
        
    public:

      typedef typename Dune::Matrix<Dune::FieldMatrix<double, 1, 1>> LocalMatrix;
    	
      ConsistentMassAssembler(double rho)
        :rho_(rho)
      {}

      template <class LocalView>
      void assemble(LocalMatrix& localMatrix, LocalView& localView) {
        
        using Element = typename LocalView::Element;
        auto element = localView.element();
        const int dim = Element::dimension;
        auto geometry = element.geometry();  
        const auto& localFE = localView.tree().child(0).finiteElement();
        int order = 2*(dim*localFE.localBasis().order()-1);
        const auto& quadRule = QuadratureRules<double, dim>::rule(element.type(), order);
        
        localMatrix.setSize(localView.size(), localView.size());
        localMatrix = 0.0;
             
        for(const auto& quadPoint : quadRule) {
    
          const auto quadPos = quadPoint.position();
          const double integrationElement = geometry.integrationElement(quadPos);
          
          std::vector<FieldVector<double, 1>> shapefunctionValues(localFE.size());
          localFE.localBasis().evaluateFunction(quadPos, shapefunctionValues);
          
          for( int i=0; i<localFE.size(); i++) {
            for( int j=0; j<localFE.size(); j++) {
              for( int k=0; k<dim; k++) {
                auto row = localView.tree().child(k).localIndex(i);
                auto col = localView.tree().child(k).localIndex(j);
                localMatrix[row][col] += quadPoint.weight()*rho_*integrationElement*(shapefunctionValues[i]*shapefunctionValues[j]);
              }
            }
          }
        }
      }
  };
}

#endif
