// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:

#ifndef HRZ_LUMPED_MASS_ASSEMBLER_HH
#define HRZ_LUMPED_MASS_ASSEMBLER_HH

#include <dune/geometry/quadraturerules.hh>

namespace Dune::Elastodynamics {

  class HRZLumpedMassAssembler {

    private:
    
      double rho_;
        
    public:

      typedef typename Dune::Matrix<Dune::FieldMatrix<double, 1, 1>> LocalMatrix;
    	
      HRZLumpedMassAssembler(double rho)
        :rho_(rho)
      {}

      template <class LocalView>
      void assemble(LocalMatrix& localMatrix, LocalView& localView) {
        
        using Element = typename LocalView::Element;
        auto element = localView.element();
        static const int dim = LocalView::GridView::dimension;
        static const int dimworld = LocalView::GridView::dimensionworld;
        auto geometry = element.geometry();  
        
        // this is cheating, but hey, it's the same in each dimension ;)
        const auto& localFE = localView.tree().child(0).finiteElement();
        int order = 2*dim*localFE.localBasis().order();
        const auto& quadRule = QuadratureRules<double, dim>::rule(element.type(), order);
        
        localMatrix.setSize(localView.size(), localView.size());
        localMatrix = 0.0;
        
        for(const auto& quadPoint : quadRule) {
    
          const auto quadPos = quadPoint.position();
          const double integrationElement = geometry.integrationElement(quadPos);
          
          std::vector<FieldVector<double, 1>> shapefunctionValues(localFE.size());
          localFE.localBasis().evaluateFunction(quadPos, shapefunctionValues);
          
          for( int i=0; i<localFE.size(); i++) {
            for( int k=0; k<dimworld; k++) {
              auto row = localView.tree().child(k).localIndex(i);
              localMatrix[row][row] += quadPoint.weight()*rho_*integrationElement*(shapefunctionValues[i]*shapefunctionValues[i]);
            }
          }
        } 
        
        std::vector<double> weight(dimworld, 0.0);
        
        for( int i=0; i<localFE.size(); i++) {
          for( int k=0; k<dimworld; k++) {
            auto row = localView.tree().child(k).localIndex(i);
            weight[k] += localMatrix[row][row];
          }
        }
        
        double elementMass = geometry.volume()*rho_; // 1D: volume = length, 2D: volume = area
                  
        for( int i=0; i<localFE.size(); i++) {
          for( int k=0; k<dimworld; k++) {
            auto row = localView.tree().child(k).localIndex(i);
            localMatrix[row][row] *= elementMass/weight[k];
          }
        }
      }
  };
}

#endif
