// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:

#ifndef NEUMANN_BOUNDARY_HH
#define NEUMANN_BOUNDARY_HH

namespace Dune::Elastodynamics {

  template<class Function>
  class NeumannBoundaryAssembler {

    private:

      const Function& neumannFunction_;

    public:

      typedef typename Dune::BlockVector<Dune::FieldVector<double, 1>> LocalVector;
      
      NeumannBoundaryAssembler(const Function& neumannFunction)
        : neumannFunction_(neumannFunction)
      {}
      
      template<class BoundaryIterator, class LocalView>
      void assemble(const BoundaryIterator& it, LocalVector& localVector, LocalView& localView) {
        
        using Element = typename LocalView::Element;  
        const int dim = Element::dimension;
        const auto& localFE = localView.tree().child(0).finiteElement();
        int quadOrder = 2*localFE.localBasis().order();
        const auto& quadRule = QuadratureRules<double, dim-1>::rule(it->geometry().type(), quadOrder);

        localVector.resize(localView.size());
        localVector = 0.0;
             
        for(const auto& quadPoint : quadRule) {

          const auto quadPos = it->geometryInInside().global(quadPoint.position());
          const auto integrationElement = it->geometry().integrationElement(quadPoint.position());
   
          std::vector<FieldVector<double, 1>> shapefunctionValues(localFE.size());
          localFE.localBasis().evaluateFunction(quadPos, shapefunctionValues);
          
          // this part is critical ! Should Neumann Value be Force or Stress, i guess Stress, but why still wrong results
          // ... with forces it wont work at all ... values are way too small
          for( int i=0; i<localFE.size(); i++) {
            for(  int k=0; k<dim; k++) {
              auto row = localView.tree().child(k).localIndex(i);
              auto neumannValue_ = neumannFunction_(localView.index(row));
              localVector[row] += quadPoint.weight()*neumannValue_*integrationElement*shapefunctionValues[i];
            }
          }
        }
      }
  };
}     

#endif
