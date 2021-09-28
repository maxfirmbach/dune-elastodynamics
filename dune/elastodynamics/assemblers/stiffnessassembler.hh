// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:

#ifndef STIFFNESS_ASSEMBLER_HH
#define STIFFNESS_ASSEMBLER_HH

#include <dune/elastodynamics/assemblers/hooketensor.hh>
#include <dune/elastodynamics/assemblers/symmetrictensor.hh>
#include <dune/geometry/quadraturerules.hh>

namespace Dune::Elastodynamics {

  class StiffnessAssembler {

    private:

      double E_, nu_;
      
      template <class DeformationGradient, class Strain>
      void computeStrain(DeformationGradient& gradient, Strain& strain) {
        for (int i=0; i<gradient.N() ; ++i) {
          strain(i,i) = gradient[i][i];
          for (int j=i+1; j<gradient.M(); ++j) {
            strain(i,j) = 0.5*(gradient[i][j] + gradient[j][i]);
          }
        }
      }
      
    public:
    
      typedef typename Dune::Matrix<Dune::FieldMatrix<double, 1, 1>> LocalMatrix;
				
      StiffnessAssembler(double E, double nu)
        : E_(E), nu_(nu)
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
          const auto invJacobian = geometry.jacobianInverseTransposed(quadPos);
          
          std::vector<FieldMatrix<double, 1, dim>> referenceGradients(localFE.size());
          localFE.localBasis().evaluateJacobian(quadPos, referenceGradients);

          std::vector<FieldVector<double, dim>> gradients(localFE.size());
          for( int i=0; i<gradients.size(); i++) {
		    invJacobian.mv(referenceGradients[i][0], gradients[i]);
		  }

		  std::vector<std::array<SymmetricTensor<dim>, dim>> strain(localFE.size());  		  
		  for( int i=0; i<localFE.size(); i++) {
		    for( int k=0; k<dim; k++ ) {	      
		      Dune::FieldMatrix<double, dim, dim> deformationGradient(0);
		      deformationGradient[k] = gradients[i];
		      computeStrain(deformationGradient, strain[i][k]);
		    }            
          }

          Dune::Elastodynamics::HookeTensor<dim> hookeTensor(E_, nu_);

          for( int i=0; i<localFE.size(); i++) {
            for( int k=0; k<dim; k++) {
              auto row = localView.tree().child(k).localIndex(i);  
              SymmetricTensor<dim> stress;
              
              hookeTensor.C.mv(strain[i][k], stress);  // maybe change to lambda*tr(e)I+2*nu*e

              for (int j=0; j<localFE.size(); j++) {
                for( int l=0; l<dim; l++) {
                  auto col = localView.tree().child(l).localIndex(j);                
                  localMatrix[row][col] += quadPoint.weight()*integrationElement*(stress*strain[j][l]);
                }
              }
            }
          }       
        }
      }
  };
}

#endif
