// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:

#ifndef STIFFNESS_ASSEMBLER_TRUSS_HH
#define STIFFNESS_ASSEMBLER_TRUSS_HH

#include <assert.h>

#include <dune/geometry/quadraturerules.hh>

namespace Dune::Elastodynamics {

  class StiffnessAssemblerTrussSimple {

    private:
    
      double E_, A_;

    public:

      typedef typename Dune::Matrix<Dune::FieldMatrix<double, 1, 1>> LocalMatrix;
      StiffnessAssemblerTrussSimple(double E, double A)
        : E_(E), A_(A)
      {}

      template <class LocalView>
      void assemble(LocalMatrix& localMatrix, LocalView& localView) {
        
        using Element = typename LocalView::Element;
        auto element = localView.element();
        static const int dim = LocalView::GridView::dimension;
        static const int dimworld = LocalView::GridView::dimensionworld;
        auto geometry = element.geometry();
        const auto& localFE = localView.tree().child(0).finiteElement();

        if(dimworld!=2 && localFE.localBasis().order()!=1) {
          std::cout << "Dimension has to be 2 and polynomial order 1!" << std::endl;
          exit(0);
        }

        localMatrix.setSize(localView.size(), localView.size());
        localMatrix = 0.0;
        
        // get coordinates from point 0 and 1
        auto x0 = geometry.corner(0)[0];
        auto y0 = geometry.corner(0)[1];
        auto x1 = geometry.corner(1)[0];
        auto y1 = geometry.corner(1)[1];
        
        double Theta, L;

        // calculate length of element
        L = sqrt(((x1-x0)*(x1-x0))+((y1-y0)*(y1-y0)));
        
        // calculate angle theta    
        Theta = atan((y1-y0)/(x1-x0));
        
        int x=0;
        for( int i=0; i<localFE.size(); i++) {
          for( int k=0; k<dimworld; k++) {
            auto row = localView.tree().child(k).localIndex(i);
            int y=0;
            for( int j=0; j<localFE.size(); j++) {
              for( int l=0; l<dimworld; l++) {
                auto col = localView.tree().child(l).localIndex(j);
                double value;
                if(x==0 && y==0) value =  pow(cos(Theta), 2.0);
                if(x==0 && y==1) value =  cos(Theta)*sin(Theta);
                if(x==0 && y==2) value = -pow(cos(Theta), 2.0);
                if(x==0 && y==3) value = -cos(Theta)*sin(Theta);
                if(x==1 && y==0) value =  cos(Theta)*sin(Theta);
                if(x==1 && y==1) value =  pow(sin(Theta), 2.0);
                if(x==1 && y==2) value = -cos(Theta)*sin(Theta);
                if(x==1 && y==3) value = -pow(sin(Theta), 2.0);
                if(x==2 && y==0) value = -pow(cos(Theta), 2.0);
                if(x==2 && y==1) value = -cos(Theta)*sin(Theta);
                if(x==2 && y==2) value =  pow(cos(Theta), 2.0);
                if(x==2 && y==3) value =  cos(Theta)*sin(Theta);
                if(x==3 && y==0) value = -cos(Theta)*sin(Theta);
                if(x==3 && y==1) value = -pow(sin(Theta), 2.0);
                if(x==3 && y==2) value =  cos(Theta)*sin(Theta);
                if(x==3 && y==3) value =  pow(sin(Theta), 2.0);
                localMatrix[row][col] = value;
                y++;
              }
            }
            x++;
          }        
        }
        localMatrix *= A_*E_/L;
    }
  };

  class StiffnessAssemblerTrussQuadrature {

    private:
    
      double E_, A_;

    public:

      typedef typename Dune::Matrix<Dune::FieldMatrix<double, 1, 1>> LocalMatrix;
      StiffnessAssemblerTrussQuadrature(double E, double A)
        : E_(E), A_(A)
      {}

      template <class LocalView>
      void assemble(LocalMatrix& localMatrix, LocalView& localView) {
        
        using Element = typename LocalView::Element;
        auto element = localView.element();
        static const int dim = LocalView::GridView::dimension;
        static const int dimworld = LocalView::GridView::dimensionworld;
        auto geometry = element.geometry();
        const auto& localFE = localView.tree().child(0).finiteElement();

        localMatrix.setSize(localView.size(), localView.size());
        localMatrix = 0.0;

        int order = 2*(dimworld*localFE.localBasis().order()-1);;
        const auto& quadRule = QuadratureRules<double, dim>::rule(element.type(), order);

        for(const auto& quadPoint : quadRule) {
    
          const auto quadPos = quadPoint.position();
          const double integrationElement = geometry.integrationElement(quadPos);
          const auto invJacobian = geometry.jacobianInverseTransposed(quadPos);

          std::vector<FieldMatrix<double, 1, dim>> referenceGradients(localFE.size());
          localFE.localBasis().evaluateJacobian(quadPos, referenceGradients);

          std::vector<FieldVector<double, dimworld>> gradients(localFE.size());
          for( int i=0; i<gradients.size(); i++) {
		    invJacobian.mv(referenceGradients[i][0], gradients[i]);
		  }

          for( int i=0; i<localFE.size(); i++) {
            for( int k=0; k<dimworld; k++) {
              auto row = localView.tree().child(k).localIndex(i);
              for( int j=0; j<localFE.size(); j++) {
                for( int l=0; l<dimworld; l++) {
                  auto col = localView.tree().child(l).localIndex(j);
                  localMatrix[row][col] += quadPoint.weight()*A_*E_*integrationElement*(gradients[i][k]*gradients[j][l]);
                }
              }
            }
          }
        }
      }
  };
}

#endif
