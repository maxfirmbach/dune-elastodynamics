// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:

/*
This file is taken from the dune-fufem module!
Title: symmetrictensor.hh
Authors: dune-fufem team dune-fufem@lists.spline.inf.fu-berlin.de
Date: 2021
Version: -
Availability: https://git.imp.fu-berlin.de/agnumpde/dune-fufem/-/blob/master/dune/fufem/quadraturerules/lumpingquadraturerule.hh
*/

#ifndef LUMPINGQUADRATURERULE_HH
#define LUMPINGQUADRATURERULE_HH

#include <dune/geometry/referenceelements.hh>
#include <dune/geometry/quadraturerules.hh>

template<class ct, int dim>
class LumpingQuadratureRule : public Dune::QuadratureRule<ct, dim> {

  using Base = Dune::QuadratureRule<ct, dim>;
  using QuadPoint = Dune::QuadraturePoint<ct, dim>;

  public:

    LumpingQuadratureRule(const Dune::GeometryType& gt) : Base(gt, 1) {
  
      auto refElement = Dune::ReferenceElements<ct, dim>::general(this->type());
      int size = refElement.size(dim);
            
      ct weight = refElement.volume()/size;
      this->reserve(size);
      for(int i=0; i<size; ++i) {
        this->push_back(QuadPoint(refElement.position(i,dim), weight));
      }
    }
    
};

#endif
