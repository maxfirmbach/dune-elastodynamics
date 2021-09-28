// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:

/*
This file is taken from the dune-fufem module, but modified in some parts!
Title: symmetrictensor.hh
Authors: dune-fufem team dune-fufem@lists.spline.inf.fu-berlin.de
Date: 2021
Version: -
Availability: https://git.imp.fu-berlin.de/agnumpde/dune-fufem/-/blob/master/dune/fufem/symmetrictensor.hh
*/

#ifndef HOOKETENSOR_HH
#define HOOKETENSOR_HH

namespace Dune::Elastodynamics {

  template <int dim>
  class HookeTensor {
  
    HookeTensor(double E, double nu) {
      std::cout << "You broke the universe!" << std::endl;
    }
  };

  template <>
  class HookeTensor<2> {
  
    private:
      static const int dim = 2;
      
    public: 
    
      Dune::FieldMatrix<double, (dim+1)*dim/2, (dim+1)*dim/2> C = 0.0;
    
      HookeTensor(double E, double nu) {
        
        // Here be careful which to choose, i guess      
        // plane strain
        /*
        C[0][0] = 1.0-nu; C[0][1] = nu;
        C[1][0] = nu;     C[1][1] = 1.0-nu;
        
        C[2][2] = 1.0-2.0*nu;
        
        C *= E/((1.0+nu)*(1.0-2.0*nu));
        */
        
        // plane stress
        C[0][0] = 1.0; C[0][1] = nu;
        C[1][0] = nu;  C[1][1] = 1.0;
        
        C[2][2] = 1.0-nu;
        
        C *= E/(1.0-nu*nu);
      }
  };
     
  template <>
  class HookeTensor<3> {
  
    private:
      static const int dim = 3;
      
    public:
    
      Dune::FieldMatrix<double, (dim+1)*dim/2, (dim+1)*dim/2> C = 0.0;
    
      HookeTensor<3>(double E, double nu) {
      
        C[0][0] = 1.0-nu; C[0][1] = nu;     C[0][2] = nu;
        C[1][0] = nu;     C[1][1] = 1.0-nu; C[1][2] = nu;
        C[2][0] = nu;     C[2][1] = nu;     C[2][2] = 1.0-nu;
        
        C[3][3] = 1.0-2.0*nu;
        C[4][4] = 1.0-2.0*nu;
        C[5][5] = 1.0-2.0*nu;
        
        C *= E/((1.0+nu)*(1.0-2.0*nu));
      
      }
  };
}

#endif
