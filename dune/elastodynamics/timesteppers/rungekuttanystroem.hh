// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:

#ifndef RUNGE_KUTTA_NYSTROEM_HH
#define RUNGE_KUTTA_NYSTROEM_HH

#include "coefficients.hh"
#include "timestepcontroller.hh"

namespace Dune {
  
  template <typename MatrixType, typename VectorType>	
  class RungeKuttaNystroem {
  
    private:
	  
      TimeStepController fixed_;
      double dt_;
	
	  MatrixType lumpedmass_, stiffness_;
	
	  int stages_, order_;
	  Dune::Matrix<Dune::FieldMatrix<double, 1, 1>> A_;
	  Dune::BlockVector<Dune::FieldVector<double, 1>> b_, b_bar_, c_;
	  Dune::BlockVector<VectorType> k;
		
    public:
	
	  RungeKuttaNystroem(MatrixType& lumpedmass,
	   	                 MatrixType& stiffness,
					     RKNCoefficients& coefficients,
					     TimeStepController& fixed)
      : lumpedmass_(lumpedmass)
	  , stiffness_(stiffness)
	  , A_(coefficients.A())
	  , b_(coefficients.b())
	  , b_bar_(coefficients.b_bar())
	  , c_(coefficients.c())
	  , stages_(coefficients.stages())
	  , order_(coefficients.order())
	  , fixed_(fixed)
      {
        // set storage for stages
        k.resize(stages_);
      }
	
      void initialize(VectorType& load) 
      {
        // initialize stages	                    
        for(int i=0; i<stages_; i++) {
	      k[i].resize(load.size());
		  k[i] = 0.0;
	    }
      }
	
	  void step(VectorType& displacement,
                VectorType& velocity,
                VectorType& acceleration,
                VectorType& load)
      {
        // get fixed timestep size
	    dt_ = fixed_.deltaT();
		
	    // calculate function evaluation vectors k
	    for(int i=0; i<stages_; i++) 
	    {
	      k[i] = displacement;
		  k[i].axpy(dt_*c_[i], velocity);
		  for (int j=0; j<i ; j++) {
		    k[i].axpy(dt_*dt_*A_[i][j], k[j]);
		  }
	
	      // function evaluation			
		  VectorType loadupdate;
	      loadupdate = load;
		  stiffness_.mmv(k[i], loadupdate);
		  lumpedmass_.mv(loadupdate, k[i]);
	    }
	    
	    // perform update
	    displacement.axpy(dt_, velocity);
	    for(int i=0; i<stages_; i++) 
	    {
	      displacement.axpy(dt_*dt_*b_bar_[i], k[i]);
		  velocity.axpy(dt_*b_[i], k[i]);
	    }	 	    
	  }	
  };
}

#endif
