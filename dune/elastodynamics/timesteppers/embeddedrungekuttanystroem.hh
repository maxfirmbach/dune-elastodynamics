// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:

#ifndef EMBEDDED_RUNGE_KUTTA_NYSTROEM_HH
#define EMBEDDED_RUNGE_KUTTA_NYSTROEM_HH

#include "coefficients.hh"
#include "timestepcontroller.hh"

namespace Dune {
  
  template <typename MatrixType, typename VectorType>	
  class EmbeddedRungeKuttaNystroem {
  
    private:
	  
      AdaptiveStepController *adaptive_;
      double dt_;
	
	  MatrixType lumpedmass_, stiffness_;
	
	  int stages_, order_;
	  Dune::Matrix<Dune::FieldMatrix<double, 1, 1>> A_;
	  Dune::BlockVector<Dune::FieldVector<double, 1>> b_, b_bar_, b_tilde_, b_bar_tilde_, c_;
	  Dune::BlockVector<VectorType> k;
		
    public:
	
	  EmbeddedRungeKuttaNystroem(MatrixType& lumpedmass,
	   	                         MatrixType& stiffness,
					             EmbeddedRKNCoefficients& coefficients,
					             AdaptiveStepController* adaptive)
      : lumpedmass_(lumpedmass)
	  , stiffness_(stiffness)
	  , A_(coefficients.A())
	  , b_(coefficients.b())
	  , b_bar_(coefficients.b_bar())
	  , b_tilde_(coefficients.b_tilde())
	  , b_bar_tilde_(coefficients.b_bar_tilde())
	  , c_(coefficients.c())
	  , stages_(coefficients.stages())
	  , order_(coefficients.order())
	  , adaptive_(adaptive)
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
      
        while(1)
        {
          // store values of last timestep
          VectorType displacement_ = displacement;
          VectorType displacement_tilde_ = displacement;
          VectorType velocity_ = velocity;
          VectorType velocity_tilde_ = velocity;
        
          // get fixed timestep size
	      dt_ = adaptive_->deltaT();
		
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
	      displacement_.axpy(dt_, velocity_);
	      displacement_tilde_.axpy(dt_, velocity_tilde_);
	      for(int i=0; i<stages_; i++) 
	      {
	        displacement_.axpy(dt_*dt_*b_bar_[i], k[i]);
		    displacement_tilde_.axpy(dt_*dt_*b_bar_tilde_[i], k[i]);
		    velocity_.axpy(dt_*b_[i], k[i]);
		    velocity_tilde_.axpy(dt_*b_tilde_[i], k[i]);
	      }
	      
	      // calculate error
	      displacement_ -= displacement_tilde_;
	      velocity_ -= velocity_tilde_;	   
	      double error_ = std::max(displacement_.infinity_norm(), velocity_.infinity_norm());
	       
	      // get new timestep
	      bool accepted = adaptive_->timeStepValid(dt_, error_, order_);

	      if(accepted)
	      {
	        displacement = displacement_tilde_;
	        velocity = velocity_tilde_;
	        break;
	      }	 
	    }	    
	  }		  
  };
}

#endif
