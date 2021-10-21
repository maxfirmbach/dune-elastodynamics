// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:

#ifndef NEWMARK_HH
#define NEWMARK_HH

#include "coefficients.hh"
#include "timestepcontroller.hh"

#include <dune/istl/umfpack.hh>

namespace Dune {

  template <typename MatrixType, typename VectorType>	
  class Newmark {
	
    private:
	
	  TimeStepController fixed_;
	  double dt_;
	
	  MatrixType efficient_mass_, mass_, stiffness_;	
	  double beta_, gamma_;

    public:
	  
      Newmark(MatrixType& mass,
              MatrixType& stiffness,
              NewmarkCoefficients& coefficients,
	          TimeStepController& fixed)
	  : mass_(mass)
	  , stiffness_(stiffness)
	  , beta_(coefficients.beta())
	  , gamma_(coefficients.gamma())
      , fixed_(fixed)
	  {}

                
      void initialize(VectorType& acceleration,
	                  VectorType load) // we only want a copy and not work on the memory here!
	  {
        // initial value calculation for acceleration
        UMFPack<MatrixType> solver(mass_);
        solver.setVerbosity(1);

        InverseOperatorResult statistics;   
        solver.apply(acceleration, load, statistics);
        
        // calculate efficient mass matrix
        dt_ = fixed_.deltaT();
        efficient_mass_ = mass_;
        efficient_mass_.axpy(beta_*dt_*dt_, stiffness_);
      }

	
      void step(VectorType& displacement,
                VectorType& velocity,
                VectorType& acceleration,
                VectorType load) // we only want a copy and not work on the memory here!
      {
        // get fixed timestepsize
        dt_ = fixed_.deltaT();
    
        // predictor      
        displacement.axpy(dt_, velocity);
        displacement.axpy((0.5-beta_)*dt_*dt_, acceleration);
        velocity.axpy((1.0-gamma_)*dt_, acceleration);
      
        // solve
        stiffness_.mmv(displacement, load);
        
        UMFPack<MatrixType> solver(efficient_mass_);
        solver.setVerbosity(1);
 
        InverseOperatorResult statistics;
        solver.apply(acceleration, load, statistics);
              
        // corrector
        velocity.axpy(gamma_*dt_, acceleration);
        displacement.axpy(beta_*dt_*dt_, acceleration);
      }
  };
}

#endif
