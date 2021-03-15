// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:

#ifndef NEWMARK_HH
#define NEWMARK_HH

#include "coefficients.hh"
#include "timestepcontroller.hh"

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
	                  VectorType& load)
	  {
        // initial value calculation for acceleration
        MatrixAdapter<MatrixType, VectorType, VectorType> massOperator(mass_);
        SeqILU0<MatrixType, VectorType, VectorType> preconditioner(mass_, 1.0);
        CGSolver<VectorType> cg(massOperator, preconditioner, 1e-10, 15000, 1);
        InverseOperatorResult statistics;
        cg.apply(acceleration, load, statistics);
    
        // calculate efficient mass matrix
        dt_ = fixed_.deltaT();
        efficient_mass_ = mass_;
        efficient_mass_.axpy(beta_*dt_*dt_, stiffness_);
      }

	
      void step(VectorType& displacement,
                VectorType& velocity,
                VectorType& acceleration,
                VectorType& load)
      {
        // get fixed timestepsize
        dt_ = fixed_.deltaT();
    
        // predictor      
        displacement.axpy(dt_, velocity);
        displacement.axpy((0.5-beta_)*dt_*dt_, acceleration);   
        velocity.axpy((1.0-gamma_)*dt_, acceleration);
      
        // solve
        stiffness_.mmv(displacement, load);
      
        MatrixAdapter<MatrixType, VectorType, VectorType> efficient_massOperator(efficient_mass_);
        SeqILU0<MatrixType, VectorType, VectorType> preconditioner(efficient_mass_, 1.0);
        CGSolver<VectorType> cg(efficient_massOperator, preconditioner, 1e-8, 15000, 1);
        InverseOperatorResult statistics;
        cg.apply(acceleration, load, statistics);
      
        // corrector
        velocity.axpy(gamma_*dt_, acceleration);
        displacement.axpy(beta_*dt_*dt_, acceleration);
      }
      
  };
}

#endif
