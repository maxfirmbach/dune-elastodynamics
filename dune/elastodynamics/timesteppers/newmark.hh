// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:

#ifndef NEWMARK_HH
#define NEWMARK_HH

#include "coefficients.hh"
#include "timestepcontroller.hh"

#include <dune/elastodynamics/preconditioners/distributedscalarproduct.hh>
#include <dune/elastodynamics/preconditioners/distributedpreconditioner.hh>

namespace Dune {

  template <typename GridView, typename MatrixType, typename VectorType>	
  class Newmark {
  
    using Jacobi          = SeqJac<MatrixType, VectorType, VectorType>;
    using Preconditioner  = DistributedPreconditioner<GridView, Jacobi>;
    using ScalarProduct   = DistributedScalarProduct<GridView, VectorType>;
    using LinearOperator  = MatrixAdapter<MatrixType, VectorType, VectorType>;
	
    private:
	
	  TimeStepController fixed_;
	  double dt_;
	
	  MatrixType efficient_mass_, mass_, stiffness_;
	  const GridView& gridView_;
	
	  double beta_, gamma_;

    public:
	  
      Newmark(GridView& gridView,
              MatrixType& mass,
              MatrixType& stiffness,
              NewmarkCoefficients& coefficients,
	          TimeStepController& fixed)
	  : gridView_(gridView)
	  , mass_(mass)
	  , stiffness_(stiffness)
	  , beta_(coefficients.beta())
	  , gamma_(coefficients.gamma())
      , fixed_(fixed)
	  {}

                
      void initialize(VectorType& acceleration,
	                  VectorType load) // we only want a copy and not work on the memory here!
	  {
        // initial value calculation for acceleration
        Jacobi         jacobi(mass_, 1, 1.0);      
        Preconditioner preconditioner(gridView_, jacobi);
        ScalarProduct  scalarProduct(gridView_);
        LinearOperator linearOperator(mass_);
  
        CGSolver<VectorType> solver(linearOperator, scalarProduct, preconditioner, 1e-8, 15000,
                                    (gridView_.comm().rank()==0) ? 1 : 0);
  
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
        
        Jacobi         jacobi(efficient_mass_, 1, 1.0);      
        Preconditioner preconditioner(gridView_, jacobi);
        ScalarProduct  scalarProduct(gridView_);
        LinearOperator linearOperator(efficient_mass_);
  
        CGSolver<VectorType> solver(linearOperator, scalarProduct, preconditioner, 1e-8, 15000,
                                    (gridView_.comm().rank()==0) ? 1 : 0);
  
        InverseOperatorResult statistics;   
        solver.apply(acceleration, load, statistics);
              
        // corrector
        velocity.axpy(gamma_*dt_, acceleration);
        displacement.axpy(beta_*dt_*dt_, acceleration);
      }
  };
}

#endif
