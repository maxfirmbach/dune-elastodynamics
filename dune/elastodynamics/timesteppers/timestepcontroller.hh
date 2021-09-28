#ifndef TIME_STEP_CONTROLLER_HH
#define TIME_STEP_CONTROLLER_HH

#include <math.h>

namespace Dune {
	class TimeStepController 
	{
		private:
		
		public:
			
			double deltaT() {
				return dt_;
			}

		protected:
			double dt_;
			double time_;
			
			TimeStepController(double time, double dt)
			: time_(time)
			, dt_(dt)
			{}
	};


	class FixedStepController : public TimeStepController
	{
		private:
		
		public:
			FixedStepController(double time, double dt)
			: TimeStepController(time, dt)
			{}
		
	};

	class AdaptiveStepController : public TimeStepController
	{
		private:
			double tol_;

		public:
			AdaptiveStepController(double time, double dt, double tol)
			: TimeStepController(time, dt)
			, tol_(tol)
			{}

			bool timeStepValid(double dt, double error, int p) {
				// check if error is bigger than given tolerance
				dt_ = 0.9*dt*std::pow(tol_/error, 1.0/(p+1.0));
				if(error > tol_) {
					// timestep not valid, calculate new dt and recalculate
					return false;
				} else {
					// error is below tolerance, timestep accepted with new dt
					return true;
				}
			}
	};
}

#endif
