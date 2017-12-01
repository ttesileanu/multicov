/** @file cppsimmaxent.cc
    @brief Perform an MCMC simulation to draw samples from a maximum entropy
           model.

    Goals:
      - should work with arbitrary combinations of alphabets
      - should feature caching of results to avoid expensive recalculations
      - should have plenty of options for tracking the trajectory
**/

#include "simobject.h"

void mexFunction(int nlhs, mxArray* plhs[], int nrhs, const mxArray* prhs[])
{
  SimulationObject worker(nlhs, plhs, nrhs, prhs);
  worker.init();
  worker.run();
}
