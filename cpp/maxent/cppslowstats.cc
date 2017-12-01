/** @file cppslowstats.cc
    @brief Calculate the expected single-site frequencies and pairwise
           correlations given the maximum entropy model parameters.
 
    This is done slowly, by summing over all the states.
**/

#include "slowstatsobject.h"

void mexFunction(int nlhs, mxArray* plhs[], int nrhs, const mxArray* prhs[])
{
  SlowStatsObject worker(nlhs, plhs, nrhs, prhs);
  worker.init();
  worker.run();
}
