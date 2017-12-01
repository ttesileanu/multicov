/** @file cppevalmaxent.cc
    @brief Calculate maximum entropy energies for a set of sequences.
**/

#include "evalobject.h"

void mexFunction(int nlhs, mxArray* plhs[], int nrhs, const mxArray* prhs[])
{
  EvaluationObject worker(nlhs, plhs, nrhs, prhs);
  worker.init();
  worker.run();
}
