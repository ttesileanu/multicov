#include "matlab_random.h"

void MatlabRandom::set_size_() const
{
  if (rand_arg_) mxDestroyArray(rand_arg_);

  rand_arg_ = mxCreateDoubleMatrix(1, 2, mxREAL);
  mxGetPr(rand_arg_)[0] = buffer_size_;
  mxGetPr(rand_arg_)[1] = 1;
}

void MatlabRandom::replenish_() const
{
  if (buffer_) mxDestroyArray(buffer_);

  // call Matlab's 'rand'
  mexCallMATLAB(1, &buffer_, 1, &rand_arg_, "rand");

  // reset the index
  idx_ = 0;
}

