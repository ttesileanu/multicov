/** @file calcobject.h
 *  @brief Defines the class that performs the work for calcseqw_cpp.
 *
 *  @author Tiberiu Tesileanu
 */
#ifndef CALCOBJECT_H_
#define CALCOBJECT_H_

#include "mexobject.h"

#include "matrix.h"

typedef Matrix<mxChar> Alignment;

class CalculateObject : public MexObject {
 public:
  /// Constructor from MEX command line.
  CalculateObject(int nlhs, mxArray* plhs[], int nrhs, const mxArray* prhs[])
    : MexObject("perturb_cpp", nlhs, plhs, nrhs, prhs), maxSeqid_(1) {}

  /// Override input checking.
  virtual void checkInputs();

  /// Override input wrappers initialization.
  virtual void initInputs();

  /// Override function run.
  virtual void run();

 private:
  /// The alignment.
  Alignment         alignment_;
  /// The threshold for calculating sequence weights.
  double            maxSeqid_;
};

#endif
