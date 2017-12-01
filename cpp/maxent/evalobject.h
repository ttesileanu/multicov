/** @file evalobject.h
 *  @brief Defines the class that performs the work for cppevalmaxent.
 */
#ifndef POTTS_EVALOBJECT_H_
#define POTTS_EVALOBJECT_H_

#include "mexobject.h"

#include "evaluator.h"
#include "matrix.h"

/// A sequence alignment.
typedef Matrix<mxChar> Alignment;

/// The object that does all the evaluation work.
class EvaluationObject : public MexObject {
 public:
  /// Constructor from MEX command line.
  EvaluationObject(int nlhs, mxArray* plhs[], int nrhs, const mxArray* prhs[])
    : MexObject("cppevalmaxent", nlhs, plhs, nrhs, prhs) {}

  /// Override input checking.
  virtual void checkInputs();

  /// Override input wrappers initialization.
  virtual void initInputs();

  /// Override function run.
  virtual void run();

 private:
  /// The alignment.
  Alignment           alignment_;
  /// The alphabets.
  std::vector<std::string>
                      alphabets_;
  /// The alphabet widths.
  std::vector<size_t> alphawidths_;
  /// Couplings.
  Couplings           couplings_;
};

#endif
