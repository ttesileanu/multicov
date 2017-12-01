/** @file mexobject.h
 *  @brief Defines an interface for wrapping a MEX job inside a class.
 *
 *  @author Tiberiu Tesileanu
 */
#ifndef MEXOBJECT_H_
#define MEXOBJECT_H_

#include "mex.h"

#include <string>

/// Interface for wrapping a MEX job inside a class.
class MexObject {
 public:
  /// Initialize with MEX command line and program name
  MexObject(
      const std::string& name,
      int nlhs,
      mxArray* plhs[],
      int nrhs,
      const mxArray* prhs[])
    : nlhs_(nlhs), plhs_(plhs), nrhs_(nrhs), prhs_(prhs), name_(name) {}

  /// Check the inputs; exit with error if one of the inputs is invalid.
  virtual void checkInputs() {}

  /// Initialize input wrappers.
  virtual void initInputs() {}

  /// Initialize helper structures.
  virtual void initHelpers() {}

  /// Initialize everything.
  void init() { checkInputs(); initInputs(); initHelpers(); }

  /// Run the function.
  virtual void run() {}

  /// Issue an error and exit.
  void error(const std::string& id, const std::string& text)
    { mexErrMsgIdAndTxt((name_ + ":" + id).c_str(), text.c_str()); }

 protected:
  int             nlhs_;
  int             nrhs_;
  mxArray**       plhs_;
  const mxArray** prhs_;
  std::string     name_;
};

#endif
