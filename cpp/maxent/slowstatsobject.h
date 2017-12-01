/** @file slowstatsobject.h
 *  @brief Defines the class that calculates the statistics given the maximum
 *         entropy model.
 *
 *  @author Tiberiu Tesileanu
 */
#ifndef POTTS_SLOWSTATSOBJECT_H_
#define POTTS_SLOWSTATSOBJECT_H_

#include "mexobject.h"

#include "evaluator.h"

class SlowStatsObject : public MexObject {
 public:
  /// Constructor from MEX command line.
  SlowStatsObject(int nlhs, mxArray* plhs[], int nrhs, const mxArray* prhs[])
    : MexObject("cppslowstats", nlhs, plhs, nrhs, prhs) {}

  /// Override input checking.
  virtual void checkInputs();

  /// Override input wrappers initialization.
  virtual void initInputs();

  /// Override function run.
  virtual void run();

 private:
  /// The alphabets.
  std::vector<std::string>
                      alphabets_;
  /// The alphabet widths.
  std::vector<size_t> alphawidths_;
  /// Couplings.
  Couplings           couplings_;
  /// Should we output progress information?
  bool                verbose_;
  /// How often to output progress information (in seconds).
  double              disp_interval_;
};

#endif
