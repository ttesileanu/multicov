/** @file simobject.h
 *  @brief Defines the class that performs the MCMC simulation.
 *
 *  @author Tiberiu Tesileanu
 */
#ifndef POTTS_SIMOBJECT_H_
#define POTTS_SIMOBJECT_H_

#include "mexobject.h"

#include "evaluator.h"

class SimulationObject : public MexObject {
 public:
  /// Constructor from MEX command line.
  SimulationObject(int nlhs, mxArray* plhs[], int nrhs, const mxArray* prhs[])
    : MexObject("cppsimmaxent", nlhs, plhs, nrhs, prhs) {}

  /// Override input checking.
  virtual void checkInputs();

  /// Override input wrappers initialization.
  virtual void initInputs();

  /// Override function run.
  virtual void run();

 private:
  /// Initial sequence.
  Sequence            seq0_;
  /// The alphabets.
  std::vector<std::string>
                      alphabets_;
  /// The alphabet widths.
  std::vector<size_t> alphawidths_;
  /// Couplings.
  Couplings           couplings_;
  /// Number of steps in the simulation.
  size_t              nsteps_;
  /// Burn-in period, during which sequences are not recorded.
  size_t              burnin_;
  /// How often to record sequences (number of steps).
  size_t              rec_step_;
  /// Should we output progress information?
  bool                verbose_;
  /// How often to output progress information (in seconds).
  double              disp_interval_;
};

#endif
