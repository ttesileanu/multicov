/** @file mcmc.h
 *  @brief The classes handling the actual Monte Carlo code.
 */
#ifndef POTTS_MCMC_H_
#define POTTS_MCMC_H_

#include "evaluator.h"
#include "matlab_random.h"

/// Markov-chain Monte Carlo class
class Mcmc {
 public:
  /// Construct with initial sequence and couplings.
  explicit Mcmc(const Sequence& seq, const Couplings& couplings);

  /// Make a step in the simulation. Return @a true if accepted.
  bool step();

  /// Get the current energy.
  double get_energy() const { return energy_; }

  /// Get the current sequence.
  const Sequence& get_sequence() const { return seq_; }

  /// Get inverse temperature.
  double get_beta() const { return beta_; }

  /// Set inverse temperature.
  void set_beta(double b) { beta_ = b; }

 private:
  /// Current sequence.
  Sequence      seq_;
  /// Evaluator.
  Evaluator     evaluator_;
  /// Current energy.
  double        energy_;
  /// Random number generator.
  MatlabRandom  random_;
  /// Inverse temperature.
  double        beta_;
};

#endif
