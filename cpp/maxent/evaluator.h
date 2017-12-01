/** @file evaluator.h
 *  @brief Define classes used to evaluate Potts model energies.
 */
#ifndef POTTS_EVALUATOR_H_
#define POTTS_EVALUATOR_H_

#include <vector>

#include "couplings.h"

/// A sequence is a vector of @a int.
typedef std::vector<int> Sequence;

/// Keeps track of a single-point mutation.
class Mutation {
 public:
  /// Position where mutation happened.
  size_t    pos;
  /// Initial character at that position.
  int       initial;
  /// Final character at that position.
  int       final;

  Mutation() : pos(0), initial(0), final(0) {}
  Mutation(size_t k, int i, int f) : pos(k), initial(i), final(f) {}
};

/** @brief Class that calculates or updates the Potts energy for a given
           sequence.
 *
 *  The Potts energy of a sequence @f$[a_1, \cdots, a_n]@f$ is given by
 *    @f[E(a_1, \cdots, a_n) = -\sum_{i<j} J_{ij}(a_i, a_j) - \sum_i h_i(a_i)@f]
 *  where @f$J_{ij}(a, b)@f$ is the coupling between character @a a at site
 *  @a i and character @a b at site @a j, and @f$h_i(a)@f$ is the value of the
 *  field at character @a a at site @a i.
 *
 *  Here we assume that the fields are included on the diagonal of the coupling
 *  matrix,
 *    @f[J_{ii}(a, a) = 2 h_i(a)@f]
 *  so that we can write
 *    @f[E(a_1, \cdots, a_n) = -\frac 12 \sum_{i,j} J_{ij}(a_i, a_j)@f]
 */
class Evaluator {
 public:
  /** @brief Constructor with initialization.
   *
   *  Note that this class assumes that the coupling matrix is symmetric.
   *  Results are undefined (and will likely be wrong) if it is not.
   */
  explicit Evaluator(const Couplings& couplings) : couplings_(couplings) {}

  /** @brief Fresh evaluation.
   *
   *  This is O(n^2), where n is the length of the sequence.
   */
  double get(const Sequence& s) const;
  /** @brief Update evaluation.
   *
   *  The sequence parameter @a s should be the sequence *after* the
   *  mutation.
   *
   *  This runs in O(n) normally, O(1) if the perturbation is trivial (target
   *  character is equal to source). When using this repeatedly to handle
   *  multi-character mutations, it becomes more efficient to use @a get()
   *  when more than about a quarter of the characters are mutated.
   */
  double get_change(const Sequence& s, const Mutation& pert) const;

  /// Get access to the couplings.
  const Couplings getCouplings() const { return couplings_; }

 private:
  Couplings     couplings_;
};

#endif
