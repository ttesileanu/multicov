#include "mcmc.h"

#include <cmath>

#include "matlab_random.h"

Mcmc::Mcmc(const Sequence& seq, const Couplings& couplings)
    : seq_(seq), evaluator_(couplings), energy_(evaluator_.get(seq_)), beta_(1)
{
}

bool Mcmc::step()
{
  // decide on a perturbation
  Mutation mut;
  // need the couplings to get sequence and alphabet size information
  const Couplings& couplings = evaluator_.getCouplings();
  const size_t n = couplings.nletters.size();
  mut.pos = random_.get_int(n);
  mut.initial = seq_[mut.pos];

  const size_t L = couplings.nletters[mut.pos];
  // make sure we don't get the same number we started with
  mut.final = random_.get_int(L - 1);
  if (mut.final >= mut.initial)
    ++mut.final;
  // perturb
  seq_[mut.pos] = mut.final;

  // estimate change in energy
  double dE = evaluator_.get_change(seq_, mut);

  // accept/reject
  if (dE <= 0 || random_.get_double() < exp(-beta_*dE)) {
    // accept step
    energy_ += dE;
    return true;
  } else {
    // reject step
    seq_[mut.pos] = mut.initial;
    return false;
  }
}
