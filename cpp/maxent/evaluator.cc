#include "evaluator.h"

double Evaluator::get(const Sequence& s) const
{
  size_t n = s.size();

  double energy = 0;

  size_t jidx = 0;
  for (size_t j = 0; j < n; ++j) {
    const int aj = s[j];
    const int Lj = couplings_.nletters[j];

    // get pointer to beginning of column
    const Couplings::Data::elem_type* p = &couplings_.data(0, jidx + aj);
    for (size_t i = 0; i < n; ++i) {
      const int ai = s[i];
      const int Li = couplings_.nletters[i];

      energy += p[ai];

      // jump to the next position
      p += Li;
    }

    jidx += Lj;
  }

  // take care of the prefactor
  energy *= -0.5;

  return energy;
}

double Evaluator::get_change(const Sequence& s, const Mutation& pert) const
{
  // handle the trivial case
  if (pert.initial == pert.final)
    return 0;

  size_t n = s.size();

  // get index for the initial and final character in modified position
  const size_t idx0 = couplings_.startidx[pert.pos];
  const size_t idx_ak1 = idx0 + pert.initial;
  const size_t idx_ak2 = idx0 + pert.final;

  const Couplings::Data& data = couplings_.data;
  const std::vector<size_t>& nletters = couplings_.nletters;

  // get pointers to the initial and final column
  const Couplings::Data::elem_type* p1 = &data(0, idx_ak1);
  const Couplings::Data::elem_type* p2 = &data(0, idx_ak2);

  double diff = 0.5*(data(idx_ak1, idx_ak1) + data(idx_ak2, idx_ak2));

  // XXX the code for single alphabet is about 2.5x faster...
  // XXX perhaps this could be accelerated by treating the alphabets one
  // XXX by one (that way Li would be constant in the inner loop)

  for (size_t i = 0; i < n; ++i) {
      const int ai = s[i];
      const int Li = nletters[i];

      diff += p1[ai] - p2[ai];

      // jump to the next position
      p1 += Li;
      p2 += Li;
  }

  return diff;
}
