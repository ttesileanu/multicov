#ifndef POTTS_COUPLINGS_H_
#define POTTS_COUPLINGS_H_

#include "mex.h"

#include <string>
#include <vector>

#include "matrix.h"

/** @brief Storage for couplings and fields (fields on diagonal).
 *
 *  Apart from the numeric values of the couplings and fields, this also
 *  stores the structure of the alignment---the alphabets that are used
 *  and where they are used.
 */
struct Couplings {
  /// Datatype used to store the actual couplings and fields.
  typedef Matrix<double> Data;

  /// The actual coupling data.
  Data                  data;
  /// The number of letters in the alphabet at each alignment position.
  std::vector<size_t>   nletters;
  /// The row indices in the coupling matrix where alignment positions start.
  std::vector<size_t>   startidx;

  /// Set the coupling data and the nletters and startidx fields.
  void set(const mxArray* dataorigin,
           const std::vector<std::string>& alphabets_,
           const std::vector<size_t>& alphawidths_);
};

#endif
