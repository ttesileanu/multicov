#include "evalobject.h"

#include <algorithm>
#include <string>

#include "utils.h"
#include "vector.h"

void EvaluationObject::checkInputs()
{
  // cppevalmaxent(sequences, alphabetletters, alphawidths, couplings)
  //   - sequences should be an array in which each sequences is a *column*
  //   - this is opposite to the usual sequences in MultiCOV, but makes it
  //     easier to access each sequence
  //   - alphabetletters is a cell array of the results of alphagetletters
  //   - they should include the gap, or be for 'gapped' alphabets
  //   - couplings should include the gaps

  // check the number of input/output arguments
  if (nrhs_ != 4)
    error("nargin", "Need four arguments.");
  if (nlhs_ > 1)
    error("nargout", "Returns at most one argument.");

  // check the input arguments
  // sequences
  if (!mxIsChar(prhs_[0]))
    error("invseqs", "The sequences should be a character matrix.");

  // alphabet letters
  if (!mxIsCell(prhs_[1]) || (mxGetM(prhs_[1]) != 1 && mxGetN(prhs_[1]) != 1))
    error("invalpha","The alphabetletters should be a cell array of strings.");
  const size_t nalphas = std::max(mxGetM(prhs_[1]), mxGetN(prhs_[1]));

  // alphawidths
  if (!mxIsDouble(prhs_[2]) || (mxGetM(prhs_[2])!=1 && mxGetN(prhs_[2])!=1))
    error("invwidth", "The alphawidths should be a numeric vector.");
  if (nalphas != std::max(mxGetM(prhs_[2]), mxGetN(prhs_[2])))
    error("missalpha", "Alphabetletters and alphawidths should have matching "
      "lengths.");

  // coupling matrix
  if (!mxIsDouble(prhs_[3]) || (mxGetM(prhs_[3]) != mxGetN(prhs_[3])) ||
      mxIsSparse(prhs_[3]) || mxIsComplex(prhs_[3]))
    error("invJ", "The coupling matrix should be non-sparse, real, square, and"
                  " contain doubles.");
}

void EvaluationObject::initInputs()
{
  // get convenient access to the alignment data
  alignment_.setStorage(mxGetChars(prhs_[0]),
    mxGetM(prhs_[0]), mxGetN(prhs_[0]));

  // store the alphabets
  const size_t nalphas = std::max(mxGetM(prhs_[1]), mxGetN(prhs_[1]));
  double* alphaw_p = mxGetPr(prhs_[2]);
  size_t len = 0;
  for (size_t i = 0; i < nalphas; ++i) {
    mxArray* crtelem = mxGetCell(prhs_[1], i);
    if (!mxIsChar(crtelem))
      error("invalpha", "The alphabet letters should be a cell array of "
        "strings.");
    
    alphabets_.push_back(to_string(crtelem));
    alphawidths_.push_back(alphaw_p[i]);
    len = len + alphawidths_.back();
  }
  if (len != alignment_.nRows())
    error("badwidth", "Alignment width does not match alphabet information.");

  // store the couplings
  couplings_.set(prhs_[3], alphabets_, alphawidths_);
}

void EvaluationObject::run()
{
  // create the output vector
  mxArray* m_energies = mxCreateDoubleMatrix(alignment_.nCols(), 1, mxREAL);
  Vector<double> energies(mxGetPr(m_energies), alignment_.nCols());

  Evaluator evaluator(couplings_);

  for (size_t k = 0; k < alignment_.nCols(); ++k) {
    Sequence seq = to_seq(&alignment_(0, k), alphabets_, alphawidths_);
    energies(k) = evaluator.get(seq);
  }

  plhs_[0] = m_energies;
}
