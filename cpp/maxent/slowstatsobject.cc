#include "slowstatsobject.h"

#include <list>
#include <stdexcept>
#include <string>

#include <cmath>

#include <stdint.h>

#include "utils.h"
#include "vector.h"

// some useful internal functions
namespace {

// call Matlab's tic command
uint64_t tic() {
  mxArray* mx_result;

  mexCallMATLAB(1, &mx_result, 0, 0, "tic");

  // read result
  if (!mxIsUint64(mx_result))
    throw std::runtime_error("Internal: Matlab's tic did not return uint64!");
  uint64_t res = *static_cast<const uint64_t*>(mxGetData(mx_result));

  // free result memory
  mxDestroyArray(mx_result);

  return res;
}

// call Matlab's toc command, using an identifier returned from tic
double toc(uint64_t id) {
  mxArray* mx_result;

  // set up the arguments
  mxArray* mx_id = mxCreateNumericMatrix(1, 1, mxUINT64_CLASS, mxREAL);
  *static_cast<uint64_t*>(mxGetData(mx_id)) = id;

  mexCallMATLAB(1, &mx_result, 1, &mx_id, "toc");

  // read result
  if (!mxIsDouble(mx_result))
    throw std::runtime_error("Internal: Matlab's toc did not return double!");
  double res = mxGetScalar(mx_result);

  // free result memory
  mxDestroyArray(mx_result);

  return res;
}

} // anonymous namespace

void SlowStatsObject::checkInputs()
{
  // cppslowstats(alphaletters, alphawidths, couplings, params)
  //   - alphabetletters is a cell array of the results of alphagetletters
  //   - they should include the gap, or be for 'gapped' alphabets
  //   - couplings should include the gaps
  //   - params is a structure -- see initInputs()

  // check the number of input/output arguments
  if (nrhs_ != 4)
    error("nargin", "Need four arguments.");
  if (nlhs_ > 3)
    error("nargout", "Returns at most three arguments.");

  // check the input arguments
  // initial sequence
  // alphabet letters
  if (!mxIsCell(prhs_[0]) || (mxGetM(prhs_[0]) != 1 && mxGetN(prhs_[0]) != 1))
    error("invalpha","The alphabetletters should be a cell array of strings.");
  const size_t nalphas = std::max(mxGetM(prhs_[0]), mxGetN(prhs_[0]));

  // alphawidths
  if (!mxIsDouble(prhs_[1]) || (mxGetM(prhs_[1])!=1 && mxGetN(prhs_[1])!=1))
    error("invwidth", "The alphawidths should be a numeric vector.");
  if (nalphas != std::max(mxGetM(prhs_[1]), mxGetN(prhs_[1])))
    error("missalpha", "Alphabetletters and alphawidths should have matching "
      "lengths.");

  // coupling matrix
  if (!mxIsDouble(prhs_[2]) || (mxGetM(prhs_[2]) != mxGetN(prhs_[2])) ||
      mxIsSparse(prhs_[2]) || mxIsComplex(prhs_[2]))
    error("invJ", "The coupling matrix should be non-sparse, real, square, and"
                  " contain doubles.");

  // additional options
  if (!mxIsEmpty(prhs_[3]) && !mxIsStruct(prhs_[3]))
    error("invopts", "Additional options should be provided in a structure.");
}

void SlowStatsObject::initInputs()
{
  // default options
  verbose_ = true;        // params.verbose
  disp_interval_ = 10;    // params.dispinterval

  // store the alphabets
  const size_t nalphas = std::max(mxGetM(prhs_[0]), mxGetN(prhs_[0]));
  double* alphaw_p = mxGetPr(prhs_[1]);
  size_t len = 0;
  for (size_t i = 0; i < nalphas; ++i) {
    mxArray* crtelem = mxGetCell(prhs_[0], i);
    if (!mxIsChar(crtelem))
      error("invalpha", "The alphabet letters should be a cell array of "
        "strings.");
    
    alphabets_.push_back(to_string(crtelem));
    alphawidths_.push_back(alphaw_p[i]);
    len += alphawidths_.back();
  }

  // store the couplings
  couplings_.set(prhs_[2], alphabets_, alphawidths_);

  // store additional options
  if (!mxIsEmpty(prhs_[3])) {
    const mxArray* tmp = 0;

    // XXX not much error checking here...
    tmp = mxGetField(prhs_[3], 0, "verbose");
    if (tmp)
      verbose_ = (mxGetScalar(tmp) > 1e-6);

    tmp = mxGetField(prhs_[3], 0, "dispinterval");
    if (tmp)
      disp_interval_ = mxGetScalar(tmp);
  }
}

size_t intpow(size_t x, size_t n)
{
  switch (n) {
    case 0:
      return 1;
    case 1:
      return x;
    case 2:
      return x*x;
    default:;
  };

  const size_t y = intpow(x, n/2);
  if (n % 2 == 0)
    return y*y;
  else
    return x*y*y;
}

void SlowStatsObject::run()
{
  // XXX hopefully size_t is 64-bit!
  size_t nsteps = 1;
  size_t n = 0;
  for (size_t k = 0; k < alphawidths_.size(); ++k) {
    nsteps *= intpow(alphabets_[k].size(), alphawidths_[k]);
    n += alphawidths_[k];
  }

  Sequence state(n, 0);

  // create the binary alignment map, and a similar map for the couplings
  std::vector<size_t> binaryMap;
  std::vector<size_t> couplingMap;
  std::vector<size_t> sizesPerPos;

  binaryMap.resize(n);
  couplingMap.resize(n);
  sizesPerPos.resize(n);

  size_t a = 0; // this tracks the position in the coupling matrix
  size_t b = 0; // this tracks the position in the binary alignment
  size_t alphaidx = 0;
  size_t crt = 0;
  for (size_t k = 0; k < n; ++k) {
    couplingMap[k] = a;
    binaryMap[k] = b;

    const size_t size = alphabets_[alphaidx].size();
    sizesPerPos[k] = size;
    a += size;
    b += size - 1; // gap ignored in binary alignment
    ++crt;

    if (crt >= alphawidths_[alphaidx]) {
      crt = 0;
      ++alphaidx;
    }
  }

  // create the output matrices
  const size_t binSize = b;
  mxArray* mxfi = mxCreateDoubleMatrix(binSize, 1, mxREAL);
  mxArray* mxfij = mxCreateDoubleMatrix(binSize, binSize, mxREAL);

  // get convenient access to outputs
  Vector<double> fi(mxGetPr(mxfi), binSize);
  Matrix<double> fij(mxGetPr(mxfij), binSize, binSize);

  Evaluator evaluator(couplings_);
  double energy = evaluator.get(state);
  
  // partition function
  double Z = 0;

  const size_t check_step = 100000;
  size_t i = 0;
  uint64_t tic_id = tic();
  while (state[0] < sizesPerPos[0]) {
    // test!
/*    if (std::abs(evaluator.get(state) - energy) > 1e-6)
      throw std::runtime_error("bad energy!");*/

    // check if we need to look at the time and/or CTRL+C
    if (i % check_step == 0) {
      // make sure we don't fall out of sync
      energy = evaluator.get(state);
      if (verbose_) {
        double dt = toc(tic_id);
        if (dt > disp_interval_) {
          mexPrintf("%f\% done\n", 100*(double)i / nsteps);
          mexEvalString("drawnow;");
          tic_id = tic();
        }
      }
    }

    const double exp_energy = std::exp(-energy);
    Z += exp_energy;

    // update fi, fij
    for (size_t x = 0; x < n; ++x) {
      int nx = state[x];
      if (nx < 1) // skip gaps
        continue;
      size_t idxx = binaryMap[x] + nx - 1;

      fi(idxx) += exp_energy;
      fij(idxx, idxx) += exp_energy;

      for (size_t y = x + 1; y < n; ++y) {
        int ny = state[y];
        if (ny < 1) // skip gaps
          continue;
        size_t idxy = binaryMap[y] + ny - 1;

        fij(idxx, idxy) += exp_energy;
        fij(idxy, idxx) += exp_energy;
      }
    }

    ++i;

    // iterate
    int oldValue = state.back();
    ++state.back();
    size_t k = n - 1;
    while (state[k] >= sizesPerPos[k] && k > 0) {
      state[k] = 0;
      energy += evaluator.get_change(state, Mutation(k, oldValue, 0));

      --k;
      oldValue = state[k];
      ++state[k];
    }

    if (state[k] < sizesPerPos[k])
      energy += evaluator.get_change(state, Mutation(k, oldValue, state[k]));
  }

  // normalize
  for (size_t x = 0; x < binSize; ++x) {
    fi(x) /= Z;
    for (size_t y = 0; y < binSize; ++y) {
      fij(x, y) /= Z;
    }
  }

  // send the output to Matlab
  if (nlhs_ > 0)
    plhs_[0] = mxfi;
  else
    mxDestroyArray(mxfi);

  if (nlhs_ > 1)
    plhs_[1] = mxfij;
  else
    mxDestroyArray(mxfij);

  if (nlhs_ > 2)
    plhs_[2] = mxCreateDoubleScalar(Z);
}
