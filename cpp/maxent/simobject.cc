#include "simobject.h"

#include <list>
#include <stdexcept>
#include <string>

#include <stdint.h>

#include "mcmc.h"
#include "utils.h"
#include "vector.h"

// XXX this is not a public API function from Matlab, but it should work to
// detect CTRL+C
extern "C" bool utIsInterruptPending();
extern "C" bool utSetInterruptPending(bool);

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

// one element in the simulation's history
struct HistoryElement {
  Sequence      sequence;
  double        energy;
  double        acceptance;
};

// the history
typedef std::list<HistoryElement> History;

void SimulationObject::checkInputs()
{
  // cppsimmaxent(iniseq, alphaletters, alphawidths, couplings, nsteps, params)
  //   - alphabetletters is a cell array of the results of alphagetletters
  //   - they should include the gap, or be for 'gapped' alphabets
  //   - couplings should include the gaps
  //   - nsteps == 0 means run forever (or until CTRL+C is pressed)
  //   - params is a structure -- see initInputs()

  // check the number of input/output arguments
  if (nrhs_ != 6)
    error("nargin", "Need six arguments.");
  if (nlhs_ > 1)
    error("nargout", "Returns at most one argument.");

  // check the input arguments
  // initial sequence
  if (!mxIsChar(prhs_[0]) || (mxGetM(prhs_[0]) != 1 && mxGetN(prhs_[0]) != 1))
    error("invseq0", "The initial sequence should be a character vector.");
  const size_t seqLen = mxGetM(prhs_[0])*mxGetN(prhs_[0]);

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

  // number of steps
  if (!mxIsDouble(prhs_[4]) || mxGetM(prhs_[4]) != 1 || mxGetN(prhs_[4]) != 1)
    error("invnsteps", "Number of steps should be a scalar.");

  // additional options
  if (!mxIsEmpty(prhs_[5]) && !mxIsStruct(prhs_[5]))
    error("invopts", "Additional options should be provided in a structure.");
}

void SimulationObject::initInputs()
{
  // default options
  burnin_ = 0;            // params.burnin
  rec_step_ = 1;          // params.recstep
  verbose_ = true;        // params.verbose
  disp_interval_ = 10;    // params.dispinterval

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
  const size_t iniseq_len = mxGetM(prhs_[0])*mxGetN(prhs_[0]);
  if (len != iniseq_len)
    error("badwidth", "Initial sequence width does not match alphabet "
      "information.");

  // store the couplings
  couplings_.set(prhs_[3], alphabets_, alphawidths_);

  // store the intial sequence
  seq0_ = to_seq(mxGetChars(prhs_[0]), alphabets_, alphawidths_);

  // store the number of steps
  // we're automatically rounding down
  nsteps_ = mxGetScalar(prhs_[4]);
  // store additional options
  if (!mxIsEmpty(prhs_[5])) {
    const mxArray* tmp = 0;

    // XXX not much error checking here...
    tmp = mxGetField(prhs_[5], 0, "burnin");
    if (tmp) burnin_ = mxGetScalar(tmp);

    tmp = mxGetField(prhs_[5], 0, "recstep");
    if (tmp) rec_step_ = mxGetScalar(tmp);

    tmp = mxGetField(prhs_[5], 0, "verbose");
    if (tmp) verbose_ = (mxGetScalar(tmp) > 1e-6);

    tmp = mxGetField(prhs_[5], 0, "dispinterval");
    if (tmp) disp_interval_ = mxGetScalar(tmp);
  }
}

void SimulationObject::run()
{
  Mcmc mc(seq0_, couplings_);
  History history;

  size_t check_step = 100000;
  size_t i = 0;
  uint64_t tic_id = tic();
  size_t n_accepted = 0;
  size_t n_total = 0;
  while (nsteps_ == 0 || i < nsteps_) {
    // check if we need to look at the time and/or CTRL+C
    if (i % check_step == 0) {
     if (utIsInterruptPending()) { // CRTL+C -- stop early
        utSetInterruptPending(false);
        break;
      }

      if (verbose_) {
        double dt = toc(tic_id);
        if (dt > disp_interval_) {
          mexPrintf("Step: %llu, energy: %f\n", (unsigned long long)i,
            mc.get_energy());
          mexEvalString("drawnow;");
          tic_id = tic();
        }
      }
    }

    // check if we need to store
    if (i >= burnin_ && (i - burnin_) % rec_step_ == 0) {
      HistoryElement element;
      element.sequence = mc.get_sequence();
      element.energy = mc.get_energy();
      element.acceptance = (n_total > 0)?((double)n_accepted / n_total):0;
      history.push_back(element);
      n_accepted = 0;
      n_total = 0;
    }

    // make the step
    if (mc.step()) ++n_accepted;
    ++n_total;

    ++i;
  }

  // send the output to Matlab
  if (nlhs_ > 0) {
    // make a structure for the output
    const char* ret_names[] = {"nsteps", "sequences", "energies",
                               "acceptances"};
    const unsigned n_retfields = sizeof(ret_names)/sizeof(*ret_names);
    plhs_[0] = mxCreateStructMatrix(1, 1, n_retfields, ret_names);

    // fill the output structure
    mxArray* m_nsteps = mxCreateDoubleScalar(i);
    mxSetField(plhs_[0], 0, "nsteps", m_nsteps);

    mwSize s_sizes[2] = {history.size(), seq0_.size()};
    mxArray* m_sequences = mxCreateCharArray(2, s_sizes);
    Matrix<mxChar> sequences(mxGetChars(m_sequences), s_sizes[0], s_sizes[1]);
    mxSetField(plhs_[0], 0, "sequences", m_sequences);

    mxArray* m_energies = mxCreateDoubleMatrix(s_sizes[0], 1, mxREAL);
    Vector<double> energies(mxGetPr(m_energies), s_sizes[0]);
    mxSetField(plhs_[0], 0, "energies", m_energies);

    mxArray* m_acceptances = mxCreateDoubleMatrix(s_sizes[0], 1, mxREAL);
    Vector<double> acceptances(mxGetPr(m_acceptances), s_sizes[0]);
    mxSetField(plhs_[0], 0, "acceptances", m_acceptances);

    size_t j = 0;
    const size_t nalphas = alphabets_.size();
    for (History::const_iterator i = history.begin(); i != history.end();
         ++i, ++j)
    {
      size_t crtidx = 0;
      for (size_t a = 0; a < nalphas; ++a) {
        const std::string& letters = alphabets_[a];
        const size_t w = alphawidths_[a];
        for (size_t p = 0; p < w; ++p, ++crtidx) {
          sequences(j, crtidx) = letters[i -> sequence[crtidx]];
        }
      }

      energies(j) = i -> energy;
      acceptances(j) = i -> acceptance;
    }
  }
}
