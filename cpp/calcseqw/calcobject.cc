#include "calcobject.h"

#include <cmath>

// implementation of CalculateObject methods
void CalculateObject::checkInputs()
{
  // check the number of inputs
  if (nrhs_ != 2)
    error("nargin", "Need 2 input arguments.");
  // check the number of outputs
  if (nlhs_ > 1)
    error("nargout", "Can have at most 1 output arguments.");

  // calcseqw_cpp(alignmatrix, threshold)

  // check the input arguments
  if (!mxIsChar(prhs_[0]))
    error("invalign", "The alignment data should be a character matrix.");
  if (!mxIsDouble(prhs_[1]) || mxIsSparse(prhs_[1]) || mxIsComplex(prhs_[1]) ||
      mxGetM(prhs_[1]) != 1 || mxGetN(prhs_[1]) != 1)
    error("invmaxseqid", "Max_seqid should be a real scalar.");
}

void CalculateObject::initInputs()
{
  // get convenient access to the alignment data
  alignment_.setStorage(mxGetChars(prhs_[0]),
    mxGetM(prhs_[0]), mxGetN(prhs_[0]));

  maxSeqid_ = *mxGetPr(prhs_[1]);
}

void CalculateObject::run()
{
  const size_t nSeqs = alignment_.nRows();
  const size_t nPos = alignment_.nCols();

  // create the sparse output matrix
  mxArray* mxAdjacency = mxCreateSparseLogicalMatrix(nSeqs, nSeqs, nSeqs);

//  const size_t minDiffs = round(nPos*(1 - maxSeqid_));
  const size_t minDiffs = floor(nPos*(1 - maxSeqid_));
  const double growthFactor = 2;
  size_t nElems = 0;
  mxGetJc(mxAdjacency)[0] = 0; // first column starts at 0
  for (size_t j = 0; j < nSeqs; ++j) {
    for (size_t i = j + 1; i < nSeqs; ++i) {
      size_t nDiffs = 0;
      for (size_t k = 0; k < nPos; ++k)
        if (alignment_(i, k) != alignment_(j, k))
          ++nDiffs;

       if (nDiffs <= minDiffs) {
        // put a 1 in the matrix
        // but first decide whether we need to grow the storage
        const size_t maxStorage = mxGetNzmax(mxAdjacency);
        if (nElems == maxStorage) {
          // need reallocation
          const size_t newStorage = growthFactor*maxStorage;
          mxSetNzmax(mxAdjacency, newStorage);
          mxSetIr(mxAdjacency, (mwIndex*)mxRealloc(mxGetIr(mxAdjacency),
            sizeof(mwIndex)*newStorage));
          // XXX this is very sketchy, but Matlab doesn't have mxSetLogicals...
          mxSetPr(mxAdjacency, (double*)mxRealloc(mxGetLogicals(mxAdjacency),
            sizeof(mxLogical)*newStorage));
        }

        mxGetLogicals(mxAdjacency)[nElems] = 1;
        mxGetIr(mxAdjacency)[nElems++] = i;
      }
    }
    // mark the end of the column
    mxGetJc(mxAdjacency)[j + 1] = nElems;
  }

  if (nlhs_ > 0)
    plhs_[0] = mxAdjacency;
  else
    mxDestroyArray(mxAdjacency);
}
