/** @file calcseqw_cpp.cc
 *  @brief Calculate adjacency matrix for an alignment, in the style of DCA --
 *         putting an edge between two sequences if their sequence identity is
 *         above a given threshold.
 *
 *  @author Tiberiu Tesileanu
 */
#include "calcobject.h"

void mexFunction(int nlhs, mxArray* plhs[], int nrhs, const mxArray* prhs[])
{
  CalculateObject worker(nlhs, plhs, nrhs, prhs);
  worker.init();
  worker.run();
}
