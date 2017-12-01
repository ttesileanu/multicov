/** @file utils.h
 *  @brief Some useful functions for the maximum entropy model evaluator and
 *         MCMC simulation.
 */
#ifndef POTTS_UTILS_H_
#define POTTS_UTILS_H_

#include "mex.h"

#include <string>
#include <vector>

#include "evaluator.h"

/** @brief Convert a Matlab string to a C++ string.
 *
 *  Note that this does not check that the input is a string.
 */
std::string to_string(const mxArray* v0);

/** @brief Convert a Matlab character vector to a numeric sequence.
 *
 *  This uses @a alphabets, a vector of strings each of which is a list of the
 *  letters in that alphabet, and @a alphawidths, a vector of number of
 *  columns taken by each alphabet, to perform the conversion.
 */
Sequence to_seq(const mxChar* s, const std::vector<std::string>& alphabets,
    const std::vector<size_t> alphawidths);

#endif
