/** @file matlab_random.h
 *  @brief Generate random numbers by calling Matlab functions.
 */
#ifndef POTTS_MATLAB_RANDOM_H_
#define POTTS_MATLAB_RANDOM_H_

#include "mex.h"

/** @brief Convenience class for calling Matlab's random generator functions.
 *
 *  Random numbers are obtained in chunks, to avoid the overhead of repeated
 *  Matlab function calls. These numbers are then returned to the user when
 *  requested.
 */
class MatlabRandom {
 public:
  MatlabRandom() : buffer_(0), idx_(0), buffer_size_(1000000), rand_arg_(0)
    { set_size_(); replenish_(); }
  ~MatlabRandom() {
    if (buffer_) mxDestroyArray(buffer_);
    if (rand_arg_) mxDestroyArray(rand_arg_);
  }

  /// Return a floating point number in the interval [0, 1).
  double get_double() const {
    // replenish buffer if needed
    if (idx_ >= buffer_size_) replenish_();
    return mxGetPr(buffer_)[idx_++];
  }
  /// Return a floating point number in the interval [0, @a x).
  double get_double(double x) const { return x*get_double(); }
  /// Return a positive integer in the interval [0, n).
  size_t get_int(size_t n) const { return n*get_double(); }

  /// Get the buffer size.
  size_t get_buffer_size() const { return buffer_size_; }
  /// Set the buffer size. This discards current buffer.
  void set_buffer_size(size_t n)
    { buffer_size_ = n; set_size_(); replenish_(); }

 private:
  /// Create the array that is used as argument for rand.
  void set_size_() const;
  /// Replenish the buffer.
  void replenish_() const;

  /// Random number buffer.
  mutable mxArray*  buffer_;
  /// Index in buffer.
  mutable size_t    idx_;
  /// Buffer size.
  size_t            buffer_size_;
  /// Array needed as argument for rand.
  mutable mxArray*  rand_arg_;
};

#endif
