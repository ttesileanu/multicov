/** @file matrix.h
 *  @brief Defines a matrix class that can use either its own, managed storage
 *         space, or storage space provided and managed by the user.
 *
 *  Another feature is the ability to use either row-major or column-major
 *  ordering. This is useful to interface with Matlab.
 *
 *  @author Tiberiu Tesileanu
 */
#ifndef MATRIX_H_
#define MATRIX_H_

/** @brief This is a matrix-like container that is capable of either managing
 *         its own space or using storage space provided by the user.
 *
 *  Another feature of this class is that, apart from having the possibility of
 *  containing any value type, it can also choose between using column-major
 *  ordering (the default; Matlab style), or row-major ordering (C style).
 */
template <class T, bool col = true>
class Matrix {
 public:
  /// Type of values contained by the matrix.
  typedef T         elem_type;
  /// Type used for parameters that look like values of the matrix.
  // XXX this is essentially tweaked to PODs, to avoid depending on something
  // that provides call_traits (like BOOST)...
  typedef T         param_type;

  /// Is the matrix column-order?
  static const bool columnOrder = col;

  /// Empty constructor.
  Matrix() : data_(0), data_owner_(false), rows_(0), cols_(0) {}

  /** @brief Constructor with initialization of matrix.
   *
   *  Note that this turns on data owning by the matrix, and allocates the
   *  necessary space (unless the size is zero).
   */
  Matrix(size_t rows, size_t cols)
      : data_owner_(true), rows_(rows), cols_(cols) { allocate_(); }

  /** @brief Constructor with initialization of matrix. This also fills the
   *         matrix with the given value.
   *
   *  Note that this turns on data owning by the matrix, and allocates the
   *  necessary space (unless the size is zero).
   */
  Matrix(size_t rows, size_t cols, param_type elem)
      : data_owner_(true), rows_(rows), cols_(cols)
    { allocate_(); fill(elem); }

  /** @brief Constructor using external storage.
   *
   *  This means that the matrix does not own the data, and will not deallocate
   *  it upon destruction.
   */
  Matrix(elem_type* space, size_t rows = 0, size_t cols = 0)
      : data_(space), data_owner_(false), rows_(rows), cols_(cols) {}
  
  /// Copy constructor.
  Matrix(const Matrix<T>& old) {
    rows_ = old.rows_; cols_ = old.cols_;
    data_owner_ = old.data_owner_;
    if (data_owner_) {
      allocate_();
      std::copy(old.data_, old.data_ + numElements(), data_);
    } else {
      data_ = old.data_;
    }
  }

  /// Destructor.
  ~Matrix() { free_(); }

  /// Get number of elements.
  size_t numElements() const { return rows_*cols_; }

  /// Fill the matrix with the given element.
  void fill(param_type elem) { std::fill_n(data_, numElements(), elem); }

  /** @brief Resize the matrix.
   *
   *  If the matrix owns its data, this destroys its contents. Otherwise,
   *  the contents are unaffected.
   */
  void resize(size_t rows, size_t cols) {
    // if the matrix was already unallocated, allocate space
    rows_ = rows; cols_ = cols;
    if (data_) {
      free_();
    } else {
      data_owner_ = true;
    }
    if (data_owner_)
      allocate_();
  }
  
  /** @brief Tell the matrix to use the given memory area for storage. The
   *         initial size of the matrix is 0x0.
   *
   *  If the matrix is in owner mode, it will take ownership of the memory,
   *  and free it with @a delete[] upon destruction.
   */
  void setStorage(elem_type* p) { data_ = p; rows_ = 0; cols_ = 0; }

  /** @brief Tell the matrix to use the given memory area for storage.
   *
   *  This sets the matrix to a non-owner mode.
   */
  void setStorage(elem_type* p, size_t rows, size_t cols)
    { data_owner_ = false; data_ = p; rows_ = rows; cols_ = cols; }

  /** @brief Set the matrix's dimensions to zero.
   *
   *  If the matrix owns its data, this deletes all data. Otherwise, the data
   *  is unaffected.
   */
  void clear() { resize(0, 0); }

  /** @brief Turn the matrix into the owner of its data. This means that it will
   *         deallocate this data (with @a delete[]) upon destruction!
   *
   *  This does nothing if the matrix is already a data owner.
   */
  void ownData() { data_owner_ = true; }

  /// Make the matrix not own its data anymore.
  void disownData() { data_owner_ = false; }

  /// Get direct access to the matrix's data.
  elem_type* getData() { return data_; }

  /// Get direct read-only access to the matrix's data.
  const elem_type* getData() const { return data_; }

  /** @brief This creates a new data storage, that the matrix manages.
   *
   *  This function does not carry over any previous contents of the matrix.
   */
  void createData() { free_(); ownData(); allocate_(); }

  /// Get number of rows in array.
  size_t nRows() const { return rows_; }

  /// Get number of columns in array.
  size_t nCols() const { return cols_; }

  /// Is the matrix empty? (i.e., has zero columns and zero rows)
  bool isEmpty() const { return (rows_ == 0 && cols_ == 0); }

  /// Access the matrix.
  elem_type& operator()(size_t i, size_t j) {
    return (col?data_[i + j*rows_]:data_[i*cols_ + j]);
  }

  /// Access the matrix read-only-ly.
  param_type operator()(size_t i, size_t j) const {
    return (col?data_[i + j*rows_]:data_[i*cols_ + j]);
  }

  /// Assignment operator.
  Matrix<T>& operator=(Matrix<T> old) { swap(*this, old); return *this; }

  /// Swap.
  friend void swap(Matrix<T>& first, Matrix<T>& second) {
    using std::swap;

    swap(first.data_, second.data_);
    swap(first.data_owner_, second.data_owner_);
    swap(first.rows_, second.rows_);
    swap(first.cols_, second.cols_);
  }

 private:
  /// Allocate own storage.
  void allocate_() {
    size_t n = numElements();
    if (n) data_ = new elem_type[n];
      else data_ = 0;
  }
  /// Free any allocated storage, if we own it.
  void free_() {
    if (data_owner_) {
      delete[] data_;
      data_ = 0;
    }
  }

  /// Pointer to data storage.
  elem_type*  data_;
  /// Is data owned by the matrix?
  bool        data_owner_;
  /// Number of rows in the matrix.
  size_t      rows_;
  /// Number of columns in the matrix.
  size_t      cols_;
};

#endif
