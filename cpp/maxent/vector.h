/** @file vector.h
 *  @brief Defines a vector class that can use either its own, managed storage
 *         space, or storage space provided and managed by the user.
 *
 *  @author Tiberiu Tesileanu
 */
#ifndef VECTOR_H_
#define VECTOR_H_

/** @brief This is a vector-like container that is capable of either managing
 *         its own space or using storage space provided by the user.
 */
template <class T>
class Vector {
 public:
  /// Type of values contained by the vector.
  typedef T         elem_type;
  /// Type used for parameters that look like values of the vector.
  // XXX this is essentially tweaked to PODs, to avoid depending on something
  // that provides call_traits (like BOOST)...
  typedef T         param_type;

  /// Empty constructor.
  Vector() : data_(0), data_owner_(false), n_(0) {}

  /** @brief Constructor with initialization of vector.
   *
   *  Note that this turns on data owning by the vector, and allocates the
   *  necessary space (unless the size is zero).
   */
  Vector(size_t n) : data_owner_(true), n_(n) { allocate_(); }

  /** @brief Constructor with initialization of vector. This also fills the
   *         vector with the given value.
   *
   *  Note that this turns on data owning by the vector, and allocates the
   *  necessary space (unless the size is zero).
   */
  Vector(size_t n, param_type elem) : data_owner_(true), n_(n)
    { allocate_(); fill(elem); }

  /** @brief Constructor using external storage.
   *
   *  This means that the vector does not own the data, and will not deallocate
   *  it upon destruction.
   */
  Vector(elem_type* space, size_t n = 0) : data_(space), data_owner_(false),
      n_(n) {}
  
  /// Copy constructor.
  Vector(const Vector<T>& old) {
    n_ = old.n_;
    data_owner_ = old.data_owner_;
    if (data_owner_) {
      allocate_();
      std::copy(old.data_, old.data_ + numElements(), data_);
    } else {
      data_ = old.data_;
    }
  }

  /// Destructor.
  ~Vector() { free_(); }

  /// Get number of elements.
  size_t numElements() const { return n_; }

  /// Fill the vector with the given element.
  void fill(param_type elem) { std::fill_n(data_, numElements(), elem); }

  /** @brief Resize the vector.
   *
   *  If the vector owns its data, this destroys its contents. Otherwise,
   *  the contents are unaffected.
   */
  void resize(size_t n) {
    // if the vector was already unallocated, allocate space
    n_ = n;
    if (data_) {
      free_();
    } else {
      data_owner_ = true;
    }
    if (data_owner_)
      allocate_();
  }
  
  /** @brief Tell the vector to use the given memory area for storage. The
   *         initial size is 0.
   *
   *  If the vector is in owner mode, it will take ownership of the memory,
   *  and free it with @a delete[] upon destruction.
   */
  void setStorage(elem_type* p) { data_ = p; n_ = 0; }

  /** @brief Tell the vector to use the given memory area for storage.
   *
   *  This sets the vector to non-owner mode.
   */
  void setStorage(elem_type* p, size_t n)
    { data_owner_ = false; data_ = p; n_ = n; }

  /** @brief Set the vector's dimensions to zero.
   *
   *  If the vector owns its data, this deletes all data. Otherwise, the data
   *  is unaffected.
   */
  void clear() { resize(0); }

  /** @brief Turn the vector into the owner of its data. This means that it will
   *         deallocate this data (with @a delete[]) upon destruction!
   *
   *  This does nothing if the vector is already a data owner.
   */
  void ownData() { data_owner_ = true; }

  /// Get access to last element.
  elem_type& back() { return data_[n_ - 1]; }

  /// Get last element (read-only).
  param_type back() const { return data_[n_ - 1]; }

  /// Make the vector not own its data anymore.
  void disownData() { data_owner_ = false; }

  /// Get direct access to the vector's data.
  elem_type* getData() { return data_; }

  /// Get direct read-only access to the vector's data.
  const elem_type* getData() const { return data_; }

  /** @brief This creates a new data storage, that the vector manages.
   *
   *  This function does not carry over any previous contents of the vector.
   */
  void createData() { free_(); ownData(); allocate_(); }

  /// Get number of elements in the vector.
  size_t getN() const { return n_; }

  /// Is the vector empty?
  bool isEmpty() const { return (n_ == 0); }

  /// Access the vector.
  elem_type& operator()(size_t i) {
    return data_[i];
  }

  /// Access the vector read-only-ly.
  const elem_type& operator()(size_t i) const {
    return data_[i];
  }

  /// Assignment operator.
  Vector<T>& operator=(Vector<T> old) { swap(*this, old); return *this; }

  /// Swap.
  friend void swap(Vector<T>& first, Vector<T>& second) {
    using std::swap;

    swap(first.data_, second.data_);
    swap(first.data_owner_, second.data_owner_);
    swap(first.n_, second.n_);
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
  /// Is data owned by the vector?
  bool        data_owner_;
  /// Number of elements in the vector.
  size_t      n_;
};

#endif
