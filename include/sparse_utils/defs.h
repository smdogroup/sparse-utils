#ifndef SPARSE_UTILS_DEFS_H
#define SPARSE_UTILS_DEFS_H

#include <complex>
#include <cstddef>
#include <limits>

namespace SparseUtils {

template <typename T>
using complex_t = std::complex<T>;

// The value indicating that it is not a valid index.
static constexpr int NO_INDEX = std::numeric_limits<int>::max();

// Free heap memory and set pointer to nullptr
#define FREE_ARRAY(array_ptr) \
  if (array_ptr) {            \
    delete[] array_ptr;       \
    array_ptr = nullptr;      \
  }

/**
 * @brief Check if a type is complex.
 *
 * Usage:
 *   given an arbitrary type T, if T is a complex type, then
 *     is_complex<T>::value == true, otherwise
 *     is_complex<T>::value == false
 */
template <typename T>
struct is_complex : public std::false_type {};

template <typename T>
struct is_complex<complex_t<T>> : public std::true_type {};

// Get absolute value of the real part
inline double absreal(complex_t<double> a) {
  if (a.real() >= 0.0) {
    return a.real();
  } else {
    return -a.real();
  }
}

inline double absreal(double a) {
  if (a >= 0.0) {
    return a;
  } else {
    return -a;
  }
}

}  // namespace SparseUtils

#endif  // SPARSE_UTILS_DEFS_H
