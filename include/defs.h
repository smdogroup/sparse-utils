#ifndef SPARSE_UTILS_DEFS_H
#define SPARSE_UTILS_DEFS_H

#include <complex>
#include <cstddef>
#include <limits>

namespace SparseUtils {

/**
 * @brief The type for array indexing.
 *
 * Note: ptrdiff_t is preferable over int because int is often 32 bits and might
 * not be large enough for 64-bit environments. ptrdiff_t is also preferable
 * over size_t because it is signed.
 *
 * Caveat: If the array size is so large (greater than PTRDIFF_MAX but less than
 * SIZE_MAX), element indices may not be representable by ptrdiff_t and the
 * behavior is undefined. But usually ptrdiff_t is 64-bit which is sufficient.
 *
 * Refs:
 * - https://github.com/isocpp/CppCoreGuidelines/pull/1115
 * - https://stackoverflow.com/a/48730597
 * - https://en.cppreference.com/w/cpp/types/ptrdiff_t
 *
 */
using index_t = std::ptrdiff_t;

template <typename T>
using complex_t = std::complex<T>;

// The value indicating that it is not a valid index.
static constexpr index_t NO_INDEX = std::numeric_limits<index_t>::max();

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

}  // namespace SparseUtils

#endif  // SPARSE_UTILS_DEFS_H
