#ifndef SPARSE_UTILS_UTILS_H
#define SPARSE_UTILS_UTILS_H

#include <algorithm>

#include "defs.h"
#include "sparse_matrix.h"

namespace SparseUtils {

// Compute y = alpha * A * x + beta * y
template <typename T>
void CSRMatVec(double alpha, index_t nrows, const index_t *rowp,
               const index_t *cols, const T *Avals, const T *x, double beta,
               T *y);

// Compute A * x -> y
template <typename T>
void CSCMatVec(double alpha, index_t nrows, index_t ncols, const index_t *colp,
               const index_t *rows, const T *Avals, const T *x, double beta,
               T *y);

// Based on the pattern of A, compute A^{T}. The numerical values are optional
template <typename T>
void SparseTranspose(index_t nrows, index_t ncols, const index_t *rowp,
                     const index_t *cols, const T *Avals, index_t *colp,
                     index_t *rows, T *ATvals);

// Compute the matrix-matrix product A * A^{T}
template <typename T>
void MatMatTransNumeric(index_t nrows, index_t ncols, const index_t *rowp,
                        const index_t *cols, const T *Avals,
                        const index_t *colp, const index_t *rows,
                        const T *ATvals, const index_t *Bcolp, index_t *Brows,
                        T *Bvals, index_t *flag, T *tmp);

// Compute the result C + A * D * A^{T}, where C and D are diagonal
template <typename T>
void MatMatTransNumeric(index_t nrows, index_t ncols, const T *cvals,
                        const index_t *rowp, const index_t *cols,
                        const T *Avals, const T *dvals, const index_t *colp,
                        const index_t *rows, const T *ATvals,
                        const index_t *Bcolp, index_t *Brows, T *Bvals,
                        index_t *flag, T *tmp);

// Compute the number of entries in the matrix product A * A^{T}
inline index_t MatMatTransSymbolic(index_t nrows, index_t ncols,
                                   const index_t *rowp, const index_t *cols,
                                   const index_t *colp, const index_t *rows,
                                   index_t *Bcolp, index_t *flag);

// Sort an array of length len, then remove duplicate entries and
// entries with values -1.
inline index_t RemoveDuplicates(index_t *array, index_t len,
                                index_t exclude = -1);

// Sort and make the data structure unique - remove diagonal
inline void SortAndRemoveDuplicates(index_t nvars, index_t *rowp, index_t *cols,
                                    index_t remove_diagonal = 0);

// Convert BSRMat to an unblocked, CSR format
template <typename T, index_t M, index_t N>
CSRMat<T> *bsr_to_csr(BSRMat<T, M, N> *bsr_mat);

// Convert BSRMat to an unblocked, CSC format
template <typename T, index_t M, index_t N>
CSCMat<T> *bsr_to_csc(BSRMat<T, M, N> *bsr_mat);

}  // namespace SparseUtils

#include "detail/utils_impl.h"

#endif  // SPARSE_UTILS_UTILS_H