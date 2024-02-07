
#ifndef SPARSE_UTILS_SPARSE_MATRIX_IMPL_H
#define SPARSE_UTILS_SPARSE_MATRIX_IMPL_H

#include "sparse_matrix.h"

namespace SparseUtils {

// Zero the entries of the sparse matrix
template <typename T, index_t M, index_t N>
void BSRMat<T, M, N>::zero() {
  std::fill(vals, vals + M * N * nnz, T(0.0));
}

/**
 * @brief Find the value index of a block given indices (row, col) of a block
 *
 * @param row block row index
 * @param col block column index
 * @return index_t the value index jp such that vals[jp] gives the
 * sought block, if (row, col) isn't in the nonzero pattern, NO_INDEX will be
 * returned
 */
template <typename T, index_t M, index_t N>
index_t BSRMat<T, M, N>::find_value_index(index_t row, index_t col) {
  index_t jp_start = rowp[row];
  index_t jp_end = rowp[row + 1];

  for (index_t jp = jp_start; jp < jp_end; jp++) {
    if (cols[jp] == col) {
      return jp;
    }
  }

  return NO_INDEX;
}

/**
 * @brief add values from an element matrix mat of shape (m, n)
 *
 * @tparam Mat a matrix type whose entry can be ()-indexed
 * @param m number of rows of mat
 * @param i global row indices for each entry of mat
 * @param n number of columns of mat
 * @param j global column indices for each entry of mat
 * @param mat the element matrix
 */
template <typename T, index_t M, index_t N>
template <class Mat>
void BSRMat<T, M, N>::add_values(const index_t m, const index_t i[],
                                 const index_t n, const index_t j[], Mat &mat) {
  for (index_t ii = 0; ii < m; ii++) {
    index_t block_row = i[ii] / M;
    index_t local_row = i[ii] % M;

    for (index_t jj = 0; jj < n; jj++) {
      index_t block_col = j[jj] / N;
      index_t local_col = j[jj] % N;

      index_t jp = find_value_index(block_row, block_col);
      if (jp != NO_INDEX) {
        vals[M * N * jp + N * local_row + local_col] += mat(ii, jj);
      }
    }
  }
}

/**
 * @brief add values from an block element matrix mat of shape (m * M, n * N)
 *
 * @tparam Mat a matrix type whose entry can be ()-indexed
 * @param m number of row blocks of mat
 * @param ib row block indices for each block of mat
 * @param n number of column blocks of mat
 * @param j column block indices for each block of mat
 * @param mat the element matrix
 */
template <typename T, index_t M, index_t N>
template <class Mat>
void BSRMat<T, M, N>::add_block_values(const index_t m, const index_t ib[],
                                       const index_t n, const index_t jb[],
                                       Mat &mat) {
  for (index_t ii = 0; ii < m; ii++) {
    index_t block_row = ib[ii];

    for (index_t jj = 0; jj < n; jj++) {
      index_t block_col = jb[jj];

      index_t jp = find_value_index(block_row, block_col);
      if (jp != NO_INDEX) {
        for (index_t local_row = 0; local_row < M; local_row++) {
          for (index_t local_col = 0; local_col < N; local_col++) {
            vals[M * N * jp + N * local_row + local_col] +=
                mat(M * ii + local_row, N * jj + local_col);
          }
        }
      }
    }
  }
}

/**
 * @brief Zero out rows and set diagonal entry to one for each zeroed row
 *
 * @param nbcs number of global rows to zero-out
 * @param dof global dof indices
 */
template <typename T, index_t M, index_t N>
void BSRMat<T, M, N>::zero_rows(const index_t nbcs, const index_t dof[]) {
  for (index_t ii = 0; ii < nbcs; ii++) {
    index_t block_row = dof[ii] / M;
    index_t local_row = dof[ii] % M;

    for (index_t jp = rowp[block_row]; jp < rowp[block_row + 1]; jp++) {
      for (index_t k = 0; k < N; k++) {
        vals[M * N * jp + N * local_row + k] = 0.0;
      }

      if (cols[jp] == block_row) {
        vals[M * N * jp + N * local_row + local_row] = 1.0;
      }
    }
  }
}

/**
 * @brief Convert to a dense matrix
 *
 * @param m_ output, number of rows of the dense matrix
 * @param n_ output, number of columns of the dense matrix
 * @param A_ output, dense matrix values stored in row-major ordering
 */
template <typename T, index_t M, index_t N>
void BSRMat<T, M, N>::to_dense(index_t *m_, index_t *n_, T **A_) {
  index_t m = M * nbrows;
  index_t n = N * nbcols;
  index_t size = m * n;

  T *A = new T[size];
  std::fill(A, A + size, T(0.0));

  for (index_t i = 0; i < nbrows; i++) {
    for (index_t jp = rowp[i]; jp < rowp[i + 1]; jp++) {
      index_t j = cols[jp];

      for (index_t ii = 0; ii < M; ii++) {
        const index_t irow = M * i + ii;
        for (index_t jj = 0; jj < N; jj++) {
          const index_t jcol = N * j + jj;
          A[n * irow + jcol] = vals[M * N * jp + N * ii + jj];
        }
      }
    }
  }

  *A_ = A;
  *m_ = m;
  *n_ = n;
}

/**
 * @brief Export the matrix as mtx format
 *
 * @param mtx_name the output file
 * @param epsilon only write entries such that abs(e) >= epsilon
 */
template <typename T, index_t M, index_t N>
void BSRMat<T, M, N>::write_mtx(const std::string mtx_name, double epsilon) {
  // Open file and destroy old contents, if any
  std::FILE *fp = std::fopen(mtx_name.c_str(), "w");

  // Write header
  if (is_complex<T>::value) {
    std::fprintf(fp, "%%%%MatrixMarket matrix coordinate complex general\n");
  } else {
    std::fprintf(fp, "%%%%MatrixMarket matrix coordinate real general\n");
  }

  // Write global m, n and nnz
  std::fprintf(fp, "%15d%15d%15d\n", nbrows * M, nbcols * N, nnz * M * N);

  // Write entries
  index_t nnz_mtx = 0;
  for (index_t i = 0; i < nbrows; i++) {
    for (index_t jp = rowp[i]; jp < rowp[i + 1]; jp++) {
      index_t j = cols[jp];  // (i, j) is the block index pair

      for (index_t ii = 0; ii < M; ii++) {
        const index_t irow = M * i + ii + 1;  // convert to 1-based index
        for (index_t jj = 0; jj < N; jj++) {
          // (irow, jcol) is the entry coo
          const index_t jcol = N * j + jj + 1;  // convert to 1-based index
          T val = vals[M * N * jp + N * ii + jj];
          if constexpr (is_complex<T>::value) {
            std::fprintf(fp, "%d %d %30.20e %30.20e\n", irow, jcol, val.real(),
                         val.imag());
            nnz_mtx++;
          } else if (absfunc(val) >= epsilon) {
            std::fprintf(fp, "%d %d %30.20e\n", irow, jcol, val);
            nnz_mtx++;
          }
        }
      }
    }
  }
  std::fclose(fp);

  // Modify nnz
  fp = std::fopen(mtx_name.c_str(), "r+");
  if (is_complex<T>::value) {
    std::fprintf(fp, "%%%%MatrixMarket matrix coordinate complex general\n");
  } else {
    std::fprintf(fp, "%%%%MatrixMarket matrix coordinate real general\n");
  }
  std::fprintf(fp, "%15d%15d%15d", nbrows * M, nbcols * N, nnz_mtx);
  std::fclose(fp);
  return;
}

/**
 * @brief Convert to a dense matrix
 *
 * @param m_ output, number of rows of the dense matrix
 * @param n_ output, number of columns of the dense matrix
 * @param A_ output, dense matrix values stored in row-major ordering
 */
template <typename T>
void CSRMat<T>::to_dense(index_t *m_, index_t *n_, T **A_) {
  index_t m = nrows;
  index_t n = ncols;
  index_t size = m * n;

  T *A = new T[size];
  std::fill(A, A + size, T(0.0));

  for (index_t i = 0; i < nrows; i++) {
    for (index_t jp = rowp[i]; jp < rowp[i + 1]; jp++) {
      index_t j = cols[jp];
      A[n * i + j] = vals[jp];
    }
  }

  *A_ = A;
  *m_ = m;
  *n_ = n;
}

/**
 * @brief Export the matrix as mtx format
 *
 * @param mtx_name the output file
 * @param epsilon only write entries such that abs(e) >= epsilon
 */
template <typename T>
void CSRMat<T>::write_mtx(const std::string mtx_name, double epsilon) {
  // Open file and destroy old contents, if any
  std::FILE *fp = std::fopen(mtx_name.c_str(), "w");

  // Write header
  if (is_complex<T>::value) {
    std::fprintf(fp, "%%%%MatrixMarket matrix coordinate complex general\n");
  } else {
    std::fprintf(fp, "%%%%MatrixMarket matrix coordinate real general\n");
  }

  // Write m, n and nnz
  std::fprintf(fp, "%15d%15d%15d\n", nrows, ncols, nnz);

  // Write entries
  index_t nnz_mtx = 0;
  for (index_t i = 0; i < nrows; i++) {
    for (index_t jp = rowp[i]; jp < rowp[i + 1]; jp++) {
      if constexpr (is_complex<T>::value) {
        std::fprintf(fp, "%d %d %30.20e %30.20e\n", i + 1, cols[jp] + 1,
                     vals[jp].real(), vals[jp].imag());
        nnz_mtx++;
      } else if (absfunc(vals[jp]) >= epsilon) {
        std::fprintf(fp, "%d %d %30.20e\n", i + 1, cols[jp] + 1, vals[jp]);
        nnz_mtx++;
      }
    }
  }
  std::fclose(fp);

  // Modify nnz
  fp = std::fopen(mtx_name.c_str(), "r+");
  if (is_complex<T>::value) {
    std::fprintf(fp, "%%%%MatrixMarket matrix coordinate complex general\n");
  } else {
    std::fprintf(fp, "%%%%MatrixMarket matrix coordinate real general\n");
  }
  std::fprintf(fp, "%15d%15d%15d", nrows, ncols, nnz_mtx);
  std::fclose(fp);
  return;
}

/**
 * @brief Zero out columns and set diagonal entry to one for each zeroed
 * column
 *
 * @param nbcs number of global columns to zero-out
 * @param dof global dof indices
 */
template <typename T>
void CSCMat<T>::zero_columns(const index_t nbcs, const index_t dof[]) {
  for (index_t ii = 0; ii < nbcs; ii++) {
    index_t column = dof[ii];

    for (index_t jp = colp[column]; jp < colp[column + 1]; jp++) {
      vals[jp] = 0.0;

      if (rows[jp] == column) {
        vals[jp] = 1.0;
      }
    }
  }
}

/**
 * @brief Convert to a dense matrix
 *
 * @param m_ output, number of rows of the dense matrix
 * @param n_ output, number of columns of the dense matrix
 * @param A_ output, dense matrix values stored in row-major ordering
 */
template <typename T>
void CSCMat<T>::to_dense(index_t *m_, index_t *n_, T **A_) {
  index_t m = nrows;
  index_t n = ncols;
  index_t size = m * n;

  T *A = new T[size];
  std::fill(A, A + size, T(0.0));

  for (index_t j = 0; j < ncols; j++) {
    for (index_t ip = colp[j]; ip < colp[j + 1]; ip++) {
      index_t i = rows[ip];
      A[n * i + j] = vals[ip];
    }
  }

  *A_ = A;
  *m_ = m;
  *n_ = n;
}

/**
 * @brief Export the matrix as mtx format
 *
 * @param mtx_name the output file
 * @param epsilon only write entries such that abs(e) >= epsilon
 */
template <typename T>
void CSCMat<T>::write_mtx(const std::string mtx_name, double epsilon) {
  // Open file and destroy old contents, if any
  std::FILE *fp = std::fopen(mtx_name.c_str(), "w");

  // Write header
  if (is_complex<T>::value) {
    std::fprintf(fp, "%%%%MatrixMarket matrix coordinate complex general\n");
  } else {
    std::fprintf(fp, "%%%%MatrixMarket matrix coordinate real general\n");
  }

  // Write m, n and nnz
  std::fprintf(fp, "%15d%15d%15d\n", nrows, ncols, nnz);

  // Write entries
  index_t nnz_mtx = 0;
  for (index_t j = 0; j < ncols; j++) {
    for (index_t ip = colp[j]; ip < colp[j + 1]; ip++) {
      if constexpr (is_complex<T>::value) {
        std::fprintf(fp, "%d %d %30.20e %30.20e\n", rows[ip] + 1, j + 1,
                     vals[ip].real(), vals[ip].imag());
        nnz_mtx++;
      } else if (absfunc(vals[ip]) >= epsilon) {
        std::fprintf(fp, "%d %d %30.20e\n", rows[ip] + 1, j + 1, vals[ip]);
        nnz_mtx++;
      }
    }
  }
  std::fclose(fp);

  // Modify nnz
  fp = std::fopen(mtx_name.c_str(), "r+");
  if (is_complex<T>::value) {
    std::fprintf(fp, "%%%%MatrixMarket matrix coordinate complex general\n");
  } else {
    std::fprintf(fp, "%%%%MatrixMarket matrix coordinate real general\n");
  }
  std::fprintf(fp, "%15d%15d%15d", nrows, ncols, nnz_mtx);
  std::fclose(fp);
  return;
}

}  // namespace SparseUtils

#endif  // SPARSE_UTILS_SPARSE_MATRIX_IMPL_H