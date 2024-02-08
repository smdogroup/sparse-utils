
#ifndef SPARSE_UTILS_SPARSE_MATRIX_IMPL_H
#define SPARSE_UTILS_SPARSE_MATRIX_IMPL_H

#include "sparse_matrix.h"

namespace SparseUtils {

// Zero the entries of the sparse matrix
template <typename T, int M, int N>
void BSRMat<T, M, N>::zero() {
  std::fill(vals, vals + M * N * nnz, T(0.0));
}

/**
 * @brief Find the value index of a block given indices (row, col) of a block
 *
 * @param row block row index
 * @param col block column index
 * @return int the value index jp such that vals[jp] gives the
 * sought block, if (row, col) isn't in the nonzero pattern, NO_INDEX will be
 * returned
 */
template <typename T, int M, int N>
int BSRMat<T, M, N>::find_value_index(int row, int col) {
  int jp_start = rowp[row];
  int jp_end = rowp[row + 1];

  for (int jp = jp_start; jp < jp_end; jp++) {
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
template <typename T, int M, int N>
template <class Mat>
void BSRMat<T, M, N>::add_values(const int m, const int i[],
                                 const int n, const int j[], Mat &mat) {
  for (int ii = 0; ii < m; ii++) {
    int block_row = i[ii] / M;
    int local_row = i[ii] % M;

    for (int jj = 0; jj < n; jj++) {
      int block_col = j[jj] / N;
      int local_col = j[jj] % N;

      int jp = find_value_index(block_row, block_col);
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
template <typename T, int M, int N>
template <class Mat>
void BSRMat<T, M, N>::add_block_values(const int m, const int ib[],
                                       const int n, const int jb[],
                                       Mat &mat) {
  for (int ii = 0; ii < m; ii++) {
    int block_row = ib[ii];

    for (int jj = 0; jj < n; jj++) {
      int block_col = jb[jj];

      int jp = find_value_index(block_row, block_col);
      if (jp != NO_INDEX) {
        for (int local_row = 0; local_row < M; local_row++) {
          for (int local_col = 0; local_col < N; local_col++) {
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
template <typename T, int M, int N>
void BSRMat<T, M, N>::zero_rows(const int nbcs, const int dof[]) {
  for (int ii = 0; ii < nbcs; ii++) {
    int block_row = dof[ii] / M;
    int local_row = dof[ii] % M;

    for (int jp = rowp[block_row]; jp < rowp[block_row + 1]; jp++) {
      for (int k = 0; k < N; k++) {
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
template <typename T, int M, int N>
void BSRMat<T, M, N>::to_dense(int *m_, int *n_, T **A_) {
  int m = M * nbrows;
  int n = N * nbcols;
  int size = m * n;

  T *A = new T[size];
  std::fill(A, A + size, T(0.0));

  for (int i = 0; i < nbrows; i++) {
    for (int jp = rowp[i]; jp < rowp[i + 1]; jp++) {
      int j = cols[jp];

      for (int ii = 0; ii < M; ii++) {
        const int irow = M * i + ii;
        for (int jj = 0; jj < N; jj++) {
          const int jcol = N * j + jj;
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
template <typename T, int M, int N>
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
  int nnz_mtx = 0;
  for (int i = 0; i < nbrows; i++) {
    for (int jp = rowp[i]; jp < rowp[i + 1]; jp++) {
      int j = cols[jp];  // (i, j) is the block index pair

      for (int ii = 0; ii < M; ii++) {
        const int irow = M * i + ii + 1;  // convert to 1-based index
        for (int jj = 0; jj < N; jj++) {
          // (irow, jcol) is the entry coo
          const int jcol = N * j + jj + 1;  // convert to 1-based index
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
void CSRMat<T>::to_dense(int *m_, int *n_, T **A_) {
  int m = nrows;
  int n = ncols;
  int size = m * n;

  T *A = new T[size];
  std::fill(A, A + size, T(0.0));

  for (int i = 0; i < nrows; i++) {
    for (int jp = rowp[i]; jp < rowp[i + 1]; jp++) {
      int j = cols[jp];
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
  int nnz_mtx = 0;
  for (int i = 0; i < nrows; i++) {
    for (int jp = rowp[i]; jp < rowp[i + 1]; jp++) {
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
void CSCMat<T>::zero_columns(const int nbcs, const int dof[]) {
  for (int ii = 0; ii < nbcs; ii++) {
    int column = dof[ii];

    for (int jp = colp[column]; jp < colp[column + 1]; jp++) {
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
void CSCMat<T>::to_dense(int *m_, int *n_, T **A_) {
  int m = nrows;
  int n = ncols;
  int size = m * n;

  T *A = new T[size];
  std::fill(A, A + size, T(0.0));

  for (int j = 0; j < ncols; j++) {
    for (int ip = colp[j]; ip < colp[j + 1]; ip++) {
      int i = rows[ip];
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
  int nnz_mtx = 0;
  for (int j = 0; j < ncols; j++) {
    for (int ip = colp[j]; ip < colp[j + 1]; ip++) {
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