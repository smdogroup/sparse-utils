#ifndef SPARSE_UTILS_SPARSE_MATRIX_H
#define SPARSE_UTILS_SPARSE_MATRIX_H

#include <string>

#include "defs.h"

namespace SparseUtils {

/**
 * @brief Block Compressed sparse row matrix
 *
 * @tparam T data type
 * @tparam M number of rows for each block
 * @tparam N number of columns for each block
 *
 * Example:
 *
 * Below is an illustration of a blocked sparse matrix with 2x3 blocks
 *
 *  [x x x          [x x x
 *   x x x]          x x x]
 *          [x x x          [x x x
 *           x x x]          x x x]
 *  [x x x          [x x x
 *   x x x]          x x x]
 *          [x x x          [x x x
 *           x x x]          x x x]
 *
 * for this example, M = 2, N = 3, nbrows = nbcols = 4
 * block row and column indices are 0, 1, 2, 3
 * global row and column indices are 0, 1, ..., 7
 * local row indices for each block are 0, 1
 * local column indices for each block are 0, 1, 2
 *
 * Note: blocks are stored row-by-row.
 *
 */
template <typename T, int M, int N>
class BSRMat {
 public:
  /**
   * @brief Constructor
   *
   * @tparam VecType a vector type whose entry can be []-indexed
   * @param nbrows number of rows of blocks
   * @param nbcols number of columns of blocks
   * @param nnz number of non-zero blocks, note that global nnz = nnz * M * N
   * @param rowp_ vector of row pointers
   * @param cols_ vector of column indices
   */
  BSRMat(int nbrows, int nbcols, int nnz, const int *rowp_, const int *cols_,
         const T *vals_ = nullptr)
      : nbrows(nbrows), nbcols(nbcols), nnz(nnz) {
    rowp = new int[nbrows + 1];
    cols = new int[nnz];
    vals = new T[M * N * nnz];

    for (int i = 0; i < nbrows + 1; i++) {
      rowp[i] = rowp_[i];
    }

    for (int i = 0; i < nnz; i++) {
      cols[i] = cols_[i];
    }

    if (vals_) {
      for (int i = 0; i < M * N * nnz; i++) {
        vals[i] = vals_[i];
      }
    } else {
      zero();
    }
  }

  ~BSRMat() {
    FREE_ARRAY(rowp);
    FREE_ARRAY(cols);
    FREE_ARRAY(vals);
    FREE_ARRAY(diag);
    FREE_ARRAY(perm);
    FREE_ARRAY(iperm);
    FREE_ARRAY(color_count);
  };

  // Zero the entries of the sparse matrix
  void zero();

  // Find the value index of a block given indices (row, col) of a block
  int find_value_index(int row, int col);

  // add values from an element matrix mat of shape (m, n)
  void add_values(const int m, const int i[], const int n, const int j[],
                  T mat[]);
  void add_block_values(const int m, const int i[], const int n, const int j[],
                        T mat[]);

  // Zero out rows and set diagonal entry to one for each zeroed row
  void zero_rows(const int nbcs, const int dof[]);

  // Matrix-vector multiplication
  void axpy(T x[], T y[]);

  // Convert to a dense matrix
  void to_dense(int *m_, int *n_, T **A_);

  // Export the matrix as mtx format
  void write_mtx(const std::string mtx_name = "matrix.mtx",
                 double epsilon = 0.0);

  // Number of block rows, block columns and non-zero blocks, nnz = rowp[nbrows]
  int nbrows = 0, nbcols = 0, nnz = 0;

  // rowp and cols array
  int *rowp = nullptr;  // length: nbrows + 1
  int *cols = nullptr;  // length: nnz = rowp[nbrows]
  T *vals = nullptr;

  // Pointer to the diagonal block, this is not allocated until
  // factorization
  int *diag = nullptr;  // length: nbrows

  // row-permutation perm[new row] = old row
  // This is not allocated by default
  int *perm = nullptr;

  // Inverse row-permutation iperm[old row] = new row
  // This is not allocated by default
  int *iperm = nullptr;

  // When coloring is used, its ordering is stored in the permutation array
  int num_colors = 0;          // Number of colors
  int *color_count = nullptr;  // Number of nodes with this color, not
                               // allocated by default
};

// Compressed sparse row matrix
template <typename T>
class CSRMat {
 public:
  CSRMat(int nrows, int ncols, int nnz, const int *rowp_ = nullptr,
         const int *cols_ = nullptr)
      : nrows(nrows), ncols(ncols), nnz(nnz) {
    rowp = new int[nrows + 1];
    cols = new int[nnz];
    vals = new T[nnz];

    if (rowp_ && cols_) {
      for (int i = 0; i < nrows + 1; i++) {
        rowp[i] = rowp_[i];
      }
      for (int i = 0; i < nnz; i++) {
        cols[i] = cols_[i];
      }
    }
  }

  ~CSRMat() {
    FREE_ARRAY(rowp);
    FREE_ARRAY(cols);
    FREE_ARRAY(vals);
  }

  // Convert to a dense matrix
  void to_dense(int *m_, int *n_, T **A_);

  // Export the matrix as mtx format
  void write_mtx(const std::string mtx_name = "matrix.mtx",
                 double epsilon = 0.0);

  int nrows, ncols, nnz;  // number of rows, columns and nonzeros
  int *rowp;              // length: nrows + 1
  int *cols;              // length: nnz
  T *vals;                // length: nnz
};

/**
 * @brief Compressed sparse column matrix
 */
template <typename T>
class CSCMat {
 public:
  CSCMat(int nrows, int ncols, int nnz, const int *colp_ = nullptr,
         const int *rows_ = nullptr)
      : nrows(nrows), ncols(ncols), nnz(nnz) {
    colp = new int[ncols + 1];
    rows = new int[nnz];
    vals = new T[nnz];

    if (colp_ && rows_) {
      for (int i = 0; i < ncols + 1; i++) {
        colp[i] = colp_[i];
      }
      for (int i = 0; i < nnz; i++) {
        rows[i] = rows_[i];
      }
    }
  }

  ~CSCMat() {
    FREE_ARRAY(colp);
    FREE_ARRAY(rows);
    FREE_ARRAY(vals);
  }

  // Zero out columns and set diagonal entry to one for each zeroed column
  void zero_columns(const int nbcs, const int dof[]);

  // Convert to a dense matrix
  void to_dense(int *m_, int *n_, T **A_);

  // Export the matrix as mtx format
  void write_mtx(const std::string mtx_name = "matrix.mtx",
                 double epsilon = 0.0);

  int nrows = 0, ncols = 0,
      nnz = 0;          // number of rows, columns and nonzeros
  int *colp = nullptr;  // length: ncols + 1
  int *rows = nullptr;  // length: nnz
  T *vals = nullptr;    // length: nnz
};

}  // namespace SparseUtils

#include "detail/sparse_matrix_impl.h"

#endif  // SPARSE_UTILS_SPARSE_MATRIX_H