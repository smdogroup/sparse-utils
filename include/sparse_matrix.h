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
template <typename T, index_t M, index_t N>
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
  template <class VecType>
  BSRMat(index_t nbrows, index_t nbcols, index_t nnz, const VecType &rowp_,
         const VecType &cols_)
      : nbrows(nbrows), nbcols(nbcols), nnz(nnz) {
    rowp = new index_t[nbrows + 1];
    cols = new index_t[nnz];
    vals = new T[M * N * nnz];

    for (index_t i = 0; i < nbrows + 1; i++) {
      rowp[i] = rowp_[i];
    }

    for (index_t i = 0; i < nnz; i++) {
      cols[i] = cols_[i];
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
  index_t find_value_index(index_t row, index_t col);

  // add values from an element matrix mat of shape (m, n)
  template <class Mat>
  void add_values(const index_t m, const index_t i[], const index_t n,
                  const index_t j[], Mat &mat);

  template <class Mat>
  void add_block_values(const index_t m, const index_t i[], const index_t n,
                        const index_t j[], Mat &mat);

  // Zero out rows and set diagonal entry to one for each zeroed row
  void zero_rows(const index_t nbcs, const index_t dof[]);

  // Convert to a dense matrix
  void to_dense(index_t *m_, index_t *n_, T **A_);

  // Export the matrix as mtx format
  void write_mtx(const std::string mtx_name = "matrix.mtx",
                 double epsilon = 0.0);

  // Number of block rows, block columns and non-zero blocks, nnz = rowp[nbrows]
  index_t nbrows = 0, nbcols = 0, nnz = 0;

  // rowp and cols array
  index_t *rowp = nullptr;  // length: nbrows + 1
  index_t *cols = nullptr;  // length: nnz = rowp[nbrows]
  T *vals = nullptr;

  // Pointer to the diagonal block, this is not allocated until
  // factorization
  index_t *diag = nullptr;  // length: nbrows

  // row-permutation perm[new row] = old row
  // This is not allocated by default
  index_t *perm = nullptr;

  // Inverse row-permutation iperm[old row] = new row
  // This is not allocated by default
  index_t *iperm = nullptr;

  // When coloring is used, its ordering is stored in the permutation array
  index_t num_colors = 0;          // Number of colors
  index_t *color_count = nullptr;  // Number of nodes with this color, not
                                   // allocated by default
};

// Compressed sparse row matrix
template <typename T>
class CSRMat {
 public:
  CSRMat(index_t nrows, index_t ncols, index_t nnz,
         const index_t *rowp_ = nullptr, const index_t *cols_ = nullptr)
      : nrows(nrows), ncols(ncols), nnz(nnz) {
    rowp = new index_t[nrows + 1];
    cols = new index_t[nnz];
    vals = new T[nnz];

    if (rowp_ && cols_) {
      for (index_t i = 0; i < nrows + 1; i++) {
        rowp[i] = rowp_[i];
      }
      for (index_t i = 0; i < nnz; i++) {
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
  void to_dense(index_t *m_, index_t *n_, T **A_);

  // Export the matrix as mtx format
  void write_mtx(const std::string mtx_name = "matrix.mtx",
                 double epsilon = 0.0);

  index_t nrows, ncols, nnz;  // number of rows, columns and nonzeros
  index_t *rowp;              // length: nrows + 1
  index_t *cols;              // length: nnz
  T *vals;                    // length: nnz
};

/**
 * @brief Compressed sparse column matrix
 */
template <typename T>
class CSCMat {
 public:
  CSCMat(index_t nrows, index_t ncols, index_t nnz,
         const index_t *colp_ = nullptr, const index_t *rows_ = nullptr)
      : nrows(nrows), ncols(ncols), nnz(nnz) {
    colp = new index_t[ncols + 1];
    rows = new index_t[nnz];
    vals = new T[nnz];

    if (colp_ && rows_) {
      for (index_t i = 0; i < ncols + 1; i++) {
        colp[i] = colp_[i];
      }
      for (index_t i = 0; i < nnz; i++) {
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
  void zero_columns(const index_t nbcs, const index_t dof[]);

  // Convert to a dense matrix
  void to_dense(index_t *m_, index_t *n_, T **A_);

  // Export the matrix as mtx format
  void write_mtx(const std::string mtx_name = "matrix.mtx",
                 double epsilon = 0.0);

  index_t nrows = 0, ncols = 0,
          nnz = 0;          // number of rows, columns and nonzeros
  index_t *colp = nullptr;  // length: ncols + 1
  index_t *rows = nullptr;  // length: nnz
  T *vals = nullptr;        // length: nnz
};

}  // namespace SparseUtils

#include "detail/sparse_matrix_impl.h"

#endif  // SPARSE_UTILS_SPARSE_MATRIX_H