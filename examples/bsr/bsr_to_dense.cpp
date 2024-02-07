#include <iostream>
#include <unordered_set>

#include "defs.h"
#include "sparse_matrix.h"
#include "utils.h"

using namespace SparseUtils;

static index_t constexpr M = 2;
static index_t constexpr N = 2;
using BSRMat_t = BSRMat<double, M, N>;
using CSRMat_t = CSRMat<double>;
using CSCMat_t = CSCMat<double>;

int main() {
  index_t constexpr nbrows = 3;
  index_t constexpr nbcols = 3;

  srand(0);

  std::vector<index_t> rowp(nbrows + 1);
  for (index_t i = 0; i < nbrows; i++) {
    rowp[i] = rand() % (nbcols + 1);  // [0, nbrows]
  }

  index_t presum = 0, temp = 0;
  for (index_t i = 0; i < nbrows; i++) {
    temp = rowp[i];
    rowp[i] = presum;
    presum += temp;
  }
  index_t nnz = presum;
  rowp[nbrows] = nnz;

  printf("nnz: %ld\n", nnz);

  index_t index = 0;
  std::vector<index_t> cols(nnz);
  for (index_t i = 0; i < nbrows; i++) {
    int num = rowp[i + 1] - rowp[i];
    std::unordered_set<index_t> s;
    while (s.size() < num) {
      s.insert(rand() % nbcols);  // [0, nbcols - 1]
    }
    for (auto e : s) {
      cols[index] = e;
      index++;
    }
  }

  BSRMat_t* bsr_mat = new BSRMat_t(nbrows, nbcols, nnz, rowp, cols);

  for (index_t n = 0; n < nnz; n++) {
    for (index_t ii = 0; ii < M; ii++) {
      for (index_t jj = 0; jj < N; jj++) {
        bsr_mat->vals[M * N * n + N * ii + jj] =
            (double)rand() / RAND_MAX;  // [0, 1]
      }
    }
  }

  index_t m_bsr, n_bsr;
  double* vals_bsr;
  bsr_mat->to_dense(&m_bsr, &n_bsr, &vals_bsr);

  printf("sparse matrix:\n");
  printf("rowp: ");
  for (index_t i = 0; i < bsr_mat->nbrows + 1; i++) {
    printf("%6ld  ", bsr_mat->rowp[i]);
  }
  printf("\n");
  printf("cols: ");
  for (index_t i = 0; i < bsr_mat->nnz; i++) {
    printf("%6ld  ", bsr_mat->cols[i]);
  }
  printf("\n");
  printf("vals\n");
  for (index_t i = 0; i < bsr_mat->nnz; i++) {
    for (index_t j = 0; j < M * N; j++) {
      printf("%6.2f  ", bsr_mat->vals[M * N * i + j]);
    }
    printf("\n");
  }

  printf("dense matrix:\n");
  for (index_t i = 0; i < m_bsr; i++) {
    for (index_t j = 0; j < n_bsr; j++) {
      std::printf("%6.2f  ", vals_bsr[i * n_bsr + j]);
    }
    std::printf("\n");
  }

  CSRMat_t* csr_mat = bsr_to_csr(bsr_mat);

  return 0;
}