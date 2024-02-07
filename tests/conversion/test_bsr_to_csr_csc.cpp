#include <unordered_set>
#include <vector>

#include "defs.h"
#include "sparse_matrix.h"
#include "test_commons.h"
#include "utils.h"

using namespace SparseUtils;

class BSRMatTest : public ::testing::Test {
 protected:
  static index_t constexpr M = 2;
  static index_t constexpr N = 2;
  using BSRMat_t = BSRMat<double, M, N>;
  using CSRMat_t = CSRMat<double>;
  using CSCMat_t = CSCMat<double>;

  void SetUp() override {
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

    bsr_mat = new BSRMat_t(nbrows, nbcols, nnz, rowp, cols);

    for (index_t n = 0; n < nnz; n++) {
      for (index_t ii = 0; ii < M; ii++) {
        for (index_t jj = 0; jj < N; jj++) {
          bsr_mat->vals[M * N * n + N * ii + jj] =
              (double)rand() / RAND_MAX;  // [0, 1]
        }
      }
    }
  }

  BSRMat_t *bsr_mat = nullptr;
};

TEST_F(BSRMatTest, BSR_TO_CSR) {
  CSRMat_t *csr_mat = bsr_to_csr(bsr_mat);
  index_t m_bsr, n_bsr;
  double *vals_bsr;
  bsr_mat->to_dense(&m_bsr, &n_bsr, &vals_bsr);

  index_t m_csr, n_csr;
  double *vals_csr;
  csr_mat->to_dense(&m_csr, &n_csr, &vals_csr);

  EXPECT_EQ(m_bsr, m_csr);
  EXPECT_EQ(n_bsr, n_csr);
  EXPECT_VEC_EQ(m_bsr * n_bsr, vals_bsr, vals_csr);
}

TEST_F(BSRMatTest, BSR_TO_CSC) {
  CSCMat_t *csc_mat = bsr_to_csc(bsr_mat);

  index_t m_bsr, n_bsr;
  double *vals_bsr;
  bsr_mat->to_dense(&m_bsr, &n_bsr, &vals_bsr);

  index_t m_csc, n_csc;
  double *vals_csc;
  csc_mat->to_dense(&m_csc, &n_csc, &vals_csc);

  EXPECT_EQ(m_bsr, m_csc);
  EXPECT_EQ(n_bsr, n_csc);
  EXPECT_VEC_EQ(m_bsr * n_bsr, vals_bsr, vals_csc);
}
