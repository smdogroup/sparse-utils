#include <vector>

#include "sparse_utils/sparse_matrix.h"
#include "test_commons.h"

template <typename T, int M, int N>
void test_bsr_axpy(int nbrows, int nbcols, int nnz, std::vector<int>& rowp,
                   std::vector<int>& cols, std::vector<T>& vals,
                   std::vector<T>& b, std::vector<T>& axpy_exact) {
  using BSRMat_t = SparseUtils::BSRMat<double, M, N>;
  BSRMat_t bsr(nbrows, nbcols, nnz, rowp.data(), cols.data(), vals.data());
  std::vector<double> axpy(axpy_exact.size(), 0.0);
  bsr.axpy(b.data(), axpy.data());
  EXPECT_VEC_NEAR(axpy.size(), axpy, axpy_exact, 1e-15);
}

TEST(BSRMatTest, axpy) {
  int nbrows = 8;
  int nbcols = 6;
  int nnz = 30;

  std::vector<int> rowp = {0, 5, 9, 12, 17, 20, 23, 27, 30};
  std::vector<int> cols = {0, 1, 2, 3, 5, 1, 2, 4, 5, 0, 1, 5, 0, 1, 2,
                           3, 5, 1, 3, 4, 1, 3, 5, 0, 1, 2, 3, 0, 2, 3};
  std::vector<double> vals = {
      0.5488135039273248, 0.7151893663724195, 0.6027633760716439,
      0.5448831829968969, 0.6458941130666561, 0.8917730007820798,
      0.9636627605010293, 0.7917250380826646, 0.5288949197529045,
      0.5680445610939323, 0.925596638292661,  0.832619845547938,
      0.7781567509498505, 0.8700121482468192, 0.978618342232764,
      0.7991585642167236, 0.7805291762864555, 0.6399210213275238,
      0.9446689170495839, 0.5218483217500717, 0.7742336894342167,
      0.5684339488686485, 0.6176354970758771, 0.6120957227224214,
      0.6169339968747569, 0.9437480785146242, 0.6818202991034834,
      0.6976311959272649, 0.6667667154456677, 0.6706378696181594};
  std::vector<double> b = {0.3154283509241839, 0.3637107709426226,
                           0.5701967704178796, 0.4386015134623203,
                           0.9883738380592262, 0.1020448107480281};
  std::vector<double> axpy_exact = {1.0818238759379846, 1.7103161359990606,
                                    0.6007913606166261, 1.5500516159330469,
                                    1.162860613401067,  0.5939396197613945,
                                    1.2546274056033493, 0.8943836700535315};
  test_bsr_axpy<double, 1, 1>(nbrows, nbcols, nnz, rowp, cols, vals, b,
                              axpy_exact);
}