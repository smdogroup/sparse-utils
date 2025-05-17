#ifndef SPARSE_UTILS_SPARSE_SYMBOLIC_H
#define SPARSE_UTILS_SPARSE_SYMBOLIC_H

#include <algorithm>
#include <limits>
#include <set>
#include <vector>

#include "sparse_amd.h"
#include "sparse_matrix.h"
// #include "utils/a2dprofiler.h"
#include "a2dprofiler.h"

namespace SparseUtils {

/**
 * @brief A POD datatype for sparse matrix coordinates
 */
struct COO {
  index_t row_idx;
  index_t col_idx;

  inline bool operator==(const COO& src) const {
    return (row_idx == src.row_idx && col_idx == src.col_idx);
  }
};

/*
  Given the CSR data, sort each row
*/
template <class ArrayType>
void SortCSRData(index_t nrows, ArrayType& rowp, ArrayType& cols) {
  // Sort the cols array
  index_t* it1;
  index_t* it2 = cols.data();
  for (index_t i = 0; i < nrows; i++) {
    it1 = it2;
    it2 += rowp[i + 1] - rowp[i];
    std::sort(it1, it2);  // Note: Kokkos::sort is slow!
  }
}

/*
  Add the connectivity from a connectivity list, use Kokkos unordered set
*/
#if 0
template <class ConnArray>
void BSRMatAddConnectivity(ConnArray& conn,
                           Kokkos::UnorderedMap<COO, void>& node_set) {
  Timer t("BSRMatAddConnectivity()");
  int fail = 0;
  index_t nelems = conn.extent(0);
  index_t nnodes = conn.extent(1);
  node_set.rehash(nelems * nnodes * nnodes);
  Kokkos::parallel_reduce(
      nelems,
      KOKKOS_LAMBDA(const index_t i, int& error) {
        for (index_t j1 = 0; j1 < nnodes; j1++) {
          for (index_t j2 = 0; j2 < nnodes; j2++) {
            auto result = node_set.insert(COO{conn(i, j1), conn(i, j2)});
            error += result.failed();
          }
        }
      },
      fail);
  if (fail) {
    char msg[256];
    std::snprintf(msg, sizeof(msg),
                  "BSRMatAddConnectivity failed with %d failing inserts.",
                  fail);
    throw std::runtime_error(msg);
  }
}

/*
  Create a BSRMat from the set of node pairs, use Kokkos unordered set
*/
template <typename T, index_t M>
BSRMat<T, M, M>* BSRMatFromNodeSet(index_t nnodes,
                                   Kokkos::UnorderedMap<COO, void>& node_set) {
  Timer t("BSRMatFromNodeSet(2)");
  // Find the number of nodes referenced by other nodes
  IdxArray1D_t rowp("rowp", (nnodes + 1));

  // Loop over COO entries and count number of entries in each row
  Kokkos::parallel_for(
      node_set.capacity(), KOKKOS_LAMBDA(index_t i) {
        if (node_set.valid_at(i)) {
          auto key = node_set.key_at(i);
          Kokkos::atomic_increment(&rowp[key.row_idx + 1]);
        }
      });
  Kokkos::fence();

  // Set the pointer into the rows
  rowp[0] = 0;
  for (index_t i = 0; i < nnodes; i++) {
    rowp[i + 1] += rowp[i];
  }

  index_t nnz = rowp[nnodes];
  IdxArray1D_t cols("cols", nnz);

  // Maintain a hash map to track the offset from rowp for each row:
  // offset_tracker[row_idx] = current_offset
  Kokkos::UnorderedMap<index_t, index_t> offset_tracker(nnodes);

  // Set current_offet to 0 for all rows
  int fail = 0;
  Kokkos::parallel_reduce(
      nnodes,
      KOKKOS_LAMBDA(const index_t i, int& error) {
        auto result = offset_tracker.insert(i, index_t(0));
        error += result.failed();
      },
      fail);
  if (fail) {
    char msg[256];
    std::snprintf(msg, sizeof(msg),
                  "Populating offset_tracker failed with %d failing inserts.",
                  fail);
    throw std::runtime_error(msg);
  }

  // Loop over the nodes to increment the tracker and populate cols
  Kokkos::parallel_for(
      node_set.capacity(), KOKKOS_LAMBDA(index_t i) {
        if (node_set.valid_at(i)) {
          // Get row index
          COO node = node_set.key_at(i);
          index_t row_idx = node.row_idx;
          index_t col_idx = node.col_idx;

          // If row index already in the hash map, increment its count,
          // otherwise insert the row index to hash map
          index_t offset =
              Kokkos::atomic_fetch_add(&offset_tracker.value_at(row_idx), 1);
          cols[rowp[row_idx] + offset] = col_idx;
        }
      });
  Kokkos::fence();

  // Sort the cols array
  SortCSRData(nnodes, rowp, cols);

  BSRMat<T, M, M>* A = new BSRMat<T, M, M>(nnodes, nnodes, nnz, rowp, cols);

  return A;
}
#endif

#if 0
/*
  Compute the non-zero pattern of the matrix based on the connectivity pattern
*/
template <typename T, index_t M, class ConnArray>
BSRMat<T, M, M>* BSRMatFromConnectivity(ConnArray& conn) {
  // Set the number of elements
  index_t nelems = conn.extent(0);

  // Find the number of nodes
  index_t nnodes = 0;
  for (index_t i = 0; i < conn.extent(0); i++) {
    for (index_t j = 0; j < conn.extent(1); j++) {
      if (conn(i, j) > nnodes) {
        nnodes = conn(i, j);
      }
    }
  }
  nnodes++;

  // Insert all the nodes into the node set
  std::set<std::pair<index_t, index_t>> node_set;
  for (index_t i = 0; i < nelems; i++) {
    for (index_t j1 = 0; j1 < conn.extent(1); j1++) {
      for (index_t j2 = 0; j2 < conn.extent(1); j2++) {
        node_set.insert(std::pair<index_t, index_t>(conn(i, j1), conn(i, j2)));
      }
    }
  }

  // Find the number of nodes referenced by other nodes
  std::vector<index_t> rowp(nnodes + 1);

  typename std::set<std::pair<index_t, index_t>>::iterator it;
  for (it = node_set.begin(); it != node_set.end(); it++) {
    rowp[it->first + 1] += 1;
  }

  // Set the pointer into the rows
  rowp[0] = 0;
  for (index_t i = 0; i < nnodes; i++) {
    rowp[i + 1] += rowp[i];
  }

  index_t nnz = rowp[nnodes];
  std::vector<index_t> cols(nnz);

  for (it = node_set.begin(); it != node_set.end(); it++) {
    cols[rowp[it->first]] = it->second;
    rowp[it->first]++;
  }

  // Reset the pointer into the nodes
  for (index_t i = nnodes; i > 0; i--) {
    rowp[i] = rowp[i - 1];
  }
  rowp[0] = 0;

  // Sort the cols array
  SortCSRData(nnodes, rowp, cols);

  BSRMat<T, M, M>* A =
      new BSRMat<T, M, M>(nnodes, nnodes, nnz, rowp, cols);

  return A;
}
#endif

/*
  Compute the non-zero pattern of the matrix based on the connectivity pattern
*/
template <typename T, index_t M>
BSRMat<T, M, M>* BSRMatFromConnectivityCUDA(int nelems, int nnodes,
                                            int nodes_per_elem,
                                            const int32_t* conn) {
  // Insert all the nodes into the node set
  std::set<std::pair<index_t, index_t>> node_set;
  for (index_t i = 0; i < nelems; i++) {
    for (index_t j1 = 0; j1 < nodes_per_elem; j1++) {
      for (index_t j2 = 0; j2 < nodes_per_elem; j2++) {
        node_set.insert(std::pair<index_t, index_t>(
            conn[nodes_per_elem * i + j1], conn[nodes_per_elem * i + j2]));
      }
    }
  }

  // Find the number of nodes referenced by other nodes
  std::vector<index_t> rowp(nnodes + 1);

  typename std::set<std::pair<index_t, index_t>>::iterator it;
  for (it = node_set.begin(); it != node_set.end(); it++) {
    rowp[it->first + 1] += 1;
  }

  // Set the pointer into the rows
  rowp[0] = 0;
  for (index_t i = 0; i < nnodes; i++) {
    rowp[i + 1] += rowp[i];
  }

  index_t nnz = rowp[nnodes];
  std::vector<index_t> cols(nnz);

  for (it = node_set.begin(); it != node_set.end(); it++) {
    cols[rowp[it->first]] = it->second;
    rowp[it->first]++;
  }

  // Reset the pointer into the nodes
  for (index_t i = nnodes; i > 0; i--) {
    rowp[i] = rowp[i - 1];
  }
  rowp[0] = 0;

  // Sort the cols array
  SortCSRData(nnodes, rowp, cols);

  BSRMat<T, M, M>* A =
      new BSRMat<T, M, M>(nnodes, nnodes, nnz, rowp.data(), cols.data());

  return A;
}

#if 0
template <typename T, index_t M, class ConnArray>
BSRMat<T, M, M>* BSRMatFromConnectivityDeprecated(ConnArray& conn) {
  // Set the number of elements
  index_t nelems = conn.extent(0);

  // Find the number of nodes
  index_t nnodes = 0;
  for (index_t i = 0; i < conn.extent(0); i++) {
    for (index_t j = 0; j < conn.extent(1); j++) {
      if (conn(i, j) > nnodes) {
        nnodes = conn(i, j);
      }
    }
  }
  nnodes++;

  // Create data to store node -> element connectivity
  std::vector<index_t> node_to_elem_ptr(nnodes + 1);
  for (index_t i = 0; i < conn.extent(0); i++) {
    for (index_t j = 0; j < conn.extent(1); j++) {
      index_t node = conn(i, j);
      node_to_elem_ptr[node + 1]++;
    }
  }

  for (index_t i = 0; i < nnodes; i++) {
    node_to_elem_ptr[i + 1] += node_to_elem_ptr[i];
  }

  std::vector<index_t> node_to_elem(node_to_elem_ptr[nnodes]);
  for (int i = 0; i < conn.extent(0); i++) {
    for (index_t j = 0; j < conn.extent(1); j++) {
      index_t node = conn(i, j);
      node_to_elem[node_to_elem_ptr[node]] = i;
      node_to_elem_ptr[node]++;
    }
  }

  // Do an in-place sort of the row and column data
  SortCSRData(nnodes, node_to_elem_ptr, node_to_elem, node_to_elem);

  // Reset the element pointer
  std::vector<index_t> rowp(nnodes + 1);
  std::vector<index_t> cols;

  rowp[0] = 0;

  // The set of nodes for each
  std::set<index_t> node_set;  // The set of nodes
  for (index_t i = 0; i < nnodes; i++) {
    for (index_t j = node_to_elem_ptr[i]; j < node_to_elem_ptr[i + 1]; j++) {
      int elem = node_to_elem[j];

      for (index_t k = 0; k < conn.extent(1); k++) {
        node_set.insert(conn(elem, k));
      }
    }

    // Set the rowp indices for the next set
    rowp[i + 1] = rowp[i] + node_set.size();

    // Push the values
    typename std::set<index_t>::iterator it;
    for (it = node_set.begin(); it != node_set.end(); it++) {
      cols.push_back(*it);
    }

    node_set.clear();
  }

  index_t nnz = rowp[nnodes];

  BSRMat<T, M, M>* A =
      new BSRMat<T, M, M>(nnodes, nnodes, nnz, rowp, cols);

  return A;
}
#endif

template <class VecType>
index_t CSRFactorSymbolic(const index_t nrows, const VecType Arowp,
                          const VecType Acols, std::vector<index_t>& rowp,
                          std::vector<index_t>& cols) {
  index_t nnz = 0;

  // Column indices associated with the current row
  std::vector<index_t> rcols(nrows);

  // Row, column and diagonal index data for the new factored matrix
  std::vector<index_t> diag(nrows);
  // printf("begin CSRFactorSymbolic\n");

  rowp[0] = 0;
  for (index_t i = 0; i < nrows; i++) {
    index_t nr = 0;  // Number of entries in the current row

    // printf("i %d\n", i);

    // Add the matrix elements to the current row of the matrix.
    // These new elements are sorted.
    for (index_t jp = Arowp[i]; jp < Arowp[i + 1]; jp++) {
      rcols[nr] = Acols[jp];
      nr++;
    }

    // Now, perform the symbolic factorization-- this generates new entries
    // Loop over entries in this row, before the diagonal
    index_t j = 0;
    for (; rcols[j] < i; j++) {
      index_t p = j + 1;                    // The index into rcols
      index_t kp_end = rowp[rcols[j] + 1];  // the end of row number cols[j]

      // Start with the first entry after the diagonal in row, cols[j]
      // k is the index into cols for row cols[j]
      for (index_t kp = diag[rcols[j]] + 1; kp < kp_end; kp++) {
        // Increment p to an entry where we may have cols[k] == rcols[p]
        while (p < nr && rcols[p] < cols[kp]) {
          p++;
        }

        // Add the element into the list of new entries
        if (p >= nr || rcols[p] != cols[kp]) {
          // Insert the new entry into the list and keep the list sorted
          for (index_t n = nr; n > p; n--) {
            rcols[n] = rcols[n - 1];
          }
          rcols[p] = cols[kp];
          nr++;
        }
      }
    }

    // Make sure that we don't exceed the capacity of the vector
    if (nnz + nr > cols.size()) {
      cols.resize(2 * nnz + nr);
    }

    // Add the new elements
    for (index_t n = 0; n < nr; n++, nnz++) {
      cols[nnz] = rcols[n];
    }

    rowp[i + 1] = nnz;
    diag[i] = rowp[i] + j;
  }

  // printf("finish CSRFactorSymbolic\n");

  return nnz;
}

/*
  Find the reordering to reduce the fill in during factorization
*/
// need to add this back for re-ordering version and fix a2d-multiphysics code
// Sean : see BSRMatAMDFactorSymbolicCUDA below for corrections there are
// several mistakes in this code and it reorders the matrix wrong

#if 0
template <typename T, index_t M> BSRMat<T, M, M>*
BSRMatAMDFactorSymbolic(BSRMat<T, M, M>& A,
                                         double fill_factor = 5.0) {
  // Copy over the non-zero structure of the matrix
  int nrows = A.nbrows;
  IdxArray1D_t rowp("rowp", A.nbrows + 1);
  IdxArray1D_t cols("cols", A.nnz);
  IdxArray1D_t perm_("perm_", A.nbrows);

  // Copy the values to rowp and cols
  BLAS::copy(rowp, A.rowp);
  BLAS::copy(cols, A.cols);

  // Compute the re-ordering
  int* interface_nodes = NULL;
  int ninterface_nodes = 0;
  int ndep_vars = 0;
  int* dep_vars = NULL;
  int* indep_ptr = NULL;
  int* indep_vars = NULL;
  int use_exact_degree = 0;
  amd_order_interface(nrows, (int*)rowp.data(), (int*)cols.data(),
                      (int*)perm_.data(), interface_nodes, ninterface_nodes,
                      ndep_vars, dep_vars, indep_ptr, indep_vars,
                      use_exact_degree);

  // Set up the factorization
  // perm[new var] = old_var
  // iperm[old var] = new var

  // Set the permutation array
  IdxArray1D_t perm("perm", A.nbrows);
  IdxArray1D_t iperm("iperm", A.nbrows);
  BLAS::copy(perm, perm_);

  for (index_t i = 0; i < A.nbrows; i++) {
    iperm[perm[i]] = i;
  }

  // Allocate the new arrays for re-ordering the vector
  IdxArray1D_t Arowp("Arowp", A.nbrows + 1);
  IdxArray1D_t Acols("Acols", A.nnz);

  // Re-order the matrix
  Arowp[0] = 0;
  index_t nnz = 0;
  for (index_t i = 0; i < A.nbrows;
       i++) {  // Loop over the new rows of the matrix
    index_t iold = perm[i];

    // Find the old column numbres and convert them to new ones
    for (index_t jp = A.rowp[iold]; jp < A.rowp[iold + 1]; jp++, nnz++) {
      Acols[nnz] = iperm[A.cols[jp]];
    }

    // After copying, update the size
    Arowp[i + 1] = Arowp[i] + (A.rowp[iold + 1] - A.rowp[iold]);
  }

  // Sort the data for the permuted matrix
  SortCSRData(A.nbrows, Arowp, Acols);

  // Compute the symbolic matrix
  std::vector<index_t> Afrowp(A.nbrows + 1);
  std::vector<index_t> Afcols(index_t(fill_factor * nnz));

  index_t Afnnz = CSRFactorSymbolic(A.nbrows, Arowp, Acols, Afrowp, Afcols);

  BSRMat<T, M, M>* Afactor =
      new BSRMat<T, M, M>(A.nbrows, A.nbrows, Afnnz, Afrowp, Afcols);

  // Set up the non-zero pattern for the new matrix
  Afactor->perm = perm;
  Afactor->iperm = iperm;

  return Afactor;
}
#endif

template <typename T, index_t M>
BSRMat<T, M, M>* BSRMatAMDFactorSymbolicCUDA(BSRMat<T, M, M>& A,
                                             double fill_factor = 5.0) {
  // modified version for CUDA
  // Copy over the non-zero structure of the matrix
  int nrows = A.nbrows;
  // IdxArray1D_t rowp("rowp", A.nbrows + 1);
  // IdxArray1D_t cols("cols", A.nnz);
  // IdxArray1D_t perm_("perm_", A.nbrows);

  std::vector<index_t> rowp(A.nbrows + 1);
  std::vector<index_t> cols(A.nnz);
  std::vector<index_t> perm_(A.nbrows);

  // Copy the values to rowp and cols
  // BLAS::copy(rowp, A.rowp);
  // BLAS::copy(cols, A.cols);
  rowp.assign(A.rowp, A.rowp + A.nbrows + 1);
  cols.assign(A.cols, A.cols + A.nnz);

  // return;

  // Compute the re-ordering
  int* interface_nodes = NULL;
  int ninterface_nodes = 0;
  int ndep_vars = 0;
  int* dep_vars = NULL;
  int* indep_ptr = NULL;
  int* indep_vars = NULL;
  int use_exact_degree = 0;
  amd_order_interface(nrows, (int*)rowp.data(), (int*)cols.data(),
                      (int*)perm_.data(), interface_nodes, ninterface_nodes,
                      ndep_vars, dep_vars, indep_ptr, indep_vars,
                      use_exact_degree);

  // Set up the factorization
  // perm[new var] = old_var
  // iperm[old var] = new var

  // Set the permutation array
  std::vector<index_t> perm(A.nbrows);
  std::vector<index_t> iperm(A.nbrows);
  // IdxArray1D_t perm("perm", A.nbrows);
  // IdxArray1D_t iperm("iperm", A.nbrows);
  // BLAS::copy(perm, perm_);
  // std::copy(perm, perm_, A.nbrows);
  perm.assign(perm_.data(), perm_.data() + A.nbrows);

  for (index_t i = 0; i < A.nbrows; i++) {
    iperm[perm[i]] = i;
  }

  // Allocate the new arrays for re-ordering the vector
  // IdxArray1D_t Arowp("Arowp", A.nbrows + 1);
  // IdxArray1D_t Acols("Acols", A.nnz);
  std::vector<index_t> Arowp(A.nbrows + 1);
  std::vector<index_t> Acols(A.nnz);

  // Re-order the matrix
  Arowp[0] = 0;
  index_t nnz = 0;
  // // Sean - had at one point switched to perm on cols and rows (opposite of
  // convention here)
  // // leads to higher fillin (so switched back below)
  // for (index_t inew = 0; inew < A.nbrows;
  //      inew++) {  // Loop over the new rows of the matrix
  //   index_t iold = iperm[inew];

  //   // Find the old column numbres and convert them to new ones
  //   for (index_t jp = A.rowp[iold]; jp < A.rowp[iold + 1]; jp++, nnz++) {
  //     Acols[nnz] = perm[A.cols[jp]];
  //   }

  //   // After copying, update the size
  //   Arowp[inew + 1] = Arowp[inew] + (A.rowp[iold + 1] - A.rowp[iold]);
  // }

  // iperm is used to reorder matrix here essentially
  for (index_t i = 0; i < A.nbrows;
       i++) {  // Loop over the new rows of the matrix
    index_t iold = perm[i];

    // Find the old column numbres and convert them to new ones
    for (index_t jp = A.rowp[iold]; jp < A.rowp[iold + 1]; jp++, nnz++) {
      Acols[nnz] = iperm[A.cols[jp]];
    }

    // After copying, update the size
    Arowp[i + 1] = Arowp[i] + (A.rowp[iold + 1] - A.rowp[iold]);
  }

  // Sort the data for the permuted matrix
  SortCSRData(A.nbrows, Arowp, Acols);

  // Compute the symbolic matrix
  std::vector<index_t> Afrowp(A.nbrows + 1);
  std::vector<index_t> Afcols(index_t(fill_factor * nnz));

  index_t Afnnz = CSRFactorSymbolic(A.nbrows, Arowp, Acols, Afrowp, Afcols);

  // deep copy rowp, cols, perm, iperm (otherwise std::vector go out of scope
  // and data)
  // index_t *new_rowp, *new_cols
  index_t *new_perm, *new_iperm;
  // new_rowp = new index_t[A.nbrows + 1];
  // std::copy(Afrowp.begin(), Afrowp.end(), new_rowp);
  // new_cols = new index_t[Afnnz];
  // std::copy(Afcols.begin(), Afcols.end(), new_cols);
  new_perm = new index_t[A.nbrows];
  std::copy(perm.begin(), perm.end(), new_perm);
  new_iperm = new index_t[A.nbrows];
  std::copy(iperm.begin(), iperm.end(), new_iperm);

  BSRMat<T, M, M>* Afactor = new BSRMat<T, M, M>(A.nbrows, A.nbrows, Afnnz,
                                                 Afrowp.data(), Afcols.data());

  // Set up the non-zero pattern for the new matrix
  // deep copy
  Afactor->perm = new_perm;
  Afactor->iperm = new_iperm;

  return Afactor;
}

template <typename T, index_t M>
BSRMat<T, M, M>* BSRMatReorderSymbolicCUDA(BSRMat<T, M, M>& A,
                                           double fill_factor = 5.0) {
  // modified version for CUDA
  // Copy over the non-zero structure of the matrix
  int nrows = A.nbrows;

  std::vector<index_t> rowp(A.nbrows + 1);
  std::vector<index_t> cols(A.nnz);

  // Copy the values to rowp and cols
  rowp.assign(A.rowp, A.rowp + A.nbrows + 1);
  cols.assign(A.cols, A.cols + A.nnz);

  // Set the permutation array
  std::vector<index_t> perm(A.nbrows);
  std::vector<index_t> iperm(A.nbrows);
  perm.assign(A.perm, A.perm + A.nbrows);
  iperm.assign(A.iperm, A.iperm + A.nbrows);

  // Allocate the new arrays for re-ordering the vector
  std::vector<index_t> Arowp(A.nbrows + 1);
  std::vector<index_t> Acols(A.nnz);

  // Re-order the matrix
  Arowp[0] = 0;
  index_t nnz = 0;

  // iperm is used to reorder matrix here essentially
  for (index_t i = 0; i < A.nbrows;
       i++) {  // Loop over the new rows of the matrix
    index_t iold = perm[i];

    // Find the old column numbres and convert them to new ones
    for (index_t jp = A.rowp[iold]; jp < A.rowp[iold + 1]; jp++, nnz++) {
      Acols[nnz] = iperm[A.cols[jp]];
    }

    // After copying, update the size
    Arowp[i + 1] = Arowp[i] + (A.rowp[iold + 1] - A.rowp[iold]);
  }

  // Sort the data for the permuted matrix
  SortCSRData(A.nbrows, Arowp, Acols);

  // Compute the symbolic matrix
  std::vector<index_t> Afrowp(A.nbrows + 1);
  std::vector<index_t> Afcols(index_t(fill_factor * nnz));

  index_t Afnnz = CSRFactorSymbolic(A.nbrows, Arowp, Acols, Afrowp, Afcols);

  BSRMat<T, M, M>* Afactor = new BSRMat<T, M, M>(A.nbrows, A.nbrows, Afnnz,
                                                 Afrowp.data(), Afcols.data());

  return Afactor;
}

template <typename T, int M>
BSRMat<T, M, M>* BSRMatApplyPerm(BSRMat<T, M, M>& A, int* perm, int* iperm) {
  // accept perm, iperm from separate input so not destroyed when matrix is
  // deleted

  std::vector<index_t> Arowp(A.nbrows + 1);
  std::vector<index_t> Acols(A.nnz);

  // Re-order the matrix
  Arowp[0] = 0;
  index_t nnz = 0;

  // iperm is used to reorder matrix here essentially
  for (index_t i = 0; i < A.nbrows;
       i++) {  // Loop over the new rows of the matrix
    index_t iold = perm[i];

    // Find the old column numbres and convert them to new ones
    for (index_t jp = A.rowp[iold]; jp < A.rowp[iold + 1]; jp++, nnz++) {
      Acols[nnz] = iperm[A.cols[jp]];
    }

    // After copying, update the size
    Arowp[i + 1] = Arowp[i] + (A.rowp[iold + 1] - A.rowp[iold]);
  }

  // Sort the data for the permuted matrix
  SortCSRData(A.nbrows, Arowp, Acols);

  // printf("Method 1 rowp: ");
  // printVec<int>(A.nbrows + 1, Arowp.data());

  //
  BSRMat<T, M, M>* Afactor = new BSRMat<T, M, M>(A.nbrows, A.nbrows, A.nnz,
                                                 Arowp.data(), Acols.data());
  return Afactor;
}

/*
  Symbolic factorization stage
*/
template <typename T, int M>
BSRMat<T, M, M>* BSRMatFactorSymbolic(BSRMat<T, M, M>& A,
                                      double fill_factor = 5.0) {
  std::vector<index_t> rowp(A.nbrows + 1);
  std::vector<index_t> cols(index_t(fill_factor * A.nnz));

  index_t nnz = CSRFactorSymbolic(A.nbrows, A.rowp, A.cols, rowp, cols);

  // printf("Method 2 rowp: ");
  // printVec<int>(A.nbrows + 1, rowp.data());

  BSRMat<T, M, M>* Afactor = new BSRMat<T, M, M>(
      A.nbrows, A.nbrows, nnz, rowp.data(), cols.data(), nullptr);

  return Afactor;
}

template <typename T, int M>
BSRMat<T, M, M>* BSRMatReorderFactorSymbolic(BSRMat<T, M, M>& A, int* perm,
                                             double fill) {
  // printf("in BSRMatReorderFactorSymbolic:");
  // deep copy perm so it doesn't get modified or deleted later
  int* _perm = new int[A.nbrows];
  int* iperm = new int[A.nbrows];
  for (int inode = 0; inode < A.nbrows; inode++) {
    _perm[inode] = perm[inode];
    iperm[perm[inode]] = inode;
  }
  // printf("checkpt1\n");

  // don't want to put these in the matrix, will get deleted later
  // A.perm = _perm;
  // A.iperm = iperm;

  auto A2 = BSRMatApplyPerm<T, M>(A, perm, iperm);
  // printf("checkpt2\n");
  return BSRMatFactorSymbolic<T, M>(*A2, fill);
}

#if 0
/*
  Compute the non-zero pattern for C = A * B
*/
template <typename T, index_t M, index_t N, index_t P>
BSRMat<T, M, P>* BSRMatMatMultSymbolic(BSRMat<T, M, N>& A, BSRMat<T, N, P>& B,
                                       double fill_factor = 2.0) {
  index_t nrows = A.nbrows;
  index_t ncols = B.nbcols;
  index_t nnz = 0;

  // Column indices associated with the current row
  const index_t empty = MAX_INDEX;
  const index_t first_entry = MAX_INDEX - 1;
  std::vector<index_t> next(ncols, empty);

  // Row, column and diagonal index data for the new factored matrix
  IdxArray1D_t rowp("rowp", nrows + 1);
  IdxArray1D_t cols("cols", index_t(fill_factor * A.nnz));

  // Compute the non-zero structure of the resulting matrix C = A * B
  // one row at a time
  for (index_t i = 0; i < A.nbrows; i++) {
    index_t head = first_entry;
    index_t num_cols = 0;  // The size of the temporary cols array

    // Add the non-zero pattern to this matrix from each row of B
    // for each column of A.
    index_t jp_end = A.rowp[i + 1];
    for (index_t jp = A.rowp[i]; jp < jp_end; jp++) {
      index_t j = A.cols[jp];

      // Merge the two arrays into cols
      index_t kp_end = B.rowp[j + 1];
      for (index_t kp = B.rowp[j]; kp < kp_end; kp++) {
        index_t k = B.cols[kp];
        if (next[k] == empty) {
          next[k] = head;
          head = k;
          num_cols++;
        }
      }
    }

    if (nnz + num_cols > cols.size()) {
      Kokkos::resize(cols, nnz + num_cols + A.nnz);
    }

    // Reverse through the list
    for (index_t j = 0; j < num_cols; j++) {
      cols[nnz] = head;
      nnz++;
      index_t temp = head;
      head = next[head];
      next[temp] = empty;
    }

    rowp[i + 1] = nnz;
  }

  // Sort the CSR data
  SortCSRData(nrows, rowp, cols);

  BSRMat<T, M, P>* bsr = new BSRMat<T, M, P>(nrows, ncols, nnz, rowp, cols);

  return bsr;
}
#endif

/*
  Compute the non-zero pattern for C = S + A * B
*/
#if 0
template <typename T, index_t M, index_t N, index_t P>
BSRMat<T, M, P>* BSRMatMatMultAddSymbolic(BSRMat<T, M, P>& S,
                                          BSRMat<T, M, N>& A,
                                          BSRMat<T, N, P>& B,
                                          double fill_factor = 2.0) {
  index_t nrows = A.nbrows;
  index_t ncols = B.nbcols;
  index_t nnz = 0;

  // Column indices associated with the current row
  const index_t empty = MAX_INDEX;
  const index_t first_entry = MAX_INDEX - 1;
  std::vector<index_t> next(ncols, empty);

  // Row, column and diagonal index data for the new factored matrix
  IdxArray1D_t rowp("rowp", nrows + 1);
  IdxArray1D_t cols("cols", index_t(fill_factor * A.nnz));

  // Compute the non-zero structure of the resulting matrix C = A * B
  // one row at a time
  for (index_t i = 0; i < A.nbrows; i++) {
    int head = int(first_entry);
    index_t num_cols = 0;  // The size of the temporary cols array

    // Merge the two arrays into cols
    index_t jp_end = S.rowp[i + 1];
    for (index_t jp = S.rowp[i]; jp < jp_end; jp++) {
      index_t j = S.cols[jp];
      if (next[j] == empty) {
        next[j] = head;
        head = j;
        num_cols++;
      }
    }

    // Add the non-zero pattern to this matrix from each row of B
    // for each column of A.
    jp_end = A.rowp[i + 1];
    for (index_t jp = A.rowp[i]; jp < jp_end; jp++) {
      index_t j = A.cols[jp];

      index_t kp_end = B.rowp[j + 1];
      for (index_t kp = B.rowp[j]; kp < kp_end; kp++) {
        index_t k = B.cols[kp];
        if (next[k] == empty) {
          next[k] = head;
          head = k;
          num_cols++;
        }
      }
    }

    if (nnz + num_cols > cols.size()) {
      Kokkos::resize(cols, nnz + num_cols + A.nnz);
    }

    // Reverse through the list
    for (index_t j = 0; j < num_cols; j++) {
      cols[nnz] = head;
      nnz++;
      index_t temp = head;
      head = next[head];
      next[temp] = empty;
    }

    rowp[i + 1] = nnz;
  }

  // Sort the CSR data
  SortCSRData(nrows, rowp, cols);

  BSRMat<T, M, P>* bsr = new BSRMat<T, M, P>(nrows, ncols, nnz, rowp, cols);

  return bsr;
}
#endif

/*
  Compute the non-zero pattern of the transpose of the matrix
*/
template <typename T, index_t M, index_t N>
BSRMat<T, N, M>* BSRMatMakeTransposeSymbolic(BSRMat<T, M, N>& A) {
  // The number of rows and columns for the transposed matrix
  index_t nrows = A.nbcols;
  index_t ncols = A.nbrows;

  std::vector<index_t> rowp(nrows + 1, 0);

  // Count up the number of references
  for (index_t i = 0; i < A.nbrows; i++) {
    for (index_t jp = A.rowp[i]; jp < A.rowp[i + 1]; jp++) {
      index_t j = A.cols[jp];
      rowp[j + 1]++;
    }
  }

  for (index_t i = 0; i < nrows; i++) {
    rowp[i + 1] += rowp[i];
  }

  index_t nnz = rowp[nrows];
  std::vector<index_t> cols(nnz);
  for (index_t i = 0; i < A.nbrows; i++) {
    for (index_t jp = A.rowp[i]; jp < A.rowp[i + 1]; jp++) {
      index_t j = A.cols[jp];
      cols[rowp[j]] = i;
      rowp[j]++;
    }
  }

  // Re-set the rowp array
  for (index_t i = nrows; i > 0; i--) {
    rowp[i] = rowp[i - 1];
  }
  rowp[0] = 0;

  // Create the new BSR matrix
  BSRMat<T, N, M>* At =
      new BSRMat<T, N, M>(nrows, ncols, nnz, rowp.data(), cols.data());

  return At;
}

/*
  Make a transpose matrix
*/
template <typename T, index_t M, index_t N>
BSRMat<T, N, M>* BSRMatMakeTranspose(BSRMat<T, M, N>& A) {
  BSRMat<T, N, M>* At = BSRMatMakeTransposeSymbolic(A);
  At->zero();

  // Loop over the values in A
  for (index_t i = 0; i < A.nbrows; i++) {
    for (index_t jp = A.rowp[i]; jp < A.rowp[i + 1]; jp++) {
      index_t j = A.cols[jp];

      index_t kp = At->find_value_index(j, i);  // Find At(j, i)
      if (kp != NO_INDEX) {
        for (index_t k1 = 0; k1 < M; k1++) {
          for (index_t k2 = 0; k2 < N; k2++) {
            At->vals(kp, k2, k1) = A.vals(jp, k1, k2);
          }
        }
      } else {
        std::cerr
            << "BSRMatMakeTranspose: Non-zero pattern does not match cannot "
               "copy values"
            << std::endl;
      }
    }
  }

  return At;
}

/*
  Duplicate the non-zero pattern of A without copying the values
*/
template <typename T, index_t M, index_t N>
BSRMat<T, M, N>* BSRMatDuplicate(BSRMat<T, M, N>& A) {
  return new BSRMat<T, M, N>(A.nbrows, A.nbcols, A.nnz, A.rowp, A.cols);
}

/*
  Multicolor code for a single process using a greedy algorithm
*/
#if 0
template <class VecType>
index_t CSRMultiColorOrder(const index_t nvars, const index_t rowp[],
                           const index_t cols[], VecType colors, VecType perm) {
  // Allocate a temporary array to store the
  const index_t empty = MAX_INDEX;
  std::vector<index_t> tmp(nvars + 1);
  std::fill(tmp.begin(), tmp.begin() + nvars, empty);
  for (index_t i = 0; i < nvars; i++) {
    colors[i] = empty;
  }

  index_t num_colors = 0;
  for (index_t i = 0; i < nvars; i++) {
    // Find the minimum color that is not referred to by any adjacent
    // node.
    const index_t jp_end = rowp[i + 1];
    for (index_t jp = rowp[i]; jp < jp_end; jp++) {
      int j = cols[jp];
      if (colors[j] != empty) {
        tmp[colors[j]] = i;
      }
    }

    // Set the color for this variable if it already exists
    bool flag = true;
    for (index_t k = 0; k < num_colors; k++) {
      if (tmp[k] != i) {
        colors[i] = k;
        flag = false;
        break;
      }
    }

    // Create a new color
    if (flag) {
      colors[i] = num_colors;
      num_colors++;
    }
  }

  // Now that all the nodes have been colored, order them
  std::fill(tmp.begin(), tmp.begin() + num_colors + 1, 0);

  // Count up the number of nodes for each color
  for (index_t i = 0; i < nvars; i++) {
    tmp[colors[i] + 1]++;
  }

  // Set tmp as an offset for each color
  for (int i = 1; i < num_colors + 1; i++) {
    tmp[i] += tmp[i - 1];
  }

  // Create the new color variables
  for (int i = 0; i < nvars; i++) {
    perm[tmp[colors[i]]] = i;
    tmp[colors[i]]++;
  }

  return num_colors;
}

template <typename T, index_t M>
void BSRMatMultiColorOrder(BSRMat<T, M, M>& A) {
  A.perm = IdxArray1D_t("A.perm", A.nbrows);

  IdxArray1D_t colors("colors", A.nbrows);

  // Use A->perm as a temporary variable to avoid double allocation
  A.num_colors = CSRMultiColorOrder(A.nbrows, A.rowp.data(), A.cols.data(),
                                    colors.data(), A.perm.data());

  // Count up the number of nodes with each color
  A.color_count = IdxArray1D_t("A.color_count", A.num_colors);

  for (index_t i = 0; i < A.nbrows; i++) {
    A.color_count[colors[i]]++;
  }
}
#endif

}  // namespace SparseUtils

#endif  // SPARSE_UTILS_SPARSE_SYMBOLIC_H
