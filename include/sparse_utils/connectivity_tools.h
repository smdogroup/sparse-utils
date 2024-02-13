#ifndef SPARSE_UTILS_CONNECTIVITY_TOOLS_H
#define SPARSE_UTILS_CONNECTIVITY_TOOLS_H

namespace SparseUtils {

void NodeToElementFromConnectivity(int num_nodes, int num_elements,
                                   const int nodes_per_element,
                                   const int* element_nodes,
                                   int** node_to_elem_ptr_,
                                   int** node_to_elem_) {
  // Create data to store node -> element connectivity
  int* node_to_elem_ptr = new int[num_nodes + 1];
  for (int i = 0; i < num_nodes + 1; i++) {
    node_to_elem_ptr[i] = 0;
  }

  const int* ptr = element_nodes;
  for (int i = 0; i < num_elements; i++) {
    for (int j = 0; j < nodes_per_element; j++, ptr++) {
      node_to_elem_ptr[ptr[0] + 1]++;
    }
  }

  for (int i = 0; i < num_nodes; i++) {
    node_to_elem_ptr[i + 1] += node_to_elem_ptr[i];
  }

  // Set up the node to element data
  int* node_to_elem = new int[node_to_elem_ptr[num_nodes]];
  ptr = element_nodes;
  for (int i = 0; i < num_elements; i++) {
    for (int j = 0; j < nodes_per_element; j++, ptr++) {
      int node = ptr[0];
      node_to_elem[node_to_elem_ptr[node]] = i;
      node_to_elem_ptr[node]++;
    }
  }

  for (int i = num_nodes; i > 0; i--) {
    node_to_elem_ptr[i] = node_to_elem_ptr[i - 1];
  }
  node_to_elem_ptr[0] = 0;

  // Set the outputs
  *node_to_elem_ptr_ = node_to_elem_ptr;
  *node_to_elem_ = node_to_elem;
}

void CSRFromConnectivity(int num_nodes, int num_elements,
                         const int nodes_per_element, const int* element_nodes,
                         int** rowp_, int** cols_) {
  int* node_to_elem_ptr = nullptr;
  int* node_to_elem = nullptr;
  NodeToElementFromConnectivity(num_nodes, num_elements, nodes_per_element,
                                element_nodes, &node_to_elem_ptr,
                                &node_to_elem);

  // Set up the CSR data structure
  int* rowp = new int[num_nodes + 1];
  int* counter = new int[num_nodes];

  // Initialize the counter
  for (int i = 0; i < num_nodes; i++) {
    counter[i] = -1;
  }

  // Count up the number of non-zero entries
  rowp[0] = 0;
  for (int i = 0, nnz = 0; i < num_nodes; i++) {
    for (int j = node_to_elem_ptr[i]; j < node_to_elem_ptr[i + 1]; j++) {
      int elem = node_to_elem[j];
      for (int k = 0; k < nodes_per_element; k++) {
        int node = element_nodes[elem * nodes_per_element + k];
        if (counter[node] < i) {
          counter[node] = i;
          nnz++;
        }
      }
    }
    rowp[i + 1] = nnz;
  }

  // Allocate the column indices
  int* cols = new int[rowp[num_nodes]];

  // Reset the counter
  for (int i = 0; i < num_nodes; i++) {
    counter[i] = -1;
  }

  // Count up the number of non-zero entries
  for (int i = 0, nnz = 0; i < num_nodes; i++) {
    for (int j = node_to_elem_ptr[i]; j < node_to_elem_ptr[i + 1]; j++) {
      int elem = node_to_elem[j];
      for (int k = 0; k < nodes_per_element; k++) {
        int node = element_nodes[elem * nodes_per_element + k];
        if (counter[node] < i) {
          counter[node] = i;
          cols[nnz] = node;
          nnz++;
        }
      }
    }
  }

  *rowp_ = rowp;
  *cols_ = cols;

  // Free unused data
  delete[] node_to_elem_ptr;
  delete[] node_to_elem;
  delete[] counter;
}

}  // namespace SparseUtils

#endif  // SPARSE_UTILS_CONNECTIVITY_TOOLS_H
