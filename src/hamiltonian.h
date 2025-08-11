// File: src/hamiltonian.h
// Purpose: Declares the data structures and functions for loading and managing
//          the Hamiltonian from the input file.

#ifndef HAMILTONIAN_H
#define HAMILTONIAN_H

#include "datatypes.h"

// --- Global Hamiltonian Data ---
// These pointers will be allocated and filled by hamiltonian_load().
// The 'extern' keyword means they are defined in hamiltonian.c but
// can be accessed by any file that includes this header.

extern bitset_t** P_matrix; // Array of pointers to bitsets
extern bitset_t** cycles;   // Array of pointers to bitsets

// Diagonal Term (D0)
extern complex_t* D0_coeff;
extern bitset_t** D0_product;
extern int D0_size;

// Off-Diagonal Terms (D_k)
extern int* D_sizes;
extern complex_t** D_coeffs;
extern bitset_t*** D_products; // 3D array: [op_idx][term_idx][bitset_ptr]


// --- Public API Functions ---

/**
 * @brief Loads all simulation parameters and Hamiltonian data from the input file.
 *        This function allocates all necessary memory for the global arrays.
 * @param filename The path to the input file (e.g., "pmqmc.in").
 * @param params A pointer to the SimParams struct to be filled.
 * @return 0 on success, -1 on failure.
 */
int hamiltonian_load(const char* filename, SimParams* params);

/**
 * @brief Frees all memory that was dynamically allocated by hamiltonian_load().
 *        Should be called at the end of the program.
 */
void hamiltonian_free(const SimParams* params);

#endif // HAMILTONIAN_H