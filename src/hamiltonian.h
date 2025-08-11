// File: src/hamiltonian.h
// Purpose: Declares the data structures and functions for loading and managing
//          the Hamiltonian from the input file.

#ifndef HAMILTONIAN_H
#define HAMILTONIAN_H

#include "datatypes.h"

// --- Hamiltonian Data Structure ---
// This struct encapsulates all data loaded from the pmqmc.in file.
typedef struct {
    // Permutation operators and their cycles
    bitset_t** P_matrix;
    bitset_t** cycles;

    // Diagonal Term (D0)
    int D0_size;
    complex_t* D0_coeff;
    bitset_t** D0_product;

    // Off-Diagonal Terms (D_k)
    int* D_sizes;
    complex_t** D_coeffs;
    bitset_t*** D_products;
} Hamiltonian;

// --- Public API Functions ---

/**
 * @brief Creates and loads a Hamiltonian struct from the input file.
 *        This function allocates all necessary memory for the Hamiltonian.
 * @param filename The path to the input file (e.g., "pmqmc.in").
 * @param params A pointer to the SimParams struct to be filled.
 * @return A pointer to the newly created Hamiltonian, or NULL on failure.
 */
Hamiltonian* hamiltonian_create_and_load(const char* filename, SimParams* params);

/**
 * @brief Frees all memory that was dynamically allocated by hamiltonian_load().
 *        Should be called at the end of the program.
 * @param h A pointer to the Hamiltonian to be freed.
 * @param params A pointer to the SimParams struct.
 */
void hamiltonian_free(Hamiltonian* h, const SimParams* params);

#endif // HAMILTONIAN_H