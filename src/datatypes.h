// File: src/datatypes.h
// Purpose: Defines fundamental, general-purpose data types and the bitset API
//          used throughout the QMC simulation.

#ifndef DATATYPES_H
#define DATATYPES_H

#include <complex.h> // For standard C complex numbers

// --- Dynamic Bitset Type ---
// A flexible bitset that can handle an arbitrary number of bits (N),
// determined at runtime.
typedef struct {
    unsigned char* bytes; // Dynamically allocated array to store bits
    int num_bits;         // The total number of bits (e.g., N)
    int num_bytes;        // The allocated size of the bytes array
} bitset_t;

// --- Configuration Parameters ---
// This struct will hold all parameters read from pmqmc.in.
// It is populated once during initialization.
typedef struct {
    int N;
    int NOP;
    int NCYCLES;
    double BETA;
    long long TSTEPS;
    long long STEPS;
    int STEPS_PER_MEASUREMENT;
    long long SKIP_MEASUREMENTS;
    int QMAX;
    int NBINS;
    int WORM;  // Boolean: 1 for worm updates, 0 for composite updates
} SimParams;


// --- Worm Data Structure ---
// This struct holds the state of the worm algorithm.
typedef struct {
    int active;      // Is the worm currently active?
    int i, j;        // Positions of the worm's head and tail
    int k, l;        // Operator indices at the head and tail
    bitset_t* z_k;   // Spin configuration at the head
    bitset_t* z_l;   // Spin configuration at the tail
} Worm;


// --- Standard Complex Number Type ---
// A convenient typedef for C99's double complex type.
typedef double complex complex_t;


// --- Bitset Function Prototypes (Public API) ---

/**
 * @brief Creates and allocates a new bitset, initialized to all zeros.
 * @param num_bits The number of bits the bitset should hold (N).
 * @return A pointer to the newly created bitset, or NULL on failure.
 */
bitset_t* bitset_create(int num_bits);

/**
 * @brief Creates a bitset from a string of '0's and '1's.
 * @param str A null-terminated string (e.g., "00101100").
 * @return A pointer to the newly created bitset, or NULL on failure.
 */
bitset_t* bitset_create_from_string(const char* str);

/**
 * @brief Frees all memory associated with a bitset.
 * @param bs A pointer to the bitset to be freed.
 */
void bitset_free(bitset_t* bs);

/**
 * @brief Copies the content of a source bitset to a destination bitset.
 * @param dest The destination bitset (must be pre-allocated and same size as src).
 * @param src The source bitset.
 */
void bitset_copy(bitset_t* dest, const bitset_t* src);

/**
 * @brief Sets the bit at a given position to 1.
 * @param bs The bitset to modify.
 * @param pos The 0-indexed position of the bit to set.
 */
void bitset_set(bitset_t* bs, int pos);

/**
 * @brief Flips the bit at a given position (0->1, 1->0).
 * @param bs The bitset to modify.
 * @param pos The 0-indexed position of the bit to flip.
 */
void bitset_flip(bitset_t* bs, int pos);

/**
 * @brief Gets the value of the bit (0 or 1) at a given position.
 * @param bs The bitset to read from.
 * @param pos The 0-indexed position of the bit to get.
 * @return The integer value (0 or 1) of the bit.
 */
int bitset_get(const bitset_t* bs, int pos);

/**
 * @brief Performs a bitwise XOR operation (dest ^= src).
 * @param dest The bitset that will be modified.
 * @param src The bitset to XOR with.
 */
void bitset_xor(bitset_t* dest, const bitset_t* src);

/**
 * @brief Performs a bitwise AND operation (dest &= src).
 * @param dest The bitset that will be modified.
 * @param src The bitset to AND with.
 */
void bitset_and(bitset_t* dest, const bitset_t* src);

/**
 * @brief Counts the number of set bits (1s) in the bitset.
 * @param bs The bitset to count.
 * @return The total number of set bits.
 */
int bitset_count(const bitset_t* bs);

/**
 * @brief Prints a string representation of the bitset to stdout.
 * @param bs The bitset to print.
 */
void bitset_print(const bitset_t* bs);


#endif // DATATYPES_H