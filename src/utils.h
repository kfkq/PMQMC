// File: src/utils.h
// Purpose: Declares general-purpose utility functions, primarily for random
//          number generation needed by the simulation.

#ifndef UTILS_H
#define UTILS_H

/**
 * @brief Initializes the random number generator with a specific seed.
 *        This should be called once at the very beginning of the program.
 * @param seed The seed for the RNG. Using the same seed will produce the
 *             same sequence of random numbers, which is useful for debugging.
 */
void rng_init(unsigned int seed);

/**
 * @brief Generates a random double uniformly distributed in the range [0.0, 1.0).
 * @return A random double.
 */
double rng_uniform_double();

/**
 * @brief Generates a random integer uniformly distributed in the range [0, max-1].
 * @param max The exclusive upper bound for the random integer.
 * @return A random integer.
 */
int rng_uniform_int(int max);

#endif // UTILS_H