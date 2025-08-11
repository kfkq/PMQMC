// File: src/utils.c
// Purpose: Implements the utility functions.

#include "utils.h"
#include <stdlib.h> // For srand, rand, RAND_MAX

void rng_init(unsigned int seed) {
    srand(seed);
}

double rng_uniform_double() {
    // Standard method to get a double in [0.0, 1.0)
    // The +1.0 ensures the upper bound is exclusive.
    return (double)rand() / ((double)RAND_MAX + 1.0);
}

int rng_uniform_int(int max) {
    if (max <= 0) {
        return 0; // Safety check
    }
    // This method avoids modulo bias, which can be an issue with rand() % max
    return (int)(rng_uniform_double() * max);
}