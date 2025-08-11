// File: src/utils.c
// Purpose: Implements the utility functions using Mersenne Twister PRNG.

#include "utils.h"

// Mersenne Twister constants
#define N 624
#define M 397
#define MATRIX_A 0x9908b0dfUL
#define UPPER_MASK 0x80000000UL
#define LOWER_MASK 0x7fffffffUL

// Mersenne Twister state
static unsigned long mt[N];
static int mti = N + 1;

/**
 * Initializes the Mersenne Twister with a seed.
 */
void rng_init(unsigned int seed) {
    mt[0] = seed & 0xffffffffUL;
    for (mti = 1; mti < N; mti++) {
        mt[mti] = (1812433253UL * (mt[mti - 1] ^ (mt[mti - 1] >> 30)) + mti);
        mt[mti] &= 0xffffffffUL;
    }
}

/**
 * Generates the next N words in the sequence.
 */
static void mt_generate() {
    int i;
    unsigned long y;
    
    for (i = 0; i < N - M; i++) {
        y = (mt[i] & UPPER_MASK) | (mt[i + 1] & LOWER_MASK);
        mt[i] = mt[i + M] ^ (y >> 1) ^ ((y & 1) ? MATRIX_A : 0);
    }
    
    for (; i < N - 1; i++) {
        y = (mt[i] & UPPER_MASK) | (mt[i + 1] & LOWER_MASK);
        mt[i] = mt[i + (M - N)] ^ (y >> 1) ^ ((y & 1) ? MATRIX_A : 0);
    }
    
    y = (mt[N - 1] & UPPER_MASK) | (mt[0] & LOWER_MASK);
    mt[N - 1] = mt[M - 1] ^ (y >> 1) ^ ((y & 1) ? MATRIX_A : 0);
    
    mti = 0;
}

/**
 * Generates a random 32-bit unsigned integer.
 */
static unsigned long mt_rand_int32() {
    unsigned long y;
    
    if (mti >= N) {
        mt_generate();
    }
    
    y = mt[mti++];
    
    // Tempering
    y ^= (y >> 11);
    y ^= (y << 7) & 0x9d2c5680UL;
    y ^= (y << 15) & 0xefc60000UL;
    y ^= (y >> 18);
    
    return y;
}

double rng_uniform_double() {
    // Generate a double in [0.0, 1.0)
    return ((double)mt_rand_int32() / (double)0x100000000ULL);
}

int rng_uniform_int(int max) {
    if (max <= 0) {
        return 0; // Safety check
    }
    // This method avoids modulo bias
    return (int)(rng_uniform_double() * max);
}