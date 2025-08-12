// File: src/datatypes.c
// Purpose: Implements the logic for the dynamic bitset defined in datatypes.h.

#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>
#include "datatypes.h"

// --- Bitset Function Implementations ---

bitset_t* bitset_create(int num_bits) {
    if (num_bits <= 0) {
        fprintf(stderr, "Error: bitset_create requires num_bits > 0\n");
        return NULL;
    }

    bitset_t* bs = (bitset_t*)malloc(sizeof(bitset_t));
    if (!bs) {
        perror("Failed to allocate bitset struct");
        return NULL;
    }

    bs->num_bits = num_bits;
    // Calculate bytes needed: ceil(num_bits / 8.0)
    bs->num_bytes = (num_bits + 7) / 8; 
    
    // Use calloc to initialize all bits to zero automatically
    bs->bytes = (unsigned char*)calloc(bs->num_bytes, sizeof(unsigned char));
    if (!bs->bytes) {
        perror("Failed to allocate bytes for bitset");
        free(bs);
        return NULL;
    }

    return bs;
}

void bitset_free(bitset_t* bs) {
    if (bs) {
        free(bs->bytes); // Free the internal array first
        free(bs);        // Then free the struct itself
    }
}

bitset_t* bitset_create_from_string(const char* str) {
    int len = strlen(str);
    bitset_t* bs = bitset_create(len);
    if (!bs) return NULL;

    for (int i = 0; i < len; ++i) {
        // The string "1000" corresponds to bit 3 being set.
        // Character at string index 'i' corresponds to bit at position 'len - 1 - i'.
        if (str[i] == '1') {
            bitset_set(bs, len - 1 - i);
        }
    }
    return bs;
}

void bitset_copy(bitset_t* dest, const bitset_t* src) {
    if (dest == NULL || src == NULL || dest->num_bits != src->num_bits) {
        fprintf(stderr, "Error: bitset_copy requires valid, same-sized bitsets.\n");
        return;
    }
    memcpy(dest->bytes, src->bytes, src->num_bytes);
}

void bitset_set(bitset_t* bs, int pos) {
    if (pos < 0 || pos >= bs->num_bits) return;
    int byte_index = pos / 8;
    int bit_index_in_byte = pos % 8;
    bs->bytes[byte_index] |= (1 << bit_index_in_byte);
}

void bitset_flip(bitset_t* bs, int pos) {
    if (pos < 0 || pos >= bs->num_bits) return;
    int byte_index = pos / 8;
    int bit_index_in_byte = pos % 8;
    bs->bytes[byte_index] ^= (1 << bit_index_in_byte);
}

int bitset_get(const bitset_t* bs, int pos) {
    if (pos < 0 || pos >= bs->num_bits) return 0;
    int byte_index = pos / 8;
    int bit_index_in_byte = pos % 8;
    return (bs->bytes[byte_index] >> bit_index_in_byte) & 1;
}

void bitset_xor(bitset_t* dest, const bitset_t* src) {
    if (dest == NULL || src == NULL || dest->num_bits != src->num_bits) {
        fprintf(stderr, "Error: bitset_xor requires valid, same-sized bitsets.\n");
        return;
    }
    for (int i = 0; i < dest->num_bytes; ++i) {
        dest->bytes[i] ^= src->bytes[i];
    }
}

int bitset_count(const bitset_t* bs) {
    int count = 0;
    int i = 0;
    // Process in 64-bit chunks for speed
    for (; i + 7 < bs->num_bytes; i += 8) {
        uint64_t chunk;
        memcpy(&chunk, bs->bytes + i, 8);
        count += __builtin_popcountll(chunk); // GCC/Clang intrinsic
    }
    // Process remaining bytes
    for (; i < bs->num_bytes; ++i) {
        count += __builtin_popcount(bs->bytes[i]);
    }
    return count;
}

void bitset_print(const bitset_t* bs) {
    if (!bs) {
        printf("(null)\n");
        return;
    }
    for (int i = bs->num_bits - 1; i >= 0; --i) {
        printf("%d", bitset_get(bs, i));
    }
}