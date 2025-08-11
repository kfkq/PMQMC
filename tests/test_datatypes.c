// File: src/test_datatypes.c
// Purpose: A standalone program to test the functionality of the dynamic bitset.

#include <stdio.h>
#include <stdlib.h>
#include "datatypes.h"

// A simple macro to make testing easier
#define TEST_ASSERT(condition, message) \
    if (!(condition)) { \
        fprintf(stderr, "TEST FAILED: %s\n", message); \
        exit(1); \
    } else { \
        printf("TEST PASSED: %s\n", message); \
    }

int main() {
    printf("--- Running Tests for Dynamic Bitset ---\n\n");

    // --- Test 1: Creation and Destruction ---
    printf("--- Test 1: Creation and Destruction ---\n");
    bitset_t* bs_small = bitset_create(10);
    TEST_ASSERT(bs_small != NULL, "Small bitset creation (N=10)");
    TEST_ASSERT(bs_small->num_bits == 10, "Small bitset has correct number of bits");
    TEST_ASSERT(bs_small->num_bytes == 2, "Small bitset has correct number of bytes (ceil(10/8)=2)");
    
    bitset_t* bs_large = bitset_create(100);
    TEST_ASSERT(bs_large != NULL, "Large bitset creation (N=100)");
    TEST_ASSERT(bs_large->num_bits == 100, "Large bitset has correct number of bits");
    TEST_ASSERT(bs_large->num_bytes == 13, "Large bitset has correct number of bytes (ceil(100/8)=13)");
    
    bitset_free(bs_small);
    bitset_free(bs_large);
    printf("Creation and destruction tests seem OK.\n\n");

    // --- Test 2: Set, Get, and Flip ---
    printf("--- Test 2: Set, Get, and Flip ---\n");
    bitset_t* bs = bitset_create(70);
    TEST_ASSERT(bitset_get(bs, 5) == 0, "Initial bit is 0");
    
    bitset_set(bs, 5);
    TEST_ASSERT(bitset_get(bs, 5) == 1, "bitset_set works");
    
    bitset_flip(bs, 5);
    TEST_ASSERT(bitset_get(bs, 5) == 0, "bitset_flip (1->0) works");
    
    bitset_flip(bs, 5);
    TEST_ASSERT(bitset_get(bs, 5) == 1, "bitset_flip (0->1) works");

    // Test edge cases
    bitset_set(bs, 0);
    TEST_ASSERT(bitset_get(bs, 0) == 1, "Edge case: set bit 0");
    bitset_set(bs, 69);
    TEST_ASSERT(bitset_get(bs, 69) == 1, "Edge case: set bit 69 (last bit)");
    bitset_free(bs);
    printf("Set, Get, and Flip tests seem OK.\n\n");

    // --- Test 3: Count, Copy, and XOR ---
    printf("--- Test 3: Count, Copy, and XOR ---\n");
    bitset_t* bs1 = bitset_create_from_string("01101"); // N=5, value 13
    bitset_t* bs2 = bitset_create_from_string("10101"); // N=5, value 21
    
    TEST_ASSERT(bitset_count(bs1) == 3, "bitset_count works");
    
    bitset_xor(bs1, bs2); // 13 ^ 21 = 01101 ^ 10101 = 11000 (value 24)
    TEST_ASSERT(bitset_get(bs1, 4) == 1, "XOR result bit 4 is correct");
    TEST_ASSERT(bitset_get(bs1, 3) == 1, "XOR result bit 3 is correct");
    TEST_ASSERT(bitset_get(bs1, 2) == 0, "XOR result bit 2 is correct");
    TEST_ASSERT(bitset_count(bs1) == 2, "XOR result count is correct");

    bitset_t* bs3 = bitset_create(5);
    bitset_copy(bs3, bs1);
    TEST_ASSERT(bitset_get(bs3, 4) == 1 && bitset_get(bs3, 3) == 1, "bitset_copy works");
    
    bitset_free(bs1);
    bitset_free(bs2);
    bitset_free(bs3);
    printf("Count, Copy, and XOR tests seem OK.\n\n");
    
    printf("----------------------------------------\n");
    printf("All datatypes tests passed successfully!\n");
    printf("----------------------------------------\n");

    return 0;
}