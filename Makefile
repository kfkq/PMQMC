# File: Makefile
# Description: Builds the PM-QMC simulator

# Compiler and Flags
CC = mpicc
# TODO: Edit these paths based on your system's HDF5 installation.
# Example for Debian/Ubuntu:
HDF5_INCLUDE_DIR = /usr/include/hdf5/serial
HDF5_LIB_DIR     = /usr/lib/x86_64-linux-gnu/hdf5/serial
# Example for macOS with Homebrew:
# HDF5_INCLUDE_DIR = /opt/homebrew/include
# HDF5_LIB_DIR     = /opt/homebrew/lib
# ------------------------------------
# CFLAGS: -Wall (all warnings), -Wextra (more warnings), -g (debug symbols), 
# -std=c99 (language standard), -Isrc (include path for headers)
CFLAGS = -Wall -Wextra -O3 -std=c99 -Isrc -I$(HDF5_INCLUDE_DIR)
# VPATH tells make where to find source files from src/ and tests/
VPATH = src:tests
# Linker Flags for HDF5 and MPI ---
LDFLAGS = -L$(HDF5_LIB_DIR) -lhdf5 -lm

# =============================================================================
# --- FILE DEFINITIONS ---
# =============================================================================

# --- Main Simulation Program ---
MAIN_SRCS = main.c datatypes.c divdiff.c hamiltonian.c state.c updates.c utils.c measurements.c io.c processing.c
MAIN_OBJS = $(MAIN_SRCS:.c=.o)
MAIN_EXEC = pmqmc

# --- Unit Test Programs ---
TEST_DATATYPES_SRCS = test_datatypes.c datatypes.c
TEST_DATATYPES_OBJS = $(TEST_DATATYPES_SRCS:.c=.o)
TEST_DATATYPES_EXEC = test_datatypes

TEST_DIVDIFF_SRCS = test_divdiff.c divdiff.c datatypes.c
TEST_DIVDIFF_OBJS = $(TEST_DIVDIFF_SRCS:.c=.o)
TEST_DIVDIFF_EXEC = test_divdiff

TEST_HAMILTONIAN_SRCS = test_hamiltonian.c hamiltonian.c datatypes.c
TEST_HAMILTONIAN_OBJS = $(TEST_HAMILTONIAN_SRCS:.c=.o)
TEST_HAMILTONIAN_EXEC = test_hamiltonian

TEST_STATE_SRCS = test_state.c state.c hamiltonian.c divdiff.c datatypes.c utils.c
TEST_STATE_OBJS = $(TEST_STATE_SRCS:.c=.o)
TEST_STATE_EXEC = test_state

TEST_UTILS_SRCS = test_utils.c utils.c
TEST_UTILS_OBJS = $(TEST_UTILS_SRCS:.c=.o)
TEST_UTILS_EXEC = test_utils

TEST_UPDATES_SRCS = test_updates.c updates.c state.c hamiltonian.c divdiff.c datatypes.c utils.c
TEST_UPDATES_OBJS = $(TEST_UPDATES_SRCS:.c=.o)
TEST_UPDATES_EXEC = test_updates


# =============================================================================
# --- TARGETS ---
# =============================================================================

# Default target: Build the main program and the analyzer
all: $(MAIN_EXEC) $(ANALYZER_EXEC)

# --- Main Executable Build Rules ---
$(MAIN_EXEC): $(MAIN_OBJS)
	$(CC) $(CFLAGS) -o $@ $(MAIN_OBJS) -lm $(LDFLAGS)

# --- Unified Test Target ---
# Builds and runs all individual unit tests.
test: $(TEST_DATATYPES_EXEC) $(TEST_DIVDIFF_EXEC) $(TEST_HAMILTONIAN_EXEC) $(TEST_STATE_EXEC) $(TEST_UTILS_EXEC) $(TEST_UPDATES_EXEC)
	@echo "\n========================================"
	@echo "         RUNNING ALL UNIT TESTS"
	@echo "========================================"
	./$(TEST_DATATYPES_EXEC)
	./$(TEST_DIVDIFF_EXEC)
	./$(TEST_HAMILTONIAN_EXEC)
	./$(TEST_STATE_EXEC)
	./$(TEST_UTILS_EXEC)
	./$(TEST_UPDATES_EXEC)
	@echo "\n========================================"
	@echo "      ALL UNIT TESTS PASSED"
	@echo "========================================"

# --- Individual Test Build Rules ---
$(TEST_DATATYPES_EXEC): $(TEST_DATATYPES_OBJS)
	$(CC) $(CFLAGS) -o $@ $(TEST_DATATYPES_OBJS)

$(TEST_DIVDIFF_EXEC): $(TEST_DIVDIFF_OBJS)
	$(CC) $(CFLAGS) -o $@ $(TEST_DIVDIFF_OBJS) -lm

$(TEST_HAMILTONIAN_EXEC): $(TEST_HAMILTONIAN_OBJS)
	$(CC) $(CFLAGS) -o $@ $(TEST_HAMILTONIAN_OBJS) -lm

$(TEST_STATE_EXEC): $(TEST_STATE_OBJS)
	$(CC) $(CFLAGS) -o $@ $(TEST_STATE_OBJS) -lm

$(TEST_UTILS_EXEC): $(TEST_UTILS_OBJS)
	$(CC) $(CFLAGS) -o $@ $(TEST_UTILS_OBJS) -lm

$(TEST_UPDATES_EXEC): $(TEST_UPDATES_OBJS)
	$(CC) $(CFLAGS) -o $@ $(TEST_UPDATES_OBJS) -lm


# --- Generic and Cleanup Rules ---

# Generic rule to compile .c files from src/ or tests/ into .o object files
%.o: %.c
	$(CC) $(CFLAGS) -c $< -o $@

# Clean up all build files and generated data
clean:
	rm -f *.o
	rm -f $(MAIN_EXEC)
	rm -f $(TEST_DATATYPES_EXEC) $(TEST_DIVDIFF_EXEC) $(TEST_HAMILTONIAN_EXEC) $(TEST_STATE_EXEC) $(TEST_UTILS_EXEC) $(TEST_UPDATES_EXEC)
	rm -f raw_data.h5 results.dat hamiltonian.in

# Declare targets that are not files
.PHONY: all clean test analyze