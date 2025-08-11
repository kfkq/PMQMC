# File: Makefile

# Compiler and Flags
CC = gcc
CFLAGS = -Wall -Wextra -g -std=c99 -Isrc
VPATH = src:tests

# --- Main Program Files (currently empty until we add main.c) ---
MAIN_SRCS = 
MAIN_OBJS = $(MAIN_SRCS:.c=.o)
MAIN_EXEC = pmrqmc

# --- Test Programs ---
TEST_DATATYPES_SRCS = test_datatypes.c datatypes.c
TEST_DATATYPES_OBJS = $(TEST_DATATYPES_SRCS:.c=.o)
TEST_DATATYPES_EXEC = test_datatypes

TEST_DIVDIFF_SRCS = test_divdiff.c divdiff.c datatypes.c
TEST_DIVDIFF_OBJS = $(TEST_DIVDIFF_SRCS:.c=.o)
TEST_DIVDIFF_EXEC = test_divdiff

TEST_HAMILTONIAN_SRCS = test_hamiltonian.c hamiltonian.c datatypes.c
TEST_HAMILTONIAN_OBJS = $(TEST_HAMILTONIAN_SRCS:.c=.o)
TEST_HAMILTONIAN_EXEC = test_hamiltonian

TEST_STATE_SRCS = test_state.c state.c hamiltonian.c divdiff.c datatypes.c
TEST_STATE_OBJS = $(TEST_STATE_SRCS:.c=.o)
TEST_STATE_EXEC = test_state

TEST_UTILS_SRCS = test_utils.c utils.c
TEST_UTILS_OBJS = $(TEST_UTILS_SRCS:.c=.o)
TEST_UTILS_EXEC = test_utils

# NEW TEST FOR UPDATES
TEST_UPDATES_SRCS = test_updates.c updates.c state.c hamiltonian.c divdiff.c datatypes.c utils.c
TEST_UPDATES_OBJS = $(TEST_UPDATES_SRCS:.c=.o)
TEST_UPDATES_EXEC = test_updates

# =============================================================================
# --- TARGETS ---
# =============================================================================

# The 'all' target is currently empty. To build the final program later:
# make pmrqmc
all:
	@echo "No main program to build yet. Run 'make test' to check modules."

$(MAIN_EXEC): $(MAIN_OBJS)
	$(CC) $(CFLAGS) -o $(MAIN_EXEC) $(MAIN_OBJS) -lm

# --- Unified Test Target ---
test: $(TEST_DATATYPES_EXEC) $(TEST_DIVDIFF_EXEC) $(TEST_HAMILTONIAN_EXEC) $(TEST_STATE_EXEC) $(TEST_UTILS_EXEC) $(TEST_UPDATES_EXEC)
	@echo "\n========================================"
	@echo "         RUNNING ALL TESTS"
	@echo "========================================"
	./$(TEST_DATATYPES_EXEC)
	./$(TEST_DIVDIFF_EXEC)
	./$(TEST_HAMILTONIAN_EXEC)
	./$(TEST_STATE_EXEC)
	./$(TEST_UTILS_EXEC)
	./$(TEST_UPDATES_EXEC)
	@echo "\n========================================"
	@echo "      ALL TESTS PASSED SUCCESSFULLY"
	@echo "========================================"

# --- Individual Test Build Rules ---
$(TEST_DATATYPES_EXEC): $(TEST_DATATYPES_OBJS)
	$(CC) $(CFLAGS) -o $(TEST_DATATYPES_EXEC) $(TEST_DATATYPES_OBJS)

$(TEST_DIVDIFF_EXEC): $(TEST_DIVDIFF_OBJS)
	$(CC) $(CFLAGS) -o $(TEST_DIVDIFF_EXEC) $(TEST_DIVDIFF_OBJS) -lm

$(TEST_HAMILTONIAN_EXEC): $(TEST_HAMILTONIAN_OBJS)
	$(CC) $(CFLAGS) -o $(TEST_HAMILTONIAN_EXEC) $(TEST_HAMILTONIAN_OBJS) -lm

$(TEST_STATE_EXEC): $(TEST_STATE_OBJS)
	$(CC) $(CFLAGS) -o $(TEST_STATE_EXEC) $(TEST_STATE_OBJS) -lm

$(TEST_UTILS_EXEC): $(TEST_UTILS_OBJS)
	$(CC) $(CFLAGS) -o $(TEST_UTILS_EXEC) $(TEST_UTILS_OBJS)

$(TEST_UPDATES_EXEC): $(TEST_UPDATES_OBJS)
	$(CC) $(CFLAGS) -o $(TEST_UPDATES_EXEC) $(TEST_UPDATES_OBJS) -lm


# Generic rule to compile .c files into .o object files
%.o: %.c
	$(CC) $(CFLAGS) -c $< -o $@

# Clean up all build files
clean:
	rm -f $(MAIN_OBJS) $(TEST_DATATYPES_OBJS) $(TEST_DIVDIFF_OBJS) $(TEST_HAMILTONIAN_OBJS) $(TEST_STATE_OBJS) $(TEST_UTILS_OBJS) $(TEST_UPDATES_OBJS)
	rm -f $(MAIN_EXEC) $(TEST_DATATYPES_EXEC) $(TEST_DIVDIFF_EXEC) $(TEST_HAMILTONIAN_EXEC) $(TEST_STATE_EXEC) $(TEST_UTILS_EXEC) $(TEST_UPDATES_EXEC)

.PHONY: all clean test