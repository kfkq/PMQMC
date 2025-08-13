# PMQMC

### Introduction

Permutation Matrix Representation Quantum Monte Carlo (PMR-QMC) for spin-1/2 systems.

### Requirements

*   A C compiler (GCC or Clang)
*   GNU Make
*   Python 3 and NumPy

### Usage

The workflow is a three-step process: generate input, simulate, and analyze.

1.  **Generate Input File**
    ```bash
    # (Edit the python script to define your model and parameters)
    python3 examples/1d_chain.py
    ```
    This creates `hamiltonian.in`.

2.  **Compile & Run Simulation**
    ```bash
    make
    ./pmqmc
    ```
    This reads `hamiltonian.in` and generates `raw.dat`.

3.  **Analyze Results**
    ```bash
    ./statistics_analyzer
    ```
    This reads `raw.dat` and writes the final statistics to `results.dat`.

To remove all generated files, run `make clean`.