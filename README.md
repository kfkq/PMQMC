# PM-QMC

### Introduction

A high-performance C implementation of the Permutation Matrix Representation Quantum Monte Carlo (PMR-QMC) algorithm for simulating arbitrary spin-1/2 systems.

This project uses a two-stage workflow:
1.  `pmqmc`: A fast simulator that dumps raw measurement data to a compressed HDF5 file (`raw_data.h5`).
2.  `statistics`: A lightweight post-processor that reads the HDF5 file, performs binning analysis, and writes the final physical observables to `results.dat`.

### Requirements

*   **C Compiler:** GCC or Clang.
*   **GNU Make:** For building the executables.
*   **Python 3 & NumPy:** For generating the Hamiltonian input file.
    ```bash
    pip install numpy
    ```
*   **HDF5 C Library:** Required for compiling and running.
    *   On Debian/Ubuntu: `sudo apt-get install libhdf5-dev`
    *   On macOS (Homebrew): `brew install hdf5`

### Usage

The workflow consists of three steps:

**1. Generate Input File**

Edit a Python script (e.g., `examples/1d_chain.py`) to define your model and simulation parameters. Then, run it to generate the input file:
```bash
python3 examples/1d_chain.py```
This creates `hamiltonian.in`, which includes a `DEFAULT_MEASUREMENTS_BEGIN` block to control which observables are calculated.

**2. Compile & Run Simulation**

First, edit the `Makefile` to set the `HDF5_INCLUDE_DIR` and `HDF5_LIB_DIR` paths for your system. Then, compile:
```bash
make
```
Run the simulation. It will automatically find `hamiltonian.in`:
```bash
./pmqmc
```
This runs the simulation and generates a compressed `raw_data.h5` file.

**3. Analyze Results**

Run the analyzer. It automatically finds `hamiltonian.in` and `raw_data.h5`:
```bash
./statistics
```
This reads the raw data, performs the analysis, and writes the final results to `results.dat`.

### HDF5 Data Structure

The output data is stored in `raw_data.h5`. You can inspect this file with tools like `h5dump` or load it easily in Python (`h5py`), MATLAB, etc. The internal structure is as follows:

*   **/measurements** `(Dataset)`: A 1D array of compound data structures containing the time-series data. Each row has the following fields:
    *   `sgn` (double): The sign of the configuration weight (1.0 or -1.0).
    *   `H` (double): The instantaneous value of the energy, `<H>`.
    *   `H2` (double): The instantaneous value of the energy squared, `<H^2>`.
    *   `Z_mag` (double): The instantaneous value of the Z-magnetization.
    *   `q` (int): The length of the operator sequence.

*   **/measurements** `(Attributes)`: The dataset contains metadata attributes that store the key simulation parameters used for the run:
    *   `N` (int): Number of qubits.
    *   `BETA` (double): Inverse temperature.
    *   `STEPS` (long long): Total number of Monte Carlo updates.
    *   `SKIP_MEASUREMENTS` (long long): Number of initial measurements to discard during analysis.
    *   `NBINS` (int): Number of bins to use for error analysis.
    *   ...and other parameters from the input file.