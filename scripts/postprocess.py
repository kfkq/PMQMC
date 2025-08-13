# File: scripts/postprocess.py
# Purpose: Fully automated post-processing and statistical analysis of PMQMC data.
#
# Workflow:
# 1. Reads raw HDF5 data and simulation parameters.
# 2. Calculates the autocorrelation time of the energy.
# 3. Determines an optimal number of bins based on the autocorrelation time.
# 4. Performs Jackknife resampling on binned data to get robust error bars.
# 5. Writes the final results to 'default_obs.dat' in a clean, fixed-width format.

import numpy as np
import h5py
import sys
from pathlib import Path

# --- Constants ---
HDF5_FILENAME = "raw_data.h5"
RESULTS_FILENAME = "default_obs.dat"
MIN_BINS_FOR_RELIABLE_STATS = 20 # Warn the user if suggested nbins is below this

def calculate_autocorrelation_time(data_series: np.ndarray) -> float:
    """Calculates the integrated autocorrelation time using the FFT method."""
    try:
        series = data_series - np.mean(data_series)
        f = np.fft.fft(series, n=2 * len(series))
        acf = np.fft.ifft(f * np.conj(f))[:len(series)].real
        acf /= acf[0]
        
        tau_int = 0.0
        for val in acf:
            if val > 0.05:
                tau_int += val
            else:
                break
        return max(1.0, 2 * tau_int - 1)
    except Exception:
        return 1.0

def calculate_observables_from_bins(binned_data: dict, beta: float) -> dict:
    """Calculates final observables from a set of bin averages."""
    results = {}
    
    mean_sgn = np.mean(binned_data['sgn'])
    
    # Use np.get() to safely access keys that might not exist
    mean_H_sgn = np.mean(binned_data.get('H_sgn', 0))
    mean_H2_sgn = np.mean(binned_data.get('H2_sgn', 0))
    mean_Z_mag_sgn = np.mean(binned_data.get('Z_mag_sgn', 0))

    # Handle the sign problem with the ratio formula
    results['H'] = mean_H_sgn / mean_sgn if 'H_sgn' in binned_data else 0
    results['H2'] = mean_H2_sgn / mean_sgn if 'H2_sgn' in binned_data else 0
    results['Z_mag'] = mean_Z_mag_sgn / mean_sgn if 'Z_mag_sgn' in binned_data else 0
    results['Cv'] = (results['H2'] - results['H']**2) * beta**2 if 'H' in results and 'H2' in results else 0
    
    return results

def main():
    """Main function to drive the analysis."""
    print("--- PMQMC Python Analyzer ---")

    # 1. Load data and parameters from HDF5 file
    if not Path(HDF5_FILENAME).exists():
        print(f"Error: Data file '{HDF5_FILENAME}' not found in the current directory.")
        sys.exit(1)

    with h5py.File(HDF5_FILENAME, 'r') as f:
        dset = f['/measurements']
        data = dset[:]
        params = dict(dset.attrs)
    
    print(f"Successfully loaded {len(data)} measurements.")
    skip = int(params.get('SKIP_MEASUREMENTS', 0))
    
    # 2. Calculate autocorrelation time and determine nbins
    if 'H' not in data.dtype.names:
        print("Warning: 'H' data not found. Cannot calculate autocorrelation time. Defaulting to 100 bins.")
        nbins = 100
    else:
        print("Calculating autocorrelation time for energy <H>...")
        tau = calculate_autocorrelation_time(data['H'][skip:])
        print(f"Integrated autocorrelation time (tau) ≈ {tau:.2f} measurements")

        min_bin_length = int(np.ceil(5 * tau))
        num_points_to_analyze = len(data) - skip
        suggested_nbins = num_points_to_analyze // min_bin_length
        
        if suggested_nbins < MIN_BINS_FOR_RELIABLE_STATS:
            print(f"\nWARNING: The suggested number of bins ({suggested_nbins}) is less than {MIN_BINS_FOR_RELIABLE_STATS}.")
            print("         This indicates the simulation run may be too short compared to its")
            print("         autocorrelation time. The resulting error bars may be unreliable.")
            print("         Consider running the simulation for more steps.\n")
        
        nbins = suggested_nbins

    print(f"Using {nbins} bins for Jackknife analysis.")

    # 3. Prepare data for analysis
    data = data[skip:]
    bin_length = len(data) // nbins
    if bin_length == 0:
        print(f"Error: Not enough data for {nbins} bins after skipping. Exiting.")
        sys.exit(1)
    data = data[:nbins * bin_length]

    # 4. Create primary bins
    binned_data = {}
    measured_obs = []
    for name in data.dtype.names:
        if name in ['sgn', 'q']:
            binned_data[name] = data[name].reshape((nbins, bin_length)).mean(axis=1)
        else:
            # Observables are weighted by the sign
            binned_data[f'{name}_sgn'] = (data[name] * data['sgn']).reshape((nbins, bin_length)).mean(axis=1)
            measured_obs.append(name)

    # 5. Perform Jackknife Resampling
    jackknife_estimates = {obs: [] for obs in measured_obs + ['Cv']}
    
    for i in range(nbins):
        jk_binned_data = {key: np.delete(val, i) for key, val in binned_data.items()}
        jk_results = calculate_observables_from_bins(jk_binned_data, params['BETA'])
        
        for key in jackknife_estimates:
            # Use .get() to handle cases where an observable (like Z_mag) wasn't measured
            # and thus won't be in jk_results.
            if jk_results.get(key) is not None:
                jackknife_estimates[key].append(jk_results.get(key))

    # 6. Calculate final Jackknife mean and standard error
    final_stats = {}
    for key, estimates in jackknife_estimates.items():
        if estimates: # Only calculate if the observable was measured
            n_jk = len(estimates)
            jk_mean = np.mean(estimates)
            jk_var = ((n_jk - 1) / n_jk) * np.sum((estimates - jk_mean)**2)
            final_stats[f'mean_{key}'] = jk_mean
            final_stats[f'stdev_{key}'] = np.sqrt(jk_var)

    # Add non-Jackknife stats
    final_stats['mean_sgn'] = np.mean(binned_data['sgn'])
    final_stats['stdev_sgn'] = np.std(binned_data['sgn'], ddof=1) / np.sqrt(nbins)
    final_stats['mean_q'] = np.mean(data['q'])
    final_stats['max_q'] = np.max(data['q'])

    # 7. Write results to file in the specified format
    with open(RESULTS_FILENAME, 'w') as f:
        f.write(f"mean(q) = {final_stats['mean_q']:.4f}\n")
        f.write(f"max(q) = {final_stats['max_q']}\n")
        
        # Build header and data strings dynamically
        header = ""
        data_line = ""
        
        # Always include sign
        header += f"{'<sgn>':>18s}{'stdev(sgn)':>18s}"
        data_line += f"{final_stats['mean_sgn']:>18.8f}{final_stats['stdev_sgn']:>18.8f}"
        
        # Add other observables if they were measured
        if 'H' in measured_obs:
            header += f"{'<H>':>18s}{'stdev(H)':>18s}"
            data_line += f"{final_stats.get('mean_H', 0):>18.8f}{final_stats.get('stdev_H', 0):>18.8f}"
        if 'H2' in measured_obs:
            header += f"{'<H^2>':>18s}{'stdev(H^2)':>18s}"
            data_line += f"{final_stats.get('mean_H2', 0):>18.8f}{final_stats.get('stdev_H2', 0):>18.8f}"
        if 'Z_mag' in measured_obs:
            header += f"{'<Z_mag>':>18s}{'stdev(Z_mag)':>18s}"
            data_line += f"{final_stats.get('mean_Z_mag', 0):>18.8f}{final_stats.get('stdev_Z_mag', 0):>18.8f}"

        f.write(header + "\n")
        f.write(data_line + "\n")

    print(f"\nAnalysis complete. Results written to '{RESULTS_FILENAME}'")

if __name__ == "__main__":
    main()