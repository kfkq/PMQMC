# File: scripts/postprocess.py
# Purpose: Fully automated post-processing and statistical analysis of PMQMC data.
#
# Workflow:
# 1. Reads raw HDF5 data and simulation parameters.
# 2. Automatically detects the thermalization point using a robust converging window method.
# 3. Calculates the autocorrelation time on the equilibrated data.
# 4. Determines an optimal number of bins.
# 5. Performs Jackknife resampling to get robust error bars.
# 6. Writes the final results to 'default_obs.dat'.

import numpy as np
import h5py
import sys
from pathlib import Path
import matplotlib.pyplot as plt

# --- Constants ---
HDF5_FILENAME = "raw_data.h5"
RESULTS_FILENAME = "default_obs.dat"
PLOT_FILENAME = "equilibration_plot.png"
MIN_BINS_FOR_RELIABLE_STATS = 20

def find_equilibration_point_robust(data_series: np.ndarray) -> int:
    """
    Estimates the equilibration point by finding the first window of data from the
    start that is statistically consistent with the last half of the data.
    """
    n_total = len(data_series)
    if n_total < 200: # Need enough data for a reliable baseline
        print("Warning: Not enough data to reliably detect equilibration point. Using 0.")
        return 0

    # 1. Establish a robust baseline from the last 50% of the data
    n_ref = n_total // 2
    reference_data = data_series[-n_ref:]
    ref_mean = np.mean(reference_data)
    ref_std = np.std(reference_data)

    if ref_std == 0: # Data is flat, likely already equilibrated
        return 0

    # 2. Scan from the beginning with a sliding window
    # Use a window size that is 1% of the total, but at least 100 points.
    window_size = max(100, n_total // 100)
    
    for i in range(0, n_total - window_size, window_size):
        window_data = data_series[i : i + window_size]
        window_mean = np.mean(window_data)
        window_std = np.std(window_data)

        # 3. Statistical test: Is the window mean consistent with the reference mean?
        # We check if the difference in means is significant compared to the
        # combined standard error of the two means.
        std_error_diff = np.sqrt(ref_std**2 / n_ref + window_std**2 / window_size)
        
        # If the window mean is within ~3 standard errors of the reference mean,
        # we declare this point as the start of equilibrium.
        if np.abs(window_mean - ref_mean) < 3 * std_error_diff:
            return i

    # If no convergence was found, something is unusual.
    # Warn the user and return a safe default (e.g., discard first 10%).
    print("Warning: Could not automatically detect a clear equilibration point.")
    print("         The data may not have converged. Using first 10% as a fallback.")
    return n_total // 10

def create_equilibration_plot(data_series: np.ndarray, equilibration_point: int, output_filename: str):
    """
    Generates and saves a plot showing the full time series and the
    detected equilibration point.
    """
    print(f"\nGenerating equilibration plot and saving to '{output_filename}'...")
    
    fig, ax = plt.subplots(figsize=(12, 6))
    
    x_axis = np.arange(len(data_series))
    
    # Plot the discarded (unequilibrated) part in a different color
    ax.plot(x_axis[:equilibration_point], data_series[:equilibration_point], 
            color='salmon', label='Discarded Data (Thermalization)')
    
    # Plot the equilibrated part in a more prominent color
    ax.plot(x_axis[equilibration_point:], data_series[equilibration_point:], 
            color='steelblue', label='Equilibrated Data')

    # Add a vertical line to mark the cutoff
    ax.axvline(x=equilibration_point, color='r', linestyle='--', linewidth=2, 
               label=f'Equilibration Point = {equilibration_point}')

    ax.set_title('Energy Time Series and Detected Equilibration Point')
    ax.set_xlabel('Measurement Step')
    ax.set_ylabel('Energy <H>')
    ax.legend()
    ax.grid(True, linestyle=':', alpha=0.6)
    
    try:
        plt.savefig(output_filename, dpi=150, bbox_inches='tight')
    except Exception as e:
        print(f"Warning: Could not save plot. Error: {e}")
    plt.close(fig)

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
    
    mean_H_sgn = np.mean(binned_data.get('H_sgn', 0))
    mean_H2_sgn = np.mean(binned_data.get('H2_sgn', 0))
    mean_Z_mag_sgn = np.mean(binned_data.get('Z_mag_sgn', 0))

    results['H'] = mean_H_sgn / mean_sgn if 'H_sgn' in binned_data and mean_sgn != 0 else 0
    results['H2'] = mean_H2_sgn / mean_sgn if 'H2_sgn' in binned_data and mean_sgn != 0 else 0
    results['Z_mag'] = mean_Z_mag_sgn / mean_sgn if 'Z_mag_sgn' in binned_data and mean_sgn != 0 else 0
    results['Cv'] = (results['H2'] - results['H']**2) * beta**2 if 'H' in results and 'H2' in results else 0
    
    return results

def main():
    """Main function to drive the analysis."""
    print("--- PMQMC Python Analyzer ---")

    # 1. Load data
    if not Path(HDF5_FILENAME).exists():
        print(f"Error: Data file '{HDF5_FILENAME}' not found.")
        sys.exit(1)

    with h5py.File(HDF5_FILENAME, 'r') as f:
        dset = f['/measurements']
        data = dset[:]
        params = dict(dset.attrs)
    
    print(f"Successfully loaded {len(data)} measurements.")
    
    if 'H' not in data.dtype.names:
        print("Error: 'H' data not found in HDF5 file. Cannot perform analysis.")
        sys.exit(1)

    # 2. Automatically detect equilibration point
    print("\nDetecting equilibration point using robust converging window method...")
    skip = find_equilibration_point_robust(data['H'])
    print(f"Data appears to be equilibrated after approximately {skip} measurements.")
    
    equilibrated_data = data[skip:]
    num_points_to_analyze = len(equilibrated_data)

    # 3. Generate the plot
    create_equilibration_plot(data['H'], skip, PLOT_FILENAME)

    equilibrated_data = data[skip:]
    num_points_to_analyze = len(equilibrated_data)

    # 3. Calculate autocorrelation time on equilibrated data
    print("\nCalculating autocorrelation time on equilibrated data...")
    tau = calculate_autocorrelation_time(equilibrated_data['H'])
    print(f"Integrated autocorrelation time (tau) ≈ {tau:.2f} measurements")

    min_bin_length = int(np.ceil(3 * tau))
    suggested_nbins = num_points_to_analyze // min_bin_length if min_bin_length > 0 else 200
    
    if suggested_nbins < MIN_BINS_FOR_RELIABLE_STATS:
        print(f"\nWARNING: The suggested number of bins ({suggested_nbins}) is less than {MIN_BINS_FOR_RELIABLE_STATS}.")
        print("         This indicates the simulation run may be too short for its")
        print("         autocorrelation time. The resulting error bars may be unreliable.")
    
    nbins = max(suggested_nbins, 1)
    print(f"Using {nbins} bins for final Jackknife analysis.")

    # 4. Prepare data for analysis
    bin_length = num_points_to_analyze // nbins
    if bin_length == 0:
        print(f"Error: Not enough data for {nbins} bins after skipping. Exiting.")
        sys.exit(1)
    final_data = equilibrated_data[:nbins * bin_length]

    # 5. Create primary bins
    binned_data = {}
    measured_obs = []
    for name in final_data.dtype.names:
        if name in ['sgn', 'q']:
            binned_data[name] = final_data[name].reshape((nbins, bin_length)).mean(axis=1)
        else:
            binned_data[f'{name}_sgn'] = (final_data[name] * final_data['sgn']).reshape((nbins, bin_length)).mean(axis=1)
            measured_obs.append(name)

    # 6. Perform Jackknife Resampling
    jackknife_estimates = {obs: [] for obs in measured_obs + ['Cv'] if obs in final_data.dtype.names or obs == 'Cv'}
    
    for i in range(nbins):
        jk_binned_data = {key: np.delete(val, i) for key, val in binned_data.items()}
        jk_results = calculate_observables_from_bins(jk_binned_data, params.get('BETA', 1.0))
        
        for key in jackknife_estimates:
            if jk_results.get(key) is not None:
                jackknife_estimates[key].append(jk_results.get(key))

    # 7. Calculate final Jackknife mean and standard error
    final_stats = {}
    for key, estimates in jackknife_estimates.items():
        if estimates:
            n_jk = len(estimates)
            jk_mean = np.mean(estimates)
            jk_var = ((n_jk - 1) / n_jk) * np.sum((estimates - jk_mean)**2)
            final_stats[f'mean_{key}'] = jk_mean
            final_stats[f'stdev_{key}'] = np.sqrt(jk_var)

    # Add non-Jackknife stats
    final_stats['mean_sgn'] = np.mean(binned_data['sgn'])
    final_stats['stdev_sgn'] = np.std(binned_data['sgn'], ddof=1) / np.sqrt(nbins)
    final_stats['mean_q'] = np.mean(final_data['q'])
    final_stats['max_q'] = np.max(final_data['q'])

    # 8. Write results to file
    with open(RESULTS_FILENAME, 'w') as f:
        f.write(f"mean(q) = {final_stats['mean_q']:.4f}\n")
        f.write(f"max(q) = {final_stats['max_q']}\n")
        
        header = ""
        data_line = ""
        
        header += f"{'<sgn>':>18s}{'stdev(sgn)':>18s}"
        data_line += f"{final_stats['mean_sgn']:>18.8f}{final_stats['stdev_sgn']:>18.8f}"
        
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