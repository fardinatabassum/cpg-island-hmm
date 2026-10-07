import matplotlib.pyplot as plt
import numpy as np
import math
import os
from src.viterbi_hmm import viterbi_decode
from src.data_reader import read_fasta


def run_visualization(fasta_path, start_bp, end_bp):
    """
    Optimized visualization that slices the sequence before decoding
    """
    # Setup parameters
    states = [0, 1]
    start_p = [math.log(0.5), math.log(0.5)]

    # Adjusted Transition Probabilities
    trans_p = {
        0: {0: math.log(0.9999), 1: math.log(0.0001)},
        1: {0: math.log(0.001), 1: math.log(0.999)}
    }

    # Emission calibrated for human Chr22 GC density
    emiss_p = {
        0: {'A': math.log(0.25), 'C': math.log(0.25), 'G': math.log(0.25), 'T': math.log(0.25)},
        1: {'A': math.log(0.10), 'C': math.log(0.40), 'G': math.log(0.40), 'T': math.log(0.10)}
    }

    # Local slicing
    print(f"--- Loading {fasta_path} ---")
    full_sequence = read_fasta(fasta_path)

    # Add a 500bp buffer
    buffer = 500
    local_start = max(0, start_bp - buffer)
    local_end = min(len(full_sequence), end_bp + buffer)

    sequence = full_sequence[local_start:local_end]
    print(f"Decoding local region: {start_bp} to {end_bp} (Total: {len(sequence)} bp)")

    # Run decoder
    path = viterbi_decode(sequence, states, start_p, trans_p, emiss_p)

    # Slice the results back to the requested window
    sub_seq = sequence[buffer:-buffer] if local_start > 0 else sequence[:end_bp - start_bp]
    sub_path = path[buffer:-buffer] if local_start > 0 else path[:end_bp - start_bp]
    x = np.arange(start_bp, start_bp + len(sub_seq))

    # Calculated GC content
    gc_content = []
    window_size = 100
    for i in range(len(sub_seq)):
        # Handle edges for the sliding window
        w_start = max(0, i - window_size // 2)
        w_end = min(len(sub_seq), i + window_size // 2)
        window = sub_seq[w_start:w_end]
        gc = (window.count('G') + window.count('C')) / len(window)
        gc_content.append(gc)

    # Visualization
    fig, ax1 = plt.subplots(figsize=(12, 6))

    # Plot HMM State (Blue)
    ax1.set_xlabel('Genomic Position on Chr22 (bp)', fontsize=12)
    ax1.set_ylabel('HMM State (0=BG, 1=Island)', color='tab:blue', fontsize=12)
    ax1.fill_between(x, 0, sub_path, color='tab:blue', alpha=0.2, label='Predicted Island')
    ax1.step(x, sub_path, where='post', color='tab:blue', linewidth=1.5)
    ax1.tick_params(axis='y', labelcolor='tab:blue')
    ax1.set_ylim(-0.1, 1.1)

    # Plot GC Content (Red)
    ax2 = ax1.twinx()
    ax2.set_ylabel('GC Content (100bp window)', color='tab:red', fontsize=12)
    ax2.plot(x, gc_content, color='tab:red', alpha=0.4, label='GC %')
    ax2.tick_params(axis='y', labelcolor='tab:red')
    ax2.axhline(y=0.5, color='black', linestyle='--', alpha=0.3)  # 50% Threshold

    plt.title(f"Viterbi Path vs. GC Content: Chr22 ({start_bp} - {end_bp})", fontsize=14, fontweight='bold')
    plt.tight_layout()

    # Create results folder if it does not exist
    os.makedirs("results", exist_ok=True)
    plt.savefig("results/final_cpg_visualization.png", dpi=300)
    print(f"Successfully saved plot to results/final_cpg_visualization.png")
    plt.show()


if __name__ == "__main__":
    # Window selected from coordinate: chr22  39,349,653 - 39,350,483
    run_visualization("data/chr22.fa", 39348500, 39351500)
  
