import numpy as np
import math
import argparse
from src.data_reader import read_fasta, get_header
from src.utils import save_as_bed


def viterbi_decode(sequence, states, start_p, trans_p, emiss_p, global_offset=0):
    """
    Performs Viterbi decoding to identify CpG islands in a genomic sequence.

    This implementation uses log-space transformations to ensure numerical stability
    over long genomic sequences (such as Human Chromosome 22) and includes a
    neutral coasting logic to handle assembly gaps ('N' characters).

    Args:
        sequence (str): The DNA sequence string (A, C, G, T, N).
        states (list): List of hidden states [0, 1] where 0=Background, 1=Island.
        start_p (list): Log-probabilities for the initial state distribution.
        trans_p (dict): Nested dictionary of log-transition probabilities.
        emiss_p (dict): Nested dictionary of log-emission probabilities.
        global_offset (int): Offset to align local indices with genomic coordinates.

    Returns:
        list: A list of tuples (start, end) representing detected CpG islands.
    """
    L = len(sequence)
    S = len(states)

    # Initialization
    # V[state, index] stores the maximum log-probability of reaching a state
    # B[state, index] stores the backpointer to the most likely previous state
    V = np.full((S, L), -np.inf)
    B = np.zeros((S, L), dtype=int)

    for s in range(S):
        # Initial score = Start Probability + Emission Probability of the first base
        V[s, 0] = start_p[s] + emiss_p[s].get(sequence[0], math.log(0.25))

    # Main Recursive Loop (The Trellis)
    for i in range(1, L):
        current_obs = sequence[i]

        # Assembly Gap Handling ('N')
        # If an 'N' is encountered, we maintain the previous log-score and
        # coast through the gap to preserve genomic coordinate alignment.
        if current_obs not in ['A', 'C', 'G', 'T']:
            for s in range(S):
                V[s, i] = V[s, i - 1] + math.log(0.25)
                B[s, i] = s
            continue

        for s in range(S):
            # Calculate (Previous Score + Transition to Current State) for all possible origins
            scores = [V[prev_s, i - 1] + trans_p[prev_s][s] for prev_s in range(S)]

            # Store the maximum path score and the index of the state that produced it
            V[s, i] = max(scores) + emiss_p[s].get(current_obs, math.log(0.25))
            B[s, i] = np.argmax(scores)

    # Path Reconstruction (Backtracking)
    # Identify the best ending state and follow the backpointers (B) from right to left
    best_final_state = np.argmax(V[:, L - 1])
    path = np.zeros(L, dtype=int)
    path[L - 1] = best_final_state

    for i in range(L - 2, -1, -1):
        # Looks at the backpointer for the state we just found at index i+1
        path[i] = B[path[i + 1], i + 1]

    return path

def main():
    """
    CLI execution block.
    Allows users to specify input FASTA and output BED files via the terminal.
    """
    parser = argparse.ArgumentParser(description="HMM-based CpG Island Detector for Genomic Sequences")
    parser.add_argument("-i", "--input", required=True, help="Path to the input FASTA file")
    parser.add_argument("-o", "--output", default="predictions.bed", help="Path to save the output BED file")
    args = parser.parse_args()

    # Model calibration States 0 (Background) and 1 (CpG Island)
    states = [0, 1]
    start_p = [math.log(0.5), math.log(0.5)]

    # Transition matrix adjusted to favor state persistence
    trans_p = {
        0: {0: math.log(0.9999), 1: math.log(0.0001)},
        1: {0: math.log(0.001), 1: math.log(0.999)}
    }

    # Emission matrix weighted heavily for G/C content in the Island state
    emiss_p = {
        0: {'A': math.log(0.25), 'C': math.log(0.25), 'G': math.log(0.25), 'T': math.log(0.25)},
        1: {'A': math.log(0.10), 'C': math.log(0.40), 'G': math.log(0.40), 'T': math.log(0.10)}
    }

    try:
        print(f"--- Loading data from {args.input} ---")
        sequence = read_fasta(args.input)
        chrom_name = get_header(args.input)
        print(f"Successfully loaded {len(sequence)} base pairs from {chrom_name}.")

        print("Running Viterbi decoding")
        # results = viterbi_decode(sequence, states, start_p, trans_p, emiss_p)
        path = viterbi_decode(sequence, states, start_p, trans_p, emiss_p)

        # Extract islands and apply the 200bp biological filter
        islands, current_start = [], None
        for idx, state in enumerate(path):
            if state == 1 and current_start is None:
                current_start = idx
            elif state == 0 and current_start is not None:
                if (idx - current_start) >= 200:
                    islands.append((current_start, idx))
                current_start = None

        # Capture islands that reach the end of the sequence
        if current_start is not None and (len(path) - current_start) >= 200:
            islands.append((current_start, len(path)))

        # Output results
        if islands:
            print(f"Successfully detected {len(islands)} CpG islands.")
            # Save to the root directory
            save_as_bed(islands, args.output, chrom=chrom_name)
        else:
            print("No islands detected. Check if sequence is >200bp")

    except FileNotFoundError:
        print(f"Error: {args.input} not found. Verify the file path relative to the root directory.")

    except Exception as e:
        print(f"An unexpected error occurred: {e}")


if __name__ == "__main__":
    main()
