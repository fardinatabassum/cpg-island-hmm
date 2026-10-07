import argparse


def load_bed(filename):
    """Parses a BED file into a list of tuples."""
    intervals = []
    with open(filename, 'r') as f:
        for line in f:
            if line.startswith("#") or not line.strip(): continue
            parts = line.split()
            # BED format: Chrom, Start, End
            intervals.append((int(parts[1]), int(parts[2])))
    return intervals


def calculate_metrics(predicted, reference):
    """
    Calculates Sensitivity and Precision based on interval overlaps.
    """
    ref_found_flags = [False] * len(reference)
    true_positives = 0

    # Overlap check
    for p_start, p_end in predicted:
        has_overlap = False
        for i, (r_start, r_end) in enumerate(reference):
            # Check for any overlap between intervals
            if max(p_start, r_start) < min(p_end, r_end):
                ref_found_flags[i] = True
                has_overlap = True
        # If this prediction hit any reference it is a true positive for Precision
        if has_overlap:
            true_positives += 1


    sensitivity = (sum(ref_found_flags) / len(reference)) * 100 if reference else 0
    precision = (true_positives / len(predicted)) * 100 if predicted else 0

    return sensitivity, precision, sum(ref_found_flags)


def main():
    parser = argparse.ArgumentParser(description="HMM Performance Evaluation")
    parser.add_argument("-p", "--predicted", required=True, help="Your predictions.bed")
    parser.add_argument("-r", "--reference", required=True, help="UCSC reference.bed")
    args = parser.parse_args()

    print(f"Evaluating {args.predicted} against {args.reference}")

    pred_intervals = load_bed(args.predicted)
    ref_intervals = load_bed(args.reference)

    sens, prec, tp = calculate_metrics(pred_intervals, ref_intervals)

    print(f"Total Reference Islands: {len(ref_intervals)}")
    print(f"Total Predicted Islands: {len(pred_intervals)}")
    print(f"True Positives (Overlaps): {tp}")
    print("-" * 30)
    print(f"Sensitivity (Recall): {sens:.2f}%")
    print(f"Precision:            {prec:.2f}%")


if __name__ == "__main__":
    main()
