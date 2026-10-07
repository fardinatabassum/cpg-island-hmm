def read_fasta(file_path):
    """Parses a FASTA file, strips headers, and normalizes sequence to uppercase"""
    sequence = []
    with open(file_path, 'r') as f:
        for line in f:
            if not line.startswith(">"):
                sequence.append(line.strip().upper())
    return "".join(sequence)

def get_header(file_path):
    """Retrieves the chromosome header for the BED file."""
    with open(file_path, 'r') as f:
        for line in f:
            if line.startswith(">"):
                return line.strip()[1:].split()[0]
    return "unknown_chrom"
