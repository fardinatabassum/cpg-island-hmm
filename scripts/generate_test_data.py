import random


def generate_zone(length, gc_content):
    """Generates a DNA string with a specific GC percentage."""
    bases = ['G', 'C', 'A', 'T']
    weights = [gc_content / 2, gc_content / 2, (1 - gc_content) / 2, (1 - gc_content) / 2]
    return "".join(random.choices(bases, weights=weights, k=length))


def create_messy_fasta(filename):
    # Background with low GC
    zone1 = generate_zone(300, 0.40)

    # Real Island with high GC with an 'N' gap
    island_part1 = generate_zone(200, 0.70)
    gap = "N" * 50
    island_part2 = generate_zone(200, 0.70)

    # High GC but shorter than filter 
    tease = generate_zone(100, 0.80)
    background2 = generate_zone(200, 0.40)

    # Terminal Island with High GC at the end
    terminal = generate_zone(250, 0.75)

    full_sequence = zone1 + island_part1 + gap + island_part2 + tease + background2 + terminal

    with open(filename, "w") as f:
        f.write(">Synthetic_Chr22_Test\n")
        # Write 60 characters per line (Standard FASTA format)
        for i in range(0, len(full_sequence), 60):
            f.write(full_sequence[i:i + 60] + "\n")

    print(f"Successfully generated {len(full_sequence)}bp test file: {filename}")


if __name__ == "__main__":
    create_messy_fasta("data/prototype_genome.fa")
