# CpG Island Detector using HMM

## **Project Overview**
### Project summary:
* **Research Question:** To what extent can a first-order two-state Hidden Markov Model (HMM) accurately identify CpG islands on Human Chromosome 22 compared to the UCSC Genome Browser’s curated reference benchmarks?
* **Algorithm:** Viterbi Algorithm implemented utilizing log-transformed probabilities to prevent numerical underflow, integrated with a neutral-emission handling strategy that preserves the linear genomic scale during non-informative assembly gaps. The output was also filtered using a 200bp length constraint.
  * **Input**: 
    * `FASTA` Human Chromosome 22 (GRCh38) 
    * `FASTA` Synthetic prototype sequence used for logic validation
    * `.bed` file from the UCSC Genome Browser (cpgIslandExt track) to evaluate model performance
  * **Outputs**:  
    * `.bed` 3-column BED file containing the specific start and end coordinates of every predicted CpG island that passed the 200bp biological length filter
    * `.txt` A statistical summary report containing the model's Sensitivity and Precision.
    * `.png` A genomic track plot that overlays the HMM predicted Viterbi path with local GC-content percentages, providing visual confirmation of the model’s state transitions across the chromosome.

## **Installation and Setup**

### Requirements
* Python Version: Python 3.8 or higher is required.

### Clone the repository: 
```bash
git clone https://github.com/fardinatabassum/cpg-island-hmm.git
cd cpg-island-hmm
```

### Install dependencies:
```bash 
pip install -r requirements.txt
```

### Download the Chromosome 22 FASTA
```bash
wget https://hgdownload.soe.ucsc.edu/goldenPath/hg38/chromosomes/chr22.fa.gz
gunzip chr22.fa.gz
mkdir -p data && mv chr22.fa ./data/
```

* **System-Level Requirements:**
  * Memory: The HMM uses log-space arithmetic and NumPy vectors to ensure the Chromosome 22 sequence can be decoded within standard 8GB/16GB RAM limits.
  * No External Compilers: The implementation is written in native Python.

## **Quick Start**

### **Data Preparation** 
* The `data/` folder contains a synthetic `prototype_genome.fa`. To run the full analysis, download hg38 Chromosome 22 and place it in the directory as `data/chr22.fa`.

### **Execute Analysis using the CLI**

#### Prototype Sequence
```bash
# Verify the HMM logic on the prototype sequence
python3 -m src.viterbi_hmm -i data/prototype_genome.fa -o results/prototype_islands.bed

# Evaluate the prototype results
python3 -m src.evaluation -p results/prototype_islands.bed -r data/ucsc_reference.bed
```

**Expected Output:**
* The script will output a `prototype_islands.bed` file located in `results/prototype_islands.bed`.
* **Success Indicator**: The console should report a sensitivity of ~75.85% when compared to the reference.
* **File Sample**: `chr_prototype 118 420`

#### Real-Data (Chromosome 22)
```bash
# Run the Viterbi decoding
python3 -m src.viterbi_hmm -i data/chr22.fa -o results/chr22_islands.bed --min-len 200

# Performance Evaluation
python3 -m src.evaluation -p results/chr22_islands.bed -r data/ucsc_reference.bed
```

## **Usage and Options**

### Main Scripts
* `viterbi_hmm.py`: This is the primary driver for the analysis. It parses the FASTA sequence, performs log-space decoding, and writes the results.
  * `-i`, `--input`: Path to the input genomic FASTA file (`data/chr22.fa`)
  * `-o`, `--output`: Path to save the resulting BED file (defaults to `results/out.bed`)
  * `--min-len`: The minimum length threshold (in bp) for a predicted island to be reported (defaults to `200`)
* `evaluation.py`: This utility compares the predictions against the UCSC benchmark to calculate sensitivity and precision.
* `visualization.py`: Generates the genomic track plot overlaying GC content and predicted states.

### Configuration & Parameters
Adjusting these constants will significantly alter the model's sensitivity:
* **Transition Probabilities (`trans_p`):** These control the expected average length and frequency of islands. The current calibration favors long, contiguous segments ($P = 0.999$).
* **Emission Probabilities (`emiss_p`):** These define the base composition of a CpG Island state (set at 40% G and 40% C) versus the Background state (25% for all bases).

### Advanced Usage Examples
* To validate the model on a smaller sequence and generate the GC-content track plot:
```bash
python3 -m src.viterbi_hmm -i data/prototype_genome.fa -o results/viz_test.bed --min-len 200 --visualize
```

## **Limitations and Assumptions**

### Key Assumptions
* The algorithm assumes the current genomic state is only dependent on the immediate preceding state due to it being a First-Order HMM.
* The genomic sequence length requires log-transformation of probabilities to prevent numerical underflow.
* The parser assumes FASTA headers follow standard NCBI or UCSC formatting and that input sequences consist of {A, C, G, T, N}.

### Known Limitations
* While the model achieves high sensitivity (~75.85%), precision (~57.75%) is limited by the model's tendency to predict short, high-GC regions that may lack regulatory function.
* Intermediate GC-content regions can cause state toggling ("flicker"), occasionally splitting a single biological island into smaller segments.
* Storing the full backpointer matrix for 51 million bases may cause latency on systems with less than 8GB of available RAM.
* The neutral-emission strategy preserves coordinate alignment during assembly gaps but cannot infer the biological state of unsequenced regions.
