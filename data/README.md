# Data Directory
This directory contains the genomic sequences and reference standards required for the Hidden Markov Model (HMM) validation and final analysis.

## 1. Prototype Data
* `prototype_genome.fa`: A 780bp synthetic FASTA sequence used to validate the Viterbi logic.
  * Format: FASTA
  * Length: 780bp
  * Composition: 
      * 0-119bp: Background noise (AT-rich).
      * 120-420bp: High-GC "Island" simulating a promoter region, containing a 50bp 'N' assembly gap.
      * 421-780bp: Terminal background noise to test state-switching and termination.
#### Note: A Python script for automated prototype sequence generater is planned for a future update to enhance reproducibility.

## 2. Documented Prototype Run

### Clone the repository: 
```bash
git clone https://github.com/fardinatabassum/BINF6251-final-project.git
cd BINF6251-final-project
```
### Install dependencies:
```bash 
pip install -r requirements.txt
```
### Run the prototype:
```bash
python src/viterbi_hmm.py
```

### Results
* **Output:** predictions.bed
* **Location and Use:** in the root directory. Import it into the IGV/UCSC Genome Browser.
* **Expected bed file content:**
```
chr_prototype    118    420
```

### Expected Console Output
````
--- Loading data from data/prototype_genome.fa ---
Successfully loaded 780 base pairs.
Running Viterbi decoding...
Success! Found 1 island(s).
Genomic Coordinates: [(118, 420)]
````

## 3. Obtaining Real-World Genomic Data
```bash
# Download the compressed Chromosome 22 FASTA
wget https://hgdownload.soe.ucsc.edu/goldenPath/hg38/chromosomes/chr22.fa.gz

# Decompress the file
gunzip chr22.fa.gz

# Move to the data directory
mv chr22.fa ./data/
```
