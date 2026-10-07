## BINF6251 Final Project

# Probabilistic Identification of CpG Islands in the Human Genome using Hidden Markov Models

## Research Question
The research question for my final project is, can a first-order Hidden Markov Model be used to accurately identify CpG island boundaries in human promoter regions by differentiating between island and background genomic states?

CpG Islands are important epigenetic regulators located near the transcription start sites in most genes in humans. The methylation state in these islands controls gene expression by acting as an on/off switch. Hypermethylation of CpG islands is linked to various cancers and developmental disorders. Having a way to identify these islands can provide us with valuable insights to understand these regions and how they control responses to the environment. This can then help us understand how diseases like cancer progress.

## Algorithm and Algorithm Class
The algorithm class that I will be using is Hidden Markov Models (HMM). The specific algorithm I plan to implement will be the Viterbi algorithm to decode the most likely hidden state path from the observed nucleotide sequence. I chose this algorithm because HMMs can effectively model the relationship between observable sequence data and unobserved biological states. While the nucleotide sequence is visible, the functional identity of a CpG island versus the genomic background is a hidden property that must be decoded. By using the Viterbi algorithm, we can use the Markov Property to account for the specific transition dependencies that window-based counting methods can often overlook. This probabilistic approach allows for a more accurate identification of genomic boundaries.
   
## Data Plan
I will be using the human genomic sequence data, specifically focusing on Chromosome 22. I chose chromosome 22 because it is one of the smaller chromosomes, which allows for the HMM to run faster, making testing easier and computationally less intensive. However, chromosome 22 also has high gene density, which means it has a high density of promoter regions and, consequently, CpG islands. To validate my model, I will also be using the CpG Islands track from the UCSC Table Browser that identifies the known start and end coordinates of CpG islands on this chromosome.
  
### Sources: 
* UCSC Genome Browser
* CpG Islands track from the UCSC Table Browser
* Synthetic FASTA file for initial logic testing

### Data type(s):
- FASTA file for raw genomic sequences of chromosome 22
- BED file genomic coordinates of CpG islands to calculate sensitivity and specificity

Licensing or access considerations: Public domain is being used, so no special licenses or permissions are required for access

### Prototype data:
This data will be a small synthetic DNA sequence. This sequence will be mostly background, containing low GC content with a region in the middle containing a high frequency of C and G nucleotides. This allows me to control and know where the hidden state is, thus allowing me to verify the HMM logic. I can check for the transition and emission probabilities to see if they correctly identify the transition from background to island before I move to a much larger dataset of the actual Chromosome 22. The prototype relates to the realistic dataset as a controlled benchmark. It mimics the statistical properties of human genomic DNA, such as CG suppression in the background and high density in the islands.

## Success Criteria
* The project is successful if the Viterbi implementation correctly decoded the most likely path of hidden states as Island or Background without encountering any errors. The model should also be statistically significant and overlap with the UCSC Gold Standard annotations. Achieving a high sensitivity on Chromosome 22, meaning that the model correctly identifies the majority of known CpG islands despite the background noise.
  
* Expected outputs:
   * A BED file containing the start and end positions of every predicted CpG island on Chromosome 22
   * A visualization showing the state probability for a sample region, demonstrating how the algorithm decided a region was an island

* Verifying results: 
   * Test the HMM against my prototype_genome.fa before using real genomic data
   * Compare the coordinates predicted by my hmm.py against the UCSC CpG Island track
   * Take a known promoter region on Chromosome 22 and verify that the HMM predicts an island where the Genome Browser shows high GC density      and existing annotations.

## Pitfall Scan

### Data-related issues:

* Large sections of the genomic data in Chromosome 22 are not yet fully sequenced, so they are represented by the character "N" in FASTA files. If the HMM tries to produce this character from an Island or Background state without a predefined probability, the script will crash or produce biased results.
  
* We can mitigate this by including a data cleaning step in our file parsing script to look for N characters and split the sequences at these gaps. That way, only contiguous sequences will be processed, which can then be merged.

### Algorithmic issues:

* HMM requires multiplying long chains of probabilities (transition x emission) as it moves across a sequence. Chromosome 22 is approximately 51 million base pairs. These probabilities will become very small, dropping below the floating-point precision of Python and resulting it to be rounded to zero.

* We can mitigate this by implementing the algorithm with a log-space transformation. Logarithms turn multiplication into addition, thus keeping the numerical stability.

### Evaluation issues: 

* CpG islands are relatively rare in the human genome. If the model incorrectly predicts Background state for the entire chromosome, a standard accuracy metric would still report a high success, which is not biologically correct.

* We can mitigate this by using sensitivity and specificity to evaluate performance. I will also generate a Confusion Matrix to specifically track how many true islands were missed (false negatives) versus how many Islands were falsely identified in the background (false positives).

## Planned Repository Structure (Initial Sketch)
```
/ BINF6251-final-project
├── data/
│   ├── hg38_chr22.fa            # Chromosome 22 FASTA from UCSC
│   ├── ucsc_cpg_islands.bed     # CpG island BED file
│   └── prototype_genome.fa      # Synthetic sequence
├── hmm.py                       # HMM and Viterbi implementation
├── data_readers.py              # FASTA/BED parsing and cleaning
├── evaluate.py                  # Comparison logic and validation
├── .gitignore                   # Excludes large data files
├── LICENSE                      # Project licensing
├── PROPOSAL.md                  # Project proposal
└── README.md                    # Project overview
```

## Generative AI Disclosure (If Used)
Claude was used to understand the different algorithms and biological contexts. 
