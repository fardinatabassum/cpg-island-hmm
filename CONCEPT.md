## BINF6251 Final Project
# Conceptual Progress Report

## 1. Recap of Project
* **Research Question:** Can a first-order Hidden Markov Model (HMM) accurately identify CpG island boundaries in human promoter regions on Chromosome 22?
* **Algorithm and Algorithm Class:** Hidden Markov Model (HMM) using the Viterbi Algorithm for decoding.
* **Data Type:** Genomic sequence data (FASTA) and coordinate-based annotations (BED).
* **Project Summary:** This project implements a probabilistic framework to infer hidden biological states (CpG Islands vs. Genomic Background) from observable nucleotide sequences.

    By utilizing the Viterbi algorithm, the model accounts for transition dependencies between nucleotides that simple window-based counting methods often overlook allowing for more accurate identification of genomic boundaries.

## 2. Inputs, Outputs, and Assumptions
* **Expected Input(s):**
    * `sequence`: A string of uppercase characters (A, C, G, T) representing the DNA.
    * `transition_probs`: A dictionary of log-probabilities for state transitions.
    * `emission_probs`: A dictionary of log-probabilities for nucleotide emissions per state.
* **Expected Output(s):**
    * `path`: A list of hidden states (Background or Island) of length $L$ representing the most likely path.
    * `predictions.bed`: A standard bed format file containing start and end coordinates of only the predicted islands that that meet the $\ge$ 200bp threshold. 
    * **Confusion Matrix:** A $2 \times 2$ table visualizing the overlap between predicted island bases and known island bases to identify if the model is biased toward the "Background" state.
* **Key Assumptions:**
    * All input sequences are pre-processed to uppercase A/C/G/T.
    * Emission and transition probabilities are normalized and provided in log-space.
    * The sequence is treated as a continuous linear string.
    * **Clarity Threshold:** Any predicted island shorter than 200bp will be discarded to ensure results meet biological standards for clear CpG islands.
    * **Coordinate Mapping:** The output coordinates must include a `global_offset` to account for any sequence chunking performed during the data-reading stage.

## 3. Detailed Pseudocode


```
# Core Function: viterbi_decode(sequence, states, trans_p, emiss_p):
    """
    Purpose: Decodes a DNA sequence into hidden states (Island vs Background).
    Input: Sequence must be a pre-processed string (A,C,G,T,N) from FASTA.
    Data cleaning: Case-normalization and header stripping done by the data_reader.py.
    """
# Data Structures: Trellis matrix (V) and Backpointer matrix (B)

# STEP 1: Initialization
Initialize V matrix [states x length] with -infinity
Initialize B matrix [states x length] with Null
for each state s in states:
    V[s, 0] = log(start_p[s]) + emiss_p[s, sequence[0]]

# STEP 2: Main Logic Loop (Recursive Step)
for i from 1 to sequence_length - 1:
    current_obs = sequence[i]
    
    # Edge Case: Assign neutral probability to assembly gaps ('N')
    if current_obs not in ['A', 'C', 'G', 'T']:     
        for s in range(num_states):
          # Add neutral log-prob and carry forward previous state's probability
          V[s, i] = V[s, i-1] + math.log(0.25)
          B[s, i] = s  # Pointer points to itself at the previous step
        continue # Move to next base (i+1)
      
    # Standard Viterbi Logic (for A, C, G, T)
    for each current_state s:
        V[s, i], B[s, i] = max(V[prev_s, i-1] + trans_p[prev_s, s] for prev_s in states)
        V[s, i] += emiss_p[s, current_obs]
    
# STEP 3: Path Reconstruction (Backtracking)
best_final_state = index of max value in V at last column
path = [best_final_state]

for i from length-1 down to 1:
    # Use the first element of our list (the state at i) 
    # to find the state at i-1
    prev_state = B[path[0], i]
    path.insert(0, prev_state)
 
# STEP 4: Post-Processing & Thresholding
# Filter out unclear/short islands 
Final_Islands = []
for each contiguous segment of 'Island' states in path:
    if length(segment) >= 200:
        # Convert local index to genomic coordinates
        # Start = segment_start + global_offset
        Final_Islands.append(segment_coordinates)
    else:
        # Segments < 200 are ignored (reverted to Background)

return Final_Islands
```
## 4. Complexity and Bottlenecks
* **Time Complexity:** $O(L \times S^2)$, where $L$ is sequence length (~51M for Chromosome 22) and $S$ is the number of states (equal to 2). This is linear relative to the genome size.
* **Space Complexity:** $O(L \times S)$ for storing the Trellis and Backpointer matrices.
* **Performance Bottlenecks:** Storing two matrices for 51 million bases may exceed 8GB of RAM (personal computer limitation), causing performance issues.
* **Mitigation Strategy:** I will implement sequence chunking to process the chromosome in 1MB segments between assembly gaps to reduce memory pressure, or use NumPy arrays with specific data types (int8 for backpointers) to save space.

## 5. Validation and Testing Plan
* **Small simulated Example:**

  * **Input:** "CGCGCGCGAA"

  * **Expected Result:** The first 8 bases should be labeled "Island" due to high C/G density; the final "AA" triggers a transition to "Background".

* **Stress Test Dataset:** Running the full Chromosome 22 FASTA (~51M bp) to ensure the implementation handles realistic scales.

* **Evidence of Incorrectness:**

  * If the output path length does not match the input sequence length.

  * If the model predicts an "Island" state for a long run of Adenines (Poly-A tail).

* **Automated Tests:**

  * **Unit Test:** Verify the log-transformation function handles probabilities near 0.

  * **End-to-End:** Check that the generated BED file coordinates are valid (start < stop).

## 6. Updated Pitfall and Risk Log

### Data-related issues
* **Status:** Persistent & Refined.
* **Current Context:** While the proposal suggested a separate pre-processing cleaning step to handle "N" characters (assembly gaps), the design has shifted toward algorithmic handling.
* **Mitigation:** The pseudocode includes an explicit conditional check to assign neutral probabilities to gaps, preventing script crashes while maintaining genomic coordinates. 

### Algorithmic issues
* **Status:** Mitigated by Design.
* **Current Context:** The risk of Numerical Underflow due to multiplying long chains of probabilities across 51 million base pairs was the primary concern in the proposal.
* **Mitigation:** This is now fully addressed in the pseudocode by implementing the algorithm entirely in log-space. By converting probabilities to logarithms, the implementation uses addition rather than multiplication, ensuring numerical stability throughout the decoding process.

### Evaluation issues
* **Status:** Persistent & Expanded.
* **Current Context:** There is an imbalance problem as CpG islands are rare compared to the genomic background so it is a risk for accurate validation.
* **Mitigation:** To ensure the model doesn't simply predict background for the entire chromosome and to address the rarity of CpG islands, I will use Sensitivity and Specificity instead of raw accuracy. Additionally, I have introduced a minimum length threshold (200bp) for predicted islands to filter out stochastic noise and ambiguous regions to ensure biological relevance.

### New Risk: Memory Management
* **Status:** New. 
* **Context:** Complexity analysis identified that storing matrices for Chromosome 22 could exceed standard RAM limits. 
* **Mitigation:** The implementation will utilize sequence chunking to manage the memory footprint effectively. 

## 7. Appendix: Generative AI Usage
Claude was used to understand the algorithm, work out computational limitations and format the Markdown file.
* Prompts: 
  * Explain fully and step by step how the Viterbi Algorithm works - this helped me fully breakdown the algorithms and understand the concept step by step
  * How to format and categorize in Markdown file? -Used to format headings and titles in Markdown file
  * What are some computational limitations when running Viterbi and how to solve them - Used to understand time complexities and memory limitations and suggest mitigations.
  
