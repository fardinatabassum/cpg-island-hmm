## Final Project Reflection Document
### **What Went Right**
* **Technical and Process Successes**
  * Was successfully able to implement transition to log-space Viterbi which eliminated the numerical underflow issues.
  * The neutral-emission logic for handling assembly gaps ('N' characters) successfully preserved the linear genomic scale. This ensured that our predictions remained perfectly aligned with the UCSC reference.
  * Was able to validate using a synthetic prototype genome with a decoy segment which allowed to verify the 200bp biological filter in a controlled environment before scaling to real data.
* **Personal Learning Successes**
  * I successfully moved beyond the Jupyter Notebook environment to create a modular Python package with a functional Command Line Interface (CLI).
  * Learning to write unit tests for biological edge cases (like case normalization and multi-line FASTA parsing) made the final pipeline significantly more fault-tolerant.
  * I was also able to generate a substantial plot that allowed visual confirmation of the model’s state transitions.
### **What Went Wrong or Was Hard**
  * **Concrete Challenges** 
    * Initial attempts to parse the UCSC cpgIslandExt reference data were difficult due to header metadata. I underestimated the time needed to sanitize real-world genomic files.
    * A recurring challenge was fragmentation, where the HMM would drop out of the Island state for a single base and then jump back in. This resulted in one biological island being reported as two or three smaller segments.
    * I had had to manually adjust my model probabilities to get the outcome I was expecting.
  * **Start Over Strategy**
    * If I were to start over, I would build the evaluation.py script before the Viterbi decoder. Having a way to measure sensitivity immediately would have helped me adjust the transition probabilities much faster. 
  *  **Personal Learning Obstacles**
      * Understanding the transition from probability multiplication to log-addition was conceptually difficult at first. I addressed this by manually calculating a 5-base sequence on paper to see how the numbers transformed. If I couldn't grasp it, I would have visualized the Trellis using a smaller transition matrix to debug the scoring logic.
### **Algorithmic Lessons**
  * The HMM is a natural fit for genomic sequence analysis because biological states (like islands vs. background) are inherently sequential. It captures the transition between these states better than a simple sliding-window GC-count. It was also the least intimidating one to implement but also simple which I why I had to come up with ways to make my project innovative such as implementing a CLI and a visualization.
  * **Tradeoffs**
    * The model achieved a high 75.85% sensitivity, but the precision (57.75%) reflects the classic tradeoff of probabilistic models which are good ar finding features but can also end up picking high-GC regions that lack biological regulatory function.
    * Implementing log-space adds code complexity, but it was a non-negotiable requirement for analyzing human-scale chromosomes.
  * In the lecture, Viterbi is often presented with four static bases (A, C, G, T). The surprising thing during implementation was realizing that real-world genomes are full of 'N' characters and lowercase, which requires significant pre-processing logic not usually found in lecture examples. Working with very messy data causes extra steps to be taken before running the model.
### **Future Directions**
  * For the future I would like to move from a 2-state model to a 4-state model (adding more states) could help reduce the fragmentation issue.
  * Now that we have discussed Baum-Welch in class I would use that to learn the transition and emission probabilities directly from the data, rather than using the manually calibrated values. 
  * I would also like to expand my project and test the model on different chromosomes such as Chromosome X to see if the parameters hold up across different genomic contexts.
### **Generative AI Disclosure (If Used)**
* Appendix: Generative AI Usage
  * **Tool:** Claude
  * **Role in implementation:** Generative AI was utilized for code scaffolding, debugging logical edge cases, testing, refreshing old python logic such as unit testing, Markdown formatting and engineering synthetic test data. The tool influenced the implementation in the following key areas:
    * Synthetic Sequence Generation:**Prompt** "Generate a sequence for a FASTA file that includes an AT-rich background, a high-GC island, and a 50bp assembly gap to test my model"
    * Debugging & Edge Case Handling: **Prompt** "Explain why my precision is lower than my sensitivity in this specific CpG island detector."
    * Unit testing: **Prompt** "Remind me how to write unit tests"
    * Markdown formatting: **Prompt** "How do you create tables on Markdown, help me format letters in an equation on Markdown"