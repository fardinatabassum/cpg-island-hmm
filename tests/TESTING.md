# Testing and Validation Report
This document outlines the multi-layered testing strategy used to ensure the reliability and scientific accuracy of the HMM-based CpG island detector.

## Automated Unit Testing
The project includes a suite of automated tests located in the `tests/` directory to verify fundamental components without requiring large genomic datasets.

### `tests/test_viterbi.py`
* Verifies that AT-rich sequences correctly remain in the Background state (State 0).
* Confirms that the Viterbi algorithm correctly processes log-probabilities to prevent numerical underflow.
* Ensures that the logic preserves exact genomic coordinates even when encountering long assembly gaps (N characters).

### `tests/test_data_reader.py`
* Confirms that only the chromosome ID (such as chr22) is extracted, ignoring other metadata in the FASTA header.
* Verifies that lowercase sequences are correctly converted to uppercase for HMM compatibility.
* Ensures that multi-line FASTA sequences are correctly joined into a single continuous string and newline is removed.

**To run the test:**
```
python3 -m unittest discover -v tests
```


## Prototype Validation
Before analyzing Chromosome 22, the algorithm was validated against a synthetic messy genome where the island locations were known.

* **Test Data**: Generated via `scripts/generate_test_data.py`, creating a sequence with high-GC islands, AT-rich background, and assembly gaps.
* **Goal**: The model must detect the primary island while ignoring any synthetic decoys (high GC but less than 200bp) due to the length filter.
* **Results**:
    * Expected: (120, 420)
    * Detected: `[(118, 420)]`
* **Success**: The model successfully ignored short fragments and bridged a 60bp gap to find the correct island.


##  Biological Sanity Checks 
The algorithm satisfies several critical sanity checks to ensure the output is biologically meaningful.

* The total length of the Viterbi path always matches the input sequence length, ensuring predicted coordinates align with reference data.
* There is a minimal length filter that ensures that every entry in the final `.bed` output is $\ge$ 200bp, adhering to the standard biological definition of a CpG island.


## Performance Comparison on Chromosome 22
The model was evaluated against the **UCSC cpgIslandExt** track for Human Chromosome 22 (**GRCh38**).

| Metric         | Result                       |
|:---------------|:-----------------------------|
| Sensitivity    | **75.85%**                   |
| Precision      | **57.75%**                   |
| True Positives | 559 unique islands recovered |

The **75.85%** sensitivity was achieved using a reference-first overlap logic. This ensures that if the HMM fragments a single biological island into multiple predictions, it is only counted as a single success which provides a more realistic picture.


## Hardware
Validation was performed on a standard 8GB RAM laptop. To ensure memory efficiency, the pipeline utilizes local slicing for visualizations and log-space arithmetic for global decoding.
