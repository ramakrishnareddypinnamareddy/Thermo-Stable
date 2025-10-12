# Thermo-Stable

# Protein Classifier Setup

## Prerequisites

This repository requires the **FoldX** executable to run the `run_classifier` script.

### Steps to Setup FoldX

1. **Download FoldX**  
   Get the valid FoldX executable from the official website: [FoldX Suite](https://foldxsuite.crg.eu/).

2. **Place the Executable**  
   Copy the downloaded FoldX file into the same folder where `run_classifier.py` is located.

3. **Rename the Executable**  
   Ensure the FoldX file is renamed to `foldx` exactly (case-sensitive).

4. **Give Executable Permissions**  
   Run the following command to make it executable:
   ```bash
   chmod +x foldx
