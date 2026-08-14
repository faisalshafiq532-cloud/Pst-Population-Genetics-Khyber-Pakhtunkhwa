# STRUCTURE analysis files

This directory contains files associated with the visualization of population structure results.

Population structure was inferred using the external software STRUCTURE. Replicate alignment and processing were conducted separately using CLUMPP. The R script `12_STRUCTURE_Plotting.R` was used to process and visualize the resulting population membership data.

The files in this directory are retained to support reproducibility of the STRUCTURE visualization workflow.

## Files

- `Map.xlsx`: Population or isolate metadata used in the STRUCTURE visualization workflow.
- `K2ClumppIndFile.output`: CLUMPP-aligned individual membership data for K = 2.
- `K3ClumppIndFile.output`: CLUMPP-aligned individual membership data for K = 3.
- `K4ClumppIndFile.output`: CLUMPP-aligned individual membership data for K = 4.
- `K5ClumppIndFile.output`: CLUMPP-aligned individual membership data for K = 5.
- `K6ClumppIndFile.output`: CLUMPP-aligned individual membership data for K = 6.
- `K7ClumppIndFile.output`: CLUMPP-aligned individual membership data for K = 7.

These files are used as inputs for the R-based visualization workflow and are not final manuscript figures.