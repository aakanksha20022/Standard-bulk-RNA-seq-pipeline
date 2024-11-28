RNA-Seq Analysis Workflow
----
This repository outlines the steps performed for analyzing generic FASTQ files using an RNA-Seq pipeline. Below is a summary of the process:

Workflow Steps
1. Read Preprocessing
   
Read trimming was conducted using Sickle and Scythe.

3. Read Alignment and Assembly
   
Used the 'New Tuxedo' pipeline:

HISAT2: Read alignment

StringTie: Transcript assembly

5. Count Matrix Generation
   
Generated count matrices from the transcript assembly using the external prepDE.py script.

7. Experimental Design Table
   
Created an experiment design table based on the provided instructions.

9. Data Filtering and Differential Expression Analysis
    
Applied filtering to the data and performed differential expression (DE) analysis.

11. Visualization
    
Created the following plots for data exploration and interpretation:

Dispersion plots

rlog-transformed PCA plots

“SD versus Mean” plots

MA (MvA) plots based on null hypotheses

13. Batch Effect Correction
    
Corrected rlog-transformed gene counts using a limma-based batch removal function.

15. Revised PCA
    
Generated updated PCA plots using batch-corrected rlog values.


Outputs
Count Matrices: Gene-level and transcript-level count matrices for downstream analysis.
Visualizations: Dispersion plots, PCA plots, SD vs. Mean plots, and MA plots.
