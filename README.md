Standard-bulk-RNA-seq-pipeline
---
The assignment involved analyzing generic FASTQ files through the following steps:

1. Read Preprocessing: Performed read trimming using Sickle and Scythe.
2. Read Alignment and Assembly: Utilized the 'New Tuxedo' pipeline with HISAT2 for read alignment and StringTie for transcript assembly.
3. Count Matrix Generation: Generated count matrices from the assembly using an external prepDE Python script.
4. Experimental Design Table: Created based on provided instructions.
5. Data Filtering and Differential Expression Analysis: Applied filtering and ran differential expression analysis.
6. Visualization:
  -Dispersion plots
  -rlog-transformed PCA plots
  -“SD versus Mean” plots
  -MA (MvA) plots based on null hypotheses
7.Batch Effect Correction: Applied a limma-based batch removal function to correct rlog-transformed gene counts.
8.Revised PCA: Generated PCA plots using batch-corrected rlog values.
