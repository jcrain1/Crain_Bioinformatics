# Bioinformatics Purpose

Document for sequencing data, the progression from bcl or non-demultiplexed fastq files--to documented fastq files that can be used for analysis.

## Step 1 Demultiplex

Demultiplex all files. Scripts for demultiplexing are kept here as well as within bulk data.

### Demultiplex

Using bcl-convert, develop scripts and sample files to demultiplex data. Includes queries and script files.

### For Beocat

Directory structure was set up for each skimsequencing run. Structure is:

```
├── SkimSeq_Project_Number
   ├── 00_Demultiplex
      ├── Individual samples demultiplexed
      ├── Summary Statistics
   ├── 01_Samples
      ├── Concatenated samples into single samples.
   ├── 02_trimmed
      ├── 
   ├── 03_aligned
      ├── 
   ├── 04_filtered
      ├── 
   ├── 05_indexed
      ├── 
   ├── Scripts
      ├── Shell scripts that step through process
   ├── out_error
      ├── 
```