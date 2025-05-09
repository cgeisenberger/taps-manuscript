# Single-cell Epi2-seq

This is the repository accompanying the Manuscript _Single-cell multi-omic detection of DNA methylation and histone modifications reconstructs the dynamics of epigenetic maintenance_. 
It contains links to all the (raw) sequencing data, FACS index data as well as information about the preprocessing and downstream analyses. 

## Data availability

### Single-cell Epi2-seq data
Raw data for combined methylation and histone profiling have been deposited at the Gene Expression Omnibus (accession [GSE232637](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE232637)). 

### RPE1-hTERT WGBS data
Raw WGBS data have been deposited at the Gene Expression Omnibus (accession [GSE232637](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE232637)).  

### Whole-genome bisulfite sequencing (WGBS)
The following files were retrieved from ENCODE:

* K562
  - [ENCFF867JRG.bed.gz](https://www.encodeproject.org/files/ENCFF867JRG/@@download/ENCFF867JRG.bed.gz)
  - [ENCFF721JMB.bed.gz](https://www.encodeproject.org/files/ENCFF721JMB/@@download/ENCFF721JMB.bed.gz)
* HepG2
  - [ENCFF817LMT.bed.gz](https://www.encodeproject.org/files/ENCFF817LMT/@@download/ENCFF817LMT.bed.gz)
  - [ENCFF453UDK.bed.gz](https://www.encodeproject.org/files/ENCFF453UDK/@@download/ENCFF453UDK.bed.gz)
* H1
  - [ENCFF573YXL.bed.gz](https://www.encodeproject.org/files/ENCFF573YXL/@@download/ENCFF573YXL.bed.gz)
  - [ENCFF434CNG.bed.gz](https://www.encodeproject.org/files/ENCFF434CNG/@@download/ENCFF434CNG.bed.gz)
* GM12878
  - [ENCFF614QHA.bed.gz](https://www.encodeproject.org/files/ENCFF614QHA/@@download/ENCFF614QHA.bed.gz)
  - [ENCFF570TIL.bed.gz](https://www.encodeproject.org/files/ENCFF570TIL/@@download/ENCFF570TIL.bed.gz)

### K562 ChIP-Seq data
1. [H3K27me3 dataset](https://www.encodeproject.org/experiments/ENCSR000EWB/)
  * [H3K27me3 replicate 1](https://www.encodeproject.org/files/ENCFF190OWE/@@download/ENCFF190OWE.bam)
  * [H3K27me3 replicate 2](https://www.encodeproject.org/files/ENCFF692KQZ/@@download/ENCFF692KQZ.bam)
2. [H3K36me3 dataset](https://www.encodeproject.org/experiments/ENCSR000DWB/)
  * [H3K36me3 replicate 1](https://www.encodeproject.org/files/ENCFF639PLN/@@download/ENCFF639PLN.bam)
  * [H3K36me3 replicate 2](https://www.encodeproject.org/files/ENCFF673KBG/@@download/ENCFF673KBG.bam)
3. [H3K9me3 dataset](https://www.encodeproject.org/experiments/ENCSR000APE/)

### K562 ChIC data
[Count tables](https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE164779&format=file&file=GSE164779%5Fmetadata%5FK562%5Fk9me3%2Etxt%2Egz) for single-cell sortChIC data from [Zeller et al.](https://www.nature.com/articles/s41588-022-01260-3) were downloaded from GEO repository [GSE164779](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE164779).

## Preprocessing 
Preprocessing of fastq.gz files was performed using the Package `SingleCellMultiOmics` developed by [Buys de Barbanson](https://github.com/BuysDB). The software package is available on [github](https://github.com/BuysDB/SingleCellMultiOmics).

## Code for figures
All code to produce the figures in the manuscript can be found in the subdirectory `code`.

## Adapter sequences 
Supplementary files containing all adapter and oligo sequences used for sequencing library preparation are summarized in three supplementary tables in the subdirectory `tables`.

## Index Data
FACS plots for the experiments in RPE-1 hTERT Fucci cells where deposited in the subdirectory `facs_plots`. The subfolder `index_data` contains FACS parameters for every single cell, file names refer to the data available via GEO.
