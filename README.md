# Single-cell Epi2-seq

This is the repository accompanying the single-cell Epi2-seq manuscript by Geisenberger et al. 2023. 


## Data availability
Raw scEPI2-seq data have been deposited at the Gene Expression Omnibus (accession xxx). 

### RPE1-hTERT WGBS data
Raw WGBS data have been deposited at the Gene Expression Omnibus (accession xxx).  

### K562 whole-genome bisulfite sequencing (WGBS)
Data were downloaded from ENCODE:
* [replicate 1 bed file](https://www.encodeproject.org/files/ENCFF867JRG/@@download/ENCFF867JRG.bed.gz)
* [replicate 2 bed file](https://www.encodeproject.org/files/ENCFF721JMB/@@download/ENCFF721JMB.bed.gz)

### K562 ChIP-Seq data
BAM files for H3K27me3 and H3K36me3 were downloaded from ENCODE:
* [H3K27me3 replicate 1](https://www.encodeproject.org/files/ENCFF190OWE/@@download/ENCFF190OWE.bam)
* [H3K27me3 replicate 2](https://www.encodeproject.org/files/ENCFF692KQZ/@@download/ENCFF692KQZ.bam)
* [H3K36me3 replicate 1](https://www.encodeproject.org/files/ENCFF639PLN/@@download/ENCFF639PLN.bam)
* [H3K36me3 replicate 2](https://www.encodeproject.org/files/ENCFF673KBG/@@download/ENCFF673KBG.bam)

### K562 ChIC data
[Count tables](https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE164779&format=file&file=GSE164779%5Fmetadata%5FK562%5Fk9me3%2Etxt%2Egz) for single-cell sortChIC data from [Zeller et al.](https://www.nature.com/articles/s41588-022-01260-3) were downloaded from GEO repository [GSE164779](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE164779).

## Preprocessing 
Preprocessing of fastq.gz files was performed using the Package `SingleCellMultiOmics` developed by [Buys de Barbanson](https://github.com/BuysDB). The software package is available on [github](https://github.com/BuysDB/SingleCellMultiOmics).

## Code for figures
All code to produce the figures in the manuscript can be found in the subdirectory `code`.

## Adapter sequences 
Supplementary files containing all adapter and oligo sequences used for sequencing library preparation are summarized in three supplementary tables in the subdirectory `tables`.
