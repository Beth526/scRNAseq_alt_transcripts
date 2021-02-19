# scRNAseq_alt_transcripts


The goal was to see if alternate poly-A usage could be estimated from 3' Chromium scRNA-seq datasets, and if it would be informative for clustering or cluster markers. Work in progress!

Test data used: 10x Chromium 10,000 PBMC dataset (Single Cell 3’ v.3, CellRanger v.4.0.0). Both the filtered .h5 file and the CellRanger-processed .bam file are needed for this analysis.

Note: CellRanger v4 uses STAR to align the reads, for a read to be counted for a gene it must overlap an exon by 50% AND be consistant with an annotated transcript. CellRanger v5 has a new feature to count reads from introns and all reads from the correct strand across the length of the gene. CellRanger bam files also contain corrected UMI and barcodes.

The functions are in the sc_altpolya.py file 

```{python}
filtered_barcodes, top_genes = get_top_genes_and_barcodes_list(h5_file_name, n_genes=100)
```
Gets the top genes by counts and the filtered barcodes list from the 10X filtered .h5 file. Probably the first 100 will be mostly mitochondrial and ribosomal proteins.

![top genes head](https://github.com/Beth526/scRNAseq_alt_transcripts/blob/main/images/top_gene_table.png)

```{python}
df = gtf_top_genes(top_genes, gtf_file)
```
Gets chrom, strand and exon positions set for each gene in top_genes or any dataframe with ensembl gene ids of interest in an 'ids' column

![genes with gtf info head](https://github.com/Beth526/scRNAseq_alt_transcripts/blob/main/images/gtf_table.png)

```{python}
df = samtools_view(df, filtered_barcodes, bam_file)
```
Considering gene strand, finds all reads that start between the first exon position and the last. Finds read stop positions and only retains furthest 3' read-stop for each UMI group. Allows for read_stops to occur up to 300kbp after last exon position.

![genes with read_stop info head](https://github.com/Beth526/scRNAseq_alt_transcripts/blob/main/images/samtools_table.png)

```{python}
counts, summary = seperate_into_peaks(df, num_peaks_exon = 3, num_peaks_other = 3)
```
Uses a gaussian mixture model to find peaks (up to the number specified) of read_stops for 
1) exons: read_stops occuring inside exons. The introns are removed and exons juxtaposed, so that the peaks are detected along the coding sequence of the gene.
2) others: read_stops occuring outside exons. This includes introns, and genomic sequence outside the gene, up to 300kbp downstream.

![counts table head](https://github.com/Beth526/scRNAseq_alt_transcripts/blob/main/images/counts_table.png)

![summary table head](https://github.com/Beth526/scRNAseq_alt_transcripts/blob/main/images/summary_table.png)

```{python}
counts, summary = select_alt_transcripts(final_counts, final_summary, min_per_gene = 2, max_within_gene_correlation = 1, count_greater_than_std = True)
```
Filter the counts and summaries dataframes. It can remove peaks if less than 2 were found for a gene (a single alt transcript will be highly correlated with parent gene counts), remove peaks if counts was less than standard deviation (wide sparse peaks are more likely to be noise) and remove peaks if minimmum intra gene correlation was above a certain level (not informative).

```{python}
corr_matrix = peak_correlations_for_gene(base_gene_id, counts)
```
A correlation matrix for all alternate transcripts of a given gene

```{python}
min_corr_table = min_correlations_by_gene(counts, summaries)
```
A table of the minium intragene correlations for alternate transcripts of every gene

```{python}
graph_exons(gene_name, df, num_peaks = 3, num_bins = 100)
```
Produce 3 histogram graphs for read_stops inside exons before and after assignment by GMM algorithm. Note that the GMM model is non-deterministic so sometimes the results will difer. Running this function repeatedly with a different num_peaks may give an idea of the best num_peaks parameter for the seperate_into_peaks() function.

3 graphs
1) read_stops outside exons genomic location histogram
2) read_stops outside exons genomic location histogram, x-axis limited to gene 
3) num_peaks histograms for each of the peaks found 

![example](https://github.com/Beth526/scRNAseq_alt_transcripts/blob/main/images/GNLY.jpg)


```{python}
 graph_others(gene_name, df, num_peaks = 3, num_bins = 100)
```
Produce 3 histogram graphs for read_stops outside exons before and after assignment by GMM algorithm. 

3 graphs
1) read_stops outside exons genomic location histogram
2) read_stops outside exons genomic location histogram, x-axis limited to gene 
3) num_peaks histograms for each of the peaks found 

```{python}
new_h5(counts, summaries, h5_dset_file, new_dset_file)
```
Make a new .h5 file containing all the information from the original .h5 file but with the new alternate trascript counts added.


The genes ranked top 100 - 1100 from the PBMC dataset were put through this pipeline and then an h5 file with the alternate transcripts included was analyzed using Seurat alongside the original h5 file.

The additional correlated genes biased the principal components, giving more power to genes with many highly correlated transcripts in clustering. 

However, some of the alternate transcripts were unique cluster markers, meaning that thier parent gene and other alternate transcripts were not also a marker for the cluster. There were on average 19 alternate transcript cluster markers per cluster.

One example is RAB7A
![RAB7A](https://github.com/Beth526/scRNAseq_alt_transcripts/blob/main/images/RAB7A.jpg?raw=true)
![RAB7A_alt](https://github.com/Beth526/scRNAseq_alt_transcripts/blob/main/images/RAB7A_others_2.jpg)

The most interesting finding was a TRAC intronic read_stop (potential alternate poly A site) site highly uncorrelated with other TRAC transcripts
![TRAC_exon](https://github.com/Beth526/scRNAseq_alt_transcripts/blob/main/images/TRAC_exons.jpg)
![TRAC_intron](https://github.com/Beth526/scRNAseq_alt_transcripts/blob/main/images/TRAC_other.jpg)





