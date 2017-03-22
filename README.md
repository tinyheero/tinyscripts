# tinyscripts

Stores general (tiny) scripts that have not been packages up yet.

## summarize_rnaseq_reads_by_exon.R

This is a R script focused on generate exon raw counts and expression values. It requires as input a transcript database .sqlite file. This can be generated from the `GenomicFeatures::makeTxDbFromBiomart()` function. For example:

```r
txdb <- makeTxDbFromBiomart(dataset = "hsapiens_gene_ensembl" )
saveFeatures(txdb, "~/hsapiens_txdb.sqlite")
```

This can be then be used as input into the `summarize_rnaseq_reads_by_exon.R` script. For example:

```
Rscript summarize_rnaseq_reads_by_exon.R \
  ~/hsapiens_txdb.sqlite \
  /path/to/bam \
  /path/to/outfile
```
