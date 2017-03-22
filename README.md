# tinyscripts

Stores general (tiny) scripts that have not been packages up yet.

## summarize_rnaseq_reads_by_exon.R

This is a R script focused on generate exon raw counts and expression values. It requires as input a transcript database .sqlite file. This can be generated from the `GenomicFeatures::makeTxDbFromBiomart()` function. For example:

```r
library("GenomicFeatures")
txdb <- makeTxDbFromUCSC(genome = "hg19", tablename = "knownGene")
saveDb(txdb,  "~/hsapiens_txdb.sqlite")
```

This will take some time depending on what database you are trying to download. For reasons, not completely known to me, the `.sqlite` file has to be saved in your home directory (hence why the example saves it to `~/hsapiens_txdb.sqlite`. Once this has been generated, you can use it as input into the `summarize_rnaseq_reads_by_exon.R` script. For example:

```bash
Rscript summarize_rnaseq_reads_by_exon.R \
  ~/hsapiens_txdb.sqlite \
  /path/to/bam \
  /path/to/outfile
```
