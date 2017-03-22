#!/usr/bin/env Rscript
# Description: This script is used to generate exon raw counts and expression values 
# Authors: Fong Chun Chan <fongchunchan@gmail.com>
library("GenomicFeatures")
library('GenomicAlignments')
library("Rsamtools")
library('argparse')

parser <- ArgumentParser(description='This script is used to generate exon raw counts and expression values')
parser$add_argument('--addChr', action = 'store_true',  help = 'Set the flag to add chr as a prefix to each seqlevel (chromosome name) in the txdbFile. This is needed if you are using an aligned bam file with chr prefixed and the database doesn\'t have a chr prefix [default %(default)s]')
parser$add_argument('--chr', type = 'character',  help = 'Set the chromosome you want to load the data and summarize over. If not set, then it will summarize over the entire genome. [default %(default)s]')
parser$add_argument('--ignoreStrandFlag', action = 'store_true',  help = 'Set this if your library was prepared without any strand-specific protocol [default %(default)s]')
parser$add_argument('--singleEndFlag', action = 'store_true',  help = 'Set whether the library was single-end sequencing [default %(default)s]')
parser$add_argument('--fragmentsFlag', action = 'store_true',  help = 'If the library is paired-end, setting this option means that reads with an unmapped paired will be used in the counting. By default, only mapped paired reads are included in the counting [default %(default)s]')
parser$add_argument('txdbFile', nargs = 1, type = 'character', help = 'The transcript database file (sqlite file)')
parser$add_argument('bamFile', nargs = 1, type = 'character', help = 'The input bam file')
parser$add_argument('outFile', nargs = 1, type = 'character', help = 'The output file')
arguments <- parser$parse_args()

#For debugging
if (FALSE){
	opt <- c( 
			 '--addChr',
			 '--chr', 'chr15'
			 )
	txdbFile <- '~/ensg72.biomart.sqlite'
	bamFile <- 'RCOR1_KD/WTSS/bam/HEK293_RCOR1_clone1.bam'
	outFile <- 'tmp.txt'
	arguments <- parser$parse_args( c( txdbFile, bamFile, outFile, opt ) )
}

print(paste("Loading", arguments$txdbFile))
txdb <- loadDb(arguments$txdbFile)

if( arguments$addChr && !is.null(arguments$chr) ){ # if we need to add a chr prefix to the txdb (e.g. if input chromosome is chr and the txdb does not have chr), then we need to actually remove the chr in the chromosome argument since I can't figure out how to add chr to the txdb 
	print( paste('Restricting the transcript database to chromosome -', arguments$chr))
 	seqlevels(txdb, force = TRUE) <- gsub('chr', '', arguments$chr )
 	#txdb <- restoreSeqlevels(txdb) # restores all the seqlevels (essentially removes the chromosome filter)
} 
print('... Done')

print(paste("Reading", arguments$bamFile, "and retrieving exons... "))
si <- seqinfo(BamFile(arguments$bamFile))
gr <- GRanges(seqnames(si), IRanges(100, seqlengths(si)-100))

# allExons is read in as a GRanges object and NOT as GRangesList. This has implications on how the summarizeOverlaps() function reads counts. 
allExons <- exons(txdb, columns = c('gene_id', 'exon_id', 'exon_name'))

if ( arguments$addChr ){
	print('Prefixing chr to the chromosome names for exons...')
	newSeqNames <- paste('chr', seqlevels(allExons), sep = '')
	names(newSeqNames) <- seqlevels(allExons)
	allExons <- renameSeqlevels( allExons, newSeqNames )
}

scf <- scanBamFlag( 
					isNotPrimaryRead = FALSE, # remove non-primary reads (i.e. multi-aligned reads) 
					isNotPassingQualityControls = FALSE, # remove reads that fail quality control
					)

if ( arguments$singleEndFlag ){
	print('Reading in reads as single-end')
	reads <- readGAlignments( BamFile( arguments$bamFile, yieldSize = 1000000 ), param = ScanBamParam( which = gr, flag = scf ) )
} else{
	print('Reading in reads as paired-end...')
	reads <- readGAlignmentPairs( BamFile( arguments$bamFile, yieldSize = 1000000 ), param = ScanBamParam( which = gr, flag = scf ) );
}
print('Finished')

print('Count raw exon read counts ...')
# Due to each exon being a distinct feature (GRanges object), by default the summarizeOverlap() modes of counting should technically discard any reads that are split between two exons (e.g. split-reads) or have pairs that align to two different exons. 
# This is different from the GRangesList in wich the exons would have been related through another feature (e.g. gene) and so the counting modes would take this into account. We circument this in the summarizeOverlaps, 
# by setting the inter.feature parameter to FALSE. This ensures that reads that overlap multiple distinct exons (e.g. split-reads, and paired reads that align to two different exons) will be counted as contributing 1 read to each exon.
# fragments: If set to TRUE, will count reads have an unmapped pair.
summarizedExpt <- summarizeOverlaps(allExons, reads, ignore.strand = arguments$ignoreStrandFlag, inter.feature = FALSE, fragments = arguments$fragmentsFlag )
countsForExons <- as.numeric( assays(summarizedExpt)$counts )
names(countsForExons) <- rownames(summarizedExpt)
print('... Done')

print('Generating expression values ...')
numBases <- width(allExons)
numKBases <- numBases / 1000
millionsMapped <- length(reads) / 10^6
rpm <- countsForExons / millionsMapped
rpkm <- rpm / numKBases
print('... Done')

print('Retrieving annotation data ...')
annotDf <- values(allExons)
print('...Done')

print('Preparing output...')
exonsReadDf <- data.frame( geneID = sapply(annotDf[, 'gene_id'], '[[', 1), exonID = annotDf[, 'exon_id'], exonName = annotDf[, 'exon_name'], exonCount = countsForExons, exonRPM = rpm, exonRPKM = rpkm, stringsAsFactors = FALSE )
print('...Done')

print(paste('Writing data to', arguments$outFile))
write.table(exonsReadDf, file = arguments$outFile, sep = '\t', quote = F, row.names=F)
print('...Done')
print(proc.time())
