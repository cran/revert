## ----setup, include = FALSE---------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>",
  cache = TRUE,
  eval = TRUE
  )

## ----del_example, eval=FALSE--------------------------------------------------
#  library(revert)
#  
#  bam.file1 <- system.file('extdata', 'toy_alignments_1.bam', package = 'revert')
#  
#  reversions <- getReversions(
#  	bam.file = bam.file1,
#  	genome.version = "BSgenome.Hsapiens.UCSC.hg38",
#  	chromosome = "chr17",
#  	pathog.mut.start = 43082434,
#  	pathog.mut.type = "SNV",
#  	snv.reference.allele = "G",
#  	snv.alternative.allele = "A",
#  	flanking.window = 100,
#  	minus.strand = TRUE
#  	)

## ----ins_example, eval=FALSE--------------------------------------------------
#  bam.file2 <- system.file('extdata', 'toy_alignments_2.bam', package = 'revert')
#  
#  reversions <- getReversions(
#  	bam.file = bam.file2,
#  	genome.version = "BSgenome.Hsapiens.UCSC.hg38",
#  	chromosome = "chr13",
#  	pathog.mut.start = 32338763,
#  	pathog.mut.type = "DEL",
#  	deletion.sequence = "AT",
#  	deletion.length = 2,
#  	flanking.window = 100
#  	)

## ----snv_example, eval=FALSE--------------------------------------------------
#  bam.file3 <- system.file('extdata', 'toy_alignments_3.bam', package = 'revert')
#  
#  reversions <- getReversions(
#  	bam.file = bam.file3,
#  	genome.version = "BSgenome.Hsapiens.UCSC.hg38",
#  	chromosome = "chr17",
#  	pathog.mut.start = 43092689,
#  	pathog.mut.type = "INS",
#  	insertion.sequence = "T",
#  	flanking.window = 100,
#  	minus.strand = TRUE
#  	)	

