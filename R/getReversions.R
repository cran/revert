#' Detecting reversion mutations
#' 
#' @description
#' 	Function for detecting reversions for a given pathogenic mutation from reads alignment of NGS genomic data
#' 	
#' @param bam.file
#'  Character. 
#'  The input bam file containing aligned reads.
#' @param genome.version
#'  Character. 
#'  Genome version of alignments in bam file and the position of pathogenic mutation.
#'  The name of genome version is specified by the available BSgenome data packages: BSgenome::available.genomes().
#'  Default is "BSgenome.Hsapiens.UCSC.hg38".
#' @param chromosome
#'  Character. 
#'  Name of chromosome where pathogenic mutation is located, e.g., "chr17" or "17", "chrX" or "X".
#'  The chromosome name should be concordant with the chromosome identifiers used in bam.file and genome.version
#' @param pathog.mut.start
#'  Integer. 
#'  Genomic coordinate of the start position of pathogenic mutation.
#' @param pathog.mut.type
#'  Character. 
#'  Type of pathogenic mutation: "SNV", "DEL" or "INS".
#' @param snv.reference.allele
#'  Character.
#'  Reference allele of pathogenic mutation if it is a single nucleotide variant (SNV).
#'  Default is NULL.
#' @param snv.alternative.allele
#'  Character.
#'  Alternative allele of pathogenic mutation if it is a SNV.
#'  Default is NULL.
#' @param deletion.sequence
#'  Character. 
#'  Deleted nucleotides of pathogenic mutation if it is a deletion (DEL). 
#'  Default is NULL.
#' @param deletion.length
#'  Integer. 
#'  Number of deleted nucleotides of pathogenic mutation if it is a DEL. 
#'  Parameters deletion.sequence and deletion.length can not both be NULL when pathogenic mutation is DEL.
#'  Default is NULL. 
#' @param insertion.sequence
#'  Character. 
#'  Inserted nucleotides of pathogenic mutation if it is an insertion (INS). 
#'  Default is NULL.
#' @param flanking.window
#'  Integer. 
#'  Length of flanking regions (bp) to pathogenic mutation locus for reversion detection. 
#'  Default is 100.
#' @param minus.strand
#'  Logical. 
#'  TRUE if the gene in question is on the reverse strand. 
#'  Default is FALSE. 
#' @param check.wildtype.reads
#'  Logical. 
#'  TRUE if assume wildtype reads mapped to pathogenic mutation are restored to wildtype and alternative reversions will be detected from the wildtype reads.
#'  Only used for pathogenic mutation targeted gene editing experiment.   
#'  Default is FALSE.
#'  
#' @return A list containing two tables summarizing the reversion detection:
#'	1. The numbers of different types of reads analysed in the input bam file
#'	2. Reversion mutation table including the following columns: 
#'      rev_id: Unique ID for reversion event
#'		rev_freq: Frequency of reversion event
#'		rev_type: Type of reversion event, i.e., complement reversion to pathogenic mutation, replacement reversion of pathogenic mutation, or alternative reversion to pathogenic mutation 
#'		rev_mut_number: Index of each mutation in a reversion event
#'		mut_id: Unique ID for reversion mutation 
#'		chr: Chromosome
#'		mut_start_pos: Start position of reversion mutation
#'		mut_type: Type of reversion mutation, i.e., SNV, DEL or INS
#'		mut_seq: Sequence changes of mutation, i.e., inserted or deleted sequences for indels, or reference and alternative alleles for SNVs 
#'		mut_length: Length of mut_seq, 0 for SNV
#'		mut_hgvs: HGVS Genomic DNA ID of reversion mutation
#'		pathog_mut_hgvs: Original pathogenic reversion mutation
#'		dist_to_pathog_mut: Distance between original pathogenic mutation and reversion mutation
#'
#' @importFrom BSgenome installed.genomes
#' @importFrom BSgenome getBSgenome
#' @importFrom BSgenome seqnames
#' @importFrom BSgenome.Hsapiens.UCSC.hg38 BSgenome.Hsapiens.UCSC.hg38
#' @importFrom Rsamtools BamFile
#' @importFrom Rsamtools scanBamHeader
#' @importFrom Rsamtools ScanBamParam
#' @importFrom Rsamtools scanBamFlag
#' @importFrom Rsamtools scanBam
#' @importFrom GenomicRanges GRanges
#' @importFrom IRanges IRanges
#' @importFrom Biostrings getSeq
#' @importFrom Biostrings reverseComplement
#' @importFrom Biostrings DNAString
#' @importFrom Biostrings translate
#' @importFrom Biostrings subseq
#' @importFrom stats quantile
#' 
#' @examples
#' {
#' # To detect reversions for BRCA2 mutation "chr13:g.32338763-32338764delAT"
#' bam.file2 <- system.file('extdata', 'toy_alignments_2.bam', package = 'revert')
#' reversions <- getReversions(
#'     bam.file = bam.file2, 
#'     genome.version = "BSgenome.Hsapiens.UCSC.hg38",
#'     chromosome = "chr13", 
#'     pathog.mut.start = 32338763, 
#'     pathog.mut.type = "DEL", 
#'     deletion.sequence = "AT",
#'     deletion.length = 2,
#'     flanking.window = 100,
#'     minus.strand = FALSE )
#' 
#' }
#' @export getReversions
getReversions <- function(
		bam.file, 
		genome.version = "BSgenome.Hsapiens.UCSC.hg38",
		chromosome, 
		pathog.mut.start, 
		pathog.mut.type = "SNV", 
		snv.reference.allele = NULL, 
		snv.alternative.allele = NULL,
		deletion.sequence = NULL,
		deletion.length = NULL,
		insertion.sequence = NULL,
		flanking.window = 100,
		minus.strand = FALSE, 
		check.wildtype.reads = FALSE ) {
	
	#######################
	# parameter check
	#######################
	
	if (!file.exists(bam.file)) {
		stop("BAM file not found at ", bam.file)
	}else{
	    bamfile <- BamFile(bam.file)
	    bam.header <- scanBamHeader(bamfile)
	    bam.header.chr <- names(bam.header$targets)
	}
	
	if( genome.version %in% installed.genomes() ){
	    reference.genome <- getBSgenome(genome.version)
	}else{
	    stop("Reference genome package ",genome.version," is not currently installed.\nPlease install the package ",genome.version," before running getReversions()")
	}
    
    chromosome <- as.character(chromosome)
    if( ! (chromosome %in% seqnames(reference.genome) ) ){
        stop( "The input chromosome name \"", chromosome,"\" is not concordant with the chromosome names in reference genome ", genome.version)
    }
    if( ! (chromosome %in% bam.header.chr ) ){
        stop( "The input chromosome name \"", chromosome,"\" is not concordant with the chromosome names in bam file ", bam.file)
    }
    
	if( !(pathog.mut.type %in% c("SNV","DEL","INS")) ){
		stop("pathog.mut.type ",pathog.mut.type," is not supported.\nPlease input pathog.mut.type: SNV, DEL or INS")
	}
	
	if( pathog.mut.type=="SNV" ){
		if( is.null(snv.reference.allele) | is.null(snv.alternative.allele) ){
			stop("snv.reference.allele and snv.alternative.allele are required when pathog.mut.type is SNV")
		}else{
			if( nchar(snv.reference.allele)!=1 | nchar(snv.alternative.allele)!=1 ){
				stop("Input SNV is not single nucleotide substitution")
			}else{
				if( !(snv.reference.allele %in% c("A","T","C","G") ) ){
					stop("snv.reference.allele is not a DNA nucleotide of A,C,T,G")
				}
				if( !(snv.alternative.allele %in% c("A","T","C","G") ) ){
					stop("snv.alternative.allele is not a DNA nucleotide of A,C,T,G")
				}
				
				pathog.mut.end <- pathog.mut.start
				true.reference.allele <- as.character( getSeq(reference.genome, chromosome, pathog.mut.start, pathog.mut.end) )
				if( snv.reference.allele != true.reference.allele ){
					stop("snv.reference.allele does not match reference genome")
				}else{
					roi.start <- pathog.mut.start - flanking.window
					roi.end <- pathog.mut.end + flanking.window
					snv.genotype <- paste(snv.reference.allele, snv.alternative.allele, sep=">")
					pathog.mut.hgvs <- paste0(pathog.mut.start,snv.genotype)
				}
			}
		}
		if( !is.null(insertion.sequence) ){
			warning("insertion.sequence is not used when pathog.mut.type is SNV")
		}
		if( !is.null(deletion.sequence) | !is.null(deletion.length) ){
			warning("deletion.sequence and deletion.length are not used when pathog.mut.type is SNV")
		}
	}
	
	if( pathog.mut.type=="DEL" ){
		if( is.null(deletion.sequence) & is.null(deletion.length) ){
			stop("deletion.sequence and deletion.length can not both be NULL when pathog.mut.type is DEL")
		}
		if( !is.null(deletion.sequence) & !is.null(deletion.length) ){
			if( nchar(deletion.sequence)!=deletion.length ){
				stop("deletion.sequence and deletion.length are not concordant")
			}else{
				pathog.mut.end <- pathog.mut.start + deletion.length - 1
				deletion.reference <- as.character( getSeq(reference.genome, chromosome, pathog.mut.start, pathog.mut.end) )
			}
		}
		if( !is.null(deletion.sequence) & is.null(deletion.length) ){
			deletion.length <- nchar(deletion.sequence)
			pathog.mut.end <- pathog.mut.start + deletion.length - 1
			deletion.reference <- as.character( getSeq(reference.genome, chromosome, pathog.mut.start, pathog.mut.end) )
		}
		if( is.null(deletion.sequence) & !is.null(deletion.length) ){
			pathog.mut.end <- pathog.mut.start + deletion.length - 1
			deletion.sequence <- as.character( getSeq(reference.genome, chromosome, pathog.mut.start, pathog.mut.end) ) 
			deletion.reference <- deletion.sequence
		}
		
		if( deletion.sequence != deletion.reference ){
			stop("deletion.sequence does not match reference genome")
		}
		
		roi.start <- pathog.mut.start - flanking.window
		roi.end <- pathog.mut.end + flanking.window
		pathog.mut.hgvs <- paste0(pathog.mut.start,"_",pathog.mut.end,"del")
		
		if( !is.null(snv.reference.allele) | !is.null(snv.alternative.allele) ){
			warning("snv.reference.allele and snv.alternative.allele are not used when pathog.mut.type is DEL")
		}
		if( !is.null(insertion.sequence) ){
			warning("insertion.sequence is not used when pathog.mut.type is SNV")
		}
	}
	
	if( pathog.mut.type=="INS" ){
		if( is.null(insertion.sequence) ){
			stop("insertion.sequence is required when pathog.mut.type is INS")
		}else{
			if( nchar(insertion.sequence) > sum(unlist(gregexpr("A|T|C|G",insertion.sequence))>0) ){
				stop("insert.sequence contains non DNA nucleotides")
			}else{
				pathog.mut.end <- pathog.mut.start + 1
				roi.start <- pathog.mut.start - flanking.window
				roi.end <- pathog.mut.end + flanking.window
				pathog.mut.hgvs <- paste0(pathog.mut.start,"ins",insertion.sequence)
			}
		}
		if( !is.null(snv.reference.allele) | !is.null(snv.alternative.allele) ){
			warning("snv.reference.allele and snv.alternative.allele are not used when pathog.mut.type is INS")
		}
		if( !is.null(deletion.sequence) | !is.null(deletion.length) ){
			warning("deletion.sequence and deletion.length are not used when pathog.mut.type is SNV")
		}
	}
	
	if( flanking.window <= 0 ){
		stop("flanking.window must be an integer greater than 0. Defalut is 100.")
	}
	
	message("Reversion detection for the pathogenic mutation ", chromosome, ":g.", pathog.mut.hgvs)
	
	#######################
	# MAIN
	#######################
	message("Extracting reads from the input bam file" )
	
	# Specify chromosome and region to extract reads from
	gr <- GRanges(seqnames = chromosome, ranges = IRanges(start = pathog.mut.start, end = pathog.mut.end))

	# Here we give it that region and tell it to only extract the position, CIGAR, and sequence fields
	params <- ScanBamParam(
		which = gr, 
		flag = scanBamFlag(isUnmappedQuery=FALSE, isDuplicate=FALSE),
		what = c("strand", "pos", "cigar", "seq") )
	aln <- scanBam(bamfile, param = params)
	names(aln) <- 'x'

	# message("Extracted ", length(aln$x$pos), " reads from bam file" )

	dp.total <- 0
	dp.wt.total <- 0
	dp.mt.total <- 0
	dp.rpm.total <- 0
	dp.gapN.total <- 0
	dp.rev.total <- 0
	dp.mt.rev <- 0
	dp.rpm.rev <- 0
	dp.wt.rev <- 0
	dp.rev.uniq <- 0
	
	rev.list.full <- c()
	rev.tbl.full <- c()

	if( length(aln$x$pos) > 0 ){
		
		message("Detecting reversions from ", length(aln$x$pos), " reads")
		progress.flag <- floor( quantile(1:length(aln$x$pos), seq(0.1,1,0.1)) )
		
		for (i in 1:length(aln$x$pos)) {

			# initialize flags		
			is.reads.wt <- FALSE
			is.reads.mt <- FALSE
			is.reads.rpm <- FALSE
			
			is.rpm.del <- FALSE
			is.rpm.ins <- FALSE
			is.rpm.snv <- FALSE
			
			is.inframe <- FALSE
			is.orf <- FALSE
			
			is.reversion <- FALSE
			is.new.reversion <- FALSE
			is.nonsense.reversion <- FALSE
			
			rev.type <- NULL
			rev.event <- NULL
	
			rev.tbl <- c()
	
			# Get cigar string
			reads <- as.character(aln$x$seq[i])
			aln.pos <- aln$x$pos[i]
			cigar <- aln$x$cigar[i]
			
			op.length <- as.numeric(unlist(strsplit(cigar, "M|I|D|N|S|H|P|=|X") ) )
			op.index <- unlist(gregexpr("M|I|D|N|S|H|P|=|X", cigar) )
			op.code <- substring(cigar, op.index, op.index)
			op <- data.frame( code = op.code, length = op.length, stringsAsFactors = FALSE )
			
			op <- op[ op$code != "H", ]
			op <- op[ op$code != "P", ]
			if( op$code[1]=="S" ){
				reads <- substring( reads, op$length[1]+1 )
			}
			if( op$code[nrow(op)]=="S" ){
				reads.rev <- intToUtf8( rev( utf8ToInt(reads) ) )
				reads.rev.clipped <- substring( reads.rev, op$length[nrow(op)]+1 )
				reads <- intToUtf8( rev( utf8ToInt(reads.rev.clipped) ) )
			}
			op <- op[ op$code != "S", ]
	
			ref.start <- aln.pos
			ref.length <- sum( op$length[ op$code %in% c("M","D","N","=","X") ] )
			ref.end <- ref.start + ref.length - 1
			
			if( ref.start > pathog.mut.start | ref.end < pathog.mut.end ){
				next
			}
				
			dp.total <- dp.total + 1
			
			### decipher CIGAR to map a single reads to genome
			aln.df <- data.frame( 
				cigar = rep(op$code, op$length),
				cigar_index = 1:sum(op$length),
				stringsAsFactors = FALSE )
			
			ref <- as.character( getSeq(reference.genome, chromosome, ref.start, ref.end) )
	
			ref.df <- aln.df[ aln.df$cigar %in% c("M","D","N","=","X"), ]
			ref.df$ref_index <- 1:nrow(ref.df)
			ref.df$ref_pos <- ref.start:ref.end
			ref.df$ref_seq <- unlist( strsplit(ref, "") )
			ref.df$cigar <- NULL
			
			reads.df <- aln.df[ aln.df$cigar %in% c("M","I","=","X"), ]
			reads.df$reads_index <- 1:nrow(reads.df)
			reads.df$reads_seq <- unlist( strsplit(reads, "") )
			reads.df$cigar <- NULL
			
			aln.df <- merge(aln.df, ref.df, by="cigar_index", all.x=TRUE)
			aln.df <- merge(aln.df, reads.df, by="cigar_index", all.x=TRUE)	
			
			aln.df$ref_seq[ aln.df$cigar=="I" ] <- "-"
			aln.df$reads_seq[ aln.df$cigar=="D" ] <- "-"
			aln.df$reads_seq[ aln.df$cigar=="N" ] <- "*"
			
			aln.df$snv <- ifelse( 
				aln.df$cigar%in%c("M","=","X") & aln.df$ref_seq!=aln.df$reads_seq, 
				paste(aln.df$ref_seq,aln.df$reads_seq, sep=">"),
				NA )
	
			del.event <- ifelse(aln.df$cigar=="D", aln.df$ref_seq, ",")
			del.event <- paste(del.event, collapse="")
			del.event <- unlist( strsplit(del.event, ",") )
			del.event <- del.event[ del.event != "" ]
			del.event.index <- which(aln.df$cigar=="D") 
			del.event.pos <- aln.df$ref_pos[setdiff(del.event.index, del.event.index+1)]
			ref.del.df <- data.frame( 
				ref_pos = del.event.pos,
				del = del.event,
				stringsAsFactors = FALSE )
			aln.df <- merge(aln.df, ref.del.df, by="ref_pos", all.x=TRUE)
			aln.df <- aln.df[ order(aln.df$cigar_index) ,]
			
			ins.event <- ifelse(aln.df$cigar=="I", aln.df$reads_seq, ",")
			ins.event <- paste(ins.event, collapse="")
			ins.event <- unlist( strsplit(ins.event, ",") )
			ins.event <- ins.event[ ins.event != "" ]
			ins.event.index <- which(aln.df$cigar=="I") 
			ins.event.pos <- aln.df$ref_pos[setdiff(ins.event.index-1, ins.event.index)]
			ref.ins.df <- data.frame( 
				ref_pos = ins.event.pos,
				ins = ins.event,
				stringsAsFactors = FALSE )
			aln.df <- merge(aln.df, ref.ins.df, by="ref_pos", all.x=TRUE)
			aln.df <- aln.df[ order(aln.df$cigar_index) ,]
			
			### extract alignment for ROI
			if( roi.start < ref.start ){
				roi.start.cigar.index <- aln.df$cigar_index[aln.df$ref_pos==ref.start & !is.na(aln.df$ref_pos) ]
			}else{
				roi.start.cigar.index <- aln.df$cigar_index[ aln.df$ref_pos==roi.start & !is.na(aln.df$ref_pos) ]
				if( aln.df$cigar[roi.start.cigar.index]=="D" & is.na(aln.df$del[roi.start.cigar.index]) ){
					roi.start.cigar.index <- max( aln.df$cigar_index[aln.df$cigar_index<roi.start.cigar.index & !is.na(aln.df$del)] )
				}
			}
			
			if( roi.end > ref.end ){
				roi.end.cigar.index <- aln.df$cigar_index[aln.df$ref_pos==ref.end & !is.na(aln.df$ref_pos) ]
			}else{
				roi.end.cigar.index <- aln.df$cigar_index[ aln.df$ref_pos==roi.end & !is.na(aln.df$ref_pos) ]
				if( !is.na( aln.df$ins[roi.end.cigar.index] ) ){
					roi.end.cigar.index <- roi.end.cigar.index + nchar(aln.df$ins[roi.end.cigar.index])
				}
			}
			
			roi.aln.df <- aln.df[ aln.df$cigar_index %in% c(roi.start.cigar.index:roi.end.cigar.index), ]
			
			if( sum(roi.aln.df$cigar=="N") > 0 ){
				dp.gapN.total <- dp.gapN.total + 1
				next
			}
			
			roi.ref.df <- roi.aln.df[ !is.na(roi.aln.df$ref_index), ]
			roi.ref.df <- roi.ref.df[ order(roi.ref.df$ref_index), ]
			roi.ref <- paste(roi.ref.df$ref_seq, collapse="")
			
			roi.reads.df <- roi.aln.df[ !is.na(roi.aln.df$reads_index), ]
			roi.reads.df <- roi.reads.df[ order(roi.reads.df$reads_index), ]
			roi.reads <- paste(roi.reads.df$reads_seq, collapse="")
			
			### identify reads type: wildtype, mutant, replacement-mutant
			if( pathog.mut.type == "SNV" ){
				
				true.genotype <- roi.ref.df$snv[roi.ref.df$ref_pos==pathog.mut.start]
				
				if( snv.genotype==true.genotype & !is.na(true.genotype) ){
					is.reads.mt <- TRUE
					dp.mt.total <- dp.mt.total + 1
				}else{
					rpm.snv <- true.genotype
					rpm.ins <- roi.ref.df$ins[ roi.ref.df$ref_pos %in% c(pathog.mut.start-1,pathog.mut.start) ] 
					rpm.del <- roi.ref.df$cigar[ roi.ref.df$ref_pos %in% c( pathog.mut.start:pathog.mut.end ) ]
					
					if( !is.na(rpm.snv) ){
						is.rpm.snv <- TRUE
					}
					if( sum(is.na(rpm.ins)) < length(rpm.ins) ){
						is.rpm.ins <- TRUE
					}
					if( sum(rpm.del=="D") > 0 ){
						is.rpm.del <- TRUE
					}
					
					if( is.rpm.snv | is.rpm.ins | is.rpm.del ){
						is.reads.rpm <- TRUE
						dp.rpm.total <- dp.rpm.total + 1
					}else{
						is.reads.wt <- TRUE
						dp.wt.total <- dp.wt.total + 1
					}
				}
			}
			
			if( pathog.mut.type == "DEL" ){
				
				true.deletion <- roi.ref.df$del[ roi.ref.df$ref_pos==pathog.mut.start ]
				
				if( deletion.sequence==true.deletion & !is.na(true.deletion) ){
					is.reads.mt <- TRUE
					dp.mt.total <- dp.mt.total + 1
				}else{
					rpm.snv <- roi.ref.df$snv[ roi.ref.df$ref_pos %in% c( pathog.mut.start:pathog.mut.end ) ]
					rpm.ins <- roi.ref.df$ins[ roi.ref.df$ref_pos %in% c( (pathog.mut.start-1):pathog.mut.end ) ] 
					rpm.del <- roi.ref.df$cigar[ roi.ref.df$ref_pos %in% c( pathog.mut.start:pathog.mut.end ) ]
					
					if( sum(is.na(rpm.snv)) < length(rpm.snv) ){
						is.rpm.snv <- TRUE
					}
					if( sum(is.na(rpm.ins)) < length(rpm.ins) ){
						is.rpm.ins <- TRUE
					}
					if( sum(rpm.del=="D") > 0 ){
						is.rpm.del <- TRUE
					}
					
					if( is.rpm.snv | is.rpm.ins | is.rpm.del ){
						is.reads.rpm <- TRUE
						dp.rpm.total <- dp.rpm.total + 1
					}else{
						is.reads.wt <- TRUE
						dp.wt.total <- dp.wt.total + 1
					}
				}
			}
			
			if( pathog.mut.type == "INS" ){
				
				true.insertion <- roi.ref.df$ins[ roi.ref.df$ref_pos==pathog.mut.start ] 
				
				if( insertion.sequence==true.insertion & !is.na(true.insertion) ){
					is.reads.mt <- TRUE
					dp.mt.total <- dp.mt.total + 1
				}else{
				    rpm.snv <- roi.ref.df$snv[ roi.ref.df$ref_pos %in% c( pathog.mut.start:pathog.mut.end ) ]
					rpm.ins <- true.insertion
					rpm.del <- roi.ref.df$cigar[ roi.ref.df$ref_pos %in% c( pathog.mut.start, pathog.mut.end ) ]
					
					if( sum(is.na(rpm.snv)) < length(rpm.snv) ){
					    is.rpm.snv <- TRUE
					}
					if( !is.na(rpm.ins) ){
						is.rpm.ins <- TRUE
					}
					if( sum(rpm.del=="D") > 0 ){
						is.rpm.del <- TRUE
					}
					
					if( is.rpm.snv | is.rpm.ins | is.rpm.del ){
						is.reads.rpm <- TRUE
						dp.rpm.total <- dp.rpm.total + 1
					}else{
						is.reads.wt <- TRUE
						dp.wt.total <- dp.wt.total + 1
					}
				}
			}
			
			
			### check nonsense SNV reversions
			if( is.reads.mt & pathog.mut.type=="SNV" ){
				
				nsRev.index <- c()
				nsRev.pos <- c()
				
				snv.cigar.index <- roi.aln.df$cigar_index[ roi.aln.df$ref_pos==pathog.mut.start & !is.na(roi.aln.df$ref_pos) ]
				snv5b.df <- roi.aln.df[ roi.aln.df$cigar_index %in% c( (snv.cigar.index-2):(snv.cigar.index+2) ) ,]
				snv5b.df <- snv5b.df[ order(snv5b.df$cigar_index), ]
				
				if( nrow(snv5b.df)>=3 & sum(!is.na(snv5b.df$snv)) > 1 ){
					
					snv5b.reads.wRev <- ifelse( snv5b.df$cigar=="D" | snv5b.df$cigar=="I", "N", snv5b.df$reads_seq )
					snv5b.reads.wRev <- paste( snv5b.reads.wRev, collapse="" )
					
					snv5b.reads.woRev <- ifelse(
						snv5b.df$cigar_index!=snv.cigar.index & !is.na(snv5b.df$snv), 
						snv5b.df$ref_seq, 
						snv5b.df$reads_seq )
					snv5b.reads.woRev <- ifelse( snv5b.df$cigar=="D" | snv5b.df$cigar=="I", "N", snv5b.reads.woRev)
					snv5b.reads.woRev <- paste( snv5b.reads.woRev, collapse="" )
					
					if( minus.strand ){
						snv5b.reads.wRev <- reverseComplement( DNAString(snv5b.reads.wRev) )
						snv5b.reads.woRev <- reverseComplement( DNAString(snv5b.reads.woRev) )
					}else{
						snv5b.reads.wRev <- DNAString(snv5b.reads.wRev) 
						snv5b.reads.woRev <- DNAString(snv5b.reads.woRev) 
					}
					
					snv5b.reads.wRev.aa <- sapply(
						1:(length(snv5b.reads.wRev)-2), 
						function(pos) as.character(translate(subseq(snv5b.reads.wRev, start=pos, width=3), if.fuzzy.codon="X") ) 
						)
					snv5b.reads.woRev.aa <- sapply(
						1:(length(snv5b.reads.woRev)-2), 
						function(pos) as.character(translate(subseq(snv5b.reads.woRev, start=pos, width=3), if.fuzzy.codon="X") ) 
						)
					
					snv5b.aa.df <- data.frame( 
						woRev=snv5b.reads.woRev.aa,
						wRev=snv5b.reads.wRev.aa, 
						stringsAsFactors = F )
					
					nsRev.index <- which( snv5b.aa.df$woRev=="*" & (snv5b.aa.df$wRev!="*" & snv5b.aa.df$wRev!="X") )
					
					if( length(nsRev.index) > 0 ){

						if(minus.strand){
							nsRev.pos <- lapply( 
								nsRev.index, 
								function(a){ 
									b <- snv5b.df[ (nrow(snv5b.df)-(a-1)):(nrow(snv5b.df)-(a-1)-2),]
									b$ref_pos[!is.na(b$snv) & b$cigar_index!=snv.cigar.index] 
									} 
								)
						}else{
							nsRev.pos <- lapply( 
								nsRev.index, 
								function(a){ 
									b <- snv5b.df[ a:(a+2), ] 
									b$ref_pos[!is.na(b$snv) & b$cigar_index!=snv.cigar.index] 
									} 
								)
						}
						
						nsRev.pos <- unique( unlist(nsRev.pos) )
						
						if( length(nsRev.pos) > 0 ){
							is.nonsense.reversion <- TRUE
						}
					}
				}
			}
			
			
			### check frame shift and reading frame (stop codon) for mutant reads
			# mt reads: complement indels
			if( is.reads.mt ){
				
				rev.ref.df <- roi.ref.df
				
				if( pathog.mut.type=="SNV" & is.nonsense.reversion ){
					rev.ref.df$snv[ ! rev.ref.df$ref_pos %in% nsRev.pos  ] <- NA
				}else{
					rev.ref.df$snv <- NA
				}
				
				if( pathog.mut.type=="DEL" ){
					rev.ref.df$del[ rev.ref.df$ref_pos==pathog.mut.start ] <- NA
				}
				
				if( pathog.mut.type=="INS" ){
					rev.ref.df$ins[ rev.ref.df$ref_pos==pathog.mut.start ] <- NA
				}
	
				rev.ref.df <- rev.ref.df[ !is.na(rev.ref.df$snv) | !is.na(rev.ref.df$ins) | !is.na(rev.ref.df$del) ,]
	
				if( nrow(rev.ref.df)>0 ){
					is.inframe <- abs( nchar(roi.reads)-nchar(roi.ref) ) %% 3 == 0
					is.orf <- suppressWarnings(checkORF(roi.reads, minus.strand = minus.strand))
				
					if( is.inframe & is.orf ){
						dp.mt.rev <- dp.mt.rev + 1
						dp.rev.total <- dp.rev.total + 1
						is.reversion <- TRUE
						rev.type <- "complement"
					}
				}
			}
			
			### check frame shift and reading frame (stop codon) for replacement-mutant reads
			# rpm reads: all replacement mutations + complement indels
			if( is.reads.rpm ){
				
				rev.ref.df <- roi.ref.df
				rev.ref.df$snv[ ! rev.ref.df$ref_pos %in% c(pathog.mut.start:pathog.mut.end) ] <- NA
				rev.ref.df <- rev.ref.df[ !is.na(rev.ref.df$snv) | !is.na(rev.ref.df$ins) | !is.na(rev.ref.df$del) ,]
				
				if( nrow(rev.ref.df)>0 ){
					is.inframe <- abs( nchar(roi.reads)-nchar(roi.ref) ) %% 3 == 0
					is.orf <- suppressWarnings(checkORF(roi.reads, minus.strand = minus.strand))
					
					if( is.inframe & is.orf ){
						dp.rpm.rev <- dp.rpm.rev + 1
						dp.rev.total <- dp.rev.total + 1
						is.reversion <- TRUE
						rev.type <- "replacement"
					}
				}
			}
			
			if( check.wildtype.reads ){
				if( is.reads.wt ){
					
					rev.ref.df <- roi.ref.df
					rev.ref.df$snv <- NA
					rev.ref.df <- rev.ref.df[ !is.na(rev.ref.df$ins) | !is.na(rev.ref.df$del) ,]
					
					if( nrow(rev.ref.df)>0 ){
						is.inframe <- abs( nchar(roi.reads)-nchar(roi.ref) ) %% 3 == 0
						is.orf <- suppressWarnings(checkORF(roi.reads, minus.strand = minus.strand))
						
						if( is.inframe & is.orf ){
							dp.wt.rev <- dp.wt.rev + 1
							dp.rev.total <- dp.rev.total + 1
							is.reversion <- TRUE
							rev.type <- "alternative"
						}
					}
				}
			}
			
			### output reversion event
			if( is.reversion ){
	
				rev.tbl <- data.frame(
					mut_start_pos = rep( rev.ref.df$ref_pos, each=3 ),
					mut_type = rep( c("SNV","DEL","INS"), nrow(rev.ref.df) ),
					mut_seq = as.vector( t(rev.ref.df[,c("snv","del","ins")]) ),
					stringsAsFactors = FALSE  )
				rev.tbl <- rev.tbl[ !is.na(rev.tbl$mut_seq), ]
				
				rev.tbl$mut_length <- nchar(rev.tbl$mut_seq)
				rev.tbl$mut_length[rev.tbl$mut_type=="SNV"] <- 0
				
				rev.tbl$mut_end_pos <- rev.tbl$mut_start_pos + rev.tbl$mut_length - 1
				rev.tbl$mut_end_pos[rev.tbl$mut_type=="SNV"] <- rev.tbl$mut_start_pos[rev.tbl$mut_type=="SNV"]
				rev.tbl$mut_end_pos[rev.tbl$mut_type=="INS"] <- rev.tbl$mut_start_pos[rev.tbl$mut_type=="INS"] + 1
				
				rev.tbl$mut_hgvs <- rep(NA, nrow(rev.tbl))
				
				rev.tbl$mut_hgvs[rev.tbl$mut_type=="DEL"] = paste0( 
					rev.tbl$mut_start_pos[rev.tbl$mut_type=="DEL"], 
					"_", 
					rev.tbl$mut_end_pos[rev.tbl$mut_type=="DEL"], 
					"del" )
				
				rev.tbl$mut_hgvs[rev.tbl$mut_type=="INS"] = paste0( 
					rev.tbl$mut_start_pos[rev.tbl$mut_type=="INS"], 
					"ins",
					rev.tbl$mut_seq[rev.tbl$mut_type=="INS"] )
				
				rev.tbl$mut_hgvs[rev.tbl$mut_type=="SNV"] = paste0( 
					rev.tbl$mut_start_pos[rev.tbl$mut_type=="SNV"], 
					rev.tbl$mut_seq[rev.tbl$mut_type=="SNV"] )
				
				rev.event <- paste( 
					rev.type,
					paste(rev.tbl$mut_hgvs, collapse=":"),
					sep="|" )
				
				if( is.null(rev.list.full) ){
					is.new.reversion <- TRUE
					dp.rev.uniq <- 1
					rev.list.full <- data.frame(
						rev_id = paste0( "rev", dp.rev.uniq ),
						rev_freq = 1,
						rev_event = rev.event,
						stringsAsFactors = FALSE )
				}else{
					if( rev.event %in% rev.list.full$rev_event ){
						rev.list.full$rev_freq[ rev.list.full$rev_event==rev.event ] <- rev.list.full$rev_freq[ rev.list.full$rev_event==rev.event ] + 1  
					}else{
						is.new.reversion <- TRUE
						dp.rev.uniq <- dp.rev.uniq + 1
						rev.list.full[ nrow(rev.list.full)+1 , ] <- list( paste0("rev",dp.rev.uniq), 1, rev.event)
					}
				}
				
				if( is.new.reversion ){
					
					rev.tbl$rev_id <- paste0( "rev", rep(dp.rev.uniq, nrow(rev.tbl) ) )
					rev.tbl$rev_type <- rep( rev.type, nrow(rev.tbl) )
					rev.tbl$rev_mut_number <- 1:nrow(rev.tbl)
					rev.tbl$mut_id <- paste(rev.tbl$rev_id, rev.tbl$rev_mut_number, sep=".")
					rev.tbl$chr <- rep( chromosome, nrow(rev.tbl) )
					
					rev.tbl$pathog_mut_hgvs <- rep(pathog.mut.hgvs, nrow(rev.tbl))
					rev.tbl$pathog_mut_start <- rep(pathog.mut.start, nrow(rev.tbl))
					rev.tbl$pathog_mut_end <- rep(pathog.mut.end, nrow(rev.tbl))
					
					rev.tbl$dist_to_pathog_mut <- mapply(
						getMutationsDistance,
						rev.tbl$pathog_mut_start,
						rev.tbl$pathog_mut_end,
						rep(pathog.mut.type, nrow(rev.tbl)),
						rev.tbl$mut_start_pos,
						rev.tbl$mut_end_pos,
						rev.tbl$mut_type,
						SIMPLIFY = TRUE
					)
					
					rev.tbl.full <- rbind(rev.tbl.full, rev.tbl)
				} # end of is.new.reversion
			
			} # end of is.reversion
		
			if( i %in% progress.flag ){
				message( "===> ", names(progress.flag)[progress.flag==i], " reads have been processed" )
			}
			
		} # end of for loop
	
		message("Reversion detection completed")
	} # end of length(aln$x$pos)>0
	
	# cat("Reads aligned to pathogenic mutation:", dp.total, "\n")
	# cat("Reads with wild-type sequence:", dp.wt.total, "\n")
	# cat("Reads with pathogenic mutation:", dp.mt.total, "\n")
	# cat("Reads with replacement mutations:", dp.rpm.total, "\n")
	# cat("Reads with reversions:", dp.rev.total, "\n")
	# cat("Reads with complement reversions:", dp.mt.rev, "\n")
	# cat("Reads with replacement reversions:", dp.rpm.rev, "\n")
	# if( check.wildtype.reads ){
	# 	cat("Reads with alternative reversions:", dp.wt.rev, "\n")
	# }
	# cat("Unique reversions:", dp.rev.uniq, "\n")
	
	stats.df <- data.frame(
		Category = c(
			"reads_aligned_to_original_mutation", 
			"reads_wildtype", 
			"reads_with_original_mutation", 
			"reads_with_replacement_mutations", 
			"reads_with_reversions", 
			"reads_with_complement_reversions", 
			"reads_with_replacement_reversions", 
			"reads_with_alternative_reversions",
			"unique_reversions" ),
		Counts = c(dp.total, dp.wt.total, dp.mt.total, dp.rpm.total, dp.rev.total, dp.mt.rev, dp.rpm.rev, dp.wt.rev, dp.rev.uniq),
		stringsAsFactors = FALSE )
	
	if( !check.wildtype.reads ){
		stats.df <- stats.df[ stats.df$Category!="reads_with_alternative_reversions", ]
	}
	
	output.columns <- c(
		"rev_id",
		"rev_freq",
		"rev_type",
		"rev_mut_number",
		"mut_id",
		"chr", 
		"mut_start_pos",
		"mut_type",
		"mut_seq",
		"mut_length",
		"mut_hgvs",
		"pathog_mut_hgvs",
		"dist_to_pathog_mut" )
	
	if( is.null(rev.list.full) ){
		rev.tbl.full <- data.frame( matrix(ncol=length(output.columns), nrow=0) ) 
		colnames(rev.tbl.full) <- output.columns
	}else{
		rev.tbl.full$rev_freq <- rev.list.full$rev_freq[ match(rev.tbl.full$rev_id, rev.list.full$rev_id) ]
		rev.tbl.full <- rev.tbl.full[, output.columns]
	}
	
	reversions <- list(stats = stats.df,
					   rev = rev.tbl.full )
	
	return(reversions)
	
} 
		