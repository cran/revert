#' Translate a DNA sequence and determine if it results in one or more stop codons.
#' @description
#' 	Translate a DNA sequence and determine if it results in one or more stop codons.
#' @param x
#' 	DNAStringSet object containing the sequence to test.
#' @param minus.strand
#'	Logical. TRUE if the gene in question is on the reverse strand. Default is FALSE.
#' @return 
#'  Returns TRUE if at least one RF translation is free of stop codons.
#' @importFrom Biostrings DNAString
#' @importFrom Biostrings reverseComplement
#' @importFrom Biostrings translate
#' @importFrom Biostrings subseq
#' @examples
#' \dontrun{
#'     checkORF()
#' }
#' @keywords internal
#' @noRd 
checkORF <- function(x, minus.strand = FALSE) {
	# Translate a DNA sequence and determine if it results in one or more stop codons
	# Returns TRUE if at least one RF translation is free of stop codons
	
	if (grepl('N', x)) {
		# Can't determine this if there is an ambiguous base, so just return false
		return(FALSE)
	}else{
		x <- DNAString(x)
	}
	
	if (minus.strand) {
		# Switch seq to its reverse complement to translate it on the reverse strand
		x <- reverseComplement(x)
	}
	
	# Suppress warnings because it will warn us when it skips the last 1-2 bases
	# that don't make up a whole codon
	# Generate all three reading frames by sub-setting from different start positions
	reading.frames <- lapply(1:3, function(pos) subseq(x, start=pos))
	# Translate each RF, see ?translate
	aa.seq <- lapply(reading.frames, translate)
	
	has.ORF <- FALSE
	
	# Look for stop codon (*) in the AA sequence
	orf <- sapply( aa.seq, function(aa){!grepl('\\*',aa)} )
	
	if( sum(orf) > 0 ){
		has.ORF <- TRUE
	}
	
	# Return TRUE if we found one that has no stop codons
	return(has.ORF)
}



#' Get distance between the original mutation and a reversion mutation
#' @description
#' 	Get distance between the original mutation and a reversion mutation.
#' @param original.mut.start
#'  Integer. Genomic coordinate of the start position of original mutation.
#' @param original.mut.end
#'  Integer. Genomic coordinate of the end position of original mutation.
#' @param original.mut.type
#'  Character. Type of original mutation: SNV, DEL or INS.
#' @param secondary.mut.start
#'  Integer. Genomic coordinate of the start position of secondary mutation.
#' @param secondary.mut.end
#'  Integer. Genomic coordinate of the end position of secondary mutation.
#' @param secondary.mut.type
#'	Character. Type of secondary mutation: SNV, DEL or INS.
#' @return 
#'  Returns the distance in basepairs (bp) between two input mutations.
#' @examples
#' \dontrun{
#'     getMutationsDistance()
#' }
#' @keywords internal
#' @noRd
getMutationsDistance <- function(
		original.mut.start, 
		original.mut.end,
		original.mut.type, 
		secondary.mut.start, 
		secondary.mut.end,
		secondary.mut.type ){
	
	if( original.mut.type=="INS" & secondary.mut.type=="INS" ){
		original.mut.end <- original.mut.start
		secondary.mut.end <- secondary.mut.start
	}
	
	if( original.mut.type!="INS" & secondary.mut.type=="INS" ){
		original.mut.start <- original.mut.start - 1 
		secondary.mut.end <- secondary.mut.start
	}
	
	ovlp <- intersect( c(original.mut.start:original.mut.end), c(secondary.mut.start:secondary.mut.end) )
	
	if( length(ovlp)>0 ){
		dist.sec2orig <- 0
	}else{
		d1 <- secondary.mut.start - original.mut.end
		d2 <- secondary.mut.end - original.mut.start
		dist.sec2orig <- sign(d1) * min( abs(d1), abs(d2) )
	}
	
	return(dist.sec2orig)
}





