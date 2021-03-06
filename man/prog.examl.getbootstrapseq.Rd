% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/bigphylo.prog.R
\name{prog.examl.getbootstrapseq}
\alias{prog.examl.getbootstrapseq}
\title{Program to generate a booststrap alignment}
\usage{
prog.examl.getbootstrapseq()
}
\value{
NULL. A boostrap alignment is written to file in phylip format.
}
\description{
Input parameters are 'indir', 'infile', 'outdir', 'bootstrap', 'by', 'resume' and 'verbose', and specified via an 
\code{argv} string, see the Examples. The 'bootstrap' option specifies the boostrap iteration number, e. g. '-bootstrap=0' for the
first bootstrap iteration. The 'by' option specifies the way the boostrap alignment is created. Valid options are \code{codon} and \code{nucleotide}.
}
\examples{
require(big.phylo)
data(nz_h3n2)
#	H3N2 data is stored in 'seq' DNAbin matrix object
print( seq )							
#	write sequences to directory for processing
indir		<- getwd()
infile		<- 'nz_h3n2.R'
insignat	<- ''
save(seq, file=paste(indir, infile, sep='/'))	
#	create the command string
bs.id		<- 1
infile		<- substr(infile, 1, nchar(infile)-2)
argv		<<- cmd.examl.bsalignment(indir, infile, bs.id) 
argv		<<- unlist(strsplit(argv,' '))
#	create the bootstrap alignment
prog.examl.getbootstrapseq()
}
\seealso{
\code{\link{pipeline.ExaML.bootstrap.per.proc}}
}
