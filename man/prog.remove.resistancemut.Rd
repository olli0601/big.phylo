\name{prog.remove.resistancemut}
\alias{prog.remove.resistancemut}
\title{Remove resistance mutations}
\usage{
  prog.remove.resistancemut()
}
\value{
  NULL. A sequence alignment with resistance codons set to
  'NNN' is written to file.
}
\description{
  Input parameters are 'indir', 'infile', 'outdir',
  'outfile', 'alignment.start', 'resume' and 'verbose', and
  specified via an \code{argv} string, see the Examples.
  The 'alignment.start' option specifies the position of
  the first nucleotide relative to the HXB2 reference
  sequence.
}
\examples{
require(big.phylo)
#	specify file that holds an R DNAbin object
infile			<- 'myfile.R'
indir			<- 'mydir'
save(seq, file=paste(indir, infile, sep='/'))	
#	create the command string
infile			<- substr(infile, 1, nchar(infile)-2)
outfile			<- paste(infile,'out',sep='_')
alignment.start	<- 2253
argv			<<- cmd.rm.resistance(indir, infile, outfile, alignment.start=alignment.start) 
argv			<<- unlist(strsplit(argv,' '))
#	create the bootstrap alignment
#prog.remove.resistancemut()
}

