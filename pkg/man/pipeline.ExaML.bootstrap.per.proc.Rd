\name{pipeline.ExaML.bootstrap.per.proc}
\alias{pipeline.ExaML.bootstrap.per.proc}
\title{Compute phylogeny with bootstrap values}
\usage{
  pipeline.ExaML.bootstrap.per.proc(indir, infile,
    outdir = indir, bs.from = 0, bs.n = 500, bs.to = bs.n,
    hpc.walltime = 1, hpc.q = NA, hpc.sys = cmd.hpcsys(),
    hpc.mem = "500mb", hpc.nproc = 1, verbose = 1)
}
\arguments{
  \item{indir}{Directory of the sequences in file
  \code{infile}}

  \item{infile}{File name to an DNAbin matrix. The file
  name must end in \code{.R}}

  \item{outdir}{Output directory for the bootstrap
  calculations}

  \item{bs.from}{Iteration number of bootstrap
  calculations. Defaults to \code{0}.}

  \item{bs.n}{Total number of boostrap iterations. Defaults
  to \code{500}.}

  \item{bs.to}{Final iteration number of this call to
  \code{pipeline.ExaML.bootstrap.per.proc}. Defaults to
  \code{bs.n}.}

  \item{hpc.sys}{Server specification for PBS header.}

  \item{hpc.walltime}{Walltime specification for PBS
  header.}

  \item{hpc.mem}{RAM specification for PBS header.}

  \item{hpc.nproc}{Number of processors specification for
  PBS header.}

  \item{hpc.q}{Queue specification for PBS header.}

  \item{verbose}{Flag to run function in verbose mode.}
}
\value{
  NULL. Creates shell files in \code{outdir}, and attempts
  to submit those to a queuing system.
}
\description{
  Compute a maximum likelihood phylogeny with the
  \code{ExaML} program. Bootstrap values capture the
  uncertainty in the tree phylogeny and are computed as
  follows. Columns of the multiple sequence alignment are
  bootstrapped to generate a random sequence data set. The
  starting tree is randomly constructed. A single phylogeny
  is constructed, and referred to as a boostrap phylogeny.
  Based on a set of bootstrap phylogenies, the one with
  best likelihood is determined and designated as the final
  ExaML phylogeny. Node uncertainty is computed by
  calculating how often any particular subtree in the final
  ExaML phylogeny occurs in the set of boostrap
  phylogenies.
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
#	Write the a shell file to create the first ExaML boostrap tree
pipeline.ExaML.bootstrap.per.proc(indir, infile, bs.from=0, bs.to=0, hpc.walltime=1, hpc.q=NA, hpc.mem="500mb", hpc.nproc=1)
#	Resume and write the shell files to create two more ExaML boostrap trees
pipeline.ExaML.bootstrap.per.proc(indir, infile, bs.from=1, bs.to=2, hpc.walltime=1, hpc.q=NA, hpc.mem="500mb", hpc.nproc=1)
#
#	NOTE: if the current Unix system does not have a 'qsub' command, then a shell file will be generated that can
#	be run manually on the current Unix system. Support for an HPC system is currently implemented for the Imperial
#	server 'cx1.hpc.ic.ac.uk'. To see how the final HPC script would look like, try
#
#	pipeline.ExaML.bootstrap.per.proc(indir, infile, bs.from=1, bs.to=2, hpc.sys='cx1.hpc.ic.ac.uk', hpc.walltime=1, hpc.q=NA, hpc.mem="500mb", hpc.nproc=1)
#
#	It is not difficult to enable support for a different server.
#	See ?cmd.hpcwrapper
}

