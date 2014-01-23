\name{cmd.examl.bsalignment}
\alias{cmd.examl.bsalignment}
\title{Produce the shell command that creates the bootstrap alignment.}
\usage{
  cmd.examl.bsalignment(indir, infile, bs.id,
    outdir = indir, prog.bscreate = PR.EXAML.BSCREATE,
    opt.bootstrap.by = "codon", resume = 0, verbose = 1)
}
\arguments{
  \item{opt.bootstrap.by}{Character string that specifies
  specifies the way the boostrap alignment is created.
  Valid options are \code{codon} and \code{nucleotide}.}

  \item{indir}{Directory of the sequences in file
  \code{infile}}

  \item{infile}{File name to an DNAbin matrix. The file
  name must end in \code{.R}}

  \item{outdir}{Output directory for the bootstrap
  calculations}

  \item{verbose}{Flag to run function in verbose mode.}
}
\value{
  Character string
}
\description{
  Internal code.
}

