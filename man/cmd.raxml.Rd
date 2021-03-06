% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/bigphylo.cmd.R
\name{cmd.raxml}
\alias{cmd.raxml}
\title{Produce a single RAxML shell command.}
\usage{
cmd.raxml(infile.fasta, outfile = paste(infile.fasta, ".newick", sep =
  ""), pr = PR.RAXML, pr.args = "-m GTRCAT --HKY85 -p 42")
}
\value{
Character string
}
\description{
Produce a single RAxML shell command.
}
