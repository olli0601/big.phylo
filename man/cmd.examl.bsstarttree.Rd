% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/bigphylo.cmd.R
\name{cmd.examl.bsstarttree}
\alias{cmd.examl.bsstarttree}
\title{Produce a shell command to compute the final boostrap tree.}
\usage{
cmd.examl.bsstarttree(indir, infile, bs.from = 0, bs.to = 99,
  bs.n = bs.to - bs.from + ifelse(bs.from == 0, 1, 0), outdir = indir,
  prog.mpi = "mpiexec", prog.parser = PR.EXAML.PARSER,
  prog.starttree = PR.EXAML.STARTTREE, prog.examl = PR.EXAML.EXAML,
  args.examl = "-m GAMMA -D", prog.supportadder = PR.EXAML.BS,
  resume = 1, verbose = 1)
}
\value{
Character string.
}
\description{
Internal code.
}
