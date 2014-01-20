# This file contains parallel processing scripts that generate and submit PBS scripts to a high performance system. 
######################################################################################

#' @export 
#' @title Compute phylogeny with bootstrap values 
#' @description Compute a maximum likelihood phylogeny with the \code{ExaML} program. Bootstrap values
#' 	capture the uncertainty in the tree phylogeny and are computed as follows. Columns of the multiple
#' 	sequence alignment are bootstrapped to generate a random sequence data set. The starting tree is randomly
#'  constructed. A single phylogeny is constructed, and referred to as a boostrap phylogeny. Based on a set 
#'  of bootstrap phylogenies, the one with best likelihood is determined and designated as the final ExaML
#'  phylogeny. Node uncertainty is computed by calculating how often any particular subtree in the final ExaML
#'  phylogeny occurs in the set of boostrap phylogenies.
#' @param indir			Directory of the sequences in file \code{infile}
#' @param infile		File name to an DNAbin matrix. The file name must end in \code{.R}
#' @param outdir		Output directory for the bootstrap calculations
#' @param bs.from		Iteration number of bootstrap calculations. Defaults to \code{0}.
#' @param bs.n			Total number of boostrap iterations. Defaults to \code{500}.
#' @param bs.to			Final iteration number of this call to \code{pipeline.ExaML.bootstrap.per.proc}. Defaults to \code{bs.n}.
#' @param hpc.walltime	Walltime specification for PBS header.
#' @param hpc.mem		RAM specification for PBS header.
#' @param hpc.nproc		Number of processors specification for PBS header.
#' @param hpc.q			Queue specification for PBS header.
#' @param verbose		Flag to run function in verbose mode.
#' @return NULL. Creates shell files in \code{outdir}, and attempts to submit those to a queuing system.
#' @example example/ex.h3n2.ExaML.bootstrap.per.proc.R
pipeline.ExaML.bootstrap.per.proc<- function(indir, infile, outdir=indir, bs.from=0, bs.n= 500, bs.to= bs.n, hpc.walltime=1, hpc.q=NA, hpc.mem="500mb", hpc.nproc=1, verbose=1)
{	
	#	sense check
	if(!grepl('.R',infile))							stop("expect R infile that ends in .R")		
	file		<- paste(indir,'/',infile,sep='')
	tmp			<- load(file)
	if(verbose)	cat(paste('\nloaded file=', tmp))
	if(tmp!='seq')
		eval(parse(text=paste("seq<- ",tmp,sep='')))
	if(!"DNAbin"%in%class(seq) || !is.matrix(seq))	stop("expect R infile that contains a DNAbin matrix")
	#	generate shell commands to run a construct a single ExaML bootstrap phylogeny
	infile		<- substr(infile, 1, nchar(infile)-2)
	cmd			<- cmd.examl.bootstrap(indir, infile, bs.from=bs.from, bs.to=bs.to, bs.n=bs.n, outdir=outdir, resume=1, verbose=verbose)
	#	add HPC wrapper and submit job
	dummy		<- lapply(cmd, function(x)
			{				
				x		<- cmd.hpcwrapper(x, hpc.walltime=24, hpc.q=hpc.q, hpc.mem=hpc.mem, hpc.nproc=hpc.nproc)
				signat	<- paste(strsplit(date(),split=' ')[[1]],collapse='_',sep='')
				outfile	<- paste("exa",signat,sep='.')
				#cat(x)
				cmd.hpccaller(outdir, outfile, x)
				Sys.sleep(1)
			})	
}



