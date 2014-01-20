# This file contains parallel processing scripts that generate and submit PBS scripts to a high performance system. 
######################################################################################

#' @title Identify candidate recombinant sequences and corresponding parent sequences in parallel
#' @description A large set of sequences in file \code{infile} in directory \code{indir} is scanned 
#' for recombinants and associated parent sequences with \code{3SEQ}. The sequences are processed in \code{batch.n}
#' separate jobs. For each job, this function creates a separate shell file that can be run in a UNIX environment.
#' If an HPC system is detected, the separate shell files are submitted to the queuing system.
#' @param indir		Directory of the sequences in file \code{infile}
#' @param infile 	File name to an DNAbin matrix. The file name must end in \code{.R}
#' @param batch.n	Number of jobs to run in parallel.
#' @param hpc.walltime	Walltime specification for PBS header.
#' @param hpc.q			Queue specification for PBS header.
#' @param hpc.mem		RAM specification for PBS header.
#' @param hpc.nproc		Number of processors specification for PBS header.
#' @param verbose		Flag to run function in verbose mode.
#' @return NULL. Creates shell files, and attempts to submit those to a queuing system
#' @example example/package.neisseria.run.3seq.R
#' @example example/package.mtDNA.run.3seq.R
#' @export
pipeline.recom.run.3seq<- function(indir, infile, batch.n=100, hpc.sys= cmd.hpcsys(), hpc.walltime=35, hpc.q=NA, hpc.mem="3850mb", hpc.nproc=1, verbose=1)
{
	#load sequences into 'seq' object 
	if(!grepl('.R',infile))							stop("expect R infile that ends in .R")		
	file			<- paste(indir,'/',infile,sep='')
	tmp				<- load(file)
	if(verbose)	cat(paste('\nloaded file=', tmp))
	if(tmp!='seq')
		eval(parse(text=paste("seq<- ",tmp,sep='')))
	if(!"DNAbin"%in%class(seq) || !is.matrix(seq))	stop("expect R infile that contains a DNAbin matrix")
	#write phylip file
	infile			<- substr(infile, 1, nchar(infile)-2)	
	file			<- paste(indir,'/',infile,".phylip",sep='')
	if(verbose)	cat(paste('\nwrite file=', file))
	seq.write.dna.phylip(seq, file)
	#generate parallel calls to 3seq
	#	batch.n	for 1e4 sequences, does about 100 in 25hrs, so request 100 batches for walltime 25 expected + 10hrs grace
	batch.seq		<- round(seq.int(0,nrow(seq),len=batch.n),d=0)
	if(tail(batch.seq,1)!=nrow(seq))
		batch.seq	<- c(batch.seq, nrow(seq))
	batch.seq		<- rbind(batch.seq[-length(batch.seq)], batch.seq[-1]-1)
	#batch.seq		<- batch.seq[,1:10]	#test run
	lapply(seq_len(ncol(batch.seq)),function(j)
			{					
				cmd			<- cmd.recombination.run.3seq(infile=file, outfile=paste(indir,'/',infile,'_',batch.seq[1,j],'-',batch.seq[2,j],".3seq",sep=''), recomb.3seq.siglevel=0.1, nproc=1, recomb.3seq.testvsall.beginatseq=batch.seq[1,j], recomb.3seq.testvsall.endatseq=batch.seq[2,j], verbose=1)				
				cmd			<- cmd.hpcwrapper(cmd, hpc.sys=hpc.sys, hpc.walltime=hpc.walltime, hpc.q=hpc.q, hpc.mem=hpc.mem,  hpc.nproc=hpc.nproc)
				if(verbose)	cat(cmd)
				outdir		<- indir
				outfile		<- paste("r3seq",paste(strsplit(date(),split=' ')[[1]],collapse='_',sep=''),sep='.')									
				cmd.hpccaller(outdir, outfile, cmd)			
			})	
	NULL
}
######################################################################################
#' @export
pipeline.recom.get.phyloincongruence.for.candidates<- function(indir, infile, insignat, resume=0, verbose=1, hpc.walltime=35, hpc.q=NA, hpc.mem="600mb", hpc.nproc=1)
{	
	argv				<<-	cmd.recombination.process.3SEQ.output(indir, infile, insignat, resume=1, verbose=1) 
	argv				<<- unlist(strsplit(argv,' '))
	df.recomb			<- prog.recom.process.3SEQ.output()	
	
	triplets			<- seq_len(nrow(df.recomb))
	#triplets			<- 147:nrow(df.recomb)
	dummy	<- lapply(triplets, function(i)
			{					
				if(verbose)	cat(paste("\nprocess triplet number",i,"\n"))
				argv				<<- cmd.recombination.check.candidates(indir, infile, insignat, i, resume=resume, verbose=1, hpc.walltime=hpc.walltime, hpc.q=hpc.q, hpc.mem=hpc.mem, hpc.nproc=hpc.nproc)
				argv				<<- unlist(strsplit(argv,' '))
				prog.recom.get.incongruence()		#this starts ExaML for the ith triplet			
			})	
}
######################################################################################
#' @export
pipeline.recom.plot.phyloincongruence.for.candidates<- function(indir, infile, insignat, resume=0, verbose=1)
{	
	argv				<<- cmd.recombination.plot.incongruence(indir, infile, insignat, triplet.id=NA, opt.select="ng2", verbose=1)
	argv				<<- unlist(strsplit(argv,' '))
	prog.recom.plot.incongruence()		
	
	argv				<<- cmd.recombination.plot.incongruence(indir, infile, insignat, triplet.id=NA, opt.select="g2", verbose=1)
	argv				<<- unlist(strsplit(argv,' '))
	prog.recom.plot.incongruence()	
}



