
#' @export
#' @import ape
#' @title Program to generate a booststrap alignment
#' @description Input parameters are 'indir', 'infile', 'outdir', 'bootstrap', 'by', 'resume' and 'verbose', and specified via an 
#' \code{argv} string, see the Examples. The 'bootstrap' option specifies the boostrap iteration number, e. g. '-bootstrap=0' for the
#' first bootstrap iteration. The 'by' option specifies the way the boostrap alignment is created. Valid options are \code{codon} and \code{nucleotide}.
#' @return NULL. A boostrap alignment is written to file in phylip format.
#' @example example/ex.h3n2.ExaML.getbootstrapseq.R
#' @seealso \code{\link{pipeline.ExaML.bootstrap.per.proc}}
#' 
prog.examl.getbootstrapseq<- function()
{		
	require(big.phylo)		#need this for command line execution
	indir				<- outdir		<- infile	<- ''	
	verbose				<- resume		<- 1
	opt.bootstrap.by	<- "codon"
	bs					<- 0
	check.any.bs.identical	<- 0
	
	if(exists("argv"))
	{
		tmp<- na.omit(sapply(argv,function(arg)
						{	switch(substr(arg,2,6),
									indir= return(substr(arg,8,nchar(arg))),NA)	}))
		if(length(tmp)>0) indir<- tmp[1]
		tmp<- na.omit(sapply(argv,function(arg)
						{	switch(substr(arg,2,7),
									outdir= return(substr(arg,9,nchar(arg))),NA)	}))
		if(length(tmp)>0) outdir<- tmp[1]
		tmp<- na.omit(sapply(argv,function(arg)
						{	switch(substr(arg,2,7),
									infile= return(substr(arg,9,nchar(arg))),NA)	}))
		if(length(tmp)>0) infile<- tmp[1]				
		tmp<- na.omit(sapply(argv,function(arg)
						{	switch(substr(arg,2,7),
									resume= return(as.numeric(substr(arg,9,nchar(arg)))),NA)	}))
		if(length(tmp)>0) resume<- tmp[1]
		tmp<- na.omit(sapply(argv,function(arg)
						{	switch(substr(arg,2,2),
									v= return(as.numeric(substr(arg,4,nchar(arg)))),NA)	}))
		if(length(tmp)>0) verbose<- tmp[1]
		tmp<- na.omit(sapply(argv,function(arg)
						{	switch(substr(arg,2,10),
									bootstrap= return(as.numeric(substr(arg,12,nchar(arg)))),NA)	}))
		if(length(tmp)>0) bs<- tmp[1]
		tmp<- na.omit(sapply(argv,function(arg)
						{	switch(substr(arg,2,3),
									by= return(substr(arg,5,nchar(arg))),NA)	}))
		if(length(tmp)>0) opt.bootstrap.by<- tmp[1]
	}
	if(verbose)
	{
		print( indir ) 
		print(outdir)
		print(infile)
		print(verbose)
		print(resume)
		print(bs)
		print(opt.bootstrap.by)		
	}
	if(!opt.bootstrap.by%in%c("nucleotide","codon"))	stop("Unexpected opt.bootstrap.by")		
	pattern 	<- paste(infile,".phylip.",sprintf("%03d",bs),sep='')
	file		<- list.files(path=outdir, pattern=pattern, full.names=1)
	if(!resume || !length(file))	
	{					
		file		<- paste(outdir,"/",infile,".R",sep='')
		if(verbose) cat(paste("\nread",file))
		tmp			<- load(file)
		if(length(tmp)!=1)		
			stop("Unexpected lenght of loaded objects")
		eval(parse(text=paste("seq.PROT.RT<- ",tmp,sep='')))
		if(!"DNAbin"%in%class(seq.PROT.RT) || !is.matrix(seq.PROT.RT))	
			stop("expect R infile that contains a DNAbin matrix")
		print(seq.PROT.RT)
		#print(bs)
		if(bs)		#keep bs0 intact
		{
			dummy			<- 0
			any.eq			<- 1
			j				<- 0
			while(any.eq)
			{
				j			<- j+1
				if(opt.bootstrap.by=="codon")
				{
					bs.blocks.n	<- floor( ncol(seq.PROT.RT )/3)
					bs.blocks.s	<- sample(seq_len(bs.blocks.n),bs.blocks.n,replace=T)-1
					bs.seq.s	<- as.numeric( sapply(bs.blocks.s,function(x)		3*x+c(1,2,3)		) )
				}
				else if(opt.bootstrap.by=="nucleotide")
				{
					bs.seq.s	<- sample(seq_len(ncol(seq.PROT.RT )),ncol(seq.PROT.RT ),replace=T)
				}
				seq.BS		<- seq.PROT.RT[,bs.seq.s]				
				if(check.any.bs.identical)
				{
					if(verbose) cat(paste("\ncheck for identity proposed boostrap seq alignment no",j))
					#check no seqs identical								
					for(i1 in seq_len(nrow(seq.BS)-1))
					{
						seq1		<- seq.BS[i1,]
						tmp			<- 1-sapply(seq.int(i1+1,nrow(seq.BS)),function(i2)
													{		
														.C("hivc_dist_ambiguous_dna", seq1, seq.BS[i2,], ncol(seq1), dummy )[[4]]			
													})
						#print(tmp)
						if(any(tmp==0))
						{
							print(tmp)
							break
						}									
						if(i1==nrow(seq.BS)-1)
							any.eq	<- 0
					}
					if(verbose) cat(paste("\nchecked for identity proposed boostrap seq alignment no",j,"is any identical",any.eq))
				}
				else
					any.eq	<- 0
			}					
		}
		else
		{
			cat(paste("\nkeep boostrap seq alignment no",bs,"as original"))
			seq.BS	<- seq.PROT.RT
		}
		file		<- paste(outdir,"/",infile,".phylip.",sprintf("%03d",bs),sep='')
		cat(paste("\nsave boostrap seq alignment to",file))
		seq.write.dna.phylip(seq.BS, file=file)
	}
	else
		cat("\nfound boostrap sequence alignment")
}

#' @export 
#' @title Remove resistance mutations
#' @description Input parameters are 'indir', 'infile', 'outdir', 'outfile', 'alignment.start', 'resume' and 'verbose', and specified via an 
#' \code{argv} string, see the Examples. The 'alignment.start' option specifies the position of the first nucleotide relative to the HXB2
#' reference sequence.
#' @return NULL. A sequence alignment with resistance codons set to 'NNN' is written to file.
#' @example example/ex.remove.resistancemut.R
#' 
prog.remove.resistancemut<- function()
{
	library(big.phylo)
	#load data.table 'dr' with drug resistance mutations
	data(IAS.primarydrugresistance.201303)
	
	indir			<- paste(DATA,"tmp",sep='/')
	infile			<- "ATHENA_2013_03_CurAll+LANL_Sequences"	
	outdir			<- paste(DATA,"tmp",sep='/')
	outfile			<- "ATHENA_2013_03_NoDRAll+LANL_Sequences"	
	alignment.start	<- 2253	
	verbose			<- 1
	if(exists("argv"))
	{
		tmp<- na.omit(sapply(argv,function(arg)
						{	switch(substr(arg,2,6),
									indir= return(substr(arg,8,nchar(arg))),NA)	}))
		if(length(tmp)>0) indir<- tmp[1]
		tmp<- na.omit(sapply(argv,function(arg)
						{	switch(substr(arg,2,7),
									infile= return(substr(arg,9,nchar(arg))),NA)	}))
		if(length(tmp)>0) infile<- tmp[1]		
		tmp<- na.omit(sapply(argv,function(arg)
						{	switch(substr(arg,2,7),
									outdir= return(substr(arg,9,nchar(arg))),NA)	}))
		if(length(tmp)>0) outdir<- tmp[1]
		tmp<- na.omit(sapply(argv,function(arg)
						{	switch(substr(arg,2,8),
									outfile= return(substr(arg,10,nchar(arg))),NA)	}))
		if(length(tmp)>0) outfile<- tmp[1]		
		tmp<- na.omit(sapply(argv,function(arg)
						{	switch(substr(arg,2,16),
									alignment.start= return(substr(arg,18,nchar(arg))),NA)	}))
		if(length(tmp)>0) alignment.start<- tmp[1]	
		tmp<- na.omit(sapply(argv,function(arg)
						{	switch(substr(arg,2,2),
									v= return(as.numeric(substr(arg,4,nchar(arg)))),NA)	}))
		if(length(tmp)>0) verbose<- tmp[1]
	}
	if(verbose)
	{
		print(indir)
		print(infile)
		print(outdir)
		print(outfile)		
		print(alignment.start)
		print(verbose)
	}	
	#load alignment
	if(!grepl('.R',infile))							stop("expect R infile that ends in .R")		
	file			<- paste(indir,'/',infile,sep='')
	tmp				<- load(file)
	if(verbose)	cat(paste('\nloaded file=', tmp))
	if(tmp!='seq')
		eval(parse(text=paste("seq<- ",tmp,sep='')))
	if(!"DNAbin"%in%class(seq) || !is.matrix(seq))	stop("expect R infile that contains a DNAbin matrix")

	#modify dr table for particular alignment	
	set(dr, NULL, "Alignment.nuc.pos", dr[,Alignment.nuc.pos]-alignment.start+1)	
	#if alignment equals any of the drug resistance mutants, replace with NNN	
	seq				<- seq.rm.drugresistance(as.character(seq), dr, verbose=verbose, rtn.DNAbin=1 )		
	#save No Drug resistance alignment to file
	file			<- paste(outdir,'/',outfile,'_',".R",sep='')
	if(verbose)	cat(paste("\nwrite R file to",file))
	save(seq, file=file)	
	file			<- paste(outdir,'/',outfile,".phylip",sep='')
	if(verbose)	cat(paste("\nwrite phylip file to",file))
	seq.write.dna.phylip(seq, file=file)					
	file			<- paste(outdir,'/',outfile,".fasta",sep='')			
	if(verbose)	cat(paste("\nwrite fasta file to",file))
	write.dna(seq, file=file, format="fasta", colsep='', colw=ncol(seq), blocksep=0)	
	seq
}
