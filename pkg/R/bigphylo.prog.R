
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
prog.remove.resistancemut<- function()
{
	library(big.phylo)	
	
	#load drug resistance mutations and select unique mutants by codon
	dr		<- as.data.table( read.csv( paste( CODE.HOME,"/data/IAS_primarydrugresistance_201303.csv",sep='' ), stringsAsFactors=F ) )	
	dr[,Alignment.nuc.pos:= (Gene.codon.number-1)*3+Gene.HXB2pos ]		
	dr		<- dr[,	{	tmp<- unique(Mutant); list(Mutant=tmp, Gene.codon.number=Gene.codon.number[1], Wild.type=Wild.type[1], DR.name=DR.name[1])	}, by=Alignment.nuc.pos]
	#select nucleotide codes that are consistent with drug resistance mutants
	nt2aa	<- as.data.table( read.csv( paste( CODE.HOME,"/data/standard_nt_code.csv",sep='' ), stringsAsFactors=F ) )
	setnames(nt2aa,c("AA","NTs"),c("Mutant","Mutant.NTs"))
	nt2aa	<- subset(nt2aa, select=c(Mutant,Mutant.NTs))
	dr		<- merge(dr, nt2aa, all.x=1, by="Mutant", allow.cartesian=TRUE)
	setkey(dr, "Alignment.nuc.pos")
	#print(dr, nrows=250)
	dr		<- subset(dr, select=c(Alignment.nuc.pos, Mutant.NTs, DR.name))
	set(dr, NULL, "Mutant.NTs", tolower(dr[,Mutant.NTs]))
	
	indir			<- paste(DATA,"tmp",sep='/')
	infile			<- "ATHENA_2013_03_CurAll+LANL_Sequences"
	insignat		<- "Sat_Jun_16_17/23/46_2013"
	outdir			<- paste(DATA,"tmp",sep='/')
	outfile			<- "ATHENA_2013_03_NoDRAll+LANL_Sequences"
	outsignat		<- "Thu_Aug_01_17/05/23_2013"
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
						{	switch(substr(arg,2,9),
									insignat= return(substr(arg,11,nchar(arg))),NA)	}))
		if(length(tmp)>0) insignat<- tmp[1]		
		tmp<- na.omit(sapply(argv,function(arg)
						{	switch(substr(arg,2,7),
									outdir= return(substr(arg,9,nchar(arg))),NA)	}))
		if(length(tmp)>0) outdir<- tmp[1]
		tmp<- na.omit(sapply(argv,function(arg)
						{	switch(substr(arg,2,8),
									outfile= return(substr(arg,10,nchar(arg))),NA)	}))
		if(length(tmp)>0) outfile<- tmp[1]		
		tmp<- na.omit(sapply(argv,function(arg)
						{	switch(substr(arg,2,10),
									outsignat= return(substr(arg,12,nchar(arg))),NA)	}))
		if(length(tmp)>0) outsignat<- tmp[1]	
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
		print(insignat)
		print(outdir)
		print(outfile)		
		print(outsignat)
		print(alignment.start)
		print(verbose)
	}	
	#load alignment
	file				<- paste(indir,'/',infile,'_',gsub('/',':',insignat),".R",sep='')
	load(file)	
	#modify dr table for particular alignment	
	set(dr, NULL, "Alignment.nuc.pos", dr[,Alignment.nuc.pos]-alignment.start+1)
	
	#remove	likely.nonB.outliers
	cat("\nchange infile: remove	likely.nonB.outliers")
	likely.nonB.outliers	<- c("R03-07193","2006G206","PROT+P51_B.AU.1995.C92.AF538307","2008G084")
	likely.nonB.outliers	<- which(rownames(seq.PROT.RT) %in% likely.nonB.outliers)
	seq.PROT.RT				<- seq.PROT.RT[-likely.nonB.outliers,]
	
	#if alignment equals any of the drug resistance mutants, replace with NNN	
	seq.PROT.RT			<- seq.rm.drugresistance(as.character(seq.PROT.RT), dr, verbose=verbose, rtn.DNAbin=1 )	
	
	#save No Drug resistance alignment to file
	file								<- paste(outdir,'/',outfile,'_',gsub('/',':',outsignat),".R",sep='')
	if(verbose)	cat(paste("\nwrite R file to",file))
	save(seq.PROT.RT, file=file)	
	file								<- paste(outdir,'/',outfile,'_',gsub('/',':',outsignat),".phylip",sep='')
	if(verbose)	cat(paste("\nwrite phylip file to",file))
	seq.write.dna.phylip(seq.PROT.RT, file=file)					
	file								<- paste(outdir,'/',outfile,'_',gsub('/',':',outsignat),".fasta",sep='')			
	if(verbose)	cat(paste("\nwrite fasta file to",file))
	write.dna(seq.PROT.RT, file=file, format="fasta", colsep='', colw=ncol(seq.PROT.RT), blocksep=0)
	
	seq.PROT.RT
}
