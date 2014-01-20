######################################################################################
#create PROT+RT data set of first sequences from all patients
#' @export
prog.examl.getbootstrapseq<- function(check.any.bs.identical=0)
{	
	library(ape)
	library(data.table)
	library(hivclust)
	require(recombination.analyzer)
	
	indir				<- outdir		<- paste(DATA,"tmp",sep='/')
	infile				<- "ATHENA_2013_03_FirstCurSequences_PROTRT"
	signat.out			<- signat.in	<- "Sat_May_11_14/23/46_2013"
	verbose				<- resume		<- 1
	opt.bootstrap.by	<- "codon"
	bs					<- 0
	
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
						{	switch(substr(arg,2,9),
									insignat= return(substr(arg,11,nchar(arg))),NA)	}))
		if(length(tmp)>0) signat.in<- tmp[1]
		tmp<- na.omit(sapply(argv,function(arg)
						{	switch(substr(arg,2,10),
									outsignat= return(substr(arg,12,nchar(arg))),NA)	}))
		if(length(tmp)>0) signat.out<- tmp[1]		
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
	if(1)
	{
		print( indir ) 
		print(outdir)
		print(infile)
		print(signat.in)
		print(signat.out)
		print(verbose)
		print(resume)
		print(bs)
		print(opt.bootstrap.by)
		print(signat.in)
		print(signat.out)
	}
	if(!opt.bootstrap.by%in%c("nucleotide","codon"))	stop("Unexpected opt.bootstrap.by")		
	pattern 	<- paste(infile,"_",gsub('/',':',signat.out),".phylip.",sprintf("%03d",bs),sep='')
	file		<- list.files(path=outdir, pattern=pattern, full.names=1)
	if(!resume || !length(file))	
	{					
		file		<- paste(outdir,"/",infile,"_",gsub('/',':',signat.out),".R",sep='')
		if(verbose) cat(paste("\nread",file))
		tmp			<- load(file)
		if(length(tmp)!=1)	stop("Unexpected lenght of loaded objects")
		eval(parse(text=paste("seq.PROT.RT<- ",tmp,sep='')))		
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
		file		<- paste(outdir,"/",infile,"_",gsub('/',':',signat.out),".phylip.",sprintf("%03d",bs),sep='')
		cat(paste("\nsave boostrap seq alignment to",file))
		seq.write.dna.phylip(seq.BS, file=file)
	}
	else
		cat("\nfound boostrap sequence alignment")
}
######################################################################################
#' @title Return \code{3SEQ} output as a data.table 
#' @description Collect \code{3SEQ} output of multiple, consecutive job files for a given set of sequences 
#' in file \code{infile} in directory \code{indir}. The \code{3SEQ} output must cover all sequences in 
#' \code{infile}. Input parameters are 'indir', 'infile', 'resume' and 'verbose', and specified via an 
#' \code{argv} string, see the Examples.
#' @return Data table of potential recombinant sequences and associated parent sequences that are identified at a 
#' given p-value that is corrected for multiple comparisons.
#' @seealso \code{\link{pipeline.recom.run.3seq}}
#' @example example/package.mtDNA.process.3seq.R
#' @export
prog.recom.process.3SEQ.output<- function()
{	
	require(recombination.analyzer)
	
	verbose		<- 1
	resume		<- 1
	indir		<- getwd()		
	infile		<- "ATHENA_2013_03_NoDRAll+LANL_Sequences"			
	
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
		#		
		tmp<- na.omit(sapply(argv,function(arg)
						{	switch(substr(arg,2,7),
									resume= return(as.numeric(substr(arg,9,nchar(arg)))),NA)	}))
		if(length(tmp)>0) resume<- tmp[1]
		tmp<- na.omit(sapply(argv,function(arg)
						{	switch(substr(arg,2,2),
									v= return(as.numeric(substr(arg,4,nchar(arg)))),NA)	}))
		if(length(tmp)>0) verbose<- tmp[1]		
	}	
	if(verbose)
	{
		print(verbose)
		print(resume)
		print(indir)		
		print(infile)				
	}
	#load sequences as 'seq' object
	if(!grepl('.R',infile))							stop("expect R infile that ends in .R")
	file		<- paste(indir,'/',infile,sep='')
	if(verbose)	cat(paste("\nload file ",file))			
	tmp			<- load(file)
	if(verbose)	cat(paste('\nloaded file=', tmp))
	if(tmp!='seq')
		eval(parse(text=paste("seq<- ",tmp,sep='')))
	if(!"DNAbin"%in%class(seq) || !is.matrix(seq))	stop("expect R infile that contains a DNAbin matrix")
	infile		<- substr(infile, 1, nchar(infile)-2)
	#resume if possible
	file		<- paste(indir,'/',infile,"_3seq",".R",sep='')
	if(resume)
	{				
		options(show.error.messages = FALSE)		
		if(verbose)	cat(paste("\ntry to resume file ",file))
		readAttempt<-	try(suppressWarnings(load(file)))
		options(show.error.messages = TRUE)
		if(!inherits(readAttempt, "try-error") && verbose)	cat(paste("\nresumed file ",file))			
	}
	if(!resume || inherits(readAttempt, "try-error"))
	{
		if(verbose)	cat(paste("\ngenerate file ",file))
		#	collect 3SEQ output files
		tmp			<- list.files(indir, pattern="3seq$")
		files		<- tmp[grepl(infile, tmp, fixed=1)]
		tmp			<- sapply(strsplit(files,'.', fixed=1), function(x) x[[1]] )
		tmp			<- sapply(strsplit(tmp,'_'), function(x) tail(x,1) )		#see if filename has '_xx-xx' string at last pos  
		files		<- files[grepl('-', tmp)]		
		if(verbose)	cat(paste("\nFound 3seq output files matching infile and insignat, n=",length(files)))
		#	figure out if any consecutive files are missing	
		tmp			<- tmp[grepl('-', tmp)]
		df.seqchunks<- as.data.table( t( sapply( strsplit(tmp, '-'), as.numeric) ) )
		setnames(df.seqchunks, c("V1","V2"), c("start","end"))
		setkey(df.seqchunks, start)
		tmp			<- subset(df.seqchunks, df.seqchunks[-1,start]-1 != df.seqchunks[-nrow(df.seqchunks),end] )
		if(nrow(tmp)){		print(tmp); stop("Found missing sequence chunks")		}
		#	determine non-empty files and number of columns		
		tmp			<- sapply(files, function(x)
				{
					tmp	<- count.fields(paste(indir, '/', x, sep=''), skip=1, sep='\t')
					ifelse(length(tmp), max(tmp), 0)			
				})
		col.max		<- max(tmp)		
		files		<- files[tmp>0]
		if(verbose)	cat(paste("\nFound non-empty 3seq files, n=",length(files)))
		#	read non-empty files
		col.names		<- c("parent1","parent2","child","m","n","k","p","[p_max]","HS?","log(p)","DS(p)","DS(p)","min_rec_length")
		if(col.max>length(col.names))
			col.names	<- c(col.names, paste("bp",seq_len( col.max-length(col.names) ),sep=''))
		df.recomb		<- lapply(files, function(x)
				{	
					if(verbose)	cat(paste("\nprocess file", x))
					tmp					<- read.delim(paste(indir, '/', x, sep=''), skip=1, header=0, fill=1, col.names=col.names, stringsAsFactors=0, strip.white=T)					
					as.data.table(tmp)			
				})
		df.recomb		<- rbindlist(df.recomb)
		#	extract FASTASampleCodes from job output
		tmp				<- df.recomb[, list(m1= regexpr("-[",parent1,fixed=1), m2=regexpr("-[",parent2,fixed=1), mc=regexpr("-[",child,fixed=1)) ]
		if(any(tmp<1))	stop("Unexpected parent1 or parent2. Parsing error?")
		df.recomb		<- cbind( df.recomb, tmp )
		df.recomb[,dummy:=seq_len(nrow(df.recomb))]
		set(df.recomb, NULL, "parent1", df.recomb[,substr(parent1, 1, m1-1), by="dummy"][,V1])
		set(df.recomb, NULL, "parent2", df.recomb[,substr(parent2, 1, m2-1), by="dummy"][,V1])
		set(df.recomb, NULL, "child", df.recomb[,substr(child, 1, mc-1), by="dummy"][,V1])
		#	extract unique recombinants from job output
		setkey(df.recomb,child)
		df.recomb		<- unique(df.recomb)
		if(verbose)	cat(paste("\nFound recombinant sequences, n=",nrow(df.recomb)))
		#	order parents
		tmp				<- df.recomb[, {
					z<- sort(c(parent1, parent2))
					list(parent1=z[1], parent2=z[2])
				} ,by="dummy"]
		df.recomb		<- df.recomb[, setdiff( colnames(df.recomb), c("parent1","parent2","p","X.p_max.","HS.","DS.p.","m1","m2","mc") ), with=F]
		df.recomb		<- merge(tmp, df.recomb, by="dummy")
		#	check if triplets are unique
		tmp				<- df.recomb[, {
					z<- sort(c(parent1, parent2, child))
					list(parent1=z[1], parent2=z[2], child=z[3])
				} ,by="dummy"]
		setkeyv(tmp, c("parent1","parent2","child"))
		tmp				<- unique(tmp)
		if(nrow(tmp)<nrow(df.recomb))	warning("triplets in df.recomb not unique")
		#
		tmp				<- subset(df.recomb, select=c(parent1,parent2,child))
		setkey(tmp, parent1, parent2)
		tmp				<- unique(tmp)
		tmp[, parentpair:=seq_len(nrow(tmp))]				
		if(verbose)	cat(paste("\nNumber of unique parent pairs, n=",nrow(tmp)))
		#
		tmp				<- unique(c(df.recomb[, parent1], df.recomb[, parent2]))
		if(verbose)	cat(paste("\nNumber of unique parent sequences, n=",length(tmp)))
		#
		if(length(unique(df.recomb[,child]))!=nrow(df.recomb))	stop("Unexpected duplicate children - select parents with smallest p?")
		#
		tmp				<- subset(df.recomb, select=c(parent1,parent2))
		setkey(tmp, parent1, parent2)
		tmp				<- unique(tmp)
		tmp[, parentpair:=seq_len(nrow(tmp))]				
		if(verbose)	cat(paste("\nNumber of unique parent pairs, n=",nrow(tmp)))
		df.recomb		<- merge(df.recomb, tmp, by=c("parent1","parent2"))
		setnames(df.recomb,c("log.p.","DS.p..1"),c("logp","adjp"))
		#	candidate breakpoints bp1 etc are all overlapping, consider only bp1 for simplicity	
		#	evaluate midpoint of breapoint regions		
		tmp<- strsplit(df.recomb[,bp1], ',')
		if(any(sapply(tmp, length )!=2)) 	stop("\nunexpected bp1: missing ','")		
		df.recomb[, bp1.1:= sapply(tmp, "[", 1)]
		df.recomb[, bp1.2:= sapply(tmp, "[", 2)]
		set(df.recomb, NULL, "bp1.1", round( sapply(strsplit(df.recomb[,bp1.1], '-'), function(x)	mean(as.numeric(x))	) ) )
		set(df.recomb, NULL, "bp1.2", round( sapply(strsplit(df.recomb[,bp1.2], '-'), function(x)	min(as.numeric(x))	) ) )		
		df.recomb	<- merge(df.recomb, df.recomb[, list(child.len=max(which(as.character(seq[child,])!='-'))), by="dummy"], by="dummy" )
		if(any(df.recomb[, child.len-bp1.2]<0))	stop("\nunexpected breakpoint past end of child sequence")
		df.recomb[, child.start:= 1]
		tmp			<- which( df.recomb[, bp1.1<10] )
		set(df.recomb, tmp, "child.start", df.recomb[tmp, bp1.1])
		tmp			<- which( df.recomb[, child.len-bp1.2<25] )
		set(df.recomb, tmp, "child.len", as.integer(df.recomb[tmp, bp1.2]))
		#	save potential recombinants
		file		<- paste(indir,'/',infile,"_3seq",".R",sep='')
		if(verbose)	cat(paste("\nSave candidate triplets to",file))
		save(df.recomb, file=file)		
	}
	if(0)
	{
		setnames(df.recomb, "dummy", "triplet.id")
		#
		#	get candidate recombinant sequences
		#
		df.recombseq	<- data.table( FASTASampleCode= unique( c(df.recomb[, parent1], df.recomb[, parent2], df.recomb[, child]) ) )
		df.recombseq	<- df.recombseq[ , 	{
												tmp<- subset(df.recomb, parent1==FASTASampleCode | parent2==FASTASampleCode | child==FASTASampleCode )
												list( n.triplets=nrow(tmp), min.p=min(tmp[,adjp]), med.p=median(tmp[,adjp]), triplet.id=tmp[,triplet.id] )
											} , by="FASTASampleCode"]
		setkey(df.recombseq, n.triplets, FASTASampleCode)
		#plot( df.recombseq[,n.triplets], df.recombseq[,min.p], pch=18 )
		
		#
		#determine how many other triplet sequences there are for an m candidate
		#
		df.mrecombseq	<- subset(df.recombseq, n.triplets>2)[, {
					tmp			<- triplet.id
					tmp			<- subset(df.recomb, triplet.id%in%tmp)
					triplet.seq	<- setdiff( unique( c(tmp[, parent1], tmp[, parent2], tmp[, child]) ), FASTASampleCode)
					list( triplet.seq=triplet.seq, triplet.seq.n=length(triplet.seq)  )
				}, by="FASTASampleCode"]
		df.mrecombbp	<- subset(df.recombseq, n.triplets>2)[, {
					tmp			<- triplet.id
					tmp			<- subset(df.recomb, triplet.id%in%tmp)
					list( bp1= tmp[,bp1], bp1.1= tmp[,bp1.1], bp1.2= tmp[,bp1.2] )																	
				}, by="FASTASampleCode"]														
		setkeyv(df.mrecombbp, c("FASTASampleCode", "bp1.1", "bp1.2"))												
		#	breakpoints among mrecombinants are not necessarily the same
		print(df.mrecombbp, n=250)										
	}
	df.recomb
}
######################################################################################
#' @export
prog.recom.plot.incongruence<- function()
{
	require(RColorBrewer)
	require(ape)
	require(recombination.analyzer)
	
	verbose		<- 1
	resume		<- 1
	indir		<- paste(DATA,"tmp",sep='/')		
	infile		<- "ATHENA_2013_03_NoDRAll+LANL_Sequences"			
	insignat	<- "Thu_Aug_01_17/05/23_2013"
	
	id			<- NA		
	bs.n		<- 500
	select		<- ''
	#select		<- 'ng2'
	
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
		#		
		tmp<- na.omit(sapply(argv,function(arg)
						{	switch(substr(arg,2,2),
									v= return(as.numeric(substr(arg,4,nchar(arg)))),NA)	}))
		if(length(tmp)>0) verbose<- tmp[1]
		#
		tmp<- na.omit(sapply(argv,function(arg)
						{	switch(substr(arg,2,10),
									tripletid= return(as.numeric(substr(arg,12,nchar(arg)))),NA)	}))
		if(length(tmp)>0) id<- tmp[1]
		tmp<- na.omit(sapply(argv,function(arg)
						{	switch(substr(arg,2,7),
									select= return(substr(arg,9,nchar(arg))),NA)	}))
		if(length(tmp)>0) select<- tmp[1]
	}	
	if(verbose)
	{
		print(verbose)		
		print(indir)		
		print(infile)			
		print(insignat)
		print(id)		
		print(bs.n)		
		print(select)
	}
	cols				<- brewer.pal(6,"Paired")
	pattern				<- c("^in_parent1","^out_parent1","^in_parent2","^out_parent2","^in_child","^out_child")
	cex					<- 0.6
	thresh.nodesupport	<- 0.6
	edge.length.outliers<- 0.5
	
	if(is.na(id))
	{
		#	read candidate triplets
		argv				<<-	cmd.recombination.process.3SEQ.output(indir, infile, insignat, resume=1, verbose=1) 
		argv				<<- unlist(strsplit(argv,' '))
		df.recomb			<- prog.recom.process.3SEQ.output()
		setnames(df.recomb, "dummy", "triplet.id")
		setkey(df.recomb, triplet.id)
		#	get candidate recombinant sequences
		df.recombseq	<- data.table( FASTASampleCode= unique( c(df.recomb[, parent1], df.recomb[, parent2], df.recomb[, child]) ) )
		df.recombseq	<- df.recombseq[ , 	{
					tmp<- subset(df.recomb, parent1==FASTASampleCode | parent2==FASTASampleCode | child==FASTASampleCode )
					list( n.triplets=nrow(tmp), min.p=min(tmp[,adjp]), med.p=median(tmp[,adjp]), triplet.id=tmp[,triplet.id] )
				} , by="FASTASampleCode"]
		setkey(df.recombseq, FASTASampleCode)
		#	select unclear triplets	
		if(select=="ng2")
		{
			tmp				<- unique(subset(df.recombseq, n.triplets>2)[, triplet.id])			
			df.recomb		<- df.recomb[J(setdiff(df.recomb[,triplet.id], tmp )),]
			if(verbose)	cat(paste("\nSelected unclear triplets, n=", nrow(df.recomb)))
		}
		#	select triplets that have a given FASTASampleCode
		if(select %in% unique(subset(df.recombseq, n.triplets>2)[, FASTASampleCode]))
		{
			tmp				<- subset(df.recombseq, n.triplets>2)[J(select),]
			df.recomb		<- df.recomb[J(tmp[,triplet.id]),]
		}
		#	analyze FASTASampleCodes that occur in >2 triplets
		if(select=="g2")
		{
			tmp				<- unique(subset(df.recombseq, n.triplets>2)[, FASTASampleCode])
			dummy			<- lapply(tmp, function(x)
						{
							argv			<<- cmd.recombination.plot.incongruence(indir, infile, insignat, prog= PR.RECOMB.PLOTINCONGRUENCE, opt.select=x,verbose=1)
							argv			<<- unlist(strsplit(argv,' '))
							prog.recom.plot.incongruence()			
						})
			stop()
		}
		
		#	read available checks for triplets
		files				<- list.files(indir, pattern=".newick$")
		files				<- files[ grepl(paste('_3seqcheck_',sep=''), files) & grepl(infile, files, fixed=1) & grepl(gsub('/',':',insignat), files) ]
		tmp					<- regexpr("3seqcheck_",files)
		tmp					<- sapply(seq_along(files), function(i){		substr(files[i],tmp[i],nchar(files[i]))		})
		files.df			<- data.table(	file=files, 
											triplet.id= sapply(strsplit(tmp, '_'), function(x)		substr(x[[2]],3,nchar(x[[2]]))	),
											region= sapply(strsplit(tmp, '_'), function(x)		substr(x[[3]],2,nchar(x[[3]]))	)		)
		set(files.df, NULL, "triplet.id", as.numeric(files.df[,triplet.id]))							
		setkey(files.df, triplet.id, region)
		if(verbose)	cat(paste("\nFound files, n=", nrow(files.df)))
		files.df			<- merge(files.df, files.df[,list(region.n= length(region)), by="triplet.id"], by="triplet.id")
		files.df			<- subset(files.df, region.n>1)
		if(verbose)	cat(paste("\nFound checked triplets, n=", nrow(files.df)/2))
		if(verbose) cat(paste("\nTriplets still to check=",paste(setdiff( df.recomb[,triplet.id],unique(files.df[,triplet.id]) ), collapse=', ')))
		#	select candidate triplets for which checks available
		files.df			<- merge(files.df, df.recomb,by="triplet.id")
		if(verbose)	cat(paste("\nSelected triplets for plotting, n=", nrow(files.df)/2))
		setkey(files.df, parentpair, adjp)
		#setkey(files.df, adjp)
		
		file				<- paste(indir,'/',infile,"_3seqcheck_examlbs",bs.n,'_',select,'_',gsub('/',':',insignat),".pdf",sep='')
		if(verbose)	cat(paste("\nPlot both phylogenies to", file))		
		pdf(file, width=12, height=6)
		def.par 			<- par(no.readonly = TRUE)
		par(mar=c(2,0.5,1,0))
		layout(matrix(c(1,1,2,2), 2, 2))					
		dummy				<- lapply( unique( files.df[, triplet.id] ), function(z)
				{
					#x<- subset(files.df, triplet.id==58)
					x<- subset(files.df, triplet.id==z)
					#	read tree corresponding to 'in' recombinant region
					tmp					<- paste(indir,'/',x[1, file],sep='')
					if(verbose)	cat(paste("\nRead file", tmp))
					ph.in				<- ladderize( read.tree(tmp) )		
					tmp					<- as.numeric( ph.in$node.label )
					tmp[is.na(tmp)]		<- 0 
					ph.in$node.label	<- tmp/100
					#	read tree corresponding to 'out' recombinant region
					tmp					<- paste(indir,'/',x[2, file],sep='')		
					if(verbose)	cat(paste("\nRead file", tmp))
					ph.out				<- ladderize( read.tree(tmp) )		
					tmp					<- as.numeric( ph.out$node.label )
					tmp[is.na(tmp)]		<- 0 
					ph.out$node.label	<- tmp/100
					#	remove outlier filler sequences
					outliers			<- c()
					tmp					<- which( ph.in$edge.length>edge.length.outliers )
					if(length(tmp))
					{
						tmp				<- ph.in$edge[tmp,2]
						outliers		<- ph.in$tip.label[ tmp[tmp<Ntip(ph.in)] ]
					}
					if(length(tmp))
					{						
						tmp				<- which( ph.out$edge.length>1 )
						tmp				<- ph.out$edge[tmp,2]
						outliers		<- c(outliers, ph.out$tip.label[ tmp[tmp<Ntip(ph.out)] ])
					}
					if(length(outliers))
					{
						ph.in			<- drop.tip(ph.in, outliers)
						ph.out			<- drop.tip(ph.out, outliers)
					}
					#	show half way stable subtrees -- ideally would contain triplet sequences
					clustering.in		<- hivc.clu.clusterbythresh(ph.in, thresh.nodesupport=thresh.nodesupport, nodesupport=ph.in$node.label, retval="all")
					clustering.out		<- hivc.clu.clusterbythresh(ph.out, thresh.nodesupport=thresh.nodesupport, nodesupport=ph.out$node.label, retval="all")
					#	show filler sequences in different color
					tip.color.in		<- rep("black", Ntip(ph.in))	
					for(i in seq_along(pattern))
					{
						tmp						<- which( grepl(pattern[i],ph.in$tip) )
						if(length(tmp))	
							tip.color.in[ tmp ]	<- cols[i]
					}
					tip.color.out		<- rep("black", Ntip(ph.out))	
					for(i in seq_along(pattern))
					{
						tmp						<- which( grepl(pattern[i],ph.out$tip) )
						if(length(tmp))	
							tip.color.out[ tmp ]<- cols[i]
					}	
					#	plot both phylogenies side by side	
					dummy<- hivc.clu.plot(ph.in, clustering.in[["clu.mem"]], show.tip.label=T, cex.nodelabel=cex, cex.edge.incluster=2*cex, no.margin=F, tip.color=tip.color.in)
					mtext(paste("parent.pair=",x[,parentpair]," triplet.id=",x[1,triplet.id]," log10p=",round(log10(x[1,adjp]),d=2)," region 'in' length=",x[1,bp1.2-bp1.1]), side = 3, cex=cex)
					axisPhylo(cex=cex)
					dummy<- hivc.clu.plot(ph.out, clustering.out[["clu.mem"]], show.tip.label=T, cex.nodelabel=cex, cex.edge.incluster=2*cex, no.margin=F, tip.color=tip.color.out)		
					mtext(paste("parent.pair=",x[,parentpair]," triplet.id=",x[1,triplet.id]," log10p=",round(log10(x[1,adjp]),d=2)," region 'out' length=",x[1,child.len-child.start-bp1.2+bp1.1]), side = 3, cex=cex)
					axisPhylo(cex=cex)
				})
		par(def.par)
		dev.off()
	
		#x<- "R12-15108"
		#subset(df.recomb, parent1==x | parent2==x | child==x )
	}
	if(!is.na(id))
	{
		#	read tree corresponding to 'in' recombinant region
		file				<- paste(indir,'/',infile,"_3seqcheck_id",id,"_rIn_examlbs",bs.n,'_',gsub('/',':',insignat),".newick",sep='')
		if(verbose)	cat(paste("\nRead file", file))
		ph.in				<- ladderize( read.tree(file) )		
		tmp					<- as.numeric( ph.in$node.label )
		tmp[is.na(tmp)]		<- 0 
		ph.in$node.label	<- tmp/100
		#	read tree corresponding to 'out' recombinant region
		file				<- paste(indir,'/',infile,"_3seqcheck_id",id,"_rOut_examlbs",bs.n,'_',gsub('/',':',insignat),".newick",sep='')
		if(verbose)	cat(paste("\nRead file", file))
		ph.out				<- ladderize( read.tree(file) )		
		tmp					<- as.numeric( ph.out$node.label )
		tmp[is.na(tmp)]		<- 0 
		ph.out$node.label	<- tmp/100
		#	remove outlier filler sequences
		outliers			<- c()
		tmp					<- which( ph.in$edge.length>edge.length.outliers )
		if(length(tmp))
		{
			tmp				<- ph.in$edge[tmp,2]
			outliers		<- ph.in$tip.label[ tmp[tmp<Ntip(ph.in)] ]
		}
		if(length(tmp))
		{						
			tmp				<- which( ph.out$edge.length>1 )
			tmp				<- ph.out$edge[tmp,2]
			outliers		<- c(outliers, ph.out$tip.label[ tmp[tmp<Ntip(ph.out)] ])
		}
		if(length(outliers))
		{
			ph.in			<- drop.tip(ph.in, outliers)
			ph.out			<- drop.tip(ph.out, outliers)
		}
		#	show half way stable subtrees -- ideally would contain triplet sequences
		clustering.in		<- hivc.clu.clusterbythresh(ph.in, thresh.nodesupport=thresh.nodesupport, nodesupport=ph.in$node.label, retval="all")
		clustering.out		<- hivc.clu.clusterbythresh(ph.out, thresh.nodesupport=thresh.nodesupport, nodesupport=ph.out$node.label, retval="all")
		#	show filler sequences in different color
		tip.color.in		<- rep("black", Ntip(ph.in))	
		for(i in seq_along(pattern))
		{
			tmp						<- which( grepl(pattern[i],ph.in$tip) )
			if(length(tmp))	
				tip.color.in[ tmp ]	<- cols[i]
		}
		tip.color.out		<- rep("black", Ntip(ph.out))	
		for(i in seq_along(pattern))
		{
			tmp						<- which( grepl(pattern[i],ph.out$tip) )
			if(length(tmp))	
				tip.color.out[ tmp ]<- cols[i]
		}	
		#	plot both phylogenies side by side
		file				<- paste(indir,'/',infile,"_3seqcheck_id",id,"_examlbs",bs.n,'_',gsub('/',':',insignat),".pdf",sep='')
		if(verbose)	cat(paste("\nPlot both phylogenies to", file))
		pdf(file, width=8, height=6)
		def.par 			<- par(no.readonly = TRUE)
		layout(matrix(c(1,1,2,2), 2, 2))		
		dummy<- hivc.clu.plot(ph.in, clustering.in[["clu.mem"]], show.tip.label=T, cex.nodelabel=cex, cex.edge.incluster=2*cex, no.margin=F, tip.color=tip.color.in)
		mtext(paste("parent.pair=",x[,parentpair]," triplet.id=",x[1,triplet.id]," log10p=",round(log10(x[1,adjp]),d=2)," region 'in' length=",x[1,bp1.2-bp1.1]), side = 3, cex=cex)
		axisPhylo(cex=cex)
		dummy<- hivc.clu.plot(ph.out, clustering.out[["clu.mem"]], show.tip.label=T, cex.nodelabel=cex, cex.edge.incluster=2*cex, no.margin=F, tip.color=tip.color.out)		
		mtext(paste("parent.pair=",x[,parentpair]," triplet.id=",x[1,triplet.id]," log10p=",round(log10(x[1,adjp]),d=2)," region 'out' length=",x[1,child.len-child.start-bp1.2+bp1.1]), side = 3, cex=cex)
		axisPhylo(cex=cex)
		par(def.par)
		dev.off()
	}
}
######################################################################################
#' @export
prog.recom.get.incongruence<- function()
{	
	require(ape)
	require(data.table)
	require(recombination.analyzer)
	#default arguments
	verbose		<- 1
	resume		<- 0
	seq.select.n<- 10
	bs.from		<- 0
	bs.to		<- 499
	bs.n		<- 500	
	hpc.walltime<- 36
	hpc.mem		<- "600mb"
	hpc.nproc	<- 1		
	hpc.q		<- "pqeph"
	#default job
	indir		<- paste(DATA,"tmp",sep='/')		
	infile		<- "ATHENA_2013_03_NoDRAll+LANL_Sequences"			
	insignat	<- "Thu_Aug_01_17/05/23_2013"	
	id			<- 51	
	
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
		#		
		tmp<- na.omit(sapply(argv,function(arg)
						{	switch(substr(arg,2,7),
									resume= return(as.numeric(substr(arg,9,nchar(arg)))),NA)	}))
		if(length(tmp)>0) resume<- tmp[1]
		tmp<- na.omit(sapply(argv,function(arg)
						{	switch(substr(arg,2,2),
									v= return(as.numeric(substr(arg,4,nchar(arg)))),NA)	}))
		if(length(tmp)>0) verbose<- tmp[1]
		#
		tmp<- na.omit(sapply(argv,function(arg)
						{	switch(substr(arg,2,10),
									tripletid= return(as.numeric(substr(arg,12,nchar(arg)))),NA)	}))
		if(length(tmp)>0) id<- tmp[1]
		#
		tmp<- na.omit(sapply(argv,function(arg)
						{	switch(substr(arg,2,13),
									hpc.walltime= return(substr(arg,15,nchar(arg))),NA)	}))
		if(length(tmp)>0) hpc.walltime<- tmp[1]
		tmp<- na.omit(sapply(argv,function(arg)
						{	switch(substr(arg,2,8),
									hpc.mem= return(substr(arg,10,nchar(arg))),NA)	}))
		if(length(tmp)>0) hpc.mem<- tmp[1]
		tmp<- na.omit(sapply(argv,function(arg)
						{	switch(substr(arg,2,6),
									hpc.q= return(substr(arg,8,nchar(arg))),NA)	}))
		if(length(tmp)>0) hpc.q<- tmp[1]
		tmp<- na.omit(sapply(argv,function(arg)
						{	switch(substr(arg,2,10),
									hpc.nproc= return(as.numeric(substr(arg,12,nchar(arg)))),NA)	}))
		if(length(tmp)>0) hpc.nproc<- tmp[1]
	}	
	if(verbose)
	{
		print(verbose)
		print(resume)
		print(indir)		
		print(infile)			
		print(insignat)
		print(id)
		print(seq.select.n)
		print(bs.from)
		print(bs.to)
		print(bs.n)
		print(hpc.walltime)
		print(hpc.mem)
		print(hpc.nproc)		
		print(hpc.q)
	}
	if(resume)
	{
		tmp			<- 1
		file		<- paste(indir,'/',infile,"_3seqcheck_id",id,"_rIn_",gsub('/',':',insignat),".R",sep='')
		options(show.error.messages = FALSE)		
		if(verbose)	cat(paste("\ntry to load file ",file))
		readAttempt	<-	try(suppressWarnings(load(file)))
		options(show.error.messages = TRUE)
		if(!inherits(readAttempt, "try-error") && verbose)		cat(paste("\nloaded file=",file))
		if(inherits(readAttempt, "try-error"))					tmp		<- 0
		
		file		<- paste(indir,'/',infile,"_3seqcheck_id",id,"_rOut_",gsub('/',':',insignat),".R",sep='')
		options(show.error.messages = FALSE)		
		if(verbose)	cat(paste("\ntry to load file ",file))
		readAttempt	<-	try(suppressWarnings(load(file)))
		options(show.error.messages = TRUE)
		if(!inherits(readAttempt, "try-error") && verbose)		cat(paste("\nloaded file=",file))
		if(inherits(readAttempt, "try-error"))					tmp		<- 0		
	}
	if(!resume || tmp==0)
	{
		file		<- paste(indir,'/',infile,'_', gsub('/',':',insignat),".R",sep='')
		if(verbose)	cat(paste("\nload file ",file))			
		load(file)
		#	loaded seq.PROT.RT
		file		<- paste(indir,'/',infile,"_3seq_", gsub('/',':',insignat),".R",sep='')		
		options(show.error.messages = FALSE)		
		if(verbose)	cat(paste("\ntry to load file ",file))
		readAttempt	<-	try(suppressWarnings(load(file)))
		options(show.error.messages = TRUE)
		if(inherits(readAttempt, "try-error"))	stop(paste("\nCannot find 3SEQ file, run prog.recom.process.3SEQ.output?, file=",file))			
		#	loaded df.recomb
		
		#
		#	process triplet for dummy id	
		#	
		df.recomb	<- subset(df.recomb, dummy==id)
		if(verbose)	cat(paste("\nprocess triplet number",id))
		if(verbose)	print(df.recomb)
		#	create sequence matrices corresponding to the two breakpoint regions 
		seq.in		<- seq.PROT.RT[,seq.int(df.recomb[,bp1.1],df.recomb[,bp1.2])]
		seq.out		<- if(df.recomb[,child.start]<df.recomb[,bp1.1]-1) seq.int(df.recomb[,child.start],df.recomb[,bp1.1]-1) else numeric(0) 
		seq.out		<- if(df.recomb[,bp1.2]+1<df.recomb[,child.len]) c(seq.out,seq.int(df.recomb[,bp1.2]+1, df.recomb[,child.len]))	else 	seq.out
		seq.out		<- seq.PROT.RT[,seq.out]
		seq.select.f<- ifelse(min(ncol(seq.out),ncol(seq.in))<150, 10, 10)
		if(verbose)	cat(paste("\nsetting inflation factor to",seq.select.f))
		seq.select.n<- seq.select.n * seq.select.f
		#	select background sequences for child based on sequence similarity
		if(verbose)	cat(paste("\ncompute genetic distances for parent1 parent2 child"))
		tmp				<- which( rownames(seq.PROT.RT)==df.recomb[,child] )		
		dummy			<- 0				
		seq.dist					<- 1 - sapply(seq_len(nrow(seq.in))[-tmp],function(i){		.C("hivc_dist_ambiguous_dna", seq.in[tmp,], seq.in[i,], ncol(seq.in), dummy )[[4]]			})	
		seq.dist[is.nan(seq.dist)]	<- Inf
		seq.df			<- data.table( FASTASampleCode=rownames(seq.in)[-tmp] , dist=seq.dist, group="child", region="in" ) 		
		seq.dist					<- 1 - sapply(seq_len(nrow(seq.out))[-tmp],function(i){		.C("hivc_dist_ambiguous_dna", seq.out[tmp,], seq.out[i,], ncol(seq.out), dummy )[[4]]			})	
		seq.dist[is.nan(seq.dist)]	<- Inf
		seq.df			<- rbind(seq.df, data.table( FASTASampleCode=rownames(seq.out)[-tmp] , dist=seq.dist, group="child", region="out" )) 
		#	select background sequences for parent1 based on sequence similarity
		tmp				<- which( rownames(seq.PROT.RT)==df.recomb[,parent1] )				
		seq.dist					<- 1 - sapply(seq_len(nrow(seq.in))[-tmp],function(i){		.C("hivc_dist_ambiguous_dna", seq.in[tmp,], seq.in[i,], ncol(seq.in), dummy )[[4]]			})	
		seq.dist[is.nan(seq.dist)]	<- Inf
		seq.df			<- rbind(seq.df,data.table( FASTASampleCode=rownames(seq.in)[-tmp] , dist=seq.dist, group="parent1", region="in" )) 		
		seq.dist					<- 1 - sapply(seq_len(nrow(seq.out))[-tmp],function(i){		.C("hivc_dist_ambiguous_dna", seq.out[tmp,], seq.out[i,], ncol(seq.out), dummy )[[4]]			})	
		seq.dist[is.nan(seq.dist)]	<- Inf
		seq.df			<- rbind(seq.df, data.table( FASTASampleCode=rownames(seq.out)[-tmp] , dist=seq.dist, group="parent1", region="out" )) 
		#	select background sequences for parent2 based on sequence similarity
		tmp				<- which( rownames(seq.PROT.RT)==df.recomb[,parent2] )				
		seq.dist					<- 1 - sapply(seq_len(nrow(seq.in))[-tmp],function(i){		.C("hivc_dist_ambiguous_dna", seq.in[tmp,], seq.in[i,], ncol(seq.in), dummy )[[4]]			})	
		seq.dist[is.nan(seq.dist)]	<- Inf
		seq.df			<- rbind(seq.df,data.table( FASTASampleCode=rownames(seq.in)[-tmp] , dist=seq.dist, group="parent2", region="in" )) 		
		seq.dist					<- 1 - sapply(seq_len(nrow(seq.out))[-tmp],function(i){		.C("hivc_dist_ambiguous_dna", seq.out[tmp,], seq.out[i,], ncol(seq.out), dummy )[[4]]			})	
		seq.dist[is.nan(seq.dist)]	<- Inf
		seq.df			<- rbind(seq.df, data.table( FASTASampleCode=rownames(seq.out)[-tmp] , dist=seq.dist, group="parent2", region="out" ))
		#	first pass:
		#	select closest FASTASampleCode by group and region to verify recombination breakpoint by phylogenetic incongruence
		setkeyv(seq.df, c("group","region","dist"))
		seq.df			<- subset(seq.df, dist>0)
		if(verbose)	cat(paste("\nFound related sequences with dist>0, n=",nrow(seq.df)))
		#	for each group and region, select the n closest sequence names (non-unique)
		seq.df			<- seq.df[	,	list(FASTASampleCode=FASTASampleCode[seq_len(seq.select.n)], dist=dist[seq_len(seq.select.n)]), by=c("group","region")]
		#	keep each sequence name once, for the group it is closest to		
		seq.df			<- seq.df[, {
										tmp<- which.min(dist)
										list(dist= dist[tmp], group=group[tmp], region=region[tmp])
									}, by=c("FASTASampleCode")]
		setkeyv(seq.df, c("group","region","dist"))					
		if(verbose)	cat(paste("\ndetermined candidates for balancing filler sequences, n=",nrow(seq.df)))
		#	select unique sequences
		tmp				<- c( df.recomb[,parent1],df.recomb[,parent2],df.recomb[,child], seq.df[, FASTASampleCode] )
		tmp				<- hivc.seq.unique(seq.in[tmp,])
		seq.df			<- merge( data.table(FASTASampleCode=rownames(tmp)), seq.df, by="FASTASampleCode" )
		tmp				<- c( df.recomb[,parent1],df.recomb[,parent2],df.recomb[,child], seq.df[, FASTASampleCode] )
		tmp				<- hivc.seq.unique(seq.out[tmp,])
		seq.df			<- merge( data.table(FASTASampleCode=rownames(tmp)), seq.df, by="FASTASampleCode" )
		seq.df			<- subset(seq.df, !FASTASampleCode%in%c( df.recomb[,parent1],df.recomb[,parent2],df.recomb[,child] ) )		
		if(verbose)	cat(paste("\nfound candidates for filler sequences that are unique on both recombinant regions, n=",nrow(seq.df)))
		if(verbose)	print( seq.df[	,	list(n=length(FASTASampleCode)) ,by=c("group","region")] )
		#		
		seq.select.n	<- seq.select.n/seq.select.f 
		#	select 'seq.select.n' unique filler sequences for 'in' region, balancing by group as much as possible
		if(verbose)	cat(paste("\nSelect filler sequences for recombinant region 'in'"))
		tmp						<- subset( seq.df, region=='in' )[,FASTASampleCode]
		tmp						<- hivc.seq.unique(seq.in[tmp ,])
		seq.in.df				<- merge( data.table(FASTASampleCode=rownames(tmp)), subset(seq.df,region=="in"), by="FASTASampleCode" )
		setkey(seq.in.df, dist)
		if(nrow(seq.in.df)<seq.select.n)	cat(paste("\ncan only select less than the requested number of sequences, n=",nrow(seq.in.df)))
		tmp						<- rbind( seq.in.df, data.table(FASTASampleCode=NA, dist=NA, group=c("child","parent1","parent2"), region=NA) )
		seq.in.order			<- tmp[	,	list(n=length(na.omit(FASTASampleCode))) ,by=c("group")]		
		seq.in.order			<- seq.in.order[order(n),]
		overflow				<- 0
		ans						<- data.table(FASTASampleCode=NA, dist=NA, group=NA, region=NA)
		for(x in seq.in.order[,group])
		{			
			#print(x)
			tmp					<- subset(seq.in.df, group==x)
			#print(tmp)
			ans					<- rbind(tmp[seq_len( min(seq.select.n+overflow, nrow(tmp)) ),], ans )
			overflow			<- ifelse(seq.select.n+overflow<nrow(tmp), 0, seq.select.n+overflow-nrow(tmp))
		}
		if(overflow>0)	cat(paste("\nNot as many filler sequences as requested for recombinant region 'in', n=",nrow(seq.out.df)))
		seq.in.df				<- ans[-nrow(ans),]
		if(verbose)	cat(paste("\nSelected balancing sequences for recombinant region 'in', n=",nrow(seq.in.df)))
		if(verbose)	print( seq.in.df[	,	list(n=length(FASTASampleCode)) ,by=c("group","region")] )
		#
		#	select 'seq.select.n' unique filler sequences for 'out' region, balancing by group as much as possible
		#
		if(verbose)	cat(paste("\nSelect sequences for recombinant region 'out'"))
		tmp						<- subset( seq.df, region=='out' )[,FASTASampleCode] 
		tmp						<- hivc.seq.unique(seq.out[tmp,])		
		seq.out.df				<- merge( data.table(FASTASampleCode=rownames(tmp)), subset(seq.df,region=="out"), by="FASTASampleCode" )
		setkey(seq.out.df, dist)
		if(nrow(seq.out.df)<seq.select.n)	cat(paste("\ncan only select less than the requested number of sequences, n=",nrow(seq.out.df)))
		tmp						<- rbind( seq.out.df, data.table(FASTASampleCode=NA, dist=NA, group=c("child","parent1","parent2"), region=NA) )
		seq.out.order			<- tmp[	,	list(n=length(FASTASampleCode)) ,by=c("group")]		
		seq.out.order			<- seq.out.order[order(n),]
		overflow				<- 0
		ans						<- data.table(FASTASampleCode=NA, dist=NA, group=NA, region=NA)		
		for(x in seq.out.order[,group])
		{			
			#print(x)
			tmp					<- subset(seq.out.df, group==x)
			#print(tmp)
			ans					<- rbind(tmp[seq_len( min(seq.select.n+overflow, nrow(tmp)) ),], ans )
			overflow			<- ifelse(seq.select.n+overflow<nrow(tmp), 0, seq.select.n+overflow-nrow(tmp))
		}
		if(overflow>0)	cat(paste("\nNot as many filler sequences as requested for recombinant region 'out', n=",nrow(seq.out.df)))
		seq.out.df				<- ans[-nrow(ans),]
		if(verbose)	cat(paste("\nSelected balancing sequences for recombinant region 'out', n=",nrow(seq.out.df)))
		if(verbose)	print( seq.out.df[	,	list(n=length(FASTASampleCode)) ,by=c("group","region")] )
		#	combine unique filler sequences
		seq.df					<- rbind( seq.in.df, seq.out.df )
		if(verbose)	cat(paste("\nSelected balancing set of closest filler sequences"))
		if(verbose)	print( seq.df[	,	list(n=length(FASTASampleCode)) ,by=c("group","region")] )		
		if( length(unique( seq.df[, FASTASampleCode] ))!=nrow(seq.df) )		stop("Unexpected non-unique sequence names")
		if( nrow(hivc.seq.unique( seq.in[ seq.df[, FASTASampleCode], ] ))!=nrow(seq.df) ) stop("Unexpected non-unique 'in' sequences")
		if( nrow(hivc.seq.unique( seq.out[ seq.df[, FASTASampleCode], ] ))!=nrow(seq.df) ) stop("Unexpected non-unique 'out' sequences")
		#	could be that the triplet sequences are not unique among each other 
		#	if so, make change to one of the triplet sequences
		seq.select				<- c( df.recomb[,parent1],df.recomb[,parent2],df.recomb[,child], seq.df[, FASTASampleCode] )
		tmp						<- hivc.seq.unique( seq.in[ seq.select, ] )
		if( !length(setdiff(c( df.recomb[,parent1],df.recomb[,parent2],df.recomb[,child] ),rownames(tmp))) )
			seq.in				<- tmp
		if( length(setdiff(c( df.recomb[,parent1],df.recomb[,parent2],df.recomb[,child] ),rownames(tmp))) )
		{
			tmp					<- setdiff(c( df.recomb[,parent1],df.recomb[,parent2],df.recomb[,child] ),rownames(tmp) )		#name of sequence in triplet that is identical with one other sequence in triplet
			if(verbose)	cat(paste("\nFound identical triplet sequence for region 'in'", tmp))
			seq.in				<- as.character(seq.in)
			seq.in[tmp,1]		<- ifelse(seq.in[tmp,1]=='t','c',ifelse(seq.in[tmp,1]=='c','t',ifelse(seq.in[tmp,1]=='a','g','a')))
			seq.in				<- as.DNAbin(seq.in)
			tmp					<- hivc.seq.unique( seq.in[ seq.select, ] )
			if( length(setdiff(c( df.recomb[,parent1],df.recomb[,parent2],df.recomb[,child] ),rownames(tmp))) )	stop("Unexpected duplicate for 'in'")
			seq.in				<- tmp
		}			
		seq.out					<- seq.out[seq.select,]
		tmp						<- hivc.seq.unique(seq.out)
		if( !length(setdiff(c( df.recomb[,parent1],df.recomb[,parent2],df.recomb[,child] ),rownames(tmp))) )
			seq.out				<- tmp
		if( length(setdiff(c( df.recomb[,parent1],df.recomb[,parent2],df.recomb[,child] ),rownames(tmp))) )
		{
			tmp					<- setdiff(c( df.recomb[,parent1],df.recomb[,parent2],df.recomb[,child] ),rownames(tmp) )		#name of sequence in triplet that is identical with one other sequence in triplet
			if(verbose)	cat(paste("\nFound identical triplet sequence for region 'out'", tmp))
			seq.out				<- as.character(seq.out)
			seq.out[tmp,1]		<- ifelse(seq.out[tmp,1]=='t','c',ifelse(seq.out[tmp,1]=='c','t',ifelse(seq.out[tmp,1]=='a','g','a')))
			seq.out				<- as.DNAbin(seq.out)
			tmp					<- hivc.seq.unique(seq.out[seq.select,])
			if( length(setdiff(c( df.recomb[,parent1],df.recomb[,parent2],df.recomb[,child] ),rownames(tmp))) )	stop("Unexpected duplicate for 'out'")
			seq.out				<- tmp
		}
		if(any(rownames(seq.in)!=rownames(seq.out)))	stop("Unexpected unequal sequences selected")
		#	reset rownames
		tmp						<- c( paste(c("tparent1","tparent2","tchild"),seq.select[1:3],sep='_'), seq.df[,list(label= paste(region,group,FASTASampleCode,sep='_')), by="FASTASampleCode"][,label] )
		rownames(seq.in)		<- tmp
		rownames(seq.out)		<- tmp
		#	save
		file		<- paste(indir,'/',infile,"_3seqcheck_id",id,"_rIn_",gsub('/',':',insignat),".R",sep='')
		if(verbose) cat(paste("\nsave to ",file))
		save(seq.in, file=file)		
		file		<- paste(indir,'/',infile,"_3seqcheck_id",id,"_rOut_",gsub('/',':',insignat),".R",sep='')
		if(verbose) cat(paste("\nsave to ",file))
		save(seq.out, file=file)
	}
	if(1)
	{
		#
		#	run bootstrap ExaML for region 'in', all boostraps on one processor
		#			
		cmd				<- NULL
		file			<- paste(indir,'/',infile,"_3seqcheck_id",id,"_rIn_examlbs",bs.n,'_',gsub('/',':',insignat),".newick",sep='')		
		if(!resume || !file.exists(file))
		{
			infile.exa	<- paste(infile,"_3seqcheck_id",id,"_rIn",sep='')		
			cmd			<- cmd.examl.bootstrap.on.one.machine(indir, infile.exa, gsub('/',':',insignat),gsub('/',':',insignat), bs.from=bs.from, bs.to=bs.to, bs.n=bs.n, outdir=indir, opt.bootstrap.by="nucleotide", resume=1, verbose=1)
		}
		#
		#	run bootstrap ExaML for region 'out', all boostraps on one processor
		#						
		file			<- paste(indir,'/',infile,"_3seqcheck_id",id,"_rOut_examlbs",bs.n,'_',gsub('/',':',insignat),".newick",sep='')
		if(!resume || !file.exists(file))
		{
			infile.exa	<- paste(infile,"_3seqcheck_id",id,"_rOut",sep='')		
			cmd			<- c(cmd, cmd.examl.bootstrap.on.one.machine(indir, infile.exa, gsub('/',':',insignat),gsub('/',':',insignat), bs.from=bs.from, bs.to=bs.to, bs.n=bs.n, outdir=indir, opt.bootstrap.by="nucleotide", resume=1, verbose=1))
		}
		#
		if(verbose) cat(paste("\ncreated ExaML bootstrap runs, n=",length(cmd)))
		if(!is.null(cmd))
		{
			cmd			<- paste(cmd,collapse='\n')
			#cmd		<- paste(cmd,cmd.recombination.plot.incongruence(indir, infile, gsub('/',':',insignat), triplet.id=id, verbose=1),sep='')				
			#cat(cmd)
			if(verbose) cat(paste("\nqsub ExaML bootstrap runs, hpc.walltime=",hpc.walltime," hpc.mem=",hpc.mem," hpc.nproc=",hpc.nproc," hpc.q=",hpc.q))
			cmd			<- cmd.hpcwrapper(cmd, hpc.walltime=hpc.walltime, hpc.q=hpc.q, hpc.mem=hpc.mem, hpc.nproc=hpc.nproc)
			signat		<- paste(strsplit(date(),split=' ')[[1]],collapse='_',sep='')
			outdir		<- paste(DATA,"tmp",sep='/')
			outfile		<- paste("3sc",signat,sep='.')
			#cat(cmd)			
			cmd.hpccaller(outdir, outfile, cmd)
			Sys.sleep(1)
		}
	}
}



