
#' @export
#' @title Write a phylip file
seq.write.dna.phylip<- function(seq.DNAbin.mat, file)
{		
	tmp<- cbind( rownames(seq.DNAbin.mat), apply( as.character( seq.DNAbin.mat ), 1, function(x) paste(x,sep='',collapse='')  ) )
	tmp<- paste(t(tmp),collapse='\n',sep='')	
	tmp<- paste( paste(c(nrow(seq.DNAbin.mat),ncol(seq.DNAbin.mat)),sep='',collapse=' '),'\n',tmp,'\n',collapse='',sep='' )
	cat(tmp, file=file)
}

#' @export
#' @title Strip gaps
seq.strip.gap<- function(seq.DNAbin.mat, strip.max.len=NA, strip.pc=NA, gap.chars='-')
{
	stopifnot(	!is.na(strip.max.len)&is.na(strip.pc) | is.na(strip.max.len)&!is.na(strip.pc), is.na(strip.pc)|strip.pc<1, is.na(strip.pc)|strip.pc>0, is.na(strip.max.len)|strip.max.len>0)	
	seq.DNAbin.mat	<- as.character(seq.DNAbin.mat)
	if(length(gap.chars)==1)
		tmp			<- apply(seq.DNAbin.mat,2,function(x) length(which(x==gap.chars)))
	if(length(gap.chars)>1)
		tmp			<- apply(seq.DNAbin.mat,2,function(x) length(which(x%in%gap.chars)))
	tmp				<- tmp/nrow(seq.DNAbin.mat)
	if(!is.na(strip.pc))
		return(as.DNAbin(seq.DNAbin.mat[,tmp<strip.pc,drop=FALSE]))
	if(!is.na(strip.max.len))
		for(strip.pc in rev(seq(0.01,1.00,0.01)))
		{
			z		<- which(tmp<strip.pc)
			if(length(z)<=strip.max.len)
			{
				cat('\nStripping gaps at pc=',strip.pc)
				return(as.DNAbin(seq.DNAbin.mat[,z,drop=FALSE]))
			}
		}
}

#' @export
#' @title Find the row numbers in a sequence matrix with a particular sequence pattern
seq.find<- function(char.matrix, pos0= NA, from= c(), verbose=1)
{
	if(is.na(pos0)) 	stop("start position of token to be replaced is missing")
	if(!length(from))	stop("token to be replaced is missing")
	query.colidx	<- seq.int(pos0,pos0+length(from)-1)
	query.yes		<- which( apply(char.matrix, 1, function(x)	all(x[query.colidx]==from) ) )
	query.yes	
}

#' @export 
#' @title Compute a subset of unique sequences
seq.unique<- function(seq.DNAbin.matrix)
{
	x<- as.character(seq.DNAbin.matrix)
	x<- apply(x, 1, function(z) paste(z,collapse=''))
	seq.DNAbin.matrix[!duplicated(x),]			
}

#' @export 
#' @title Perform online BLAST search
#' @description see also package 'RFLPtools'
seq.blast<- function (x, database = "nr", hitListSize = "10", filter = "L", expect = "10", program = "blastn", organism= "HIV-1") 
{
	baseUrl <- "http://www.ncbi.nlm.nih.gov/blast/Blast.cgi"
	query <- paste("QUERY=", as.character(x), "&DATABASE=", database, "&ORGANISM=", organism,
			"&HITLIST_SIZE=", hitListSize, "&FILTER=", filter, "&EXPECT=", 
			expect, "&PROGRAM=", program, sep = "")
	url0 <- sprintf("%s?%s&CMD=Put", baseUrl, query)
	results <- tempfile()
	Sys.sleep(5)
	require(XML)
	post <- htmlTreeParse(url0, useInternalNodes = TRUE)
	x <- post[["string(//comment()[contains(., \"QBlastInfoBegin\")])"]]
	rid <- sub(".*RID = ([[:alnum:]]+).*", "\\1", x)
	rtoe <- as.integer(sub(".*RTOE = ([[:digit:]]+).*", "\\1", 
					x))
	url1 <- sprintf("%s?RID=%s&FORMAT_TYPE=XML&CMD=Get", baseUrl, 
			rid)
	Sys.sleep(rtoe)
	result <- annotate:::.tryParseResult(url1)
	qseq <- xpathApply(result, "//Hsp_qseq", xmlValue)
	hseq <- xpathApply(result, "//Hsp_hseq", xmlValue)
	require(Biostrings)
	res <- list()
	for (i in seq_len(length(qseq))) {
		res[i] <- DNAMultipleAlignment(c(hseq[[i]], qseq[[i]]), 
				rowmask = as(IRanges(), "NormalIRanges"), colmask = as(IRanges(), 
						"NormalIRanges"))
	}
	res
}

#' @export 
#' @title Read BLAST output
#' @description slight modification of read.blast() in pkg RFLPtools; expects blast was run with -outfmt 6
seq.blast.read<- function (file, sep = "\t") 
{
	require(data.table)
	x <- read.table(file = file, header = FALSE, sep = sep, quote = "\"", dec = ".", fill = TRUE, comment.char = "", stringsAsFactors = FALSE)
	if (ncol(x) != 12) 
		stop("Data in given file", basename(file), "has wrong dimension!")
	names(x) <- c("query.id", "subject.id", "identity", "alignment.length","mismatches", "gap.opens", "q.start", "q.end", "s.start","s.end", "evalue", "bit.score")
	data.table(x, key="query.id")
}

#' @export
#' @title Remove drug resistance mutations
seq.rm.drugresistance<- function(seq, outfile=NA)
{
	require(data.table)
	stopifnot( any(rownames(seq)=='HXB2') )		#	expect HXB2 as reference in alignment
	load(system.file(package="big.phylo", 'AC_drugresistance_201508.rda'))	
	load(system.file(package="big.phylo", 'refseq_hiv1_hxb2.rda'))
	hxb2			<- paste(hxb2[, HXB2.K03455], collapse='')
	seq.hxb2		<- paste(as.character(seq['HXB2',]),collapse='')
	#	check that HXB2 in alignment is HXB2
	seq.hxb2.ng		<- gsub('-+','',seq.hxb2)
	stopifnot(nchar(seq.hxb2.ng)<=nchar(hxb2))					#	expect HXB2 in alignment is part of true HXB2
	seq.hxb2.st		<- as.integer(regexpr(substr(seq.hxb2.ng,1,20),hxb2))
	stopifnot(seq.hxb2.st>0)									#	expect HXB2 in alignment is part of true HXB2
	stopifnot(nchar(seq.hxb2.ng)+seq.hxb2.st-1L<=nchar(hxb2))	#	expect HXB2 in alignment is part of true HXB2	
	stopifnot( seq.hxb2.ng==substr(hxb2, seq.hxb2.st, nchar(seq.hxb2.ng)+seq.hxb2.st-1L) )
	#	find coordinates of real HXB2 in alignment
	seq.hxb2.pos	<- c()
	tmp				<- gregexpr('-+',seq.hxb2)[[1]]
	if(tmp[1]>-1)
	{
		seq.hxb2.pos<- data.table(GP_ST= as.integer(tmp), GP_LEN= attr(tmp,'match.length'))
		seq.hxb2.pos<- seq.hxb2.pos[, list(GP_POS= seq.int(GP_ST, length=GP_LEN)), by='GP_ST'][, GP_POS]		
	}
	seq.hxb2.pos 	<- data.table(HXB2INSEQ_POS= setdiff(seq_len(nchar(seq.hxb2)), seq.hxb2.pos))
	seq.hxb2.pos[, HXB2_POS:= seq_len(nrow(seq.hxb2.pos))]
	set(seq.hxb2.pos, NULL, 'HXB2_POS', seq.hxb2.pos[, HXB2_POS+seq.hxb2.st-1L])
	stopifnot( seq.hxb2.pos[, tail(HXB2_POS,1)]<=nchar(hxb2) )
	#	merge with dr
	setnames(dr, 'HXB2.pos', 'HXB2_POS')
	dr				<- merge(dr, seq.hxb2.pos, by='HXB2_POS')
	setnames(dr, 'HXB2INSEQ_POS', 'Alignment.nuc.pos')
	tmp			<- seq.rm.drugresistance.internal(as.character(seq), dr, verbose=1, rtn.DNAbin=1 )
	dr.info		<- tmp$nodr.info
	nodr.seq	<- tmp$nodr.seq
	#
	#	save
	#	
	if(!is.na(outfile))
		save(seq, nodr.seq, dr.info, file=outfile)
	list(nodr.info=dr.info, nodr.seq=nodr.seq)	
}

seq.rm.drugresistance.internal<- function(char.matrix, dr, verbose=1, rtn.DNAbin=1)
{
	if(verbose)	cat(paste("\nchecking for potential drug resistance mutations, n=",nrow(dr)))
	nodr.info	<- dr[, {
							query.yes	<- seq.find(char.matrix, Alignment.nuc.pos, unlist(strsplit(unlist(Mutant.NTs),'')))
							if(length(query.yes))
							{
								ans		<- rownames(char.matrix)[query.yes]
								#print( char.matrix[query.yes, seq.int(dr[i,Alignment.nuc.pos]-3, length.out=9)] ); stop()
							}	
							if(!length(query.yes))
								ans		<- NA_character_				
							list(TAXA=ans)
						}, by=c('HXB2_POS','DR.name','Gene.codon.number','Wild.type','Mutant.NTs','Alignment.nuc.pos')]
	nodr.info	<- subset(nodr.info, !is.na(TAXA))			
	for(i in nodr.info[, sort(unique(Alignment.nuc.pos))])
	{
		tmp		<- subset(nodr.info, Alignment.nuc.pos==i)[, TAXA]
		stopifnot(length(tmp)>0)
		cat(paste('\nsetting at pos',i,'to nnn for taxa, n=', length(tmp)))
		char.matrix[tmp,	seq.int(i, length.out=3) ]<- matrix("n", nrow=length(tmp), ncol=3)	
	}
	if(rtn.DNAbin)
		char.matrix	<- as.DNAbin(char.matrix)	
	list(nodr.seq=char.matrix, nodr.info=nodr.info)	
}	

#' @export
#' @title Compute pairwise distance matrix
seq.dist<- function(seq.DNAbin.matrix, verbose=1)
{
	if(0)
	{
		require(ape)
		ans<- dist.dna(seq.DNAbin.matrix, model="raw", as.matrix=1)
	}
	if(1)
	{
		library(bigmemory)
		options(bigmemory.typecast.warning=FALSE)
		big.matrix.charmax<- 127
		dummy	<- 0
		ans<- big.matrix(nrow(seq.DNAbin.matrix),nrow(seq.DNAbin.matrix), dimnames=list(rownames(seq.DNAbin.matrix),c()), type="char", init=NA)
		if(nrow(seq.DNAbin.matrix)>5)
		{
			ans[1,1:6]<- 2^seq.int(6,11)-1
			if(is.na(ans[1,2]))		stop("unexpected behaviour of bigmemory")
			ans[1,1:6]<- NA
		}						
		for(i1 in seq.int(1,nrow(seq.DNAbin.matrix)-1))
		{			
			seq1<- seq.DNAbin.matrix[i1,]
			time<- system.time	(
					tmp	<- 1 - sapply(seq.int(i1+1,nrow(seq.DNAbin.matrix)),function(i2){		.C("hivc_dist_ambiguous_dna", seq1, seq.DNAbin.matrix[i2,], ncol(seq1), dummy )[[4]]			})
			)[3]
			if(verbose)	cat(paste("\ncompute distance of row",i1,"entries",nrow(seq.DNAbin.matrix)-i1,"took",time))			
			tmp												<- round(tmp*1e3,d=0)			
			tmp[tmp>big.matrix.charmax]						<- big.matrix.charmax
			ans[i1, seq.int(i1+1,nrow(seq.DNAbin.matrix))]	<- tmp
		}		
	}
	ans
}

#' @export
#' @title Replace nucleotide codes in sequence alignment
seq.replace<- function(seq.DNAbin.matrix, code.from='?', code.to='n', verbose=0)
{
	seq.DNAbin.matrix	<- as.character(seq.DNAbin.matrix)		
	seq.DNAbin.matrix	<- apply(seq.DNAbin.matrix, 2, function(col) 		gsub(code.from,code.to,col,fixed=1)			)	
	as.DNAbin( seq.DNAbin.matrix )
}

#' @export
#' @title Remove gaps in sequence alignment
seq.rmgaps<- function(seq.DNAbin.matrix, rm.only.col.gaps=1, verbose=0)
{
	seq.DNAbin.matrix		<- as.character(seq.DNAbin.matrix)		
	if(!rm.only.col.gaps)
	{	
		if(is.matrix(seq.DNAbin.matrix))
		{
			tmp					<- lapply(seq_len(nrow(seq.DNAbin.matrix)), function(i){	seq.DNAbin.matrix[i, seq.DNAbin.matrix[i,]!="-" & seq.DNAbin.matrix[i,]!="?"]	})
			names(tmp)			<- rownames(seq.DNAbin.matrix)
		}
		else
		{
			tmp					<- lapply(seq_along(seq.DNAbin.matrix), function(i){	seq.DNAbin.matrix[[i]][ seq.DNAbin.matrix[[i]]!="-" & seq.DNAbin.matrix[[i]]!="?"]	})
			names(tmp)			<- names(seq.DNAbin.matrix)
		}		
		seq.DNAbin.matrix	<- tmp
	}
	else
	{		
		nogap				<- which( !apply(seq.DNAbin.matrix,2,function(x) all(x=="-" || x=="?")) )
		if(verbose)	cat(paste("\nremove gaps, n=",ncol(seq.DNAbin.matrix)-length(nogap)))
		seq.DNAbin.matrix	<- seq.DNAbin.matrix[,nogap]	
	}
	as.DNAbin( seq.DNAbin.matrix )
}

#' @export
#' @title Read sequences from GenBank
seq.read.GenBank<- function (access.nb, seq.names = access.nb, species.names = TRUE, gene.names = FALSE, as.character = FALSE, attributes= c("origin")) 
{
	require(ape)
	N <- length(access.nb)
	nrequest <- N%/%400 + as.logical(N%%400)
	X <- character(0)
	for (i in 1:nrequest) {
		a <- (i - 1) * 400 + 1
		b <- 400 * i
		if (i == nrequest) 
			b <- N
		URL <- paste("http://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nucleotide&id=", 
				paste(access.nb[a:b], collapse = ","), "&rettype=gb&retmode=text", 
				sep = "")
		X <- c(X, scan(file = URL, what = "", sep = "\n", quiet = TRUE))
	}
	
	FI <- grep("^ {0,}ORIGIN", X) + 1
	LA <- which(X == "//") - 1
	obj <- vector("list", N)
	for (i in 1:N) {
		tmp <- gsub("[[:digit:] ]", "", X[FI[i]:LA[i]])
		obj[[i]] <- unlist(strsplit(tmp, NULL))
	}
	names(obj) <- seq.names
	if (!as.character) 
		obj <- as.DNAbin(obj)
	if(length(attributes)) 
	{
		attr.lines	<- lapply(attributes, function(attr) grep(paste("^ {0,}/",attr,sep=''), X) 	)
		attr.fields	<- lapply(seq_along(attr.lines),function(i)		gsub("\"$","",gsub(paste("^ {0,}/",attributes[i],"=\"",sep=''),"",X[attr.lines[[i]]]))		)						
		for(i in seq_along(attributes))
			attr(obj, attributes[[i]])<- attr.fields[[i]]
	}
	if (gene.names) {
		tmp <- character(N)
		sp <- grep(" +gene +<", X)
		for (i in 1:N) tmp[i] <- unlist(strsplit(X[sp[i + 1L]], 
							" +/gene=\""))[2]
		attr(obj, "gene") <- gsub("\"$", "", tmp)
	}
	obj
}

#' @export
seq.create.referencepairs<- function(dir.name= DATA)
{
	if(0)	#generate ATHENA_2013_hptn052.rda
	{
		gban				<- c( paste("JN",seq.int(247047,247075),sep=''),paste("JN",seq.int(634296,634492),sep='') )	
		tmp					<- hivc.read.GenBank(gban, as.character=0, attributes= c("isolate","country","collection_date"))		
		file				<- "ATHENA_2013_hptn052.fa"
		write.dna(tmp, paste(DATA,file,sep='/'), format = "fasta" )
		cmd					<- hiv.cmd.clustalo(paste(dir.name,"tmp",sep='/'), file, signat='', outdir=paste(dir.name,"tmp",sep='/'))
		cat(cmd)
		if(0) system(cmd)		
		file				<- paste(DATA,"tmp/ATHENA_2013_hptn052.fa.clustalo",sep='/')
		hptn052				<- read.dna(file, format = "fasta" )		
		rownames(hptn052)	<- attr(tmp,"isolate")		
		file				<- paste(DATA,"tmp/ATHENA_2013_hptn052.rda",sep='/')
		save(hptn052, file= file)
	}
	file<- paste(DATA,"tmp/ATHENA_2013_hptn052.rda",sep='/')		
	load(file)
	#select control sequences and compute pairwise distances between them
	hptn052.control		<- grep("control", rownames(hptn052))	
	hptn052.control.d	<- hivc.pwdist( hptn052[hptn052.control,] )
	hptn052.control.d	<- hptn052.control.d[ upper.tri(hptn052.control.d) ]
	#select positive sequences and compute pairwise distance between them	
	hptn052.sdc			<- c(grep("A", rownames(hptn052)), grep("B", rownames(hptn052)))
	hptn052.sdc			<- hptn052[hptn052.sdc,]
	hptn052.sdc.ind		<- sapply(strsplit(rownames(hptn052.sdc),'-'), function(x)
			{ 
				as.numeric(x[2]) + ifelse(x[3]=='I',0,0.5)
			})
	hptn052.sdc.ind		<- sapply( unique( hptn052.sdc.ind ), function(x){		hptn052.sdc[hptn052.sdc.ind==x,]	})
	hptn052.sdc.ind		<- lapply( which( sapply(hptn052.sdc.ind,nrow)==2 ), function(i)  hptn052.sdc.ind[[i]] )
	hptn052.sdc.ind.d	<- sapply( hptn052.sdc.ind, function(x) hivc.pwdist(x)[1,2] )
	
	list( gd.islink= hptn052.sdc.ind.d, gd.unlinked= hptn052.control.d )
}

#' @export
#' @title Write a NEXUS file
seq.write.dna.nexus<- function(seq.DNAbin.mat, ph=NULL, file=NULL, nexus.format="DNA",nexus.gap='-', nexus.missing='?', nexus.interleave="NO")
{		
	tmp		<- cbind( rownames(seq.DNAbin.mat), apply( as.character( seq.DNAbin.mat ), 1, function(x) paste(x,sep='',collapse='')  ) )
	tmp		<- apply(tmp, 1, function(x) paste(x, collapse='\t', sep=''))
	tmp		<- paste(tmp, collapse='\n',sep='')
	header	<- paste( "#NEXUS\nBEGIN DATA;\nDIMENSIONS NTAX=",nrow(seq.DNAbin.mat)," NCHAR=",ncol(seq.DNAbin.mat),";\nFORMAT DATATYPE=",nexus.format," MISSING=",nexus.missing," GAP=",nexus.gap," INTERLEAVE=",nexus.interleave,";\nMATRIX\n", collapse='',sep='')
	tmp		<- paste(header, tmp, "\n;\nEND;\n", sep='')
	if(!is.null(ph))
	{
		tmp	<- paste(tmp,"#BEGIN TAXA;\nTAXLABELS ", paste(ph$tip.label, sep='',collapse=' '), ';\nEND;\n\n',sep='')
		tmp	<- paste(tmp, "BEGIN TREES;\nTREE tree1 = ", write.tree(ph), "\nEND;\n", sep='')
	}
	if(!is.null(file))
		cat(tmp, file=file)
	tmp
}		

#' @export
#' @title Compute the length of sequences
seq.length<- function(seq.DNAbin.mat, exclude=c('-','?'))
{
	counts	<- apply(seq.DNAbin.mat,1,function(x) base.freq(x, freq=1, all=1))
	apply(counts[ !rownames(counts)%in%exclude, ],2,sum)
}

#' @export
#' @title Compute the proportion of ambiguous nucleotides
seq.proportion.ambiguous<- function(seq.DNAbin.mat, exclude=c('-','?'))
{
	counts	<- apply(seq.DNAbin.mat,1,function(x) base.freq(x, freq=1, all=1))
	len		<- apply(counts[ !rownames(counts)%in%exclude, ],2,sum)
	pa		<- apply(counts[c("r", "m", "w", "s", "k", "y", "v", "h", "d", "b"),],2,sum)
	pa/len
}

#' @export
#' @title Compute the GC content
seq.gc.content<- function(seq.DNAbin.mat)
{	
	rna.gc.fraction.n		<- c('a','c','g','t',	'r','m','w','s',	'k','y','v','h',		'd','b','n','-','?')
	rna.gc.fraction			<- c( 0, 1, 1, 0,		0.5, 0.5, 0, 1, 	1/2, 1/2, 2/3, 1/3,		1/3,2/3, 1/4, 0, 0)		#this fraction assumes that U is synonymous with T and that U does not occur in the code
	rna.gc.sum				<- c( T, T, T, T,       T, T, T, T,         T, T, T, T,				T, T, T, F, F )	
	counts					<- apply(seq.DNAbin.mat,1,function(x) base.freq(x, freq=1, all=1))
	apply(counts*rna.gc.fraction,2,sum) / apply(counts[rna.gc.sum,],2,sum)
}

######################################################################################
hivc.clu.polyphyletic.clusters<- function(cluphy.df, cluphy.subtrees=NULL, ph=NULL, clustering=NULL, verbose=1, plot.file=NA, pdf.scaley=25, pdf.xlim=NULL, cex.nodelabel=0.2, cex.tiplabel=0.2, adj.tiplabel= c(-0.15,0.5))
{
	if(!is.null(ph) && !is.null(clustering))
	{
		#get node corresponding to index of selected clusters		
		cluphy.cluidx				<- clustering[["clu.idx"]][ unique( cluphy.df[,cluster] ) ]		
		if(any(is.na(cluphy.cluidx)))	stop("unexcpected NA in cluphy.cluidx")
		#get selected clusters into single phylogeny
		cluphy.subtrees				<- lapply(cluphy.cluidx, function(x)		extract.clade(ph, x, root.edge= 1, interactive = FALSE) 		)
		names(cluphy.subtrees)		<- unique( cluphy.df[,cluster] )
	}
	else if(is.null(cluphy.subtrees))
		stop("expect either ('ph' and 'clustering') or 'cluphy.subtrees' as input")
	cluphy							<- eval(parse(text=paste('cluphy.subtrees[[',seq_along(cluphy.subtrees),']]', sep='',collapse='+')))
	print(cluphy)
	plot.coordinates				<- NULL
	#plot selected clusters
	if(!is.na(plot.file))
	{		
		cluphy.tiplabels						<- hivc.clu.get.tiplabels( cluphy, copy(cluphy.df) )
		if(verbose) cat(paste("\nwrite tree to file",plot.file))
		color.except.rootedge					<- rep(1, Nnode(cluphy, internal.only=F))
		color.except.rootedge[Ntip(cluphy)+1]	<- NA
		plot.coordinates						<- hivc.clu.plot(cluphy, color.except.rootedge, file=plot.file, pdf.scaley=pdf.scaley, pdf.off=0, pdf.xlim=pdf.xlim, cex.nodelabel=cex.nodelabel )		
		hivc.clu.plot.tiplabels( seq_len(Ntip(cluphy)), cluphy.tiplabels$text, cluphy.tiplabels$col, cex=cex.tiplabel, adj=adj.tiplabel, add.xinch=0, add.yinch=0 )
		dev.off()
	}	
	list(cluphy=cluphy, cluphy.subtrees=cluphy.subtrees, plot.coordinates=plot.coordinates)
}	
######################################################################################
hivc.clu.get.tiplabels<- function(ph, 	df.info, col.notmsm="#4EB3D3", col.Early="#EF9708", col.highVL="#FEE391", col.AfterTreat="#D4B9DA", col.green="#D9F0A3", col.latePres="#FA9FB5",
										select=c("CountryInfection","Trm","Sex","isAcute","lRNA.early","NegT","AnyPos_T1","PosSeqT","lRNA.hb4tr_LT","lRNA_bTS","lRNA_TS","lRNA_aTS","lRNAi_bTS","lRNAi_aTS","AnyT_T1","TrImo_bTS","TrImo_aTS","PosCD4_T1","CD4_T1","CD4_bTS","CD4_TS","CD4_aTS","Patient","RegionHospital") )
{
	require(colorspace)
	require(RColorBrewer)
	#
	#	PATIENT
	#
	#set colors CountryInfection
	tmp					<- rep("transparent",nrow(df.info))
	df.info[,CountryInfection.col:=tmp]
	set(df.info, which(df.info[,!is.na(CountryInfection) & CountryInfection!="NL"]), "CountryInfection.col", col.notmsm)
	#set colors Patient
	tmp					<- unique( df.info[,Patient] )
	#tmp2				<- diverge_hcl(length(tmp), h = c(246, 40), c = 96, l = c(85, 90))		
	tmp					<- data.table(Patient=tmp, Patient.col="transparent", key="Patient")
	df.info				<- merge(  df.info, tmp, all.x=1, by="Patient" )
	#set colors Sex
	df.info[,Sex.col:="transparent"]						
	set(df.info, which(df.info[,!is.na(Sex) & Sex!="M"]), "Sex.col", col.notmsm)
	#set colors MSM or BI	
	df.info[,Trm.col:="transparent"]						
	set(df.info, which(df.info[,!is.na(Trm) & Trm!="MSM"]), "Trm.col", col.notmsm)
	#set colors RegionHospital
	tmp					<- levels( df.info[,RegionHospital] )		
	tmp2				<- brewer.pal(length(tmp), "Dark2")
	tmp					<- data.table(RegionHospital=tmp, RegionHospital.col=tmp2, key="RegionHospital")
	df.info				<- merge(  df.info, tmp, all.x=1, by="RegionHospital" )
	#set colors isAcute
	df.info[,isAcute.col:="transparent"]						
	set(df.info, which(df.info[,isAcute%in%c("Yes","Maybe")]), "isAcute.col", col.Early)
	#
	#	TIMELINE
	#
	#	set transparent colors:		PosSeqT NegT AnyPos_T1	
	df.info[, NegT.col:="transparent"]
	df.info[, AnyPos_T1.col:="transparent"]
	#
	#	TREATMENT
	#
	#	set AnyT_T1.col				if 	PosSeqT<AnyT_T1 	orange 		else pink
	df.info[, AnyT_T1.col:="transparent"]
	select.seqb4tr		<- which(df.info[, !is.na(PosSeqT) & !is.na(AnyT_T1) & PosSeqT<=AnyT_T1])
	select.seqatr		<- which(df.info[, !is.na(PosSeqT) & !is.na(AnyT_T1) & PosSeqT>AnyT_T1])
	set(df.info, select.seqb4tr, "AnyT_T1.col", col.Early)
	set(df.info, select.seqatr, "AnyT_T1.col", col.AfterTreat)
	df.info[, TrImo_bTS.col:="transparent"]
	df.info[, TrImo_aTS.col:="transparent"]
	set(df.info, which(df.info[, TrImo_aTS>4]), "TrImo_aTS.col", col.AfterTreat)
	set(df.info, which(df.info[, TrImo_bTS>4]), "TrImo_bTS.col", col.AfterTreat)
	df.info[, PosSeqT.col:="transparent"]
	set(df.info, select.seqb4tr, "PosSeqT.col", col.Early)
	set(df.info, select.seqatr, "PosSeqT.col", col.AfterTreat)	
	#
	#	VIRAL LOAD
	#	
	df.info[, lRNA_bTS.col:="transparent"]
	set(df.info,	which(df.info[,lRNA_bTS>5]),	"lRNA_bTS.col",		col.highVL) 
	df.info[, lRNA_TS.col:="transparent"]
	set(df.info,	which(df.info[,lRNA_TS>5]),	"lRNA_TS.col",		col.highVL) 
	df.info[, lRNA_aTS.col:="transparent"]
	set(df.info,	which(df.info[,lRNA_aTS>3.5]),	"lRNA_aTS.col",		col.highVL) 
	df.info[, lRNAi_bTS.col:="transparent"]
	set(df.info,	which(df.info[,lRNAi_bTS>0.75]),	"lRNAi_bTS.col",		col.highVL) 
	df.info[, lRNAi_aTS.col:="transparent"]
	set(df.info,	which(df.info[,lRNAi_aTS>0.25]),	"lRNAi_aTS.col",		col.highVL) 
	df.info[, lRNA.hb4tr_LT.col:="transparent"]
	set(df.info,	which(df.info[,!is.na(lRNA.hb4tr_LT) & PosSeqT<=lRNA.hb4tr_LT]),	"lRNA.hb4tr_LT.col",		col.highVL) 
	df.info[, lRNA.early.col:="transparent"]
	set(df.info,	which(df.info[, lRNA.early]),	"lRNA.early.col",		col.Early) 
	#
	#	CD4
	#		
	df.info[, CD4_T1.col:="transparent"]
	set(df.info,	which(df.info[,CD4_T1>350]),	"CD4_T1.col",		col.green)	
	df.info[, PosCD4_T1.col:="transparent"]
	set(df.info,	which(df.info[,CD4_T1>350]),	"PosCD4_T1.col",		col.green)	
	df.info[, CD4_bTS.col:="transparent"]
	set(df.info,	which(df.info[,CD4_bTS>350]),	"CD4_bTS.col",		col.green)	
	df.info[, CD4_TS.col:="transparent"]
	set(df.info,	which(df.info[,CD4_TS>350]),	"CD4_TS.col",		col.green)	
	df.info[, CD4_aTS.col:="transparent"]
	set(df.info,	which(df.info[,CD4_aTS>350]),	"CD4_aTS.col",		col.green)
	#
	#	color late presenter
	#
	tmp<- which(df.info[,CD4_T1<350])
	set(df.info,	tmp,	"CD4_T1.col",			col.latePres)
	set(df.info,	tmp,	"PosCD4_T1.col",		col.latePres)
	set(df.info,	tmp,	"lRNA_bTS.col",			col.latePres)
	set(df.info,	tmp,	"lRNAi_bTS.col",		col.latePres)
	#	
	#	convert time to string	& handle inaccurate NegT or AnyPos_T1
	#
	tmp		<- which(df.info[, as.POSIXlt(NegT)$mday==1 & as.POSIXlt(NegT)$mon==0])				#since NegT were reset, these are the ones with inaccurate month
	set(df.info,NULL,	"NegT", 		as.character( df.info[,NegT], "%y.%m" ))
	set(df.info, tmp,	"NegT", 		paste(df.info[tmp,substr(NegT,1,3)],'??',sep='') )
	tmp		<- which(df.info[, as.POSIXlt(AnyPos_T1)$mday==31 & as.POSIXlt(AnyPos_T1)$mon==11])	#since AnyPos_T1 were reset, these are the ones with inaccurate month
	set(df.info,NULL,	"AnyPos_T1",	as.character( df.info[,AnyPos_T1], "%y.%m" ))
	set(df.info, tmp, 	"AnyPos_T1", 	paste(df.info[tmp,substr(AnyPos_T1,1,3)],'??',sep='') )
	set(df.info,NULL,	"PosSeqT", 		as.character( df.info[,PosSeqT], "%y.%m" ))
	set(df.info,NULL,	"lRNA.hb4tr_LT",as.character( df.info[,lRNA.hb4tr_LT], "%y.%m" ))
	set(df.info,NULL,	"PosCD4_T1", 	as.character( df.info[,PosCD4_T1], "%y.%m" ))
	set(df.info,NULL,	"AnyT_T1", 		as.character( df.info[,AnyT_T1], "%y.%m" ))
	#	
	#	set isAcute to either Y or M or N
	set(df.info,NULL,"isAcute", 		as.character( df.info[,isAcute]))
	set(df.info,which(df.info[,isAcute=="No"]),"isAcute",'N')
	set(df.info,which(df.info[,isAcute=="Maybe"]),"isAcute",'M')
	set(df.info,which(df.info[,isAcute=="Yes"]),"isAcute",'Y')
	#	set lRNA.early to either HVLE or ----
	set(df.info,NULL,"lRNA.early", 		as.character( df.info[, lRNA.early]))
	set(df.info,which(df.info[,lRNA.early=="TRUE"]),"lRNA.early","Y")
	set(df.info,which(df.info[,lRNA.early=="FALSE"]),"lRNA.early","N")
	#	set lRNAi_bTS lRNAi_aTS
	set(df.info,NULL,"lRNAi_bTS", 		as.character( round(df.info[, lRNAi_bTS],d=2)))
	set(df.info,NULL,"lRNAi_aTS", 		as.character( round(df.info[, lRNAi_aTS],d=2)))
	#	set TrImo_bTS TrImo_aTS
	set(df.info,NULL,"TrImo_bTS", 		as.character( round(df.info[, TrImo_bTS],d=1)))
	set(df.info,NULL,"TrImo_aTS", 		as.character( round(df.info[, TrImo_aTS],d=1)))
	#	set CD4_T1
	set(df.info,NULL,"CD4_T1", 		as.character( round(df.info[, CD4_T1],d=0)))
	#	set CD4_TS CD4_bTS CD4_aTS
	set(df.info,NULL,"CD4_TS", 		as.character( df.info[, CD4_TS]))
	set(df.info,NULL,"CD4_bTS", 		as.character( df.info[, CD4_bTS]))
	set(df.info,NULL,"CD4_aTS", 		as.character( df.info[, CD4_aTS]))
	#	set 'Amst' to 'A'
	set(df.info,NULL,"RegionHospital", 		as.character( df.info[,RegionHospital]))
	set(df.info,which(df.info[,RegionHospital=="Amst"]),"RegionHospital",'A')	
	#
	#	handle missing entries -- ensure that alignment is OK
	#
	setkey(df.info, Patient)
	tmp					<- which(is.na(df.info[,Patient]))
	set(df.info, tmp, "Patient", '')
	set(df.info, tmp, "Patient.col", "transparent")
	tmp					<- which(is.na(df.info[,RegionHospital]))
	set(df.info, tmp, "RegionHospital", '-')
	set(df.info, tmp, "RegionHospital.col", "transparent")		
	tmp					<- which(is.na(df.info[,CountryInfection]))
	set(df.info, tmp, "CountryInfection", "--")
	set(df.info, tmp, "CountryInfection.col", "transparent")				
	tmp					<- which(is.na(df.info[,Trm]))
	set(df.info, tmp, "Trm", '--')
	set(df.info, tmp, "Trm.col", "transparent")
	tmp					<- which(is.na(df.info[,Sex]))
	set(df.info, tmp, "Sex", '-')
	set(df.info, tmp, "Sex.col", "transparent")	
	set(df.info, which(df.info[,is.na(isAcute)]), "isAcute", 	'-')	
	tmp					<- which(is.na(df.info[,AnyPos_T1]))
	set(df.info, tmp, "AnyPos_T1", "--.--")	
	set(df.info, tmp, "AnyPos_T1.col", "transparent")
	tmp					<- which(is.na(df.info[,NegT]))
	set(df.info, tmp, "NegT", "--.--")	
	set(df.info, tmp, "NegT.col", "transparent")
	tmp					<- which(is.na(df.info[,PosSeqT]))
	set(df.info, tmp, "PosSeqT", "--.--")	
	set(df.info, tmp, "PosSeqT.col", "transparent")		
	tmp					<- which(is.na(df.info[,lRNA.hb4tr_LT]))
	set(df.info, tmp, "lRNA.hb4tr_LT", "--.--")		
	tmp					<- which(is.na(df.info[,PosCD4_T1]))
	set(df.info, tmp, "PosCD4_T1", "--.--")		
	tmp					<- which(is.na(df.info[,AnyT_T1]))
	set(df.info, tmp, "AnyT_T1", "--.--")		
	set(df.info, which(df.info[,is.na(TrImo_aTS)]), "TrImo_aTS", 	'-')
	set(df.info, which(df.info[,is.na(TrImo_bTS)]), "TrImo_bTS", 	'-')	
	set(df.info, which(df.info[,is.na(CD4_T1)]), "CD4_T1", 			'---')
	set(df.info, which(df.info[,is.na(PosCD4_T1)]), "PosCD4_T1", 	'---')
	set(df.info, which(df.info[,is.na(CD4_TS)]), "CD4_TS", 			'---')
	set(df.info, which(df.info[,is.na(CD4_bTS)]), "CD4_bTS", 		'---')
	set(df.info, which(df.info[,is.na(CD4_aTS)]), "CD4_aTS", 		'---')
	#
	#	add suffixes 
	#	
	set(df.info,NULL,"AnyPos_T1",		paste("d:",df.info[,AnyPos_T1],sep=''))
	set(df.info,NULL,"NegT", 			paste("n:",df.info[,NegT],sep=''))
	set(df.info,NULL,"PosSeqT", 		paste("s:",df.info[,PosSeqT],"   ",sep=''))
	set(df.info,NULL,"lRNA.hb4tr_LT",	paste("VLeh:",df.info[,lRNA.hb4tr_LT],sep=''))
	set(df.info,NULL,"AnyT_T1", 		paste("TR+:",df.info[,AnyT_T1],sep=''))
	set(df.info,NULL,"TrImo_bTS", 		paste("TR-:",df.info[,TrImo_bTS],sep=''))
	set(df.info,NULL,"TrImo_aTS", 		paste(":",df.info[,TrImo_aTS],"   ",sep=''))			
	set(df.info,NULL,"isAcute", 		paste("AC:",df.info[,isAcute],sep=''))
	set(df.info,NULL,"lRNA.early", 		paste(":",df.info[,lRNA.early],"   ",sep=''))
	set(df.info,NULL,"CountryInfection",paste("inf:",df.info[,CountryInfection],sep=''))
	set(df.info,NULL,"RegionHospital",	paste("",df.info[,RegionHospital],sep=''))
	set(df.info,NULL,"lRNA_bTS", 		paste("VL:",df.info[,lRNA_bTS],sep=''))
	set(df.info,NULL,"lRNA_TS", 		paste(":",df.info[,lRNA_TS],sep=''))
	set(df.info,NULL,"lRNA_aTS", 		paste(":",df.info[,lRNA_aTS],sep=''))
	set(df.info,NULL,"lRNAi_bTS", 		paste("VI:",df.info[,lRNAi_bTS],sep=''))
	set(df.info,NULL,"lRNAi_aTS", 		paste(":",df.info[,lRNAi_aTS],"   ",sep=''))	
	set(df.info,NULL,"PosCD4_T1", 		paste("CD4:",df.info[,PosCD4_T1],sep=''))
	set(df.info,NULL,"CD4_T1", 			paste(":",df.info[,CD4_T1],sep=''))	
	set(df.info,NULL,"CD4_bTS", 		paste(":",df.info[,CD4_bTS],sep=''))
	set(df.info,NULL,"CD4_TS", 			paste(":",df.info[,CD4_TS],sep=''))
	set(df.info,NULL,"CD4_aTS", 		paste(":",df.info[,CD4_aTS],"   ",sep=''))
	#
	#get df.info into order of tips
	#
	setkey(df.info, FASTASampleCode)
	df.info				<- df.info[ph$tip.label,]
	#select text and col matrix 
	text				<- t( as.matrix( df.info[,select, with=0] ) )
	colnames(text)		<- ph$tip.label
	col					<- t( as.matrix( df.info[,paste(select,".col",sep=''),with=0] ) )
	colnames(col)		<- ph$tip.label
	ans<- list(text=text, col=col)
	ans
}
######################################################################################
#prepare standard format of tip labels -- requires df.info to be sorted along the tips as they appear in a phylogeny
hivc.clu.get.tiplabels.v1<- function(ph, df.info)
{
	require(colorspace)
	require(RColorBrewer)
	#set colors CountryInfection
	tmp					<- rep("transparent",nrow(df.info))
	df.info[,CountryInfection.col:=tmp]
	set(df.info, which(df.info[,CountryInfection=="NL"]), "CountryInfection.col", "#EF9708")
	#set colors Patient
	tmp					<- unique( df.info[,Patient] )
	tmp2				<- diverge_hcl(length(tmp), h = c(246, 40), c = 96, l = c(85, 90))		
	tmp					<- data.table(Patient=tmp, Patient.col=tmp2, key="Patient")
	df.info				<- merge(  df.info, tmp, all.x=1, by="Patient" )
	#set colors Sex
	tmp					<- rep("transparent",nrow(df.info))
	df.info[,Sex.col:=tmp]						
	set(df.info, which(df.info[,Sex=="M"]), "Sex.col", "#EF9708")
	#set colors MSM or BI
	tmp					<- rep("transparent",nrow(df.info))
	df.info[,Trm.col:=tmp]						
	set(df.info, which(df.info[,Trm%in%c("MSM","BI")]), "Trm.col", "#EF9708")
	#set colors RegionHospital
	tmp					<- levels( df.info[,RegionHospital] )		
	tmp2				<- brewer.pal(length(tmp), "Dark2")
	tmp					<- data.table(RegionHospital=tmp, RegionHospital.col=tmp2, key="RegionHospital")
	df.info				<- merge(  df.info, tmp, all.x=1, by="RegionHospital" )		
	#set colors time		
	tmp					<- range( range(df.info[, NegT],na.rm=1), range(df.info[, AnyPos_T1],na.rm=1) )
	tmp					<- as.POSIXlt( seq.Date(tmp[1],tmp[2]+365,by="years") )$year
	tmp2				<- heat_hcl(length(tmp), h = c(0, -100), l = c(75, 40), c = c(40, 80), power = 1)
	yearcols			<- data.table(Year=tmp, Year.col=tmp2)
	#set colors PosSeqT
	tmp					<- data.table(PosSeqT= unique( df.info[,PosSeqT] ), key="PosSeqT" )
	tmp2				<- tmp[, as.POSIXlt(PosSeqT)$year]
	tmp[,"Year":=tmp2]
	tmp					<- subset( merge(tmp, yearcols, all.x=1, by="Year"), select=c(PosSeqT,Year.col))
	setnames(tmp,"Year.col","PosSeqT.col")
	df.info			<- merge(  df.info, tmp, all.x=1, by="PosSeqT" )
	#set colors NegT
	tmp					<- data.table(NegT= unique( df.info[,NegT] ), key="NegT" )
	tmp2				<- tmp[, as.POSIXlt(NegT)$year]
	tmp[,"Year":=tmp2]
	tmp					<- subset( merge(tmp, yearcols, all.x=1, by="Year"), select=c(NegT,Year.col))
	setnames(tmp,"Year.col","NegT.col")
	df.info			<- merge(  df.info, tmp, all.x=1, by="NegT" )
	#set colors AnyPos_T1
	tmp					<- data.table(AnyPos_T1= unique( df.info[,AnyPos_T1] ), key="AnyPos_T1" )
	tmp2				<- tmp[, as.POSIXlt(AnyPos_T1)$year]
	tmp[,"Year":=tmp2]
	tmp					<- subset( merge(tmp, yearcols, all.x=1, by="Year"), select=c(AnyPos_T1,Year.col))
	setnames(tmp,"Year.col","AnyPos_T1.col")
	df.info			<- merge(  df.info, tmp, all.x=1, by="AnyPos_T1" )									
	
	#convert time to string
	set(df.info,NULL,"AnyPos_T1",substr(as.character( df.info[,AnyPos_T1] ),1,7))
	set(df.info,NULL,"NegT",substr(as.character( df.info[,NegT] ),1,7))
	set(df.info,NULL,"PosSeqT",substr(as.character( df.info[,PosSeqT] ),1,7))
	
	#handle missing entries
	setkey(df.info, Patient)
	tmp					<- which(is.na(df.info[,Patient]))
	set(df.info, tmp, "Patient", '')
	set(df.info, tmp, "Patient.col", "transparent")
	tmp					<- which(is.na(df.info[,RegionHospital]))
	set(df.info, tmp, "RegionHospital", '-')
	set(df.info, tmp, "RegionHospital.col", "transparent")		
	tmp					<- which(is.na(df.info[,CountryInfection]))
	set(df.info, tmp, "CountryInfection", "--")
	set(df.info, tmp, "CountryInfection.col", "transparent")				
	tmp					<- which(is.na(df.info[,Trm]))
	set(df.info, tmp, "Trm", '')
	set(df.info, tmp, "Trm.col", "transparent")
	tmp					<- which(is.na(df.info[,Sex]))
	set(df.info, tmp, "Sex", '')
	set(df.info, tmp, "Sex.col", "transparent")
	tmp					<- which(is.na(df.info[,AnyPos_T1]))
	set(df.info, tmp, "AnyPos_T1", "-------")
	set(df.info, NULL, "AnyPos_T1", paste("HIV+:",df.info[,AnyPos_T1],sep=''))
	set(df.info, tmp, "AnyPos_T1.col", "transparent")
	tmp					<- which(is.na(df.info[,NegT]))
	set(df.info, tmp, "NegT", "-------")
	set(df.info, NULL, "NegT", paste("HIV-:",df.info[,NegT],sep=''))
	set(df.info, tmp, "NegT.col", "transparent")
	tmp					<- which(is.na(df.info[,PosSeqT]))
	set(df.info, tmp, "PosSeqT", "-------")
	set(df.info, NULL, "PosSeqT", paste("HIVS:",df.info[,PosSeqT],sep=''))
	set(df.info, tmp, "PosSeqT.col", "transparent")		
	#get df.info into order of tips
	setkey(df.info, FASTASampleCode)
	df.info				<- df.info[ph$tip.label,]
	#select text and col matrix 
	text				<- t( as.matrix( subset(df.info,select=c(CountryInfection, Trm, Sex, NegT, AnyPos_T1, PosSeqT, Patient, RegionHospital)) ) )
	colnames(text)		<- ph$tip.label
	col					<- t( as.matrix( subset(df.info,select=c(CountryInfection.col, Trm.col, Sex.col, NegT.col, AnyPos_T1.col, PosSeqT.col, Patient.col, RegionHospital.col)) ) )
	colnames(col)		<- ph$tip.label
	ans<- list(text=text, col=col)
	ans
}	
######################################################################################
hivc.clu.plot.tiplabels<- function (tip, text, col, xx=NULL, adj = c(-0.05, 0.5), cex=1, add.xinch= 0.03, add.yinch= 0.02) 
{		
	lastPP 			<- get("last_plot.phylo", envir = .PlotPhyloEnv)
	if(is.null(xx))
		xx 			<- lastPP$xx[tip]
	yy 				<- lastPP$yy[tip]
	if(length(tip)==1)
	{
		#assume vector of text and col of same length
		if(!is.vector(text) || !is.vector(col))	stop("expect vector text and vector col")
		if(length(text)!=length(col))			stop("expect same length of text and col")
		wh				<- sapply(text, function(x){	 c(xinch(strwidth(x, units = "inches", cex = cex)), yinch(strheight(x, units = "inches", cex = cex))) 		})
		wh.total		<- c(sum(wh[1,]),max(wh[2,]))
		coord			<- matrix(NA,6,length(text),dimnames=list(c("xl","xr","yb","yt","xx","yy"),c()))
		coord["xl",]	<- xx - wh.total[1] * adj[1] - xinch(add.xinch)
		tmp				<- cumsum(wh[1,]) + xinch(add.xinch) / ncol(wh)
		coord["xl",]	<- coord["xl",] + c(0, tmp[-ncol(coord)])
		coord["xr",]	<- coord["xl",] + wh[1,] + xinch(add.xinch)/ ncol(wh)
		coord["yb",]	<- yy - wh.total[2] * adj[2] - yinch(add.yinch)
		coord["yt",]	<- coord["yb",] + wh.total[2] + yinch(add.yinch)
		coord["xx",]	<- coord["xl",] + c( diff( coord[1,] ), diff(coord[1:2,ncol(wh)])  ) / 2
		coord["yy",]	<- rep(yy, ncol(coord))
#print(wh); print(wh.total); print(coord)		
	}
	else
	{
		#assume matrix of text and col of same ncol
		if(!is.matrix(text) || !is.matrix(col))				stop("expect matrix text and vector col")
		if(ncol(text)!=ncol(col) || nrow(text)!=nrow(col))	stop("expect same dimensions of text and col")
		if(ncol(text)!=length(tip))							stop("expect cols of text and col to correspond to tips")
		coord	<- lapply( seq_along(xx),function(i)
				{
					wh				<- sapply(text[,i], function(x){	 c(xinch(strwidth(x, units = "inches", cex = cex)), yinch(strheight(x, units = "inches", cex = cex))) 		})
					wh.total		<- c(sum(wh[1,]),max(wh[2,]))
					coord			<- matrix(NA,6,nrow(text),dimnames=list(c("xl","xr","yb","yt","xx","yy"),c()))
					coord["xl",]	<- xx[i] - wh.total[1] * adj[1] - xinch(add.xinch)
					tmp				<- cumsum(wh[1,]) + xinch(add.xinch) / ncol(wh)
					coord["xl",]	<- coord["xl",] + c(0, tmp[-ncol(coord)])
					coord["xr",]	<- coord["xl",] + wh[1,] + xinch(add.xinch)/ ncol(wh)
					coord["yb",]	<- yy[i] - wh.total[2] * adj[2] - yinch(add.yinch)
					coord["yt",]	<- coord["yb",] + wh.total[2] + yinch(add.yinch)
					coord["xx",]	<- coord["xl",] + c( diff( coord[1,] ), diff(coord[1:2,ncol(wh)])  ) / 2
					coord["yy",]	<- rep(yy[i], ncol(coord))
					coord
				})
		coord	<- do.call("cbind",coord)
		text	<- as.vector(text)
		col		<- as.vector(col)
#print(coords); print(text); print(col)
	}
	rect(coord["xl",], coord["yb",], coord["xr",], coord["yt",], col = col, border=NA)
	text(coord["xx",], coord["yy",], text, cex=cex)
}
######################################################################################
hivc.clu.plot<- function(	ph, clu, edge.col.basic="black", show.node.label= T, show.tip.label=F, file=NULL,  
							highlight.edge.of.tiplabel=NULL, highlight.edge.of.tiplabel.col="red", 
							highlight.cluster=NULL, highlight.cluster.col="red",							
							pdf.scaley=10, pdf.width= 7, pdf.height=pdf.scaley*7, pdf.off=1, pdf.xlim=NULL,
							cex.nodelabel=0.5, cex.edge.incluster=1, cex.edge.outcluster= cex.edge.incluster/3, no.margin=T, ...)
{
	require(colorspace)	
	clu.edge							<- clu[ ph$edge[,1] ]
	clu.edge							<- clu.edge+1				#set col index for non-clustering edges to 1
	clu.edge[is.na(clu.edge)]			<- 1
	#cols.n								<- length(unique(clu.edge))-1	
	#cols								<- c("black",rainbow_hcl(cols.n, start = 30, end = 300))	
	clu.edge.col						<- rep(edge.col.basic, nrow(ph$edge))	#cols[clu.edge]
	clu.edge.col[clu.edge==1]			<- "grey50" 
	if(!is.null(highlight.cluster))
	{
		if(length(highlight.cluster.col)==1)
			highlight.cluster.col<- rep(highlight.cluster.col,length(highlight.cluster))
		for(i in seq_along(highlight.cluster))
			clu.edge.col[clu.edge%in%(highlight.cluster[[i]]+1)]<- 	highlight.cluster.col[i]	
	}	
	if(!is.null(highlight.edge.of.tiplabel))
	{
		highlight.edge<- lapply(highlight.edge.of.tiplabel, function(x)		which( ph$edge[,2]%in%which( substr(ph$tip.label, 1, nchar(x))==x ) )		)					
		if(length(highlight.edge.of.tiplabel.col)==1)
			highlight.edge.of.tiplabel.col	<- rep(highlight.edge, length(highlight.edge.of.tiplabel.col) )		
		for(i in seq_along(highlight.edge))
			clu.edge.col[highlight.edge[[i]]]		<- highlight.edge.of.tiplabel.col[i]		
	}	
	clu.edge.width						<- rep(cex.edge.outcluster, length(clu.edge))
	clu.edge.width[clu.edge!=1]			<- cex.edge.incluster
	clu.edge.lty						<- rep(1, length(clu.edge))
	clu.edge.lty[clu.edge!=1]			<- 1
	if(class(file)=="character")
		pdf(file,width=pdf.width,height=pdf.height)
	if(no.margin) 	par(mar=c(0,0,0,0))
	plot.coordinates					<- plot(ph, show.tip.label=show.tip.label, show.node.label=show.node.label, cex=cex.nodelabel, edge.color=clu.edge.col, edge.width=clu.edge.width, edge.lty=clu.edge.lty, x.lim=pdf.xlim, ...)
	if(class(file)=="character" && pdf.off)
		dev.off()
	plot.coordinates
}



