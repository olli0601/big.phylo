require(big.phylo)
data(nz_h3n2)
#	H3N2 data is stored in 'seq' DNAbin matrix object
print( seq )							
#	write sequences to directory for processing
indir		<- getwd()
infile		<- 'nz_h3n2.R'
insignat	<- ''
save(seq, file=paste(indir, infile, sep='/'))	
#	create the command string
bs.id		<- 1
infile		<- substr(infile, 1, nchar(infile)-2)
argv		<<- cmd.examl.bsalignment(indir, infile, bs.id) 
argv		<<- unlist(strsplit(argv,' '))
#	create the bootstrap alignment
prog.examl.getbootstrapseq()
