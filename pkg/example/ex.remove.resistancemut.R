require(big.phylo)
#	specify file that holds an R DNAbin object
infile			<- 'myfile.R'
indir			<- 'mydir'
save(seq, file=paste(indir, infile, sep='/'))	
#	create the command string
infile			<- substr(infile, 1, nchar(infile)-2)
outfile			<- paste(infile,'out',sep='_')
alignment.start	<- 2253
argv			<<- cmd.rm.resistance(indir, infile, outfile, alignment.start=alignment.start) 
argv			<<- unlist(strsplit(argv,' '))
#	create the bootstrap alignment
#prog.remove.resistancemut()