require(big.phylo)
#	specify file that holds an R DNAbin object
infile			<- 'myfile.R'
indir			<- 'mydir'
save(seq, file=paste(indir, infile, sep='/'))	
#	create the command string
infile			<- substr(infile, 1, nchar(infile)-2)
outfile			<- paste(infile,'out',sep='_')
argv			<<- cmd.rm.resistance(indir, infile, outfile) 
argv			<<- unlist(strsplit(argv,' '))
#	create the bootstrap alignment
#prog.remove.resistancemut()
