require(big.phylo)
data(nz_h3n2)
#	H3N2 data is stored in 'seq' DNAbin matrix object
print( seq )							
#	write sequences to directory for processing
indir		<- getwd()
infile		<- 'nz_h3n2.R'
insignat	<- ''
save(seq, file=paste(indir, infile, sep='/'))	
#	Create the first ExaML boostrap tree
pipeline.ExaML.bootstrap.per.proc(indir, infile, bs.from=0, bs.to=0, hpc.walltime=1, hpc.q=NA, hpc.mem="500mb", hpc.nproc=1)