require(recombination.analyzer)
data(nz_h3n2)
#H3N2 data is stored in 'seq' DNAbin matrix object
print( seq )							
#write sequences to directory for processing
indir		<- getwd()
infile		<- 'nz_h3n2.R'
insignat	<- ''
save(seq, file=paste(indir, infile, sep='/'))	
#Run 3SEQ in a single batch job
pipeline.recom.run.3seq(indir, infile, batch.n=1, hpc.walltime=1, hpc.q=NA, hpc.mem="500mb", hpc.nproc=1)
