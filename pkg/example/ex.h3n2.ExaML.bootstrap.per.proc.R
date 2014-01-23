require(big.phylo)
data(nz_h3n2)
#	H3N2 data is stored in 'seq' DNAbin matrix object
print( seq )							
#	write sequences to directory for processing
indir		<- getwd()
infile		<- 'nz_h3n2.R'
insignat	<- ''
save(seq, file=paste(indir, infile, sep='/'))	
#	Write the a shell file to create the first ExaML boostrap tree
pipeline.ExaML.bootstrap.per.proc(indir, infile, bs.from=0, bs.to=0, hpc.walltime=1, hpc.q=NA, hpc.mem="500mb", hpc.nproc=1)
#	Resume and write the shell files to create two more ExaML boostrap trees
pipeline.ExaML.bootstrap.per.proc(indir, infile, bs.from=1, bs.to=2, hpc.walltime=1, hpc.q=NA, hpc.mem="500mb", hpc.nproc=1)
#
#	NOTE: if the current Unix system does not have a 'qsub' command, then a shell file will be generated that can
#	be run manually on the current Unix system. Support for an HPC system is currently implemented for the Imperial
#	server 'cx1.hpc.ic.ac.uk'. To see how the final HPC script would look like, try
#
#	pipeline.ExaML.bootstrap.per.proc(indir, infile, bs.from=1, bs.to=2, hpc.sys='cx1.hpc.ic.ac.uk', hpc.walltime=1, hpc.q=NA, hpc.mem="500mb", hpc.nproc=1)
#
#	It is not difficult to enable support for a different server.
#	See ?cmd.hpcwrapper
