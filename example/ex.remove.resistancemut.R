#
#	remove HIV drug resistance mutations from an alignment
#
require(ape)
require(data.table)
require(big.phylo)
#
#	using the R function: this requires an ape DNAbin object
#
load( system.file(package=PR.PACKAGE, "data", "HIV_example.rda") )
print(seq)
seq.rm.drugresistance(seq)
#
#	using Rscript: this requires fasta file name
#	The snippet below shows how to call Rscript from the command line
#
indir			<- system.file(package=PR.PACKAGE, "data")
infile			<- "HIV_example.fasta"
outdir			<- "XXX"
outfile			<- "HIV_example_NoDRM.fasta"
cat( cmd.rm.resistance(indir, infile, outfile, outdir=outdir) ) 
