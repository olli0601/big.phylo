#!/bin/sh

#######################################################
# start: create and check bootstrap alignment
#######################################################
/Users/Oliver/Library/R/2.15/library/big.phylo/misc/startme.R -exe=BOOTSTRAPSEQ -resume=0 -bootstrap=0 -by=codon -indir=/Users/Oliver/git/big.phylo/pkg/example -infile=nz_h3n2 -outdir=/Users/Oliver/git/big.phylo/pkg/example
#######################################################
# end: create and check bootstrap alignment
#######################################################
#######################################################
# start: compute ExaML tree
####################################################### 
echo 'run /Users/Oliver/Library/R/2.15/library/big.phylo/ext/ExaML-parser'
[ -s /Users/Oliver/git/big.phylo/pkg/example/ExaML_result.nz_h3n2.finaltree.000 ] && [ -s /Users/Oliver/git/big.phylo/pkg/example/ExaML_info.nz_h3n2.finaltree.000 ] && exit 1
CWDEXAML=$(pwd)
cd /Users/Oliver/git/big.phylo/pkg/example
/Users/Oliver/Library/R/2.15/library/big.phylo/ext/ExaML-parser -m DNA -s /Users/Oliver/git/big.phylo/pkg/example/nz_h3n2.phylip.000 -n nz_h3n2.phylip.examl.000 
echo 'end /Users/Oliver/Library/R/2.15/library/big.phylo/ext/ExaML-parser' 
echo 'run /Users/Oliver/Library/R/2.15/library/big.phylo/ext/ExaML-parsimonator'
/Users/Oliver/Library/R/2.15/library/big.phylo/ext/ExaML-parsimonator -p53354 -s /Users/Oliver/git/big.phylo/pkg/example/nz_h3n2.phylip.000 -n nz_h3n2.starttree.000 
echo 'end /Users/Oliver/Library/R/2.15/library/big.phylo/ext/ExaML-parsimonator' 
echo 'run /Users/Oliver/Library/R/2.15/library/big.phylo/ext/examl'
mpirun -np 1 /Users/Oliver/Library/R/2.15/library/big.phylo/ext/examl -m GAMMA -D -s nz_h3n2.phylip.examl.000.binary -t RAxML_parsimonyTree.nz_h3n2.starttree.000.0 -n nz_h3n2.finaltree.000 
echo 'end /Users/Oliver/Library/R/2.15/library/big.phylo/ext/examl' 
echo 'start cleanup'
find . -name '*phylip.000*' -delete
find . -name '*phylip.examl.000*' -delete
find . -name '*starttree.000*' -delete
find . -name 'ExaML_binaryCheckpoint.*?finaltree.000*' -delete
find . -name 'ExaML_log.*?finaltree.000*' -delete 
echo 'end cleanup'
cd $CWDEXAML
#######################################################
# end: compute ExaML tree
#######################################################
#######################################################
# start: check if all ExaML boostrap trees have been computed and if yes create ExaML bootstrap tree
#######################################################
CWD=$(pwd)
cd /Users/Oliver/git/big.phylo/pkg/example 
echo 'check if all bootstrap samples have been computed'
if [ $(find . -name 'ExaML_result.nz_h3n2.finaltree*' | wc -l) -eq 500 ]; then 
	echo 'all bootstrap samples computed -- find best tree and add bootstrap support values'
	for i in $(seq 0 499); do cat ExaML_result.nz_h3n2.finaltree.$(printf %03d $i) >> ExaML_result.nz_h3n2.bstree; done
	BSBEST=$(grep 'Likelihood of best tree' ExaML_info.* | awk '{print $5,$1;}' | sort -n | tail -1 | grep -o 'finaltree.*' | cut -d':' -f 1) 
	echo "best tree is $BSBEST"
	/Users/Oliver/Library/R/2.15/library/big.phylo/ext/ExaML-raxml -f b -m GTRCAT -s nz_h3n2.phylip.examl.binary -t ExaML_result.nz_h3n2.$BSBEST -z ExaML_result.nz_h3n2.bstree -n nz_h3n2.supporttree 
	echo 'all bootstrap samples computed -- found best tree and added bootstrap support values' 
	echo 'start cleanup' 
	rm RAxML_info*
	mv RAxML_bipartitions.nz_h3n2.supporttree nz_h3n2_examlbs500.newick
	rm RAxML_bipartitionsBranchLabels.nz_h3n2.supporttree
	rm ExaML_result.nz_h3n2.bstree
	zip nz_h3n2_examlout.zip  ExaML_result.nz_h3n2.* ExaML_info.nz_h3n2.*
	rm ExaML_result.nz_h3n2.* ExaML_info.nz_h3n2.* 
	echo 'end cleanup'
fi
cd $CWD
#######################################################
# end: check if all ExaML boostrap trees have been computed and if yes create ExaML bootstrap tree
#######################################################
