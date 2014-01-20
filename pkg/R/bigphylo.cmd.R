#' This file contains R functions that provide an interface to HPC scripting, and the 3SEQ and ExaML programs
#' Shell scripts are generated that can be either run directly, or submitted to an HPC system.
#' At present, cmd.hpccaller provides an interface to the CX1B high performance system that runs PBS (portable batch system)

PR.PACKAGE					<- "recombination.analyzer"
PR.STARTME					<- system.file(package=PR.PACKAGE, "misc", "3seq.startme.R") 
PR.EXAML.BSCREATE			<- paste(PR.STARTME,"-exeBOOTSTRAPSEQ",sep=' ')
PR.RECOMB.3SEQ				<- system.file(package=PR.PACKAGE, "ext", "3seq") 
PR.RECOMB.PROCESS3SEQOUTPUT	<- paste(PR.STARTME,"-exe=RECOMB.PROCESS3SEQOUT",sep=' ')
PR.RECOMB.CHECKCANDIDATES	<- paste(PR.STARTME,"-exe=RECOMB.CHECKCANDIDATES",sep=' ')
PR.RECOMB.PLOTINCONGRUENCE	<- paste(PR.STARTME,"-exe=RECOMB.PLOTINCONGRUENCE",sep=' ')
PR.EXAML.PARSER				<- system.file(package=PR.PACKAGE, "ext", "ExaML-parser") 
PR.EXAML.STARTTREE			<- system.file(package=PR.PACKAGE, "ext", "ExaML-parsimonator")
PR.EXAML.EXAML				<- system.file(package=PR.PACKAGE, "ext", "examl")
PR.EXAML.BS					<- system.file(package=PR.PACKAGE, "ext", "ExaML-raxml")
HPC.NPROC					<- {tmp<- c(1,4); names(tmp)<- c("debug","cx1.hpc.ic.ac.uk"); tmp}
HPC.MPIRUN					<- {tmp<- c("mpirun","mpiexec"); names(tmp)<- c("debug","cx1.hpc.ic.ac.uk"); tmp}
HPC.CX1.IMPERIAL			<- "cx1.hpc.ic.ac.uk"		#this is set to system('domainname',intern=T) for the hpc cluster of choice
HPC.MEM						<- "1750mb"
HPC.CX1.IMPERIAL.LOAD		<- "module load intel-suite mpi R/2.15"

#generate 3seq command
#' @export
cmd.recombination.run.3seq<- function(infile, outfile=paste(infile,".3s.rec",sep=''), recomb.3seq.siglevel=0.1, recomb.3seq.testvsall.beginatseq=NA, recomb.3seq.testvsall.endatseq=NA, prog= PR.RECOMB.3SEQ, nproc=1, verbose=1)
{
	cmd<- "#######################################################
# start: run 3Seq
#######################################################"
	cmd<- paste(cmd,paste("\necho \'run ",prog,"\'\n",sep=''))				
	#default commands
	cmd<- paste(cmd,prog,sep=" ")	
	cmd<- paste(cmd, " -x ",infile," -id ",outfile, sep='' )	
	if(!is.na(recomb.3seq.testvsall.beginatseq))
		cmd<- paste(cmd, " -b",recomb.3seq.testvsall.beginatseq, sep='' )
	if(!is.na(recomb.3seq.testvsall.endatseq))
		cmd<- paste(cmd, " -e",recomb.3seq.testvsall.endatseq, sep='' )	
	cmd<- paste(cmd, " -hs -t",recomb.3seq.siglevel, sep='')					
	#verbose stuff
	cmd<- paste(cmd,paste("\necho \'end ",prog,"\'\n",sep=''))
	cmd<- paste(cmd,"#######################################################
# end: run 3Seq
#######################################################\n",sep='')
	cmd			
}

#process 3seq output
#' @export
cmd.recombination.process.3SEQ.output<- function(indir, infile, insignat, prog= PR.RECOMB.PROCESS3SEQOUTPUT, resume=1, verbose=1)
{
	cmd<- "#######################################################
# start: run prog.recom.process.3SEQ.output
#######################################################"
	cmd<- paste(cmd,paste("\necho \'run ",prog,"\'\n",sep=''))
	#default commands
	cmd<- paste(cmd,prog," -v=",verbose," -resume=",resume,sep='')
	cmd<- paste(cmd," -indir=",indir," -infile=",infile," -insignat=",insignat,sep='')
	#verbose stuff
	cmd<- paste(cmd,paste("\necho \'end ",prog,"\'\n",sep=''))
	cmd<- paste(cmd,"#######################################################
# end: run prog.recom.process.3SEQ.output
#######################################################\n",sep='')
	cmd	
}

#check 3seq candidate recombinants
#' @export
cmd.recombination.check.candidates<- function(indir, infile, insignat, triplet.id, prog= PR.RECOMB.CHECKCANDIDATES, resume=1, verbose=1,hpc.walltime=NA, hpc.q=NA, hpc.mem=NA, hpc.nproc=NA)
{
	cmd<- "#######################################################
# start: prog.recom.get.incongruence
#######################################################"
	cmd<- paste(cmd,paste("\necho \'run ",prog,"\'\n",sep=''))
	#default commands
	cmd<- paste(cmd,prog," -v=",verbose," -resume=",resume,sep='')
	cmd<- paste(cmd," -indir=",indir," -infile=",infile," -insignat=",insignat," -tripletid=",triplet.id,sep='')	
	if(!is.na(hpc.walltime))	cmd<- paste(cmd," -hpc.walltime=",hpc.walltime,sep='')
	if(!is.na(hpc.q))			cmd<- paste(cmd," -hpc.q=",hpc.q,sep='')
	if(!is.na(hpc.mem))			cmd<- paste(cmd," -hpc.mem=",hpc.mem,sep='') 
	if(!is.na(hpc.nproc))		cmd<- paste(cmd," -hpc.nproc=",hpc.nproc,sep='')	
	#verbose stuff
	cmd<- paste(cmd,paste("\necho \'end ",prog,"\'\n",sep=''))
	cmd<- paste(cmd,"#######################################################
# end: prog.recom.get.incongruence
#######################################################\n",sep='')
	cmd	
}

cmd.recombination.plot.incongruence<- function(indir, infile, insignat, triplet.id=NA, prog= PR.RECOMB.PLOTINCONGRUENCE, opt.select=NA,verbose=1)
{
	cmd<- "#######################################################
# start: prog.recom.plot.incongruence
#######################################################"
	cmd<- paste(cmd,paste("\necho \'run ",prog,"\'\n",sep=''))
	#default commands
	cmd<- paste(cmd,prog," -v=",verbose,sep='')
	cmd<- paste(cmd," -indir=",indir," -infile=",infile," -insignat=",insignat, sep='')
	if(!is.na(opt.select))
		cmd<- paste(cmd," -select=",opt.select,sep='')
	if(!is.na(triplet.id))
		cmd<- paste(cmd," -tripletid=",triplet.id,sep='')
	#verbose stuff
	cmd<- paste(cmd,paste("\necho \'end ",prog,"\'\n",sep=''))
	cmd<- paste(cmd,"#######################################################
# end: prog.recom.plot.incongruence
#######################################################\n",sep='')
	cmd	
}

#' @export
cmd.examl<- function(indir, infile, signat.in, signat.out, outdir=indir, prog.parser= PR.EXAML.PARSER, prog.starttree= PR.EXAML.STARTTREE, args.starttree.seed=12345, args.starttree.bsid= NA, prog.examl= PR.EXAML.EXAML, args.examl="-m GAMMA -D", resume=1, verbose=1)
{
	if(is.na(args.starttree.bsid))
		args.starttree.bsid	<- "000"
	else
		args.starttree.bsid	<-	sprintf("%03d",args.starttree.bsid)
	cmd<- "#######################################################
# start: compute ExaML tree
#######################################################"
	cmd<- paste(cmd,paste("\necho \'run ",prog.parser,"\'\n",sep=''))
	#if output files are found and resume, don t do anything
	if(resume)
	{
		cmd		<- paste(cmd,"[ -s ",outdir,'/ExaML_result.',infile,'_',signat.in,".finaltree.",args.starttree.bsid," ] && ", sep='')	#if file non-zero
		cmd		<- paste(cmd,"[ -s ",outdir,'/ExaML_info.',infile,'_',signat.in,".finaltree.",args.starttree.bsid," ] ", sep='')		#and if file non-zero
		cmd		<- paste(cmd,"&& exit 1\n",sep='')		#then exit
	}
	#default commands for parser					
	cmd			<- paste(cmd,"CWDEXAML=$(pwd)\n",sep='')
	cmd			<- paste(cmd,"cd ",outdir,'\n',sep='')
	tmp			<- paste(indir,paste(infile,'_',signat.in,".phylip.",args.starttree.bsid,sep=''),sep='/')
	cmd			<- paste(cmd,prog.parser," -m DNA -s ",tmp,sep='')
	tmp			<- paste(infile,'_',signat.out,".phylip.examl.",args.starttree.bsid,sep='')
	cmd			<- paste(cmd," -n ",tmp,sep='')
	#verbose stuff for parser	
	cmd			<- paste(cmd,paste("\necho \'end ",prog.parser,"\'",sep=''))
		
	cmd			<- paste(cmd,paste("\necho \'run ",prog.starttree,"\'\n",sep=''))
	#default commands for start tree
	tmp			<- paste(indir,paste(infile,'_',signat.in,".phylip.",args.starttree.bsid,sep=''),sep='/')	
	cmd			<- paste(cmd,prog.starttree," -p",args.starttree.seed," -s ",tmp,sep='')	
	tmp			<- paste(infile,'_',signat.out,".starttree.",args.starttree.bsid,sep='')
	cmd			<- paste(cmd," -n ",tmp,sep='')
	#verbose stuff
	cmd			<- paste(cmd,paste("\necho \'end ",prog.starttree,"\'",sep=''))
	
	cmd			<- paste(cmd,paste("\necho \'run ",prog.examl,"\'\n",sep=''))
	#default commands for final tree
	tmp			<- cmd.hpcsys()
	if(tmp=="debug")
		cmd		<- paste(cmd,HPC.MPIRUN[tmp]," -np ",HPC.NPROC[tmp],' ',prog.examl,' ',args.examl,sep='')
	else if(tmp==HPC.CX1.IMPERIAL)
		cmd		<- paste(cmd,HPC.MPIRUN[tmp],' ',prog.examl,' ',args.examl,sep='')
	else	
		stop("unknown hpc sys")
	tmp			<- paste(infile,'_',signat.out,".phylip.examl.",args.starttree.bsid,".binary",sep='')
	cmd			<- paste(cmd," -s ",tmp,sep='')
	tmp			<- paste("RAxML_parsimonyTree.",infile,'_',signat.out,".starttree.",args.starttree.bsid,".0",sep='')
	cmd			<- paste(cmd," -t ",tmp,sep='')
	tmp			<- paste(infile,'_',signat.out,".finaltree.",args.starttree.bsid,sep='')
	cmd			<- paste(cmd," -n ",tmp,sep='')		
	cmd			<- paste(cmd,paste("\necho \'end ",prog.examl,"\'",sep=''))
	
	#delete ExaML output that is not further needed 
	cmd			<- paste(cmd,paste("\necho \'start cleanup\'",sep=''))	
	cmd			<- paste(cmd,"\nfind . -name \'*phylip.",args.starttree.bsid,"*\' -delete",sep='')
	cmd			<- paste(cmd,"\nfind . -name \'*phylip.examl.",args.starttree.bsid,"*\' -delete",sep='')
	cmd			<- paste(cmd,"\nfind . -name \'*starttree.",args.starttree.bsid,"*\' -delete",sep='')
	cmd 		<- paste(cmd,"\nfind . -name \'ExaML_binaryCheckpoint.*?finaltree.",args.starttree.bsid,"*\' -delete", sep='' )
	cmd 		<- paste(cmd,"\nfind . -name \'ExaML_log.*?finaltree.",args.starttree.bsid,"*\' -delete", sep='' )
	cmd			<- paste(cmd,paste("\necho \'end cleanup\'",sep=''))		
	cmd			<- paste(cmd,"\ncd $CWDEXAML",sep='')
	cmd			<- paste(cmd,"\n#######################################################
# end: compute ExaML tree
#######################################################\n",sep='')
	cmd
}

# 	creates a shell command to create a new bootstrap alignment over codon positions of an input alignment
#' @export 
cmd.examl.bsalignment<- function(indir, infile, signat.in, signat.out, bs.id, outdir=indir, prog.bscreate= PR.EXAML.BSCREATE, opt.bootstrap.by="codon",resume=0, verbose=1)
{
	cmd			<- paste("#######################################################
# start: create and check bootstrap alignment
#######################################################\n",sep='')
	cmd			<- paste(cmd,prog.bscreate," -resume=",resume," -bootstrap=",bs.id," -by=",opt.bootstrap.by,sep='')
	cmd			<- paste(cmd," -indir=",indir," -infile=",infile," -outdir=",outdir,sep='')
	cmd			<- paste(cmd," -insignat=",signat.in," -outsignat=",signat.out,sep='')
	cmd			<- paste(cmd,"\n#######################################################
# end: create and check bootstrap alignment
#######################################################",sep='')			
	cmd
}

#' @export
cmd.examl.bootstrap.on.one.machine<- function(indir, infile, signat.in, signat.out, bs.from=0, bs.to=99, bs.n=bs.to-bs.from+ifelse(bs.from==0,1,0), outdir=indir, prog.parser= PR.EXAML.PARSER, prog.starttree= PR.EXAML.STARTTREE, prog.examl=PR.EXAML.EXAML, opt.bootstrap.by="codon", args.examl="-m GAMMA -D", prog.supportadder=PR.EXAML.BS, tmpdir.prefix="examl", resume=1, verbose=1)
{
	hpcsys			<- cmd.hpcsys()
	hpcsys			<- "cx1.hpc.ic.ac.uk"
	#create number of seeds for the number of runs being processed, which could be less than bs.n
	bs.id			<- seq.int(bs.from,bs.to)
	bs.seeds		<- floor( runif(length(bs.id), 1e4, 1e5-1) )
	tmpdir			<- paste("$CWD/",tmpdir.prefix,'_',format(Sys.time(),"%y-%m-%d-%H-%M-%S"),sep='')
	cmd				<- sapply(seq_along(bs.seeds), function(i)
			{
				cmd			<- ''
				if(hpcsys=="debug")						#my MAC - don t use scratch
				{
					cmd		<- paste(cmd,cmd.examl.bsalignment(indir, infile, signat.in, signat.out, bs.id[i], opt.bootstrap.by=opt.bootstrap.by, outdir=indir, verbose=verbose),sep='\n')
					cmd		<- paste(cmd,cmd.examl(indir, infile, signat.in, signat.out, outdir=outdir, prog.parser= prog.parser, prog.starttree= prog.starttree, args.starttree.seed=bs.seeds[i], args.starttree.bsid= bs.id[i], prog.examl=prog.examl, args.examl=args.examl, resume=resume, verbose=verbose),sep='\n')
				}
				else if(hpcsys=="cx1.hpc.ic.ac.uk")		#imperial - use scratch directory
				{										
					if(i==1)
					{
						cmd	<- paste(cmd,"\nCWD=$(pwd)",sep='')
						cmd	<- paste(cmd,"\necho $CWD",sep='')
						cmd	<- paste(cmd,"\nmkdir -p ",tmpdir,sep='')
						tmp	<- paste(indir,'/',infile,'_',gsub('/',':',signat.in),".R",sep='')
						cmd	<- paste(cmd,"\ncp ",tmp," ",tmpdir,sep='')
					}
					cmd		<- paste(cmd,cmd.examl.bsalignment(tmpdir, infile, signat.in, signat.out, bs.id[i], opt.bootstrap.by=opt.bootstrap.by, outdir=tmpdir, verbose=verbose),sep='\n')
					cmd		<- paste(cmd,cmd.examl(tmpdir, infile, signat.in, signat.out, outdir=tmpdir, prog.parser= prog.parser, prog.starttree= prog.starttree, args.starttree.seed=bs.seeds[i], args.starttree.bsid= bs.id[i], prog.examl=prog.examl, args.examl=args.examl, resume=resume, verbose=verbose),sep='\n')
				}
				cmd
			})
	cmd				<- paste(cmd, collapse='')
	cmd				<- paste(cmd,"#######################################################
# start: check if all ExaML boostrap trees have been computed and if yes create ExaML bootstrap tree
#######################################################",sep='')	
	cmd			<- paste(cmd,"\nCWD=$(pwd)\n",sep='')
	cmd			<- paste(cmd,"cd ",tmpdir,sep='') 
	cmd			<- paste(cmd,paste("\necho \'check if all bootstrap samples have been computed\'",sep=''))
	#compute bs tree even when some errors	
	tmp			<- paste("\nif [ $(find . -name 'ExaML_result.",infile,'_',gsub('/',':',signat.in),".finaltree*' | wc -l) -gt ",round(bs.n*0.95)," ]; then",sep='')
	cmd			<- paste(cmd,tmp,sep='')
	cmd			<- paste(cmd,paste("\n\techo \'all bootstrap samples computed -- find best tree and add bootstrap support values\'",sep=''))				
	tmp			<- c(	paste("ExaML_result.",infile,'_',signat.out,".finaltree",sep=''), 	paste("ExaML_result.",infile,'_',signat.out,".bstree",sep=''))
	#add all bootstrap final trees into bs tree file
	cmd			<- paste(cmd,"\n\tfor i in $(seq 0 ",bs.n-1,"); do cat ",tmp[1],".$(printf %03d $i) >> ",tmp[2],"; done",sep='')
	#identify suffix finaltree.XXX of final tree with highest likelihood
	cmd			<- paste(cmd,"\n\tBSBEST=$(grep 'Likelihood of best tree' ExaML_info.* | awk '{print $5,$1;}' | sort -n | tail -1 | grep -o 'finaltree.*' | cut -d':' -f 1)",sep='')
	cmd			<- paste(cmd,paste("\n\techo \"best tree is $BSBEST\"",sep=''))		
	#create final tree with bootstrap support values
	tmp			<- c(	paste(infile,'_',signat.out,".phylip.examl.binary",sep=''),	paste("ExaML_result.",infile,'_',signat.out,".$BSBEST",sep=''),	paste("ExaML_result.",infile,'_',signat.out,".bstree",sep=''), paste(infile,'_',signat.out,".supporttree",sep=''))
	cmd			<- paste(cmd,"\n\t",prog.supportadder," -f b -m GTRCAT -s ",tmp[1]," -t ",tmp[2]," -z ",tmp[3]," -n ",tmp[4],sep='' )
	cmd			<- paste(cmd,paste("\n\techo \'all bootstrap samples computed -- found best tree and added bootstrap support values\'",sep=''))										
	#delete ExaML output that is not further needed
	cmd			<- paste(cmd,paste("\n\techo \'start cleanup\'",sep=''))
	cmd			<- paste(cmd,"\n\trm RAxML_info*",sep=' ')								
	cmd			<- paste(cmd,paste("\n\tmv RAxML_bipartitions.",infile,'_',signat.out,".supporttree ",infile,"_examlbs",bs.n,'_',signat.out,".newick",sep=''),sep='')				
	cmd			<- paste(cmd,paste("\n\trm RAxML_bipartitionsBranchLabels.",infile,'_',signat.out,".supporttree",sep=''),sep='')
	cmd			<- paste(cmd,paste("\n\trm ExaML_result.",infile,'_',signat.out,".bstree",sep=''),sep='')									
	cmd			<- paste(cmd,paste("\n\techo \'end cleanup\'",sep=''))
	#copy bstree to outdir
	cmd			<- paste(cmd,paste("\n\tcp -f ",infile,"_examlbs",bs.n,'_',signat.out,".newick",' ',outdir,sep=''),sep='')
	cmd			<- paste(cmd,"\nfi",sep='')
	cmd			<- paste(cmd,"\n#######################################################
# end: check if all ExaML boostrap trees have been computed and if yes create ExaML bootstrap tree
#######################################################",sep='')
cmd			<- paste(cmd,"\n#######################################################
# start: zip and copy ExaML output
#######################################################",sep='')	
	#outside if:	zip and copy ExaML output to outdir just in case something went wrong
	cmd			<- paste(cmd,paste("\nzip ",infile,'_examlout_',signat.out,".zip  ExaML_result.",infile,'_',signat.out,".* ExaML_info.",infile,'_',signat.out,".*",sep=''),sep='')
	cmd			<- paste(cmd,paste("\ncp -f ",infile,'_examlout_',signat.out,".zip",' ',outdir,sep=''),sep='')
	cmd			<- paste(cmd,paste("\nrm ExaML_result.",infile,'_',signat.out,".* ExaML_info.",infile,'_',signat.out,".*",sep=''),sep='')
	cmd			<- paste(cmd,"\ncd $CWD",sep='')
	cmd			<- paste(cmd,"\n#######################################################
# end: zip and copy ExaML output
#######################################################\n",sep='')				
	cmd
}

#' @export
cmd.examl.bootstrap<- function(indir, infile, signat.in, signat.out, bs.from=0, bs.to=99, bs.n=bs.to-bs.from+ifelse(bs.from==0,1,0), outdir=indir, prog.parser= PR.EXAML.PARSER, prog.starttree= PR.EXAML.STARTTREE, prog.examl=PR.EXAML.EXAML, opt.bootstrap.by="codon", args.examl="-m GAMMA -D", prog.supportadder=PR.EXAML.BS, tmpdir.prefix="examl", resume=1, verbose=1)
{
	hpcsys			<- cmd.hpcsys()
	#hpcsys			<- "cx1.hpc.ic.ac.uk"
	#create number of seeds for the number of runs being processed, which could be less than bs.n
	bs.id			<- seq.int(bs.from,bs.to)
	bs.seeds		<- floor( runif(length(bs.id), 1e4, 1e5-1) )
	tmpdir.prefix	<- paste(tmpdir.prefix,'_',format(Sys.time(),"%y-%m-%d-%H-%M-%S"),sep='')
	lapply(seq_along(bs.seeds), function(i)
			{
				cmd			<- ''				
				if(hpcsys=="debug")						#my MAC - don t use scratch
				{
					cmd		<- paste(cmd,cmd.examl.bsalignment(indir, infile, signat.in, signat.out, bs.id[i], opt.bootstrap.by=opt.bootstrap.by, outdir=indir, verbose=verbose),sep='\n')
					cmd		<- paste(cmd,cmd.examl(indir, infile, signat.in, signat.out, outdir=outdir, prog.parser= prog.parser, prog.starttree= prog.starttree, args.starttree.seed=bs.seeds[i], args.starttree.bsid= bs.id[i], prog.examl=prog.examl, args.examl=args.examl, resume=resume, verbose=verbose),sep='\n')
				}
				else if(hpcsys=="cx1.hpc.ic.ac.uk")		#imperial - use scratch directory
				{										
					if(resume)
					{
						cmd	<- paste(cmd,"\nif ", sep='')
						cmd	<- paste(cmd,"[ ! -s ",outdir,'/ExaML_result.',infile,'_',signat.in,".finaltree.",sprintf("%03d",bs.id[i])," ]", sep='')
						cmd	<- paste(cmd," || ", sep='')
						cmd	<- paste(cmd,"[ ! -s ",outdir,'/ExaML_info.',infile,'_',signat.in,".finaltree.",sprintf("%03d",bs.id[i])," ];", sep='')
						cmd	<- paste(cmd," then\n",sep='')
						cmd	<- paste(cmd,"#######################################################
# start: not indented if statement -- don t do anything if ExaML output exists already
#######################################################",sep='')
					}					
					cmd		<- paste(cmd,"CWD=$(pwd)\n",sep='\n')
					cmd		<- paste(cmd,"echo $CWD\n",sep='')
					tmpdir	<- paste("$CWD/",tmpdir.prefix,sep='')
					cmd		<- paste(cmd,"mkdir -p ",tmpdir,'\n',sep='')
					tmp		<- paste(indir,'/',infile,'_',gsub('/',':',signat.in),".R",sep='')
					cmd		<- paste(cmd,"cp ",tmp," ",tmpdir,sep='')
					cmd		<- paste(cmd,cmd.examl.bsalignment(tmpdir, infile, signat.in, signat.out, bs.id[i], opt.bootstrap.by=opt.bootstrap.by, outdir=tmpdir, verbose=verbose),sep='\n')
					cmd		<- paste(cmd,cmd.examl(tmpdir, infile, signat.in, signat.out, outdir=tmpdir, prog.parser= prog.parser, prog.starttree= prog.starttree, args.starttree.seed=bs.seeds[i], args.starttree.bsid= bs.id[i], prog.examl=prog.examl, args.examl=args.examl, resume=resume, verbose=verbose),sep='\n')
					cmd		<- paste(cmd,"cp -f ",tmpdir,"/ExaML_result.",infile,'_',gsub('/',':',signat.in),".finaltree.",sprintf("%03d",bs.id[i]),' ', outdir,'\n',sep='')
					cmd		<- paste(cmd,"cp -f ",tmpdir,"/ExaML_info.",infile,'_',gsub('/',':',signat.in),".finaltree.",sprintf("%03d",bs.id[i]),' ', outdir,'\n',sep='')
					if(resume)
					{
					cmd		<- paste(cmd,"#######################################################
# end: not indented if statement -- don t do anything if ExaML output exists already
#######################################################\nelse\n\techo 'resumed bootstrap number ",bs.id[i],"'\nfi\n",sep='')					
					}
				}				
				cmd			<- paste(cmd,"#######################################################
# start: check if all ExaML boostrap trees have been computed and if yes create ExaML bootstrap tree
#######################################################",sep='')	
				cmd			<- paste(cmd,"\nCWD=$(pwd)\n",sep='')
				cmd			<- paste(cmd,"cd ",outdir,sep='') 
				cmd			<- paste(cmd,paste("\necho \'check if all bootstrap samples have been computed\'",sep=''))
				tmp			<- paste("\nif [ $(find . -name 'ExaML_result.",infile,'_',gsub('/',':',signat.in),".finaltree*' | wc -l) -eq ",bs.n," ]; then",sep='')
				cmd			<- paste(cmd,tmp,sep='')
				cmd			<- paste(cmd,paste("\n\techo \'all bootstrap samples computed -- find best tree and add bootstrap support values\'",sep=''))				
				tmp			<- c(	paste("ExaML_result.",infile,'_',signat.out,".finaltree",sep=''), 	paste("ExaML_result.",infile,'_',signat.out,".bstree",sep=''))
				#add all bootstrap final trees into bs tree file
				cmd			<- paste(cmd,"\n\tfor i in $(seq 0 ",bs.n-1,"); do cat ",tmp[1],".$(printf %03d $i) >> ",tmp[2],"; done",sep='')
				#identify suffix finaltree.XXX of final tree with highest likelihood
				cmd			<- paste(cmd,"\n\tBSBEST=$(grep 'Likelihood of best tree' ExaML_info.* | awk '{print $5,$1;}' | sort -n | tail -1 | grep -o 'finaltree.*' | cut -d':' -f 1)",sep='')
				cmd			<- paste(cmd,paste("\n\techo \"best tree is $BSBEST\"",sep=''))		
				#create final tree with bootstrap support values
				tmp			<- c(	paste(infile,'_',signat.out,".phylip.examl.binary",sep=''),	paste("ExaML_result.",infile,'_',signat.out,".$BSBEST",sep=''),	paste("ExaML_result.",infile,'_',signat.out,".bstree",sep=''), paste(infile,'_',signat.out,".supporttree",sep=''))
				cmd			<- paste(cmd,"\n\t",prog.supportadder," -f b -m GTRCAT -s ",tmp[1]," -t ",tmp[2]," -z ",tmp[3]," -n ",tmp[4],sep='' )
				cmd			<- paste(cmd,paste("\n\techo \'all bootstrap samples computed -- found best tree and added bootstrap support values\'",sep=''))										
				#delete ExaML output that is not further needed
				cmd			<- paste(cmd,paste("\n\techo \'start cleanup\'",sep=''))
				cmd			<- paste(cmd,"\n\trm RAxML_info*",sep=' ')								
				cmd			<- paste(cmd,paste("\n\tmv RAxML_bipartitions.",infile,'_',signat.out,".supporttree ",infile,"_examlbs",bs.n,'_',signat.out,".newick",sep=''),sep='')				
				cmd			<- paste(cmd,paste("\n\trm RAxML_bipartitionsBranchLabels.",infile,'_',signat.out,".supporttree",sep=''),sep='')
				cmd			<- paste(cmd,paste("\n\trm ExaML_result.",infile,'_',signat.out,".bstree",sep=''),sep='')								
				cmd			<- paste(cmd,paste("\n\tzip ",infile,'_examlout_',signat.out,".zip  ExaML_result.",infile,'_',signat.out,".* ExaML_info.",infile,'_',signat.out,".*",sep=''),sep='')
				cmd			<- paste(cmd,paste("\n\trm ExaML_result.",infile,'_',signat.out,".* ExaML_info.",infile,'_',signat.out,".*",sep=''),sep='')
				cmd			<- paste(cmd,paste("\n\techo \'end cleanup\'",sep=''))
				cmd			<- paste(cmd,"\nfi",sep='')
				cmd			<- paste(cmd,"\ncd $CWD",sep='')
				cmd			<- paste(cmd,"\n#######################################################
# end: check if all ExaML boostrap trees have been computed and if yes create ExaML bootstrap tree
#######################################################\n",sep='')				
			})
}

#' @export
cmd.examl.bsstarttree<- function(indir, infile, signat.in, signat.out, bs.from=0, bs.to=99, bs.n=bs.to-bs.from+ifelse(bs.from==0,1,0),outdir=indir, prog.parser= PR.EXAML.PARSER, prog.starttree= PR.EXAML.STARTTREE, prog.examl=PR.EXAML.EXAML, args.examl="-m GAMMA -D", prog.supportadder=PR.EXAML.BS, resume=1, verbose=1)
{
	#create number of seeds for the number of runs being processed, which could be less than bs.n
	bs.id	<- seq.int(bs.from,bs.to)
	bs.seeds<- floor( runif(length(bs.id), 1e4, 1e5-1) )
	lapply(seq_along(bs.seeds), function(i)
			{
				cmd			<- paste("cp ",paste(indir,paste(infile,'_',signat.in,".phylip",sep=''),sep='/'),sep='')
				cmd			<- paste(cmd,' ',paste(indir,paste(infile,'_',signat.in,".phylip.",sprintf("%03d",bs.id[i]),sep=''),sep='/'),"\n",sep='')
				cmd			<- paste(cmd, cmd.examl(indir, infile, signat.in, signat.out, outdir=outdir, prog.parser= prog.parser, prog.starttree= prog.starttree, args.starttree.seed=bs.seeds[i], args.starttree.bsid= bs.id[i], prog.examl=prog.examl, args.examl=args.examl, resume=resume, verbose=verbose),sep='')
				curr.dir	<- getwd()
				cmd			<- paste(cmd,"#######################################################
# start: check if all ExaML boostrap trees have been computed and if yes create ExaML bootstrap tree
#######################################################",sep='')				
				cmd			<- paste(cmd,"\ncd ",outdir,sep='')
				cmd			<- paste(cmd,paste("\necho \'check if all bootstrap samples have been computed\'",sep=''))
				tmp			<- paste("\nif [ $(find . -name 'ExaML_result*' | wc -l) -eq ",bs.n," ]; then",sep='')
				cmd			<- paste(cmd,tmp,sep='')
				cmd			<- paste(cmd,paste("\n\techo \'all bootstrap samples computed -- find best tree and add bootstrap support values\'",sep=''))				
				tmp			<- c(	paste("ExaML_result.",infile,'_',signat.out,".finaltree",sep=''), 	paste("ExaML_result.",infile,'_',signat.out,".bstree",sep=''))
				#add all bootstrap final trees into bs tree file
				cmd			<- paste(cmd,"\n\tfor i in $(seq 0 ",bs.n-1,"); do cat ",tmp[1],".$(printf %03d $i) >> ",tmp[2],"; done",sep='')
				#identify suffix finaltree.XXX of final tree with highest likelihood
				cmd			<- paste(cmd,"\n\tBSBEST=$(grep 'Likelihood of best tree' ExaML_info.* | awk '{print $5,$1;}' | sort -n | tail -1 | grep -o 'finaltree.*' | cut -d':' -f 1)",sep='')
				cmd			<- paste(cmd,paste("\n\techo \"best tree is $BSBEST\"",sep=''))		
				#create final tree with bootstrap support values
				tmp			<- c(	paste(infile,'_',signat.out,".phylip.examl.binary",sep=''),	paste("ExaML_result.",infile,'_',signat.out,".$BSBEST",sep=''),	paste("ExaML_result.",infile,'_',signat.out,".bstree",sep=''), paste(infile,'_',signat.out,".supporttree",sep=''))
				cmd			<- paste(cmd,"\n\t",prog.supportadder," -f b -m GTRCAT -s ",tmp[1]," -t ",tmp[2]," -z ",tmp[3]," -n ",tmp[4],sep='' )
				cmd			<- paste(cmd,paste("\n\techo \'all bootstrap samples computed -- found best tree and added bootstrap support values\'",sep=''))										
				#delete ExaML output that is not further needed
				cmd			<- paste(cmd,paste("\n\techo \'start cleanup\'",sep=''))
				cmd			<- paste(cmd,"\n\trm RAxML_info*",sep=' ')
				cmd			<- paste(cmd,paste("\n\trm ",infile,'_',signat.out,".phylip.examl.binary",sep=''),sep='')				
				cmd			<- paste(cmd,paste("\n\tmv RAxML_bipartitions.",infile,'_',signat.out,".supporttree ",infile,"_examlbs",bs.n,'_',signat.out,".newick",sep=''),sep='')				
				cmd			<- paste(cmd,paste("\n\trm RAxML_bipartitionsBranchLabels.",infile,'_',signat.out,".supporttree",sep=''),sep='')
				cmd			<- paste(cmd,paste("\n\trm ExaML_result.",infile,'_',signat.out,".bstree",sep=''),sep='')								
				cmd			<- paste(cmd,paste("\n\ttar -zcf ",infile,'_examlout_',signat.out,".tar.gz  ExaML_result.",infile,'_',signat.out,".* ExaML_info.",infile,'_',signat.out,".*",sep=''),sep='')
				cmd			<- paste(cmd,paste("\n\trm ExaML_result.",infile,'_',signat.out,".* ExaML_info.",infile,'_',signat.out,".*",sep=''),sep='')
				cmd			<- paste(cmd,paste("\n\techo \'end cleanup\'",sep=''))
				cmd			<- paste(cmd,"\nfi",sep='')
				cmd			<- paste(cmd,"\ncd ",curr.dir,sep='')
				cmd			<- paste(cmd,"\n#######################################################
# end: check if all ExaML boostrap trees have been computed and if yes create ExaML bootstrap tree
#######################################################\n",sep='')
				#cat(cmd)
				#stop()
				#if [ $(find -E . -name 'ExaML_result*' | wc -l)==2 ]; then echo 'hello'; fi
			})
}

#' @export
cmd.examl.cleanup<- function(outdir, prog= PR.EXAML.EXAML)
{
	cmd<- "#######################################################
# clean up after ExaML tree
#######################################################"
	cmd<- paste(cmd,paste("\necho \'clean after ",prog,"\'\n",sep=''))	
	
	tmp<- list.files(outdir, full.names=1)
	#tmp<- tmp[c(grep("examlstarttree",tmp), grep("examlbin",tmp),grep("examl.binary",tmp),grep("phylip",tmp))]
	#cmd<- paste(cmd, "\nrm ", paste(tmp,collapse=' ',sep=''), '\n', sep='')	
	cmd<- paste(cmd,"rm ",paste(paste(outdir,c("RAxML_info*","RAxML_parsimonyTree*","ExaML_binaryCheckpoint*"),sep='/'),collapse=' ',sep=''),sep=' ')
	
	cmd<- paste(cmd,paste("\necho \'cleaned up after ",prog,"\'\n",sep=''))
	cmd
}

######################################################################################
cmd.hpcsys<- function()
{
	tmp<- system('domainname',intern=T)
	if(!nchar(tmp))	tmp<- "debug"
	tmp
}

#' @export
cmd.hpcwrapper.cx1.ic.ac.uk<- function(hpc.walltime=24, hpc.mem=HPC.MEM, hpc.nproc=1, hpc.q=NA)
{
	wrap<- "#!/bin/sh"
	tmp	<- paste("#PBS -l walltime=",hpc.walltime,":59:59,pcput=",hpc.walltime,":45:00",sep='')
	wrap<- paste(wrap, tmp, sep='\n')		
	tmp	<- paste("#PBS -l select=1:ncpus=",hpc.nproc,":mem=",hpc.mem,sep='')
	wrap<- paste(wrap, tmp, sep='\n')
	wrap<- paste(wrap, "#PBS -j oe", sep='\n')
	if(!is.na(hpc.q))
		wrap<- paste(wrap, paste("#PBS -q",hpc.q), sep='\n\n')
	wrap<- paste(wrap, HPC.CX1.IMPERIAL.LOAD, sep='\n')
	wrap
}

#add additional high performance computing information 
#' @export
cmd.hpcwrapper<- function(cmd, hpc.sys= cmd.hpcsys(), hpc.walltime=24, hpc.mem=HPC.MEM, hpc.nproc=1, hpc.q=NA)
{	
	#hpc.sys<- HPC.CX1.IMPERIAL
	if(hpc.sys==HPC.CX1.IMPERIAL)
		wrap<- cmd.hpcwrapper.cx1.ic.ac.uk(hpc.walltime=hpc.walltime, hpc.mem=hpc.mem, hpc.nproc=hpc.nproc, hpc.q=hpc.q)
	else
	{
		wrap<- "#!/bin/sh"
		cat(paste("\ndetected no HPC system and no hpcwrapper generated, domain name is",hpc.sys))
	}
	cmd<- lapply(seq_along(cmd),function(i){	paste(wrap,cmd[[i]],sep='\n')	})
	if(length(cmd)==1)
		cmd<- unlist(cmd)
	cmd	
}

#create high performance computing qsub file and submit
#' @export
cmd.hpccaller<- function(outdir, outfile, cmd)
{
	if( nchar( Sys.which("qsub") ) )
	{
		file	<- paste(outdir,'/',outfile,'.qsub',sep='')
		cat(paste("\nwrite HPC script to",file,"\n"))
		cat(cmd,file=file)
		cmd		<- paste("qsub",file)
		cat( cmd )
		cat( system(cmd, intern=TRUE) )
		Sys.sleep(1)
	}
	else
	{
		file	<- paste(outdir,'/',outfile,'.sh',sep='')
		cat(paste("\nwrite Shell script to\n",file,"\nNo 'qsub' function detected to submit the shell script automatically.\nStart this shell file manually\n"))
		cat(cmd,file=file)
		Sys.chmod(file, mode = "777")		
	}
	
}