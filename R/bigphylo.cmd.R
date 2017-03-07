PR.PACKAGE					<- "big.phylo" 
PR.EXAML.BSCREATE			<- paste('Rscript', system.file(package=PR.PACKAGE, "create.bootstrapalignment.Rscript"))
PR.RM.RESISTANCE			<- paste('Rscript', system.file(package=PR.PACKAGE, "rm.drm.Rscript"))
PR.STRIP.GAPS				<- paste('Rscript', system.file(package=PR.PACKAGE, "strip.gaps.Rscript"))
PR.TREEDATER				<- paste('Rscript', system.file(package=PR.PACKAGE, "treedater.Rscript"))
PR.EXAML.PARSER				<- system.file(package=PR.PACKAGE, "ext", "ExaML-parser") 
PR.EXAML.STARTTREE			<- system.file(package=PR.PACKAGE, "ext", "ExaML-parsimonator")
PR.EXAML.EXAML				<- system.file(package=PR.PACKAGE, "ext", "examl")
PR.FASTTREE					<- system.file(package=PR.PACKAGE, "ext", "FastTree")
PR.LSD						<- system.file(package=PR.PACKAGE, "ext", "lsd")
PR.LSDDATES					<- paste('Rscript', system.file(package=PR.PACKAGE, "lsd.dates.Rscript"))
PR.MVR						<- paste('Rscript',system.file(package=PR.PACKAGE, "big.mvr.Rscript"),sep=' ')
PR.PHYD						<- paste('java -jar ', system.file(package=PR.PACKAGE, "ext", "PhyDstar.jar"), sep='')
PR.JMODELTEST				<- paste('java -jar ', system.file(package=PR.PACKAGE, "ext", "jmodeltest-2.1.10", "jModelTest.jar"), sep='')
PR.EXAML.BS					<- system.file(package=PR.PACKAGE, "ext", "ExaML-raxml")
HPC.MEM						<- "1750mb"
HPC.CX1.IMPERIAL.LOAD		<- "module load intel-suite mpi R/3.2.0"

#'	@export 
cmd.examl<- function(indir, infile, outdir=indir, prog.mpi='mpiexec', prog.parser= PR.EXAML.PARSER, args.parser="-m DNA",prog.starttree= PR.EXAML.STARTTREE, prog.rndstarttree=PR.EXAML.BS, args.starttree.type='parsimony', args.starttree.seed=12345, args.starttree.bsid= NA, prog.examl= PR.EXAML.EXAML, args.examl="-m GAMMA -D", resume=0, verbose=1)
{
	#"mpirun -np 4"
	if(is.na(args.starttree.bsid))
		args.starttree.bsid	<- "000"
	else
		args.starttree.bsid	<-	sprintf("%03d",args.starttree.bsid)
	cmd<- "#######################################################
# start: compute ExaML tree
#######################################################\n"	
	#if output files are found and resume, don t do anything
	if(resume)
	{
		cmd		<- paste(cmd,"[ -s ",outdir,'/ExaML_result.',infile,".finaltree.",args.starttree.bsid," ] && ", sep='')	#if file non-zero
		cmd		<- paste(cmd,"[ -s ",outdir,'/ExaML_info.',infile,".finaltree.",args.starttree.bsid," ] ", sep='')		#and if file non-zero
		cmd		<- paste(cmd,"&& exit 1\n",sep='')		#then exit
	}
	#default commands for parser					
	cmd			<- paste(cmd,"CWDEXAML=$(pwd)\n",sep='')
	cmd			<- paste(cmd,"cd ",outdir,'\n',sep='')
	cmd			<- paste(cmd,"echo \'run ",prog.parser,"\'\n",sep='')
	tmp			<- paste(indir,paste(infile,".phylip.",args.starttree.bsid,sep=''),sep='/')
	cmd			<- paste(cmd,prog.parser,' ',args.parser,' -s ',tmp,sep='')
	tmp			<- paste(infile,".phylip.examl.",args.starttree.bsid,sep='')
	cmd			<- paste(cmd," -n ",tmp,sep='')
	#verbose stuff for parser	
	cmd			<- paste(cmd,paste("\necho \'end ",prog.parser,"\'",sep=''))	
	cmd			<- paste(cmd,paste("\necho \'run ",prog.starttree,"\'\n",sep=''))
	#default commands for start tree
	if(args.starttree.type!='random')
	{
		tmp			<- paste(indir,paste(infile,".phylip.",args.starttree.bsid,sep=''),sep='/')	
		cmd			<- paste(cmd,prog.starttree," -p",args.starttree.seed," -s ",tmp,sep='')	
		tmp			<- paste(infile,".starttree.",args.starttree.bsid,sep='')
		cmd			<- paste(cmd," -n ",tmp,sep='')		
	}
	if(args.starttree.type=='random')
	{
		tmp			<- paste(indir,paste(infile,".phylip.",args.starttree.bsid,sep=''),sep='/')	
		cmd			<- paste(cmd,prog.rndstarttree," -y -d -m GTRCAT -p",args.starttree.seed," -s ",tmp,sep='')	
		tmp			<- paste(infile,".starttree.",args.starttree.bsid,sep='')
		cmd			<- paste(cmd," -n ",tmp,sep='')		
	}
	#verbose stuff
	cmd			<- paste(cmd,paste("\necho \'end ",prog.starttree,"\'",sep=''))	
	cmd			<- paste(cmd,paste("\necho \'run ",prog.examl,"\'\n",sep=''))
	#default commands for final tree
	cmd			<- paste(cmd, prog.mpi,' ',prog.examl,' ',args.examl,sep='')
	tmp			<- paste(infile,".phylip.examl.",args.starttree.bsid,".binary",sep='')
	cmd			<- paste(cmd," -s ",tmp,sep='')
	if(args.starttree.type!='random')
		tmp			<- paste("RAxML_parsimonyTree.",infile,".starttree.",args.starttree.bsid, ".0", sep='')
	if(args.starttree.type=='random')
		tmp			<- paste("RAxML_randomTree.",infile,".starttree.",args.starttree.bsid, sep='')	
	cmd			<- paste(cmd," -t ",tmp,sep='')
	tmp			<- paste(infile,".finaltree.",args.starttree.bsid,sep='')
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

#' @export
#' @title Produce a single PhyD* shell command. 
#' @return	Character string
cmd.phydstar<- function(infile.d, infile.v=NA, outfile=NA, pr=PR.PHYD, method='BioNJ', fs=15, binary=TRUE, negative.branch.length=FALSE, lower.triangular=TRUE)
{
	stopifnot(method%in%c('BioNJ','MVR','NJ','UNJ'))	
	stopifnot(negative.branch.length%in%c(TRUE,FALSE))
	stopifnot(lower.triangular%in%c(TRUE,FALSE))
	stopifnot(binary%in%c(TRUE,FALSE))
	
	cmd				<- paste("#######################################################
# start: PhyD*
#######################################################\n",sep='')
	cmd	<- paste(cmd, pr, " -d ",method," -p ",fs," -n ",as.character(factor(negative.branch.length, levels=c(TRUE,FALSE),labels=c('Y','N')))," -b ",as.character(factor(binary, levels=c(TRUE,FALSE),labels=c('Y','N'))),sep='')	
	if(lower.triangular)
		cmd	<- paste(cmd, ' -l',sep='')
	cmd	<- paste(cmd, ' -i ', infile.d, sep='')
	if(!is.na(infile.v))
		cmd	<- paste(cmd, ' -v ', infile.v, sep='')
	if(!is.na(outfile))		
		cmd	<- paste(cmd,'\n','mv ',infile.d,'_',tolower(method),'.t',' ',outfile,'\n',sep='')
	cmd	<- paste(cmd, "#######################################################
# end: PhyD*
#######################################################\n",sep='')	
	cmd
}

#' @export
#' @title Produce a single shell command to run MVR. 
#' @return	Character string
cmd.mvr<- function(infile, outfile, prog=PR.MVR, method='MVR', complete.distance.matrix=FALSE, name.gd.col='GD')
{
	cmd<- "#######################################################
# start: run big.mvr.Rscript 
#######################################################"
	cmd		<- paste(cmd, paste("\necho \'run ",prog,"\'\n",sep=''))
	cmd		<- paste(cmd, paste(prog,' -infile=', infile,' -outfile=',outfile,' -method=',method,' -complete.distance.matrix=',as.integer(complete.distance.matrix),' -name.gd.col=',name.gd.col,' \n', sep=''))
	cmd		<- paste(cmd, paste("echo \'end ",prog,"\'\n",sep=''))
	cmd		<- paste(cmd,"#######################################################
# end: run big.mvr.Rscript
#######################################################\n",sep='')
	cmd
}


#' @export cmd.treedater
#' @title Generate commands to date trees with TreeDater. 
#' @return	Character string
cmd.treedater<- function(infile.tree, infile.dates, outfile, pr=PR.TREEDATER, root=NA, ali.len=NA, omega0=NA)
{
	cmd	<- paste(pr, ' --infile.dates "',infile.dates,'" --infile.tree "',infile.tree,'" --outfile "', outfile,'"', sep='')	
	if(!is.na(ali.len))
		cmd	<- paste(cmd, ' --ali.len ', ali.len, sep='')
	if(!is.na(root))	
		cmd	<- paste(cmd, ' --root "', root, '"', sep='')	
	if(!is.na(omega0))
		cmd	<- paste(cmd, ' --omega0 ', omega0, sep='')
	cmd
}

#' @export
#' @title Generate the LSD Dates file. 
#' @return	Character string
cmd.lsd.dates<- function(infile.dates, infile.tree, outfile.lsd.dates, pr=PR.LSDDATES, run.lsd=FALSE, outfile.lsd=NA, root=NA, exclude.missing.dates=FALSE, outfile.tree=NA, ali.len=NA)
{
	cmd	<- paste(pr, ' --infile.dates "',infile.dates,'" --infile.tree "',infile.tree,'" --outfile.lsd.dates "', outfile.lsd.dates,'"', sep='')	
	if(!is.na(outfile.lsd))
		cmd	<- paste(cmd, ' --outfile.lsd "', outfile.lsd, '"', sep='')
	if(!is.na(ali.len))
		cmd	<- paste(cmd, ' --ali.len ', ali.len, sep='')
	if(!is.na(root))
	{
		stopifnot(!is.na(outfile.tree))
		cmd	<- paste(cmd, ' --root "', root, '"', sep='')	
	}				
	if(run.lsd)
		cmd	<- paste(cmd, ' --run.lsd', sep='')
	if(exclude.missing.dates)
	{
		stopifnot(!is.na(outfile.tree))
		cmd	<- paste(cmd, ' --exclude.missing.dates', sep='')
	}
	if(!is.na(outfile.tree))
		cmd	<- paste(cmd, ' --outfile.tree "', outfile.tree, '" ', sep='')			
	cmd
}


#' @export
#' @title Produce a single LSD shell command. 
#' @return	Character string
cmd.jmodeltest<- function(infile.fasta, outfile=paste0(infile.fasta,'.jmodeltest'), pr=PR.JMODELTEST, pr.args='-f -i -g 4 -s 3 -DT -S NNI -t ML', nproc=1)
{		
	cmd				<- paste("#######################################################
# start: jModelTest
#######################################################\n",sep='')												
	cmd				<- paste(cmd,"CWD=$(pwd)\n",sep='')
	cmd				<- paste(cmd,"echo $CWD\n",sep='')	
	tmpdir.prefix	<- paste('jm_',format(Sys.time(),"%y-%m-%d-%H-%M-%S"),sep='')
	tmpdir			<- paste("$CWD/",tmpdir.prefix,sep='')
	tmp.fasta		<- file.path(tmpdir, basename(infile.fasta))
	tmp.out			<- file.path(tmpdir, 'out')
	cmd				<- paste(cmd,"mkdir -p ",tmpdir,'\n',sep='')
	cmd				<- paste(cmd,'cp "',infile.fasta,'" ',tmpdir,'\n', sep='')	
	cmd				<- paste(cmd, pr,' -d "', tmp.fasta,'" -o "', tmp.out,'" -tr ',nproc, ' ', pr.args,'\n', sep='')
	cmd				<- paste(cmd, "mv ", tmp.out,' "',outfile,'"\n',sep='')
	cmd				<- paste(cmd, "#######################################################
# end: jModelTest
#######################################################\n",sep='')
	cmd
}

#' @export
#' @title Produce a single LSD shell command. 
#' @return	Character string
cmd.lsd<- function(infile.tree, infile.dates, ali.nrow, outfile=infile.tree, pr=PR.LSD, pr.args='-v 2 -c -b 10 -r as')
{		
	cmd				<- paste("#######################################################
# start: LSD
#######################################################\n",sep='')												
	cmd				<- paste(cmd,"CWD=$(pwd)\n",sep='')
	cmd				<- paste(cmd,"echo $CWD\n",sep='')	
	tmpdir.prefix	<- paste('lsd_',format(Sys.time(),"%y-%m-%d-%H-%M-%S"),sep='')
	tmpdir			<- paste("$CWD/",tmpdir.prefix,sep='')
	tmp.tree		<- file.path(tmpdir, basename(infile.tree))
	tmp.dates		<- file.path(tmpdir, basename(infile.dates))
	tmp.out			<- file.path(tmpdir, basename(outfile))	
	cmd				<- paste(cmd,"mkdir -p ",tmpdir,'\n',sep='')
	cmd				<- paste(cmd,'cp "',infile.tree,'" ',tmpdir,'\n', sep='')
	cmd				<- paste(cmd,'cp "',infile.dates,'" ',tmpdir,'\n', sep='')
	cmd				<- paste(cmd, pr,' -i ', tmp.tree,' -d ', tmp.dates,' -s ',ali.nrow, ' -o ', tmp.out, ' ', pr.args,'\n', sep='')
	cmd				<- paste(cmd, "mv ", tmp.out,'* "',dirname(outfile),'"\n',sep='')
	cmd				<- paste(cmd, "#######################################################
# end: LSD
#######################################################\n",sep='')
	cmd
}

#' @export
#' @title Produce a single FastTree shell command. 
#' @return	Character string
cmd.fasttree<- function(infile.fasta, outfile=paste(infile.fasta,'.newick',sep=''), pr=PR.FASTTREE, pr.args='-nt -gtr -gamma')
{		
	cmd				<- paste("#######################################################
# start: FASTTREE
#######################################################\n",sep='')												
	cmd				<- paste(cmd,"CWD=$(pwd)\n",sep='')
	cmd				<- paste(cmd,"echo $CWD\n",sep='')	
	tmpdir.prefix	<- paste('ft_',format(Sys.time(),"%y-%m-%d-%H-%M-%S"),sep='')
	tmpdir			<- paste("$CWD/",tmpdir.prefix,sep='')
	tmp.in			<- file.path(tmpdir, basename(infile.fasta))
	tmp.out			<- file.path(tmpdir, basename(outfile))
	tmp.log			<- file.path(tmpdir, paste(basename(outfile),'.log',sep=''))
	cmd				<- paste(cmd,"mkdir -p ",tmpdir,'\n',sep='')
	cmd				<- paste(cmd,'cp "',infile.fasta,'" ',tmp.in,'\n', sep='')	
	cmd				<- paste(cmd, pr,' ',pr.args,' -log ',tmp.log,' < ', tmp.in,' > ', tmp.out,'\n', sep='')
	cmd				<- paste(cmd, "mv ", tmp.out,'* "',dirname(outfile),'"\n',sep='')
	cmd				<- paste(cmd, "#######################################################
# end: FASTTREE
#######################################################\n",sep='')
	cmd
}

#' @export
#' @title Produce a single ExaML shell command.
#' @description Internal code. 
#' @inheritParams cmd.examl.bootstrap
#' @return	Character string
cmd.examl.single<- function(indir, infile, outdir=indir, outfile=infile, prog.mpi='mpiexec', prog.bscreate=PR.EXAML.BSCREATE, opt.bootstrap.by="nucleotide",prog.parser= PR.EXAML.PARSER, args.parser="-m DNA",prog.starttree= PR.EXAML.STARTTREE, prog.rndstarttree=PR.EXAML.BS, args.starttree.seed=12345, args.starttree.type='parsimony', prog.examl= PR.EXAML.EXAML, args.examl="-m GAMMA -D", bs.seed=floor(runif(1, 1e4, 1e5-1)), verbose=1)
{
	infile			<- gsub('\\.rda|\\.R|\\.fasta|\\.fa','',infile)
	tmpdir.prefix	<- paste('exa_',format(Sys.time(),"%y-%m-%d-%H-%M-%S"),sep='')	
	cmd				<- paste("#######################################################
# start: single examl
#######################################################\n",sep='')												
	cmd				<- paste(cmd,"CWD=$(pwd)\n",sep='')
	cmd				<- paste(cmd,"echo $CWD\n",sep='')
	tmpdir			<- paste("$CWD/",tmpdir.prefix,sep='')
	cmd				<- paste(cmd,"mkdir -p ",tmpdir,'\n',sep='')
	cmd				<- paste(cmd,'cp ',indir,'/',infile,'* ',tmpdir,sep='')
	tmp				<- gsub('-q ','',regmatches(args.parser,regexpr('-q .*',args.parser)))
	if(length(tmp))
		cmd			<- paste(cmd,"\ncp ",indir,'/',tmp," ",tmpdir,sep='')
	cmd				<- paste(cmd,cmd.examl.bsalignment(tmpdir, infile, 0, prog.bscreate=prog.bscreate, opt.bootstrap.by=opt.bootstrap.by, outdir=tmpdir, verbose=verbose),sep='\n')
	cmd				<- paste(cmd,cmd.examl(tmpdir, infile, outdir=tmpdir, prog.mpi=prog.mpi, prog.parser= prog.parser, args.parser=args.parser, prog.starttree= prog.starttree, prog.rndstarttree=prog.rndstarttree, args.starttree.seed=bs.seed, args.starttree.type=args.starttree.type, args.starttree.bsid=0, prog.examl=prog.examl, args.examl=args.examl, resume=0, verbose=verbose),sep='\n')
	cmd				<- paste(cmd,"mv ",tmpdir,"/ExaML_result.",infile,".finaltree.",sprintf("%03d",0),' ', outdir,'/',outfile,'_examl.newick\n',sep='')
	cmd				<- paste(cmd,"mv ",tmpdir,"/ExaML_info.",infile,".finaltree.",sprintf("%03d",0),' ', outdir,'/',outfile,'_examl.txt\n',sep='')
	cmd				<- paste(cmd,"mv ",tmpdir,"/ExaML_modelFile.",infile,".finaltree.",sprintf("%03d",0),' ', outdir,'/',outfile,'_examl.params\n',sep='')
	cmd			<- paste(cmd,"#######################################################
# end: single examl
#######################################################\n",sep='')
	cmd
}

#' @export 
#' @title Produce the shell command that creates the bootstrap alignment.
#' @description Internal code.
#' @inheritParams pipeline.ExaML.bootstrap.per.proc
#' @param opt.bootstrap.by	Character string that specifies specifies the way the boostrap alignment is created. Valid options are \code{codon} and \code{nucleotide}.
#' @return	Character string
cmd.examl.bsalignment<- function(indir, infile, bs.id, outdir=indir, prog.bscreate= PR.EXAML.BSCREATE, opt.bootstrap.by="codon",resume=0, verbose=1)
{
	cmd			<- paste("#######################################################
# start: create and check bootstrap alignment
#######################################################\n",sep='')
	cmd			<- paste(cmd,prog.bscreate," -resume=",resume," -bootstrap=",bs.id," -by=",opt.bootstrap.by,sep='')
	cmd			<- paste(cmd," -indir=",indir," -infile=",infile," -outdir=",outdir,sep='')	
	cmd			<- paste(cmd,"\n#######################################################
# end: create and check bootstrap alignment
#######################################################",sep='')			
	cmd
}

#' @export 
#' @title Produce the shell command that masks resistance mutations.
#' @description Internal code.
#' @inheritParams pipeline.ExaML.bootstrap.per.proc
#' @param alignment.start Number that specifies the position of the first nucleotide relative to the HXB2 reference sequence
#' @return	Character string
cmd.rm.resistance<- function(indir, infile, outfile, outdir=indir, prog= PR.RM.RESISTANCE, verbose=1)
{
	cmd			<- paste("#######################################################
# start: mask codons with resistance mutations with NNN
#######################################################\n",sep='')	
	cmd			<- paste(cmd,prog," -indir=",indir," -infile=",infile," -outdir=",outdir," -outfile=",outfile,' ',sep='')	
	cmd			<- paste(cmd,"\n#######################################################
# end: mask codons with resistance mutations with NNN
#######################################################",sep='')			
	cmd
}


#' @export
#' @title Produce the ExaML boostrap shell command, all boostraps on a single processor.
#' @description Internal code.
#' @inheritParams pipeline.ExaML.bootstrap.per.proc
#' @return	Character string
cmd.examl.bootstrap.on.one.machine<- function(indir, infile, bs.from=0, bs.to=99, bs.n=bs.to-bs.from+ifelse(bs.from==0,1,0), outdir=indir, prog.mpi='mpiexec', prog.parser= PR.EXAML.PARSER, prog.starttree= PR.EXAML.STARTTREE, prog.examl=PR.EXAML.EXAML, opt.bootstrap.by="codon", args.examl="-m GAMMA -D", prog.supportadder=PR.EXAML.BS, tmpdir.prefix="examl", resume=1, verbose=1)
{
	infile			<- gsub('\\.rda|\\.R|\\.fasta|\\.fa','',infile)
	#create number of seeds for the number of runs being processed, which could be less than bs.n
	bs.id			<- seq.int(bs.from,bs.to)
	bs.seeds		<- floor( runif(length(bs.id), 1e4, 1e5-1) )
	tmpdir			<- paste("$CWD/",tmpdir.prefix,'_',format(Sys.time(),"%y-%m-%d-%H-%M-%S"),sep='')
	cmd				<- sapply(seq_along(bs.seeds), function(i)
			{
				cmd			<- ''
				if(i==1)
				{
					cmd	<- paste(cmd,"\nCWD=$(pwd)",sep='')
					cmd	<- paste(cmd,"\necho $CWD",sep='')
					cmd	<- paste(cmd,"\nmkdir -p ",tmpdir,sep='')
					tmp	<- paste(indir,'/',infile,sep='')
					cmd	<- paste(cmd,'\ncp "',tmp,'"* ',tmpdir,sep='')
				}
				cmd		<- paste(cmd,cmd.examl.bsalignment(tmpdir, infile, bs.id[i], opt.bootstrap.by=opt.bootstrap.by, outdir=tmpdir, verbose=verbose),sep='\n')
				cmd		<- paste(cmd,cmd.examl(tmpdir, infile, outdir=tmpdir, prog.mpi=prog.mpi, prog.parser= prog.parser, prog.starttree= prog.starttree, args.starttree.seed=bs.seeds[i], args.starttree.bsid= bs.id[i], prog.examl=prog.examl, args.examl=args.examl, resume=resume, verbose=verbose),sep='\n')				
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
	tmp			<- paste("\nif [ $(find . -name 'ExaML_result.",infile,".finaltree*' | wc -l) -ge ",round(bs.n*0.95)," ]; then",sep='')
	cmd			<- paste(cmd,tmp,sep='')
	cmd			<- paste(cmd,paste("\n\techo \'all bootstrap samples computed -- find best tree and add bootstrap support values\'",sep=''))				
	tmp			<- c(	paste("ExaML_result.",infile,".finaltree",sep=''), 	paste("ExaML_result.",infile,".bstree",sep=''))
	#add all bootstrap final trees into bs tree file
	cmd			<- paste(cmd,"\n\tfor i in $(seq 0 ",bs.n-1,"); do cat ",tmp[1],".$(printf %03d $i) >> ",tmp[2],"; done",sep='')
	#identify suffix finaltree.XXX of final tree with highest likelihood
	cmd			<- paste(cmd,"\n\tBSBEST=$(grep 'Likelihood of best tree' ExaML_info.* | awk '{print $5,$1;}' | sort -n | tail -1 | grep -o 'finaltree.*' | cut -d':' -f 1)",sep='')
	cmd			<- paste(cmd,paste("\n\techo \"best tree is $BSBEST\"",sep=''))		
	#create final tree with bootstrap support values
	tmp			<- c(	paste(infile,".phylip.examl.binary",sep=''),	paste("ExaML_result.",infile,".$BSBEST",sep=''),	paste("ExaML_result.",infile,".bstree",sep=''), paste(infile,".supporttree",sep=''))
	cmd			<- paste(cmd,"\n\t",prog.supportadder," -f b -m GTRCAT -s ",tmp[1]," -t ",tmp[2]," -z ",tmp[3]," -n ",tmp[4],sep='' )
	#create data tree with bootstrap support values
	tmp			<- c(	paste(infile,".phylip.examl.binary",sep=''),	paste("ExaML_result.",infile,".finaltree.000",sep=''),	paste("ExaML_result.",infile,".bstree",sep=''), paste(infile,".datatree",sep=''))	
	cmd			<- paste(cmd,"\n\t",prog.supportadder," -f b -m GTRCAT -s ",tmp[1]," -t ",tmp[2]," -z ",tmp[3]," -n ",tmp[4],sep='' )
	cmd			<- paste(cmd,paste("\n\techo \'all bootstrap samples computed -- found best tree and added bootstrap support values\'",sep=''))										
	#delete ExaML output that is not further needed
	cmd			<- paste(cmd,paste("\n\techo \'start cleanup\'",sep=''))
	cmd			<- paste(cmd,"\n\trm RAxML_info*",sep=' ')								
	cmd			<- paste(cmd,paste("\n\tmv RAxML_bipartitions.",infile,".supporttree ",infile,"_examl_mxbs",bs.n,".newick",sep=''),sep='')
	cmd			<- paste(cmd,paste("\n\tmv RAxML_bipartitions.",infile,".datatree ",infile,"_examl_dtbs",bs.n,".newick",sep=''),sep='')
	cmd			<- paste(cmd,paste("\n\trm RAxML_bipartitionsBranchLabels.",infile,".supporttree",sep=''),sep='')						
	cmd			<- paste(cmd,paste("\n\trm RAxML_bipartitionsBranchLabels.",infile,".datatree",sep=''),sep='')
	cmd			<- paste(cmd,paste("\n\trm ExaML_result.",infile,".bstree",sep=''),sep='')									
	cmd			<- paste(cmd,paste("\n\techo \'end cleanup\'",sep=''))
	#copy ML bstree to outdir
	cmd			<- paste(cmd,paste("\n\tcp -f ",infile,"_examl_mxbs",bs.n,".newick",' "',outdir,'"',sep=''),sep='')
	#copy data bstree to outdir
	cmd			<- paste(cmd,paste("\n\tcp -f ",infile,"_examl_dtbs",bs.n,".newick",' "',outdir,'"',sep=''),sep='')
	cmd			<- paste(cmd,"\nfi",sep='')
	cmd			<- paste(cmd,"\n#######################################################
# end: check if all ExaML boostrap trees have been computed and if yes create ExaML bootstrap tree
#######################################################",sep='')
cmd			<- paste(cmd,"\n#######################################################
# start: zip and copy ExaML output
#######################################################",sep='')	
	#outside if:	zip and copy ExaML output to outdir just in case something went wrong
	cmd			<- paste(cmd,paste("\nzip ",infile,'_examlout',".zip  ExaML_result.",infile,".* ExaML_info.",infile,".* ExaML_modelFile.",infile,".* ",  sep=''),sep='')
	cmd			<- paste(cmd,paste("\ncp -f ",infile,'_examlout',".zip",' "',outdir,'"',sep=''),sep='')
	cmd			<- paste(cmd,paste("\nrm ExaML_result.",infile,".* ExaML_info.",infile,".*",sep=''),sep='')
	cmd			<- paste(cmd,"\ncd $CWD",sep='')
	cmd			<- paste(cmd,"\n#######################################################
# end: zip and copy ExaML output
#######################################################\n",sep='')				
	cmd
}

#' @export
#' @title Generate the ExaML boostrap shell command, one bootstrap per processor.
#' @description Internal code.
#' @inheritParams pipeline.ExaML.bootstrap.per.proc
#' @return	Character string
cmd.examl.bootstrap<- function(indir, infile, bs.from=0, bs.to=99, bs.n=bs.to-bs.from+ifelse(bs.from==0,1,0), outdir=indir, prog.mpi='mpiexec', prog.bscreate=PR.EXAML.BSCREATE, prog.parser= PR.EXAML.PARSER, prog.starttree= PR.EXAML.STARTTREE, prog.examl=PR.EXAML.EXAML, opt.bootstrap.by="codon", args.examl="-m GAMMA -D", args.parser="-m DNA", prog.supportadder=PR.EXAML.BS, tmpdir.prefix="examl", resume=1, verbose=1)
{
	infile			<- gsub('\\.rda|\\.R|\\.fasta|\\.fa','',infile)
	#create number of seeds for the number of runs being processed, which could be less than bs.n
	bs.id			<- seq.int(bs.from,bs.to)
	bs.seeds		<- floor( runif(length(bs.id), 1e4, 1e5-1) )
	tmpdir.prefix	<- paste(tmpdir.prefix,'_',format(Sys.time(),"%y-%m-%d-%H-%M-%S"),sep='')
	lapply(seq_along(bs.seeds), function(i)
			{
				cmd			<- ''				
				if(resume)
				{
					cmd	<- paste(cmd,"\nif ", sep='')
					cmd	<- paste(cmd,"[ ! -s ",outdir,'/ExaML_result.',infile,".finaltree.",sprintf("%03d",bs.id[i])," ]", sep='')
					cmd	<- paste(cmd," || ", sep='')
					cmd	<- paste(cmd,"[ ! -s ",outdir,'/ExaML_info.',infile,".finaltree.",sprintf("%03d",bs.id[i])," ];", sep='')
					cmd	<- paste(cmd," then\n",sep='')
					cmd	<- paste(cmd,"#######################################################	
# start: not indented if statement -- don t do anything if ExaML output exists already
#######################################################",sep='')
				}					
				cmd		<- paste(cmd,"CWD=$(pwd)\n",sep='\n')
				cmd		<- paste(cmd,"echo $CWD\n",sep='')
				tmpdir	<- paste("$CWD/",tmpdir.prefix,sep='')
				cmd		<- paste(cmd,"mkdir -p ",tmpdir,'\n',sep='')
				tmp		<- paste(indir,'/',infile,".*",sep='')
				cmd		<- paste(cmd,"cp ",tmp," ",tmpdir,sep='')
				tmp		<- gsub('-q ','',regmatches(args.parser,regexpr('-q .*',args.parser)))
				if(length(tmp))
				{
					cmd	<- paste(cmd,"\ncp ",indir,'/',tmp," ",tmpdir,sep='')
				}
				cmd		<- paste(cmd,cmd.examl.bsalignment(tmpdir, infile, bs.id[i], prog.bscreate=prog.bscreate, opt.bootstrap.by=opt.bootstrap.by, outdir=tmpdir, verbose=verbose),sep='\n')
				cmd		<- paste(cmd,cmd.examl(tmpdir, infile, outdir=tmpdir, prog.mpi=prog.mpi, prog.parser= prog.parser, args.parser=args.parser, prog.starttree= prog.starttree, args.starttree.seed=bs.seeds[i], args.starttree.bsid= bs.id[i], prog.examl=prog.examl, args.examl=args.examl, resume=resume, verbose=verbose),sep='\n')
				cmd		<- paste(cmd,"cp -f ",tmpdir,"/ExaML_result.",infile,".finaltree.",sprintf("%03d",bs.id[i]),' ', outdir,'\n',sep='')
				cmd		<- paste(cmd,"cp -f ",tmpdir,"/ExaML_info.",infile,".finaltree.",sprintf("%03d",bs.id[i]),' ', outdir,'\n',sep='')
				cmd		<- paste(cmd,"cp -f ",tmpdir,"/ExaML_modelFile.",infile,".finaltree.",sprintf("%03d",bs.id[i]),' ', outdir,'\n',sep='')	
				if(resume)
				{
						cmd		<- paste(cmd,"#######################################################
# end: not indented if statement -- don t do anything if ExaML output exists already
#######################################################\nelse\n\techo 'resumed bootstrap number ",bs.id[i],"'\nfi\n",sep='')					
				}
								
				cmd			<- paste(cmd,"#######################################################
# start: check if all ExaML boostrap trees have been computed and if yes create ExaML bootstrap tree
#######################################################",sep='')	
				cmd			<- paste(cmd,"\nCWD=$(pwd)\n",sep='')
				cmd			<- paste(cmd,"cd ",outdir,sep='') 
				cmd			<- paste(cmd,paste("\necho \'check if all bootstrap samples have been computed\'",sep=''))
				tmp			<- paste("\nif [ $(find . -name 'ExaML_result.",infile,".finaltree*' | wc -l) -eq ",bs.n," ]; then",sep='')
				cmd			<- paste(cmd,tmp,sep='')
				cmd			<- paste(cmd,paste("\n\techo \'all bootstrap samples computed -- find best tree and add bootstrap support values\'",sep=''))				
				tmp			<- c(	paste("ExaML_result.",infile,".finaltree",sep=''), 	paste("ExaML_result.",infile,".bstree",sep=''))
				#add all bootstrap final trees into bs tree file
				cmd			<- paste(cmd,"\n\tfor i in $(seq 0 ",bs.n-1,"); do cat ",tmp[1],".$(printf %03d $i) >> ",tmp[2],"; done",sep='')
				#identify suffix finaltree.XXX of final tree with highest likelihood
				cmd			<- paste(cmd,"\n\tBSBEST=$(grep 'Likelihood of best tree' ExaML_info.* | awk '{print $5,$1;}' | sort -n | tail -1 | grep -o 'finaltree.*' | cut -d':' -f 1)",sep='')
				cmd			<- paste(cmd,paste("\n\techo \"best tree is $BSBEST\"",sep=''))		
				#create final tree with bootstrap support values
				tmp			<- c(	paste(infile,".phylip.examl.binary",sep=''),	paste("ExaML_result.",infile,".$BSBEST",sep=''),	paste("ExaML_result.",infile,".bstree",sep=''), paste(infile,".supporttree",sep=''))
				cmd			<- paste(cmd,"\n\t",prog.supportadder," -f b -m GTRCAT -s ",tmp[1]," -t ",tmp[2]," -z ",tmp[3]," -n ",tmp[4],sep='' )
				#create data tree with bootstrap support values
				tmp			<- c(	paste(infile,".phylip.examl.binary",sep=''),	paste("ExaML_result.",infile,".finaltree.000",sep=''),	paste("ExaML_result.",infile,".bstree",sep=''), paste(infile,".datatree",sep=''))	
				cmd			<- paste(cmd,"\n\t",prog.supportadder," -f b -m GTRCAT -s ",tmp[1]," -t ",tmp[2]," -z ",tmp[3]," -n ",tmp[4],sep='' )				
				cmd			<- paste(cmd,paste("\n\techo \'all bootstrap samples computed -- found best tree and added bootstrap support values\'",sep=''))										
				#delete ExaML output that is not further needed
				cmd			<- paste(cmd,paste("\n\techo \'start cleanup\'",sep=''))
				cmd			<- paste(cmd,"\n\trm RAxML_info*",sep=' ')								
				cmd			<- paste(cmd,paste("\n\tmv RAxML_bipartitions.",infile,".supporttree ",infile,"_examl_mxbs",bs.n,".newick",sep=''),sep='')
				cmd			<- paste(cmd,paste("\n\tmv RAxML_bipartitions.",infile,".datatree ",infile,"_examl_dtbs",bs.n,".newick",sep=''),sep='')
				cmd			<- paste(cmd,paste("\n\trm RAxML_bipartitionsBranchLabels.",infile,".supporttree",sep=''),sep='')						
				cmd			<- paste(cmd,paste("\n\trm RAxML_bipartitionsBranchLabels.",infile,".datatree",sep=''),sep='')
				cmd			<- paste(cmd,paste("\n\trm ExaML_result.",infile,".bstree",sep=''),sep='')									
				cmd			<- paste(cmd,paste("\n\techo \'end cleanup\'",sep=''))				
				cmd			<- paste(cmd,"\nfi",sep='')
				cmd			<- paste(cmd,"\ncd $CWD",sep='')
				cmd			<- paste(cmd,"\n#######################################################
# end: check if all ExaML boostrap trees have been computed and if yes create ExaML bootstrap tree
#######################################################\n",sep='')				
			})
}

#' @export
#' @title Produce a shell command to compute the final boostrap tree.
#' @description Internal code.
#' @return Character string.
cmd.examl.bsstarttree<- function(indir, infile, bs.from=0, bs.to=99, bs.n=bs.to-bs.from+ifelse(bs.from==0,1,0),outdir=indir, prog.mpi='mpiexec', prog.parser= PR.EXAML.PARSER, prog.starttree= PR.EXAML.STARTTREE, prog.examl=PR.EXAML.EXAML, args.examl="-m GAMMA -D", prog.supportadder=PR.EXAML.BS, resume=1, verbose=1)
{
	#create number of seeds for the number of runs being processed, which could be less than bs.n
	bs.id	<- seq.int(bs.from,bs.to)
	bs.seeds<- floor( runif(length(bs.id), 1e4, 1e5-1) )
	lapply(seq_along(bs.seeds), function(i)
			{
				cmd			<- paste("cp ",paste(indir,paste(infile,".phylip",sep=''),sep='/'),sep='')
				cmd			<- paste(cmd,' ',paste(indir,paste(infile,".phylip.",sprintf("%03d",bs.id[i]),sep=''),sep='/'),"\n",sep='')
				cmd			<- paste(cmd, cmd.examl(indir, infile, outdir=outdir, prog.mpi=prog.mpi, prog.parser= prog.parser, prog.starttree= prog.starttree, args.starttree.seed=bs.seeds[i], args.starttree.bsid= bs.id[i], prog.examl=prog.examl, args.examl=args.examl, resume=resume, verbose=verbose),sep='')
				curr.dir	<- getwd()
				cmd			<- paste(cmd,"#######################################################
# start: check if all ExaML boostrap trees have been computed and if yes create ExaML bootstrap tree
#######################################################",sep='')				
				cmd			<- paste(cmd,"\ncd ",outdir,sep='')
				cmd			<- paste(cmd,paste("\necho \'check if all bootstrap samples have been computed\'",sep=''))
				tmp			<- paste("\nif [ $(find . -name 'ExaML_result*' | wc -l) -eq ",bs.n," ]; then",sep='')
				cmd			<- paste(cmd,tmp,sep='')
				cmd			<- paste(cmd,paste("\n\techo \'all bootstrap samples computed -- find best tree and add bootstrap support values\'",sep=''))				
				tmp			<- c(	paste("ExaML_result.",infile,".finaltree",sep=''), 	paste("ExaML_result.",infile,".bstree",sep=''))
				#add all bootstrap final trees into bs tree file
				cmd			<- paste(cmd,"\n\tfor i in $(seq 0 ",bs.n-1,"); do cat ",tmp[1],".$(printf %03d $i) >> ",tmp[2],"; done",sep='')
				#identify suffix finaltree.XXX of final tree with highest likelihood
				cmd			<- paste(cmd,"\n\tBSBEST=$(grep 'Likelihood of best tree' ExaML_info.* | awk '{print $5,$1;}' | sort -n | tail -1 | grep -o 'finaltree.*' | cut -d':' -f 1)",sep='')
				cmd			<- paste(cmd,paste("\n\techo \"best tree is $BSBEST\"",sep=''))		
				#create final tree with bootstrap support values
				tmp			<- c(	paste(infile,".phylip.examl.binary",sep=''),	paste("ExaML_result.",infile,".$BSBEST",sep=''),	paste("ExaML_result.",infile,".bstree",sep=''), paste(infile,".supporttree",sep=''))
				cmd			<- paste(cmd,"\n\t",prog.supportadder," -f b -m GTRCAT -s ",tmp[1]," -t ",tmp[2]," -z ",tmp[3]," -n ",tmp[4],sep='' )
				cmd			<- paste(cmd,paste("\n\techo \'all bootstrap samples computed -- found best tree and added bootstrap support values\'",sep=''))										
				#delete ExaML output that is not further needed
				cmd			<- paste(cmd,paste("\n\techo \'start cleanup\'",sep=''))
				cmd			<- paste(cmd,"\n\trm RAxML_info*",sep=' ')
				cmd			<- paste(cmd,paste("\n\trm ",infile,".phylip.examl.binary",sep=''),sep='')				
				cmd			<- paste(cmd,paste("\n\tmv RAxML_bipartitions.",infile,".supporttree ",infile,"_examlbs",bs.n,".newick",sep=''),sep='')				
				cmd			<- paste(cmd,paste("\n\trm RAxML_bipartitionsBranchLabels.",infile,".supporttree",sep=''),sep='')
				cmd			<- paste(cmd,paste("\n\trm ExaML_result.",infile,".bstree",sep=''),sep='')								
				cmd			<- paste(cmd,paste("\n\ttar -zcf ",infile,'_examlout',".tar.gz  ExaML_result.",infile,".* ExaML_info.",infile,".*",sep=''),sep='')
				cmd			<- paste(cmd,paste("\n\trm ExaML_result.",infile,".* ExaML_info.",infile,".*",sep=''),sep='')
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
#' @title Procude a shell command to strip gaps from alignment.
#' @description Internal code.
#' @return Character string
cmd.strip.gaps	<- function(infile, outfile=infile, strip.max.len=NA, strip.pc=NA, prog=PR.STRIP.GAPS)
{
	cmd		<- "#######################################################
# start: strip.gaps
#######################################################\n"
	cmd		<- paste(cmd, prog, ' -infile=',infile,' ', ' -outfile=',outfile,' ',sep='')
	if(!is.na(strip.max.len))
		cmd	<- paste(cmd, ' -strip.max.len=',strip.max.len, sep='')
	if(!is.na(strip.pc))
		cmd	<- paste(cmd, ' -strip.pc=',strip.pc, sep='')
	cmd		<- paste(cmd,"\n#######################################################
# end: strip.gaps
#######################################################\n",sep='')
	cmd
}

#' @export
#' @title Procude a shell command to clean up after an ExaML run.
#' @description Internal code.
#' @return Character string
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

#' @export
cmd.hpcwrapper.cx1.ic.ac.uk<- function(hpc.walltime=24, hpc.mem=HPC.MEM, hpc.nproc=1, hpc.q=NA, hpc.load=HPC.CX1.IMPERIAL.LOAD)
{
	wrap<- "#!/bin/sh"
	tmp	<- paste("#PBS -l walltime=",hpc.walltime,":59:59,pcput=",hpc.walltime,":45:00",sep='')
	wrap<- paste(wrap, tmp, sep='\n')		
	tmp	<- paste("#PBS -l select=1:ncpus=",hpc.nproc,":mem=",hpc.mem,sep='')
	wrap<- paste(wrap, tmp, sep='\n')
	wrap<- paste(wrap, "#PBS -j oe", sep='\n')
	if(!is.na(hpc.q))
		wrap<- paste(wrap, paste("#PBS -q",hpc.q), sep='\n\n')
	wrap<- paste(wrap, hpc.load, sep='\n')
	wrap
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