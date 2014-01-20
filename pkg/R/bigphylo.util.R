###############################################################################
#	some utility functions
.ls.objects <- function (pos = 1, pattern, order.by,
		decreasing=FALSE, head=FALSE, n=5) {
	napply <- function(names, fn) sapply(names, function(x)
					fn(get(x, pos = pos)))
	names <- ls(pos = pos, pattern = pattern)
	obj.class <- napply(names, function(x) as.character(class(x))[1])
	obj.mode <- napply(names, mode)
	obj.type <- ifelse(is.na(obj.class), obj.mode, obj.class)
	obj.prettysize <- napply(names, function(x) {
				capture.output(print(object.size(x), units = "auto")) })
	obj.size <- napply(names, object.size)
	obj.dim <- t(napply(names, function(x)
						as.numeric(dim(x))[1:2]))
	vec <- is.na(obj.dim)[, 1] & (obj.type != "function")
	obj.dim[vec, 1] <- napply(names, length)[vec]
	out <- data.frame(obj.type, obj.size, obj.prettysize, obj.dim)
	names(out) <- c("Type", "Size", "PrettySize", "Rows", "Columns")
	if (!missing(order.by))
		out <- out[order(out[[order.by]], decreasing=decreasing), ]
	if (head)
		out <- head(out, n)
	out
}

# from: http://stackoverflow.com/questions/1358003/tricks-to-manage-the-available-memory-in-an-r-session
lsos <- function(..., n=10) {
	.ls.objects(..., order.by="Size", decreasing=TRUE, head=TRUE, n=n)
}

my.mkdir<-function(root,data.name)
{
	if(length(dir(root,pattern=paste('^',data.name,'$',sep='')))==0)
		system(paste("mkdir ",paste(root,data.name,sep='/'),sep=''))
}

package.roxygenize<- function()
{
	require(roxygen2)		
	roxygenize(CODE.HOME)
}

package.generate.rdafiles<- function()
{
	require(ape)
	files	<- c('mtDNA','den2','neisseria'	)
	dummy	<- lapply(files, function(x)
			{
				cat(paste('\nprocess',x))
				file			<- paste(CODE.HOME,'/data/',x,'.phylip',sep='')				
				seq				<- read.dna(file, format='sequential')
				rownames(seq)	<- gsub('\\s','',rownames(seq))
				file			<- paste(CODE.HOME,'/data/',x,'.rda',sep='')
				cat(paste('\nsave seq to file=',file))
				save(seq, file=file)				
			})		
	files	<- c('nz_h1n1','nz_h3n2')
	dummy	<- lapply(files, function(x)
			{
				cat(paste('\nprocess',x))
				file			<- paste(CODE.HOME,'/data/',x,'.phylip',sep='')				
				seq				<- read.dna(file, format='interleaved')
				rownames(seq)	<- gsub('\\s','',rownames(seq))
				file			<- paste(CODE.HOME,'/data/',x,'.rda',sep='')
				cat(paste('\nsave seq to file=',file))
				save(seq, file=file)				
			})
}

my.fade.col<-function(col,alpha=0.5)
{
	return(rgb(col2rgb(col)[1]/255,col2rgb(col)[2]/255,col2rgb(col)[3]/255,alpha))
}

my.dumpframes<- function()
{
	geterrmessage()
	dump.frames()
	cat(paste("\nmy.dumpframes dump 'last.dump' to file",paste(DATA,paste("debug_",paste(strsplit(date(),' ')[[1]],collapse='_'),".rda\n",sep=''),sep='/')))
	save(last.dump, file=paste(DATA,paste("debug_",paste(strsplit(date(),' ')[[1]],collapse='_'),".rda",sep=''),sep='/'))
	q()
}

my.aggregate<- function(x, bins)
{
	bins	<- cbind(bins[-length(bins)],bins[-1])			
	ans		<- numeric(length(x))
	for(i in seq_len(nrow(bins)))
		ans[which(x<bins[i,2] & x>=bins[i,1])]<- i
	factor(ans, levels=seq_len(nrow(bins)), labels=apply(bins,1,function(z) paste(z,collapse=',')))		
}

my.intersect.n<- function( list.to.intersect )
{
	if(length(list.to.intersect)==2)
		ans<- intersect(list.to.intersect[[1]],list.to.intersect[[2]])
	else
	{
		ans	<- table(unlist(list.to.intersect))														
		ans <- as.numeric(attr(ans,"dimnames")[[1]])[ans==length(list.to.intersect)]	
		
	}	
	ans
}

print.v<- function(x,cut=3,digits=4,prefix= "simu_",print.char= TRUE, as.R= TRUE)
{
	if(as.R)
	{
		tmp<- paste("c(",paste(c(x,recursive=T),collapse=',',sep=''),')',sep='')
		if(!is.null(names(x)))
			tmp<- paste("{tmp<-", tmp, "; names(tmp)<- ", paste('c("',paste(c(names(x),recursive=T),collapse='", "',sep=''),'")',sep=''), "; tmp}", sep= '', collapse= '')
	}
	else
	{
		if(!is.null(names(x)))
		{
			m<- matrix(NA,nrow=2,ncol=length(x))
			m[1,]<- substr(names(x),1,cut)
			m[2,]<- round( x, digits=digits )
			if(cut==0)		m<- m[2,]
			tmp<- gsub('.',',',paste(prefix,paste(as.vector(m), collapse='_',sep=''),sep=''),fixed=T)
		}
		else
			tmp<- gsub('.',',',paste(prefix,paste(round( x, digits=digits ), collapse='_',sep=''),sep=''),fixed=T)
	}
	if(print.char) print(tmp)
	tmp
}

my.barplot.table<- function(x, xh, xlab, ylab, cols, x.baseline=0, xax=as.numeric(colnames(x)), legend.loc=NULL)
{
	z		<- rbind( rep(x.baseline, ncol(x)), apply(x, 2, cumsum) )
	xlim	<- range(xax)
	ylim	<- range(z)
	
	plot(1,1,type='n',bty='n',xlab=xlab,ylab=ylab, xlim=xlim, ylim=ylim)
	lapply(2:nrow(z), function(i)
			{
				lapply(seq_len(ncol(z)), function(j)
						{
							polygon(	c(xax[j]-xh,xax[j]+xh,xax[j]+xh,xax[j]-xh),		c(rep(z[i-1,j],2),rep(z[i,j],2)), border=NA, col=cols[i-1] 	)
						})
			} )
	if(!is.null(legend.loc))
		legend(legend.loc, fill= cols, legend=rownames(x), border=NA, bty='n')
}

plot.2D.dens<- function(x,y,xlab,ylab,xlim=NA,ylim=NA,nbin=NA,width.infl=2,n.hists=5,method="gauss", palette= "topo", persp.theta= -30, persp.phi= 30, zero.abline=TRUE, ...)
{
	if(!method%in%c("gauss","ash","persp"))	stop("plot.2D.dens: exception 1a")
	if(!palette%in%c("topo","heat","gray"))	stop("plot.2D.dens: exception 1b")
	switch(method,
			gauss=
					{
						require(KernSmooth)
						require(fields)
						x.bw<- width.infl*diff(summary(x)[c(2,5)])
						y.bw<- width.infl*diff(summary(y)[c(2,5)])
						if(!x.bw) x.bw<- EPS
						if(!y.bw) y.bw<- EPS
						if(any(is.na(xlim)))	xlim<- range(x)+c(-1.5,1.5)*x.bw
						if(any(is.na(ylim)))	ylim<- range(y)+c(-1.5,1.5)*y.bw
						dens <- bkde2D(cbind(x, y), range.x=list(xlim,ylim),bandwidth=c(x.bw,y.bw))
						contour(dens$x1, dens$x2, dens$fhat,xlab=xlab,ylab=ylab)
						if(zero.abline) abline(v=0,col="black",lty=3,lwd=1.5)
						if(zero.abline) abline(h=0,col="black",lty=3,lwd=1.5)
					},
			ash={
				require(ash)
				if(any(is.na(xlim))) xlim<- range(x)*1.05
				if(any(is.na(ylim))) ylim<- range(y)*1.05
				if(any(is.na(nbin))) nbin<- 2*c(nclass.Sturges(x),nclass.Sturges(y))
				bins<- bin2(cbind(x, y), ab=rbind(xlim,ylim),nbin=nbin)
				f <- ash2(bins,rep(n.hists,2))
				#image(f$x,f$y,f$z, col=rainbow(50,start= 4/6,end=0),xlab=xlab,ylab=ylab)
				if(palette=="topo")					image(f$x,f$y,f$z, col=tail( topo.colors(trunc(50*1.4)), 50 ),xlab=xlab,ylab=ylab,...)
				else if(palette=="gray")			image(f$x,f$y,f$z, col=head( rev(gray(seq(0,.95,len=trunc(50*1.4)))), 50),xlab=xlab,ylab=ylab,...)					
				else								image(f$x,f$y,f$z, col=heat.colors( 50 ),xlab=xlab,ylab=ylab,...)
				contour(f$x,f$y,f$z,add=TRUE, nlevels= 5)
				if(zero.abline) abline(v=0,col="black",lty=2,lwd=2)
				if(zero.abline) abline(h=0,col="black",lty=2,lwd=2)
			},
			persp={
				require(ash)
				require(MASS)
				if(any(is.na(xlim))) xlim<- range(x)*1.05
				if(any(is.na(ylim))) ylim<- range(y)*1.05
				bins<- bin2(cbind(x, y), ab=rbind(xlim,ylim),nbin=2*c(nclass.Sturges(x),nclass.Sturges(y)))
				f <- ash2(bins,rep(n.hists,2))
				
				nrz <- nrow(f$z)
				ncz <- ncol(f$z)
				col<- tail(topo.colors(trunc(1.4 * 50)),50)
				fcol      <- col[trunc(f$z / max(f$z)*(50-1))+1]
				dim(fcol) <- c(nrz,ncz)
				fcol      <- fcol[-nrz,-ncz]
				par(mar=c(1/2,1/2,1/2,1/2))
				persp(x=f$x,y=f$y,z=f$z,col= fcol,theta=persp.theta,phi=persp.phi,xlab=xlab,ylab=ylab,zlab='', ticktype= "detailed" )
			})
}

plot.persplocfit<- function(x, pv, theta= 30, phi= 20, palette= "gray",tcl=-0.05,...)
{	
	d <- x$mi["d"]
	ev <- x$mi["ev"]
	where <- "grid"
	
	pv <- match(pv, x$vnames)
	tv <- (1:d)[-pv]
	vrs <- c(pv, tv)
	if (any(duplicated(vrs))) 
		warning("Duplicated variables in pv, tv")
	if (any((vrs <= 0) | (vrs > d))) 
		stop("Invalid variable numbers in pv, tv")
	m <- ifelse(d == 1, 100, 40) 
	m <- rep(m, d)
	m[tv] <- mtv<- 6		
	xl <- x$box
	marg <- lfmarg(xl, m)
	pred <- locfit:::preplot.locfit(x, marg, band = "none", tr = NULL, what = "coef", get.data = 0, f3d = 0)
	z<- matrix(pred$fit, nrow = length(marg[[1]]))
	nbcol <- 100
	if(palette=="gray")	
		color<- head( rev(gray(seq(0,0.95,len=trunc(nbcol*1.4)))), nbcol)
	else
	{
		jet.colors <- colorRampPalette( c("blue", "green") )
		color <- jet.colors(nbcol)
	}
	# Compute the z-value at the facet centres
	nrz <- nrow(z)
	ncz <- ncol(z)
	zfacet <- z[-1, -1] + z[-1, -ncz] + z[-nrz, -1] + z[-nrz, -ncz]
	# Recode facet z-values into color indices
	facetcol <- cut(zfacet, nbcol)		
	par(mar=c(0,2.5,0,0), tcl=tcl)
	pmat<- persp(marg[[1]], marg[[2]], z, zlim= range(z)*1.1, col = color[facetcol], shade= 0.1, border=NA, ticktype = "detailed", ltheta = 120,theta = theta, phi = phi, expand = 0.75, box=1, ... )
	
	list(pmat=pmat, x= marg[[1]], y= marg[[2]], z= z)
}
