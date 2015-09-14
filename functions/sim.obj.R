###############################################################
#
#	createRandomClusters:
#		ni:		Vector of size 1*k, containing the number
#				of observations per each group (the i-th's)
#		p:		Number of dimensions
#		k:		Number of clusters
#		cent:	Matrix (k*p) of centres
#
#	return value:
#		SimObj containing:
#		n:				Number of observations
#		p:				Number of dimensions (input p)
#		real.k:			Number f clusters (input k)
#		real.cent: 		Matrix of centres (input cent)
#		dat:			Matrix of simulated data (n*p)
#		real.assign:	Matrix of indicator assignments
#						(cell i,j=1 if observation i
#						belongs to cluster j) - n*k
#
###############################################################
createRandomClusters <- function(ni, p, k, cent, sdi)
{
	# Create the object containing simulation
	simObj = NULL;
	n = sum(ni);
	simObj$n = n;
	simObj$p = p;
	simObj$real.k = k;
	simObj$real.cent = cent;

	# Create the skeleton for the data and assignment matrix
	dat = matrix( rep(0, n*p), nrow= n, ncol= p );
	assgn = matrix( rep(0, n*k), nrow= n, ncol= k );

	hi = 1;
	for (i in 1:k)
	{
		curr.ni = ni[i];	# num. observations in i-th cluster
		hi.from = hi;
		hi.to = hi.from + curr.ni - 1;

		for (j in 1:p)
		{
			curr.cent = cent[i, j];
			# Hypothesis - each cluster comes from
			# different population
			x = rnorm( curr.ni, mean= curr.cent, sd= sdi );
			dat[hi.from:hi.to, j] = x;
		}

		assgn[hi.from:hi.to, i] = 1;
		hi = hi + curr.ni;
	}

	# Shuffle them (both data and grouping vector)
	prmt = sample(1:n, replace=FALSE);
	simObj$dat = dat[prmt,];
	simObj$real.assign = matrix( assgn[prmt,], nrow=n );

	return (simObj);
}



###############################################################
#
#	createRandomCenters:
#		p:		Number of dimensions
#		k:		Number of clusters
#		vec.min:	Vector of size 1*p for lower bound
#		vec.max:	Vector of size 1*p for upper bound
#
#	return value:
#		Matrix of size k*p of random number between
#		[vec.min:vec.max]
#
###############################################################
createRandomCenters <- function(p, k, vec.min, vec.max)
{
	cent = matrix( rep(0, k*p), nrow= k, ncol= p);
	for (j in 1:p)
	{
		cent[,j] = runif(k, min= vec.min[j], max= vec.max[j]);
	}

	return ( cent );
}



###############################################################
#
#	write.simObject:
#		simObj:	Object of type SimObj (as created
#				via createRandomClusters)
#		filename:	name of file to save to
#
#	return value:
#		None. Creates "filename" containing SimObj for
#		later load (via read.simObject)
#
###############################################################
write.simObject <- function( simObj, filename )
{
	write( simObj$n, file=filename, append=FALSE );
	write( simObj$p, file=filename, append=TRUE );
	write( simObj$real.k, file=filename, append=TRUE );
	write.table( simObj$real.cent, file=filename, append=TRUE, row.name=FALSE, col.name=FALSE, sep="\t" );
	write.table( simObj$dat, file=filename, append=TRUE, row.name=FALSE, col.name=FALSE, sep="\t" );
	write.table( simObj$real.assign, file=filename, append=TRUE, row.name=FALSE, col.name=FALSE, sep="\t" );
}



###############################################################
#
#	read.simObject:
#		filename:	Name of file to be loaded. The file
#				must have been created via 
#				write.simObject
#
#	return value:
#		simObject:	The reloaded object created beforehand
#				via write.simObject.
#
###############################################################
read.simObject <- function( filename )
{
	ww = readLines( filename);
	o = NULL;

	idx = 1;
	n = as.numeric( ww[idx] );
	o$n = n;

	idx = 2;
	p = as.numeric( ww[idx] );
	o$p = p;

	idx = 3;
	real.k = as.numeric( ww[idx] );
	o$real.k = real.k;

	idx = (idx+1):(idx+real.k);
	real.cent = unlist( strsplit( ww[idx], "\t" ) );
	real.cent = as.numeric( real.cent );
	o$real.cent = t( matrix( real.cent, ncol=real.k, nrow=p ) );

	idx = max(idx);
	idx = (idx+1):(idx+n);
	dat = unlist( strsplit( ww[idx], "\t" ) );
	dat = as.numeric( dat );
	o$dat = t( matrix( dat, ncol=n, nrow=p ) );

	idx = max(idx);
	idx = (idx+1):(idx+n);
	real.assign = unlist( strsplit( ww[idx], "\t" ) );
	real.assign = as.numeric( real.assign );
	o$real.assign = t( matrix( real.assign, ncol=n, nrow=real.k) );

	return( o );
}



###############################################################
#
#	draw.clustering.state:
#		dat:		Data matrix
#		cent:	Center matrix
#		X:		Assignment matrix
#		ttl:		Title of the plot
#		method:	Colors or Shapes
#
#	activity:
#		Plots the centers and the clusters
#
#	return value:
#		None - A procedure
#
###############################################################
draw.clustering.state <- function( dat, cent, X, ttl, method="Shapes", lims="diff")
{
	k = dim(X)[2];

	if (lims == "diff")
	{
		x.lim = c( min(dat[,1]), max(dat[,1]) );
		y.lim = c( min(dat[,2]), max(dat[,2]) );
	}
	else if (lims == "Sqaures")
	{
		x.lim = c( min(dat), max(dat) );
		y.lim = x.lim;
	}

	grp = X %*% 1:k;
	plot( cent, col="yellow", pch=19, cex=5, xlim=x.lim, ylim=y.lim, xlab="", ylab="", main=ttl);

	if (method == "Colors")
	{
		points(dat, col=grp, pch=19);
	}
	else if (method == "Shapes")
	{
		points(dat, col="black", pch=grp);
	}
}
