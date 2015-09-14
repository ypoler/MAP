library(mclust);

###############################################################
#
#	calculateRSS
#
###############################################################
calculateRSS <- function(n, ni, p, k, X, cent, dat, method)
{
	#########################################################
	#
	#	distanceRegular:
	#		parameters are equal to main's without 
	#		"method"
	#
	#	return value:
	#		basic Euclidean distance between each
	#		observations and given center.
	#
	#########################################################
	distanceRegular <- function(n, ni, p, k, X, cent, dat)
	{
		cnt = X %*% cent;
		dst = dat-cnt;
		dst = dst^2;
		rss = sum( dst )
		return ( rss );
	}

	#########################################################
	#
	#	distanceMAP:
	#		parameters are equal to main's without 
	#		"method"
	#
	#	return value:
	#		MAP's first two term - log(Wk) +
	#		1/n*sum( log(ni) )
	#
	#########################################################
	distanceMAP <- function(n, ni, p, k, X, cent, dat)
	{
		Wk = distanceRegular(n, ni, p, k, X, cent, dat);
		lnNj = log(ni);
		return ( log(Wk) + sum(lnNj)/n );
	}

	
	#########################################################
	#
	#	distanceWeight:
	#		parameters are equal to main's without 
	#		"method"
	#
	#	return value:
	#		distance including the ni weight per cluster.
	#
	#########################################################
	distanceWeight <- function(n, ni, p, k, X, cent, dat)
	{
		rss = 0;
		cnt = X %*% cent;
		dst = dat-cnt;
		dst = dst^2;
		for (i in 1:k)
		{
			idx = (X[,i] == 1)/n;
			rssi = sum( dst[idx, ] )*sum(idx);
			rss = rss + rssi;
		}


		return (rss);
	}

	# Chose between methods
	rss = switch (method, 
			"Regular" = distanceRegular(n, ni, p, k, X, cent, dat),
			"MAP" = distanceMAP(n, ni, p, k, X, cent, dat),
			"Weight" = distanceWeight(n, ni, p, k, X, cent, dat)
		  );

	return (rss);
}



###############################################################
#
#	calculateNumberOfMembers
#
###############################################################
calculateNumberOfMembers <- function(X)
{
	return( apply( X, 2, sum ) );
}



###############################################################
#
#	calculateCenters
#
###############################################################
calculateCenters <- function(ni, X, dat)
{
	sum.cluster = t(dat) %*% X;		# sum of all observations in each cluster
	new.cent = t(sum.cluster)/ni;	# mean for each cluster

	return (new.cent);
}



###############################################################
#
#	calculateAssignment
#
###############################################################
calculateAssignment <- function(n, p, k, cent, dat)
{
	# Calculate distance between each observation and each center
	dist = NULL;
	for (j in 1:k)
	{
		curr.cent = cent[j,];
		clust.diff = apply( dat, 1, "-", curr.cent );	# diff between each observation and center j
		clust.dist = apply( clust.diff^2, 2, sum );		# euclidean distance
		dist = cbind( dist, clust.dist );
	}

	# create assignment matrix as minimum of distances
	X = matrix( rep(0, n*k), nrow=n, ncol=k );
	for (i in 1:n)
	{
		idx = which.min( dist[i,] );
		X[i,idx] = 1;
	}

	return (X);
}



###############################################################
#
#	findClusters
#
###############################################################
findClusters <- function(n, p, k, init.X, init.cent, dat, method)
{
	retObj = NULL;

	# Set initialized members as first step
	init.cent = matrix( init.cent, ncol=p);
	curr.cent = init.cent;
	curr.X = init.X;

	# Check if all clusters have members
	curr.ni = calculateNumberOfMembers(curr.X);
	if ( min(curr.ni) == 0)
	{
		return ( NULL );
	}

	# Calculate RSS and initialize the
	curr.RSS = calculateRSS(n, curr.ni, p, k, curr.X, curr.cent, dat, method);

	# Loop and recalculate until converge
	iter = 1;
	cont = TRUE;
	while (cont)
	{
		cont = FALSE;
		new.cent = calculateCenters(curr.ni, curr.X, dat);
		new.X = calculateAssignment(n, p, k, curr.cent, dat);

		# Check that all clusters have members
		new.ni = calculateNumberOfMembers(new.X);
		if ( min(new.ni) == 0)
		{
			return ( NULL );
		}

		new.RSS = calculateRSS(n, new.ni, p, k, new.X, new.cent, dat, method);
		if (new.RSS < curr.RSS)
		{
			curr.cent = new.cent;
			curr.X = new.X;
			curr.RSS = new.RSS;
			curr.ni = new.ni;
			cont = TRUE;
			iter = iter+1;
		}
	}

	retObj$final.RSS = curr.RSS;
	retObj$final.cent = curr.cent;
	retObj$final.assign = curr.X;
	retObj$final.ni = curr.ni;
	retObj$final.iter = iter;

	return (retObj);
}



###############################################################
#
#	kmeansIteration
#
###############################################################
kmeansIteration <- function( simObj, k, method, startpntmet )
{
	start.time = Sys.time();

	retObj = simObj;
	n = simObj$n;
	p = simObj$p;
	dat = simObj$dat;

	# Set centers according to chosen method
	init.cent = switch (startpntmet, 
			"Random" = randomStartPoints( dat, n, k ),
			"Fraley" = fraleyStartPoints( dat, n, k )
		  );
	# in case k=1 than set at matrix instead of vector
	init.cent = matrix( init.cent, ncol=p );

	init.X = calculateAssignment(n, p, k, init.cent, dat);
	
	# Execute the algorithm
	kmeans.ret = findClusters( n, p, k, init.X, init.cent, dat, method );
	if ( is.null(kmeans.ret) )
	{
		retObj = NULL;
	}
	else
	{
		retObj$final.RSS = kmeans.ret$final.RSS;
		retObj$final.iter = kmeans.ret$final.iter;
		retObj$final.cent = kmeans.ret$final.cent;
		retObj$final.assign = kmeans.ret$final.assign;
		retObj$final.ni = kmeans.ret$final.ni;
		retObj$final.k = k;
		retObj$final.time = (Sys.time() - start.time);
	}

	return( retObj );
}



###############################################################
#
#	fraleyStartPoints:
#		dat:	data matrix (N*p)
#		n:	Number of observations (not in use)
#		k:	desired k
#
#	return value:
#		centres table (p*k)
#
###############################################################
fraleyStartPoints <- function( dat, n, k )
{
	hcTree = hc(modelName="VVV", data=dat);
	cl = hclass( hcTree, k );
#	hc = hclust( dist(dat^2), "cen" );
#	memb = cutree(hc, k);
	cent = NULL;
	for (i in 1:k)
	{
#		cent <- rbind(cent, colMeans(dat[memb == i, , drop = FALSE])) 
		cent <- rbind(cent, colMeans(dat[cl == i, , drop = FALSE])) 
	}
	return ( cent );
}



###############################################################
#
#	randomStartPoints:
#		dat:	data matrix (N*p)
#		n:	Number of observations
#		k:	desired k
#
#	return value:
#		centres table (p*k)
#
###############################################################
randomStartPoints <- function( dat, n, k)
{
	smp.cent = sample(n, k);
	cent = dat[smp.cent, ];
	return ( cent );
}


















###############################################################
#
#	kmeans.step:
#		n:			Number of observations
#		p:			Size of dimensions
#		k:			Number of clusters
#		curr.X:		Current assignment matrix
#		curr.cent:	Current centers matrix
#		dat:			Data matrix
#		method:		Choose "Regular" or "Weighted"
#
#	activity:
#		applies a single step of k-means.
#
#	return value:
#		Object containing:
#			RSS:		RSS value of new configuration
#			cent:	Centres of new configuration
#			X:		Assignment matrix of new configuration
#			ni:		Number of observations per cluster
#
###############################################################
kmeans.step <- function(n, p, k, curr.X, curr.cent, dat, method)
{
	retObj = NULL;

	curr.cent = matrix( curr.cent, ncol=p);		# Fix curr.cent matrix

	# Check if all clusters have members
	curr.ni = calculateNumberOfMembers(curr.X);
	if ( min(curr.ni) == 0)
	{
		return ( NULL );
	}

	# Calculate RSS and initialize the
	curr.RSS = calculateRSS(n, curr.ni, p, k, curr.X, curr.cent, dat, method);

#	new.cent = calculateCenters(curr.ni, curr.X, dat);
	new.X = calculateAssignment(n, p, k, curr.cent, dat);

	# Check that all clusters have members
	new.ni = calculateNumberOfMembers(new.X);
	if ( min(new.ni) == 0)
	{
		return ( NULL );
	}
	new.cent = calculateCenters(new.ni, new.X, dat);

	new.RSS = calculateRSS(n, new.ni, p, k, new.X, new.cent, dat, method);
	if (new.RSS < curr.RSS)
	{
		retObj$RSS = new.RSS;
		retObj$cent = new.cent;
		retObj$X = new.X;
		retObj$ni = new.ni;
	}

	return (retObj);
}
