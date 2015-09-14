
##############################################################
#
#	MAP.stat:
#		Set a score table for the correct K using the
#		the MAP method over a number of times
#
#	parameters:
#		simObj:		Simulated object
#		min.k:		minimum k
#		max.k:		maximum k
#		strt:   	"Random" or "Fraley" starting points
#		n:			number of repeats (optional=100)
#		approach:	"Binom" or "Unif"
#		q:			In case of binom approach, estimated q
#
##############################################################
MAP.stat <- function( simObj, min.k, max.k, strt="Random", n=100, approach, q=NA, a=NA )
{
	ans = rep(0, max.k-min.k+1);

	for (i in 1:n)
	{
		if ((i %% 10) == 0)
		{
			cat("MAP", approach,"q=",q,"a=",a,"at iteration", i, "/", n, "\n");
		}
		k = MAP.approach(simObj, min.k, max.k, strt, approach, q, a );
		ans[k-min.k+1] = ans[k-min.k+1]+1;
	}

	mat  = rbind( min.k:max.k, ans );
	return ( mat);
}



##############################################################
#
#	MAP.approach:
#		Calculates the correct K given a range
#		using the MAP rule (assumptions...)
#
#	parameters:
#		simObj:	Simulated object
#		min.k:	minimum k
#		max.k:	maximum k
#
##############################################################
MAP.approach <- function(simObj, min.k, max.k, strt, approach, q, a )
{
	maps = rep(Inf, max.k);
	p = simObj$p;
	n = simObj$n;
	dat = simObj$dat;
	kmeans.method = "MAP";
	if (approach == "Binder")	kmeans.method = "Regular";

	for (j in min.k:max.k)
	{
		
		ansObj = kmeansIteration( simObj, j, kmeans.method, strt );
		if (!is.null(ansObj))
		{
			Wk = MAP.Wk(dat, ansObj$final.cent, ansObj$final.assign);
			Pen.Nj = MAP.pen.nj(n, ansObj$final.ni );

			Pen.priors = switch (approach, 
				"Binom" = MAP.penalty.priors.binom(n, p, j, q),
				"No2" = MAP.penalty.priors.2(n, p, j, a),
				"Unif" = MAP.penalty.priors.unif(n, p, j, min.k, max.k),
				"Binder" = MAP.penalty.binder(n, p, j, a) - Pen.Nj  # cancel
			  );
		#	Pen.priors = MAP.penalty.priors(n, p, j);
			maps[j] = log(Wk) + Pen.Nj + Pen.priors;
		}
	}
	return ( which.min(maps) );
}


##############################################################
#
#	MAP.Wk:
#		Calculates the within-cluster distances
#
#	parameters:
#		dat:	data
#		cent:	matrix of centroids
#		X:	assignment matrix (which observation 
#			belongs to which cluster)
#
##############################################################
MAP.Wk <- function(dat, cent, X)
{
	cnt = X %*% cent;
	dst = dat-cnt;
	dst = dst^2;
	Wk = sum(dst);
	return(Wk);
}


##############################################################
#
#	MAP.pen.nj:
#		Calculates the MAP penalty of the sizes of 
#		the clusters
#
#	parameters:
#		n:	total number of observations
#		nj:	size of the j-th cluster 
#
##############################################################
MAP.pen.nj <- function(n, nj)
{
	lnNj = log(nj);
	return( sum(lnNj)/n );
}


##############################################################
#
#	MAP.penalty.priors:
#		Calculates the penalty of the MAP rule
#		for the binomial priors (assumptions....)
#
#	parameters:
#		n:	total number of observations
#		p:	dimension of the data
#		k:	number of the clusters
#
##############################################################
MAP.penalty.priors.binom <- function(n, p, k, q)
{
#	q = 0.02;
	term1 = 2/p*log(k);
	term2 = (1-q)/(q*k);
	term2 = (2*k)/(n*p)*log(term2);
	return (term1 + term2);
}



##############################################################
#
#	MAP.penalty.priors:
#		Calculates the penalty of the MAP rule
#		for the #2 priors (assumptions....)
#
#	parameters:
#		n:	total number of observations
#		p:	dimension of the data
#		k:	number of the clusters
#
##############################################################
MAP.penalty.priors.2 <- function(n, p, k, a)
{
	term1 = 2/p*log(k);
	term2 = 2*(1-1/p)*k/n*log(k);
	term3 = 2*k/(n*p)*log(a);
	return ( term1 + term2 + term3 );
}



##############################################################
#
#	MAP.penalty.priors.unif:
#		Calculates the penalty when the uniform
#		prior.
#
#	parameters:
#		n:	total number of observations
#		p:	dimension of the data
#		k:	number of the clusters
#		kmin:	minimum tested number of clusters
#		kmax: maximal tested number of clusters
#
##############################################################
MAP.penalty.priors.unif <- function(n, p, k, kmin, kmax)
{
	term1 = 2/p*log(k);
	term2 = 2/(n*p) * log(kmax-kmin);
	term3 = 2/(n*p) * log( factorial(n)/(factorial(k)*factorial(n-k)) );
	term4 = -2*k/(n*p) * log(k);
	return (term1 + term2 + term3 + term4 );
}




##############################################################
#
#	MAP.penalty.binder:
#		Calculates the penalty for binder's prior
#
#	parameters:
#		n:	total number of observations
#		p:	dimension of the data
#		k:	number of the clusters
#		a:	binder's parameter
#
##############################################################
MAP.penalty.binder <- function(n, p, k, a)
{
	j = 1:(k-1);
	res = (k-1)*log(a) - sum(log(j));
	return( -2/(n*p) * res );
}
