##############################################################
#
#	Validity.approach:
#		Calculate the correct K given a range
#		using the validity (Compactness-seperation) method
#
#	parameters:
#		simObj:	Simulated object
#		min.k:	minimum k
#		max.k:	maximum k
#
##############################################################
Validity.approach <- function( simObj, min.k, max.k, strt ) 
{
	intras = rep(Inf, max.k);
	inters = rep(-Inf, max.k);
	p = simObj$p;

	if (min.k == 1)
	{
		min.k = 2;
	}

	for (j in min.k:max.k)
	{
		ansObj = kmeansIteration( simObj, j, "Regular", strt );
		if (!is.null(ansObj))
		{
			intras[j] = Dintra(simObj$dat, ansObj$final.assign, ansObj$final.cent, simObj$n);
			inters[j] = Dinter(ansObj$final.cent, j);
		}
	}
	validity = intras/inters;
	return ( which.min(validity) );
}



##############################################################
#
#	Dintra:
#		Calculate the within-clusters distances
#
#	parameters:
#		dat:	data
#		X:	assignment
#		cent:	centers
#		n:	number of observations
#
##############################################################
Dintra <- function(dat, X, cent, n)
{
	cnt = X %*% cent;
	dst = dat-cnt;
	dst = dst^2;
	dst = apply(dst, 1, sum);
	dst = sqrt(dst);
	rss = sum( dst );
	return ( rss/n );
}



##############################################################
#
#	Dinter:
#		Calculate the between-cluster minimal distance
#
#	patameters:
#		cent:	centers
#		k:	number of clustres
#
##############################################################
Dinter <- function( cent, k )
{
	dst = dist(cent);
	return( min(dst) );
}


##############################################################
#
#	Validity.stat:
#		Set a score table for the correct K using the
#		the K-L method over a number of times
#
#	parameters:
#		simObj:	Simulated object
#		min.k:	minimum k
#		max.k:	maximum k
#		strt:   "Random" or "Fraley" starting points
#		n:		number of repeats (optional=100)
#
##############################################################
Validity.stat <- function( simObj, min.k, max.k, strt="Random", n=100 )
{
	ans = rep(0, max.k-min.k+1);

	for (i in 1:n)
	{
		if ((i %% 10) == 0)
		{
			cat("Validity at iteration", i, "/", n, "\n");
		}
		k = Validity.approach(simObj, min.k, max.k, strt);
		ans[k-min.k+1] = ans[k-min.k+1]+1;
	}

	mat  = rbind( min.k:max.k, ans );
	return ( mat);
}


