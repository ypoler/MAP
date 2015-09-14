##############################################################
#
#	Jump.approach:
#		Calculate the correct K given a range
#		using the Jump method
#
#	parameters:
#		simObj:	Simulated object
#		min.k:	minimum k
#		max.k:	maximum k
#
##############################################################
Jump.approach <- function( simObj, min.k, max.k, strt ) 
{
	cov.x = cov(simObj$dat);
	p = simObj$p;
	Y = p/2;
	Dj = rep(NA, max.k);

	for (j in min.k:max.k)
	{
		ansObj = kmeansIteration( simObj, j, "Regular", strt);
		if (!is.null(ansObj))
		{
			Dj[j] = Jump.Wk(simObj$dat, ansObj$final.assign, ansObj$final.cent, simObj$n, j)/p;
		}
	}

	J = diff( Dj^(-Y) );
	return( which.max(J)+1 );
}



##############################################################
#
#	Jump.Wk:
#		Calculate Wk (eqvivalent to dk in the
#		jump article).
#
#	parameters:
#		dat:	dataset
#		X:	assignment matrix
#		cent:	centers
#		n:	number of observations
#		k:	number of clusters
#
##############################################################
Jump.Wk <- function(dat, X, cent, n, k)
{
	W = 0;
	for (i in 1:n)
	{
		for (j in 1:k)
		{
			xi = dat[i,];
			mj = cent[j,];
			zij = X[i,j];
			W = W + zij * t(xi-mj) %*% (xi-mj);
		}
	}
	return(W);
}



##############################################################
#
#	Jump.stat:
#		Set a score table for the correct K using the
#		the Jump method over a number of times
#
#	parameters:
#		simObj:	Simulated object
#		min.k:	minimum k
#		max.k:	maximum k
#		strt:   "Random" or "Fraley" starting points
#		n:		number of repeats (optional=100)
#
##############################################################
Jump.stat <- function( simObj, min.k, max.k, strt="Random", n=100 )
{
	ans = rep(0, max.k-min.k+1);

	for (i in 1:n)
	{
		if ((i %% 10) == 0)
		{
			cat("Jump at iteration", i, "/", n, "\n");
		}
		k = Jump.approach(simObj, min.k, max.k, strt);
		ans[k-min.k+1] = ans[k-min.k+1]+1;
	}

	mat  = rbind( min.k:max.k, ans );
	return ( mat);
}


