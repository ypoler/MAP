##############################################################
#
#	KL.approach:
#		Calculate the correct K given a range
#		using the K-L method
#
#	parameters:
#		simObj:	Simulated object
#		min.k:	minimum k
#		max.k:	maximum k
#
##############################################################
KL.approach <- function( simObj, min.k, max.k, strt ) 
{
	DIFF = rep(0, max.k);
	j = min.k;
	p = simObj$p;
	ansObj = kmeansIteration( simObj, j, "Regular", strt);
	Sj = ansObj$final.RSS;
	for (j in (min.k+1):max.k)
	{
		Sj1 = Sj;
		ansObj = kmeansIteration( simObj, j, "Regular", "Random");
		if (is.null(ansObj))
		{
			Sj = Inf;
		}
		else
		{
			Sj = ansObj$final.RSS;
		}
		DIFF[j] = (j-1)^(p/2)*Sj1 - j^(p/2)*Sj;
	}
	KL = rep(0, max.k-1);
	for (j in 1:(max.k-1))
	{
		KL[j] = DIFF[j]/DIFF[j+1];
	}

	return ( which.max(abs(KL)) );
}



##############################################################
#
#	KL.stat:
#		Set a score table for the correct K using the
#		the K-L method over a number of times
#
#	parameters:
#		simObj:	Simulated object
#		min.k:	minimum k
#		max.k:	maximum k
#		n:		number of repeats (optional=100)
#
##############################################################
KL.stat <- function( simObj, min.k, max.k, strt="Random", n=100 )
{
	ans = rep(0, max.k-min.k+1);

	for (i in 1:n)
	{
		if ((i %% 10) == 0)
		{
			cat("KL at iteration", i, "/", n, "\n");
		}
		k = KL.approach(simObj, min.k, max.k, strt);
		ans[k-min.k+1] = ans[k-min.k+1]+1;
	}

	mat  = rbind( min.k:max.k, ans );
	return ( mat);
}


