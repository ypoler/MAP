##############################################################
#
#	Gap.approach:
#		Calculate the correct K given a range
#		using the Gap method
#
#	parameters:
#		simObj:	Simulated object
#		min.k:	minimum k
#		max.k:	maximum k
#		B:		number of repeats
#
##############################################################
Gap.approach <- function( simObj, min.k, max.k, strt, B ) 
{
	rng = rbind(apply(simObj$dat, 2, min),
				apply(simObj$dat, 2, max) );

	simB = simObj;
	Wbj = rep(NA, B);
	GAP = rep(NA, max.k);
	Sj = rep(NA, max.k);
	for (j in min.k:max.k)
	{
		# calculate the ref. data
		for (i in 1:B)
		{
			sim.dat = createSimDat(rng, simObj$n, simObj$p);
			simB$dat = sim.dat;
			simB.ans = kmeansIteration(simB, j, "Regular", strt);
			if (!is.null(simB.ans))
			{
				Wbj[i] = Wk(sim.dat, simB.ans$final.assign, simObj$p, j); 
			}
		}
		
		# Calculates the Wk got after true clustering and then
		# the GAP and Sj for this J
		ansObj = kmeansIteration(simObj, j, "Regular", "Random");
		if (!is.null(ansObj))
		{
			Wj = Wk(simObj$dat, ansObj$final.assign, simObj$p, j);
			GAP[j] = mean(log(Wbj), na.rm=TRUE) - log(Wj);
			cnt = sum( !is.na(Wbj) );
			Sj[j] = sd(log(Wbj), na.rm=TRUE) * sqrt(1+1/cnt);
		}
	}

	# Find the 1st J that's bigger than previous (after outliers)
	for (j in min.k:(max.k-1))
	{
		if (GAP[j] >= GAP[j+1] - Sj[j+1])
		{
			return (j);
		}
	}
	return(j);
}


##############################################################
#
#	createSimDat:
#		create a refference dataset of size n from a
#		uniform distribution ranging between original
#		data range.
#
#	parameters:
#		rng:	matrix of minimal and maximal ranges
#			(per each dimension)
#		n:	number of observations
#		p:	number of dimension
#
##############################################################
createSimDat <- function(rng, n, p)
{
	for (j in 1:p)
	{
		vec.dat = runif(n, min=rng[1, j], max=rng[2, j]);
		if (j == 1)
		{
			sim.dat = vec.dat;
		}
		else
		{
			sim.dat = cbind(sim.dat, vec.dat);
		}
	}
	return (sim.dat);
}



##############################################################
#
#	Wk:
#		calculates the Wk measure for GAP
#
#	parameters:
#		dat:	complete data set
#		X:	assignment matrix
#		p:	number of dimensions
#		k:	number of clusters to be checked
#
##############################################################
Wk <- function(dat, X, p, k)
{
	W = 0;
	for (j in 1:k)
	{
		idx = (X[,j] == 1);
		Cj = matrix( dat[idx,], ncol=p);	# in case |Cj|=1 so it won't become vector
		nj = dim(Cj)[1];
		Dj = sum( dist(Cj, diag=FALSE, upper = FALSE) );	#use the distance matrix
		#Dj = DClusterj(Cj, nj);
		W = W + Dj/nj;
	}
	return ( W );
}



##############################################################
#
#	DClusterj:
#		Calculates the overall distance between each
#		two observations belonging to the same cluster.
#
#	parameters:
#		Cj:	a cluster - set of observations
#		nj:	size of the cluster
#
##############################################################
##DClusterj <- function(Cj, nj)
##{
##	Dj = 0;
##	for (a in 1:nj)
##	{
##		for (b in a:nj)
##		{
##			x.a = Cj[a,];
##			x.b = Cj[b,];
##			dst = sum( (x.a - x.b)^2 );
##			dst = sqrt(dst);
##			Dj = Dj + dst;
##		}
##	}
##	return (Dj);
##} 


##############################################################
#
#	Gap.stat:
#		Set a score table for the correct K using the
#		the Gap method over a number of times
#
#	parameters:
#		simObj:	Simulated object
#		min.k:	minimum k
#		max.k:	maximum k
#		strt:   "Random" or "Fraley" starting points
#		n:		number of repeats (optional=100)
#		B:		size of ref. data (optional=50)
#
##############################################################
Gap.stat <- function( simObj, min.k, max.k, strt="Random", n=100, B=50 )
{
	ans = rep(0, max.k-min.k+1);

	for (i in 1:n)
	{
		if ((i %% 10) == 0)
		{
			cat("Gap at iteration", i, "/", n, "\n");
		}
		k = Gap.approach(simObj, min.k, max.k, strt, B);
		ans[k-min.k+1] = ans[k-min.k+1]+1;
	}

	mat  = rbind( min.k:max.k, ans );
	return ( mat);
}


