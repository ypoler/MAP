#################################################
#												#
#	get.dist.min:								#
#		Get the pair (x,y) with minimum 		#
#		distance between them, where x,y		#
#		are from two different groups: 			#
#		Used for single-link.					#
#												#
#												#
#	input:										#
#		obj:	simObj							#
#												#
#	return:										#
#		pair of (x,y): 							#
#												#
#################################################
get.dist.min <- function( obj )
{
	pair = NULL;
	# translate assignment matrix to vector of group indexes
	grp = apply(obj$real.assign %*% (1:obj$real.k), 1, sum);

	# distance matrix
	dists = as.matrix(dist(obj$dat));
	n = obj$n;
	last.min = Inf;

	# Go over matrix and find minimum pair
	for (i in 1:n)
	{
		for (j in i:n)
		{
			if (grp[i] == grp[j])
			{
				dists[i, j] = Inf;
				dists[j, i] = Inf;
			}
			else
			{
				if (dists[i, j] < last.min)
				{
					last.min = dists[i, j];
					pair$x = i;
					pair$y = j;
				}
			}
		}
	}

	return (pair);
}




#################################################
#												#
#	get.dist.max:								#
#		Get the pair (x,y) with maximum 		#
#		distance between them, where x,y		#
#		are from two different groups: 			#
#		Used for single-link.					#
#												#
#												#
#	input:										#
#		obj:	simObj							#
#												#
#	return:										#
#		pair of (x,y): 							#
#												#
#################################################
get.dist.max <- function( obj )
{
	pair = NULL;
	# translate assignment matrix to vector of group indexes
	grp = apply(obj$real.assign %*% (1:obj$real.k), 1, sum);

	# distance matrix
	dists = as.matrix(dist(obj$dat));
	n = obj$n;
	last.max = -Inf;

	# Go over matrix and find minimum pair
	for (i in 1:n)
	{
		for (j in i:n)
		{
			if (grp[i] == grp[j])
			{
				dists[i, j] = -Inf;
				dists[j, i] = -Inf;
			}
			else
			{
				if (dists[i, j] > last.max)
				{
					last.max = dists[i, j];
					pair$x = i;
					pair$y = j;
				}
			}
		}
	}

	return (pair);
}




#################################################
#												#
#	get.dist.mean:								#
#		Get the pair (x,y) with mean	 		#
#		distance between them, where x,y		#
#		are from two different groups: 			#
#		Used for single-link.					#
#												#
#												#
#	input:										#
#		obj:	simObj							#
#												#
#	return:										#
#		pair of (x,y): 							#
#												#
#################################################
get.dist.mean <- function( obj )
{
	pair = NULL;
	# translate assignment matrix to vector of group indexes
	grp = apply(obj$real.assign %*% (1:obj$real.k), 1, sum);

	pair = NULL;
	for (j in 1:real.k)
	{
		idx = (grp == j);
		mat = obj$dat[idx, ];
		mat.mean = apply(mat, 2, mean);
		pair = rbind(pair, mat.mean);
	}

	return (pair);
}
