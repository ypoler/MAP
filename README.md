1. About:
===========

This library contains the functions used in my Thesis:

"Bayesian Approach to Clustering 
Finding the number of clusters: the MAP rule

Submitted in partial fulfillemnt of graduation requirments for the degree of M.Sc. in Statistics

prepared under the supervision of rofessor Felix Abramovich 

January 2010

Tel Aviv University
The Raymond and Beverly Sackler Faculty of Exact Sciences
School of Mathematical Sciences
Department of Statistics and Operation Research

By Yinat Trompoler"


The functions are both for the MAP rule described in the Thesis and other approches compared to: 

KL: KRZANOWSKI, W. J., AND LAI, Y. T. A criterion for determining the number
    of groups in a data set using sum-of-squares clustering. Biometrics 44, 1 (1988), 23–34.

Validity: RAY, S., AND TURI, R. Determination of number of clusters in k-means clustering
          and application in colour image segmentation. In Proceedings of the 4th
          International Conference on Advances in Pattern Recognition and Digital Techniques (1999), pp. 137–143.


Jump: SUGAR, C. A., AND JAMES, G. M. Finding the number of clusters in a dataset:
      An information-theoretic approach. Journal of the American Statistical Association 93, 463 (2003), 750–763.

GAP: TIBSHIRANI, R., WALTHER, G., AND HASTIE, T. Estimating the number of
     clusters in a data set via the gap statistic. J. R. Statist. Soc. B 63, 2 (2001), 411–423.


The MAP functions are for Binom, Uniform priors. 




2. Required packages:
=======================

mclust (from c-ran repository)




3. Installation:
==================

The current version is of non-package format, therefore, all function files should be included manually to the R scripts.

1. Extract the files to a directory
2. Open R and set the working direcory to the root of the extracted folder (above directories "approaches" and "functions")
3. In your script which needs to use the criteria, include the following files:
source("functions\\sim.obj.R");
source("functions\\kmeans.functions.R");
source("approaches\\KL.R");
source("approaches\\Validity.R");
source("approaches\\Gap.R");
source("approaches\\Jump.R");
source("approaches\\MAP.R");




4. Data object:
==================

The object which represents the data is called "simObj". It isn't of strong-typed form, and contains the following fields:
n:		Number of observations
p:		Number of dimensions (input p)
real.k:		Number f clusters (input k)
real.cent: 	Matrix of centres (input cent)
dat:		Matrix of simulated data (n*p)
real.assign:	Matrix of indicator assignments
		(cell i,j=1 if observation i
		belongs to cluster j) - n*k

The clustering functions require only n, p and dat to be part of given simObj - the rest of the fields are created in case the simObj is simulated via createRandomClusters function.



5. API:
=========

Each approach, has an API of its own which takes a simObj, range number of clusters, number of repeats to run the k-means iteration:

KL:
#	KL.stat:
#		Set a score table for the correct K using the
#		the K-L method over a number of times
#
#	parameters:
#		simObj:	Simulated object
#		min.k:	minimum k
#		max.k:	maximum k
#		strt:   "Random" or "Fraley" starting points
#		n:	number of repeats (optional=100)

Validity:
#	Validity.stat:
#		Set a score table for the correct K using the
#		the K-L method over a number of times
#
#	parameters:
#		simObj:	Simulated object
#		min.k:	minimum k
#		max.k:	maximum k
#		strt:   "Random" or "Fraley" starting points
#		n:	number of repeats (optional=100)

Jump:
#	Jump.stat:
#		Set a score table for the correct K using the
#		the Jump method over a number of times
#
#	parameters:
#		simObj:	Simulated object
#		min.k:	minimum k
#		max.k:	maximum k
#		strt:   "Random" or "Fraley" starting points
#		n:	number of repeats (optional=100)

Gap:
#	Gap.stat:
#		Set a score table for the correct K using the
#		the Gap method over a number of times
#
#	parameters:
#		simObj:	Simulated object
#		min.k:	minimum k
#		max.k:	maximum k
#		strt:   "Random" or "Fraley" starting points
#		n:	number of repeats (optional=100)
#		B:	size of ref. data (optional=50)
NOTE: Gap approach is very time consuming, as its calculation differs from other approaches

MAP:
#	MAP.stat:
#		Set a score table for the correct K using the
#		the MAP method over a number of times
#
#	parameters:
#		simObj:		Simulated object
#		min.k:		minimum k
#		max.k:		maximum k
#		strt:   	"Random" or "Fraley" starting points
#		n:		number of repeats (optional=100)
#		approach:	"Binom" or "Unif"
#		q:		In case of binom approach, estimated q



6. Example:
==============

a. simulate 2D 3 clusters dataset:
------------------------------------

k = 3;
ni = rep(50, k);
n = sum(ni);
p = 10;
cent = matrix( runif(p*k, min=-10, max=10), ncol=p, nrow=k)
simObj = createRandomClusters( ni, p, k, cent, sdi=.5 );

# Write the simObject to a file
write.simObject(simObj, "datasets\\sim.3.clusters.2D.txt");

# Read the simObject from the file
simObj = read.simObject("datasets\\sim.3.clusters.2D.txt");

# Cluster number range
k.min = 1;
k.max = 6;

# Run KL approach on the data
res = KL.stat( simObj, k.min, k.max, strt="Random", n=50 );
plot(res[2,]~res[1,], xlab="num. clusters", ylab="freq.")	# plot the frq. of the hist

# Run Validity approach on the data
res = Validity.stat( simObj, k.min, k.max, strt="Random", n=50 );
plot(res[2,]~res[1,], xlab="num. clusters", ylab="freq.")	# plot the frq. of the hist

# Run Jump approach on the data
res = Jump.stat( simObj, k.min, k.max, strt="Random", n=50 );
plot(res[2,]~res[1,], xlab="num. clusters", ylab="freq.")	# plot the frq. of the hist

# Run MAP with binom prior (q=0.05) on the data
res = MAP.stat( simObj, k.min, k.max, strt="Random", n=50, approach="Binom", q=0.05 );
plot(res[2,]~res[1,], xlab="num. clusters", ylab="freq.")	# plot the frq. of the hist

# Run MAP with uniform prior on the data
res = MAP.stat( simObj, k.min, k.max, strt="Random", n=50, approach="Unif");
plot(res[2,]~res[1,], xlab="num. clusters", ylab="freq.")	# plot the frq. of the hist



b. Import Iris data:
----------------------

irisObj = NULL;
irisObj$p = 4;
irisObj$n = 150;
irisObj$dat = matrix(nrow=irisObj$n, ncol=irisObj$p);
for (i in 1:irisObj$p)
{
	irisObj$dat[,i] = iris[,1+i];
}

# Cluster number range
k.min = 1;
k.max = 7;

# Run KL approach on the data
res = KL.stat( irisObj, k.min, k.max, strt="Random", n=50 );
plot(res[2,]~res[1,], xlab="num. clusters", ylab="freq.")	# plot the frq. of the hist

# Run Validity approach on the data
res = Validity.stat( irisObj, k.min, k.max, strt="Random", n=50 );
plot(res[2,]~res[1,], xlab="num. clusters", ylab="freq.")	# plot the frq. of the hist

# Run Jump approach on the data
res = Jump.stat( irisObj, k.min, k.max, strt="Random", n=50 );
plot(res[2,]~res[1,], xlab="num. clusters", ylab="freq.")	# plot the frq. of the hist

# Run MAP with binom prior (q=0.05) on the data
res = MAP.stat( irisObj, k.min, k.max, strt="Random", n=50, approach="Binom", q=0.05 );
plot(res[2,]~res[1,], xlab="num. clusters", ylab="freq.")	# plot the frq. of the hist

# Run MAP with uniform prior on the data
res = MAP.stat( irisObj, k.min, k.max, strt="Random", n=50, approach="Unif");
plot(res[2,]~res[1,], xlab="num. clusters", ylab="freq.")	# plot the frq. of the hist

NOTE: All approaches perform poorly on the Iris example and tend to assign the largest number of clusters (7) instead of the correct number
