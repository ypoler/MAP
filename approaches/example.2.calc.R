func <- function(n, k, a, p) 
{
	res1 = factorial(n)/(factorial(k)*factorial(n-k));
	res2 = (a/k)^(k*p);
	return( res1 * res2);
}


func.binder <- function(n, k, a)
{
	res1 = a^(k-1);
	res2 = factorial(k-1);
	return( res1/res2 );
}


n=250;
k = 1:n;
a = 15;
p = 2;

#res = func(n, k, a, p);
res = func.binder(n, k, a);
plot(res~k)
which.max(res);
