#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#define RND ((double)rand()/RAND_MAX)

void randomx(int dim, double * a, double * b, double * x) {
	for(int i=0;i<dim;i++)x[i]=a[i]+RND*(b[i]-a[i]);
}

double vandercorput(int n, int base) {
	double q = 0, bk=(double)1/base;
	while(n>0){q += (n % base)*bk; n /= base; bk /= base;}
	return q;
}

void halton(int n, int dim, double * a, double * b, double * x, int errorestimate) {
	int base[] = {2,3,5,7,11,13,17,19,23,29,31,37,41,43,47,53,59,61,67,71,73,79,83,89,97};
	if (errorestimate == 0) {
		for(int i = 0; i<dim; i++) x[i]=a[i] + vandercorput(n,base[i]) * (b[i]-a[i]);
	} else {
		for(int i = 0; i<dim; i++) x[i]=a[i] + vandercorput(n,base[dim + i]) * (b[i]-a[i]);
	}
}

void plainmc(int dim, double (*f)(int dim, double * x), double * a, double * b, int N, double * res, double * err) {
	double V = 1.0; for(int i=0;i<dim;i++)V*=b[i]-a[i]; // auxilliary box volume
	double sum = 0.0, sum2 = 0.0, fx, x[dim];
	for(int i=0;i<N;i++){randomx(dim,a,b,x); fx=f(dim,x); sum+=fx; sum2+=pow(fx,2);}
	double favg = sum/N, f2avg = sum2 / N, var = f2avg - pow(favg,2);
	*res = favg*V; *err = sqrt(var/N)*V;
}

void haltonmc(int dim, double (*f)(int dim, double * x), double * a, double * b, int N, double * res, double * err) {
	double V = 1.0; for(int i=0;i<dim;i++)V*=b[i]-a[i]; // auxilliary box volume
	double sum = 0.0, sumaux = 0.0, fx, fxaux, x[dim], xaux[dim];
	for(int i=0;i<N;i++){halton(i,dim,a,b,x,0); fx=f(dim,x); sum+=fx;}
	for(int i=0;i<N;i++){halton(i,dim,a,b,xaux,1); fxaux=f(dim,xaux); sumaux+=fxaux;}
	double favg = sum/N, favgaux = sumaux/N;
	*res = favg*V; *err = fabs(favg-favgaux); // error estimate
}
