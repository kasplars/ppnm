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

int main(){
	// Volume of sphere of radius 0.5

	double a1[]={-0.5, -0.5, -0.5}, b1[]={0.5,0.5,0.5};
	double f(int dim, double * x) {
		double val = 1.0, sum = 0.0, norm; for (int i = 0;i<dim;i++) sum+=x[i]*x[i]; norm = sqrt(sum); if (norm>0.5) val*=0.0; 
		return val;
	}
	double res, err; int dim = 3, N = 1e6;
	plainmc(dim,f,a1,b1,N,&res,&err);
	printf("Volume of sphere of radius 0.5 using pseudo-Monte Carlo.\nReal volume = %.8g, calculated volume = %.8g, calculated error = %.8g. Number of points = %i\n\n\n",4.0*M_PI*pow(0.5,3)/3.0,res,err,N);
	
	// Integral given in homework 
	
	double a2[]={0.0, 0.0, 0.0}, b2[]={M_PI,M_PI,M_PI};
	dim = 3, N = 1e7;
	double g(int dim, double * x) {double aux = 1-cos(x[0])*cos(x[1])*cos(x[2]); return pow(aux,-1.0)/pow(M_PI,3);}
	plainmc(dim,g,a2,b2,N,&res,&err);
	printf("Integral given in homework.\nCalculated volume = %.8g, calculated error = %.8g. Number of points = %i\n\n\n",res,err,N);
	
	// Qausi monte carlo, volume of sphere of radius 0.5
	
	dim = 3;
	int Nmax = 1e5, numsteps = 100; double resquasi, errquasi;
	double h = 1.0 * Nmax / numsteps, step;
	double scaling1 = 0.5, scaling2 = 3e-7;
	FILE * data = fopen("comparison.txt","w");
	for (int i = 1; i<=numsteps; i++) {
		step = h*i;
		plainmc(dim,f,a1,b1,step,&res,&err);
		haltonmc(dim,f,a1,b1,step,&resquasi,&errquasi);
		fprintf(data,"%10g %10g %10g %10g %10g\n",step,res,err*sqrt(step)/scaling1,resquasi,errquasi*pow(log(N),dim)/N/scaling2);
	}
	fclose(data);
	
return 0;
}








