#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#define RND ((double)rand()/RAND_MAX)
#include "functions.h"


int main(){
	// Part A
	printf("Part A\n\n");

	// Volume of sphere of radius 0.5

	double a1[]={-0.5, -0.5, -0.5}, b1[]={0.5,0.5,0.5};
	double f(int dim, double * x) {
		double val = 1.0, sum = 0.0, norm; for (int i = 0;i<dim;i++) sum+=x[i]*x[i]; norm = sqrt(sum); if (norm>0.5) val*=0.0; 
		return val;
	}
	double res, err; int dim = 3, N = 1e6;
	plainmc(dim,f,a1,b1,N,&res,&err);
	printf("Volume of sphere of radius 0.5 using pseudo-Monte Carlo:\nReal volume = %.8g, calculated volume = %.8g, calculated error = %.8g. Number of points = %i\n\n\n",4.0*M_PI*pow(0.5,3)/3.0,res,err,N);
	
	// Integral given in homework 
	
	double a2[]={0.0, 0.0, 0.0}, b2[]={M_PI,M_PI,M_PI};
	dim = 3, N = 1e7;
	double g(int dim, double * x) {double aux = 1-cos(x[0])*cos(x[1])*cos(x[2]); return pow(aux,-1.0)/pow(M_PI,3);}
	plainmc(dim,g,a2,b2,N,&res,&err);
	printf("Integral given in homework calculated using pseudo-Monte Carlo:\nCalculated volume = %.8g, calculated error = %.8g. Number of points = %i\n\n\n",res,err,N);
	
	// Part B
	printf("Part B\n\n");
	// Qausi monte carlo, volume of sphere of radius 0.5
	
	dim = 3;
	int Nmax = 1e5, numsteps = 80; double resquasi, errquasi;
	double h = 1.0 * Nmax / numsteps, step;
	FILE * data = fopen("comparison.txt","w");
	for (int i = 1; i<=numsteps; i++) {
		step = h*i;
		plainmc(dim,f,a1,b1,step,&res,&err);
		haltonmc(dim,f,a1,b1,step,&resquasi,&errquasi);
		fprintf(data,"%5.10g %5.10g %5.10g %5.10g %5.10g\n",step,res,err,resquasi,errquasi);
	}
	fclose(data);
	
	haltonmc(dim,f,a1,b1,5000,&resquasi,&errquasi);
	
	printf("Volume of sphere of radius 0.5 using quasi-random Monte Carlo:\nCalculated volume = %.8g, calculated error = %.8g. Number of points = %i\n\n",resquasi,errquasi,5000);
	
	printf("The scalings of the errors between the pseudo- and quasi-random Monte Carlo are seen in the figure. We see that the error drops significantly\n");
	printf("faster with the number of points when using the quasi-random sequence than the pseudo-random one.");
	
return 0;
}








