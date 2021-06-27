#include <stdio.h>
#include <assert.h>
#include <math.h>
#include <gsl/gsl_integration.h>
#include "functions.h"

double f(double x) {return 1.0/sqrt(x);}
double p(double x) {return 4.0*sqrt(1-pow(x,2));}
double gauss(double x) {return exp(-x*x);}
double piintegral(double x) {return pow((x*x+1.0),-1);}

int main() {
	// Part A and B
	
	printf("Part A and B\n\n");

	double a = 0.0, b = 1.0, delta = 0.01, epsilon = 0.01, errorestimate; int numevals = 0;
	double function(double x){numevals++; return f(x);}
	double integral = recadapter(function,a,b,delta,epsilon,&errorestimate);
	
	printf("Integral of 1/sqrt(x) from %g to %g is %.10g not using Clenshaw-Curtis. Number of evaluations: %i.\n",a,b,integral,numevals);
	numevals = 0;
	
	integral = clenshaw_curtis(function,a,b,delta,epsilon,&errorestimate);
	
	printf("Integral of 1/sqrt(x) from %g to %g is %.10g using Clenshaw-Curtis. Number of evaluations: %i.\n\n",a,b,integral,numevals);
	numevals = 0;
	
	printf("The real value is 2. Relative and absolute error is set to 0.01. \nTo conclude, there is a lot of improvement in numbers of evaluations.\n\n");
	
	a = 0.0, b = 1.0, delta = 0.000000001, epsilon = 0.000000001;
	double function_new(double x){numevals++; return p(x);}
	double integral1 = clenshaw_curtis(function_new,a,b,delta,epsilon,&errorestimate); int numevals_aux = numevals; numevals = 0;
	double integral2 = recadapter(function_new,a,b,delta,epsilon,&errorestimate);
	printf("The integral of p(x) (see code) from %g to %g with and without Clenshaw-Curtis, respectively:\n\n",a,b);
	printf("%.20g; number of evaluations: %i; error: %.20g\n",integral1,numevals_aux,fabs(M_PI-integral1));
	printf("%.20g; number of evaluations: %i; error: %.20g\n\n",integral2,numevals,fabs(M_PI-integral2));
	printf("We see that Clenshaw-Curtis takes more evaluations but is more accurate.\n\n");
	
	int limits = 1e6; double result, error; numevals = 0;
	gsl_integration_workspace * workspace = gsl_integration_workspace_alloc(limits);
	gsl_function F; 
	double function_newest(double x,void * params){numevals++; return p(x);}
	F.function = &function_newest;
	
	gsl_integration_qags(&F,a,b,delta,epsilon,limits,workspace,&result,&error);
	
	printf("Using QAGS from GSL library, the result, error and number of function evaluations are found to be: \n%.15g, %.15g and %i, respectively.\n",result,error,numevals);
	
	
	
	// Part C 
	
	printf("\n\nPart C\n\n");
	
	double functionC(double x) {numevals++; return gauss(x);}
	double piinteg(double x) {numevals++; return piintegral(x);}
	double functionC2(double x,void * params) {numevals++; return gauss(x);}
	double functionC3(double x, void * params) {numevals++; return piinteg(x);}
	
	numevals = 0; delta = 0.000001; epsilon = 0.000001;
	
	printf("Using own routine:\n\n");
	
	double val = recadapter(functionC,-INFINITY, INFINITY, delta,epsilon,&errorestimate);
	printf("The gaussian function integrated from -infinity to infity is equal to %.10g with error estimate %.10g. The real value is sqrt(pi). The number of evaluations were %i.\n\n",val,errorestimate,numevals);
	
	numevals = 0; 
	
	val = recadapter(piinteg,0.0,INFINITY,delta,epsilon,&errorestimate);
	printf("The function 1/(x*x+1) integrated from 0 to infity is equal to %.10g with error estimate %.10g. The real value is pi/2. The number of evaluations were %i.\n\n",val,errorestimate,numevals);
	
	
	
	printf("Using GSL routine:\n\n");
	
	numevals = 0;
	
	gsl_integration_workspace * workspace2 = gsl_integration_workspace_alloc(limits);
	gsl_function G; 
	G.function = &functionC2;
	gsl_integration_qagi(&G,delta,epsilon,limits,workspace2,&result,&error);
	printf("The gaussian function integrated from -infinity to infity is equal to %.10g with error estimate %.10g. The real value is sqrt(pi). The number of evaluations were %i.\n\n",result,error,numevals);
	
	numevals = 0;
	
	gsl_integration_workspace * workspace3 = gsl_integration_workspace_alloc(limits);
	gsl_function H;
	H.function = &functionC3;
	gsl_integration_qagiu(&H,0,delta,epsilon,limits,workspace3,&result,&error);
	printf("The function 1/(x*x+1) integrated from 0 to infity is equal to %.10g with error estimate %.10g. The real value is pi/2. The number of evaluations were %i.\n\n",result,error,numevals);
	
	printf("All in all, GSL performs much better.\n");
	
	gsl_integration_workspace_free (workspace);
	gsl_integration_workspace_free (workspace2);
	gsl_integration_workspace_free (workspace3);
return 0;
}








