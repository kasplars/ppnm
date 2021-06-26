#include <stdio.h>
#include <assert.h>
#include <math.h>
#include <gsl/gsl_integration.h>

double f(double x) {return 1.0/sqrt(x);}
double p(double x) {return 4.0*sqrt(1-pow(x,2));}

double integrate(double (*f)(double), double a, double b, double f2, double f3, double delta, 
	double epsilon, int numrecs)
{
	// assert(numrecs<1E5);
	TRACE("1.2\n");
	double f1 = f(a+(b-a)/6), f4 = f(a+5.0*(b-a)/6); 
	double Q = (b-a)/6*(2.0*f1+f2+f3+2.0*f4), q = (b-a)/4*(f1+f2+f3+f4), err=fabs(Q-q);
	if (err < delta+epsilon*fabs(Q)) return Q;
	else return integrate(f,a,(a+b)/2,f1,f2,delta/sqrt(2),epsilon,numrecs+1)+
		    integrate(f,(a+b)/2,b,f3,f4,delta/sqrt(2),epsilon,numrecs+1);
}

double recadapter(double (*f)(double), double a, double b, double delta, double epsilon) {
	double f2 = f(a+2.0*(b-a)/6), f3 = f(a+4.0*(b-a)/6); int numrecs = 0;
	TRACE("1.1\n");
	return integrate(f,a,b,f2,f3,delta,epsilon,numrecs);
}

double clenshaw_curtis(double (*f)(double),double a,double b,double delta,double epsilon){
	// variable transformation x -> (a+b)/2+(a-b)/2*cos(t) makes \int_a^b f(x)dx on the right form.
	double g(double t){return f((a+b)/2+(a-b)/2*cos(t))*sin(t)*(b-a)/2;}
	return recadapter(g,0,M_PI,2*delta,2*epsilon);
}

int main() {
	double a = 0.0, b = 1.0, delta = 0.01, epsilon = 0.01; int numevals = 0;
	double function(double x){numevals++; return f(x);}
	double integral = recadapter(function,a,b,delta,epsilon);
	
	printf("Integral of f(x) from %g to %g is %g not using Clenshaw-Curtis. Number of evaluations: %i.\n",a,b,integral,numevals);
	numevals = 0;
	
	integral = clenshaw_curtis(function,a,b,delta,epsilon);
	
	printf("Integral of f(x) from %g to %g is %g using Clenshaw-Curtis. Number of evaluations: %i.\n",a,b,integral,numevals);
	numevals = 0;
	
	printf("So there is a lot of improvement in numbers of evaluations.\n");
	
	a = 0.0, b = 1.0, delta = 0.000000001, epsilon = 0.000000001;
	double function_new(double x){numevals++; return p(x);}
	double integral1 = clenshaw_curtis(function_new,a,b,delta,epsilon); int numevals_aux = numevals; numevals = 0;
	double integral2 = recadapter(function_new,a,b,delta,epsilon);
	printf("The integral of p(x) from %g to %g without and with Clenshaw-Curtis, respectively:\n",a,b);
	printf("%.15g; number of evaluations: %i; error: %.15g\n",integral1,numevals_aux,fabs(M_PI-integral1));
	printf("%.15g; number of evaluations: %i; error: %.15g\n",integral2,numevals,fabs(M_PI-integral2));
	printf("We see that Clenshaw-Curtis takes more evaluations but is more accurate.\n");
	
	int limits = 1e6; double result, error; numevals = 0;
	gsl_integration_workspace * workspace = gsl_integration_workspace_alloc(limits);
	gsl_function F; 
	double function_newest(double x,void * params){numevals++; return p(x);}
	F.function = &function_newest;
	
	gsl_integration_qags(&F,a,b,delta,epsilon,limits,workspace,&result,&error);
	
	printf("QAGS from GSL library, the result, error and number of function evaluations are: \n%.15g\n%.15g\n%i,\nrespectively.\n",result,error,numevals);
	
	gsl_integration_workspace_free (workspace);
	
	
return 0;
}








