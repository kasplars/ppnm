#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_linalg.h>
#include <float.h>
#include "library.h"

#ifndef NDEBUG
#define TRACE(args...) fprintf(stderr,args)
#else
#define TRACE(...)
#endif

static double epsilon;

void f(gsl_vector * x, gsl_vector * fx) {
	double scale = pow(0.001,-3);
	double t, y;
	for (int i = 0; i < x->size; i++) {
		t = gsl_vector_get(x,i);
		y = t*t*t;
		gsl_vector_set(fx,i,y);
	}
	gsl_vector_scale(fx,scale);	// remember to rescale problem!
}

void rosenbrockgrad(gsl_vector * t, gsl_vector * tx) {
	double scale = 1000;
	double x, y, gradx, grady;
	x = gsl_vector_get(t,0);
	y = gsl_vector_get(t,1);
	gradx = 400*x*x*x-400*x*y+2*x-2;
	grady = 200*(y-x*x);
	gsl_vector_set(tx,0,gradx); gsl_vector_set(tx,1,grady);
	gsl_vector_scale(tx,scale);	// remember to rescale problem!
}

void schroedinger(double r,gsl_vector * y, gsl_vector * dydt) {
	gsl_vector_set(dydt,0,gsl_vector_get(y,1));
	gsl_vector_set(dydt,1,-2.0*gsl_vector_get(y,0)*(epsilon+1.0/r));
}

void mfunc(gsl_vector * x, gsl_vector * tx) {
	double rinit = 0.001, rmax = 8.0;
	int n = 2;
	double h = 0.001, acc = 0.001, relacc = 0.001;
	
	
	gsl_vector * yinit = gsl_vector_alloc(n);
	gsl_vector * yend = gsl_vector_alloc(n);
	gsl_vector_set(yinit,0,rinit - rinit * rinit);
	gsl_vector_set(yinit,1,1.0 - 2*rinit);
	epsilon = gsl_vector_get(x,0);
	
	FILE * odeoutput = fopen("odeoutput.txt","w");
	fprintf(odeoutput,"#numerical solution\n");
	
	driver(odeoutput,schroedinger,n,rinit,yinit,rmax,yend,h,acc,relacc);
	gsl_vector_set(tx,0,gsl_vector_get(yend,0));
	
	fprintf(odeoutput,"\n\n#exact solution\n");
	double rval, ryval, numsteps = 100;
	for (int i = 0; i < numsteps; i++) {
		rval = rinit + (rmax - rinit) * i / (numsteps - 1);
		ryval = rval * exp(-rval);
		fprintf(odeoutput,"%10g %10g %10g\n",rval,ryval,0.0);
	}
	
	fclose(odeoutput);
	gsl_vector_free(yinit);
	gsl_vector_free(yend);
}

void newton(void (*f)(gsl_vector * x, gsl_vector * fx), gsl_vector * x, double eps) {
	int n = x->size; double dx = sqrt(DBL_EPSILON);
	gsl_vector * deltax = gsl_vector_alloc(n);
	gsl_vector * y = gsl_vector_alloc(n);
	gsl_vector * fx = gsl_vector_alloc(n);
	gsl_vector * fy = gsl_vector_alloc(n);
	gsl_vector * df = gsl_vector_alloc(n);
	gsl_matrix * J = gsl_matrix_alloc(n,n);
	gsl_matrix * R = gsl_matrix_alloc(n,n);
	while (1) {
		f(x,fx);
		for (int i = 0; i < n; i++) {
			gsl_vector_set(x,i,gsl_vector_get(x,i)+dx);
			f(x,df);
			gsl_vector_sub(df,fx);
			for (int j = 0; j < n; j++) gsl_matrix_set(J,j,i,1.0*gsl_vector_get(df,j)/dx);
			gsl_vector_set(x,i,gsl_vector_get(x,i)-dx);
		}
		gsl_vector_scale(fx,-1.0);
		GS_decomp(J,R);
		GS_solve(J,R,fx,deltax);
		double lambda = 2.0;
		while (1) {
			lambda /= 2.0;
			gsl_vector_memcpy(y,x);
			gsl_blas_daxpy(lambda,deltax,y);
			f(y,fy);
			if ((gsl_blas_dnrm2(fy) < (1.0 - lambda / 2) * gsl_blas_dnrm2(fx)) || (lambda < 1.0/64)) break;
		}
		gsl_vector_memcpy(x,y); gsl_vector_memcpy(fx,fy);
		if (gsl_blas_dnrm2(fx) < eps || gsl_blas_dnrm2(deltax) < dx) break;
	}
	gsl_vector_free(deltax);
	gsl_vector_free(y);
	gsl_vector_free(fx);
	gsl_vector_free(fy);
	gsl_vector_free(df);
	gsl_matrix_free(J);
	gsl_matrix_free(R);
}

int main() {
	// Simply finding root of x^3
	double eps = 0.0001;
	gsl_vector * x = gsl_vector_alloc(1);
	double xguess = 10;
	gsl_vector_set(x,0,xguess);
	newton(f,x,eps);
	printf("Root of x^3 = %10g within tolerance %10g\n\n",gsl_vector_get(x,0),eps);
	
	// Rosenbrock function
	gsl_vector * rosenbrockvector = gsl_vector_alloc(2);
	eps = 0.001;
	double rosenbrockx = 1.4, rosenbrocky = 0.9; gsl_vector_set(rosenbrockvector,0,rosenbrockx); gsl_vector_set(rosenbrockvector,1,rosenbrocky);
	newton(rosenbrockgrad,rosenbrockvector,eps);
	vector_print("Roots of the gradient of Rosenbrock function = ",rosenbrockvector);
	printf("within tolerance %10g\n\n",eps);
	gsl_vector_free(x);
	gsl_vector_free(rosenbrockvector);
	
	// Bound states of hydrogen atom with shooting method for b.v. problems
	gsl_vector * epsvector = gsl_vector_alloc(1);
	double accuracy = 0.001;
	double epsinit = -2.0; // must be negative guess
	gsl_vector_set(epsvector,0,epsinit);
	newton(mfunc,epsvector,accuracy);
	vector_print("Root of the M(epsilon) = ",epsvector);
	printf("within tolerance %10g. This agrees well with the expected energy -1/2.\n\n",accuracy);
	
	gsl_vector_free(epsvector);
	
		
return 0;
}










