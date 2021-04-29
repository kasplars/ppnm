#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_linalg.h>
#include <float.h>
#include "library.h"

double DELTA=sqrt(DBL_EPSILON);

#ifndef NDEBUG
#define TRACE(args...) fprintf(stderr,args)
#else
#define TRACE(...)
#endif

void sr1(gsl_matrix * B, gsl_vector * u, gsl_vector * y) {
	int n = B->size1;
	gsl_matrix * dB = gsl_matrix_alloc(n,n);
	double alpha; gsl_blas_ddot(u,y,&alpha);
	gsl_blas_dger(1.0/alpha,u,u,dB);
	gsl_matrix_add(B,dB);
	gsl_matrix_free(dB);
}

void armijo(void (*f)(gsl_vector * x, double * fx),gsl_vector * x, gsl_vector * gradfx, gsl_vector * dx) {
	double lambda = 1.0, inner; int num = 0, max = 1e3;
	double alpha = 1e-6;
	double fx, fxs;
	gsl_vector * xs = gsl_vector_alloc(x->size);
	while (num < max) {
		num++;
		gsl_vector_memcpy(xs,x);
		gsl_vector_add(xs,dx);
		f(xs,&fxs);
		f(x,&fx);
		gsl_blas_ddot(dx,gradfx,&inner);
		if (lambda < DELTA) break;
		if (fxs < fx + alpha*inner) break;
		lambda /= 2.0;
		gsl_vector_scale(dx,lambda);
	}
	gsl_vector_free(xs);
}

void gradient(void (*f)(gsl_vector * x, double * fx), gsl_vector * x, gsl_vector * grad) {
	double fx, faux; int n = x->size;
	f(x,&fx);
	for (int i = 0; i < n; i++) {
		gsl_vector_set(x,i,gsl_vector_get(x,i)+DELTA);
		f(x,&faux);
		faux -= fx;
		gsl_vector_set(grad,i,faux/DELTA);
		gsl_vector_set(x,i,gsl_vector_get(x,i)-DELTA);
	}
}

void testfunc(gsl_vector * vars, double * fx) {
	double x = gsl_vector_get(vars,0);
	double y = gsl_vector_get(vars,1);
	*fx = x*x+y*y + 1.0;
}

void rosenbrock(gsl_vector * vars, double * fx) {
	double x = gsl_vector_get(vars,0);
	double y = gsl_vector_get(vars,1);
	* fx = pow(1-x,2)+100*pow(y-pow(x,2),2);
}

void himmelblau(gsl_vector * vars, double * fx) {
	double x = gsl_vector_get(vars,0);
	double y = gsl_vector_get(vars,1);
	* fx = pow(x*x+y-11,2)+pow(x+y*y-7,2);
}

void higgs(gsl_vector * vars, double * fx) {
	int m = 30, n = 3;
	double breitwigner; *fx = 0;
	
	FILE * data = fopen("higgsdata.txt","r");
	gsl_matrix * datmat = gsl_matrix_alloc(m,n);
	gsl_matrix_fscanf(data,datmat);
	
	double mass = gsl_vector_get(vars,0);
	double Gamma = gsl_vector_get(vars,1);
	double A = gsl_vector_get(vars,2);
	double E, sigmai, dsigmai;
	
	for (int i = 0; i < m; i++) {
		E = gsl_matrix_get(datmat,i,0);
		sigmai = gsl_matrix_get(datmat,i,1);
		dsigmai = gsl_matrix_get(datmat,i,2);
		breitwigner = A / (pow(E-mass,2)+pow(Gamma,2)/4);
		* fx += pow(breitwigner-sigmai,2)/pow(dsigmai,2);
	}
	
	gsl_matrix_free(datmat);
	fclose(data);
}

void higgsplot(gsl_vector * x) {
	FILE * data = fopen("higgsdata.txt","r");
	FILE * plotfile = fopen("plotdata.txt","w");
	gsl_matrix * datmat = gsl_matrix_alloc(30,3);
	gsl_matrix_fscanf(data,datmat);
	
	double mass = gsl_vector_get(x,0);
	double Gamma = gsl_vector_get(x,1);
	double A = gsl_vector_get(x,2);
	double E, sigmai, dsigmai;
	
	fprintf(plotfile,"#index 0\n");
	
	for (int i = 0; i < 30; i++) {
		E = gsl_matrix_get(datmat,i,0);
		sigmai = gsl_matrix_get(datmat,i,1);
		dsigmai = gsl_matrix_get(datmat,i,2);
		fprintf(plotfile,"%10g %10g %10g\n",E,sigmai,dsigmai);
	}
	
	fprintf(plotfile,"\n\n\n#index 1\n");
	double Ein = gsl_matrix_get(datmat,0,0);
	double Efin = gsl_matrix_get(datmat,29,0);
	int numpoints = 300;
	double h = (Efin - Ein) / (numpoints - 1.0);
	for (int i = 0; i < numpoints; i++) {
		fprintf(plotfile,"%10g %10g\n",Ein+h*i,A/(pow((Ein+h*i)-mass,2)+pow(Gamma,2)/4));
	}
	
	gsl_matrix_free(datmat);
	fclose(plotfile);
	fclose(data);
}

int qnewton(void (*f)(gsl_vector * x, double * fx), gsl_vector * x, double acc) {
	int n = x->size, numsteps = 0; double alpha;
	double epsilon = 1e-6;
	gsl_matrix * B = gsl_matrix_alloc(n,n);
	gsl_vector * gradfx = gsl_vector_alloc(n);
	gsl_vector * gradfxs = gsl_vector_alloc(n);
	gsl_vector * y = gsl_vector_alloc(n);
	gsl_vector * u = gsl_vector_alloc(n);
	gsl_vector * s = gsl_vector_alloc(n);
	gsl_vector * xs = gsl_vector_alloc(n);
	gsl_matrix_set_identity(B);
	gradient(f,x,gradfx);
	while (numsteps < 1e5) {
		numsteps++;
		if (gsl_blas_dnrm2(gradfx) < acc) break;
		gsl_blas_dgemv(CblasNoTrans,-1,B,gradfx,0,s);
		armijo(f,x,gradfx,s);
		gsl_vector_memcpy(xs,x);
		gsl_vector_add(xs,s);
		gradient(f,xs,y); 
		gsl_vector_memcpy(gradfxs,y);
		gsl_vector_sub(y,gradfx);
		gsl_vector_memcpy(u,s);
		gsl_blas_dgemv(CblasNoTrans,-1,B,y,1,u);
		gsl_blas_ddot(u,y,&alpha);
		if (alpha > epsilon) sr1(B,u,y);
		gsl_vector_memcpy(x,xs);
		gsl_vector_memcpy(gradfx,gradfxs);
	}
	
	
	gsl_vector_free(gradfx);
	gsl_vector_free(gradfxs);
	gsl_matrix_free(B);
	gsl_vector_free(s);
	gsl_vector_free(u);
	gsl_vector_free(y);
	gsl_vector_free(xs);

return numsteps;
}

int main() {
	double acc = 0.00001;
	gsl_vector * x = gsl_vector_alloc(2);
	gsl_vector_set(x,0,1.2);
	gsl_vector_set(x,1,-1.2);
	int steps = qnewton(rosenbrock,x,acc);
	vector_print("Minimum of Rosenbrock's valley function = ",x);
	printf("within accuracy %10g\n",acc);
	printf("number of steps %i\n\n\n",steps);
	
	gsl_vector * y = gsl_vector_alloc(2);
	gsl_vector_set(y,0,-3.0);
	gsl_vector_set(y,1,3);
	steps = qnewton(himmelblau,y,acc);
	vector_print("Minimum of Himmelblau's function = ",y);
	printf("within accuracy %10g\n",acc);
	printf("number of steps %i\n\n\n",steps);
	
	acc = 0.00001;
	gsl_vector * z = gsl_vector_alloc(3);
	double mass = 125.0, width = 2.09, scale = 2.0;
	gsl_vector_set(z,0,mass);
	gsl_vector_set(z,1,width);
	gsl_vector_set(z,2,scale);
	steps = qnewton(higgs,z,acc);
	vector_print("Minimum of Higgs in the format [E,Gamma,A]^T = ",z);
	printf("within accuracy %10g\n",acc);
	printf("number of steps %i\n",steps);
	higgsplot(z);
	
	
	gsl_vector_free(x);
	gsl_vector_free(y);
	gsl_vector_free(z);
		
return 0;
}










