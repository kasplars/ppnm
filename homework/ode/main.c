#include <stdio.h>
#include <math.h>
#include <time.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_linalg.h>
#include "library.h"

#ifndef NDEBUG
#define TRACE(args...) fprintf(stderr,args)
#else
#define TRACE(...)
#endif

void shm(double t, gsl_vector * y, gsl_vector * dydt) {
	gsl_vector_set(dydt,0,gsl_vector_get(y,1));
	gsl_vector_set(dydt,1,-gsl_vector_get(y,0));
}

void sir(double t,gsl_vector * y, gsl_vector * dydt) {
	double Tc = 3.5, Tr = 20.0;
	double N;

	double S = gsl_vector_get(y,0), I = gsl_vector_get(y,1), R = gsl_vector_get(y,2);
	N = S + I + R;
	gsl_vector_set(dydt,0,-I*S/N/Tc);
	gsl_vector_set(dydt,1,I*(S/N/Tc - 1.0/Tr));
	gsl_vector_set(dydt,2,I/Tr);
}

void rkstep45(void (*f)(double, gsl_vector *, gsl_vector *), /* the f from dy/dt=f(t,y) */
	int n,				/* step number	*/
	double t,              		/* the current value of the variable */
	gsl_vector * yt,            	/* the current value y(t) of the sought function */
	double h,              		/* the step to be taken */
	gsl_vector * yh,             	/* output: y(t+h) */
	gsl_vector * dy             	/* output: error estimate */
){
	int i;
	gsl_vector * k1 = gsl_vector_alloc(n);
	gsl_vector * k2 = gsl_vector_alloc(n);
	gsl_vector * k3 = gsl_vector_alloc(n);
	gsl_vector * k4 = gsl_vector_alloc(n);
	gsl_vector * k5 = gsl_vector_alloc(n);
	gsl_vector * k6 = gsl_vector_alloc(n);
	gsl_vector * yn = gsl_vector_alloc(n);

	double b[7] = {0.0, 16.0/135, 0.0, 6656.0/12825, 28561.0/56430, -9.0/50, 2.0/55};
	double bstar[7] = {0.0, 25.0/216, 0.0, 1408.0/2565, 2197.0/4104, -1.0/5, 0.0};
	double c[7] = {0.0, 0.0, 1.0/4, 3.0/8, 12.0/13, 1.0, 1.0/2};
	double a1[7] = {0.0, 0.0, 1.0/4, 3.0/32, 1932.0/2197, 439.0/216, -8.0/27};
	double a2[7] = {0.0, 0.0, 0.0, 9.0/32, -7200.0/2197, -8.0, 2.0};
	double a3[7] = {0.0, 0.0, 0.0, 0.0, 7296.0/2197, 3680.0/513, -3544.0/2565};
	double a4[7] = {0.0, 0.0, 0.0, 0.0, 0.0, -845.0/4104, 1859.0/4104};
	double a5[7] = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, -11.0/40};

	TRACE("After defining coefficients\n");

	f(t,yt,k1);
	TRACE("After using function\n");
	for(i=0;i<n;i++) {gsl_vector_set(yn,i,gsl_vector_get(yt,i)+(a1[2]*gsl_vector_get(k1,i))*h);}
	f(t+c[2]*h,yn,k2);
	TRACE("After using one for loop\n");
	for(i=0;i<n;i++) {gsl_vector_set(yn,i,gsl_vector_get(yn,i)+(a1[3]*gsl_vector_get(k1,i)+a2[3]*gsl_vector_get(k2,i))*h);}
	f(t+c[3]*h,yn,k3);
	for(i=0;i<n;i++) {gsl_vector_set(yn,i,gsl_vector_get(yn,i)+(a1[4]*gsl_vector_get(k1,i)+a2[4]*gsl_vector_get(k2,i)+a3[4]*gsl_vector_get(k3,i))*h);}
	f(t+c[4]*h,yn,k4);
	for(i=0;i<n;i++) {gsl_vector_set(yn,i,gsl_vector_get(yn,i)+(a1[5]*gsl_vector_get(k1,i)+a2[5]*gsl_vector_get(k2,i)+a3[5]*gsl_vector_get(k3,i)+a4[5]*gsl_vector_get(k4,i))*h);}
	f(t+c[5]*h,yn,k5);
	for(i=0;i<n;i++) {gsl_vector_set(yh,i,gsl_vector_get(yn,i)+(a1[6]*gsl_vector_get(k1,i)+a2[6]*gsl_vector_get(k2,i)+a3[6]*gsl_vector_get(k3,i)+a4[6]*gsl_vector_get(k4,i)+a5[6]*gsl_vector_get(k5,i))*h);}
	f(t+c[6]*h,yh,k6);
	for(i=0;i<n;i++) {gsl_vector_set(yh,i,gsl_vector_get(yt,i)+(b[1]*gsl_vector_get(k1,i)
				+b[2]*gsl_vector_get(k2,i)+b[3]*gsl_vector_get(k3,i)
				+b[4]*gsl_vector_get(k4,i)+b[5]*gsl_vector_get(k5,i)+b[6]*gsl_vector_get(k6,i))*h);
			gsl_vector_set(yn,i,gsl_vector_get(yt,i)+(bstar[1]*gsl_vector_get(k1,i)
				+bstar[2]*gsl_vector_get(k2,i)+bstar[3]*gsl_vector_get(k3,i)
				+bstar[4]*gsl_vector_get(k4,i)+bstar[5]*gsl_vector_get(k5,i)+bstar[6]*gsl_vector_get(k6,i))*h);
			gsl_vector_set(dy,i,gsl_vector_get(yh,i)-gsl_vector_get(yn,i));
	}

	gsl_vector_free(k1);
	gsl_vector_free(k2);
	gsl_vector_free(k3);
	gsl_vector_free(k4);
	gsl_vector_free(k5);
	gsl_vector_free(k6);
	gsl_vector_free(yn);
}

void odeprint(FILE * output,double x, gsl_vector * y) {
	int n = y->size;
	if (n == 1) fprintf(output,"%10g %10g\n",x,gsl_vector_get(y,0));
	if (n == 2) fprintf(output,"%10g %10g %10g\n",x,gsl_vector_get(y,0),gsl_vector_get(y,1));
	if (n == 3) fprintf(output,"%10g %10g %10g %10g\n",x,gsl_vector_get(y,0),gsl_vector_get(y,1),gsl_vector_get(y,2));
	if (n == 4) fprintf(output,"%10g %10g %10g %10g %10g\n",x,gsl_vector_get(y,0),gsl_vector_get(y,1),gsl_vector_get(y,2),
								gsl_vector_get(y,3));
	if (n == 5) fprintf(output,"%10g %10g %10g %10g %10g %10g\n",x,gsl_vector_get(y,0),gsl_vector_get(y,1),
								gsl_vector_get(y,2),gsl_vector_get(y,3),gsl_vector_get(y,4));
	if (n == 6) fprintf(output,"%10g %10g %10g %10g %10g %10g %10g\n",x,gsl_vector_get(y,0),gsl_vector_get(y,1),
								gsl_vector_get(y,2),gsl_vector_get(y,3),gsl_vector_get(y,4),
								gsl_vector_get(y,5));
}

int driver(FILE * outstream,void (*f)(double t, gsl_vector * y,gsl_vector * dydt), /* right-hand-side of dy/dt=f(t,y) */
	int n,
	double a,                     /* the start-point a */
	gsl_vector * ya,                     /* y(a) */
	double b,                     /* the end-point of the integration */
	gsl_vector * yb,                     /* y(b) to be calculated */
	double h,                     /* initial step-size */
	double acc,                   /* absolute accuracy goal */
	double eps                    /* relative accuracy goal */
){
	double x = a, ei, normy, s, taui;
	if (x+h>b) h = b-x;

	int i, k = 0;
	TRACE("Before vector init\n");
	gsl_vector * dy = gsl_vector_alloc(n);
	TRACE("After vector init\n");

	while(x<b) {
		TRACE("Before rkstep\n");
		rkstep45(f,n,x,ya,h,yb,dy); k++;
		TRACE("After rkstep\n");
		s=0; for (i=0; i<n; i++) s+=pow(gsl_vector_get(dy,i),2); ei = sqrt(s);
		s=0; for (i=0; i<n; i++) s+=pow(gsl_vector_get(yb,i),2); normy = sqrt(s);
		taui = (normy*eps+acc)*sqrt(h/(b-a));
		if (ei<taui){
			x += h; for(i=0;i<n;i++) {gsl_vector_memcpy(ya,yb); odeprint(outstream,x,yb);}
		} if (ei>0) h*=pow(taui/ei,0.25)*0.95; else h*=2;
	}

	gsl_vector_free(dy);
	return k;
}

int main(){
	int n = 2;
	gsl_vector * ya = gsl_vector_alloc(2);
	gsl_vector * yb = gsl_vector_alloc(2);

	gsl_vector_set(ya,0,1); gsl_vector_set(ya,1,0);

	double a = 0, b = 10;
	double h = 0.002;
	double acc = 0.00001;
	double eps = 0.00001;

	FILE* outstream = fopen("ode-out.txt","w");

	int steps = driver(outstream,shm,n,a,ya,b,yb,h,acc,eps);
	vector_print("yb=",yb); printf("Number of steps: %i\n",steps);

	fclose(outstream);

	gsl_vector_free(ya);
	gsl_vector_free(yb);

	// SIR

	n = 3;

	gsl_vector * sirvec = gsl_vector_alloc(n);
	gsl_vector * dsir = gsl_vector_alloc(n);

	double sirvals[3] = {10000.0, 10.0, 100.0};

	for (int i = 0; i<n; i++) gsl_vector_set(sirvec,i,sirvals[i]);

	FILE * sirstream = fopen("sir.txt","w");

	a = 0; b = 200;
	int sirsteps = driver(sirstream,sir,n,a,sirvec,b,dsir,h,acc,eps);
	vector_print("dsir = ",dsir); printf("Number of steps: %i\n",sirsteps);

	fclose(sirstream);

	gsl_vector_free(sirvec);
	gsl_vector_free(dsir);
return 0;
}
