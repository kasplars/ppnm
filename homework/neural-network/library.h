#ifndef LIBRARY_H
#define LIBRARY_H

void backsub(gsl_matrix * U, gsl_vector * c);
void vector_print(char s[],gsl_vector* v);
void matrix_print(char s[],gsl_matrix* A);
void randomizer_matrix(gsl_matrix * A, int modulo);
void randomizer_vector(gsl_vector * A, int modulo);
double normColumn(gsl_matrix * A, double j);
void GS_decomp(gsl_matrix * A, gsl_matrix * R);
void GS_solve(gsl_matrix* Q, gsl_matrix* R, gsl_vector* b, gsl_vector* x);
void GS_inverse(gsl_matrix * Q, gsl_matrix * R, gsl_matrix * B);
void timesJ(gsl_matrix* A, int p, int q, double theta);
void Jtimes(gsl_matrix* A, int p, int q, double theta);
double largestoffdiag(gsl_matrix * A);
void jacobi_diag(gsl_matrix * A,gsl_matrix * V);
void gen_rel_sym(gsl_matrix * S,double range);
void rkstep45(void (*f)(double, gsl_vector *, gsl_vector *), /* the f from dy/dt=f(t,y) */
	int n,				/* step number	*/
	double t,              		/* the current value of the variable */
	gsl_vector * yt,            	/* the current value y(t) of the sought function */
	double h,              		/* the step to be taken */
	gsl_vector * yh,             	/* output: y(t+h) */
	gsl_vector * dy             	/* output: error estimate */
);
void odeprint(FILE * output,double x, gsl_vector * y); // is  used in driver
int driver(FILE * outstream,void (*f)(double t, gsl_vector * y,gsl_vector * dydt), /* right-hand-side of dy/dt=f(t,y) */
	int n,				/* ya and yb dimension */
	double a,                     /* the start-point a */
	gsl_vector * ya,                     /* y(a) */
	double b,                     /* the end-point of the integration */
	gsl_vector * yb,                     /* y(b) to be calculated */
	double h,                     /* initial step-size */
	double acc,                   /* absolute accuracy goal */
	double eps                    /* relative accuracy goal */
);

// MINIMIZATION
int qnewton(void (*f)(gsl_vector * x, double * fx), gsl_vector * x, double acc);
void gradient(void (*f)(gsl_vector * x, double * fx), gsl_vector * x, gsl_vector * grad);
void armijo(void (*f)(gsl_vector * x, double * fx),gsl_vector * x, gsl_vector * gradfx, gsl_vector * dx);
void sr1(gsl_matrix * B, gsl_vector * u, gsl_vector * y);

#endif
