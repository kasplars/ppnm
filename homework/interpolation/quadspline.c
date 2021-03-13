#include <stdio.h>
#include <math.h>
#include <assert.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_vector.h>

int binsearch(int n, double * x, double z){/* locates the interval for z by bisection */ 
	assert(n>1 && x[0] <=z && z<=x[n-1]);
	int i=0, j=n-1;
	while(j-i>1){
		int mid=(i+j)/2;
		if(z>x[mid]) i=mid; else j=mid;
		}
	return i;
}

typedef struct {int n; double *x, *y, *b, *c;} qspline;

qspline * qspline_alloc(int n, double *x, double *y){
	qspline * s = (qspline*) malloc(sizeof(qspline));
	s->c = (double*) malloc((n-1)*sizeof(double));
	s->b = (double*) malloc((n-1)*sizeof(double));
	s->x = (double*) malloc(n*sizeof(double));
	s->y = (double*) malloc(n*sizeof(double));
	s->n = n;

	int i;
	for (i=0;i<n;i++){s->x[i]=x[i],s->y[i]=y[i];}
	double dx[n-1], p[n-1];
	for (i=0;i<n-1;i++) {dx[i]=x[i+1]-x[i]; p[i]=(y[i+1]-y[i])/dx[i];}
	// Forward and backward recursion as in (1.11) and (1.12):
	s->c[0]=0;
	for (i=0;i<n-2;i++) {s->c[i+1] = 1/dx[i+1] * (p[i+1]-p[i]-s->c[i]*dx[i]);}
	s->c[n-2]/=2;
	for (i=n-3;i>=0;i--) {s->c[i] = 1/dx[i] * (p[i+1]-p[i]- s->c[i+1]*dx[i+1]);}
	for(i=0;i<n-1;i++) {s->b[i]=p[i]-s->c[i]*dx[i];}
	return s;
}

double qspline_eval(qspline * s, double z) {
	assert(z>=s->x[0] && z<=s->x[s->n-1]);
	int i = binsearch(s->n,s->x,z);
	double h=z-s->x[i];
	return s->y[i]+h*(s->b[i]+h*s->c[i]);
}

double qspline_deriv(qspline * s, double z) {
	assert(z>=s->x[0] && z<=s->x[s->n-1]);
	int i = binsearch(s->n,s->x,z); double dz = z-s->x[i];
	return s->b[i]+2*s->c[i]*dz;			// diff. (1.7)
}

double qspline_integ(qspline * s, double z) {
	assert(z>=s->x[0] && z<=s->x[s->n-1]);
	int j = binsearch(s->n,s->x,z);
	double sum = 0, dx;
	for (int i = 0; i<j;i++) {
		dx = s->x[i+1]-s->x[i];
		sum += s->y[i]*dx+s->b[i]* 1./2 * pow(dx,2) + s->c[i]*1./3*pow(dx,3);
	}
	double dz = z-s->x[j];
	sum += s->y[j]*dz+s->b[j]/2 * pow(dz,2) + s->c[j]/3*pow(dz,3);
	return sum;
}

void qspline_free(qspline * s) {
	free(s->x); free(s->y); free(s->b); free(s->c); free(s);
}

int main() {
	int n = 5;
	double xs[n], ya[n], yb[n], yc[n];
	for (int i = 0; i<n; i++) {
		xs[i]=i+1;  ya[i]=1; yb[i]=i+1; yc[i]=pow(i+1,2);
	}

	/* Chose which dataset:
		1:	{xi=i,yi=1}
		2: 	{xi=i,yi=i}
		3:	{xi=i,yi=i^2}
	*/

	int dataset = 3;

	//printf("%g\n",ya[2]);

	qspline * sa = qspline_alloc(n,xs,ya);
	qspline * sb = qspline_alloc(n,xs,yb);
	qspline * sc = qspline_alloc(n,xs,yc);

	int numPoints = 500; double z;

	FILE* outquadspline = fopen("outquadspline.txt","w");

	switch (dataset) {
		case 1:
			for (int i = 0; i<numPoints; i++) {
				z = xs[0]+(xs[n-1]-xs[0])*i/(numPoints-1);
				fprintf(outquadspline,"%10g %10g %10g %10g\n",
				z,qspline_eval(sa,z),qspline_deriv(sa,z),
				qspline_integ(sa,z));
			}
			break;
		case 2:
			for (int i = 0; i<numPoints; i++) {
				z = xs[0]+(xs[n-1]-xs[0])*i/(numPoints-1);
				fprintf(outquadspline,"%10g %10g %10g %10g\n",
				z,qspline_eval(sb,z),qspline_deriv(sb,z),
				qspline_integ(sb,z));
			}
			break;
		case 3:
			for (int i = 0; i<numPoints; i++) {
				z = xs[0]+(xs[n-1]-xs[0])*i/(numPoints-1);
				fprintf(outquadspline,"%10g %10g %10g %10g\n",
				z,qspline_eval(sc,z),qspline_deriv(sc,z),
				qspline_integ(sc,z));
			}
			break;
		default:
			printf("Choose a proper dataset.\n");
			break;
	}

	fclose(outquadspline);

	qspline_free(sa);
	qspline_free(sb);
	qspline_free(sc);
return 0;
}
