#include <stdio.h>
#include <math.h>
#include <assert.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_interp.h>

int binsearch(int n, gsl_vector * x, double z){/* locates the interval for z by bisection */ 
	assert(n>1 && gsl_vector_get(x,0)<=z && z<=gsl_vector_get(x,n-1));
	int i=0, j=n-1;
	while(j-i>1){
		int mid=(i+j)/2;
		if(z>gsl_vector_get(x,mid)) i=mid; else j=mid;
		}
	return i;
}

double linterp(int n, gsl_vector * x, gsl_vector * y, double z) {
	assert(n>1 && z>gsl_vector_get(x,0) && z<gsl_vector_get(x,n-1));
	int i = binsearch(n,x,z);
	double dy=gsl_vector_get(y,i+1)-gsl_vector_get(y,i), dx=gsl_vector_get(x,i+1)-gsl_vector_get(x,i);
	assert(dx>0);
	return gsl_vector_get(y,i)+dy/dx*(z-gsl_vector_get(x,i));
}

double func(double x) {return pow(x,2);}

double linterp_integ(int n, gsl_vector * x, gsl_vector * y, double z) {
	assert(n>1 && z>gsl_vector_get(x,0) && z<gsl_vector_get(x,n-1));
	int j = binsearch(n,x,z);
	// printf("%i",j);
	double sum = 0, dy, dx, x1, x2, y1, y2;
	for (int i = 0; i<j; i++) {
		x1 = gsl_vector_get(x,i), x2 = gsl_vector_get(x,i+1);
		y1 = gsl_vector_get(y,i), y2 = gsl_vector_get(y,i+1);
		dy=y2-y1, dx=x2-x1;
		sum += dy/dx * 0.5 * pow(dx,2) + y1 * dx;
		/*
		The integral is as follows, with a=dy/dx and b=y1:
		\int_{x_{i}}^{x_{i+1}} a*(x-x1) + b dx = \int_{0}^{x_{i+1}-x_{i}} ax+b dx
		= [a/2 * x^2 + b * x]^{x_{i+1}-x_{i}}_{0},
		by substitution.
		*/
	}
	dy=gsl_vector_get(y,j+1)-gsl_vector_get(y,j), dx=gsl_vector_get(x,j+1)-gsl_vector_get(x,j); 
	double dz = z - gsl_vector_get(x,j);
	assert(dx>0);
	sum += dy/dx * 0.5 * pow(dz,2) + gsl_vector_get(y,j) * dz;
	return sum;
}

int main(){
	// Construct test-vectors based on interval and
	int n = 10000;
	double xmin = -10, xmax = 10, xs[n], ys[n];
	gsl_vector * x = gsl_vector_alloc(n);
	gsl_vector * y = gsl_vector_alloc(n);
	for (int i=0; i<n; i++) {
		gsl_vector_set(x,i,xmin + (xmax-xmin)*i/(n-1)); gsl_vector_set(y,i,func(gsl_vector_get(x,i)));
		xs[i] = xmin + (xmax-xmin)*i/(n-1); ys[i]=func(xs[i]);
		// printf("%10g %10g\n",gsl_vector_get(x,i),gsl_vector_get(y,i));
	}

	// Point to be evaluated, print and free memory of vectors
	double z = 5;
	printf("Interpolation of f(z) on the interval [%g,%g]. Check: f(z=%g)=%g\n",xmin,xmax,z,linterp(n,x,y,z));
	printf("Integration of f(z) on interval [%g,%g]. Check: %g\n",xmin,z,linterp_integ(n,x,y,z));

	FILE* outlspline = fopen("outlspline.txt","w");
	FILE* outinteg = fopen("outlinteg.txt","w");
	gsl_interp * linear = gsl_interp_alloc(gsl_interp_linear,n);
	gsl_interp_init(linear,xs,ys,n);

	double nxmin = xmin*0.9, nxmax = xmax*0.9;

	for (int i=0; i<n; i++) {
		z = nxmin + (nxmax-nxmin)*i/(n-1);
		// printf("%10g\n",z);
		double interp_l=gsl_interp_eval(linear,xs,ys,z,NULL);
		double integ_l=gsl_interp_eval_integ(linear,xs,ys,gsl_vector_get(x,0),z,NULL);
		fprintf(outlspline,"%10g %10g %10g\n",z,linterp(n,x,y,z),interp_l);
		fprintf(outinteg,"%10g %10g %10g\n",z,linterp_integ(n,x,y,z),integ_l);
	}

	fclose(outlspline);
	fclose(outinteg);

	gsl_vector_free(x);
	gsl_vector_free(y);
	gsl_interp_free(linear);

	return 0;
}


