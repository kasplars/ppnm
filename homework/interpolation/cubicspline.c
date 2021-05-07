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

double func(double x) {if (x<0) return -1.0; else return 1.0;}

void cubicspline(gsl_vector * x, gsl_vector * y, gsl_vector * b, gsl_vector * c, gsl_vector * d) {
	int i, n = x->size;
	gsl_vector * h = gsl_vector_alloc(n-1);
	gsl_vector * p = gsl_vector_alloc(n-1);
	gsl_vector * D = gsl_vector_alloc(n);
	gsl_vector * B = gsl_vector_alloc(n);
	gsl_vector * Q = gsl_vector_alloc(n-1);
	
	for(i=0;i<n-1;i++) gsl_vector_set(h,i,gsl_vector_get(x,i+1)-gsl_vector_get(x,i));
	for(i=0;i<n-1;i++) gsl_vector_set(p,i,
			1.0/gsl_vector_get(h,i)*(gsl_vector_get(y,i+1)-gsl_vector_get(y,i)));
	
	// Making the tridiagonal system:
	gsl_vector_set(D,0,2); gsl_vector_set(D,n-1,2); gsl_vector_set(Q,0,1);
	gsl_vector_set(B,0,3*gsl_vector_get(p,0)); gsl_vector_set(B,n-1,3*gsl_vector_get(p,n-2));
	for(i=0;i<n-2;i++) gsl_vector_set(D,i+1,2*gsl_vector_get(h,i)/gsl_vector_get(h,i+1)+2);
	for(i=0;i<n-2;i++) gsl_vector_set(Q,i+1,gsl_vector_get(h,i)/gsl_vector_get(h,i+1));
	for(i=0;i<n-2;i++) gsl_vector_set(B,i+1,3*(gsl_vector_get(p,i)+gsl_vector_get(p,i+1)*gsl_vector_get(Q,i+1)));
	
	// Gauss elimination:
	for(i=1;i<n;i++){gsl_vector_set(D,i,gsl_vector_get(D,i)-gsl_vector_get(Q,i-1)/gsl_vector_get(D,i-1)); 
				gsl_vector_set(B,i,gsl_vector_get(B,i)-gsl_vector_get(B,i-1)/gsl_vector_get(D,i-1));}
	
	// Backsub:
	gsl_vector_set(b,n-1,gsl_vector_get(B,n-1)/gsl_vector_get(D,n-1));
	for(i=n-2;i>=0;i--) gsl_vector_set(b,i,(gsl_vector_get(B,i)-gsl_vector_get(Q,i)*gsl_vector_get(b,i+1))/gsl_vector_get(D,i));
	
	// eq. 20
	for(i=0;i<n-1;i++){
		gsl_vector_set(c,i,(-2*gsl_vector_get(b,i)-gsl_vector_get(b,i+1)+3*gsl_vector_get(p,i))/gsl_vector_get(h,i));
		gsl_vector_set(d,i,(gsl_vector_get(b,i)+gsl_vector_get(b,i+1)-2*gsl_vector_get(p,i))/pow(gsl_vector_get(h,i),2));
	}
	
	gsl_vector_free(h);
	gsl_vector_free(p);
	gsl_vector_free(D);
	gsl_vector_free(B);
	gsl_vector_free(Q);
}

double cubiceval(gsl_vector * x, gsl_vector * y, double z) {
	int n = x->size;
	gsl_vector * b = gsl_vector_alloc(n);
	gsl_vector * c = gsl_vector_alloc(n-1);
	gsl_vector * d = gsl_vector_alloc(n-1);
	cubicspline(x,y,b,c,d);
	int j = binsearch(n,x,z);
	double h = z-gsl_vector_get(x,j);
	double val = gsl_vector_get(y,j)+h*(gsl_vector_get(b,j)+h*(gsl_vector_get(c,j)+h*gsl_vector_get(d,j)));
	gsl_vector_free(b);
	gsl_vector_free(c);
	gsl_vector_free(d);
	
	return val;
}

double cubicinteg(gsl_vector * x, gsl_vector * y, double z) {
	int n = x->size;
	int j = binsearch(n,x,z);
	gsl_vector * b = gsl_vector_alloc(n);
	gsl_vector * c = gsl_vector_alloc(n-1);
	gsl_vector * d = gsl_vector_alloc(n-1);
	cubicspline(x,y,b,c,d);
	
	double sum = 0, dx;
	for (int i = 0; i<j;i++) {
		dx = gsl_vector_get(x,i+1)-gsl_vector_get(x,i);
		sum += gsl_vector_get(y,i)*dx+1.0/2*gsl_vector_get(b,i)*pow(dx,2)+1.0/3*gsl_vector_get(c,i)*pow(dx,3)
			+1.0/4*gsl_vector_get(d,i)*pow(dx,4);
	}
	dx = z - gsl_vector_get(x,j);
	sum += gsl_vector_get(y,j)*dx+1.0/2*gsl_vector_get(b,j)*pow(dx,2)+1.0/3*gsl_vector_get(c,j)*pow(dx,3)
			+1.0/4*gsl_vector_get(d,j)*pow(dx,4);
	
	gsl_vector_free(b);
	gsl_vector_free(c);
	gsl_vector_free(d);
	
	return sum;
}

double cubicderiv(gsl_vector * x, gsl_vector * y, double z) {
	int n = x->size, j = binsearch(n,x,z);
	gsl_vector * b = gsl_vector_alloc(n);
	gsl_vector * c = gsl_vector_alloc(n-1);
	gsl_vector * d = gsl_vector_alloc(n-1);
	cubicspline(x,y,b,c,d);
	
	double dx = z-gsl_vector_get(x,j);
	double val = gsl_vector_get(b,j)+2.0*gsl_vector_get(c,j)*dx+3.0*gsl_vector_get(d,j)*pow(dx,2);
	
	gsl_vector_free(b);
	gsl_vector_free(c);
	gsl_vector_free(d);
	
	return val;
}

int main(){
	int n = 6;
	double a = -3.5, b = 3.5;
	double dx = (b-a)/(n-1);
	gsl_vector * xs = gsl_vector_alloc(n);
	gsl_vector * ys = gsl_vector_alloc(n);
	
	for(int i = 0; i<n; i++) {
		gsl_vector_set(xs,i,a+i*dx);
		gsl_vector_set(ys,i,func(gsl_vector_get(xs,i)));
	}
	
	FILE* outcubicspline = fopen("outcubicspline.txt","w");
	
	int numpoints = 200; double z;
	
	for(int i = 0; i<numpoints; i++) {
		z = a+i*(b-a)/(numpoints-1);
		fprintf(outcubicspline,"%10g %10g %10g %10g\n",
				z,cubiceval(xs,ys,z),cubicderiv(xs,ys,z),cubicinteg(xs,ys,z));
	}
	
	fclose(outcubicspline);
	
	gsl_vector_free(xs);
	gsl_vector_free(ys);
return 0;
}
