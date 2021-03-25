#include <stdio.h>
#include <math.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_linalg.h>
#include "auxfile.h"
#include <assert.h>

#ifndef NDEBUG
#define TRACE(args...) fprintf(stderr,args)
#else
#define TRACE(...)
#endif

void least_sq(gsl_vector * x, gsl_vector * y, gsl_vector * dy, gsl_vector * c,double (*f)(double, int),int m) {
	TRACE("0LS\n");
	int n = x->size;
	gsl_matrix * A = gsl_matrix_alloc(n,m);
	gsl_matrix * R = gsl_matrix_alloc(m,m);
	gsl_vector * b = gsl_vector_alloc(n);
	TRACE("1LS\n");
	for (int i=0; i<n; i++) {
		gsl_vector_set(b,i,gsl_vector_get(y,i)/gsl_vector_get(dy,i));
		for (int j=0; j<m; j++) {
			gsl_matrix_set(A,i,j,f(gsl_vector_get(x,i),j)/gsl_vector_get(dy,i));
		}
	}

	GS_decomp(A,R);
	TRACE("2LS\n");
	gsl_blas_dgemv(CblasTrans,1,A,b,0,c);
	TRACE("3LS\n");
	backsub(R,c);
	TRACE("4LS\n");

	gsl_matrix_free(A);
	gsl_matrix_free(R);
	gsl_vector_free(b);
}

void vector_transform(gsl_vector * x, double (*f)(double)) {
	int n = x->size;
	for (int i = 0; i < n; i++) gsl_vector_set(x,i,f(gsl_vector_get(x,i)));
}

int numlines(char * name) {
   	char ch;
   	int linesCount=0;
   	//open file in read more
   	FILE *fp=fopen(name,"r");
   	//read character by character and check for new line
   	while((ch=fgetc(fp))!=EOF) {
      		if(ch=='\n')
         	linesCount++;
   	}
   	//close the file
   	fclose(fp);
return linesCount;
}

double f(double x, int i) {
	switch(i) {
		case 0:
			return 1;
			break;
		case 1:
			return x;
			break;
		case 2: return x*x;
			break;
		default:
			return NAN;
			break;
	}
}

int main() {
	TRACE("0\n");
	FILE* inputstream = fopen("input.txt","r");
	int lof = numlines("input.txt");

	int m=2;

	double doubley, doubledy, doublex;
	gsl_vector * x=gsl_vector_alloc(lof);
	gsl_vector * y=gsl_vector_alloc(lof);
	gsl_vector * dy=gsl_vector_alloc(lof);
	gsl_vector * c = gsl_vector_alloc(m);

	for (int index = 0; index < lof; index++) {
		if(fscanf(inputstream,"%lg %lg %lg",&doublex,&doubley,&doubledy)); // to ignore warning;
		gsl_vector_set(x,index,doublex); gsl_vector_set(y,index,doubley); gsl_vector_set(dy,index,doubledy);
	}

	vector_print("x =",x);
	vector_print("y =",y);
	vector_print("dy =",dy);

	gsl_vector_div(dy,y); // takes dy -> dy / y = dln(y)
	vector_transform(y,log);

	vector_print("log(y)=",y);
	vector_print("dlog(y)=",dy);

	least_sq(x,y,dy,c,f,m);
	vector_print("c =",c);

	int datapoints = 200;
	double xmin = 0, xmax = 16;

	gsl_vector * xs = gsl_vector_alloc(datapoints);
	gsl_vector * ys = gsl_vector_alloc(datapoints);

	FILE* plottingdata = fopen("plottingdata.txt","w");

	for(int i=0; i<lof; i++) {
		gsl_vector_set(xs,i,(xmax-xmin)*i/(datapoints-1));
		gsl_vector_set(ys,i,gsl_vector_get(c,0)+gsl_vector_get(xs,i)*gsl_vector_get(c,1));
		fprintf(plottingdata,"%10g %10g %10g %10g %10g\n",gsl_vector_get(xs,i),gsl_vector_get(ys,i),
			gsl_vector_get(x,i),gsl_vector_get(y,i),gsl_vector_get(dy,i));
	}
	for(int i=lof; i<datapoints; i++) {
		gsl_vector_set(xs,i,(xmax-xmin)*i/(datapoints-1));
		gsl_vector_set(ys,i,gsl_vector_get(c,0)+gsl_vector_get(xs,i)*gsl_vector_get(c,1));
		fprintf(plottingdata,"%10g %10g\n",gsl_vector_get(xs,i),gsl_vector_get(ys,i));
	}

	printf("The halflife of ThX is %g days. The modern value is approximately 3.6 days.\n",-log(2)/gsl_vector_get(c,1));

	gsl_vector_free(x);
	gsl_vector_free(y);
	gsl_vector_free(dy);
	gsl_vector_free(c);
	gsl_vector_free(xs);
	gsl_vector_free(ys);

	fclose(plottingdata);
	fclose(inputstream);
return 0;
}
