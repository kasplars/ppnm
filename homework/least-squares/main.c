#include <stdio.h>
#include <math.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_linalg.h>
#include "auxfile.h"
#include <assert.h>

void least_sq(gsl_vector * x, gsl_vector * y, gsl_vector * dy,
	gsl_vector * c,double (*f)(double, int),int m,gsl_matrix * S) {
	int n = x->size;
	gsl_matrix * A = gsl_matrix_alloc(n,m);
	gsl_matrix * R = gsl_matrix_alloc(m,m);
	gsl_vector * b = gsl_vector_alloc(n);
	for (int i=0; i<n; i++) {
		gsl_vector_set(b,i,gsl_vector_get(y,i)/gsl_vector_get(dy,i));
		for (int j=0; j<m; j++) {
			gsl_matrix_set(A,i,j,f(gsl_vector_get(x,i),j)/gsl_vector_get(dy,i));
		}
	}

	GS_decomp(A,R);
	gsl_blas_dgemv(CblasTrans,1,A,b,0,c);
	backsub(R,c);

	// Covariance matrix
	gsl_matrix * R_inv = gsl_matrix_alloc(m,m);
	gsl_matrix * I = gsl_matrix_alloc(m,m);
	gsl_matrix_set_identity(I);
	GS_inverse(I,R,R_inv);
	gsl_blas_dgemm(CblasNoTrans,CblasTrans,1,R_inv,R_inv,0,S);

	gsl_matrix_free(I);
	gsl_matrix_free(R_inv);
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
	FILE* inputstream = fopen("input.txt","r");
	int lof = numlines("input.txt");

	int m=2;

	double doubley, doubledy, doublex;
	gsl_vector * x=gsl_vector_alloc(lof);
	gsl_vector * y=gsl_vector_alloc(lof);
	gsl_vector * dy=gsl_vector_alloc(lof);
	gsl_vector * c = gsl_vector_alloc(m);
	gsl_matrix * S = gsl_matrix_alloc(m,m);

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

	least_sq(x,y,dy,c,f,m,S);
	vector_print("c =",c);
	matrix_print("The covariance matrix S =",S);

	int datapoints = 200;
	double xmin = 0, xmax = 16;

	gsl_vector * xs = gsl_vector_alloc(datapoints);
	gsl_vector * ys = gsl_vector_alloc(datapoints);
	gsl_vector * yspc = gsl_vector_alloc(datapoints);
	gsl_vector * ysmc = gsl_vector_alloc(datapoints);


	FILE* plottingdata = fopen("plottingdata.txt","w");

	for(int i=0; i<lof; i++) {
		gsl_vector_set(xs,i,(xmax-xmin)*i/(datapoints-1));
		gsl_vector_set(ys,i,gsl_vector_get(c,0)+gsl_vector_get(xs,i)*gsl_vector_get(c,1));
		gsl_vector_set(yspc,i,gsl_vector_get(c,0)+sqrt(gsl_matrix_get(S,0,0))+gsl_vector_get(xs,i)*(gsl_vector_get(c,1)+sqrt(gsl_matrix_get(S,1,1))));
		gsl_vector_set(ysmc,i,gsl_vector_get(c,0)-sqrt(gsl_matrix_get(S,0,0))+gsl_vector_get(xs,i)*(gsl_vector_get(c,1)-sqrt(gsl_matrix_get(S,1,1))));
		fprintf(plottingdata,"%10g %10g %10g %10g %10g %10g %10g\n",gsl_vector_get(xs,i),gsl_vector_get(ys,i),gsl_vector_get(yspc,i),
		gsl_vector_get(ysmc,i),gsl_vector_get(x,i),gsl_vector_get(y,i),gsl_vector_get(dy,i));
	}
	for(int i=lof; i<datapoints; i++) {
		gsl_vector_set(xs,i,(xmax-xmin)*i/(datapoints-1));
		gsl_vector_set(ys,i,gsl_vector_get(c,0)+gsl_vector_get(xs,i)*gsl_vector_get(c,1));
		gsl_vector_set(yspc,i,gsl_vector_get(c,0)+sqrt(gsl_matrix_get(S,0,0))+gsl_vector_get(xs,i)*(gsl_vector_get(c,1)+sqrt(gsl_matrix_get(S,1,1))));
		gsl_vector_set(ysmc,i,gsl_vector_get(c,0)-sqrt(gsl_matrix_get(S,0,0))+gsl_vector_get(xs,i)*(gsl_vector_get(c,1)-sqrt(gsl_matrix_get(S,1,1))));
		fprintf(plottingdata,"%10g %10g %10g %10g\n",gsl_vector_get(xs,i),gsl_vector_get(ys,i),gsl_vector_get(yspc,i),
		gsl_vector_get(ysmc,i));
	}

	printf("The halflife of ThX is %g +- %g days. The modern value is approximately 3.6 days, so it is somewhat outside the expectation.\n",-log(2)/gsl_vector_get(c,1),log(2)/pow(gsl_vector_get(c,1),2)*sqrt(gsl_matrix_get(S,1,1)));

	

	gsl_vector_free(x);
	gsl_vector_free(y);
	gsl_vector_free(dy);
	gsl_vector_free(c);
	gsl_vector_free(xs);
	gsl_vector_free(ys);
	gsl_vector_free(yspc);
	gsl_vector_free(ysmc);
	gsl_matrix_free(S);

	fclose(plottingdata);
	fclose(inputstream);
return 0;
}
