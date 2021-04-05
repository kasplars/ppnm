#include <stdio.h>
#include <stdlib.h>
#include <float.h>
#include <math.h>
#include <time.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_linalg.h>
#include "auxfile.h"

#ifndef NDEBUG
#define TRACE(args...) fprintf(stderr,args)
#else
#define TRACE(...)
#endif

void timesJ(gsl_matrix* A, int p, int q, double theta){
	double c=cos(theta),s=sin(theta);
	for(int i=0;i<A->size1;i++){
		double new_aip=c*gsl_matrix_get(A,i,p)-s*gsl_matrix_get(A,i,q);
		double new_aiq=s*gsl_matrix_get(A,i,p)+c*gsl_matrix_get(A,i,q);
		gsl_matrix_set(A,i,p,new_aip);
		gsl_matrix_set(A,i,q,new_aiq);
		}
}

void Jtimes(gsl_matrix* A, int p, int q, double theta){
	double c=cos(theta),s=sin(theta);
	for(int j=0;j<A->size2;j++){
		double new_apj= c*gsl_matrix_get(A,p,j)+s*gsl_matrix_get(A,q,j);
		double new_aqj=-s*gsl_matrix_get(A,p,j)+c*gsl_matrix_get(A,q,j);
		gsl_matrix_set(A,p,j,new_apj);
		gsl_matrix_set(A,q,j,new_aqj);
		}
}


double largestoffdiag(gsl_matrix * A) {
	double num = 0;
	int n = A->size1;
	for (int i = 0; i < n; i++) {
		for (int j = 0; j<n; j++) {
			if (gsl_matrix_get(A,i,j)>num && i!=j) {
				num = fabs(gsl_matrix_get(A,i,j));
			}
		}
	}
return num;
}

void jacobi_diag(gsl_matrix * A,gsl_matrix * V) {
	int n = A->size1; gsl_matrix_set_identity(V);
	double threshold = 10000*LDBL_EPSILON;
	do {
		for(int p=0;p<n-1;p++) {
			for(int q=p+1;q<n;q++){
			double apq=gsl_matrix_get(A,p,q);
			double app=gsl_matrix_get(A,p,p);
			double aqq=gsl_matrix_get(A,q,q);
			double theta=0.5*atan2(2*apq,aqq-app);
			//double c=cos(theta),s=sin(theta);
			//double new_app=c*c*app-2*s*c*apq+s*s*aqq;
			//double new_aqq=s*s*app+2*s*c*apq+c*c*aqq;
			timesJ(A,p,q, theta);
			Jtimes(A,p,q,-theta); // A<-J^T*A*J
			timesJ(V,p,q, theta); // V<-V*J
			}
		}
	}while(largestoffdiag(A)>threshold); //while(changed!=0);
}

void gen_rel_sym(gsl_matrix * S,double range) {
	srand(time(NULL));
	for (int i = 0; i < S->size1; i++) {
		for (int j = i; j < S->size2; j++) {
			double div = RAND_MAX / (2*range);
			gsl_matrix_set(S,i,j,-range + (rand() / div));
		}
	}
	for (int i = 1; i < S->size1; i++) {
		for (int j = 0; j<i; j++) {
			gsl_matrix_set(S,i,j,gsl_matrix_get(S,j,i));
		}
	}
}


int main() {
	int n = 5;
	gsl_matrix * S = gsl_matrix_alloc(n,n);
	gsl_matrix * V = gsl_matrix_alloc(n,n);
	gsl_matrix * vdvt = gsl_matrix_alloc(n,n);
	gsl_matrix * D = gsl_matrix_alloc(n,n);
	gsl_matrix * vtsv = gsl_matrix_alloc(n,n);
	gsl_matrix * temp = gsl_matrix_alloc(n,n);

	gen_rel_sym(S,1);
	matrix_print("S=",S);
	gsl_matrix_memcpy(D,S);
	jacobi_diag(D,V);
	matrix_print("D=",D);
	matrix_print("V=",V);
	gsl_blas_dgemm(CblasTrans,CblasNoTrans,1,V,S,0,temp);
	gsl_blas_dgemm(CblasNoTrans,CblasNoTrans,1,temp,V,0,vtsv);
	matrix_print("Check that V^T*S*V = D: ",vtsv);
	gsl_blas_dgemm(CblasNoTrans,CblasNoTrans,1,V,D,0,temp);
	gsl_blas_dgemm(CblasNoTrans,CblasTrans,1,temp,V,0,vdvt);
	matrix_print("Check that V*D*V^T = S: ",vdvt);
	gsl_blas_dgemm(CblasTrans,CblasNoTrans,1,V,V,0,temp);
	matrix_print("Check that V^T*V = 1: ",temp);

	gsl_matrix_free(S);
	gsl_matrix_free(V);
	gsl_matrix_free(vdvt);
	gsl_matrix_free(D);
	gsl_matrix_free(vtsv);
	gsl_matrix_free(temp);
return 0;
}
