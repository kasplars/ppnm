#include <stdio.h>
#include <math.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_linalg.h>
#include "auxfile.h"
#include <assert.h>

void backsub(gsl_matrix * U, gsl_vector * c) {
	for(int i=c->size-1; i>=0; i--) {
		double s=gsl_vector_get(c,i);
		for(int k=i+1; k<c->size; k++) s-=gsl_matrix_get(U,i,k)*gsl_vector_get(c,k);
		gsl_vector_set(c,i,s/gsl_matrix_get(U,i,i));
	}
}

void vector_print(char s[],gsl_vector* v){
	printf("%s\n",s);
	for(int i=0;i < v->size;i++) printf("[%10g] \n",gsl_vector_get(v,i));
}

void matrix_print(char s[],gsl_matrix* A){
	printf("%s\n",s);
	for(int i=0; i < A->size1; i++) {
		printf("[");
		for(int j=0; j < A->size2; j++) {
			printf("%10g ",gsl_matrix_get(A,i,j));
		}
	printf("]\n");
	}
}

void randomizer_matrix(gsl_matrix * A, int modulo) {
	for (int i = 0; i < A->size1; i++) {
		for (int j = 0; j < A->size2; j++) {
			gsl_matrix_set(A,i,j,rand() % modulo - modulo/2);
		}
	}
}

void randomizer_vector(gsl_vector * A, int modulo) {
	for (int i = 0; i < A->size; i++) {
		gsl_vector_set(A,i,rand() % modulo - modulo/2);
	}
}

double normColumn(gsl_matrix * A, double j) { // jth column of A
	double normSq = 0;
	for (int i = 0; i< A->size1; i++) { // ith row of A
		normSq += pow(gsl_matrix_get(A,i,j),2);
	}
	return sqrt(normSq);
}

void GS_decomp(gsl_matrix * A, gsl_matrix * R) {
	for (int i = 0; i < A->size2; i++) { // for each column in A
		gsl_matrix_set(R,i,i,normColumn(A,i));	// Set Rii
		for (int k = 0; k<A->size1; k++) {	// Set qi=ai/Rii
			gsl_matrix_set(A,k,i,gsl_matrix_get(A,k,i)/gsl_matrix_get(R,i,i));
		}
		if (i+1<A->size2) {
			for (int j = i+1; j < A->size2; j++) {
				double Rij = 0;
				for (int index = 0; index < A->size1; index++) {
					Rij += gsl_matrix_get(A,index,i)*gsl_matrix_get(A,index,j);
				}
				gsl_matrix_set(R,i,j,Rij);
				for (int index = 0; index < A->size1; index++) {
					gsl_matrix_set(A,index,j,gsl_matrix_get(A,index,j)-gsl_matrix_get(A,index,i)*Rij);
				}
			}
		}
	}
}

void GS_solve(gsl_matrix* Q, gsl_matrix* R, gsl_vector* b, gsl_vector* x) {
	gsl_blas_dgemv(CblasTrans,1,Q,b,0,x);
	backsub(R,x);
}

void GS_inverse(gsl_matrix * Q, gsl_matrix * R, gsl_matrix * B) {
	gsl_vector * x = gsl_vector_alloc(B->size1);
	gsl_vector * b = gsl_vector_alloc(B->size1);

	for (int i = 0;  i<B->size1; i++) {
		gsl_vector_set_basis(b,i);
		GS_solve(Q,R,b,x);
		for (int j = 0; j<B->size1; j++) {
			gsl_matrix_set(B,j,i,gsl_vector_get(x,j));
		}
	}
	gsl_vector_free(x);
	gsl_vector_free(b);
}
