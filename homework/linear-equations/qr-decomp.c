#include <stdio.h>
#include <math.h>
#include <time.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_linalg.h>

#ifndef NDEBUG
#define TRACE(args...) fprintf(stderr,args)
#else
#define TRACE(...)
#endif

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
		TRACE("2\n");
		gsl_matrix_set(R,i,i,normColumn(A,i));	// Set Rii
		TRACE("3\n");
		for (int k = 0; k<A->size1; k++) {	// Set qi=ai/Rii
			gsl_matrix_set(A,k,i,gsl_matrix_get(A,k,i)/gsl_matrix_get(R,i,i));
		}
		TRACE("4\n");
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

void backsub(gsl_matrix * U, gsl_vector * c) {
	for(int i=c->size-1; i>=0; i--) {
		double s=gsl_vector_get(c,i);
		for(int k=i+1; k<c->size; k++) s-=gsl_matrix_get(U,i,k)*gsl_vector_get(c,k);
		gsl_vector_set(c,i,s/gsl_matrix_get(U,i,i));
	}
}

void GS_solve(gsl_matrix* Q, gsl_matrix* R, gsl_vector* b, gsl_vector* x) {
	gsl_blas_dgemv(CblasTrans,1,Q,b,0,x);
	backsub(R,x);
}

int main() {
	srand(time(NULL));
	// PART 1

	int n = 8, m = 6;

	printf("------------------------\nPart 1\n\n");

	gsl_matrix * A = gsl_matrix_alloc(n,m);
	gsl_matrix * R = gsl_matrix_alloc(m,m);

	randomizer_matrix(A,50);
	matrix_print("The random matrix is A = ",A);
	printf("\n");
	GS_decomp(A,R);

	matrix_print("New orthogonal matrix Q from QR-decomposition: ",A);
	printf("\n");
	matrix_print("Check that R is upper triangular: ",R);
	printf("\n");

	gsl_matrix * temp = gsl_matrix_alloc(m,m);
	gsl_blas_dgemm(CblasTrans, CblasNoTrans,1,A,A,0,temp);
	matrix_print("Check that Q^T * Q = 1: ",temp);
	printf("\n");
	gsl_matrix_free(temp);


	gsl_matrix * temp2 = gsl_matrix_alloc(n,m);
	gsl_blas_dgemm(CblasNoTrans, CblasNoTrans,1,A,R,0,temp2);
	matrix_print("Check that QR=A: ",temp2);
	printf("\n");
	gsl_matrix_free(temp2);

	gsl_matrix_free(A);
	gsl_matrix_free(R);

	// PART 2
	printf("\n\n------------------------\nPart 2\n\n");

	n = 8;

	gsl_matrix * S = gsl_matrix_alloc(n,n);
	gsl_matrix * S_orthogonal = gsl_matrix_alloc(n,n);
	gsl_vector * b = gsl_vector_alloc(n);
	gsl_matrix * T = gsl_matrix_alloc(n,n);
	gsl_vector * x = gsl_vector_alloc(n);

	// Set random matrix and vector
	randomizer_matrix(S,30);
	gsl_matrix_memcpy(S_orthogonal,S);
	randomizer_vector(b,30);
	matrix_print("The random matrix S = ",S);
	printf("\n");
	vector_print("The random vetor b = ",b);
	printf("\n");

	// QR-decomposition
	GS_decomp(S_orthogonal,T);

	matrix_print("From QR-decomposition, the new orthogonal matrix is: ",S_orthogonal);
	printf("\n");
	matrix_print("and the triangular matrix is: ",T);
	printf("\n");

	// Solve for x in QRx=b, and calculate Sx (should be = b)
	GS_solve(S_orthogonal,T,b,x);
	gsl_blas_dgemv(CblasNoTrans,1,S,x,0,b);

	vector_print("Check if Sx=b for solution of QRx=b",b);
	printf("\n");

	gsl_matrix_free(S);
	gsl_matrix_free(S_orthogonal);
	gsl_vector_free(x);
	gsl_vector_free(b);
	gsl_matrix_free(T);
return 0;
}