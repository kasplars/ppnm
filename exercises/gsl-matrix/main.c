#include <stdio.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_eigen.h>

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

// Part A

int main(){
	int n = 3;

	gsl_matrix* A = gsl_matrix_alloc(n,n);
	gsl_matrix* Acopy = gsl_matrix_alloc(n,n);
	gsl_vector* b = gsl_vector_alloc(n);
	gsl_vector* x = gsl_vector_alloc(n);
	gsl_vector* y = gsl_vector_calloc(n);

	double matarray[9] = {6.13, -2.90, 5.86, 8.08, -6.31,
				-3.89, -4.36, 1.00, 0.19};
	double vecarray[3] = {6.23, 5.37, 2.29};

	// insert values in matrix
	gsl_matrix_set(A,0,0,matarray[0]);
	gsl_matrix_set(A,0,1,matarray[1]);
	gsl_matrix_set(A,0,2,matarray[2]);
	gsl_matrix_set(A,1,0,matarray[3]);
	gsl_matrix_set(A,1,1,matarray[4]);
	gsl_matrix_set(A,1,2,matarray[5]);
	gsl_matrix_set(A,2,0,matarray[6]);
	gsl_matrix_set(A,2,1,matarray[7]);
	gsl_matrix_set(A,2,2,matarray[8]);

	gsl_matrix_memcpy(Acopy,A);

	// insert values in vector
	gsl_vector_set(b,0,vecarray[0]);
	gsl_vector_set(b,1,vecarray[1]);
	gsl_vector_set(b,2,vecarray[2]);

	// solve
	gsl_linalg_HH_solve(Acopy,b,x);
	gsl_blas_dgemv(CblasNoTrans,1,A,x,0,y);

	vector_print("Solution to Ax=b is found using Householder solver: ",x);
	vector_print("Check solution by Ax. This should be b: ",y);

	gsl_matrix_free(A);
	gsl_matrix_free(Acopy);
	gsl_vector_free(b);
	gsl_vector_free(x);

// Part B

	n = 4;
	gsl_matrix* H = gsl_matrix_alloc(n,n);

	int i, j;
	for (i=0;i<n;i++) {
		for (j=0;j<n;j++) {
			double Hij = 1.0/(i+j+1);
			gsl_matrix_set(H,i,j,Hij);
		}
	}
	matrix_print("H = ",H);

	gsl_eigen_symmv_workspace * w = gsl_eigen_symmv_alloc(n);
	gsl_vector * eigenvalues = gsl_vector_alloc(n);
	gsl_matrix * eigenvectors = gsl_matrix_alloc(n,n);

	gsl_eigen_symmv(H,eigenvalues,eigenvectors,w);

	gsl_eigen_symmv_free(w);

	vector_print("The eigenvalues of H are: ",eigenvalues);
	matrix_print("The corresponding eigenvectors of H are",eigenvectors);

	gsl_vector_free(eigenvalues);
	gsl_matrix_free(eigenvectors);
	gsl_matrix_free(H);

return 0;
}


