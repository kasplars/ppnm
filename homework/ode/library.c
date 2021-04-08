#include <stdio.h>
#include <math.h>
#include <time.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_linalg.h>

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
	// double threshold = 10000*LDBL_EPSILON;
	double threshold = 10000000000*LDBL_EPSILON;
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
