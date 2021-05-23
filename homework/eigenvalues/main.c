#include <stdio.h>
#include <stdlib.h>
#include <float.h>
#include <math.h>
#include <time.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_eigen.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_linalg.h>
#include "matrix.h"

#ifndef NDEBUG
#define TRACE(args...) fprintf(stderr,args)
#else
#define TRACE(...)
#endif

void measuretime(int Nmax){
	FILE * outmeasuretime = fopen("outtime.txt","w");

	for (int n = 20; n <= Nmax; n++) {
		// own implementation
		gsl_matrix * S = gsl_matrix_alloc(n,n);
		gsl_matrix * Sgsl = gsl_matrix_alloc(n,n);
		gsl_vector * A = gsl_vector_alloc(n);
		gsl_matrix * V = gsl_matrix_alloc(n,n);
		gen_rel_sym(S,1);
		
		clock_t tic = clock();
			jacobi_diag(S,V);
		clock_t toc = clock();

		double myTime = (double)(toc-tic) / CLOCKS_PER_SEC;
		
		
		// gsl
		gen_rel_sym(Sgsl,1);
		gsl_eigen_symmv_workspace * w = gsl_eigen_symmv_alloc(n);
		
		clock_t gsltic = clock();
			gsl_eigen_symmv(Sgsl,A,V,w);
		clock_t gsltoc = clock();

		double gslTime = (double)(gsltoc-gsltic) / CLOCKS_PER_SEC;
		
		// own optimized implementation
		gsl_matrix * Soptim = gsl_matrix_alloc(n,n);
		gsl_matrix * Voptim = gsl_matrix_alloc(n,n);
		gen_rel_sym(Soptim,1);
		
		clock_t tico = clock();
			jacobi_diag_optim(Soptim,Voptim);
		clock_t toco = clock();

		double myTimeo = (double)(toco-tico) / CLOCKS_PER_SEC;






		fprintf(outmeasuretime,"%i %10g %10g %10g\n",n,myTime/pow(n,3),gslTime/pow(n,3),myTimeo/pow(n,3));

		gsl_matrix_free(S);
		gsl_matrix_free(Sgsl);
		gsl_vector_free(A);
		gsl_matrix_free(V);
		gsl_eigen_symmv_free(w);
		gsl_matrix_free(Soptim);
		gsl_matrix_free(Voptim);
	}

	fclose(outmeasuretime);
}

int main() {
	// Part A
	
	printf("Part A\n\n");

	int n = 5;	// size of square matrix
	
	// allocating matrices
	gsl_matrix * S = gsl_matrix_alloc(n,n);
	gsl_matrix * V = gsl_matrix_alloc(n,n);
	gsl_matrix * vdvt = gsl_matrix_alloc(n,n);
	gsl_matrix * D = gsl_matrix_alloc(n,n);
	gsl_matrix * vtsv = gsl_matrix_alloc(n,n);
	gsl_matrix * temp = gsl_matrix_alloc(n,n);

	gen_rel_sym(S,1);		// generates real symmetric matrix with entry values between -1 and 1
	matrix_print("We will work with the symmetric matrix S=",S);
	gsl_matrix_memcpy(D,S);		
	jacobi_diag(D,V);		// D will contain eigenvalues, V will contain corresponding eigenvectors
	printf("\nJacobi diagonalization...\n\n");
	matrix_print("Matrix of eigenvalues of S, D=",D);
	matrix_print("\nMatrix of eigenvectors of S, V=",V);
	gsl_blas_dgemm(CblasTrans,CblasNoTrans,1,V,S,0,temp);
	gsl_blas_dgemm(CblasNoTrans,CblasNoTrans,1,temp,V,0,vtsv);
	matrix_print("\nCheck that V^T*S*V = D: ",vtsv);
	gsl_blas_dgemm(CblasNoTrans,CblasNoTrans,1,V,D,0,temp);
	gsl_blas_dgemm(CblasNoTrans,CblasTrans,1,temp,V,0,vdvt);
	matrix_print("\nCheck that V*D*V^T = S: ",vdvt);
	gsl_blas_dgemm(CblasTrans,CblasNoTrans,1,V,V,0,temp);
	matrix_print("\nCheck that V^T*V = 1: ",temp);

	gsl_matrix_free(S);
	gsl_matrix_free(vdvt);
	gsl_matrix_free(D);
	gsl_matrix_free(V);
	gsl_matrix_free(vtsv);
	gsl_matrix_free(temp);

	// Part B

	printf("\n\nPart B\n\nThe results are given in the png files.\n\n\n");

	n=200;
	double s=1.0/(n+1);
	gsl_matrix * H = gsl_matrix_alloc(n,n);
	gsl_matrix * B = gsl_matrix_alloc(n,n);
	for(int i=0;i<n-1;i++){
  		gsl_matrix_set(H,i,i,-2);
		gsl_matrix_set(H,i,i+1,1);
		gsl_matrix_set(H,i+1,i,1);
  	}
	gsl_matrix_set(H,n-1,n-1,-2);
	gsl_matrix_scale(H,-pow(s,-2));

	// matrix_print("H = ",H);

	jacobi_diag(H,B);

	FILE * quantumenergyout = fopen("quantumenergyout.txt","w");

	for (int k=0; k < n; k++){
    		double exact = M_PI*M_PI*(k+1)*(k+1);
    		double calculated = gsl_matrix_get(H,k,k);
    		fprintf(quantumenergyout,"%i %g %g\n",k,calculated,exact);
	}

	fclose(quantumenergyout);

	FILE * eigenfuncs = fopen("eigenfuncs.txt","w");

	int flip[3] = {-1, -1, 1};	// we flip the sign of the eigenfunctions, since the jacobi procedure
					// has flipped their signs
	for(int k=0;k<3;k++){
		fprintf(eigenfuncs,"#INDEX %i\n",k);
  		fprintf(eigenfuncs,"%g %g %g\n",0.0,0.0,0.0);
  		for(int i=0;i<n;i++) {
			double xis = (i+1.0)/(n+1);
			double c = 1.0/sqrt(s);		// For every eigenvector u=(u0,u2,...), u^T*u = 1. But in our case, we must
							// normalize the eigenvectors such that the Riemann sum |u0|^2*s+|u1|^2*s+...=1.
							// Therefore we must multiply the vectors by c.
			fprintf(eigenfuncs,"%g %g %g\n",xis,c*flip[k]*gsl_matrix_get(B,i,k),
			sqrt(2)*sin((k+1)*M_PI*xis));
		}
  		fprintf(eigenfuncs,"%g %g %g\n\n\n",1.0,0.0,0.0);
  	}

	fclose(eigenfuncs);

	gsl_matrix_free(H);
	gsl_matrix_free(B);
	
	// Part C
	
	printf("Part C\n\n");
	int nmax = 100;
	measuretime(nmax);
	
	printf("Result is in outtime.png\n\n");

return 0;
}
