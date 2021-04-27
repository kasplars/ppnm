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

void rkstep45(void (*f)(double, gsl_vector *, gsl_vector *), /* the f from dy/dt=f(t,y) */
	int n,				/* step number	*/
	double t,              		/* the current value of the variable */
	gsl_vector * yt,            	/* the current value y(t) of the sought function */
	double h,              		/* the step to be taken */
	gsl_vector * yh,             	/* output: y(t+h) */
	gsl_vector * dy             	/* output: error estimate */
){
	int i;
	gsl_vector * k1 = gsl_vector_alloc(n);
	gsl_vector * k2 = gsl_vector_alloc(n);
	gsl_vector * k3 = gsl_vector_alloc(n);
	gsl_vector * k4 = gsl_vector_alloc(n);
	gsl_vector * k5 = gsl_vector_alloc(n);
	gsl_vector * k6 = gsl_vector_alloc(n);
	gsl_vector * yn = gsl_vector_alloc(n);

	double b[7] = {0.0, 16.0/135, 0.0, 6656.0/12825, 28561.0/56430, -9.0/50, 2.0/55};
	double bstar[7] = {0.0, 25.0/216, 0.0, 1408.0/2565, 2197.0/4104, -1.0/5, 0.0};
	double c[7] = {0.0, 0.0, 1.0/4, 3.0/8, 12.0/13, 1.0, 1.0/2};
	double a1[7] = {0.0, 0.0, 1.0/4, 3.0/32, 1932.0/2197, 439.0/216, -8.0/27};
	double a2[7] = {0.0, 0.0, 0.0, 9.0/32, -7200.0/2197, -8.0, 2.0};
	double a3[7] = {0.0, 0.0, 0.0, 0.0, 7296.0/2197, 3680.0/513, -3544.0/2565};
	double a4[7] = {0.0, 0.0, 0.0, 0.0, 0.0, -845.0/4104, 1859.0/4104};
	double a5[7] = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, -11.0/40};

	f(t,yt,k1);
	for(i=0;i<n;i++) {gsl_vector_set(yn,i,gsl_vector_get(yt,i)+(a1[2]*gsl_vector_get(k1,i))*h);}
	f(t+c[2]*h,yn,k2);
	for(i=0;i<n;i++) {gsl_vector_set(yn,i,gsl_vector_get(yn,i)+(a1[3]*gsl_vector_get(k1,i)+a2[3]*gsl_vector_get(k2,i))*h);}
	f(t+c[3]*h,yn,k3);
	for(i=0;i<n;i++) {gsl_vector_set(yn,i,gsl_vector_get(yn,i)+(a1[4]*gsl_vector_get(k1,i)+a2[4]*gsl_vector_get(k2,i)+a3[4]*gsl_vector_get(k3,i))*h);}
	f(t+c[4]*h,yn,k4);
	for(i=0;i<n;i++) {gsl_vector_set(yn,i,gsl_vector_get(yn,i)+(a1[5]*gsl_vector_get(k1,i)+a2[5]*gsl_vector_get(k2,i)+a3[5]*gsl_vector_get(k3,i)+a4[5]*gsl_vector_get(k4,i))*h);}
	f(t+c[5]*h,yn,k5);
	for(i=0;i<n;i++) {gsl_vector_set(yh,i,gsl_vector_get(yn,i)+(a1[6]*gsl_vector_get(k1,i)+a2[6]*gsl_vector_get(k2,i)+a3[6]*gsl_vector_get(k3,i)+a4[6]*gsl_vector_get(k4,i)+a5[6]*gsl_vector_get(k5,i))*h);}
	f(t+c[6]*h,yh,k6);
	for(i=0;i<n;i++) {gsl_vector_set(yh,i,gsl_vector_get(yt,i)+(b[1]*gsl_vector_get(k1,i)
				+b[2]*gsl_vector_get(k2,i)+b[3]*gsl_vector_get(k3,i)
				+b[4]*gsl_vector_get(k4,i)+b[5]*gsl_vector_get(k5,i)+b[6]*gsl_vector_get(k6,i))*h);
			gsl_vector_set(yn,i,gsl_vector_get(yt,i)+(bstar[1]*gsl_vector_get(k1,i)
				+bstar[2]*gsl_vector_get(k2,i)+bstar[3]*gsl_vector_get(k3,i)
				+bstar[4]*gsl_vector_get(k4,i)+bstar[5]*gsl_vector_get(k5,i)+bstar[6]*gsl_vector_get(k6,i))*h);
			gsl_vector_set(dy,i,gsl_vector_get(yh,i)-gsl_vector_get(yn,i));
	}

	gsl_vector_free(k1);
	gsl_vector_free(k2);
	gsl_vector_free(k3);
	gsl_vector_free(k4);
	gsl_vector_free(k5);
	gsl_vector_free(k6);
	gsl_vector_free(yn);
}

void odeprint(FILE * output,double x, gsl_vector * y) {
	int n = y->size;
	if (n == 1) fprintf(output,"%10g %10g\n",x,gsl_vector_get(y,0));
	if (n == 2) fprintf(output,"%10g %10g %10g\n",x,gsl_vector_get(y,0),gsl_vector_get(y,1));
	if (n == 3) fprintf(output,"%10g %10g %10g %10g\n",x,gsl_vector_get(y,0),gsl_vector_get(y,1),gsl_vector_get(y,2));
	if (n == 4) fprintf(output,"%10g %10g %10g %10g %10g\n",x,gsl_vector_get(y,0),gsl_vector_get(y,1),gsl_vector_get(y,2),
								gsl_vector_get(y,3));
	if (n == 5) fprintf(output,"%10g %10g %10g %10g %10g %10g\n",x,gsl_vector_get(y,0),gsl_vector_get(y,1),
								gsl_vector_get(y,2),gsl_vector_get(y,3),gsl_vector_get(y,4));
	if (n == 6) fprintf(output,"%10g %10g %10g %10g %10g %10g %10g\n",x,gsl_vector_get(y,0),gsl_vector_get(y,1),
								gsl_vector_get(y,2),gsl_vector_get(y,3),gsl_vector_get(y,4),
								gsl_vector_get(y,5));
	if (n == 7) fprintf(output,"%10g %10g %10g %10g %10g %10g %10g %10g\n",x,gsl_vector_get(y,0),gsl_vector_get(y,1),
								gsl_vector_get(y,2),gsl_vector_get(y,3),gsl_vector_get(y,4),
								gsl_vector_get(y,5),gsl_vector_get(y,6));
	if (n == 8) fprintf(output,"%10g %10g %10g %10g %10g %10g %10g %10g %10g\n",x,gsl_vector_get(y,0),gsl_vector_get(y,1),
								gsl_vector_get(y,2),gsl_vector_get(y,3),gsl_vector_get(y,4),
								gsl_vector_get(y,5),gsl_vector_get(y,6),gsl_vector_get(y,7));
	if (n == 9) fprintf(output,"%10g %10g %10g %10g %10g %10g %10g %10g %10g %10g\n",x,gsl_vector_get(y,0),gsl_vector_get(y,1),
								gsl_vector_get(y,2),gsl_vector_get(y,3),gsl_vector_get(y,4),
								gsl_vector_get(y,5),gsl_vector_get(y,6),gsl_vector_get(y,7),gsl_vector_get(y,8));
	if (n == 10) fprintf(output,"%10g %10g %10g %10g %10g %10g %10g %10g %10g %10g %10g\n",x,gsl_vector_get(y,0),gsl_vector_get(y,1),
								gsl_vector_get(y,2),gsl_vector_get(y,3),gsl_vector_get(y,4),
								gsl_vector_get(y,5),gsl_vector_get(y,6),gsl_vector_get(y,7),gsl_vector_get(y,8),
								gsl_vector_get(y,9));
	if (n == 11) fprintf(output,"%10g %10g %10g %10g %10g %10g %10g %10g %10g %10g %10g %10g\n",x,gsl_vector_get(y,0),gsl_vector_get(y,1),
								gsl_vector_get(y,2),gsl_vector_get(y,3),gsl_vector_get(y,4),
								gsl_vector_get(y,5),gsl_vector_get(y,6),gsl_vector_get(y,7),gsl_vector_get(y,8),
								gsl_vector_get(y,9),gsl_vector_get(y,10));
	if (n == 12) fprintf(output,"%10g %10g %10g %10g %10g %10g %10g %10g %10g %10g %10g %10g %10g\n",x,gsl_vector_get(y,0),gsl_vector_get(y,1),
								gsl_vector_get(y,2),gsl_vector_get(y,3),gsl_vector_get(y,4),
								gsl_vector_get(y,5),gsl_vector_get(y,6),gsl_vector_get(y,7),gsl_vector_get(y,8),
								gsl_vector_get(y,9),gsl_vector_get(y,10),gsl_vector_get(y,11));
}

int driver(FILE * outstream,void (*f)(double t, gsl_vector * y,gsl_vector * dydt), /* right-hand-side of dy/dt=f(t,y) */
	int n,
	double a,                     /* the start-point a */
	gsl_vector * ya,                     /* y(a) */
	double b,                     /* the end-point of the integration */
	gsl_vector * yb,                     /* y(b) to be calculated */
	double h,                     /* initial step-size */
	double acc,                   /* absolute accuracy goal */
	double eps                    /* relative accuracy goal */
){
	double x = a, ei, normy, s, taui;

	int i, k = 0;
	int max = 1e6;

	gsl_vector * dy = gsl_vector_alloc(n);

	while(x<b) {
		if (x+h>b) h = b-x;
		rkstep45(f,n,x,ya,h,yb,dy);
		s=0; for (i=0; i<n; i++) s+=pow(gsl_vector_get(dy,i),2); ei = sqrt(s);
		s=0; for (i=0; i<n; i++) s+=pow(gsl_vector_get(yb,i),2); normy = sqrt(s);
		taui = (normy*eps+acc)*sqrt(h/(b-a));
		if (ei<taui){
			x += h; k++; if (k>max-1) return k; for(i=0;i<n;i++) {gsl_vector_memcpy(ya,yb);} odeprint(outstream,x,yb);
		} if (ei>0) h*=pow(taui/ei,0.25)*0.95; else h*=2;
	}

	gsl_vector_free(dy);
	return k;
}
