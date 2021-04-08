#ifndef LIBRARY_H
#define LIBRARY_H

void backsub(gsl_matrix * U, gsl_vector * c);
void vector_print(char s[],gsl_vector* v);
void matrix_print(char s[],gsl_matrix* A);
void randomizer_matrix(gsl_matrix * A, int modulo);
void randomizer_vector(gsl_vector * A, int modulo);
double normColumn(gsl_matrix * A, double j);
void GS_decomp(gsl_matrix * A, gsl_matrix * R);
void GS_solve(gsl_matrix* Q, gsl_matrix* R, gsl_vector* b, gsl_vector* x);
void GS_inverse(gsl_matrix * Q, gsl_matrix * R, gsl_matrix * B);
void timesJ(gsl_matrix* A, int p, int q, double theta);
void Jtimes(gsl_matrix* A, int p, int q, double theta);
double largestoffdiag(gsl_matrix * A);
void jacobi_diag(gsl_matrix * A,gsl_matrix * V);
void gen_rel_sym(gsl_matrix * S,double range);

#endif
