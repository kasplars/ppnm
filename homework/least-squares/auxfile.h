#ifndef AUXFILE_H
#define AUXFILE_H

void backsub(gsl_matrix * U, gsl_vector * c);
void vector_print(char s[],gsl_vector* v);
void matrix_print(char s[],gsl_matrix* A);
void randomizer_matrix(gsl_matrix * A, int modulo);
void randomizer_vector(gsl_vector * A, int modulo);
double normColumn(gsl_matrix * A, double j);
void GS_decomp(gsl_matrix * A, gsl_matrix * R);

#endif
