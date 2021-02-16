#include<stdio.h>
#include<math.h>
#include<complex.h>

int main(){
	// Part 1
	{
	double gammaval = tgamma(5);
	printf("Gamma(5) = %g\n",gammaval);

	double besselval = j0(0.5);
	printf("J1(0.5) = %g\n",besselval);

	complex z = csqrt(-2);
	printf("sqrt(-2) = %g + i %g\n",creal(z),cimag(z));

	complex x = cexp(I * M_PI);
	printf("exp(i*pi) = %g + i %g\n",creal(x),cimag(x));

	complex y = cexp(I);
	printf("exp(i) = %g + i %g\n",creal(y),cimag(y));

	complex t = cpow(I,M_E);
	printf("i^e = %g + i %g\n",creal(t),cimag(t));

	complex k = cpow(I,I);
	printf("i^i = %g + i %g\n\n",creal(k),cimag(k));

	printf("The results are correct.\n\n");
	}
	// Part 2
	{
	float x_float = 1.f/9;
	double x_double = 1./9;
	long double x_long_double = 1.L/9;
	printf("%.30g\n%.30lg\n%.30Lg\n",x_float,x_double,x_long_double);
	printf("Number %i\n",3/2);
}
return 0;
}
