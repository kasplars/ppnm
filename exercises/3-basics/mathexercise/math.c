#include<stdio.h>
#include<math.h>
#include<complex.h>

int main(){
	double gammaval = tgamma(5);
	printf("Gamma(5) = %g\n",gammaval);

	double besselval = j0(0.5);
	printf("J1(0.5) = %g\n",besselval);

	complex z = csqrt(-2); 
return 0;
}
