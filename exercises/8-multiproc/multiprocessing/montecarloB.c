#include <stdio.h>
#include <stdlib.h>
#include <omp.h>
#include <math.h>

double randNumber(void) {unsigned int seed; double output = 1.0*rand_r(&seed)/RAND_MAX; return output;}

int main() {
	int Ntot = 1e8;
	double count1 = 0, count2 = 0, count3 = 0;
#pragma omp parallel sections
	{
	#pragma omp section
		{
		for(int i=1; i<Ntot; i++) {
			double x = randNumber(), y = randNumber();
			if (pow(x,2)+pow(y,2)<1) count1++;}
		}
	#pragma omp section
		{
		for(int i=1; i<Ntot; i++) {
			double x = randNumber(), y = randNumber();
			if (pow(x,2)+pow(y,2)<1) count2++;}
		}
	#pragma omp section
		{
		for(int i=1; i<Ntot; i++) {
			double x = randNumber(), y = randNumber();
			if (pow(x,2)+pow(y,2)<1) count3++;}
		}
	}
	printf("B: Pi value is %10g\n",4*(count1+count2+count3)/(3*Ntot));
return 0;
}
