#include<limits.h>
#include<float.h>
#include<stdio.h>
#include "part3.h"

void name_digit(int i){
	switch (i) {
		case 1: printf("one"); break;
		case 2: printf("two"); break;
		case 3: printf("three"); break;
		case 4: printf("four"); break;
		case 5: printf("five"); break;
		case 6: printf("six"); break;
		case 7: printf("seven"); break;
		case 8: printf("eight"); break;
		case 9: printf("nine"); break;
		default: printf("not a digit"); break;}
}

int main(){
	printf("Part 1\n\n");
	int i;
	// 1.i
	printf("Max integer is determined with while loop to be:\n");
	i=1; while(i+1>i) {i++;}
	printf("Max integer = %i\n\n",i);

	printf("Max integer is determined with for loop to be:\n");
	for (i = INT_MAX - 10;1<i;i++)
	printf("%i\n",i);
	printf("Max integer ^\n\n");

	printf("Max integer is determined with a do while loop to be:\n");
	i=1; do {i++;} while (i+1>1);
	printf("Max integer = %i\n\n",i);

	// 1.ii
	printf("Min integer is determined with while loop to be:\n");
	i=-1; while(i-1<1) {i--;}
	printf("Min integer = %i\n\n",i);

	printf("Min integer is determined with for loop to be:\n");
	for (i = INT_MIN + 10;1>i;i--)
	printf("%i\n",i);
	printf("Min integer ^\n\n");

	printf("Min integer is determined with a do while loop to be:\n");
	i=0; do {i--;} while (i-1<0);
	printf("Min integer = %i\n\n",i);

	// 1.iii
	// while loops
	printf("Calculating the machine epsilon for different types using while loops:\n\n");
	{double x=1; while(1+x!=1){x/=2;} x*=2;
	printf("my machine epsilon = %5.10g for double\n",x);
	}

	{long double x=1; while(1+x!=1){x/=2;} x*=2;
	printf("my machine epsilon = %5.10Lg for long double\n",x);
	}

	{float x=1.0; while(1+x!=1){x/=2;} x*=2;
	printf("my machine epsilon = %5.10f for float\n\n",x);
	}
	// for loops
	printf("Calculating the machine epsilon for different types using for loops\n\n");
	{double e; for(e=1; 1+e!=1; e/=2){} e*=2;
	printf("my machine epislon = %5.10g for double\n",e);
	}

	{long double e; for(e=1; 1+e!=1; e/=2){} e*=2;
	printf("my machine epislon = %5.10Lg for long double\n",e);
	}

	{float e; for(e=1; 1+e!=1; e/=2){} e*=2;
	printf("my machine epislon = %5.10f for float\n\n",e);
	}
	// do while
	printf("Calculating the machine epsilon for different types using do while loops\n\n");
	{double x=1; do {x/=2;} while(1+x!=1); x*=2;
	printf("my machine epsilon = %5.10g for double\n",x);
	}

	{long double x=1; do {x/=2;} while(1+x!=1); x*=2;
	printf("my machine epsilon = %5.10Lg for long double\n",x);
	}

	{float x=1.0; do {x/=2;} while(1+x!=1); x*=2;
	printf("my machine epsilon = %5.10f for float\n\n",x);
	}
	// Compared with built-in epsilons
	printf("Comparing with the built-in values:\n\n");
	printf("machine double epsilon = %5.10g\n",DBL_EPSILON);
	printf("machine long double epsilon = %5.10Lg\n",LDBL_EPSILON);
	printf("machine float epsilon = %5.10f\n\n",FLT_EPSILON);

	// Part 2
	printf("Part 2\n\n");
	// 2.i
	int max = INT_MAX / 2;

	float sum_up_float = 0.0f;
	for (int i=1; i<max; i++) {sum_up_float += 1.0f/i;}
	printf("The sum up float is %g\n",sum_up_float);

	float sum_down_float = 0.0f;
	for (int i=max; 0<i; i--) {sum_down_float += 1.0f/i;}
	printf("The sum down float is %g\n",sum_down_float);

	// 2.ii
	printf("Explanation for ^ is commented out.\n");
	/*
	The difference comes from the fact that FLT_EPSILON is relatively large, such that the precision
	is not quite limited. In the "up" case, we quickly arrive at a point where small numbers become
	'insignificant' such that addition by small numbers is basically addition with 0 (rounding down).
	In the "down" case, the numbers added to sum_down_float are growing larger and larger and are thus
	not becoming insignificant. Therefore sum_down_float > sum_up_float.
	*/

	// 2.iii
	/*
	Mathematically, it does not converge, since it is the harmonic series. Computationally however, the series will converge because of the finite precision, depending on max.
	*/

	// 2.iv

	double sum_up_double = 0.0;
	for (int i=1; i<max; i++) {sum_up_double += 1.0/i;}
	printf("The sum up double is %lg\n\n",sum_up_double);

	double sum_down_double = 0.0;
	for (int i=max; 0<i; i--) {sum_down_double += 1.0/i;}
	printf("The sum down double is %lg\n",sum_down_double);
	printf("Explanation for ^ is commented out.\n\n");
	/*
	Better precision with double.
	*/

	// Part 3
	printf("Part 3\n\n");
	printf("%i\n\n",equals(1,2,3,0));

	// Part 4
	printf("Part 4\n\n");
	name_digit(3); printf("\n");
return 0;}

