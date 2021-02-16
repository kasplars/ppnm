#include<limits.h>
#include<float.h>
#include<stdio.h>

int equals(double a, double b, double tau, double epsilon);

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
	printf("Part 1\n");
	/*
	// 1.i
	int i=1; while(i+1>i) {i++;}
	printf("my max int = %i\n",i);

	for (int i = INT_MAX - 10;1<i;i++)
	printf("my max int = %i\n",i);

	i=1; do {i++;} while (i>1)
	print("my max int = %i\n",i);

	// 1.ii
	int i=1; while(i+1>i) {i--;}
	printf("my min int = %i\n",i);

	for (int i = INT_MIN + 10;1>i;i--)
	printf("my min int = %i\n",i);

	i=1; do {i--;} while (i<1)
	print("my min int = %i\n",i);	
	*/

	// 1.iii
	// while loops
	printf("while loops\n");
	{double x=1; while(1+x!=1){x/=2;} x*=2;
	printf("my machine epsilon = %g for double\n",x);
	}

	{long double x=1; while(1+x!=1){x/=2;} x*=2;
	printf("my machine epsilon = %Lg for long double\n",x);
	}

	{float x=1.0; while(1+x!=1){x/=2;} x*=2;
	printf("my machine epsilon = %f for float\n\n",x);
	}
	// for loops
	printf("for loops\n");
	{double e; for(e=1; 1+e!=1; e/=2){} e*=2;
	printf("my machine epislon = %g for double\n",e);
	}

	{long double e; for(e=1; 1+e!=1; e/=2){} e*=2;
	printf("my machine epislon = %Lg for long double\n",e);
	}

	{float e; for(e=1; 1+e!=1; e/=2){} e*=2;
	printf("my machine epislon = %f for float\n\n",e);
	}
	// do while
	printf("do while loops\n");
	{double x=1; do {x/=2;} while(1+x!=1); x*=2;
	printf("my machine epsilon = %g for double\n",x);
	}

	{long double x=1; do {x/=2;} while(1+x!=1); x*=2;
	printf("my machine epsilon = %Lg for long double\n",x);
	}

	{float x=1.0; do {x/=2;} while(1+x!=1); x*=2;
	printf("my machine epsilon = %f for float\n\n",x);
	}
	// Compared with built-in epsilons
	printf("machine double epsilon = %g\n",DBL_EPSILON);
	printf("machine long double epsilon = %Lg\n",LDBL_EPSILON);
	printf("machine float epsilon = %f\n\n",FLT_EPSILON);

	// Part 2
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
	In the "up" case, there is a point where we add 0 all the time, so yes. If max did not have an
	upper limit, then the "down" case would not have a limit. I'm not really sure how to answer this.
	*/

	// 2.iv

	double sum_up_double = 0.0;
	for (int i=1; i<max; i++) {sum_up_double += 1.0/i;}
	printf("The sum up double is %lg\n",sum_up_double);

	double sum_down_double = 0.0;
	for (int i=max; 0<i; i--) {sum_down_double += 1.0/i;}
	printf("The sum down double is %lg\n",sum_down_double);
	printf("Explanation for ^ is commented out.\n");
	/*
	Better precision with double.
	*/

	// Part 3
	printf("%i\n",equals(1,2,3,0));

	// Part 4
	name_digit(3); printf("\n");
return 0;}

