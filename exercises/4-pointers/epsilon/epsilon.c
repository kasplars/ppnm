#include<limits.h>
#include<float.h>
#include<stdio.h>

int main(){
	// 1.i
	for (int i = INT_MAX - 10;1<i;i++)
	printf("My max int = %i\n",i);
	// 1.ii


	// 2.i
	int max = INT_MAX / 3;
	float sum_up_float = 0.0f;
	for (int i=1; i<max; i++) {sum_up_float += 1.0f/i;}
	printf("The sum up float is %g\n",sum_up_float);

	float sum_down_float = 0.0f;
	for (int i=max; 0<i; i--) {sum_up_float += 1.0f/i;}
	printf("The sum down float is %g\n",sum_down_float);


return 0;}
