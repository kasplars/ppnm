#include<stdio.h>
#include<math.h>

int main(){
	double x;
	int items;
	items = fscanf(stdin,"%lg",&x);
	while (items!=EOF) {	
		printf("x = %lg, sin(x) = %lg, cos(x) = %lg\n",x,sin(x),cos(x));
		items = fscanf(stdin,"%lg",&x);
	}
return 0;
}


