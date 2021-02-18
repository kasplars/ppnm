#include<stdio.h>
#include<stdlib.h>
#include<math.h>

int main(int argc, char** argv){
	int items = 1;
	double x;
	while (items!=EOF) {
		items = scanf("%lg",&x);
		printf("%lg\n",x);
	}
return 0;
}


