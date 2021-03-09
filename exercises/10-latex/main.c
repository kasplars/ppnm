#include <stdio.h>
#include <math.h>

double ex(double x){
	if(x<0)return 1/ex(-x);
	if(x>1./8)return pow(ex(x/2),2);
return 1+x*(1+x/2*(1+x/3*(1+x/4*(1+x/5*
		(1+x/6*(1+x/7*(1+x/8*
		(1+x/9*(1+x/10)))))))));
}

int main(){
	int datapoints = 200;
	double xmax = 5;
	for (int i=0;i<datapoints;i++) {
		printf("%10g %10g %10g\n",
		xmax*i/datapoints,ex(xmax*i/
		datapoints),exp(xmax*i/datapoints));
	}
	return 0;
}
