#include<stdio.h>
#include<math.h>
#include<gsl/gsl_sf_gamma.h>

double Gamma();

int main(){
	double xmin=0.1,xmax=4.5;
	for(double x=xmin;x<=xmax;x+=1.0/100) {
		printf("%10g %10g %10g %10g\n",x,tgamma(x),gsl_sf_gamma(x),Gamma(x));
	}
return 0;
}
