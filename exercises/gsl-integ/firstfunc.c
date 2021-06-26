#include<math.h>
#include<stdio.h>
#include<gsl/gsl_integration.h>

double f(double x,void* params) {
	double f = log(x) / sqrt(x);
return f;
}

int main() {
	int limit = 999;
	gsl_integration_workspace* w = gsl_integration_workspace_alloc(limit);
	double a = 0, b = 1;
	double result, error, acc = 1e-3, eps = 1e-3;
	gsl_function F;
	F.function = &f;
	gsl_integration_qags(&F,a,b,acc,eps,limit,w,&result,&error);
	gsl_integration_workspace_free(w);

	printf("result: %10g, error: %10g\n",result,error);
return 0;
}
