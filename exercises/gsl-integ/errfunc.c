#include<math.h>
#include<stdio.h>
#include<gsl/gsl_integration.h>

double f(double x, void* params) {
	// double z = *(double*) params;
	double f = 2 / sqrt(M_PI) * exp(-pow(x,2));
return f;
}

double erf(double z) {
	gsl_function F;
	F.function = &f;
	int limit = 1000; double result, error;
	gsl_integration_workspace* w = gsl_integration_workspace_alloc(limit);
  	gsl_integration_qags (&F, 0, z, 0, 1e-7, 1000, w, &result, &error);
  	gsl_integration_workspace_free (w);
return result;
}

int main() {
	double z0 = 1.4;
	for (double z=-z0;z <= z0;z += 1.0 / 20)
		printf("%10g %10g\n",z,erf(z));
return 0;
}
